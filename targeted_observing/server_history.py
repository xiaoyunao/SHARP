"""Actual exposure/detection ledger and cadence-aware TJN sampling state."""
from collections import defaultdict
from datetime import date, datetime
import json
import csv
import io
from pathlib import Path
import sqlite3

import numpy as np
from astropy.table import Table

from .catalog import BASE
from .server_common import runtime, write_json, utc_now, atomic_write
from .science_manifest import inspect


def connect(config):
    path = runtime(config)/'history.sqlite'
    path.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path, timeout=60)
    db.execute('PRAGMA journal_mode=WAL')
    db.executescript('''
      CREATE TABLE IF NOT EXISTS exposures (
        uid TEXT PRIMARY KEY, night TEXT, file TEXT, mid REAL, band TEXT,
        seconds REAL, ra REAL, dec REAL, object TEXT);
      CREATE TABLE IF NOT EXISTS detections (
        uid TEXT, object_id TEXT, mag_error REAL, good INTEGER,
        PRIMARY KEY(uid,object_id));
      CREATE TABLE IF NOT EXISTS ingested (night TEXT PRIMARY KEY, signature TEXT);
    ''')
    return db


def ingest(config, nights=None):
    """Read actual L2 and known detections; never add planned observations."""
    root = Path(config['processed_root'])
    if nights is None:
        nights = [p.name for p in root.iterdir() if p.is_dir() and p.name.isdigit()
                  and len(p.name)==8 and p.name >= config['history_start']]
    db = connect(config); updated = []
    for night in sorted(nights):
        l2 = root/night/'L2'
        if not l2.is_dir(): continue
        files = sorted(p for p in l2.iterdir() if p.name.endswith(('.fits','.fits.gz')))
        matched = root/night/'L4'/f'{night}_matched_asteroids.fits'
        signature = str([(p.name,p.stat().st_size,p.stat().st_mtime_ns) for p in files]
                        +[(str(matched),matched.stat().st_mtime_ns if matched.exists() else 0)])
        previous = db.execute('SELECT signature FROM ingested WHERE night=?',(night,)).fetchone()
        if previous and previous[0] == signature: continue
        frame_map = {}
        # Replacement/reprocessing of a night must remove its superseded state.
        # Deduplication by exposure UID also rejects the same exposure in two directories.
        old = [r[0] for r in db.execute('SELECT uid FROM exposures WHERE night=?',(night,))]
        db.executemany('DELETE FROM detections WHERE uid=?',[(x,) for x in old])
        db.execute('DELETE FROM exposures WHERE night=?',(night,))
        for p in files:
            try: r = inspect(p,night)
            except Exception: continue
            db.execute('INSERT OR REPLACE INTO exposures VALUES (?,?,?,?,?,?,?,?,?)',
                       (r['exposure_uid'],night,p.name,r['mid_unix'],r['filter'],r['exposure_s'],
                        r['ra_deg'],r['dec_deg'],r['object']))
            frame_map[p.name] = r
        if matched.exists():
            tab = Table.read(matched)
            for row in tab:
                if 'source_file' not in tab.colnames or str(row['source_file']).strip() not in frame_map: continue
                if 'number' not in tab.colnames: continue
                try:
                    number=row['number']
                    key=(str(int(number)) if not np.ma.is_masked(number) and np.isfinite(float(number))
                         else str(row['name']).strip().strip('()') if 'name' in tab.colnames else '')
                    if not key:continue
                    err = float(row['MagErr_PSF']) if not np.ma.is_masked(row['MagErr_PSF']) else float('nan')
                    good = np.isfinite(err) and 0 < err <= config['history_quality_error_mag']
                    if 'Flag' in tab.colnames: good &= not np.ma.is_masked(row['Flag']) and int(row['Flag']) == 0
                    if 'Flag_ISO_Num' in tab.colnames: good &= not np.ma.is_masked(row['Flag_ISO_Num']) and int(row['Flag_ISO_Num']) == 0
                    r = frame_map[str(row['source_file']).strip()]
                    if 'epoch' in tab.colnames and not np.ma.is_masked(row['epoch']):
                        if abs(float(row['epoch'])-r['mid_mjd'])*86400>max(60,r['exposure_s']/2):
                            continue
                    # Colour-dependent samples are retained, but the r queue
                    # never counts a g/i/z measurement as an r-band point.
                    db.execute('INSERT OR REPLACE INTO detections VALUES (?,?,?,?)',
                               (r['exposure_uid'],key,err if np.isfinite(err) else None,int(good)))
                except (ValueError,TypeError,OverflowError,np.ma.MaskError): continue
        db.execute('INSERT OR REPLACE INTO ingested VALUES (?,?)',(night,signature))
        db.commit(); updated.append(dict(night=night,science_exposures=len(frame_map)))
    db.close()
    return updated


def sampling(rows, config, periods=()):
    """Sampling readiness only: no claim of photometric calibration or unique period."""
    selected = sorted((r for r in rows if r['good'] and r['band']=='r'),key=lambda r:r['mid'])
    times = np.array([r['mid'] for r in selected]); nightly = defaultdict(list)
    for r in selected: nightly[r['night']].append(r['mid'])
    good_nights, spans, gaps = [], [], []
    for night, ts in sorted(nightly.items()):
        ts = np.unique(ts); span = float(np.ptp(ts)); delta = np.diff(ts)
        spans.append(span)
        if len(delta): gaps.extend(delta.tolist())
        # A burst of frames never satisfies a night spanning at least 30 minutes.
        segments=np.split(ts,np.flatnonzero(delta>config['tjn_max_sampling_gap_s'])+1)
        if any(len(part)>=config['tjn_min_points_per_night'] and np.ptp(part)>=config['tjn_min_night_span_s']
               for part in segments):
            good_nights.append(night)
    baseline = ((datetime.strptime(good_nights[-1],'%Y%m%d')-
                 datetime.strptime(good_nights[0],'%Y%m%d')).days if good_nights else 0)
    coverage = [len(set(np.floor((times/(3600*p)%1)*config['tjn_phase_bins']).astype(int)))/config['tjn_phase_bins']
                for p in periods if p>0] if len(times) else []
    median = float(np.median([r['error'] for r in selected])) if selected else None
    required_baseline=max(config['tjn_required_baseline_days'],3*max(periods)/24 if periods else 0)
    phase_bins={str(p):sorted(set(np.floor((times/(3600*p)%1)*config['tjn_phase_bins']).astype(int).tolist()))
                for p in periods if p>0}
    ready = (len(times)>=config['tjn_required_points'] and len(good_nights)>=config['tjn_required_nights']
             and baseline>=required_baseline and median is not None
             and median<=config['history_quality_median_mag'] and
             (not coverage or min(coverage)>=config['tjn_phase_coverage']))
    return dict(good_points=len(times),good_nights=len(good_nights),baseline_days=baseline,
                night_spans_s=spans,median_gap_s=float(np.median(gaps)) if gaps else None,
                largest_gap_s=float(max(gaps)) if gaps else None,median_error_mag=median,
                min_prior_phase_coverage=min(coverage) if coverage else None,
                required_baseline_days=required_baseline,prior_phase_bins=phase_bins,
                sampling_ready=bool(ready),last_mid_unix=float(times[-1]) if len(times) else None,
                calibrated_lightcurve_complete=None,period_solution_complete=None)


def states(config, periods=None):
    db=connect(config); objects=defaultdict(list)
    measurements=[]
    for mid,night,band,key,error,good,uid,file in db.execute('''
      SELECT e.mid,e.night,e.band,d.object_id,d.mag_error,d.good,e.uid,e.file
      FROM exposures e JOIN detections d ON e.uid=d.uid'''):
        if periods is not None and key not in periods:continue
        objects[key].append(dict(mid=mid,night=night,band=band,error=error,good=good))
        measurements.append((uid,night,file,key,mid,band,error,good))
    stream=io.StringIO();writer=csv.writer(stream);writer.writerow(['exposure_uid','night','file','object_id','mid_unix','filter','mag_error','good'])
    writer.writerows(measurements);atomic_write(runtime(config)/'history_tjn_measurements.csv',stream.getvalue())
    stream=io.StringIO();writer=csv.writer(stream);writer.writerow(['exposure_uid','night','file','mid_unix','filter','exposure_s','ra_deg','dec_deg','object'])
    writer.writerows(db.execute('SELECT * FROM exposures ORDER BY mid'));atomic_write(runtime(config)/'history_exposures.csv',stream.getvalue())
    db.close()
    result={k:sampling(rows,config,(periods or {}).get(k,())) for k,rows in objects.items()}
    write_json(runtime(config)/'tjn_sampling_status.json',dict(updated_at=utc_now(),objects=result,
        note='Actual matched detections and r-band sampling; sampling_ready is not a validated rotation period or calibrated lightcurve.'))
    return result
