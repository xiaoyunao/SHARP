#!/usr/bin/env python3
"""Read-only server audit. Run with server Python; JSON is written to stdout.

ssh <existing connection options> <host> <python> - < this_file.py > evidence.json
No pipeline imports, writes, fitting, or reference catalog updates occur.
"""
import ast
import hashlib
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.time import Time

warnings.filterwarnings('ignore')
BASE = Path('/home/admin/pipeline/smtpipeline')
REF = Path('/pipeline/ref/healpix')
SAMPLES = [('20251121', 'OBJ_MP_0004_0135'),
           ('20260218', 'OBJ_MP_1528_0435'),
           ('20260508', 'OBJ_MP_1528_0049')]
FILES = ['utils.py', 'main.py', 'parameters.py', 'astrometry/cal_pos.py',
         'astrometry/config.scamp', 'calib_flux.py', 'photometry/phot.py',
         'photometry/config1.sex', 'ccd9k.py']


def fingerprint(p):
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
            'bytes': p.stat().st_size,
            'mtime_utc': datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat()}


def function_source(p, name):
    s = p.read_text()
    for n in ast.walk(ast.parse(s)):
        if isinstance(n, ast.FunctionDef) and n.name == name:
            return ast.get_source_segment(s, n)
    raise ValueError(name)


def run():
    source = function_source(BASE/'utils.py', 'propagate_gaiacat')
    namespace = dict(np=np, Time=Time, SkyCoord=SkyCoord, u=u)
    exec(compile(source, '<audited historical function>', 'exec'), namespace)
    report = {'code': [fingerprint(BASE/f) for f in FILES],
              'function_source': source,
              'call_chain_sources': {f'{f}:{name}': function_source(BASE/f, name)
                                     for f, name in [('calib_flux.py', 'gaia_mags'),
                                                     ('calib_flux.py', 'read_gaiacat'),
                                                     ('astrometry/cal_pos.py', 'get_gaiacat'),
                                                     ('astrometry/cal_pos.py', 'run_scamp')]},
              'alias_target': str(Path('/pipeline/smtpipeline').resolve()),
              'samples': [], 'reference_tiles': {}}
    for night, stem in SAMPLES:
        p = Path('/processed1')/night/'L2'/f'{stem}_cat.fits.gz'
        meta = fingerprint(p)
        with fits.open(p) as h:
            header = h[1].header
            cat = Table(h[1].data)
        keys = ['DATE-OBS', 'MJD', 'EXPOSURE', 'EXPTIME', 'CALI_REF', 'CRVAL1', 'CRVAL2',
                'CCDZP', 'APER08', 'APER18', 'MLIM', 'NAXIS1', 'NAXIS2']
        meta['header'] = {k: header.get(k) for k in keys}
        # Frozen L2 exposes DATE-OBS or its equivalent OBS_DATE in primary header.
        l1 = p.parent.parent/'L1'/f'{stem}.fits.gz'
        l1_header = fits.getheader(l1, 1)
        meta['l1_date_obs'] = l1_header.get('DATE-OBS')
        datestr = l1_header['DATE-OBS']
        t = Time(datestr, format='isot', scale='utc')
        meta['date_obs_utc'] = t.isot
        config = Path('/processed1')/night/'parameter.json'
        conf = json.loads(config.read_text())
        meta['night_configuration'] = {k: conf[k] for k in ['ASTROMETRY', 'FLUX CALIBRATION']}
        meta['night_configuration_sha256'] = fingerprint(config)['sha256']
        vec = hp.ang2vec(header['CRVAL1'], header['CRVAL2'], lonlat=True)
        tiles = hp.query_disc(32, vec, np.radians(3.0))
        refs = []
        for tile in tiles:
            q = REF/f'healpix-{tile:05d}.fits'
            if not q.exists():
                continue
            with fits.open(q, memmap=True) as h:
                tab = h[1].data
                use = (tab['phot_g_mean_mag'] >= 12) & (tab['phot_g_mean_mag'] < 17)
                refs.append(Table(tab[use]))
            if q.name not in report['reference_tiles']:
                report['reference_tiles'][q.name] = fingerprint(q)
        ref = vstack(refs)
        meta['tile_ref_epochs'] = np.unique(ref['ref_epoch']).tolist()
        # Execute the exact extracted function on the same columns and UTC start time.
        namespace['propagate_gaiacat'](ref, t)
        cols = cat.colnames
        ra_col = next(x for x in cols if x.lower() == 'ra_win')
        dec_col = next(x for x in cols if x.lower() == 'dec_win')
        obs = SkyCoord(np.asarray(cat[ra_col])*u.deg, np.asarray(cat[dec_col])*u.deg)
        moved = SkyCoord(ref['X_WORLD']*u.deg, ref['Y_WORLD']*u.deg)
        idx, sep, _ = moved.match_to_catalog_sky(obs)
        pm = np.hypot(ref['pmra'], ref['pmdec'])
        finite = np.isfinite(pm)
        use = (sep.arcsec < 1) & finite
        chosen = np.flatnonzero(use)[np.argsort(np.asarray(pm[use]))[-4:]]
        missing = np.flatnonzero((sep.arcsec < 1) & ~finite)[:2]
        chosen = np.r_[chosen, missing]
        old = SkyCoord(ref['ra']*u.deg, ref['dec']*u.deg)
        rows = []
        for i in chosen:
            row = {name: (int(ref[name][i]) if name == 'source_id' else float(ref[name][i]))
                   for name in ['source_id', 'ra', 'dec', 'ref_epoch', 'pmra', 'pmdec', 'X_WORLD', 'Y_WORLD']}
            row.update(night=night, exposure=stem, date_obs_utc=t.isot,
                       target_jyear=float(t.jyear), reference_tile=f'healpix-{hp.ang2pix(32, row["ra"], row["dec"], lonlat=True):05d}.fits',
                       shift_arcsec=float(old[i].separation(moved[i]).arcsec),
                       l2_ra=float(obs[idx[i]].ra.deg), l2_dec=float(obs[idx[i]].dec.deg),
                       l2_to_propagated_arcsec=float(sep[i].arcsec),
                       l2_to_original_arcsec=float(obs[idx[i]].separation(old[i]).arcsec),
                       missing_pm_fallback=bool(not finite[i]))
            rows.append(row)
        meta['reference_examples'] = rows
        meta['bright_reference_matches_within_1arcsec'] = int(use.sum())
        meta['median_l2_to_propagated_arcsec'] = float(np.median(sep.arcsec[use]))
        meta['median_l2_to_original_arcsec'] = float(np.median(obs[idx[use]].separation(old[use]).arcsec))
        report['samples'].append(meta)
    def clean(value):
        if isinstance(value, dict):
            return {k: clean(v) for k, v in value.items()}
        if isinstance(value, list):
            return [clean(v) for v in value]
        return None if isinstance(value, float) and not np.isfinite(value) else value
    print(json.dumps(clean(report), indent=2, allow_nan=False))


if __name__ == '__main__':
    run()
