"""A common, content-validated science input set for known and unknown processing."""
from datetime import datetime, timedelta
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning

from .server_common import TZ, write_json, atomic_write, utc_now
from .unknown_pipeline.exposure_timing import header_mid_mjd


def inspect(path, night):
    path = Path(path)
    with fits.open(path, memmap=False, lazy_load_hdus=True) as hdus:
        if len(hdus) < 2 or not isinstance(hdus[1], fits.BinTableHDU):
            raise ValueError('not an L2 source catalog')
        header = hdus[0].header.copy()
        header.extend(hdus[1].header, update=True)
        name = str(header.get('OBJECT', ''))+' '+path.name
        kind = str(header.get('IMAGETYP', header.get('OBSTYPE', ''))).lower()
        if any(x in name.lower()+' '+kind for x in ('focus', 'bias', 'dark', 'flat','calib')):
            raise ValueError('calibration or focus exposure')
        # L2 source-table structure, a real celestial WCS and exposure timing
        # establish eligibility; a program-name prefix is not required.
        cols = set(hdus[1].columns.names)
        if not any({a,b} <= cols for a,b in (('RA_Win','DEC_Win'), ('RA_PSF','DEC_PSF'), ('A','B'))):
            raise ValueError('no recognized source sky coordinates')
        if 'objID' not in cols:
            raise ValueError('missing source identity objID')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',FITSFixedWarning)
            wcs = WCS(header)
        if not wcs.has_celestial or not np.isfinite(wcs.wcs.crval).all():
            raise ValueError('invalid celestial WCS')
        mid = header_mid_mjd(header)
        if mid is None:
            raise ValueError('missing or invalid actual exposure time')
        midtime = Time(mid, format='mjd', scale='utc')
        local = midtime.to_datetime(timezone=TZ)
        if local.hour < 12:
            local -= timedelta(days=1)
        if local.strftime('%Y%m%d') != night:
            raise ValueError('wrong_night:'+local.strftime('%Y%m%d'))
        seconds = None
        if header.get('EXPSTA') and header.get('EXPEND'):
            seconds = float((Time(header['EXPEND'])-Time(header['EXPSTA'])).sec)
        else:
            for key in ('EXPOSURE', 'EXPTIME'):
                try:
                    seconds = float(header[key]); break
                except (KeyError, ValueError, TypeError):
                    pass
        if seconds is None or not np.isfinite(seconds) or seconds <= 0:
            raise ValueError('invalid exposure duration')
        band = str(header.get('FILTER', header.get('CHANNEL', ''))).strip().lower()
        identity='|'.join((f'{mid:.10f}',f'{seconds:.4f}',band,
                           f'{wcs.wcs.crval[0]:.7f}',f'{wcs.wcs.crval[1]:.7f}'))
        uid = hashlib.sha256(identity.encode()).hexdigest()
        stat = path.stat()
        return dict(file=path.name, path=str(path.resolve()), exposure_uid=uid, mid_mjd=mid,
                    mid_unix=float(midtime.unix), exposure_s=seconds, filter=band,
                    ra_deg=float(wcs.wcs.crval[0]), dec_deg=float(wcs.wcs.crval[1]),
                    source_count=int(header.get('NAXIS2', 0)), bytes=stat.st_size,
                    mtime_ns=stat.st_mtime_ns, object=str(header.get('OBJECT','')))


def build(night, processed_root, output):
    from .server_quality import exclusion
    invalid = exclusion(processed_root, night)
    if invalid:
        result = dict(night=night, generated_at=utc_now(), accepted=[], rejected=[],
                      status='excluded_invalid_night', reason=invalid['reason'])
        write_json(output, result);atomic_write(Path(output).with_suffix('.txt'), '')
        return result
    folder = Path(processed_root)/night/'L2'
    accepted, rejected, seen = [], [], set()
    for path in sorted(folder.iterdir()) if folder.is_dir() else []:
        if not path.is_file() or not path.name.endswith(('.fits', '.fits.gz')):
            continue
        try:
            row = inspect(path, night)
            if row['exposure_uid'] in seen:
                raise ValueError('duplicate exposure identity')
            seen.add(row['exposure_uid']); accepted.append(row)
        except Exception as exc:
            rejected.append(dict(file=path.name, reason=str(exc)))
    signature = hashlib.sha256(json.dumps(accepted, sort_keys=True).encode()).hexdigest()
    result = dict(night=night, generated_at=utc_now(), input_signature=signature,
                  input_dir=str(folder), accepted=accepted, rejected=rejected,
                  status='ready' if accepted else 'no_valid_science_inputs')
    output = Path(output)
    write_json(output, result)
    atomic_write(output.with_suffix('.txt'), ''.join(r['file']+'\n' for r in accepted))
    return result


def load(path, verify=True):
    value = json.loads(Path(path).read_text())
    if verify:
        for row in value['accepted']:
            stat = Path(row['path']).stat()
            if (stat.st_size, stat.st_mtime_ns) != (row['bytes'], row['mtime_ns']):
                raise RuntimeError('Input changed since manifest: '+row['path'])
    return value
