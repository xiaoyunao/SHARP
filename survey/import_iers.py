#!/usr/bin/env python3
"""Validate a locally downloaded IERS-A table before updating Astropy's cache."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from astropy.time import Time
from astropy.utils import iers
from astropy.utils.data import cache_contents, import_file_to_cache


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file', type=Path)
    parser.add_argument('--date', required=True, help='Observing date, YYYY-MM-DD')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--backup-dir', type=Path,
                        default=Path('/pipeline/xiaoyunao/survey/runtime/iers_backups'))
    args = parser.parse_args()
    table = iers.IERS_A.read(args.file)
    times = Time([args.date + 'T00:00:00', args.date + 'T23:59:59'])
    prediction_mjd = float(table.meta['predictive_mjd'])
    age = float(max(times.mjd) - prediction_mjd)
    if age > 30 or prediction_mjd > min(times.mjd):
        raise ValueError(f'IERS prediction age is unsuitable for {args.date}: {age:.3f} days')
    _, ut1_status = table.ut1_utc(times, return_status=True)
    _, _, pm_status = table.pm_xy(times, return_status=True)
    if np.any(ut1_status < 0) or np.any(pm_status < 0):
        raise ValueError('IERS table does not cover the complete observing date')
    result = dict(date=args.date, rows=len(table), prediction_mjd=prediction_mjd,
                  prediction_age_days=age,
                  last_mjd=float(table['MJD'][-1].value),
                  sha256=hashlib.sha256(args.file.read_bytes()).hexdigest(),
                  dry_run=args.dry_run)
    if not args.dry_run:
        backup = args.backup_dir / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        backup.mkdir(parents=True)
        old = cache_contents()
        url = iers.IERS_A_URL
        if url in old:
            shutil.copy2(old[url], backup / 'finals2000A_previous.all')
        (backup / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
        import_file_to_cache(url, str(args.file), replace=True)
        cached = Path(cache_contents()[url])
        if hashlib.sha256(cached.read_bytes()).hexdigest() != result['sha256']:
            raise RuntimeError('Imported IERS cache checksum mismatch')
        result.update(cache=str(cached), backup=str(backup))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
