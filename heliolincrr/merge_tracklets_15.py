#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Merge nightly ALL tracklet FITS into a rolling window (default 15 nights).

Example:
  python merge_tracklets_window15.py 20251206 \
    --root /pipeline/xiaoyunao/data/heliolincrr \
    --ndays 15 \
    --nightly-subpath "tracklets_linreproj" \
    --nightly-name-template "tracklets_{night}_ALL.fits" \
    --outdir "/pipeline/xiaoyunao/data/heliolincrr/20251206/tracklets_windows_15"
"""

from __future__ import annotations

import argparse
import os
import json
import shutil
from pathlib import Path
from datetime import datetime, timedelta
from typing import List, Tuple

from astropy.table import Table, vstack

import warnings
from astropy.io.fits.verify import VerifyWarning

warnings.filterwarnings(
    "ignore",
    category=VerifyWarning
)


def yyyymmdd_to_date(s: str) -> datetime:
    return datetime.strptime(s, "%Y%m%d")


def date_to_yyyymmdd(d: datetime) -> str:
    return d.strftime("%Y%m%d")


def iter_nights(end_night: str, ndays: int) -> List[str]:
    end = yyyymmdd_to_date(end_night)
    start = end - timedelta(days=ndays - 1)
    nights = []
    d = start
    while d <= end:
        nights.append(date_to_yyyymmdd(d))
        d += timedelta(days=1)
    return nights


def merge_tables(tabs: List[Table]) -> Table:
    if len(tabs) == 1:
        return tabs[0]
    return vstack(tabs, metadata_conflicts="silent")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("end_night", help="YYYYMMDD (this is day 15 / last day)")
    ap.add_argument("--ndays", type=int, default=15, help="window length in nights (default 15)")
    ap.add_argument("--root", default="/pipeline/xiaoyunao/data/heliolincrr", help="data root")
    ap.add_argument("--nightly-subpath", default="tracklets_linreproj",
                    help='subpath under <root>/<night>/ (default: "tracklets_linreproj")')
    ap.add_argument("--nightly-name-template", default="tracklets_{night}_ALL.fits",
                    help='filename template under nightly-subpath')
    ap.add_argument("--outdir", default=None,
                    help="output directory")
    ap.add_argument("--overwrite", action="store_true", help="overwrite output file if exists")
    ap.add_argument("--write-index", action="store_true", help="write <out>.inputs.txt index")
    ap.add_argument("--require-all", action="store_true",
                    help="if set, require ALL nights exist; otherwise skip missing nights")
    args = ap.parse_args()

    end_night = str(args.end_night)
    nights = iter_nights(end_night, args.ndays)
    start_night = nights[0]

    outdir = Path(args.outdir) if args.outdir else (Path(args.root) / end_night / f"tracklets_windows_{args.ndays}")
    outdir.mkdir(parents=True, exist_ok=True)

    out = outdir / f"tracklets_{start_night}_{end_night}_W{args.ndays}.fits"

    tabs: List[Table] = []
    used_files: List[str] = []
    missing: List[str] = []

    for n in nights:
        fn = Path(args.root) / n / args.nightly_subpath / args.nightly_name_template.format(night=n)
        if fn.exists():
            t = Table.read(fn)
            # keep provenance
            t.meta["night"] = n
            tabs.append(t)
            used_files.append(str(fn))
        else:
            missing.append(n)
            if args.require_all:
                raise FileNotFoundError(f"Missing nightly ALL: {fn}")

    if len(tabs) == 0:
        raise RuntimeError(f"No nightly ALL tracklet files found in window {start_night}..{end_night}")

    merged = merge_tables(tabs)
    merged.meta["window_end"] = end_night
    merged.meta["window_start"] = start_night
    merged.meta["window_ndays"] = int(args.ndays)
    merged.meta["n_nights_found"] = int(len(tabs))
    merged.meta["n_nights_missing"] = int(len(missing))

    tmp_out = Path(str(out) + ".tmp")
    merged.write(tmp_out, overwrite=True)
    shutil.move(str(tmp_out), str(out))

    meta = {
        "window_start": start_night,
        "window_end": end_night,
        "window_ndays": int(args.ndays),
        "n_nights_found": int(len(tabs)),
        "n_nights_missing": int(len(missing)),
        "n_tracklets": int(len(merged)),
        "out_fits": str(out),
    }

    meta_path = str(out) + ".meta.json"
    tmp_meta = meta_path + ".tmp"
    with open(tmp_meta, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)
    shutil.move(tmp_meta, meta_path)

    print(f"[meta] wrote {meta_path}")
    print(f"[meta] n_nights_found={meta['n_nights_found']} n_tracklets={meta['n_tracklets']}")

    if args.write_index:
        idx = str(out) + ".inputs.txt"
        tmp_idx = idx + ".tmp"
        with open(tmp_idx, "w", encoding="utf-8") as f:
            for x in used_files:
                f.write(x + "\n")
            if missing:
                f.write("# missing_nights:\n")
                for n in missing:
                    f.write(n + "\n")
        shutil.move(tmp_idx, idx)

    print(f"Wrote {out}  n_tracklets={len(merged)}  nights={len(tabs)}  missing={len(missing)}")


if __name__ == "__main__":
    main()
