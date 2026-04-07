#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Merge per-group 2-point tracklet FITS files into a single nightly FITS.

Designed to work with outputs from make_tracklet_linreproj.py:
  tracklets_<NIGHT>_gXXX.fits

Example:
  python merge_tracklets_night.py 20251206 \
    --tracklet-dir /pipeline/xiaoyunao/data/heliolincrr/20251206/tracklets_linreproj
"""

import os
import glob
import argparse
import shutil
from typing import List

import numpy as np
from astropy.table import Table, vstack


def read_one(fn: str) -> Table:
    """Read a single tracklet file. Returns an Astropy Table."""
    t = Table.read(fn)
    # ensure plain columns (avoid masked scalars downstream)
    return t


def merge_tables(tabs: List[Table]) -> Table:
    if len(tabs) == 1:
        return tabs[0]
    # vstack with exact column match
    return vstack(tabs, metadata_conflicts="silent")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("night", help="e.g. 20251206")
    ap.add_argument(
        "--tracklet-dir",
        default=None,
        help=(
            "Directory that contains tracklets_<night>_g*.fits. "
            "If not set, defaults to <root>/<night>/tracklets_linreproj"
        ),
    )
    ap.add_argument(
        "--root",
        default="/pipeline/xiaoyunao/data/heliolincrr",
        help="Root dir that contains <night>/tracklets_linreproj",
    )
    ap.add_argument(
        "--pattern",
        default=None,
        help="Glob pattern (default: tracklets_<night>_g*.fits)",
    )
    ap.add_argument(
        "--out",
        default=None,
        help="Output FITS path (default: <tracklet-dir>/tracklets_<night>_ALL.fits)",
    )
    ap.add_argument(
        "--add-source-col",
        action="store_true",
        help="Add a 'srcfile' column with the basename of the input file.",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output if exists.",
    )
    ap.add_argument(
        "--write-index",
        action="store_true",
        help=(
            "Also write a small text index listing all input files used: "
            "<out>.inputs.txt"
        ),
    )

    args = ap.parse_args()

    night = str(args.night)

    tracklet_dir = (
        args.tracklet_dir
        if args.tracklet_dir is not None
        else os.path.join(args.root, night, "tracklets_linreproj")
    )
    if args.pattern is None:
        pattern = os.path.join(tracklet_dir, f"tracklets_{night}_g*.fits")
    else:
        # allow relative or absolute patterns
        pattern = args.pattern
        if not os.path.isabs(pattern):
            pattern = os.path.join(tracklet_dir, pattern)

    out = (
        args.out
        if args.out is not None
        else os.path.join(tracklet_dir, f"tracklets_{night}_ALL.fits")
    )
    os.makedirs(os.path.dirname(out), exist_ok=True)

    fns = sorted(glob.glob(pattern))
    if len(fns) == 0:
        raise FileNotFoundError(
            f"No tracklet files matched pattern: {pattern}. "
            f"Check --tracklet-dir / --pattern."
        )

    tabs = []
    for fn in fns:
        t = read_one(fn)
        if args.add_source_col:
            # create as fixed-length string for FITS compatibility
            src = os.path.basename(fn)
            t["srcfile"] = np.array([src] * len(t), dtype=f"U{max(1, len(src))}")
        tabs.append(t)

    merged = merge_tables(tabs)

    # record provenance in metadata
    merged.meta["night"] = night
    merged.meta["n_inputs"] = int(len(fns))
    merged.meta["pattern"] = os.path.basename(pattern)

    tmp_out = out + ".tmp"
    merged.write(tmp_out, overwrite=True)
    shutil.move(tmp_out, out)

    if args.write_index:
        idx = out + ".inputs.txt"
        tmp_idx = idx + ".tmp"
        with open(tmp_idx, "w", encoding="utf-8") as f:
            for fn in fns:
                f.write(fn + "\n")
        shutil.move(tmp_idx, idx)

    print(f"Merged {len(fns)} files -> {out} (n_tracklets={len(merged)})")


if __name__ == "__main__":
    main()
