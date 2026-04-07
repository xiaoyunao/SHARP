#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import os
from pathlib import Path

from astropy.table import Table, vstack


def collect_tables(parts_dir: Path, suffix: str) -> list[Table]:
    tables: list[Table] = []
    for path in sorted(parts_dir.glob(f"*{suffix}")):
        tables.append(Table.read(path))
    return tables


def main() -> None:
    ap = argparse.ArgumentParser(description="Merge per-file known-asteroid FITS outputs into night-level FITS files.")
    ap.add_argument("night", help="Night name YYYYMMDD")
    ap.add_argument("--parts-dir", required=True, help="Directory holding per-file FITS parts")
    ap.add_argument("--outdir", required=True, help="Night L4 directory for merged outputs")
    args = ap.parse_args()

    parts_dir = Path(args.parts_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_tables = collect_tables(parts_dir, "_all_asteroids.fits")
    matched_tables = collect_tables(parts_dir, "_matched_asteroids.fits")

    out_all = outdir / f"{args.night}_all_asteroids.fits"
    out_matched = outdir / f"{args.night}_matched_asteroids.fits"

    if all_tables:
        vstack(all_tables).write(out_all, overwrite=True)
        print(f"[WRITE] {out_all} (parts={len(all_tables)})")
    else:
        print("[WRITE] skipped all_asteroids (no per-file results)")

    if matched_tables:
        vstack(matched_tables).write(out_matched, overwrite=True)
        print(f"[WRITE] {out_matched} (parts={len(matched_tables)})")
    else:
        print("[WRITE] skipped matched_asteroids (no per-file matches)")


if __name__ == "__main__":
    main()
