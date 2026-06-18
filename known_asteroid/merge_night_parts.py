#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import os
from pathlib import Path

from astropy.table import Table, vstack


def remove_stale(path: Path) -> None:
    if path.exists():
        path.unlink()
        print(f"[REMOVE] stale {path}")


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
    ap.add_argument("--mask-matched-suffix", default="_mask15", help="Optional wider matched output suffix to merge")
    args = ap.parse_args()

    parts_dir = Path(args.parts_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_tables = collect_tables(parts_dir, "_all_asteroids.fits")
    matched_tables = collect_tables(parts_dir, "_matched_asteroids.fits")
    mask_matched_tables = collect_tables(parts_dir, f"_matched_asteroids{args.mask_matched_suffix}.fits")

    out_all = outdir / f"{args.night}_all_asteroids.fits"
    out_matched = outdir / f"{args.night}_matched_asteroids.fits"
    out_mask_matched = outdir / f"{args.night}_matched_asteroids{args.mask_matched_suffix}.fits"

    if all_tables:
        vstack(all_tables).write(out_all, overwrite=True)
        print(f"[WRITE] {out_all} (parts={len(all_tables)})")
    else:
        remove_stale(out_all)
        print("[WRITE] skipped all_asteroids (no per-file results)")

    if matched_tables:
        vstack(matched_tables).write(out_matched, overwrite=True)
        print(f"[WRITE] {out_matched} (parts={len(matched_tables)})")
    else:
        remove_stale(out_matched)
        print("[WRITE] skipped matched_asteroids (no per-file matches)")

    if mask_matched_tables:
        vstack(mask_matched_tables).write(out_mask_matched, overwrite=True)
        print(f"[WRITE] {out_mask_matched} (parts={len(mask_matched_tables)})")
    else:
        remove_stale(out_mask_matched)
        print("[WRITE] skipped mask matched_asteroids (no per-file matches)")


if __name__ == "__main__":
    main()
