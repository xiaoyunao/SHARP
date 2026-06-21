#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from datetime import datetime, timezone

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
    out_status = outdir / f"{args.night}_known_asteroid_status.json"

    all_rows = 0
    matched_rows = 0
    mask_matched_rows = 0
    if all_tables:
        all_table = vstack(all_tables)
        all_rows = len(all_table)
        all_table.write(out_all, overwrite=True)
        print(f"[WRITE] {out_all} (parts={len(all_tables)})")
    else:
        remove_stale(out_all)
        print("[WRITE] skipped all_asteroids (no per-file results)")

    if matched_tables:
        matched_table = vstack(matched_tables)
        matched_rows = len(matched_table)
        matched_table.write(out_matched, overwrite=True)
        print(f"[WRITE] {out_matched} (parts={len(matched_tables)})")
    else:
        remove_stale(out_matched)
        print("[WRITE] skipped matched_asteroids (no per-file matches)")

    if mask_matched_tables:
        mask_matched_table = vstack(mask_matched_tables)
        mask_matched_rows = len(mask_matched_table)
        mask_matched_table.write(out_mask_matched, overwrite=True)
        print(f"[WRITE] {out_mask_matched} (parts={len(mask_matched_tables)})")
    else:
        remove_stale(out_mask_matched)
        print("[WRITE] skipped mask matched_asteroids (no per-file matches)")

    status = {
        "night": args.night,
        "known_complete": True,
        "updated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "parts_dir": str(parts_dir),
        "all_asteroids": {
            "path": str(out_all),
            "exists": out_all.exists(),
            "parts": len(all_tables),
            "rows": all_rows,
        },
        "matched_asteroids": {
            "path": str(out_matched),
            "exists": out_matched.exists(),
            "parts": len(matched_tables),
            "rows": matched_rows,
            "sep_arcsec": 1.0,
        },
        "mask_matched_asteroids": {
            "path": str(out_mask_matched),
            "exists": out_mask_matched.exists(),
            "parts": len(mask_matched_tables),
            "rows": mask_matched_rows,
            "sep_arcsec": 1.5,
            "suffix": args.mask_matched_suffix,
        },
    }
    out_status.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[WRITE] {out_status}")


if __name__ == "__main__":
    main()
