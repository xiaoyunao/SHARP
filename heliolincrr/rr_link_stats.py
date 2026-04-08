#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.table import Table


def det_key(file_value, obj_id) -> str:
    return f"{Path(str(file_value)).name}:{int(obj_id)}"


def summarize_int_array(arr) -> dict[str, object]:
    arr = np.asarray(arr, dtype=np.int64)
    if arr.size == 0:
        return {
            "min": 0,
            "median": 0.0,
            "max": 0,
            "p90": 0.0,
            "counts": {},
        }
    vals, cnts = np.unique(arr, return_counts=True)
    return {
        "min": int(arr.min()),
        "median": float(np.median(arr)),
        "max": int(arr.max()),
        "p90": float(np.percentile(arr, 90)),
        "counts": {str(int(v)): int(c) for v, c in zip(vals, cnts)},
    }


def build_known_coverage_summary(tracklets_path: Path, rr_dir: Path, matched_path: Path, mask_dir: Path) -> dict[str, object]:
    tracklets = Table.read(tracklets_path)
    members = Table.read(rr_dir / "linkage_members.fits")
    matched = Table.read(matched_path)

    mask_keys = set()
    for fn in sorted(mask_dir.glob("OBJ_MP_*_cat.fits.gz")):
        try:
            t = Table.read(fn)
        except Exception:
            continue
        if "objID" not in t.colnames:
            continue
        for obj_id in np.asarray(t["objID"], dtype=np.int64):
            mask_keys.add(f"{fn.name}:{int(obj_id)}")

    key_to_tracklets = defaultdict(list)
    for row in tracklets:
        tid = str(row["tracklet_id"])
        key_to_tracklets[det_key(row["file1"], row["objID1"])].append(tid)
        key_to_tracklets[det_key(row["file2"], row["objID2"])].append(tid)

    tracklet_to_links = defaultdict(list)
    for lid, tid in zip(np.asarray(members["linkage_id"], dtype=int), np.asarray(members["tracklet_id"]).astype("U64")):
        tracklet_to_links[str(tid)].append(int(lid))

    keys = np.asarray([det_key(f, o) for f, o in zip(matched["source_file"], matched["objID"])], dtype="U256")
    survives_mask = np.asarray([k in mask_keys for k in keys], dtype=bool)
    in_tracklet = np.zeros(len(keys), dtype=bool)
    in_rr_link = np.zeros(len(keys), dtype=bool)
    n_rr_links = np.zeros(len(keys), dtype=np.int64)

    for i, key in enumerate(keys):
        tids = sorted(set(key_to_tracklets.get(str(key), [])))
        if not tids:
            continue
        in_tracklet[i] = True
        lids = sorted(set(lid for tid in tids for lid in tracklet_to_links.get(str(tid), [])))
        if lids:
            in_rr_link[i] = True
            n_rr_links[i] = len(lids)

    n_mask = int(np.sum(survives_mask))
    n_track = int(np.sum(in_tracklet))
    n_rr = int(np.sum(in_rr_link))
    return {
        "matched_detections_total": int(len(keys)),
        "survives_mask_gaia": n_mask,
        "in_tracklet": n_track,
        "in_rr_link": n_rr,
        "tracklet_given_mask": float(n_track / n_mask) if n_mask else 0.0,
        "rr_given_mask": float(n_rr / n_mask) if n_mask else 0.0,
        "rr_given_tracklet": float(n_rr / n_track) if n_track else 0.0,
        "n_rr_links_all": summarize_int_array(n_rr_links),
        "n_rr_links_hits_only": summarize_int_array(n_rr_links[in_rr_link]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("night", help="YYYYMMDD")
    ap.add_argument("--root-raw", default="/processed1")
    ap.add_argument("--root-out", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--rr-subdir", default="rr_links")
    ap.add_argument("--tracklets", default=None)
    ap.add_argument("--out", default=None, help="optional JSON summary output path")
    args = ap.parse_args()

    night = str(args.night)
    raw_root = Path(args.root_raw)
    out_root = Path(args.root_out)
    rr_dir = out_root / night / args.rr_subdir
    links_path = rr_dir / "links_tracklets.fits"
    members_path = rr_dir / "linkage_members.fits"
    tracklets_path = Path(args.tracklets) if args.tracklets else (
        out_root / night / "tracklets_linreproj" / f"tracklets_{night}_ALL.fits"
    )

    links = Table.read(links_path)
    members = Table.read(members_path)
    summary = {
        "night": night,
        "paths": {
            "rr_dir": str(rr_dir),
            "links_path": str(links_path),
            "members_path": str(members_path),
            "tracklets_path": str(tracklets_path),
        },
        "link_counts": {
            "n_links": int(len(links)),
            "n_member_rows": int(len(members)),
            "n_tracklets_per_link": summarize_int_array(links["n_tracklets"] if len(links) else []),
            "n_nights_per_link": summarize_int_array(links["n_nights"] if len(links) else []),
        },
    }

    matched_path = raw_root / night / "L4" / f"{night}_matched_asteroids.fits"
    mask_dir = out_root / night / "mask_gaia"
    if tracklets_path.exists() and matched_path.exists() and mask_dir.exists():
        summary["known_asteroid_coverage"] = build_known_coverage_summary(
            tracklets_path=tracklets_path,
            rr_dir=rr_dir,
            matched_path=matched_path,
            mask_dir=mask_dir,
        )

    out_path = Path(args.out) if args.out else (out_root / night / "analysis" / f"{night}_rr_link_stats.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"[write] {out_path}")


if __name__ == "__main__":
    main()
