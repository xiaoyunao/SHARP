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


def build_key_array(tab: Table, file_col: str, obj_col: str) -> np.ndarray:
    return np.asarray([det_key(f, o) for f, o in zip(tab[file_col], tab[obj_col])], dtype="U256")


def string_list_column(values_by_row, width: int = 2048) -> np.ndarray:
    return np.asarray(
        [";".join(map(str, vals)) if vals else "" for vals in values_by_row],
        dtype=f"U{width}",
    )


def summarize_int_array(arr: np.ndarray) -> dict[str, object]:
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


def load_tracklet_detection_members(tracklets_path: Path):
    if not tracklets_path.exists():
        return {}, {}

    t = Table.read(tracklets_path)
    if not {"tracklet_id", "file1", "file2", "objID1", "objID2"}.issubset(t.colnames):
        return {}, {}

    key_to_tracklets = defaultdict(list)
    key_to_det_rows = defaultdict(list)
    tids = np.asarray(t["tracklet_id"]).astype("U64")
    k1 = build_key_array(t, "file1", "objID1")
    k2 = build_key_array(t, "file2", "objID2")

    for i, tid in enumerate(tids):
        key_to_tracklets[str(k1[i])].append(str(tid))
        key_to_tracklets[str(k2[i])].append(str(tid))
        key_to_det_rows[str(k1[i])].append((str(tid), 1))
        key_to_det_rows[str(k2[i])].append((str(tid), 2))

    return key_to_tracklets, key_to_det_rows


def load_rr_members(rr_dir: Path):
    links_path = rr_dir / "links_tracklets.fits"
    members_path = rr_dir / "linkage_members.fits"
    if not links_path.exists() or not members_path.exists():
        return {}, {}

    members = Table.read(members_path)
    if len(members) == 0:
        return {}, {}

    tracklet_to_links = defaultdict(list)
    for lid, tid in zip(np.asarray(members["linkage_id"], dtype=int),
                        np.asarray(members["tracklet_id"]).astype("U64")):
        tracklet_to_links[str(tid)].append(int(lid))

    link_summary = {}
    orbit_path = rr_dir / "orbit_confirm" / "orbit_links.fits"
    if orbit_path.exists():
        orbit = Table.read(orbit_path)
        if len(orbit) > 0:
            for row in orbit:
                link_summary[int(row["linkage_id"])] = {
                    "fit_ok": bool(row["fit_ok"]) if "fit_ok" in orbit.colnames else False,
                    "is_good": bool(row["is_good"]) if "is_good" in orbit.colnames else False,
                }

    return tracklet_to_links, link_summary


def load_mask_gaia_keys(mask_dir: Path):
    if not mask_dir.exists():
        return set()

    keys = set()
    for fn in sorted(mask_dir.glob("OBJ_MP_*_cat.fits.gz")):
        try:
            t = Table.read(fn)
        except Exception:
            continue
        if "objID" not in t.colnames:
            continue
        base = fn.name
        for obj_id in np.asarray(t["objID"], dtype=np.int64):
            keys.add(f"{base}:{int(obj_id)}")
    return keys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("night", help="YYYYMMDD")
    ap.add_argument("--root-raw", default="/processed1")
    ap.add_argument("--root-out", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--rr-subdir", default="rr_links", help="rr subdir under <root-out>/<night>/")
    ap.add_argument("--tracklets", default=None, help="optional explicit tracklets ALL fits")
    ap.add_argument("--out", default=None)
    ap.add_argument("--summary-out", default=None, help="optional JSON summary output path")
    args = ap.parse_args()

    night = str(args.night)
    raw_root = Path(args.root_raw)
    out_root = Path(args.root_out)

    l4_path = raw_root / night / "L4" / f"{night}_matched_asteroids.fits"
    if not l4_path.exists():
        raise FileNotFoundError(f"L4 matched asteroids file not found: {l4_path}")

    l4 = Table.read(l4_path)
    need = {"source_file", "objID"}
    if not need.issubset(l4.colnames):
        raise RuntimeError(f"L4 table missing required columns: {sorted(need - set(l4.colnames))}")

    l4_keys = build_key_array(l4, "source_file", "objID")

    tracklets_path = Path(args.tracklets) if args.tracklets else (
        out_root / night / "tracklets_linreproj" / f"tracklets_{night}_ALL.fits"
    )
    rr_dir = out_root / night / args.rr_subdir
    mask_dir = out_root / night / "mask_gaia"

    key_to_tracklets, key_to_det_rows = load_tracklet_detection_members(tracklets_path)
    tracklet_to_links, link_summary = load_rr_members(rr_dir)
    mask_keys = load_mask_gaia_keys(mask_dir)

    survives_mask = np.asarray([k in mask_keys for k in l4_keys], dtype=bool)

    tracklet_lists = []
    tracklet_endpoint_lists = []
    link_lists = []
    fit_ok_any = []
    is_good_any = []

    for key in l4_keys:
        tids = sorted(set(key_to_tracklets.get(str(key), [])))
        trk_rows = key_to_det_rows.get(str(key), [])
        endpoints = [f"{tid}:{ep}" for tid, ep in trk_rows]

        lids = sorted(set(
            lid
            for tid in tids
            for lid in tracklet_to_links.get(str(tid), [])
        ))

        tracklet_lists.append(tids)
        tracklet_endpoint_lists.append(endpoints)
        link_lists.append(lids)
        fit_ok_any.append(any(link_summary.get(lid, {}).get("fit_ok", False) for lid in lids))
        is_good_any.append(any(link_summary.get(lid, {}).get("is_good", False) for lid in lids))

    out = l4.copy(copy_data=True)
    out["det_key"] = np.asarray(l4_keys, dtype="U256")
    out["survives_mask_gaia"] = survives_mask
    out["in_tracklet"] = np.asarray([len(x) > 0 for x in tracklet_lists], dtype=bool)
    out["n_tracklets_found"] = np.asarray([len(x) for x in tracklet_lists], dtype=np.int32)
    out["tracklet_ids"] = string_list_column(tracklet_lists)
    out["tracklet_endpoints"] = string_list_column(tracklet_endpoint_lists)
    out["in_rr_link"] = np.asarray([len(x) > 0 for x in link_lists], dtype=bool)
    out["n_rr_links"] = np.asarray([len(x) for x in link_lists], dtype=np.int32)
    out["linkage_ids"] = string_list_column(link_lists, width=1024)
    out["fit_ok_any"] = np.asarray(fit_ok_any, dtype=bool)
    out["is_good_any"] = np.asarray(is_good_any, dtype=bool)

    out.meta["night"] = night
    out.meta["l4_file"] = l4_path.name
    out.meta["track_file"] = tracklets_path.name
    out.meta["rr_dir"] = rr_dir.name

    out_path = Path(args.out) if args.out else (out_root / night / "analysis" / f"{night}_known_asteroid_comparison.fits")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.write(out_path, overwrite=True)

    n_total = len(out)
    n_mask = int(np.sum(np.asarray(out["survives_mask_gaia"], dtype=bool)))
    n_track = int(np.sum(np.asarray(out["in_tracklet"], dtype=bool)))
    n_rr = int(np.sum(np.asarray(out["in_rr_link"], dtype=bool)))
    n_good = int(np.sum(np.asarray(out["is_good_any"], dtype=bool)))
    n_fit_ok = int(np.sum(np.asarray(out["fit_ok_any"], dtype=bool)))
    hit_mask = np.asarray(out["in_rr_link"], dtype=bool)
    n_rr_links_arr = np.asarray(out["n_rr_links"], dtype=np.int64)

    summary = {
        "night": night,
        "paths": {
            "l4_path": str(l4_path),
            "tracklets_path": str(tracklets_path),
            "rr_dir": str(rr_dir),
            "mask_dir": str(mask_dir),
            "comparison_fits": str(out_path),
        },
        "counts": {
            "total": n_total,
            "survives_mask_gaia": n_mask,
            "in_tracklet": n_track,
            "in_rr_link": n_rr,
            "fit_ok_any": n_fit_ok,
            "is_good_any": n_good,
        },
        "fractions": {
            "tracklet_given_mask": float(n_track / n_mask) if n_mask else 0.0,
            "rr_given_mask": float(n_rr / n_mask) if n_mask else 0.0,
            "rr_given_tracklet": float(n_rr / n_track) if n_track else 0.0,
            "fit_ok_given_rr": float(n_fit_ok / n_rr) if n_rr else 0.0,
            "is_good_given_rr": float(n_good / n_rr) if n_rr else 0.0,
        },
        "n_rr_links_all": summarize_int_array(n_rr_links_arr),
        "n_rr_links_hits_only": summarize_int_array(n_rr_links_arr[hit_mask]),
    }
    summary_path = Path(args.summary_out) if args.summary_out else (
        out_root / night / "analysis" / f"{night}_known_asteroid_comparison_summary.json"
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"[write] {out_path}")
    print(f"[write] {summary_path}")
    print(
        f"[summary] total={n_total} survives_mask_gaia={n_mask} "
        f"in_tracklet={n_track} in_rr_link={n_rr} fit_ok_any={n_fit_ok} is_good_any={n_good}"
    )


if __name__ == "__main__":
    main()
