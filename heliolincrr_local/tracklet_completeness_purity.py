#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
from astropy.table import Table


def det_key(file_value, obj_id) -> str:
    return f"{Path(str(file_value)).name}:{int(obj_id)}"


def asteroid_key(row) -> str:
    number = row["number"]
    if np.isfinite(number) and int(number) > 0:
        return f"num:{int(number)}"
    return f"name:{str(row['name']).strip()}"


def asteroid_label(row) -> str:
    number = row["number"]
    if np.isfinite(number) and int(number) > 0:
        return f"{int(number)} {str(row['name']).strip()}"
    return str(row["name"]).strip()


def build_known_detection_maps(matched: Table):
    det_to_obj = {}
    obj_det_counts = Counter()
    obj_labels = {}

    for row in matched:
        key = det_key(row["source_file"], row["objID"])
        obj = asteroid_key(row)
        det_to_obj[key] = obj
        obj_det_counts[obj] += 1
        obj_labels[obj] = asteroid_label(row)

    return det_to_obj, obj_det_counts, obj_labels


def summarize_tracklets(tracklets: Table, det_to_obj: dict[str, str]):
    summary = {
        "tracklets_total": int(len(tracklets)),
        "tracklets_both_known": 0,
        "tracklets_one_known": 0,
        "tracklets_zero_known": 0,
        "tracklets_both_known_same_object": 0,
        "tracklets_both_known_diff_object": 0,
    }
    recovered_objects = set()
    tracklets_by_object = Counter()

    for row in tracklets:
        k1 = det_key(row["file1"], row["objID1"])
        k2 = det_key(row["file2"], row["objID2"])
        o1 = det_to_obj.get(k1)
        o2 = det_to_obj.get(k2)

        n_known = int(o1 is not None) + int(o2 is not None)
        if n_known == 2:
            summary["tracklets_both_known"] += 1
            if o1 == o2:
                summary["tracklets_both_known_same_object"] += 1
                recovered_objects.add(o1)
                tracklets_by_object[o1] += 1
            else:
                summary["tracklets_both_known_diff_object"] += 1
        elif n_known == 1:
            summary["tracklets_one_known"] += 1
        else:
            summary["tracklets_zero_known"] += 1

    return summary, recovered_objects, tracklets_by_object


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("night", help="YYYYMMDD")
    project_root = Path(__file__).resolve().parent.parent
    ap.add_argument("--root-raw", default=str(project_root / "local_data" / "processed1"))
    ap.add_argument("--root-out", default=str(project_root / "local_data" / "heliolincrr"))
    ap.add_argument("--tracklets", default=None, help="optional explicit nightly tracklets FITS")
    ap.add_argument("--out", default=None, help="optional JSON summary output path")
    args = ap.parse_args()

    night = str(args.night)
    raw_root = Path(args.root_raw)
    out_root = Path(args.root_out)

    matched_path = raw_root / night / "L4" / f"{night}_matched_asteroids.fits"
    tracklets_path = Path(args.tracklets) if args.tracklets else (
        out_root / night / "tracklets_linreproj" / f"tracklets_{night}_ALL.fits"
    )

    matched = Table.read(matched_path)
    tracklets = Table.read(tracklets_path)

    need_matched = {"name", "number", "source_file", "objID"}
    need_tracklets = {"file1", "file2", "objID1", "objID2", "tracklet_id"}
    missing_matched = need_matched - set(matched.colnames)
    missing_tracklets = need_tracklets - set(tracklets.colnames)
    if missing_matched:
        raise RuntimeError(f"matched_asteroids missing columns: {sorted(missing_matched)}")
    if missing_tracklets:
        raise RuntimeError(f"tracklets missing columns: {sorted(missing_tracklets)}")

    det_to_obj, obj_det_counts, obj_labels = build_known_detection_maps(matched)
    tracklet_summary, recovered_objects, tracklets_by_object = summarize_tracklets(tracklets, det_to_obj)

    eligible_objects = sorted([obj for obj, n in obj_det_counts.items() if n >= 2])
    recovered_eligible = sorted([obj for obj in eligible_objects if obj in recovered_objects])

    eligible_detection_total = int(sum(obj_det_counts[obj] for obj in eligible_objects))
    recovered_detection_total = int(sum(obj_det_counts[obj] for obj in recovered_eligible))

    summary = {
        "night": night,
        "matched_detections_total": int(len(matched)),
        "known_objects_total": int(len(obj_det_counts)),
        "eligible_known_objects_gte2": int(len(eligible_objects)),
        "recovered_known_objects_gte2": int(len(recovered_eligible)),
        "completeness_object_fraction": (
            float(len(recovered_eligible) / len(eligible_objects)) if eligible_objects else 0.0
        ),
        "eligible_known_detections_total": eligible_detection_total,
        "recovered_known_detections_total": recovered_detection_total,
        "completeness_detection_fraction": (
            float(recovered_detection_total / eligible_detection_total) if eligible_detection_total else 0.0
        ),
        **tracklet_summary,
        "purity_both_known_fraction": (
            float(tracklet_summary["tracklets_both_known"] / tracklet_summary["tracklets_total"])
            if tracklet_summary["tracklets_total"] else 0.0
        ),
        "purity_both_known_same_object_fraction": (
            float(tracklet_summary["tracklets_both_known_same_object"] / tracklet_summary["tracklets_total"])
            if tracklet_summary["tracklets_total"] else 0.0
        ),
        "purity_one_known_fraction": (
            float(tracklet_summary["tracklets_one_known"] / tracklet_summary["tracklets_total"])
            if tracklet_summary["tracklets_total"] else 0.0
        ),
        "purity_zero_known_fraction": (
            float(tracklet_summary["tracklets_zero_known"] / tracklet_summary["tracklets_total"])
            if tracklet_summary["tracklets_total"] else 0.0
        ),
        "top_recovered_objects": [
            {
                "object_key": obj,
                "label": obj_labels.get(obj, obj),
                "matched_detections": int(obj_det_counts[obj]),
                "same_object_tracklets": int(tracklets_by_object[obj]),
            }
            for obj in sorted(tracklets_by_object, key=lambda x: (-tracklets_by_object[x], -obj_det_counts[x], x))[:20]
        ],
    }

    print(json.dumps(summary, indent=2, sort_keys=True))

    if args.out:
        out_path = Path(args.out)
    else:
        out_path = out_root / night / "analysis" / f"{night}_tracklet_completeness_purity.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[write] {out_path}")


if __name__ == "__main__":
    main()
