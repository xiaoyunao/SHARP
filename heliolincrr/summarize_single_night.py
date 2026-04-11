#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.table import Table


def det_key(file_value, obj_id) -> str:
    return f"{Path(str(file_value)).name}:{int(obj_id)}"


def object_key(row) -> str:
    number = row["number"] if "number" in row.colnames else np.nan
    try:
        if number is not np.ma.masked and not math.isnan(float(number)):
            return f"num:{int(number)}"
    except Exception:
        pass
    return f"name:{str(row['name']).strip()}"


def format_float(value) -> float | None:
    try:
        if value is np.ma.masked:
            return None
        value = float(value)
    except Exception:
        return None
    return value if np.isfinite(value) else None


def semicolon_join(values) -> str:
    return ";".join(str(v) for v in values)


def summarize_count_map(counter: Counter) -> dict[str, int]:
    return {str(k): int(v) for k, v in sorted(counter.items(), key=lambda item: str(item[0]))}


def build_known_maps(matched: Table) -> tuple[dict[str, str], dict[str, set[str]]]:
    known_by_det: dict[str, str] = {}
    object_to_detkeys: dict[str, set[str]] = defaultdict(set)
    for row in matched:
        key = det_key(row["source_file"], row["objID"])
        obj = object_key(row)
        known_by_det[key] = obj
        object_to_detkeys[obj].add(key)
    return known_by_det, object_to_detkeys


def load_mask_keys(mask_dir: Path) -> set[str]:
    keys: set[str] = set()
    if not mask_dir.exists():
        return keys
    for fn in sorted(mask_dir.glob("OBJ_MP_*_cat.fits.gz")):
        try:
            t = Table.read(fn)
        except Exception:
            continue
        if "objID" not in t.colnames:
            continue
        for obj_id in np.asarray(t["objID"], dtype=np.int64):
            keys.add(f"{fn.name}:{int(obj_id)}")
    return keys


def classify_tracklets(tracklets: Table, known_by_det: dict[str, str]) -> dict[str, dict[str, object]]:
    tracklet_info: dict[str, dict[str, object]] = {}
    object_to_same2_tracklets: dict[str, list[str]] = defaultdict(list)
    before_counts = Counter()
    for row in tracklets:
        tid = str(row["tracklet_id"])
        k1 = det_key(row["file1"], row["objID1"])
        k2 = det_key(row["file2"], row["objID2"])
        o1 = known_by_det.get(k1)
        o2 = known_by_det.get(k2)
        if o1 and o2:
            cls = "2_same" if o1 == o2 else "2_diff"
        elif o1 or o2:
            cls = "1"
        else:
            cls = "0"
        info = {
            "tracklet_id": tid,
            "det_keys": [k1, k2],
            "endpoint_objects": [o1, o2],
            "class": cls,
            "group": int(row["group"]) if "group" in tracklets.colnames else -1,
            "exp_i": int(row["exp_i"]) if "exp_i" in tracklets.colnames else -1,
            "exp_j": int(row["exp_j"]) if "exp_j" in tracklets.colnames else -1,
            "v_arcsec_hr": float(row["v_arcsec_hr"]) if "v_arcsec_hr" in tracklets.colnames else None,
            "mjd1": float(row["mjd1"]),
            "mjd2": float(row["mjd2"]),
            "ra1": float(row["ra1"]),
            "ra2": float(row["ra2"]),
            "dec1": float(row["dec1"]),
            "dec2": float(row["dec2"]),
            "file1": Path(str(row["file1"])).name,
            "file2": Path(str(row["file2"])).name,
            "objID1": int(row["objID1"]),
            "objID2": int(row["objID2"]),
        }
        tracklet_info[tid] = info
        before_counts[cls] += 1
        if cls == "2_same":
            object_to_same2_tracklets[o1].append(tid)
    return {
        "tracklet_info": tracklet_info,
        "object_to_same2_tracklets": object_to_same2_tracklets,
        "before_counts": before_counts,
    }


def classify_links(
    members: Table,
    orbit_links: Table,
    tracklet_info: dict[str, dict[str, object]],
) -> tuple[dict[int, list[str]], dict[int, dict[str, object]], Counter, Counter, Counter]:
    link_to_tids: dict[int, list[str]] = defaultdict(list)
    for lid, tid in zip(np.asarray(members["linkage_id"], dtype=int), np.asarray(members["tracklet_id"]).astype("U64")):
        link_to_tids[int(lid)].append(str(tid))

    link_info: dict[int, dict[str, object]] = {}
    link_class_counts = Counter()
    fit_ok_class_counts = Counter()
    is_good_class_counts = Counter()
    orbit_by_lid = {int(row["linkage_id"]): row for row in orbit_links}

    for lid, tids in link_to_tids.items():
        objects: list[str] = []
        known_det_count = 0
        total_det_count = 0
        for tid in tids:
            for obj in tracklet_info[tid]["endpoint_objects"]:
                total_det_count += 1
                if obj:
                    known_det_count += 1
                    objects.append(obj)
        uniq = sorted(set(objects))
        if known_det_count == 0:
            link_class = "all_non_asteroid"
        elif known_det_count == total_det_count:
            link_class = "all_same_asteroid" if len(uniq) == 1 else "all_asteroid_mixed_objects"
        else:
            link_class = "mixed_with_non_asteroid"

        row = orbit_by_lid.get(lid)
        fit_ok = bool(row["fit_ok"]) if row is not None else False
        is_good = bool(row["is_good"]) if row is not None else False
        link_class_counts[link_class] += 1
        if fit_ok:
            fit_ok_class_counts[link_class] += 1
        if is_good:
            is_good_class_counts[link_class] += 1
        link_info[lid] = {
            "linkage_id": lid,
            "tracklet_ids": list(tids),
            "link_class": link_class,
            "known_det_count": known_det_count,
            "total_det_count": total_det_count,
            "objects": uniq,
            "fit_ok": fit_ok,
            "is_good": is_good,
            "orbit_row": row,
        }
    return link_to_tids, link_info, link_class_counts, fit_ok_class_counts, is_good_class_counts


def build_unknown_catalog_rows(
    link_info: dict[int, dict[str, object]],
    tracklet_info: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for lid in sorted(link_info):
        info = link_info[lid]
        if not (info["fit_ok"] and info["link_class"] == "all_non_asteroid"):
            continue
        orbit_row = info["orbit_row"]
        tids = info["tracklet_ids"]
        det_files = []
        det_objids = []
        det_mjds = []
        det_ras = []
        det_decs = []
        det_groups = []
        det_exp_pairs = []
        det_keys_seen = set()
        for tid in tids:
            tr = tracklet_info[tid]
            for endpoint in (1, 2):
                file_name = tr[f"file{endpoint}"]
                obj_id = tr[f"objID{endpoint}"]
                key = (file_name, obj_id)
                if key in det_keys_seen:
                    continue
                det_keys_seen.add(key)
                det_files.append(file_name.replace("_cat.fits.gz", ".fits.gz"))
                det_objids.append(obj_id)
                det_mjds.append(tr[f"mjd{endpoint}"])
                det_ras.append(tr[f"ra{endpoint}"])
                det_decs.append(tr[f"dec{endpoint}"])
                det_groups.append(tr["group"])
                det_exp_pairs.append(f"{tr['exp_i']}->{tr['exp_j']}")

        rows.append(
            {
                "linkage_id": lid,
                "n_tracklets": int(orbit_row["n_tracklets"]),
                "n_obs": int(orbit_row["n_obs"]),
                "fit_ok": bool(orbit_row["fit_ok"]),
                "is_good": bool(orbit_row["is_good"]),
                "rms_arcsec": format_float(orbit_row["rms_arcsec"]),
                "med_arcsec": format_float(orbit_row["med_arcsec"]),
                "max_arcsec": format_float(orbit_row["max_arcsec"]),
                "a_au": format_float(orbit_row["a_au"]),
                "ecc": format_float(orbit_row["ecc"]),
                "inc_deg": format_float(orbit_row["inc_deg"]),
                "raan_deg": format_float(orbit_row["raan_deg"]),
                "argp_deg": format_float(orbit_row["argp_deg"]),
                "nu_deg": format_float(orbit_row["nu_deg"]),
                "best_v1_kms": format_float(orbit_row["best_v1_kms"]),
                "lin_rms_arcsec": format_float(orbit_row["lin_rms_arcsec"]),
                "lin_speed_arcsec_per_day": format_float(orbit_row["lin_speed_arcsec_per_day"]),
                "lin_dir_deg": format_float(orbit_row["lin_dir_deg"]),
                "tracklet_ids": semicolon_join(tids),
                "image_names": semicolon_join(det_files),
                "objids": semicolon_join(det_objids),
                "mjds": semicolon_join(f"{x:.8f}" for x in det_mjds),
                "ras_deg": semicolon_join(f"{x:.8f}" for x in det_ras),
                "decs_deg": semicolon_join(f"{x:.8f}" for x in det_decs),
                "groups": semicolon_join(det_groups),
                "exp_pairs": semicolon_join(det_exp_pairs),
            }
        )
    return rows


def build_text_summary(summary: dict[str, object]) -> str:
    t_before = summary["tracklets_before_link"]
    t_after = summary["tracklets_after_link"]
    same2 = summary["same_asteroid_2endpoint_tracklets"]
    link_counts = summary["link_class_counts"]
    fit_counts = summary["fit_ok_link_class_counts"]
    lines = [
        f"Night: {summary['night']}",
        "",
        "Tracklets before link:",
        f"  0 endpoint asteroid: {t_before['0_endpoint_asteroid']}",
        f"  1 endpoint asteroid: {t_before['1_endpoint_asteroid']}",
        f"  2 endpoints same asteroid: {t_before['2_endpoints_same_asteroid']}",
        f"  2 endpoints different asteroids: {t_before['2_endpoints_different_asteroids']}",
        "",
        "Tracklets after link:",
        f"  0 endpoint asteroid: {t_after['0_endpoint_asteroid']}",
        f"  1 endpoint asteroid: {t_after['1_endpoint_asteroid']}",
        f"  2 endpoints same asteroid: {t_after['2_endpoints_same_asteroid']}",
        f"  2 endpoints different asteroids: {t_after['2_endpoints_different_asteroids']}",
        "",
        "Same-object 2-endpoint tracklets:",
        f"  total: {same2['total_tracklets']}",
        f"  singleton-object tracklets: {same2['singleton_object_tracklets']}",
        f"  multi-tracklet-object tracklets: {same2['multi_tracklet_object_tracklets']}",
        f"  multi-tracklet-object linked: {same2['multi_tracklet_object_tracklets_linked']}",
        f"  linked fraction excluding singletons: {same2['linked_fraction_excluding_singletons']:.6f}" if same2["linked_fraction_excluding_singletons"] is not None else "  linked fraction excluding singletons: null",
        "",
        "Link class counts:",
        f"  all_same_asteroid: {link_counts.get('all_same_asteroid', 0)}",
        f"  all_asteroid_mixed_objects: {link_counts.get('all_asteroid_mixed_objects', 0)}",
        f"  mixed_with_non_asteroid: {link_counts.get('mixed_with_non_asteroid', 0)}",
        f"  all_non_asteroid: {link_counts.get('all_non_asteroid', 0)}",
        "",
        "fit_ok link class counts:",
        f"  all_same_asteroid: {fit_counts.get('all_same_asteroid', 0)}",
        f"  all_asteroid_mixed_objects: {fit_counts.get('all_asteroid_mixed_objects', 0)}",
        f"  mixed_with_non_asteroid: {fit_counts.get('mixed_with_non_asteroid', 0)}",
        f"  all_non_asteroid: {fit_counts.get('all_non_asteroid', 0)}",
        "",
        f"Unknown fit_ok catalog rows: {summary['unknown_fit_ok_catalog']['count']}",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize the full single-night heliolincrr pipeline.")
    ap.add_argument("night", help="Target night in YYYYMMDD format")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--root-out", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--rr-subdir", default="rr_links")
    ap.add_argument("--tracklets", default=None)
    ap.add_argument("--analysis-outdir", default=None)
    ap.add_argument("--summary-json", default=None)
    ap.add_argument("--summary-txt", default=None)
    ap.add_argument("--unknown-fits", default=None)
    ap.add_argument("--unknown-json", default=None)
    args = ap.parse_args()

    night = str(args.night)
    processed_root = Path(args.processed_root)
    root_out = Path(args.root_out)
    rr_dir = root_out / night / args.rr_subdir
    tracklets_path = Path(args.tracklets) if args.tracklets else (root_out / night / "tracklets_linreproj" / f"tracklets_{night}_ALL.fits")
    analysis_outdir = Path(args.analysis_outdir) if args.analysis_outdir else (root_out / night / "analysis")
    analysis_outdir.mkdir(parents=True, exist_ok=True)

    summary_json = Path(args.summary_json) if args.summary_json else (analysis_outdir / f"{night}_single_night_summary.json")
    summary_txt = Path(args.summary_txt) if args.summary_txt else (analysis_outdir / f"{night}_single_night_summary.txt")
    unknown_fits = Path(args.unknown_fits) if args.unknown_fits else (processed_root / night / "L4" / f"{night}_unknown_links.fits")
    unknown_json = Path(args.unknown_json) if args.unknown_json else (processed_root / night / "L4" / f"{night}_unknown_links.json")

    links = Table.read(rr_dir / "links_tracklets.fits")
    members = Table.read(rr_dir / "linkage_members.fits")
    tracklets = Table.read(tracklets_path)
    orbit_links = Table.read(rr_dir / "orbit_confirm" / "orbit_links.fits")
    matched_path = processed_root / night / "L4" / f"{night}_matched_asteroids.fits"
    matched = Table.read(matched_path) if matched_path.exists() else Table()
    mask_dir = root_out / night / "mask_gaia"
    mask_keys = load_mask_keys(mask_dir)

    if len(matched) > 0:
        known_by_det, object_to_detkeys = build_known_maps(matched)
    else:
        known_by_det, object_to_detkeys = {}, defaultdict(set)

    trk = classify_tracklets(tracklets, known_by_det)
    tracklet_info = trk["tracklet_info"]
    object_to_same2_tracklets = trk["object_to_same2_tracklets"]
    linked_tids = {str(t) for t in np.asarray(members["tracklet_id"]).astype("U64")}
    after_counts = Counter(tracklet_info[tid]["class"] for tid in linked_tids)

    same2_singleton_tracklets = set()
    same2_multi_tracklets = set()
    for obj, tids in object_to_same2_tracklets.items():
        if len(tids) == 1:
            same2_singleton_tracklets.update(tids)
        else:
            same2_multi_tracklets.update(tids)

    link_to_tids, link_info, link_class_counts, fit_ok_class_counts, is_good_class_counts = classify_links(
        members, orbit_links, tracklet_info
    )

    matched_total = int(len(matched))
    survives_mask = int(sum(1 for key in known_by_det if key in mask_keys))

    failed_links = []
    for lid in sorted(link_info):
        info = link_info[lid]
        if info["fit_ok"]:
            continue
        orbit_row = info["orbit_row"]
        failed_links.append(
            {
                "linkage_id": lid,
                "link_class": info["link_class"],
                "tracklet_classes": [tracklet_info[tid]["class"] for tid in info["tracklet_ids"]],
                "tracklet_ids": info["tracklet_ids"],
                "fail_reason": str(orbit_row["fail_reason"]),
                "fail_counts": str(orbit_row["fail_counts"]),
                "min_rejected_max_v_kms": format_float(orbit_row["min_rejected_max_v_kms"]),
                "best_rejected_max_v_kms": format_float(orbit_row["best_rejected_max_v_kms"]),
            }
        )

    unknown_rows = build_unknown_catalog_rows(link_info, tracklet_info)
    unknown_fits.parent.mkdir(parents=True, exist_ok=True)
    unknown_json.parent.mkdir(parents=True, exist_ok=True)
    if unknown_rows:
        Table(rows=unknown_rows).write(unknown_fits, overwrite=True)
    else:
        if unknown_fits.exists():
            unknown_fits.unlink()
    unknown_json.write_text(json.dumps(unknown_rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    summary = {
        "night": night,
        "paths": {
            "tracklets": str(tracklets_path),
            "rr_dir": str(rr_dir),
            "matched_asteroids": str(matched_path),
            "summary_json": str(summary_json),
            "summary_txt": str(summary_txt),
            "unknown_fits": str(unknown_fits),
            "unknown_json": str(unknown_json),
        },
        "counts": {
            "matched_detections_total": matched_total,
            "mask_gaia_known_detections": survives_mask,
            "tracklets_total": int(len(tracklets)),
            "links_total": int(len(links)),
            "member_rows_total": int(len(members)),
            "orbit_fit_ok": int(np.sum(np.asarray(orbit_links["fit_ok"], dtype=bool))),
            "orbit_is_good": int(np.sum(np.asarray(orbit_links["is_good"], dtype=bool))),
        },
        "tracklets_before_link": {
            "0_endpoint_asteroid": int(trk["before_counts"]["0"]),
            "1_endpoint_asteroid": int(trk["before_counts"]["1"]),
            "2_endpoints_same_asteroid": int(trk["before_counts"]["2_same"]),
            "2_endpoints_different_asteroids": int(trk["before_counts"]["2_diff"]),
        },
        "tracklets_after_link": {
            "0_endpoint_asteroid": int(after_counts["0"]),
            "1_endpoint_asteroid": int(after_counts["1"]),
            "2_endpoints_same_asteroid": int(after_counts["2_same"]),
            "2_endpoints_different_asteroids": int(after_counts["2_diff"]),
        },
        "same_asteroid_2endpoint_tracklets": {
            "total_tracklets": int(trk["before_counts"]["2_same"]),
            "singleton_object_tracklets": int(len(same2_singleton_tracklets)),
            "singleton_objects": int(sum(1 for tids in object_to_same2_tracklets.values() if len(tids) == 1)),
            "multi_tracklet_object_tracklets": int(len(same2_multi_tracklets)),
            "multi_tracklet_objects": int(sum(1 for tids in object_to_same2_tracklets.values() if len(tids) > 1)),
            "multi_tracklet_object_tracklets_linked": int(sum(1 for tid in same2_multi_tracklets if tid in linked_tids)),
            "linked_fraction_excluding_singletons": (
                float(sum(1 for tid in same2_multi_tracklets if tid in linked_tids) / len(same2_multi_tracklets))
                if same2_multi_tracklets
                else None
            ),
        },
        "link_class_counts": summarize_count_map(link_class_counts),
        "fit_ok_link_class_counts": summarize_count_map(fit_ok_class_counts),
        "is_good_link_class_counts": summarize_count_map(is_good_class_counts),
        "failed_link_summary": {
            "count": int(len(failed_links)),
            "by_link_class": summarize_count_map(Counter(row["link_class"] for row in failed_links)),
            "by_tracklet_class_pattern": summarize_count_map(Counter(tuple(row["tracklet_classes"]) for row in failed_links)),
            "rows": failed_links,
        },
        "unknown_fit_ok_catalog": {
            "count": int(len(unknown_rows)),
            "rows": unknown_rows,
        },
    }

    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    summary_txt.write_text(build_text_summary(summary), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"[write] {summary_json}")
    print(f"[write] {summary_txt}")
    print(f"[write] {unknown_json}")
    if unknown_rows:
        print(f"[write] {unknown_fits}")


if __name__ == "__main__":
    main()
