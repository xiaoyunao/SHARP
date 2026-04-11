#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table


def yyyymmdd_iter(start: str, end: str) -> list[str]:
    start_dt = datetime.strptime(start, "%Y%m%d")
    end_dt = datetime.strptime(end, "%Y%m%d")
    nights: list[str] = []
    dt = start_dt
    while dt <= end_dt:
        nights.append(dt.strftime("%Y%m%d"))
        dt += timedelta(days=1)
    return nights


def det_key(file_value: object, obj_id: object) -> str:
    return f"{Path(str(file_value)).name}:{int(obj_id)}"


def object_key(row, colnames: list[str]) -> str:
    if "number" in colnames:
        number = row["number"]
        try:
            if number is not np.ma.masked and np.isfinite(float(number)):
                return f"num:{int(number)}"
        except Exception:
            pass
    return f"name:{str(row['name']).strip()}"


def build_known_maps(matched: Table) -> dict[str, str]:
    known_by_det: dict[str, str] = {}
    colnames = list(matched.colnames)
    for row in matched:
        known_by_det[det_key(row["source_file"], row["objID"])] = object_key(row, colnames)
    return known_by_det


def classify_link_classes(
    members: Table,
    tracklets: Table,
    matched: Table,
) -> dict[int, str]:
    known_by_det = build_known_maps(matched) if len(matched) else {}
    tracklet_objects: dict[str, tuple[str | None, str | None]] = {}
    for row in tracklets:
        tid = str(row["tracklet_id"])
        k1 = det_key(row["file1"], row["objID1"])
        k2 = det_key(row["file2"], row["objID2"])
        tracklet_objects[tid] = (known_by_det.get(k1), known_by_det.get(k2))

    link_to_tids: dict[int, list[str]] = defaultdict(list)
    for lid, tid in zip(np.asarray(members["linkage_id"], dtype=int), np.asarray(members["tracklet_id"]).astype("U64")):
        link_to_tids[int(lid)].append(str(tid))

    out: dict[int, str] = {}
    for lid, tids in link_to_tids.items():
        objects: list[str] = []
        known_det_count = 0
        total_det_count = 0
        for tid in tids:
            o1, o2 = tracklet_objects[tid]
            total_det_count += 2
            if o1:
                known_det_count += 1
                objects.append(o1)
            if o2:
                known_det_count += 1
                objects.append(o2)
        uniq = sorted(set(objects))
        if known_det_count == 0:
            out[lid] = "all_non_asteroid"
        elif known_det_count == total_det_count:
            out[lid] = "all_same_asteroid" if len(uniq) == 1 else "all_asteroid_mixed_objects"
        else:
            out[lid] = "mixed_with_non_asteroid"
    return out


def finite_list(values: list[object]) -> list[float]:
    out: list[float] = []
    for value in values:
        try:
            x = float(value)
        except Exception:
            continue
        if np.isfinite(x):
            out.append(x)
    return out


def percentile_or_none(values: list[float], q: float) -> float | None:
    if not values:
        return None
    return float(np.percentile(np.asarray(values, dtype=float), q))


def float_or_none(value: object) -> float | None:
    try:
        x = float(value)
    except Exception:
        return None
    return x if np.isfinite(x) else None


def plot_series(
    nights: list[str],
    key: str,
    label: str,
    rows: list[dict[str, object]],
    output_path: Path,
    title: str,
    ylabel: str,
) -> None:
    fig, ax = plt.subplots(figsize=(16, 5))
    x = np.arange(len(nights))
    y = [float_or_none(row.get(key)) for row in rows]
    y_plot = [np.nan if v is None else v for v in y]
    ax.bar(x, y_plot, color="#4C72B0", alpha=0.9, width=0.8)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Night")
    ax.set_xticks(x)
    ax.set_xticklabels(nights, rotation=90, fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.legend([label], loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_histogram(values: list[float], output_path: Path, title: str, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    if values:
        bins = min(40, max(10, int(np.sqrt(len(values)))))
        ax.hist(values, bins=bins, color="#4C72B0", alpha=0.85, edgecolor="white")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_boxplot_by_night(
    nights: list[str],
    per_night_values: list[list[float]],
    output_path: Path,
    title: str,
    ylabel: str,
) -> None:
    kept = [(night, vals) for night, vals in zip(nights, per_night_values) if vals]
    fig, ax = plt.subplots(figsize=(16, 5))
    if kept:
        ax.boxplot([vals for _, vals in kept], tick_labels=[night for night, _ in kept], showfliers=False)
        ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_scatter(
    x_values: list[float],
    y_values: list[float],
    output_path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    *,
    xlog: bool = False,
    ylog: bool = False,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    if x_values and y_values:
        ax.scatter(x_values, y_values, s=8, alpha=0.28, c="#4C72B0", edgecolors="none")
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Aggregate nightly heliolincrr metrics and plots.")
    ap.add_argument("--start", default="20251116")
    ap.add_argument("--end", default="20260116")
    ap.add_argument("--root-out", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--outdir", default="/pipeline/xiaoyunao/stats/heliolincrr")
    args = ap.parse_args()

    root_out = Path(args.root_out)
    processed_root = Path(args.processed_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    masked_plot_nights = {"20251226", "20260111", "20260119", "20260123"}

    nights = yyyymmdd_iter(args.start, args.end)
    rows: list[dict[str, object]] = []
    known_rms_all: list[float] = []
    known_med_all: list[float] = []
    known_max_all: list[float] = []
    known_a_all: list[float] = []
    known_ecc_all: list[float] = []
    known_inc_all: list[float] = []
    per_night_known_rms: list[list[float]] = []
    available_nights: list[str] = []
    missing_nights: list[str] = []

    for night in nights:
        summary_path = root_out / night / "analysis" / f"{night}_single_night_summary.json"
        if not summary_path.exists():
            missing_nights.append(night)
            continue

        summary = json.loads(summary_path.read_text())
        req = summary["requested_metrics"]
        row = {
            "night": night,
            "mask_known_delete_fraction": req["mask_gaia"]["known_detection_delete_fraction"],
            "mask_known_fraction_in_catalog": req["mask_gaia"]["known_detection_fraction_in_masked_catalog"],
            "tracklets_two_same_known": req["tracklets"]["two_endpoints_same_known"],
            "tracklets_two_diff_known": req["tracklets"]["two_endpoints_different_known"],
            "tracklets_one_known": req["tracklets"]["one_endpoint_known"],
            "tracklets_zero_known": req["tracklets"]["all_non_asteroid"],
            "links_all_same_known": req["links"]["all_same_known"],
            "links_all_diff_known": req["links"]["all_different_known"],
            "links_known_plus_non_known": req["links"]["known_plus_non_known"],
            "links_all_non_known": req["links"]["all_non_known"],
            "known_ge3_total": req["known_objects_with_ge_3_detections"]["total_objects"],
            "known_ge3_linked": req["known_objects_with_ge_3_detections"]["linked_objects"],
            "known_ge3_linked_fraction": req["known_objects_with_ge_3_detections"]["linked_fraction"],
            "links_total": summary["counts"]["links_total"],
            "orbit_fit_ok_links": req["orbit_fit"]["fit_ok_links"],
            "orbit_fit_ok_fraction": req["orbit_fit"]["fit_ok_fraction"],
            "orbit_is_good_links": req["orbit_fit"]["is_good_links"],
            "orbit_is_good_fraction": req["orbit_fit"]["is_good_fraction"],
            "orbit_all_same_known_links": req["orbit_fit"]["all_same_known_links"],
            "orbit_fit_ok_all_non_known_links": req["orbit_fit"]["fit_ok_all_non_known_links"],
            "unknown_fit_ok_links": summary["unknown_fit_ok_catalog"]["count"],
        }

        orbit_path = root_out / night / "rr_links" / "orbit_confirm" / "orbit_links.fits"
        members_path = root_out / night / "rr_links" / "linkage_members.fits"
        tracklets_path = root_out / night / "tracklets_linreproj" / f"tracklets_{night}_ALL.fits"
        matched_path = processed_root / night / "L4" / f"{night}_matched_asteroids.fits"

        known_rms: list[float] = []
        known_med: list[float] = []
        known_max: list[float] = []
        known_a: list[float] = []
        known_ecc: list[float] = []
        known_inc: list[float] = []
        if orbit_path.exists() and members_path.exists() and tracklets_path.exists() and matched_path.exists():
            orbit_links = Table.read(orbit_path)
            members = Table.read(members_path)
            tracklets = Table.read(tracklets_path)
            matched = Table.read(matched_path)
            link_classes = classify_link_classes(members, tracklets, matched)
            for orbit_row in orbit_links:
                lid = int(orbit_row["linkage_id"])
                if not bool(orbit_row["fit_ok"]):
                    continue
                if link_classes.get(lid) != "all_same_asteroid":
                    continue
                rms = float_or_none(orbit_row["rms_arcsec"])
                med = float_or_none(orbit_row["med_arcsec"])
                mx = float_or_none(orbit_row["max_arcsec"])
                if rms is not None:
                    known_rms.append(rms)
                if med is not None:
                    known_med.append(med)
                if mx is not None:
                    known_max.append(mx)
                a_au = float_or_none(orbit_row["a_au"])
                ecc = float_or_none(orbit_row["ecc"])
                inc = float_or_none(orbit_row["inc_deg"])
                if a_au is not None:
                    known_a.append(a_au)
                if ecc is not None:
                    known_ecc.append(ecc)
                if inc is not None:
                    known_inc.append(inc)

        row["known_fit_ok_all_same_known_rms_median"] = percentile_or_none(known_rms, 50)
        row["known_fit_ok_all_same_known_rms_p90"] = percentile_or_none(known_rms, 90)
        row["known_fit_ok_all_same_known_rms_max"] = max(known_rms) if known_rms else None
        row["known_fit_ok_all_same_known_med_median"] = percentile_or_none(known_med, 50)
        row["known_fit_ok_all_same_known_max_median"] = percentile_or_none(known_max, 50)
        row["known_fit_ok_all_same_known_a_au_median"] = percentile_or_none(known_a, 50)
        row["known_fit_ok_all_same_known_ecc_median"] = percentile_or_none(known_ecc, 50)
        row["known_fit_ok_all_same_known_inc_deg_median"] = percentile_or_none(known_inc, 50)
        row["known_fit_ok_all_same_known_count"] = len(known_rms)
        row["masked_in_plots"] = int(night in masked_plot_nights)

        rows.append(row)
        available_nights.append(night)
        per_night_known_rms.append(known_rms)
        known_rms_all.extend(known_rms)
        known_med_all.extend(known_med)
        known_max_all.extend(known_max)
        known_a_all.extend(known_a)
        known_ecc_all.extend(known_ecc)
        known_inc_all.extend(known_inc)

    if not rows:
        raise SystemExit("No nightly summaries found in requested range.")

    csv_path = outdir / f"nightly_metrics_{args.start}_{args.end}.csv"
    json_path = outdir / f"nightly_metrics_{args.start}_{args.end}.json"
    manifest_path = outdir / f"nightly_metrics_{args.start}_{args.end}_manifest.json"

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    json_path.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    manifest = {
        "start": args.start,
        "end": args.end,
        "available_nights": available_nights,
        "missing_nights": missing_nights,
        "masked_plot_nights": sorted(masked_plot_nights),
        "nights_with_summary": len(available_nights),
        "nights_missing_summary": len(missing_nights),
        "known_fit_ok_all_same_known_rms_samples": len(known_rms_all),
        "output_csv": str(csv_path),
        "output_json": str(json_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    for old_png in outdir.glob("*.png"):
        old_png.unlink()

    plot_rows = [row for row in rows if row["night"] not in masked_plot_nights]
    plot_nights = [row["night"] for row in plot_rows]
    plot_boxplot_rows = [row for row in rows if row["night"] not in masked_plot_nights]
    plot_boxplot_nights = [row["night"] for row in plot_boxplot_rows]
    plot_boxplot_values = [vals for night, vals in zip(available_nights, per_night_known_rms) if night not in masked_plot_nights]

    single_metric_plots = [
        ("01_mask_known_delete_fraction.png", "mask_known_delete_fraction", "Known delete fraction", "Mask Gaia known-deletion fraction by night", "Fraction"),
        ("02_mask_known_fraction_in_catalog.png", "mask_known_fraction_in_catalog", "Known fraction in masked catalog", "Known detection fraction in masked catalog by night", "Fraction"),
        ("03_tracklets_two_same_known.png", "tracklets_two_same_known", "2 endpoints same known", "Tracklets with two endpoints on same known asteroid", "Tracklet count"),
        ("04_tracklets_two_diff_known.png", "tracklets_two_diff_known", "2 endpoints different known", "Tracklets with two endpoints on different known asteroids", "Tracklet count"),
        ("05_tracklets_one_known.png", "tracklets_one_known", "1 endpoint known", "Tracklets with one known-asteroid endpoint", "Tracklet count"),
        ("06_tracklets_zero_known.png", "tracklets_zero_known", "0 endpoint known", "Tracklets with zero known-asteroid endpoints", "Tracklet count"),
        ("07_links_all_same_known.png", "links_all_same_known", "All same known", "Links whose members are all the same known asteroid", "Link count"),
        ("08_links_all_diff_known.png", "links_all_diff_known", "All different known", "Links containing only known asteroids but mixed objects", "Link count"),
        ("09_links_known_plus_non_known.png", "links_known_plus_non_known", "Known + non-known", "Links mixing known and non-known detections", "Link count"),
        ("10_links_all_non_known.png", "links_all_non_known", "All non-known", "Links made entirely of non-known detections", "Link count"),
        ("11_known_ge3_total.png", "known_ge3_total", "Known objects >=3 detections", "Known asteroids with >=3 detections in one night", "Object count"),
        ("12_known_ge3_linked.png", "known_ge3_linked", "Known objects >=3 detections linked", "Known asteroids with >=3 detections that were linked", "Object count"),
        ("13_known_ge3_linked_fraction.png", "known_ge3_linked_fraction", "Known >=3 linked fraction", "Linked fraction for known asteroids with >=3 detections", "Fraction"),
        ("14_orbit_fit_ok_links.png", "orbit_fit_ok_links", "fit_ok links", "Orbit-fit-ok link count by night", "Link count"),
        ("15_orbit_fit_ok_fraction.png", "orbit_fit_ok_fraction", "fit_ok fraction", "Orbit-fit-ok fraction by night", "Fraction"),
        ("16_orbit_is_good_links.png", "orbit_is_good_links", "is_good links", "Orbit-fit-good link count by night", "Link count"),
        ("17_orbit_is_good_fraction.png", "orbit_is_good_fraction", "is_good fraction", "Orbit-fit-good fraction by night", "Fraction"),
        ("18_orbit_all_same_known_links.png", "orbit_all_same_known_links", "All same known links", "All-same-known link count by night", "Link count"),
        ("19_orbit_fit_ok_all_non_known_links.png", "orbit_fit_ok_all_non_known_links", "fit_ok all non-known links", "fit_ok links that are all non-known by night", "Link count"),
        ("20_unknown_fit_ok_links.png", "unknown_fit_ok_links", "Unknown fit_ok links", "Unknown fit_ok link count by night", "Link count"),
        ("21_known_fit_rms_median.png", "known_fit_ok_all_same_known_rms_median", "Median RMS", "Median RMS of fit_ok all-same-known links", "Arcsec"),
        ("22_known_fit_rms_p90.png", "known_fit_ok_all_same_known_rms_p90", "P90 RMS", "P90 RMS of fit_ok all-same-known links", "Arcsec"),
        ("23_known_fit_rms_max.png", "known_fit_ok_all_same_known_rms_max", "Max RMS", "Max RMS of fit_ok all-same-known links", "Arcsec"),
    ]

    for filename, key, label, title, ylabel in single_metric_plots:
        plot_series(plot_nights, key, label, plot_rows, outdir / filename, title, ylabel)

    plot_histogram(
        known_rms_all,
        outdir / "27_known_fit_rms_histogram.png",
        "Known All-Same-Asteroid fit_ok RMS Distribution",
        "RMS (arcsec)",
    )
    plot_boxplot_by_night(
        plot_boxplot_nights,
        plot_boxplot_values,
        outdir / "28_known_fit_rms_boxplot_by_night.png",
        "Known All-Same-Asteroid fit_ok RMS by Night",
        "RMS (arcsec)",
    )
    plot_histogram(
        known_a_all,
        outdir / "29_known_fit_a_au_histogram.png",
        "Known All-Same-Asteroid fit_ok semi-major axis distribution",
        "a (AU)",
    )
    plot_histogram(
        known_ecc_all,
        outdir / "30_known_fit_ecc_histogram.png",
        "Known All-Same-Asteroid fit_ok eccentricity distribution",
        "ecc",
    )
    plot_histogram(
        known_inc_all,
        outdir / "31_known_fit_inc_deg_histogram.png",
        "Known All-Same-Asteroid fit_ok inclination distribution",
        "inc (deg)",
    )

    a_ecc_pairs = [
        (a, ecc)
        for a, ecc in zip(known_a_all, known_ecc_all)
        if 0.0 < a < 1000.0 and 0.0 <= ecc < 1.5
    ]
    a_inc_pairs = [
        (a, inc)
        for a, inc in zip(known_a_all, known_inc_all)
        if 0.0 < a < 1000.0 and 0.0 <= inc < 180.0
    ]
    ecc_inc_pairs = [
        (ecc, inc)
        for ecc, inc in zip(known_ecc_all, known_inc_all)
        if 0.0 <= ecc < 1.5 and 0.0 <= inc < 180.0
    ]

    plot_scatter(
        [x for x, _ in a_ecc_pairs],
        [y for _, y in a_ecc_pairs],
        outdir / "24_known_fit_a_vs_ecc_scatter.png",
        "Known fit_ok all-same-known links: semimajor axis vs eccentricity",
        "Semimajor axis a [AU]",
        "Eccentricity",
        xlog=True,
    )
    plot_scatter(
        [x for x, _ in a_inc_pairs],
        [y for _, y in a_inc_pairs],
        outdir / "25_known_fit_a_vs_inc_scatter.png",
        "Known fit_ok all-same-known links: semimajor axis vs inclination",
        "Semimajor axis a [AU]",
        "Inclination [deg]",
        xlog=True,
        ylog=True,
    )
    plot_scatter(
        [x for x, _ in ecc_inc_pairs],
        [y for _, y in ecc_inc_pairs],
        outdir / "26_known_fit_ecc_vs_inc_scatter.png",
        "Known fit_ok all-same-known links: eccentricity vs inclination",
        "Eccentricity",
        "Inclination [deg]",
        ylog=True,
    )

    print(f"[write] {csv_path}")
    print(f"[write] {json_path}")
    print(f"[write] {manifest_path}")
    for name in sorted(p.name for p in outdir.glob("*.png")):
        print(f"[write] {outdir / name}")


if __name__ == "__main__":
    main()
