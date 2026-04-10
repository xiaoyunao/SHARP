#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.table import Table


def detect_mode_from_rr_dir(rr_dir: Path) -> str:
    p = str(rr_dir)
    if p.endswith("/rr_links_15") or "/rr_links_15/" in p:
        return "w15"
    if p.endswith("/rr_links") or "/rr_links/" in p:
        return "night"
    return "unknown"


def extract_night_from_rr_dir(rr_dir: Path) -> str:
    if rr_dir.name in {"rr_links", "rr_links_15"}:
        return rr_dir.parent.name
    return rr_dir.name


def summarize_numeric(arr: np.ndarray) -> dict[str, object]:
    arr = np.asarray(arr, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {
            "count": 0,
            "min": None,
            "p10": None,
            "median": None,
            "mean": None,
            "p90": None,
            "p95": None,
            "max": None,
        }
    return {
        "count": int(arr.size),
        "min": float(np.min(arr)),
        "p10": float(np.percentile(arr, 10)),
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "p90": float(np.percentile(arr, 90)),
        "p95": float(np.percentile(arr, 95)),
        "max": float(np.max(arr)),
    }


def summarize_bool(mask: np.ndarray) -> dict[str, object]:
    mask = np.asarray(mask, dtype=bool)
    total = int(mask.size)
    true_count = int(np.sum(mask))
    return {
        "count": true_count,
        "total": total,
        "fraction": float(true_count / total) if total else 0.0,
    }


def bucket_counts(values: np.ndarray, bucket_specs: list[tuple[str, int | None, int | None]]) -> dict[str, object]:
    values = np.asarray(values, dtype=np.int64)
    total = int(values.size)
    out: dict[str, object] = {}
    for label, lo, hi in bucket_specs:
        mask = np.ones(values.shape, dtype=bool)
        if lo is not None:
            mask &= values >= int(lo)
        if hi is not None:
            mask &= values <= int(hi)
        count = int(np.sum(mask))
        out[label] = {
            "count": count,
            "fraction": float(count / total) if total else 0.0,
        }
    return out


def bucket_metric_by_mask(
    values: np.ndarray,
    metric_mask: np.ndarray,
    bucket_specs: list[tuple[str, int | None, int | None]],
) -> dict[str, object]:
    values = np.asarray(values, dtype=np.int64)
    metric_mask = np.asarray(metric_mask, dtype=bool)
    out: dict[str, object] = {}
    for label, lo, hi in bucket_specs:
        bucket_mask = np.ones(values.shape, dtype=bool)
        if lo is not None:
            bucket_mask &= values >= int(lo)
        if hi is not None:
            bucket_mask &= values <= int(hi)
        bucket_total = int(np.sum(bucket_mask))
        metric_count = int(np.sum(bucket_mask & metric_mask))
        out[label] = {
            "bucket_count": bucket_total,
            "metric_count": metric_count,
            "metric_fraction_within_bucket": float(metric_count / bucket_total) if bucket_total else 0.0,
        }
    return out


def parse_linkage_ids(values) -> list[list[int]]:
    out: list[list[int]] = []
    for raw in values:
        s = str(raw).strip()
        if not s:
            out.append([])
            continue
        out.append([int(x) for x in s.split(";") if x.strip()])
    return out


def summarize_known_pipeline(comp: Table) -> dict[str, object]:
    survives_mask = np.asarray(comp["survives_mask_gaia"], dtype=bool)
    in_tracklet = np.asarray(comp["in_tracklet"], dtype=bool)
    in_rr_link = np.asarray(comp["in_rr_link"], dtype=bool)
    fit_ok_any = np.asarray(comp["fit_ok_any"], dtype=bool)
    is_good_any = np.asarray(comp["is_good_any"], dtype=bool)

    total = len(comp)
    n_mask = int(np.sum(survives_mask))
    n_track = int(np.sum(in_tracklet))
    n_rr = int(np.sum(in_rr_link))
    n_fit_ok = int(np.sum(fit_ok_any))
    n_good = int(np.sum(is_good_any))

    return {
        "counts": {
            "total": int(total),
            "survives_mask_gaia": n_mask,
            "in_tracklet": n_track,
            "in_rr_link": n_rr,
            "fit_ok_any": n_fit_ok,
            "is_good_any": n_good,
        },
        "fractions": {
            "mask_given_total": float(n_mask / total) if total else 0.0,
            "tracklet_given_total": float(n_track / total) if total else 0.0,
            "rr_given_total": float(n_rr / total) if total else 0.0,
            "fit_ok_given_total": float(n_fit_ok / total) if total else 0.0,
            "is_good_given_total": float(n_good / total) if total else 0.0,
            "tracklet_given_mask": float(n_track / n_mask) if n_mask else 0.0,
            "rr_given_mask": float(n_rr / n_mask) if n_mask else 0.0,
            "rr_given_tracklet": float(n_rr / n_track) if n_track else 0.0,
            "fit_ok_given_rr": float(n_fit_ok / n_rr) if n_rr else 0.0,
            "is_good_given_rr": float(n_good / n_rr) if n_rr else 0.0,
            "is_good_given_fit_ok": float(n_good / n_fit_ok) if n_fit_ok else 0.0,
        },
        "n_rr_links_all": summarize_numeric(np.asarray(comp["n_rr_links"], dtype=np.int64)),
        "n_rr_links_hits_only": summarize_numeric(np.asarray(comp["n_rr_links"], dtype=np.int64)[in_rr_link]),
    }


def summarize_known_link_outcomes(comp: Table, orbit_links: Table) -> dict[str, object]:
    linkage_lists = parse_linkage_ids(np.asarray(comp["linkage_ids"]).astype("U2048"))
    known_link_ids = sorted({lid for lids in linkage_lists for lid in lids})
    orbit_by_lid = {
        int(row["linkage_id"]): row
        for row in orbit_links
    }
    rows = [orbit_by_lid[lid] for lid in known_link_ids if lid in orbit_by_lid]
    if not rows:
        return {
            "counts": {
                "known_hit_links": 0,
                "fit_ok_links": 0,
                "is_good_links": 0,
            },
            "fractions": {
                "fit_ok_given_known_hit_link": 0.0,
                "is_good_given_known_hit_link": 0.0,
                "is_good_given_fit_ok_link": 0.0,
            },
        }

    fit_ok = np.asarray([bool(r["fit_ok"]) for r in rows], dtype=bool)
    is_good = np.asarray([bool(r["is_good"]) for r in rows], dtype=bool)
    n_total = len(rows)
    n_fit = int(np.sum(fit_ok))
    n_good = int(np.sum(is_good))
    return {
        "counts": {
            "known_hit_links": int(n_total),
            "fit_ok_links": n_fit,
            "is_good_links": n_good,
        },
        "fractions": {
            "fit_ok_given_known_hit_link": float(n_fit / n_total) if n_total else 0.0,
            "is_good_given_known_hit_link": float(n_good / n_total) if n_total else 0.0,
            "is_good_given_fit_ok_link": float(n_good / n_fit) if n_fit else 0.0,
        },
    }


def summarize_residual_rows(resid: Table, orbit_links: Table) -> dict[str, object]:
    if len(resid) == 0:
        return {
            "all_rows": summarize_numeric(np.array([], dtype=np.float64)),
            "used_rows": summarize_numeric(np.array([], dtype=np.float64)),
            "fit_ok_used_rows": summarize_numeric(np.array([], dtype=np.float64)),
            "is_good_used_rows": summarize_numeric(np.array([], dtype=np.float64)),
        }

    fit_ok_map = {int(r["linkage_id"]): bool(r["fit_ok"]) for r in orbit_links}
    good_map = {int(r["linkage_id"]): bool(r["is_good"]) for r in orbit_links}

    lid = np.asarray(resid["linkage_id"], dtype=np.int64)
    values = np.asarray(resid["resid_arcsec"], dtype=np.float64)
    used = np.asarray(resid["used"], dtype=bool)
    fit_ok = np.asarray([fit_ok_map.get(int(x), False) for x in lid], dtype=bool)
    is_good = np.asarray([good_map.get(int(x), False) for x in lid], dtype=bool)

    return {
        "all_rows": summarize_numeric(values),
        "used_rows": summarize_numeric(values[used]),
        "fit_ok_used_rows": summarize_numeric(values[used & fit_ok]),
        "is_good_used_rows": summarize_numeric(values[used & is_good]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rr-dir", required=True)
    ap.add_argument("--comparison-fits", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rr_dir = Path(args.rr_dir).expanduser()
    orbit_dir = rr_dir / "orbit_confirm"
    orbit_links_path = orbit_dir / "orbit_links.fits"
    orbit_resid_path = orbit_dir / "orbit_obs_residuals.fits"
    if not orbit_links_path.exists():
        raise FileNotFoundError(f"orbit_links.fits not found: {orbit_links_path}")

    night = extract_night_from_rr_dir(rr_dir)
    mode = detect_mode_from_rr_dir(rr_dir)

    comparison_fits = (
        Path(args.comparison_fits).expanduser()
        if args.comparison_fits
        else rr_dir.parent / "analysis" / f"{night}_rr_known_asteroid_comparison.fits"
    )

    orbit_links = Table.read(orbit_links_path)
    orbit_resid = Table.read(orbit_resid_path) if orbit_resid_path.exists() else Table()
    comp = Table.read(comparison_fits) if comparison_fits.exists() else None

    fit_ok = np.asarray(orbit_links["fit_ok"], dtype=bool)
    is_good = np.asarray(orbit_links["is_good"], dtype=bool)
    n_tracklets = np.asarray(orbit_links["n_tracklets"], dtype=np.int64)
    n_obs = np.asarray(orbit_links["n_obs"], dtype=np.int64)
    n_nights = np.asarray(orbit_links["n_nights"], dtype=np.int64)

    bucket_specs_tracklets = [
        ("2", 2, 2),
        ("3", 3, 3),
        ("4", 4, 4),
        ("5-6", 5, 6),
        ("7-9", 7, 9),
        ("10+", 10, None),
    ]
    bucket_specs_obs = [
        ("4", 4, 4),
        ("5-6", 5, 6),
        ("7-8", 7, 8),
        ("9-12", 9, 12),
        ("13+", 13, None),
    ]
    bucket_specs_nights = [
        ("1", 1, 1),
        ("2", 2, 2),
        ("3+", 3, None),
    ]

    summary = {
        "night": night,
        "mode": mode,
        "paths": {
            "rr_dir": str(rr_dir),
            "orbit_links": str(orbit_links_path),
            "orbit_obs_residuals": str(orbit_resid_path) if orbit_resid_path.exists() else None,
            "comparison_fits": str(comparison_fits) if comparison_fits.exists() else None,
        },
        "orbit_link_counts": {
            "total_links": int(len(orbit_links)),
            "fit_ok": int(np.sum(fit_ok)),
            "is_good": int(np.sum(is_good)),
        },
        "orbit_link_fractions": {
            "fit_ok_given_all_links": float(np.mean(fit_ok)) if len(orbit_links) else 0.0,
            "is_good_given_all_links": float(np.mean(is_good)) if len(orbit_links) else 0.0,
            "is_good_given_fit_ok": float(np.sum(is_good) / np.sum(fit_ok)) if np.sum(fit_ok) else 0.0,
        },
        "bucketed_link_counts": {
            "n_tracklets_all": bucket_counts(n_tracklets, bucket_specs_tracklets),
            "n_obs_all": bucket_counts(n_obs, bucket_specs_obs),
            "n_nights_all": bucket_counts(n_nights, bucket_specs_nights),
        },
        "bucketed_orbit_success": {
            "fit_ok_by_n_tracklets": bucket_metric_by_mask(n_tracklets, fit_ok, bucket_specs_tracklets),
            "is_good_by_n_tracklets": bucket_metric_by_mask(n_tracklets, is_good, bucket_specs_tracklets),
            "fit_ok_by_n_obs": bucket_metric_by_mask(n_obs, fit_ok, bucket_specs_obs),
            "is_good_by_n_obs": bucket_metric_by_mask(n_obs, is_good, bucket_specs_obs),
            "fit_ok_by_n_nights": bucket_metric_by_mask(n_nights, fit_ok, bucket_specs_nights),
            "is_good_by_n_nights": bucket_metric_by_mask(n_nights, is_good, bucket_specs_nights),
        },
        "orbit_metric_distributions": {
            "rms_arcsec_all": summarize_numeric(np.asarray(orbit_links["rms_arcsec"], dtype=np.float64)),
            "rms_arcsec_fit_ok": summarize_numeric(np.asarray(orbit_links["rms_arcsec"], dtype=np.float64)[fit_ok]),
            "rms_arcsec_is_good": summarize_numeric(np.asarray(orbit_links["rms_arcsec"], dtype=np.float64)[is_good]),
            "max_arcsec_all": summarize_numeric(np.asarray(orbit_links["max_arcsec"], dtype=np.float64)),
            "max_arcsec_fit_ok": summarize_numeric(np.asarray(orbit_links["max_arcsec"], dtype=np.float64)[fit_ok]),
            "max_arcsec_is_good": summarize_numeric(np.asarray(orbit_links["max_arcsec"], dtype=np.float64)[is_good]),
            "med_arcsec_all": summarize_numeric(np.asarray(orbit_links["med_arcsec"], dtype=np.float64)),
            "lin_rms_arcsec_all": summarize_numeric(np.asarray(orbit_links["lin_rms_arcsec"], dtype=np.float64)),
            "lin_rms_arcsec_fit_ok": summarize_numeric(np.asarray(orbit_links["lin_rms_arcsec"], dtype=np.float64)[fit_ok]),
            "lin_rms_arcsec_is_good": summarize_numeric(np.asarray(orbit_links["lin_rms_arcsec"], dtype=np.float64)[is_good]),
            "a_au_fit_ok": summarize_numeric(np.asarray(orbit_links["a_au"], dtype=np.float64)[fit_ok]),
            "a_au_is_good": summarize_numeric(np.asarray(orbit_links["a_au"], dtype=np.float64)[is_good]),
            "ecc_fit_ok": summarize_numeric(np.asarray(orbit_links["ecc"], dtype=np.float64)[fit_ok]),
            "ecc_is_good": summarize_numeric(np.asarray(orbit_links["ecc"], dtype=np.float64)[is_good]),
        },
        "residual_row_distributions": summarize_residual_rows(orbit_resid, orbit_links),
    }

    if comp is not None:
        summary["known_asteroid_pipeline"] = summarize_known_pipeline(comp)
        summary["known_hit_link_outcomes"] = summarize_known_link_outcomes(comp, orbit_links)
    else:
        summary["known_asteroid_pipeline"] = None
        summary["known_hit_link_outcomes"] = None

    out_path = (
        Path(args.out).expanduser()
        if args.out
        else rr_dir.parent / "analysis" / f"{night}_orbit_fit_stats.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"[write] {out_path}")
    print(
        f"[summary] total_links={summary['orbit_link_counts']['total_links']} "
        f"fit_ok={summary['orbit_link_counts']['fit_ok']} "
        f"is_good={summary['orbit_link_counts']['is_good']}"
    )
    if summary["known_asteroid_pipeline"] is not None:
        kp = summary["known_asteroid_pipeline"]["counts"]
        print(
            f"[known] total={kp['total']} in_rr_link={kp['in_rr_link']} "
            f"fit_ok_any={kp['fit_ok_any']} is_good_any={kp['is_good_any']}"
        )


if __name__ == "__main__":
    main()
