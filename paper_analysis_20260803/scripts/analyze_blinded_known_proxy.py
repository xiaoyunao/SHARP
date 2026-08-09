#!/usr/bin/env python3
"""Measure identity-blind unknown-link survival using known detections.

The unknown linker does not use known-object identity until its post-orbit
classification/subtraction step.  Therefore, matching frozen 1.5-arcsec known
detections to frozen orbit-link observation keys provides a retrospective
link/orbit survival proxy conditional on L2 detection.  It is not an image-
level injection test and does not measure L2 detection completeness.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


V_BINS = np.array([8, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 22.5], dtype=float)
MEASURED_MAG_BINS = np.array([10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], dtype=float)
RATE_BINS = np.array([3, 5, 10, 20, 30, 40, 50, 63], dtype=float)
DIRECTION_BINS = np.arange(0, 361, 45, dtype=float)
EDGE_BINS = np.array([0, 100, 300, 500, 1000, 2000, 4000, np.inf], dtype=float)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def object_key(frame: pd.DataFrame) -> pd.Series:
    number = frame["asteroid_number"].astype("string").str.strip()
    name = frame["asteroid_name"].astype("string").str.strip()
    valid_number = number.notna() & ~number.isin(["", "None", "nan", "--"])
    fallback = "NAME:" + name.fillna("UNKNOWN")
    return fallback.mask(valid_number, "NUMBER:" + number)


def normalize_obj_id(value) -> str:
    try:
        numeric = float(value)
        if np.isfinite(numeric) and numeric.is_integer():
            return str(int(numeric))
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def finite_median(values: pd.Series) -> float:
    array = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    array = array[np.isfinite(array)]
    return float(np.median(array)) if len(array) else np.nan


def wilson(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return np.nan, np.nan
    fraction = success / total
    denominator = 1 + z * z / total
    center = (fraction + z * z / (2 * total)) / denominator
    half = z * math.sqrt(fraction * (1 - fraction) / total + z * z / (4 * total * total)) / denominator
    return center - half, center + half


def binned_rows(frame: pd.DataFrame, column: str, edges: np.ndarray, label: str) -> list[dict]:
    values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
    rows = []
    for lower, upper in zip(edges[:-1], edges[1:]):
        selection = np.isfinite(values) & (values >= lower) & (values < upper)
        sample = frame.loc[selection]
        denominator = int(len(sample))
        for stage in ["strict_linked", "fit_ok_survived", "is_good_survived"]:
            success = int(sample[stage].fillna(False).sum())
            low, high = wilson(success, denominator)
            rows.append(
                {
                    "metric": label,
                    "lower": lower,
                    "upper": upper,
                    "stage": stage,
                    "eligible_object_n": denominator,
                    "survived_object_n": success,
                    "survival_fraction": success / denominator if denominator else np.nan,
                    "ci95_low": low,
                    "ci95_high": high,
                    "denominator_definition": "object-night with >=3 distinct 1.5-arcsec known detections and 3<=first-last rate<=63 arcsec/h",
                }
            )
    return rows


def load_conditions(path: Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    frame = pd.read_csv(path, dtype={"night": "string"})
    frame["night"] = frame["night"].str.zfill(8)
    frame["source_file"] = frame["file_name"].astype(str)
    keep = [
        "night",
        "source_file",
        "n_sources",
        "seeing_arcsec",
        "sky_adu",
        "sky_rms_adu",
        "sky_mag",
        "limiting_mag",
    ]
    keep = [column for column in keep if column in frame]
    frame = frame[keep]
    for column in keep:
        if column not in {"night", "source_file"}:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame.drop_duplicates(["night", "source_file"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--products-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--l2-manifest", type=Path)
    args = parser.parse_args()
    products = args.products_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)

    mask_path = products / "known_mask15.parquet"
    obs_path = products / "orbit_obs_residuals.parquet"
    link_path = products / "orbit_links.parquet"
    detections = pd.read_parquet(mask_path)
    detections["night"] = detections["night"].astype("string").str.zfill(8)
    detections["object_key"] = object_key(detections)
    detections["obs_key"] = (
        detections["source_file"].astype(str) + ":" + detections["obj_id"].map(normalize_obj_id)
    )
    source_detection_rows = int(len(detections))
    detections = detections.drop_duplicates(["night", "object_key", "obs_key"]).copy()
    conditions = load_conditions(args.l2_manifest.resolve() if args.l2_manifest else None)
    if not conditions.empty:
        detections = detections.merge(
            conditions, on=["night", "source_file"], how="left", validate="many_to_one"
        )

    detections["detector_edge_distance_px"] = np.minimum.reduce(
        [
            pd.to_numeric(detections["x_win_px"], errors="coerce"),
            pd.to_numeric(detections["y_win_px"], errors="coerce"),
            9216.0 - pd.to_numeric(detections["x_win_px"], errors="coerce"),
            9232.0 - pd.to_numeric(detections["y_win_px"], errors="coerce"),
        ]
    )
    detections = detections.sort_values(["night", "object_key", "epoch_mjd", "obs_key"])
    grouped = detections.groupby(["night", "object_key"], sort=False)
    first = grouped.first().reset_index()
    last = grouped.last().reset_index()
    object_frame = grouped.agg(
        detected_1p5arcsec_n=("obs_key", "nunique"),
        first_mjd=("epoch_mjd", "min"),
        last_mjd=("epoch_mjd", "max"),
        median_pred_v_mag=("pred_v_mag", finite_median),
        median_measured_mag=("mag_aper4", finite_median),
        median_detector_edge_px=("detector_edge_distance_px", finite_median),
    ).reset_index()
    for column in ["n_sources", "seeing_arcsec", "sky_adu", "sky_rms_adu", "sky_mag", "limiting_mag"]:
        if column in detections:
            extra = grouped[column].apply(finite_median).rename(f"median_{column}").reset_index()
            object_frame = object_frame.merge(extra, on=["night", "object_key"], how="left", validate="one_to_one")

    endpoints = first[["night", "object_key", "pred_ra_deg", "pred_dec_deg"]].merge(
        last[["night", "object_key", "pred_ra_deg", "pred_dec_deg"]],
        on=["night", "object_key"],
        suffixes=("_first", "_last"),
        validate="one_to_one",
    )
    object_frame = object_frame.merge(endpoints, on=["night", "object_key"], how="left", validate="one_to_one")
    dra = (
        (object_frame["pred_ra_deg_last"] - object_frame["pred_ra_deg_first"] + 180.0) % 360.0
    ) - 180.0
    mean_dec = (object_frame["pred_dec_deg_first"] + object_frame["pred_dec_deg_last"]) / 2.0
    dra_cosdec = dra * np.cos(np.deg2rad(mean_dec))
    ddec = object_frame["pred_dec_deg_last"] - object_frame["pred_dec_deg_first"]
    span_hours = (object_frame["last_mjd"] - object_frame["first_mjd"]) * 24.0
    object_frame["first_last_rate_arcsec_per_hour"] = np.divide(
        np.hypot(dra_cosdec, ddec) * 3600.0,
        span_hours,
        out=np.full(len(object_frame), np.nan),
        where=span_hours > 0,
    )
    object_frame["first_last_direction_deg"] = np.mod(np.rad2deg(np.arctan2(ddec, dra_cosdec)), 360.0)
    object_frame["eligible_ge3_detections"] = object_frame["detected_1p5arcsec_n"] >= 3
    object_frame["eligible_ge3_and_rate"] = (
        object_frame["eligible_ge3_detections"]
        & object_frame["first_last_rate_arcsec_per_hour"].between(3.0, 63.0, inclusive="both")
    )

    observation_map = detections[["night", "obs_key", "object_key"]].drop_duplicates()
    association_n = observation_map.groupby(["night", "obs_key"])["object_key"].transform("nunique")
    ambiguous_obs_key_n = int(observation_map.loc[association_n > 1, ["night", "obs_key"]].drop_duplicates().shape[0])
    observation_map = observation_map.loc[association_n == 1]
    orbit_obs = pd.read_parquet(obs_path, columns=["night", "linkage_id", "obs_key", "used"])
    orbit_obs["night"] = orbit_obs["night"].astype("string").str.zfill(8)
    mapped = orbit_obs.merge(observation_map, on=["night", "obs_key"], how="left", validate="many_to_one")
    classified = mapped.groupby(["night", "linkage_id"], as_index=False).agg(
        link_obs_n=("obs_key", "nunique"),
        mapped_obs_n=("object_key", "count"),
        mapped_object_n=("object_key", "nunique"),
        mapped_object_key=("object_key", "first"),
        used_obs_n=("used", lambda values: int(pd.Series(values).fillna(False).astype(bool).sum())),
    )
    classified["strict_single_known_object"] = (
        (classified["mapped_object_n"] == 1)
        & (classified["mapped_obs_n"] == classified["link_obs_n"])
        & (classified["link_obs_n"] >= 3)
    )
    orbit_links = pd.read_parquet(link_path, columns=["night", "linkage_id", "fit_ok", "is_good"])
    orbit_links["night"] = orbit_links["night"].astype("string").str.zfill(8)
    classified = classified.merge(
        orbit_links, on=["night", "linkage_id"], how="left", validate="one_to_one"
    )
    classified.to_csv(output / "blinded_known_link_classification.csv", index=False)

    strict = classified.loc[classified["strict_single_known_object"]].copy()
    object_survival = strict.groupby(["night", "mapped_object_key"], as_index=False).agg(
        strict_link_n=("linkage_id", "nunique"),
        strict_linked=("linkage_id", lambda values: len(values) > 0),
        fit_ok_survived=("fit_ok", lambda values: bool(pd.Series(values).fillna(False).astype(bool).any())),
        is_good_survived=("is_good", lambda values: bool(pd.Series(values).fillna(False).astype(bool).any())),
    ).rename(columns={"mapped_object_key": "object_key"})
    object_frame = object_frame.merge(
        object_survival, on=["night", "object_key"], how="left", validate="one_to_one"
    )
    for column in ["strict_linked", "fit_ok_survived", "is_good_survived"]:
        object_frame[column] = object_frame[column].astype("boolean").fillna(False).astype(bool)
    object_frame["strict_link_n"] = object_frame["strict_link_n"].fillna(0).astype(int)
    object_frame.to_csv(output / "blinded_known_recovery_by_object.csv", index=False)

    eligible = object_frame.loc[object_frame["eligible_ge3_and_rate"]].copy()
    binned: list[dict] = []
    for column, edges, label in [
        ("median_pred_v_mag", V_BINS, "predicted_v_mag"),
        ("median_measured_mag", MEASURED_MAG_BINS, "measured_aperture_magnitude"),
        ("first_last_rate_arcsec_per_hour", RATE_BINS, "first_last_rate_arcsec_per_hour"),
        ("first_last_direction_deg", DIRECTION_BINS, "first_last_direction_deg_east_through_north"),
        ("median_detector_edge_px", EDGE_BINS, "median_detector_edge_distance_px"),
    ]:
        binned.extend(binned_rows(eligible, column, edges, label))
    for condition_column in ["median_n_sources", "median_seeing_arcsec", "median_sky_adu", "median_limiting_mag"]:
        if condition_column not in eligible:
            continue
        values = pd.to_numeric(eligible[condition_column], errors="coerce")
        finite_values = values[np.isfinite(values)]
        if finite_values.nunique() < 4:
            continue
        edges = np.unique(np.nanpercentile(finite_values, [0, 20, 40, 60, 80, 100]))
        if len(edges) >= 3:
            edges[-1] = np.nextafter(edges[-1], np.inf)
            binned.extend(binned_rows(eligible, condition_column, edges, condition_column))
    pd.DataFrame(binned).to_csv(output / "blinded_known_recovery_binned.csv", index=False)

    def count_survival(frame: pd.DataFrame) -> dict[str, int | float | None]:
        denominator = int(len(frame))
        result: dict[str, int | float | None] = {"eligible_object_n": denominator}
        for stage in ["strict_linked", "fit_ok_survived", "is_good_survived"]:
            success = int(frame[stage].sum())
            result[f"{stage}_n"] = success
            result[f"{stage}_fraction"] = success / denominator if denominator else None
        return result

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": "retrospective identity-blind linkage/orbit survival proxy conditional on known 1.5-arcsec detections",
        "guardrail": "not L2 detection completeness; not image injection; intermediate Gaia/static/tracklet losses are not separately identifiable from the frozen final link products",
        "source_known_mask15_rows": source_detection_rows,
        "deduplicated_known_detection_rows": int(len(detections)),
        "object_night_rows": int(len(object_frame)),
        "eligible_ge3_detections_n": int(object_frame["eligible_ge3_detections"].sum()),
        "eligible_ge3_and_rate_n": int(object_frame["eligible_ge3_and_rate"].sum()),
        "ambiguous_detection_keys_excluded_from_link_mapping_n": ambiguous_obs_key_n,
        "strict_single_known_object_links_n": int(classified["strict_single_known_object"].sum()),
        "survival_for_ge3_and_rate": count_survival(eligible),
        "survival_for_ge3_all_rates": count_survival(object_frame.loc[object_frame["eligible_ge3_detections"]]),
        "input_hashes": {
            "known_mask15.parquet": sha256(mask_path),
            "orbit_obs_residuals.parquet": sha256(obs_path),
            "orbit_links.parquet": sha256(link_path),
        },
    }
    if args.l2_manifest:
        summary["input_hashes"]["l2_manifest.csv"] = sha256(args.l2_manifest.resolve())
    (output / "blinded_known_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
