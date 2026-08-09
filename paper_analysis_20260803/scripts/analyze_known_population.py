#!/usr/bin/env python3
"""Build nominal known-object recovery and astrometric-residual products.

The prediction denominator is the production ``all_asteroids`` table: V < 22.5
inside the nominal WCS-projected detector rectangle.  It is deliberately called
``nominal`` rather than ``detectable`` because bad pixels, saturation, missing
catalog regions, and a formal limiting-magnitude selection are not encoded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord, get_body, get_sun
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers


iers.conf.auto_download = False


V_BINS = np.array([8, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 22.5], dtype=float)
RATE_BINS = np.array([0, 3, 5, 10, 20, 40, 80, 160, 320, 640, np.inf], dtype=float)
MEASURED_MAG_BINS = np.array([10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], dtype=float)
EDGE_BINS = np.array([0, 100, 300, 500, 1000, 2000, 4000, np.inf], dtype=float)
SEEING_BINS = np.array([0, 1, 1.5, 2, 2.5, 3, 4, 6, np.inf], dtype=float)
SOURCE_COUNT_BINS = np.array([0, 5000, 10000, 20000, 40000, 80000, 160000, np.inf], dtype=float)
SKY_ADU_BINS = np.array([0, 100, 300, 1000, 3000, 10000, 30000, np.inf], dtype=float)
LIMITING_MAG_BINS = np.array([12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25], dtype=float)
MOON_SEPARATION_BINS = np.array([0, 30, 60, 90, 120, 150, 180.001], dtype=float)
MOON_ILLUMINATION_BINS = np.array([0, 0.25, 0.5, 0.75, 1.001], dtype=float)


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
    valid_name = name.notna() & ~name.isin(["", "None", "nan", "--"])
    fallback = pd.Series("UNKNOWN", index=frame.index, dtype="string")
    fallback.loc[valid_name] = "NAME:" + name.loc[valid_name]
    fallback.loc[valid_number] = "NUMBER:" + number.loc[valid_number]
    return fallback


def key_series(frame: pd.DataFrame) -> pd.Series:
    return frame["source_file"].astype("string") + "||" + frame["object_key"].astype("string")


def angular_speed(frame: pd.DataFrame) -> np.ndarray:
    """Central/one-sided finite-difference speed in arcsec per hour."""
    if frame.empty:
        return np.empty(0)
    work = frame[["object_key", "epoch_mjd", "pred_ra_deg", "pred_dec_deg"]].copy()
    work["_original"] = np.arange(len(work))
    work = work.sort_values(["object_key", "epoch_mjd", "_original"])
    grouped = work.groupby("object_key", sort=False)
    previous = grouped[["epoch_mjd", "pred_ra_deg", "pred_dec_deg"]].shift(1)
    following = grouped[["epoch_mjd", "pred_ra_deg", "pred_dec_deg"]].shift(-1)
    has_previous = previous["epoch_mjd"].notna()
    has_following = following["epoch_mjd"].notna()
    t1 = np.where(has_previous, previous["epoch_mjd"], work["epoch_mjd"])
    ra1 = np.where(has_previous, previous["pred_ra_deg"], work["pred_ra_deg"])
    dec1 = np.where(has_previous, previous["pred_dec_deg"], work["pred_dec_deg"])
    t2 = np.where(has_following, following["epoch_mjd"], work["epoch_mjd"])
    ra2 = np.where(has_following, following["pred_ra_deg"], work["pred_ra_deg"])
    dec2 = np.where(has_following, following["pred_dec_deg"], work["pred_dec_deg"])
    dra = ((ra2 - ra1 + 180.0) % 360.0) - 180.0
    dra *= np.cos(np.deg2rad((dec1 + dec2) / 2.0))
    separation_arcsec = np.hypot(dra, dec2 - dec1) * 3600.0
    delta_hours = (t2 - t1) * 24.0
    speed = np.divide(
        separation_arcsec,
        delta_hours,
        out=np.full_like(separation_arcsec, np.nan, dtype=float),
        where=delta_hours > 0,
    )
    result = np.full(len(work), np.nan)
    result[work["_original"].to_numpy()] = speed
    return result


def wilson_interval(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return np.nan, np.nan
    p = success / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return center - half, center + half


def binned_recovery(frame: pd.DataFrame, variable: str, bins: np.ndarray, label: str) -> list[dict]:
    values = pd.to_numeric(frame[variable], errors="coerce").to_numpy(dtype=float)
    matched = frame["matched_1arcsec"].to_numpy(dtype=bool)
    rows = []
    for lower, upper in zip(bins[:-1], bins[1:]):
        selection = np.isfinite(values) & (values >= lower) & (values < upper)
        total = int(selection.sum())
        success = int(matched[selection].sum())
        low, high = wilson_interval(success, total)
        rows.append(
            {
                "metric": label,
                "lower": lower,
                "upper": upper,
                "denominator_prediction_n": total,
                "matched_1arcsec_n": success,
                "recovery_fraction": success / total if total else np.nan,
                "ci95_low": low,
                "ci95_high": high,
                "denominator_definition": "V<22.5 and inside nominal WCS detector rectangle",
            }
        )
    return rows


def residual_bins(frame: pd.DataFrame, variable: str, bins: np.ndarray, label: str) -> list[dict]:
    values = pd.to_numeric(frame[variable], errors="coerce").to_numpy(dtype=float)
    rows = []
    for lower, upper in zip(bins[:-1], bins[1:]):
        selection = np.isfinite(values) & (values >= lower) & (values < upper)
        sample = frame.loc[selection]
        radial = sample["radial_residual_arcsec"].to_numpy(dtype=float)
        dra = sample["dra_cosdec_arcsec"].to_numpy(dtype=float)
        ddec = sample["ddec_arcsec"].to_numpy(dtype=float)
        rows.append(
            {
                "metric": label,
                "lower": lower,
                "upper": upper,
                "matched_n": int(len(sample)),
                "median_radial_arcsec": float(np.nanmedian(radial)) if len(sample) else np.nan,
                "p84_radial_arcsec": float(np.nanpercentile(radial, 84)) if len(sample) else np.nan,
                "median_dra_cosdec_arcsec": float(np.nanmedian(dra)) if len(sample) else np.nan,
                "median_ddec_arcsec": float(np.nanmedian(ddec)) if len(sample) else np.nan,
                "robust_sigma_dra_arcsec": float(0.7413 * (np.nanpercentile(dra, 75) - np.nanpercentile(dra, 25))) if len(sample) else np.nan,
                "robust_sigma_ddec_arcsec": float(0.7413 * (np.nanpercentile(ddec, 75) - np.nanpercentile(ddec, 25))) if len(sample) else np.nan,
            }
        )
    return rows


def conditions_table(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path, dtype={"night": "string"})
    frame["night"] = frame["night"].str.zfill(8)
    keep = [
        "night",
        "file_name",
        "mjd",
        "exposure_s",
        "seeing_arcsec",
        "sky_adu",
        "sky_rms_adu",
        "sky_mag",
        "limiting_mag",
        "n_sources",
        "center_ra_deg",
        "center_dec_deg",
    ]
    keep = [column for column in keep if column in frame]
    frame = frame[keep].rename(columns={"file_name": "source_file"})
    for column in frame.columns:
        if column not in {"night", "source_file"}:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    required = {"mjd", "center_ra_deg", "center_dec_deg"}
    if required.issubset(frame.columns):
        valid = frame[list(required)].notna().all(axis=1)
        frame["moon_separation_deg"] = np.nan
        frame["moon_illumination_fraction"] = np.nan
        if valid.any():
            midpoint = frame.loc[valid, "mjd"].to_numpy(dtype=float)
            if "exposure_s" in frame:
                midpoint += frame.loc[valid, "exposure_s"].fillna(30.0).to_numpy(dtype=float) / (2 * 86400.0)
            times = Time(midpoint, format="mjd", scale="utc")
            site = EarthLocation.from_geodetic(117.575 * u.deg, 40.393 * u.deg, 960.0 * u.m)
            pointing = SkyCoord(
                ra=frame.loc[valid, "center_ra_deg"].to_numpy(dtype=float) * u.deg,
                dec=frame.loc[valid, "center_dec_deg"].to_numpy(dtype=float) * u.deg,
                frame="icrs",
            )
            moon = get_body("moon", times, location=site)
            apparent_pointing = pointing.transform_to(moon.frame)
            geocentric_moon = get_body("moon", times)
            sun = get_sun(times)
            frame.loc[valid, "moon_separation_deg"] = apparent_pointing.separation(moon).deg
            frame.loc[valid, "moon_illumination_fraction"] = (
                1.0 - np.cos(sun.separation(geocentric_moon).rad)
            ) / 2.0
    return frame.drop_duplicates(["night", "source_file"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--products-dir", type=Path, required=True)
    parser.add_argument("--file-status", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--l2-manifest", type=Path)
    parser.add_argument("--random-shift", type=Path)
    parser.add_argument("--write-detection-csv", action="store_true")
    args = parser.parse_args()
    products = args.products_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    matched = pd.read_parquet(products / "known_matched.parquet")
    matched["night"] = matched["night"].astype(str).str.zfill(8)
    matched["object_key"] = object_key(matched)
    matched["prediction_key"] = key_series(matched)
    duplicate_matched = int(matched.duplicated(["night", "prediction_key"]).sum())
    if duplicate_matched:
        raise ValueError(f"matched prediction key is not unique: {duplicate_matched} duplicate rows")

    dra_deg = ((matched["obs_ra_win_deg"] - matched["pred_ra_deg"] + 180.0) % 360.0) - 180.0
    matched["dra_cosdec_arcsec"] = dra_deg * np.cos(np.deg2rad(matched["pred_dec_deg"])) * 3600.0
    matched["ddec_arcsec"] = (matched["obs_dec_win_deg"] - matched["pred_dec_deg"]) * 3600.0
    matched["radial_residual_arcsec"] = np.hypot(
        matched["dra_cosdec_arcsec"], matched["ddec_arcsec"]
    )
    matched["detector_edge_distance_px"] = np.minimum.reduce(
        [
            matched["x_win_px"],
            matched["y_win_px"],
            9216.0 - matched["x_win_px"],
            9232.0 - matched["y_win_px"],
        ]
    )

    condition = conditions_table(args.l2_manifest)
    if not condition.empty:
        matched = matched.merge(condition, on=["night", "source_file"], how="left", validate="many_to_one")

    file_status = pd.read_csv(args.file_status, dtype={"night": "string"})
    nights = sorted(
        file_status.loc[
            (file_status["product"] == "known_all")
            & file_status["status"].astype(str).str.startswith("ok"),
            "night",
        ].astype(str).str.zfill(8).unique()
    )
    all_dataset = ds.dataset(products / "known_all.parquet", format="parquet")
    parquet_temp = output / "known_recovery_by_detection.parquet.inprogress"
    parquet_final = output / "known_recovery_by_detection.parquet"
    csv_temp = output / "known_recovery_by_detection.csv.inprogress"
    csv_final = output / "known_recovery_by_detection.csv"
    parquet_writer = None
    csv_header = True
    binned_rows: list[dict] = []
    nightly_rows: list[dict] = []
    matched_rate_parts: list[pd.DataFrame] = []
    condition_finite_counts: dict[str, int] = {}
    all_duplicate_keys = 0
    prediction_rows = 0
    source_prediction_rows = 0
    try:
        for index, night in enumerate(nights, 1):
            table = all_dataset.to_table(filter=ds.field("night") == night)
            frame = table.to_pandas()
            frame["night"] = frame["night"].astype(str).str.zfill(8)
            frame["object_key"] = object_key(frame)
            frame["prediction_key"] = key_series(frame)
            duplicates = int(frame.duplicated("prediction_key").sum())
            all_duplicate_keys += duplicates
            source_prediction_rows += len(frame)
            frame["prediction_duplicate_rank"] = frame.groupby(
                "prediction_key", sort=False
            ).cumcount()
            frame["prediction_duplicate_n"] = frame.groupby(
                "prediction_key", sort=False
            )["prediction_key"].transform("size")
            frame = frame.loc[frame["prediction_duplicate_rank"] == 0].copy()
            frame["rate_arcsec_per_hour"] = angular_speed(frame)
            matched_keys = set(matched.loc[matched["night"] == night, "prediction_key"])
            frame["matched_1arcsec"] = frame["prediction_key"].isin(matched_keys)
            if not condition.empty:
                frame = frame.merge(
                    condition, on=["night", "source_file"], how="left", validate="many_to_one"
                )
            write_columns = [
                "night",
                "source_file",
                "object_key",
                "asteroid_name",
                "asteroid_number",
                "epoch_mjd",
                "pred_ra_deg",
                "pred_dec_deg",
                "pred_v_mag",
                "rate_arcsec_per_hour",
                "matched_1arcsec",
                "prediction_duplicate_n",
            ]
            write_columns += [
                column
                for column in (
                    "seeing_arcsec",
                    "sky_adu",
                    "sky_rms_adu",
                    "sky_mag",
                    "limiting_mag",
                    "n_sources",
                    "center_ra_deg",
                    "center_dec_deg",
                    "moon_separation_deg",
                    "moon_illumination_fraction",
                )
                if column in frame
            ]
            write_frame = frame[write_columns].copy()
            arrow = pa.Table.from_pandas(write_frame, preserve_index=False)
            if parquet_writer is None:
                parquet_writer = pq.ParquetWriter(
                    parquet_temp, arrow.schema, compression="zstd", compression_level=6
                )
            parquet_writer.write_table(arrow)
            if args.write_detection_csv:
                write_frame.to_csv(csv_temp, mode="w" if csv_header else "a", header=csv_header, index=False)
                csv_header = False

            binned_rows.extend(binned_recovery(frame, "pred_v_mag", V_BINS, "predicted_v_mag"))
            binned_rows.extend(
                binned_recovery(frame, "rate_arcsec_per_hour", RATE_BINS, "sky_rate_arcsec_per_hour")
            )
            condition_specs = [
                ("seeing_arcsec", SEEING_BINS, "seeing_arcsec"),
                ("n_sources", SOURCE_COUNT_BINS, "l2_catalog_source_count"),
                ("sky_adu", SKY_ADU_BINS, "sky_background_adu"),
                ("limiting_mag", LIMITING_MAG_BINS, "header_limiting_magnitude"),
                ("moon_separation_deg", MOON_SEPARATION_BINS, "moon_separation_deg"),
                ("moon_illumination_fraction", MOON_ILLUMINATION_BINS, "moon_illumination_fraction"),
            ]
            for variable, bins, label in condition_specs:
                if variable in frame:
                    binned_rows.extend(binned_recovery(frame, variable, bins, label))
                    condition_finite_counts[variable] = condition_finite_counts.get(variable, 0) + int(
                        pd.to_numeric(frame[variable], errors="coerce").notna().sum()
                    )
            nightly_rows.append(
                {
                    "night": night,
                    "source_prediction_row_n": int(
                        len(frame) + int((frame["prediction_duplicate_n"] - 1).sum())
                    ),
                    "predicted_nominal_n": len(frame),
                    "duplicate_prediction_row_n": int(
                        (frame["prediction_duplicate_n"] - 1).sum()
                    ),
                    "matched_1arcsec_n": int(frame["matched_1arcsec"].sum()),
                    "recovery_fraction": float(frame["matched_1arcsec"].mean()),
                    "predicted_rate_finite_n": int(frame["rate_arcsec_per_hour"].notna().sum()),
                }
            )
            matched_rate_parts.append(
                frame.loc[frame["matched_1arcsec"], ["night", "prediction_key", "rate_arcsec_per_hour"]]
            )
            prediction_rows += len(frame)
            if index % 10 == 0 or index == len(nights):
                print(f"[known] nights={index}/{len(nights)} predictions={prediction_rows}", flush=True)
        parquet_writer.close()
        parquet_writer = None
        os.replace(parquet_temp, parquet_final)
        if args.write_detection_csv:
            os.replace(csv_temp, csv_final)
    finally:
        if parquet_writer is not None:
            parquet_writer.close()
        if parquet_temp.exists():
            parquet_temp.unlink()
        if csv_temp.exists():
            csv_temp.unlink()

    matched_rates = pd.concat(matched_rate_parts, ignore_index=True)
    matched = matched.merge(
        matched_rates, on=["night", "prediction_key"], how="left", validate="one_to_one"
    )
    residual_columns = [
        "night",
        "source_file",
        "object_key",
        "epoch_mjd",
        "pred_ra_deg",
        "pred_dec_deg",
        "obs_ra_win_deg",
        "obs_dec_win_deg",
        "dra_cosdec_arcsec",
        "ddec_arcsec",
        "radial_residual_arcsec",
        "pred_v_mag",
        "mag_kron",
        "mag_aper4",
        "mag_psf",
        "rate_arcsec_per_hour",
        "x_win_px",
        "y_win_px",
        "detector_edge_distance_px",
        "fwhm_px",
        "flag",
        "flag_iso_num",
    ]
    residual_columns += [
        column
        for column in (
            "seeing_arcsec",
            "sky_adu",
            "sky_rms_adu",
            "sky_mag",
            "limiting_mag",
            "n_sources",
            "moon_separation_deg",
            "moon_illumination_fraction",
        )
        if column in matched
    ]
    residual = matched[residual_columns].copy()
    residual.to_parquet(output / "known_astrometric_residuals.parquet", index=False)
    Table.from_pandas(residual).write(
        output / "known_astrometric_residuals.fits", overwrite=True
    )
    fits_matched = matched.copy()
    for column in fits_matched.columns:
        if isinstance(fits_matched[column].dtype, pd.StringDtype) or fits_matched[column].dtype == object:
            fits_matched[column] = fits_matched[column].fillna("").astype(str)
    Table.from_pandas(fits_matched).write(
        output / "known_matches_quality_selected.fits", overwrite=True
    )

    binned = pd.DataFrame(binned_rows)
    binned = (
        binned.groupby(["metric", "lower", "upper"], dropna=False, as_index=False)
        .agg(
            denominator_prediction_n=("denominator_prediction_n", "sum"),
            matched_1arcsec_n=("matched_1arcsec_n", "sum"),
            denominator_definition=("denominator_definition", "first"),
        )
    )
    intervals = [
        wilson_interval(int(row.matched_1arcsec_n), int(row.denominator_prediction_n))
        for row in binned.itertuples(index=False)
    ]
    binned["recovery_fraction"] = binned["matched_1arcsec_n"] / binned["denominator_prediction_n"]
    binned["ci95_low"] = [item[0] for item in intervals]
    binned["ci95_high"] = [item[1] for item in intervals]
    if args.random_shift and args.random_shift.exists():
        random_shift = pd.read_csv(args.random_shift)
        random_shift.to_csv(output / "known_random_shift_test.csv", index=False)
    binned.to_csv(output / "known_recovery_binned.csv", index=False)
    pd.DataFrame(nightly_rows).to_csv(output / "known_recovery_by_night.csv", index=False)

    residual_binned = []
    residual_binned.extend(residual_bins(residual, "pred_v_mag", V_BINS, "predicted_v_mag"))
    residual_binned.extend(
        residual_bins(residual, "mag_aper4", MEASURED_MAG_BINS, "measured_mag_aper4")
    )
    residual_binned.extend(
        residual_bins(residual, "rate_arcsec_per_hour", RATE_BINS, "sky_rate_arcsec_per_hour")
    )
    residual_binned.extend(
        residual_bins(residual, "detector_edge_distance_px", EDGE_BINS, "detector_edge_distance_px")
    )
    for variable, bins, label in [
        ("seeing_arcsec", SEEING_BINS, "seeing_arcsec"),
        ("n_sources", SOURCE_COUNT_BINS, "l2_catalog_source_count"),
        ("sky_adu", SKY_ADU_BINS, "sky_background_adu"),
        ("limiting_mag", LIMITING_MAG_BINS, "header_limiting_magnitude"),
        ("moon_separation_deg", MOON_SEPARATION_BINS, "moon_separation_deg"),
        ("moon_illumination_fraction", MOON_ILLUMINATION_BINS, "moon_illumination_fraction"),
    ]:
        if variable in residual:
            residual_binned.extend(residual_bins(residual, variable, bins, label))
    pd.DataFrame(residual_binned).to_csv(output / "known_astrometric_residuals_binned.csv", index=False)

    radial = residual["radial_residual_arcsec"].to_numpy(dtype=float)
    dra = residual["dra_cosdec_arcsec"].to_numpy(dtype=float)
    ddec = residual["ddec_arcsec"].to_numpy(dtype=float)
    unique_objects = matched["object_key"].nunique()
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "prediction_denominator_definition": "V<22.5 and nominal WCS-projected detector rectangle; not a true detectable denominator",
        "source_prediction_rows_including_duplicates": int(source_prediction_rows),
        "predicted_nominal_n": int(prediction_rows),
        "matched_1arcsec_n": int(len(matched)),
        "matched_unique_identity_keys": int(unique_objects),
        "matched_prediction_key_duplicates": duplicate_matched,
        "prediction_key_duplicates": all_duplicate_keys,
        "match_fraction_nominal": float(len(matched) / prediction_rows),
        "median_radial_residual_arcsec": float(np.nanmedian(radial)),
        "p84_radial_residual_arcsec": float(np.nanpercentile(radial, 84)),
        "p95_radial_residual_arcsec": float(np.nanpercentile(radial, 95)),
        "median_dra_cosdec_arcsec": float(np.nanmedian(dra)),
        "median_ddec_arcsec": float(np.nanmedian(ddec)),
        "robust_sigma_dra_arcsec": float(0.7413 * (np.nanpercentile(dra, 75) - np.nanpercentile(dra, 25))),
        "robust_sigma_ddec_arcsec": float(0.7413 * (np.nanpercentile(ddec, 75) - np.nanpercentile(ddec, 25))),
        "radial_distribution_truncated_at_arcsec": 1.0,
        "condition_finite_prediction_counts": condition_finite_counts,
        "moon_context_site_status": "provisional_scheduler_known_coordinates_117.575E_40.393N_960m",
        "iers_status": "bundled_table_auto_download_disabled; out-of-range UT1 uses Astropy degraded-accuracy fallback",
        "input_hashes": {
            "known_all.parquet": sha256(products / "known_all.parquet"),
            "known_matched.parquet": sha256(products / "known_matched.parquet"),
        },
        "outputs": {
            "known_recovery_by_detection.parquet_sha256": sha256(parquet_final),
            "known_astrometric_residuals.fits_sha256": sha256(output / "known_astrometric_residuals.fits"),
        },
    }
    if args.l2_manifest and args.l2_manifest.exists():
        summary["input_hashes"]["l2_manifest.csv"] = sha256(args.l2_manifest)
    if args.write_detection_csv:
        summary["outputs"]["known_recovery_by_detection.csv_sha256"] = sha256(csv_final)
    (output / "known_population_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
