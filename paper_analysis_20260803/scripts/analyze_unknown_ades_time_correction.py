#!/usr/bin/env python3
"""Audit the unknown-ADES exposure-start to midpoint time correction.

This script is analysis-only.  It reads the frozen retained-detection table and
the final frozen L2 manifest, joins them exactly on ``night`` plus
``catalog_name == file_name``, and writes a correction audit table and JSON
summary.  It does not create ADES/MPC payloads and does not contact any service.

Exposure fallback policy is explicit: by default no fallback is permitted.  A
fallback is used only when ``--fallback-exposure-s`` is supplied, and every
fallback row is labelled in the CSV and counted in the JSON summary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DETECTIONS = (
    PROJECT_ROOT / "snapshot" / "review_sample" / "unknown_high_confidence_detections.csv"
)
DEFAULT_L2_MANIFEST = PROJECT_ROOT / "snapshot" / "inventory" / "l2_manifest.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "snapshot" / "derived_unknown"
OUTPUT_CSV_NAME = "unknown_ades_midpoint_correction_179.csv"
OUTPUT_JSON_NAME = "unknown_ades_midpoint_correction_summary.json"


DETECTION_REQUIRED = {
    "night",
    "trk_sub",
    "linkage_id",
    "detection_index",
    "image_name",
    "catalog_name",
    "catalog_path",
    "MJD",
    "lin_speed_arcsec_per_day",
    "detection_key",
    "final_paper_status",
}
MANIFEST_REQUIRED = {
    "night",
    "file_name",
    "path",
    "strict_standard_catalog",
    "mjd",
    "exposure_s",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_bool(values: pd.Series, name: str) -> pd.Series:
    normalized = values.astype("string").str.strip().str.lower()
    mapping = {
        "true": True,
        "1": True,
        "yes": True,
        "y": True,
        "t": True,
        "false": False,
        "0": False,
        "no": False,
        "n": False,
        "f": False,
    }
    invalid = normalized.isna() | ~normalized.isin(mapping)
    if invalid.any():
        raise ValueError(f"{name} contains invalid booleans: {values[invalid].head().tolist()}")
    return normalized.map(mapping).astype(bool)


def normalize_night(values: pd.Series, name: str) -> pd.Series:
    nights = values.astype("string").str.strip().str.replace(r"\.0$", "", regex=True).str.zfill(8)
    invalid = nights.isna() | ~nights.str.fullmatch(r"\d{8}", na=False)
    if invalid.any():
        raise ValueError(f"{name} contains invalid night keys: {values[invalid].head().tolist()}")
    return nights


def validate_exact_names(values: pd.Series, name: str) -> pd.Series:
    names = values.astype("string")
    invalid = names.isna() | names.eq("") | names.ne(names.str.strip())
    if invalid.any():
        raise ValueError(
            f"{name} must be non-empty and have no surrounding whitespace: "
            f"{values[invalid].head().tolist()}"
        )
    return names


def numeric(values: pd.Series, name: str, *, positive: bool = False) -> pd.Series:
    parsed = pd.to_numeric(values, errors="coerce")
    invalid = ~np.isfinite(parsed)
    if positive:
        invalid |= parsed <= 0
    if invalid.any():
        raise ValueError(f"{name} contains invalid numeric values: {values[invalid].head().tolist()}")
    return parsed.astype(float)


def summary_stats(values: pd.Series | np.ndarray) -> dict[str, float | None]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if not len(array):
        return {"p50": None, "p90": None, "max": None}
    p50, p90 = np.percentile(array, [50, 90])
    return {"p50": float(p50), "p90": float(p90), "max": float(np.max(array))}


def load_inputs(detections_path: Path, l2_manifest_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not detections_path.is_file():
        raise FileNotFoundError(f"Frozen retained-detection table not found: {detections_path}")
    if not l2_manifest_path.is_file():
        raise FileNotFoundError(
            f"Final frozen L2 manifest not found: {l2_manifest_path}\n"
            "The 30 s fallback cannot replace a missing manifest because exact "
            "night+file association must be established first. Re-run after the "
            "final inventory/l2_manifest.csv is frozen."
        )
    detections = pd.read_csv(
        detections_path,
        dtype={"night": "string", "trk_sub": "string"},
        low_memory=False,
    )
    manifest = pd.read_csv(
        l2_manifest_path,
        dtype={"night": "string", "file_name": "string"},
        low_memory=False,
    )
    missing_detection = sorted(DETECTION_REQUIRED - set(detections.columns))
    missing_manifest = sorted(MANIFEST_REQUIRED - set(manifest.columns))
    if missing_detection:
        raise ValueError(
            f"{detections_path} lacks required columns: {', '.join(missing_detection)}"
        )
    if missing_manifest:
        raise ValueError(
            f"{l2_manifest_path} lacks required columns: {', '.join(missing_manifest)}"
        )
    return detections, manifest


def analyze(
    detections: pd.DataFrame,
    manifest: pd.DataFrame,
    *,
    detections_path: Path,
    manifest_path: Path,
    expected_rows: int,
    fallback_exposure_s: float | None,
    max_start_disagreement_s: float,
    minimum_speed_arcsec_per_hour: float,
    maximum_speed_arcsec_per_hour: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    if expected_rows <= 0:
        raise ValueError("expected_rows must be positive")
    if len(detections) != expected_rows:
        raise ValueError(
            f"Frozen retained-detection row count is {len(detections):,}; "
            f"expected {expected_rows:,}. No output was written."
        )
    if fallback_exposure_s is not None and (
        not np.isfinite(fallback_exposure_s) or fallback_exposure_s <= 0
    ):
        raise ValueError("fallback_exposure_s must be positive when supplied")
    if max_start_disagreement_s < 0:
        raise ValueError("max_start_disagreement_s must be non-negative")
    if minimum_speed_arcsec_per_hour <= 0 or maximum_speed_arcsec_per_hour <= 0:
        raise ValueError("speed bounds must be positive")
    if minimum_speed_arcsec_per_hour > maximum_speed_arcsec_per_hour:
        raise ValueError("minimum speed cannot exceed maximum speed")

    detections = detections.copy()
    manifest = manifest.copy()
    detections["night"] = normalize_night(detections["night"], "detection night")
    manifest["night"] = normalize_night(manifest["night"], "manifest night")
    detections["catalog_name"] = validate_exact_names(
        detections["catalog_name"], "detection catalog_name"
    )
    manifest["file_name"] = validate_exact_names(manifest["file_name"], "manifest file_name")
    detections["detection_key"] = validate_exact_names(
        detections["detection_key"], "detection_key"
    )
    if detections["detection_key"].duplicated().any():
        duplicates = detections.loc[
            detections["detection_key"].duplicated(False), "detection_key"
        ].head().tolist()
        raise ValueError(f"retained detection_key is not unique: {duplicates}")
    retained_status = "retained_after_posthoc_audit"
    unexpected_status = ~detections["final_paper_status"].astype("string").eq(retained_status)
    if unexpected_status.any():
        counts = detections.loc[unexpected_status, "final_paper_status"].value_counts(
            dropna=False
        ).to_dict()
        raise ValueError(
            "High-confidence detection input contains rows outside the retained "
            f"post-audit population: {counts}"
        )

    strict = parse_bool(manifest["strict_standard_catalog"], "strict_standard_catalog")
    strict_manifest = manifest.loc[strict].copy()
    duplicate_manifest = strict_manifest.duplicated(["night", "file_name"], keep=False)
    if duplicate_manifest.any():
        examples = strict_manifest.loc[
            duplicate_manifest, ["night", "file_name", "path"]
        ].head().to_dict("records")
        raise ValueError(f"strict L2 manifest keys are not unique: {examples}")

    manifest_columns = strict_manifest[
        ["night", "file_name", "path", "mjd", "exposure_s"]
    ].rename(
        columns={
            "file_name": "manifest_file_name",
            "path": "manifest_path",
            "mjd": "manifest_start_mjd_utc",
            "exposure_s": "manifest_exposure_s",
        }
    )
    joined = detections.merge(
        manifest_columns,
        left_on=["night", "catalog_name"],
        right_on=["night", "manifest_file_name"],
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    unmatched = joined.loc[joined["_merge"].ne("both"), ["night", "catalog_name"]]
    unmatched_files = (
        unmatched.drop_duplicates().sort_values(["night", "catalog_name"]).to_dict("records")
    )
    if unmatched_files:
        preview = ", ".join(
            f"{row['night']}::{row['catalog_name']}" for row in unmatched_files[:8]
        )
        raise ValueError(
            f"Exact night+catalog/file_name association failed for {len(unmatched_files)} "
            f"unique L2 files ({len(unmatched)} detection rows): {preview}. "
            "No output was written."
        )
    joined.drop(columns="_merge", inplace=True)

    production_start = numeric(joined["MJD"], "retained detection MJD", positive=True)
    manifest_start = numeric(
        joined["manifest_start_mjd_utc"], "manifest L2 start MJD", positive=True
    )
    joined["production_ades_start_mjd_utc"] = production_start
    joined["l2_manifest_start_mjd_utc"] = manifest_start
    joined["start_mjd_disagreement_s"] = (production_start - manifest_start) * 86400.0
    disagreement = joined["start_mjd_disagreement_s"].abs() > max_start_disagreement_s
    if disagreement.any():
        examples = joined.loc[
            disagreement,
            ["night", "catalog_name", "start_mjd_disagreement_s"],
        ].head().to_dict("records")
        raise ValueError(
            "Detection MJD and exact-matched L2 manifest MJD disagree beyond "
            f"{max_start_disagreement_s:g} s: {examples}. No output was written."
        )

    exposure = pd.to_numeric(joined["manifest_exposure_s"], errors="coerce")
    exposure_missing = ~np.isfinite(exposure) | (exposure <= 0)
    fallback_count = int(exposure_missing.sum())
    if fallback_count and fallback_exposure_s is None:
        examples = joined.loc[
            exposure_missing, ["night", "catalog_name", "manifest_exposure_s"]
        ].head().to_dict("records")
        raise ValueError(
            f"{fallback_count} matched detections lack a positive L2 exposure_s: {examples}. "
            "No fallback was used and no output was written. To authorize an explicit "
            "30 s fallback, rerun with --fallback-exposure-s 30."
        )
    joined["exposure_source"] = "l2_manifest"
    if fallback_count:
        exposure.loc[exposure_missing] = float(fallback_exposure_s)
        joined.loc[exposure_missing, "exposure_source"] = "explicit_cli_fallback"
    joined["exposure_s"] = exposure.astype(float)
    joined["midpoint_correction_s"] = joined["exposure_s"] / 2.0
    joined["corrected_midpoint_mjd_utc"] = (
        joined["production_ades_start_mjd_utc"]
        + joined["midpoint_correction_s"] / 86400.0
    )
    corrected_time = Time(
        joined["corrected_midpoint_mjd_utc"].to_numpy(),
        format="mjd",
        scale="utc",
        precision=6,
    )
    joined["corrected_midpoint_isot_utc"] = corrected_time.isot

    link_speed = numeric(
        joined["lin_speed_arcsec_per_day"],
        "lin_speed_arcsec_per_day",
        positive=True,
    ) / 24.0
    joined["link_linear_speed_arcsec_per_hour"] = link_speed
    joined["link_speed_within_production_3_to_63_arcsec_per_hour"] = link_speed.between(
        minimum_speed_arcsec_per_hour,
        maximum_speed_arcsec_per_hour,
        inclusive="both",
    )
    hours = joined["midpoint_correction_s"] / 3600.0
    joined["equivalent_shift_at_min_speed_arcsec"] = hours * minimum_speed_arcsec_per_hour
    joined["equivalent_shift_at_max_speed_arcsec"] = hours * maximum_speed_arcsec_per_hour
    joined["equivalent_shift_at_link_speed_arcsec"] = hours * link_speed

    identity_columns = [
        "night",
        "trk_sub",
        "linkage_id",
        "detection_index",
        "image_name",
        "catalog_name",
        "catalog_path",
        "manifest_file_name",
        "manifest_path",
        "detection_key",
        "final_paper_status",
    ]
    audit_columns = [
        "production_ades_start_mjd_utc",
        "l2_manifest_start_mjd_utc",
        "start_mjd_disagreement_s",
        "manifest_exposure_s",
        "exposure_s",
        "exposure_source",
        "midpoint_correction_s",
        "corrected_midpoint_mjd_utc",
        "corrected_midpoint_isot_utc",
        "lin_speed_arcsec_per_day",
        "link_linear_speed_arcsec_per_hour",
        "link_speed_within_production_3_to_63_arcsec_per_hour",
        "equivalent_shift_at_min_speed_arcsec",
        "equivalent_shift_at_max_speed_arcsec",
        "equivalent_shift_at_link_speed_arcsec",
    ]
    output = joined[identity_columns + audit_columns].copy()
    output.sort_values(
        ["night", "trk_sub", "linkage_id", "detection_index"], inplace=True
    )
    output.reset_index(drop=True, inplace=True)

    summary: dict[str, object] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "unknown ADES exposure-start to exposure-midpoint correction",
        "correction_direction": "add half the matched L2 exposure duration to production ADES MJD",
        "row_count": int(len(output)),
        "expected_row_count": int(expected_rows),
        "unique_detection_count": int(output["detection_key"].nunique()),
        "unique_l2_file_count": int(output[["night", "catalog_name"]].drop_duplicates().shape[0]),
        "unique_single_night_linkage_count": int(
            output[["night", "trk_sub", "linkage_id"]].drop_duplicates().shape[0]
        ),
        "input_profile": {
            "detection_rows": int(len(detections)),
            "detection_columns": int(len(detections.columns)),
            "l2_manifest_rows": int(len(manifest)),
            "strict_l2_manifest_rows": int(len(strict_manifest)),
            "retained_status_required": retained_status,
        },
        "exact_join": {
            "left_keys": ["night", "catalog_name"],
            "right_keys": ["night", "file_name"],
            "unmatched_detection_rows": 0,
            "unmatched_unique_files": [],
            "duplicate_strict_manifest_keys": 0,
        },
        "fallback_exposure": {
            "enabled": fallback_exposure_s is not None,
            "requested_value_s": (
                None if fallback_exposure_s is None else float(fallback_exposure_s)
            ),
            "used_row_count": fallback_count,
            "policy": "disabled unless explicitly supplied on CLI",
        },
        "midpoint_correction_s": summary_stats(output["midpoint_correction_s"]),
        "start_mjd_disagreement_abs_s": summary_stats(
            output["start_mjd_disagreement_s"].abs()
        ),
        "production_speed_bounds_arcsec_per_hour": {
            "minimum": float(minimum_speed_arcsec_per_hour),
            "maximum": float(maximum_speed_arcsec_per_hour),
            "rows_outside_bounds": int(
                (~output["link_speed_within_production_3_to_63_arcsec_per_hour"]).sum()
            ),
        },
        "equivalent_along_track_shift_arcsec": {
            f"at_{minimum_speed_arcsec_per_hour:g}_arcsec_per_hour": summary_stats(
                output["equivalent_shift_at_min_speed_arcsec"]
            ),
            f"at_{maximum_speed_arcsec_per_hour:g}_arcsec_per_hour": summary_stats(
                output["equivalent_shift_at_max_speed_arcsec"]
            ),
            "at_each_link_linear_speed": summary_stats(
                output["equivalent_shift_at_link_speed_arcsec"]
            ),
        },
        "input_hashes": {
            "unknown_high_confidence_detections.csv": sha256(detections_path),
            "l2_manifest.csv": sha256(manifest_path),
        },
        "safety": {
            "mpc_or_ades_file_generated": False,
            "network_access_performed": False,
            "production_files_modified": False,
        },
    }
    return output, summary


def output_paths(output_dir: Path) -> tuple[Path, Path]:
    return output_dir / OUTPUT_CSV_NAME, output_dir / OUTPUT_JSON_NAME


def refuse_existing(csv_path: Path, json_path: Path) -> None:
    existing = [str(path) for path in [csv_path, json_path] if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing audit output(s): " + ", ".join(existing)
        )


def write_outputs(
    output: pd.DataFrame,
    summary: dict[str, object],
    output_dir: Path,
) -> tuple[Path, Path]:
    csv_path, json_path = output_paths(output_dir)
    refuse_existing(csv_path, json_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Exclusive creation preserves the refusal-to-overwrite contract even if
    # another process creates an output after the preflight check.
    with csv_path.open("x", encoding="utf-8", newline="") as handle:
        output.to_csv(handle, index=False)
    try:
        with json_path.open("x", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except Exception:
        # Leave the CSV in place rather than silently deleting an audit artifact;
        # the next run will refuse overwrite and expose the interrupted state.
        raise
    return csv_path, json_path


def run_self_test() -> None:
    """Synthetic smoke test including explicit fallback and overwrite refusal."""

    with tempfile.TemporaryDirectory(prefix="unknown_ades_time_smoke.") as tmp:
        directory = Path(tmp)
        detections_path = directory / "detections.csv"
        manifest_path = directory / "l2_manifest.csv"
        output_dir = directory / "out"
        base_mjd = 61000.25
        detections = pd.DataFrame(
            {
                "night": ["20251121", "20251121", "20251122"],
                "trk_sub": ["0001", "0001", "0002"],
                "linkage_id": [1, 1, 2],
                "detection_index": [0, 1, 0],
                "image_name": ["A.fits.gz", "B.fits.gz", "C.fits.gz"],
                "catalog_name": ["A_cat.fits.gz", "B_cat.fits.gz", "C_cat.fits.gz"],
                "catalog_path": ["/frozen/A", "/frozen/B", "/frozen/C"],
                "MJD": [base_mjd, base_mjd + 0.01, base_mjd + 1.0],
                "lin_speed_arcsec_per_day": [72.0, 720.0, 1512.0],
                "detection_key": ["A::1", "B::2", "C::3"],
                "final_paper_status": ["retained_after_posthoc_audit"] * 3,
            }
        )
        manifest = pd.DataFrame(
            {
                "night": ["20251121", "20251121", "20251122"],
                "file_name": ["A_cat.fits.gz", "B_cat.fits.gz", "C_cat.fits.gz"],
                "path": ["/frozen/A", "/frozen/B", "/frozen/C"],
                "strict_standard_catalog": [True, True, True],
                "mjd": [base_mjd, base_mjd + 0.01, base_mjd + 1.0],
                "exposure_s": [30.0, 32.0, np.nan],
            }
        )
        detections.to_csv(detections_path, index=False)
        manifest.to_csv(manifest_path, index=False)
        loaded_detections, loaded_manifest = load_inputs(detections_path, manifest_path)
        try:
            analyze(
                loaded_detections,
                loaded_manifest,
                detections_path=detections_path,
                manifest_path=manifest_path,
                expected_rows=3,
                fallback_exposure_s=None,
                max_start_disagreement_s=0.1,
                minimum_speed_arcsec_per_hour=3.0,
                maximum_speed_arcsec_per_hour=63.0,
            )
        except ValueError as exc:
            if "No fallback was used" not in str(exc):
                raise AssertionError("self-test did not exercise no-fallback refusal") from exc
        else:
            raise AssertionError("self-test expected missing-exposure refusal")

        output, summary = analyze(
            loaded_detections,
            loaded_manifest,
            detections_path=detections_path,
            manifest_path=manifest_path,
            expected_rows=3,
            fallback_exposure_s=30.0,
            max_start_disagreement_s=0.1,
            minimum_speed_arcsec_per_hour=3.0,
            maximum_speed_arcsec_per_hour=63.0,
        )
        if output["midpoint_correction_s"].tolist() != [15.0, 16.0, 15.0]:
            raise AssertionError("midpoint offsets are incorrect")
        if summary["fallback_exposure"]["used_row_count"] != 1:
            raise AssertionError("fallback use was not recorded")
        csv_path, json_path = write_outputs(output, summary, output_dir)
        if not csv_path.is_file() or not json_path.is_file():
            raise AssertionError("self-test outputs were not written")
        try:
            write_outputs(output, summary, output_dir)
        except FileExistsError:
            pass
        else:
            raise AssertionError("overwrite refusal did not trigger")
    print("synthetic self-test passed (no fallback by default; explicit 30 s fallback exercised)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--detections", type=Path, default=DEFAULT_DETECTIONS)
    parser.add_argument("--l2-manifest", type=Path, default=DEFAULT_L2_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--expected-rows", type=int, default=179)
    parser.add_argument(
        "--fallback-exposure-s",
        type=float,
        help="Explicit fallback for matched rows with missing exposure_s; disabled by default",
    )
    parser.add_argument("--max-start-disagreement-s", type=float, default=0.1)
    parser.add_argument("--minimum-speed-arcsec-per-hour", type=float, default=3.0)
    parser.add_argument("--maximum-speed-arcsec-per-hour", type=float, default=63.0)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_test()
        return
    output_dir = args.output_dir.expanduser().resolve(strict=False)
    csv_path, json_path = output_paths(output_dir)
    refuse_existing(csv_path, json_path)
    detections_path = args.detections.expanduser().resolve(strict=False)
    manifest_path = args.l2_manifest.expanduser().resolve(strict=False)
    detections, manifest = load_inputs(detections_path, manifest_path)
    output, summary = analyze(
        detections,
        manifest,
        detections_path=detections_path,
        manifest_path=manifest_path,
        expected_rows=args.expected_rows,
        fallback_exposure_s=args.fallback_exposure_s,
        max_start_disagreement_s=args.max_start_disagreement_s,
        minimum_speed_arcsec_per_hour=args.minimum_speed_arcsec_per_hour,
        maximum_speed_arcsec_per_hour=args.maximum_speed_arcsec_per_hour,
    )
    csv_path, json_path = write_outputs(output, summary, output_dir)
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(
        "fallback exposure: "
        f"enabled={summary['fallback_exposure']['enabled']}, "
        f"used_rows={summary['fallback_exposure']['used_row_count']}"
    )
    print("No ADES/MPC submission file was generated.")


if __name__ == "__main__":
    main()
