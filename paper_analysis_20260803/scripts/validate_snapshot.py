#!/usr/bin/env python3
"""Validate closure of the frozen 2026-08-03 PASP analysis snapshot.

The expected headline values in this module are acceptance targets only.  Every
reported observation is recomputed from a row-level manifest, Parquet metadata,
or a row-level derived table; expected values are never substituted for missing
or failed computations.  An unavailable or not-yet-approved artifact is marked
``SKIP``.  A present, readable artifact that contradicts an expectation is
marked ``FAIL``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import struct
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    import pandas as pd
except ImportError:  # pragma: no cover - exercised only in a deficient runtime
    pd = None

try:
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover - exercised only in a deficient runtime
    pq = None


EXPECTED_COUNTS: dict[str, int] = {
    "strict_raw_exposures": 41_074,
    "strict_raw_nights": 134,
    "strict_raw_fields": 1_430,
    "engineering_raw_exposures": 41_152,
    "engineering_raw_nights": 136,
    "strict_l2_catalogs": 40_399,
    "strict_l2_nights": 131,
    "known_prediction_rows": 13_311_871,
    "known_matched_1arcsec_rows": 534_780,
    "known_mask15_rows": 563_612,
    "unknown_catalog_linkages": 4_762,
    "unknown_memberships": 14_299,
    "unknown_unique_detections": 14_159,
    "unknown_nonempty_nights": 116,
    "unknown_true_zero_nights": 8,
    "initial_review_linkages": 67,
    "retained_review_linkages": 58,
    "rejected_review_linkages": 9,
    "retained_detection_memberships": 179,
}

FIGURE_STEMS = (
    "fig01_system_architecture",
    "fig02_footprint_coverage",
    "fig03_data_accounting",
    "fig04_scheduler_example",
    "fig05_known_method_and_residuals",
    "fig06_unknown_pipeline_examples",
    "fig07_nightly_exposure_history",
    "fig08_known_results",
    "fig09_unknown_funnel",
    "fig10_high_confidence_distributions",
    "fig11_twilight_context",
    "fig12_operations_timeline",
)

EXPECTED_STAGE_UNITS = {
    "l2_detection_n": "detection",
    "mag_flag_prefilter_detection_n": "detection",
    "gaia_survivor_n": "detection",
    "tracklet_n": "tracklet",
    "link_n": "linkage",
    "orbit_fit_n": "linkage",
    "orbit_is_good_n": "linkage",
    "orbit_fit_non_known_n": "linkage",
    "orbit_is_good_non_known_n": "linkage",
    "unknown_n": "linkage",
    "review_real_n": "linkage",
    "audit_real_n": "linkage",
}

RETAINED_STATUS = "retained_after_posthoc_audit"
REJECTED_STATUS = "rejected_posthoc"
VALID_STATUSES = {"PASS", "FAIL", "SKIP"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_value(value: Any) -> Any:
    """Convert common scalar/path containers to JSON-safe values."""

    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_value(item) for item in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


@dataclass
class CheckBook:
    checks: list[dict[str, Any]] = field(default_factory=list)

    def add(
        self,
        check_id: str,
        category: str,
        status: str,
        description: str,
        *,
        expected: Any = None,
        observed: Any = None,
        sources: Iterable[Path | str] = (),
        detail: str = "",
    ) -> None:
        if status not in VALID_STATUSES:
            raise ValueError(f"invalid check status: {status}")
        self.checks.append(
            {
                "id": check_id,
                "category": category,
                "status": status,
                "description": description,
                "expected": json_value(expected),
                "observed": json_value(observed),
                "sources": [str(source) for source in sources],
                "detail": detail,
            }
        )

    def metric(
        self,
        check_id: str,
        category: str,
        description: str,
        expected: int | float,
        observed: int | float,
        source: Path,
        *,
        detail: str = "Directly recomputed from the cited artifact.",
    ) -> None:
        self.add(
            check_id,
            category,
            "PASS" if observed == expected else "FAIL",
            description,
            expected=expected,
            observed=observed,
            sources=[source],
            detail=detail,
        )

    def skipped_metric(
        self,
        check_id: str,
        category: str,
        description: str,
        expected: Any,
        source: Path,
        detail: str,
    ) -> None:
        self.add(
            check_id,
            category,
            "SKIP",
            description,
            expected=expected,
            observed=None,
            sources=[source],
            detail=detail,
        )


def parse_bool_series(series: Any) -> Any:
    if pd is None:
        raise RuntimeError("pandas is unavailable")
    if str(series.dtype).lower() in {"bool", "boolean"}:
        return series.fillna(False).astype(bool)
    return (
        series.astype("string")
        .str.strip()
        .str.lower()
        .isin({"1", "true", "yes", "y", "t"})
    )


def parse_date_key(value: Any) -> date:
    text = str(value).strip().replace("-", "")
    return datetime.strptime(text, "%Y%m%d").date()


def normalized_link_keys(
    frame: Any,
    *,
    night_column: str = "night",
) -> tuple[set[str], int]:
    required = {night_column, "trk_sub", "linkage_id"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"missing link-key columns: {sorted(missing)}")
    values = frame[[night_column, "trk_sub", "linkage_id"]].copy()
    missing_rows = int(values.isna().any(axis=1).sum())
    values[night_column] = (
        values[night_column].astype("string").str.strip().str.zfill(8)
    )
    values["trk_sub"] = values["trk_sub"].astype("string").str.strip()
    numeric_id = pd.to_numeric(values["linkage_id"], errors="coerce")
    missing_rows += int(numeric_id.isna().sum())
    valid = values[night_column].notna() & values["trk_sub"].notna() & numeric_id.notna()
    keys = {
        f"{night}::{trk_sub}::{int(linkage_id)}"
        for night, trk_sub, linkage_id in zip(
            values.loc[valid, night_column],
            values.loc[valid, "trk_sub"],
            numeric_id.loc[valid],
        )
    }
    return keys, missing_rows


def normalized_detection_keys(frame: Any) -> tuple[set[str], int]:
    required = {"night", "image_name", "obj_id"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"missing detection-key columns: {sorted(missing)}")
    values = frame[["night", "image_name", "obj_id"]].copy()
    missing_rows = int(values.isna().any(axis=1).sum())
    for column in values:
        values[column] = values[column].astype("string").str.strip()
    values["night"] = values["night"].str.zfill(8)
    valid = values.notna().all(axis=1)
    keys = {
        f"{night}::{image_name}::{obj_id}"
        for night, image_name, obj_id in values.loc[valid].itertuples(index=False, name=None)
    }
    return keys, missing_rows


def load_csv(path: Path, **kwargs: Any) -> Any:
    if pd is None:
        raise RuntimeError("pandas is unavailable")
    return pd.read_csv(path, **kwargs)


def add_unavailable_metrics(
    book: CheckBook,
    metrics: Iterable[tuple[str, str, int]],
    category: str,
    source: Path,
    detail: str,
) -> None:
    for check_id, description, expected in metrics:
        book.skipped_metric(check_id, category, description, expected, source, detail)


def validate_config(book: CheckBook, path: Path) -> dict[str, Any] | None:
    if not path.exists():
        book.add(
            "config.snapshot",
            "configuration",
            "SKIP",
            "Snapshot configuration is present and readable.",
            expected="readable JSON configuration",
            observed=None,
            sources=[path],
            detail="Configuration has not been staged yet.",
        )
        return None
    try:
        config = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        book.add(
            "config.snapshot",
            "configuration",
            "FAIL",
            "Snapshot configuration is present and readable.",
            expected="readable JSON configuration",
            observed=type(exc).__name__,
            sources=[path],
            detail=str(exc),
        )
        return None
    identity_observed = {
        "snapshot_label": config.get("snapshot_label"),
        "schema_version": config.get("schema_version"),
    }
    identity_expected = {"snapshot_label": "2026-08-03", "schema_version": "1.0"}
    book.add(
        "config.snapshot",
        "configuration",
        "PASS" if identity_observed == identity_expected else "FAIL",
        "Snapshot configuration is present and readable.",
        expected=identity_expected,
        observed=identity_observed,
        sources=[path],
    )

    rules = config.get("quality_rules")
    required_rule_maps = (
        "exclude_from_primary_science",
        "quarantine_unknown",
        "raw_without_l2",
    )
    errors: list[str] = []
    observed: dict[str, Any] = {}
    if not isinstance(rules, dict):
        errors.append("quality_rules is not an object")
    else:
        try:
            start = parse_date_key(config["observation_start"])
            end = parse_date_key(config["observation_end"])
        except Exception as exc:
            errors.append(f"invalid observation interval: {exc}")
            start = end = None
        for name in required_rule_maps:
            mapping = rules.get(name)
            if not isinstance(mapping, dict):
                errors.append(f"{name} is not an object")
                continue
            observed[name] = sorted(str(key) for key in mapping)
            for night, reason in mapping.items():
                try:
                    parsed = parse_date_key(night)
                    if start is not None and not start <= parsed <= end:
                        errors.append(f"{name}:{night} lies outside the snapshot interval")
                except Exception:
                    errors.append(f"{name}:{night} is not YYYYMMDD")
                if not str(reason).strip():
                    errors.append(f"{name}:{night} has an empty reason")
        observed["status"] = rules.get("status")
    book.add(
        "quality.config_structure",
        "quality mask",
        "FAIL" if errors else "PASS",
        "Quality-mask categories, dates, and reasons are structurally valid.",
        expected=list(required_rule_maps),
        observed=observed,
        sources=[path],
        detail="; ".join(errors),
    )

    status = str(rules.get("status", "")) if isinstance(rules, dict) else ""
    signed_off = bool(status) and not any(
        token in status.lower() for token in ("provisional", "pending", "requires")
    )
    book.add(
        "quality.author_signoff",
        "quality mask",
        "PASS" if signed_off else "SKIP",
        "The frozen quality mask has author sign-off.",
        expected="author-approved frozen mask",
        observed=status or None,
        sources=[path],
        detail="A provisional mask is usable for engineering QA but is not final paper evidence."
        if not signed_off
        else "",
    )

    unknown = config.get("unknown", {})
    known = config.get("known", {})
    expected_units = {
        "known.official_match_radius_arcsec": 1.0,
        "known.unknown_mask_radius_arcsec": 1.5,
        "unknown.gaia_match_radius_arcsec": 1.5,
        "unknown.minimum_speed_arcsec_per_hour": 3.0,
        "unknown.maximum_speed_arcsec_per_hour": 63.0,
        "unknown.maximum_magnitude_difference": 1.0,
        "unknown.maximum_speed_difference_arcsec_per_hour": 5.0,
        "unknown.maximum_direction_difference_deg": 10.0,
        "unknown.standard_exposure_midpoint_offset_s": 15.0,
    }
    observed_units = {
        "known.official_match_radius_arcsec": known.get("official_match_radius_arcsec"),
        "known.unknown_mask_radius_arcsec": known.get("unknown_mask_radius_arcsec"),
        "unknown.gaia_match_radius_arcsec": unknown.get("gaia_match_radius_arcsec"),
        "unknown.minimum_speed_arcsec_per_hour": unknown.get(
            "minimum_speed_arcsec_per_hour"
        ),
        "unknown.maximum_speed_arcsec_per_hour": unknown.get(
            "maximum_speed_arcsec_per_hour"
        ),
        "unknown.maximum_magnitude_difference": unknown.get(
            "maximum_magnitude_difference"
        ),
        "unknown.maximum_speed_difference_arcsec_per_hour": unknown.get(
            "maximum_speed_difference_arcsec_per_hour"
        ),
        "unknown.maximum_direction_difference_deg": unknown.get(
            "maximum_direction_difference_deg"
        ),
        "unknown.standard_exposure_midpoint_offset_s": unknown.get(
            "standard_exposure_midpoint_offset_s"
        ),
    }
    unit_match = observed_units == expected_units
    book.add(
        "units.config_definitions",
        "units",
        "PASS" if unit_match else "FAIL",
        "Configured thresholds retain explicit physical units and snapshot values.",
        expected=expected_units,
        observed=observed_units,
        sources=[path],
        detail="The 15 s value is the standard 30 s exposure midpoint offset; it is a definition, not a correction applied here.",
    )
    time_expected = {
        "unknown_ades_time_as_executed": "exposure_start",
        "orbit_and_tracklet_time_reference": "exposure_midpoint",
        "standard_exposure_midpoint_offset_s": 15.0,
    }
    time_observed = {key: unknown.get(key) for key in time_expected}
    book.add(
        "units.unknown_time_reference",
        "units",
        "PASS" if time_observed == time_expected else "FAIL",
        "Unknown ADES and orbit/tracklet time references are explicitly distinguished.",
        expected=time_expected,
        observed=time_observed,
        sources=[path],
    )
    return config


def validate_inventory(
    book: CheckBook,
    inventory: Path,
    computed: dict[str, Any],
) -> None:
    raw_metrics = (
        ("inventory.strict_raw_exposures", "Strict standard-science raw exposures", EXPECTED_COUNTS["strict_raw_exposures"]),
        ("inventory.strict_raw_nights", "Nights containing strict standard-science raw exposures", EXPECTED_COUNTS["strict_raw_nights"]),
        ("inventory.strict_raw_fields", "Unique strict standard-science field identifiers", EXPECTED_COUNTS["strict_raw_fields"]),
    )
    raw_path = inventory / "raw_manifest.csv"
    if not raw_path.exists():
        add_unavailable_metrics(
            book,
            raw_metrics,
            "inventory",
            raw_path,
            "The strict raw manifest has not been staged.",
        )
    elif pd is None:
        add_unavailable_metrics(
            book, raw_metrics, "inventory", raw_path, "pandas is unavailable."
        )
    else:
        try:
            raw = load_csv(raw_path, dtype={"night": "string", "field_id": "string"})
            required = {"night", "field_id", "strict_standard_science"}
            missing = required.difference(raw.columns)
            if missing:
                raise ValueError(f"missing columns: {sorted(missing)}")
            strict = raw.loc[parse_bool_series(raw["strict_standard_science"])].copy()
            observed = {
                "strict_raw_exposures": int(len(strict)),
                "strict_raw_nights": int(strict["night"].astype("string").str.zfill(8).nunique()),
                "strict_raw_fields": int(strict["field_id"].astype("string").str.strip().nunique()),
            }
            for (check_id, description, expected), key in zip(raw_metrics, observed):
                book.metric(check_id, "inventory", description, expected, observed[key], raw_path)
                computed[key] = observed[key]
        except Exception as exc:
            book.add(
                "inventory.raw_manifest_readability",
                "inventory",
                "FAIL",
                "The raw manifest can be parsed at the required grain.",
                expected="night, field_id, strict_standard_science",
                observed=type(exc).__name__,
                sources=[raw_path],
                detail=str(exc),
            )
            add_unavailable_metrics(
                book,
                raw_metrics,
                "inventory",
                raw_path,
                "Metrics were not computed because the present manifest is invalid.",
            )

    engineering_metrics = (
        ("inventory.engineering_raw_exposures", "All engineering MP raw exposures", EXPECTED_COUNTS["engineering_raw_exposures"]),
        ("inventory.engineering_raw_nights", "Nights containing any engineering MP exposure", EXPECTED_COUNTS["engineering_raw_nights"]),
    )
    engineering_path = inventory / "raw_engineering_manifest.csv"
    if not engineering_path.exists():
        add_unavailable_metrics(
            book,
            engineering_metrics,
            "inventory",
            engineering_path,
            "The engineering raw manifest has not been staged.",
        )
    elif pd is None:
        add_unavailable_metrics(
            book,
            engineering_metrics,
            "inventory",
            engineering_path,
            "pandas is unavailable.",
        )
    else:
        try:
            engineering = load_csv(engineering_path, dtype={"night": "string"})
            if "night" not in engineering:
                raise ValueError("missing night column")
            observed = {
                "engineering_raw_exposures": int(len(engineering)),
                "engineering_raw_nights": int(
                    engineering["night"].astype("string").str.zfill(8).nunique()
                ),
            }
            for (check_id, description, expected), key in zip(
                engineering_metrics, observed
            ):
                book.metric(
                    check_id,
                    "inventory",
                    description,
                    expected,
                    observed[key],
                    engineering_path,
                )
                computed[key] = observed[key]
        except Exception as exc:
            book.add(
                "inventory.engineering_manifest_readability",
                "inventory",
                "FAIL",
                "The engineering raw manifest can be parsed at the required grain.",
                expected="night column",
                observed=type(exc).__name__,
                sources=[engineering_path],
                detail=str(exc),
            )
            add_unavailable_metrics(
                book,
                engineering_metrics,
                "inventory",
                engineering_path,
                "Metrics were not computed because the present manifest is invalid.",
            )

    l2_metrics = (
        ("inventory.strict_l2_catalogs", "Strict standard L2 catalogs", EXPECTED_COUNTS["strict_l2_catalogs"]),
        ("inventory.strict_l2_nights", "Nights containing strict standard L2 catalogs", EXPECTED_COUNTS["strict_l2_nights"]),
    )
    l2_path = inventory / "l2_manifest.csv"
    if not l2_path.exists():
        add_unavailable_metrics(
            book,
            l2_metrics,
            "inventory",
            l2_path,
            "The L2 manifest has not been staged; absence is not interpreted as zero catalogs.",
        )
    elif pd is None:
        add_unavailable_metrics(
            book, l2_metrics, "inventory", l2_path, "pandas is unavailable."
        )
    else:
        try:
            l2 = load_csv(l2_path, dtype={"night": "string"})
            required = {"night", "strict_standard_catalog"}
            missing = required.difference(l2.columns)
            if missing:
                raise ValueError(f"missing columns: {sorted(missing)}")
            strict = l2.loc[parse_bool_series(l2["strict_standard_catalog"])]
            observed = {
                "strict_l2_catalogs": int(len(strict)),
                "strict_l2_nights": int(
                    strict["night"].astype("string").str.zfill(8).nunique()
                ),
            }
            for (check_id, description, expected), key in zip(l2_metrics, observed):
                book.metric(check_id, "inventory", description, expected, observed[key], l2_path)
                computed[key] = observed[key]
        except Exception as exc:
            book.add(
                "inventory.l2_manifest_readability",
                "inventory",
                "FAIL",
                "The L2 manifest can be parsed at the strict-catalog grain.",
                expected="night, strict_standard_catalog",
                observed=type(exc).__name__,
                sources=[l2_path],
                detail=str(exc),
            )
            add_unavailable_metrics(
                book,
                l2_metrics,
                "inventory",
                l2_path,
                "Metrics were not computed because the present manifest is invalid.",
            )


def parquet_row_count(path: Path) -> int:
    if pq is None:
        raise RuntimeError("pyarrow is unavailable")
    return int(pq.ParquetFile(path).metadata.num_rows)


def inferred_unit(column: str) -> str | None:
    exact = {"mjd": "d", "mjds": "d", "epoch_mjd": "d"}
    if column in exact:
        return exact[column]
    suffixes = (
        ("_arcsec_per_day", "arcsec/day"),
        ("_au_day2", "AU/day2"),
        ("_au_day", "AU/day"),
        ("_kms", "km/s"),
        ("_arcsec", "arcsec"),
        ("_mjd", "d"),
        ("_deg", "deg"),
        ("_au", "AU"),
        ("_px", "pix"),
    )
    for suffix, unit in suffixes:
        if column.endswith(suffix):
            return unit
    if column == "pred_v_mag" or column.startswith("mag_") or column.startswith("magerr_"):
        return "mag"
    return None


def validate_parquet_units(book: CheckBook, frozen: Path) -> None:
    core_products = (
        "known_all.parquet",
        "known_matched.parquet",
        "known_mask15.parquet",
        "unknown_catalog.parquet",
        "orbit_links.parquet",
        "orbit_obs_residuals.parquet",
    )
    present = [frozen / name for name in core_products if (frozen / name).exists()]
    missing_products = [name for name in core_products if not (frozen / name).exists()]
    if not present:
        book.add(
            "units.parquet_metadata",
            "units",
            "SKIP",
            "Core Parquet physical columns carry unit metadata.",
            expected="unit metadata on all unit-bearing core columns",
            observed=None,
            sources=[frozen],
            detail="No core Parquet products are available.",
        )
        return
    if pq is None:
        book.add(
            "units.parquet_metadata",
            "units",
            "SKIP",
            "Core Parquet physical columns carry unit metadata.",
            expected="unit metadata on all unit-bearing core columns",
            observed=None,
            sources=present,
            detail="pyarrow is unavailable.",
        )
        return
    wrong: list[dict[str, Any]] = []
    absent: list[dict[str, Any]] = []
    checked = 0
    try:
        for path in present:
            schema = pq.read_schema(path)
            for item in schema:
                expected = inferred_unit(item.name)
                if expected is None:
                    continue
                checked += 1
                metadata = item.metadata or {}
                raw_observed = metadata.get(b"unit")
                if raw_observed is None:
                    absent.append({"file": path.name, "column": item.name, "expected": expected})
                    continue
                observed = raw_observed.decode("utf-8", errors="replace")
                if observed != expected:
                    wrong.append(
                        {
                            "file": path.name,
                            "column": item.name,
                            "expected": expected,
                            "observed": observed,
                        }
                    )
    except Exception as exc:
        book.add(
            "units.parquet_metadata",
            "units",
            "FAIL",
            "Core Parquet physical columns carry unit metadata.",
            expected="readable schemas with consistent units",
            observed=type(exc).__name__,
            sources=present,
            detail=str(exc),
        )
        return
    status = "FAIL" if wrong else ("SKIP" if absent or missing_products else "PASS")
    book.add(
        "units.parquet_metadata",
        "units",
        status,
        "Core Parquet physical columns carry unit metadata.",
        expected="unit metadata on all unit-bearing core columns",
        observed={
            "columns_checked": checked,
            "wrong": wrong,
            "missing_metadata": absent,
            "missing_products": missing_products,
        },
        sources=present,
        detail="Contradictory units fail; absent metadata or not-yet-produced files remain incomplete."
        if status != "PASS"
        else "",
    )

    review_path = frozen / "unknown_review_detections.parquet"
    if not review_path.exists() or pq is None:
        book.add(
            "units.unknown_membership_metadata",
            "units",
            "SKIP",
            "Unknown membership physical columns carry unit metadata.",
            expected="explicit units for time, coordinates, photometry, and motion",
            observed=None,
            sources=[review_path],
            detail="Artifact or Parquet schema reader is unavailable.",
        )
    else:
        schema = pq.read_schema(review_path)
        missing = []
        wrong = []
        for item in schema:
            expected = inferred_unit(item.name)
            if expected is None:
                continue
            raw = (item.metadata or {}).get(b"unit")
            if raw is None:
                missing.append({"column": item.name, "expected": expected})
            elif raw.decode("utf-8", errors="replace") != expected:
                wrong.append(
                    {
                        "column": item.name,
                        "expected": expected,
                        "observed": raw.decode("utf-8", errors="replace"),
                    }
                )
        book.add(
            "units.unknown_membership_metadata",
            "units",
            "FAIL" if wrong else ("SKIP" if missing else "PASS"),
            "Unknown membership physical columns carry unit metadata.",
            expected="explicit units for time, coordinates, photometry, and motion",
            observed={"missing_metadata": missing, "wrong": wrong},
            sources=[review_path],
            detail="Missing metadata is an unfinished definition, not a numerical mismatch."
            if missing and not wrong
            else "",
        )


def validate_frozen_products(
    book: CheckBook,
    frozen: Path,
    computed: dict[str, Any],
) -> None:
    products = (
        (
            "known_all.parquet",
            "known.known_prediction_rows",
            "Frozen known-object prediction rows",
            "known_prediction_rows",
        ),
        (
            "known_matched.parquet",
            "known.known_matched_1arcsec_rows",
            "Frozen 1 arcsec known-object matched rows",
            "known_matched_1arcsec_rows",
        ),
        (
            "known_mask15.parquet",
            "known.known_mask15_rows",
            "Frozen 1.5 arcsec known-mask rows",
            "known_mask15_rows",
        ),
        (
            "unknown_catalog.parquet",
            "unknown.catalog_linkages",
            "Formal unknown-catalog linkages",
            "unknown_catalog_linkages",
        ),
    )
    for filename, check_id, description, key in products:
        path = frozen / filename
        expected = EXPECTED_COUNTS[key]
        if not path.exists():
            book.skipped_metric(
                check_id,
                "frozen products",
                description,
                expected,
                path,
                "The Parquet product has not been staged.",
            )
            continue
        if pq is None:
            book.skipped_metric(
                check_id,
                "frozen products",
                description,
                expected,
                path,
                "pyarrow is unavailable; the expected value was not substituted.",
            )
            continue
        try:
            observed = parquet_row_count(path)
        except Exception as exc:
            book.add(
                f"{check_id}.readability",
                "frozen products",
                "FAIL",
                f"{description} Parquet metadata is readable.",
                expected="valid Parquet metadata",
                observed=type(exc).__name__,
                sources=[path],
                detail=str(exc),
            )
            book.skipped_metric(
                check_id,
                "frozen products",
                description,
                expected,
                path,
                "The present Parquet product could not be read.",
            )
            continue
        book.metric(check_id, "frozen products", description, expected, observed, path)
        computed[key] = observed

    catalog_path = frozen / "unknown_catalog.parquet"
    if not catalog_path.exists():
        book.add(
            "unknown.formal_quality_flags",
            "unknown population",
            "SKIP",
            "Every formal unknown-catalog row is both fit_ok and is_good.",
            expected={"rows": 4_762, "fit_ok_true": 4_762, "is_good_true": 4_762},
            observed=None,
            sources=[catalog_path],
            detail="The formal catalog has not been staged.",
        )
    elif pq is None or pd is None:
        book.add(
            "unknown.formal_quality_flags",
            "unknown population",
            "SKIP",
            "Every formal unknown-catalog row is both fit_ok and is_good.",
            expected={"rows": 4_762, "fit_ok_true": 4_762, "is_good_true": 4_762},
            observed=None,
            sources=[catalog_path],
            detail="pandas or pyarrow is unavailable.",
        )
    else:
        try:
            flags = pq.read_table(catalog_path, columns=["fit_ok", "is_good"]).to_pandas()
            fit = parse_bool_series(flags["fit_ok"])
            good = parse_bool_series(flags["is_good"])
            observed = {
                "rows": int(len(flags)),
                "fit_ok_true": int(fit.sum()),
                "is_good_true": int(good.sum()),
                "fit_ok_false_or_null": int((~fit).sum()),
                "is_good_false_or_null": int((~good).sum()),
            }
            expected = {
                "rows": EXPECTED_COUNTS["unknown_catalog_linkages"],
                "fit_ok_true": EXPECTED_COUNTS["unknown_catalog_linkages"],
                "is_good_true": EXPECTED_COUNTS["unknown_catalog_linkages"],
                "fit_ok_false_or_null": 0,
                "is_good_false_or_null": 0,
            }
            book.add(
                "unknown.formal_quality_flags",
                "unknown population",
                "PASS" if observed == expected else "FAIL",
                "Every formal unknown-catalog row is both fit_ok and is_good.",
                expected=expected,
                observed=observed,
                sources=[catalog_path],
                detail="Flags are read from the formal Parquet rows, not from the derived summary JSON.",
            )
        except Exception as exc:
            book.add(
                "unknown.formal_quality_flags",
                "unknown population",
                "FAIL",
                "Every formal unknown-catalog row is both fit_ok and is_good.",
                expected={"fit_ok_false_or_null": 0, "is_good_false_or_null": 0},
                observed=type(exc).__name__,
                sources=[catalog_path],
                detail=str(exc),
            )

    status_path = frozen / "file_status.csv"
    zero_metrics = (
        (
            "unknown.nonempty_nights",
            "Nights with a present, nonempty formal unknown catalog",
            EXPECTED_COUNTS["unknown_nonempty_nights"],
        ),
        (
            "unknown.true_zero_nights",
            "Nights with a present formal unknown catalog containing zero rows",
            EXPECTED_COUNTS["unknown_true_zero_nights"],
        ),
    )
    if not status_path.exists():
        add_unavailable_metrics(
            book,
            zero_metrics,
            "unknown population",
            status_path,
            "The file-status ledger has not been staged.",
        )
    elif pd is None:
        add_unavailable_metrics(
            book,
            zero_metrics,
            "unknown population",
            status_path,
            "pandas is unavailable.",
        )
    else:
        try:
            status = load_csv(status_path, dtype={"night": "string"})
            required = {"product", "exists", "status", "rows_written", "night"}
            missing = required.difference(status.columns)
            if missing:
                raise ValueError(f"missing columns: {sorted(missing)}")
            unknown = status.loc[status["product"].eq("unknown_catalog")].copy()
            present_ok = parse_bool_series(unknown["exists"]) & (
                unknown["status"].astype("string").str.startswith("ok")
            )
            rows = pd.to_numeric(unknown["rows_written"], errors="coerce")
            nonempty = int((present_ok & rows.gt(0)).sum())
            true_zero = int((present_ok & rows.eq(0)).sum())
            for (check_id, description, expected), observed in zip(
                zero_metrics, (nonempty, true_zero)
            ):
                book.metric(
                    check_id,
                    "unknown population",
                    description,
                    expected,
                    observed,
                    status_path,
                    detail="Computed from exists/status/rows_written, so a missing file is not a true zero.",
                )
            computed["unknown_nonempty_nights"] = nonempty
            computed["unknown_true_zero_nights"] = true_zero
            book.add(
                "unknown.catalog_night_accounting",
                "unknown population",
                "PASS" if nonempty + true_zero == 124 else "FAIL",
                "Present formal unknown catalogs partition into nonempty and true-zero nights.",
                expected={"present_catalog_nights": 124, "nonempty_plus_zero": "116 + 8"},
                observed={
                    "present_catalog_nights": int(present_ok.sum()),
                    "nonempty_nights": nonempty,
                    "true_zero_nights": true_zero,
                },
                sources=[status_path],
            )
        except Exception as exc:
            book.add(
                "unknown.file_status_readability",
                "unknown population",
                "FAIL",
                "The formal unknown file-status ledger can distinguish missing and true-zero nights.",
                expected="product, exists, status, rows_written, night",
                observed=type(exc).__name__,
                sources=[status_path],
                detail=str(exc),
            )
            add_unavailable_metrics(
                book,
                zero_metrics,
                "unknown population",
                status_path,
                "Counts were not computed because the present ledger is invalid.",
            )
    validate_parquet_units(book, frozen)


def validate_unknown_memberships(
    book: CheckBook,
    derived_unknown: Path,
    computed: dict[str, Any],
) -> tuple[Any | None, set[str] | None]:
    path = derived_unknown / "unknown_all_link_memberships.csv"
    metric_specs = (
        (
            "unknown.membership_rows",
            "Unknown linkage-detection membership rows",
            EXPECTED_COUNTS["unknown_memberships"],
        ),
        (
            "unknown.globally_unique_detections",
            "Globally unique (night, image, objID) unknown detections",
            EXPECTED_COUNTS["unknown_unique_detections"],
        ),
    )
    if not path.exists():
        add_unavailable_metrics(
            book,
            metric_specs,
            "unknown population",
            path,
            "The derived membership table has not been produced.",
        )
        return None, None
    if pd is None:
        add_unavailable_metrics(
            book,
            metric_specs,
            "unknown population",
            path,
            "pandas is unavailable.",
        )
        return None, None
    try:
        columns = set(load_csv(path, nrows=0).columns)
        wanted = {
            "night",
            "image_name",
            "obj_id",
            "trk_sub",
            "linkage_id",
            "final_paper_status",
            "lin_speed_arcsec_per_day",
        }
        memberships = load_csv(
            path,
            usecols=sorted(columns.intersection(wanted)),
            dtype={
                "night": "string",
                "image_name": "string",
                "obj_id": "string",
                "trk_sub": "string",
                "linkage_id": "string",
                "final_paper_status": "string",
            },
        )
        detection_keys, missing_key_rows = normalized_detection_keys(memberships)
        observed_memberships = int(len(memberships))
        observed_unique = int(len(detection_keys))
        book.metric(
            metric_specs[0][0],
            "unknown population",
            metric_specs[0][1],
            metric_specs[0][2],
            observed_memberships,
            path,
            detail="Every table row is one linkage-detection membership.",
        )
        book.metric(
            metric_specs[1][0],
            "unknown population",
            metric_specs[1][1],
            metric_specs[1][2],
            observed_unique,
            path,
            detail="Keys are independently reconstructed from night, image_name, and obj_id.",
        )
        computed["unknown_memberships"] = observed_memberships
        computed["unknown_unique_detections"] = observed_unique
        book.add(
            "unknown.detection_key_completeness",
            "unknown population",
            "PASS" if missing_key_rows == 0 else "FAIL",
            "Every unknown membership has a complete detection identity.",
            expected=0,
            observed=missing_key_rows,
            sources=[path],
        )
        book.add(
            "unknown.membership_unique_distinction",
            "unknown population",
            "PASS" if observed_memberships >= observed_unique else "FAIL",
            "Membership and globally unique detection denominators remain distinct.",
            expected={
                "memberships": EXPECTED_COUNTS["unknown_memberships"],
                "unique": EXPECTED_COUNTS["unknown_unique_detections"],
                "difference": 140,
            },
            observed={
                "memberships": observed_memberships,
                "unique": observed_unique,
                "difference": observed_memberships - observed_unique,
            },
            sources=[path],
        )
    except Exception as exc:
        book.add(
            "unknown.membership_readability",
            "unknown population",
            "FAIL",
            "The derived membership table is readable at membership and detection grain.",
            expected="night, image_name, obj_id",
            observed=type(exc).__name__,
            sources=[path],
            detail=str(exc),
        )
        add_unavailable_metrics(
            book,
            metric_specs,
            "unknown population",
            path,
            "Metrics were not computed because the present table is invalid.",
        )
        return None, None

    unique_path = derived_unknown / "unknown_unique_detections.csv"
    if not unique_path.exists():
        book.add(
            "unknown.unique_table_reconciliation",
            "unknown population",
            "SKIP",
            "The frozen unique-detection table equals the membership-derived key set.",
            expected=len(detection_keys),
            observed=None,
            sources=[path, unique_path],
            detail="The unique-detection derivative has not been produced.",
        )
    else:
        try:
            columns = set(load_csv(unique_path, nrows=0).columns)
            unique_frame = load_csv(
                unique_path,
                usecols=sorted(columns.intersection({"night", "image_name", "obj_id"})),
                dtype={"night": "string", "image_name": "string", "obj_id": "string"},
            )
            frozen_keys, missing_rows = normalized_detection_keys(unique_frame)
            match = frozen_keys == detection_keys and missing_rows == 0 and len(unique_frame) == len(frozen_keys)
            book.add(
                "unknown.unique_table_reconciliation",
                "unknown population",
                "PASS" if match else "FAIL",
                "The frozen unique-detection table equals the membership-derived key set.",
                expected={"rows": len(detection_keys), "keys": len(detection_keys)},
                observed={
                    "rows": int(len(unique_frame)),
                    "keys": len(frozen_keys),
                    "missing_key_rows": missing_rows,
                    "symmetric_difference": len(frozen_keys.symmetric_difference(detection_keys)),
                },
                sources=[path, unique_path],
            )
        except Exception as exc:
            book.add(
                "unknown.unique_table_reconciliation",
                "unknown population",
                "FAIL",
                "The frozen unique-detection table equals the membership-derived key set.",
                expected=len(detection_keys),
                observed=type(exc).__name__,
                sources=[path, unique_path],
                detail=str(exc),
            )
    return memberships, detection_keys


def read_link_key_source(
    path: Path,
    *,
    night_column: str = "night",
    retained_only: bool = False,
) -> tuple[set[str], int, int]:
    columns = set(load_csv(path, nrows=0).columns)
    wanted = {night_column, "trk_sub", "linkage_id", "final_paper_status"}
    frame = load_csv(
        path,
        usecols=sorted(columns.intersection(wanted)),
        dtype={
            night_column: "string",
            "trk_sub": "string",
            "linkage_id": "string",
            "final_paper_status": "string",
        },
    )
    if retained_only:
        if "final_paper_status" not in frame:
            raise ValueError("missing final_paper_status")
        frame = frame.loc[frame["final_paper_status"].eq(RETAINED_STATUS)]
    keys, missing_rows = normalized_link_keys(frame, night_column=night_column)
    return keys, int(len(frame)), missing_rows


def validate_review_sample(
    book: CheckBook,
    review: Path,
    derived_unknown: Path,
    frozen: Path,
    memberships: Any | None,
    computed: dict[str, Any],
) -> None:
    status_path = review / "review_and_mpc_status.csv"
    metrics = (
        ("review.initial_linkages", "Initially human-selected linkages", "initial_review_linkages"),
        ("review.retained_linkages", "Post-audit retained linkages", "retained_review_linkages"),
        ("review.rejected_linkages", "Post-audit rejected linkages", "rejected_review_linkages"),
    )
    retained_keys: set[str] | None = None
    if not status_path.exists():
        for check_id, description, key in metrics:
            book.skipped_metric(
                check_id,
                "review sample",
                description,
                EXPECTED_COUNTS[key],
                status_path,
                "The review-status ledger has not been staged.",
            )
    elif pd is None:
        for check_id, description, key in metrics:
            book.skipped_metric(
                check_id,
                "review sample",
                description,
                EXPECTED_COUNTS[key],
                status_path,
                "pandas is unavailable.",
            )
    else:
        try:
            status = load_csv(
                status_path,
                dtype={"origin_night": "string", "trk_sub": "string", "linkage_id": "string"},
            )
            required = {"origin_night", "trk_sub", "linkage_id", "final_paper_status"}
            missing = required.difference(status.columns)
            if missing:
                raise ValueError(f"missing columns: {sorted(missing)}")
            retained = status.loc[status["final_paper_status"].eq(RETAINED_STATUS)]
            rejected = status.loc[status["final_paper_status"].eq(REJECTED_STATUS)]
            initial = pd.concat([retained, rejected], ignore_index=True)
            observed = {
                "initial_review_linkages": int(len(initial)),
                "retained_review_linkages": int(len(retained)),
                "rejected_review_linkages": int(len(rejected)),
            }
            for check_id, description, key in metrics:
                book.metric(
                    check_id,
                    "review sample",
                    description,
                    EXPECTED_COUNTS[key],
                    observed[key],
                    status_path,
                )
                computed[key] = observed[key]
            retained_keys, missing_key_rows = normalized_link_keys(
                retained, night_column="origin_night"
            )
            unknown_statuses = sorted(
                set(status["final_paper_status"].dropna().astype(str))
                .difference({RETAINED_STATUS, REJECTED_STATUS})
            )
            exact_partition = (
                observed["initial_review_linkages"]
                == observed["retained_review_linkages"] + observed["rejected_review_linkages"]
                and len(status) == observed["initial_review_linkages"]
                and missing_key_rows == 0
                and not unknown_statuses
            )
            book.add(
                "review.audit_partition",
                "review sample",
                "PASS" if exact_partition else "FAIL",
                "The 67-link initial sample partitions exactly into 58 retained and 9 rejected links.",
                expected={"initial": 67, "retained": 58, "rejected": 9},
                observed={
                    "ledger_rows": int(len(status)),
                    "initial": observed["initial_review_linkages"],
                    "retained": observed["retained_review_linkages"],
                    "rejected": observed["rejected_review_linkages"],
                    "retained_unique_keys": len(retained_keys),
                    "missing_key_rows": missing_key_rows,
                    "unrecognized_statuses": unknown_statuses,
                },
                sources=[status_path],
            )
        except Exception as exc:
            book.add(
                "review.status_readability",
                "review sample",
                "FAIL",
                "The review-status ledger is readable at linkage grain.",
                expected="origin_night, trk_sub, linkage_id, final_paper_status",
                observed=type(exc).__name__,
                sources=[status_path],
                detail=str(exc),
            )

    retained_member_keys: set[str] | None = None
    if memberships is None or "final_paper_status" not in memberships:
        book.skipped_metric(
            "review.retained_detection_memberships",
            "review sample",
            "Detection memberships in the retained 58-link sample",
            EXPECTED_COUNTS["retained_detection_memberships"],
            derived_unknown / "unknown_all_link_memberships.csv",
            "The membership table or its final status column is unavailable.",
        )
    else:
        retained_memberships = memberships.loc[
            memberships["final_paper_status"].eq(RETAINED_STATUS)
        ].copy()
        observed_rows = int(len(retained_memberships))
        book.metric(
            "review.retained_detection_memberships",
            "review sample",
            "Detection memberships in the retained 58-link sample",
            EXPECTED_COUNTS["retained_detection_memberships"],
            observed_rows,
            derived_unknown / "unknown_all_link_memberships.csv",
            detail="Rows are filtered by final_paper_status from the full membership table.",
        )
        computed["retained_detection_memberships"] = observed_rows
        try:
            retained_detection_keys, missing_detection_rows = normalized_detection_keys(
                retained_memberships
            )
            retained_member_keys, missing_link_rows = normalized_link_keys(retained_memberships)
            observed = {
                "membership_rows": observed_rows,
                "unique_detection_keys": len(retained_detection_keys),
                "unique_link_keys": len(retained_member_keys),
                "missing_detection_key_rows": missing_detection_rows,
                "missing_link_key_rows": missing_link_rows,
            }
            expected = {
                "membership_rows": 179,
                "unique_detection_keys": 179,
                "unique_link_keys": 58,
                "missing_detection_key_rows": 0,
                "missing_link_key_rows": 0,
            }
            book.add(
                "review.retained_detection_identity",
                "review sample",
                "PASS" if observed == expected else "FAIL",
                "The retained sample has 179 complete, globally unique detection keys across 58 links.",
                expected=expected,
                observed=observed,
                sources=[derived_unknown / "unknown_all_link_memberships.csv"],
            )
        except Exception as exc:
            book.add(
                "review.retained_detection_identity",
                "review sample",
                "FAIL",
                "The retained sample has 179 complete, globally unique detection keys across 58 links.",
                expected={"unique_detection_keys": 179, "unique_link_keys": 58},
                observed=type(exc).__name__,
                sources=[derived_unknown / "unknown_all_link_memberships.csv"],
                detail=str(exc),
            )

    expected_key_sources = {
        "review_status_retained": status_path,
        "review_links": review / "unknown_high_confidence_links.csv",
        "review_detections": review / "unknown_high_confidence_detections.csv",
        "derived_links": derived_unknown / "unknown_high_confidence_links_recomputed.csv",
        "derived_detections": derived_unknown
        / "unknown_high_confidence_detections_recomputed.csv",
    }
    key_sets: dict[str, set[str]] = {}
    source_rows: dict[str, int] = {}
    missing_sources: list[str] = []
    errors: list[str] = []
    if retained_keys is not None:
        key_sets["review_status_retained"] = retained_keys
        source_rows["review_status_retained"] = len(retained_keys)
    else:
        missing_sources.append("review_status_retained")
    for name, path in expected_key_sources.items():
        if name == "review_status_retained":
            continue
        if not path.exists():
            missing_sources.append(name)
            continue
        try:
            keys, rows, missing_rows = read_link_key_source(path)
            key_sets[name] = keys
            source_rows[name] = rows
            if missing_rows:
                errors.append(f"{name} has {missing_rows} incomplete key rows")
        except Exception as exc:
            errors.append(f"{name}: {type(exc).__name__}: {exc}")
    if retained_member_keys is not None:
        key_sets["full_membership_retained"] = retained_member_keys
        source_rows["full_membership_retained"] = int(
            memberships["final_paper_status"].eq(RETAINED_STATUS).sum()
        )
    else:
        missing_sources.append("full_membership_retained")

    reference = retained_keys or (next(iter(key_sets.values())) if key_sets else set())
    symmetric_differences = {
        name: len(keys.symmetric_difference(reference)) for name, keys in key_sets.items()
    }
    wrong_sizes = {name: len(keys) for name, keys in key_sets.items() if len(keys) != 58}

    catalog_path = frozen / "unknown_catalog.parquet"
    catalog_missing_retained: int | None = None
    if reference and catalog_path.exists() and pq is not None and pd is not None:
        try:
            catalog = pq.read_table(
                catalog_path, columns=["night", "trk_sub", "linkage_id"]
            ).to_pandas()
            catalog_keys, catalog_missing_rows = normalized_link_keys(catalog)
            catalog_missing_retained = len(reference.difference(catalog_keys))
            if catalog_missing_rows:
                errors.append(
                    f"formal catalog has {catalog_missing_rows} incomplete linkage keys"
                )
        except Exception as exc:
            errors.append(f"formal catalog key read: {type(exc).__name__}: {exc}")
    elif reference:
        missing_sources.append("formal_unknown_catalog")

    mismatch = (
        bool(errors)
        or bool(wrong_sizes)
        or any(symmetric_differences.values())
        or catalog_missing_retained not in (None, 0)
    )
    status = "FAIL" if mismatch else ("SKIP" if missing_sources or not key_sets else "PASS")
    book.add(
        "review.retained_key_consistency",
        "review sample",
        status,
        "The same 58 retained linkage keys propagate through review and derived products.",
        expected={"unique_link_keys_per_source": 58, "symmetric_difference": 0},
        observed={
            "source_rows": source_rows,
            "source_unique_keys": {name: len(keys) for name, keys in key_sets.items()},
            "symmetric_differences": symmetric_differences,
            "wrong_sizes": wrong_sizes,
            "formal_catalog_missing_retained": catalog_missing_retained,
            "missing_sources": sorted(set(missing_sources)),
            "errors": errors,
        },
        sources=list(expected_key_sources.values())
        + [derived_unknown / "unknown_all_link_memberships.csv", catalog_path],
        detail="A missing derivative is SKIP; a present source with a different key set is FAIL.",
    )


def validate_speed_units(book: CheckBook, derived_unknown: Path) -> None:
    path = derived_unknown / "unknown_all_links.csv"
    if not path.exists():
        book.add(
            "units.unknown_speed_conversion",
            "units",
            "SKIP",
            "Unknown-link arcsec/day values convert to arcsec/hour by division by 24.",
            expected="speed_arcsec_per_hour = lin_speed_arcsec_per_day / 24",
            observed=None,
            sources=[path],
            detail="The derived link table has not been produced.",
        )
        return
    if pd is None:
        book.add(
            "units.unknown_speed_conversion",
            "units",
            "SKIP",
            "Unknown-link arcsec/day values convert to arcsec/hour by division by 24.",
            expected="speed_arcsec_per_hour = lin_speed_arcsec_per_day / 24",
            observed=None,
            sources=[path],
            detail="pandas is unavailable.",
        )
        return
    try:
        frame = load_csv(
            path,
            usecols=["lin_speed_arcsec_per_day", "speed_arcsec_per_hour"],
        )
        per_day = pd.to_numeric(frame["lin_speed_arcsec_per_day"], errors="coerce")
        per_hour = pd.to_numeric(frame["speed_arcsec_per_hour"], errors="coerce")
        valid = per_day.notna() & per_hour.notna()
        error = (per_hour.loc[valid] - per_day.loc[valid] / 24.0).abs()
        max_error = float(error.max()) if len(error) else None
        bad_rows = int((error > 1e-10).sum())
        book.add(
            "units.unknown_speed_conversion",
            "units",
            "PASS" if bad_rows == 0 and int(valid.sum()) > 0 else "FAIL",
            "Unknown-link arcsec/day values convert to arcsec/hour by division by 24.",
            expected={"bad_rows": 0, "formula": "per_day / 24"},
            observed={
                "rows_checked": int(valid.sum()),
                "bad_rows": bad_rows,
                "max_absolute_error_arcsec_per_hour": max_error,
            },
            sources=[path],
        )
    except Exception as exc:
        book.add(
            "units.unknown_speed_conversion",
            "units",
            "FAIL",
            "Unknown-link arcsec/day values convert to arcsec/hour by division by 24.",
            expected="speed_arcsec_per_hour = lin_speed_arcsec_per_day / 24",
            observed=type(exc).__name__,
            sources=[path],
            detail=str(exc),
        )


HASH_LINE = re.compile(r"^([0-9a-fA-F]{64})\s+[ *]?(.+?)\s*$")


def validate_hash_manifest(
    book: CheckBook,
    check_id: str,
    description: str,
    manifest: Path,
) -> None:
    if not manifest.exists():
        book.add(
            check_id,
            "provenance",
            "SKIP",
            description,
            expected="hashes.sha256 with resolvable entries",
            observed=None,
            sources=[manifest],
            detail="The hash manifest has not been produced.",
        )
        return
    malformed: list[int] = []
    missing: list[str] = []
    mismatched: list[str] = []
    verified: list[str] = []
    try:
        for line_number, raw_line in enumerate(
            manifest.read_text(encoding="utf-8").splitlines(), start=1
        ):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            match = HASH_LINE.match(line)
            if not match:
                malformed.append(line_number)
                continue
            expected, relative = match.groups()
            path = manifest.parent / relative
            if not path.exists():
                missing.append(relative)
                continue
            if sha256(path).lower() != expected.lower():
                mismatched.append(relative)
            else:
                verified.append(relative)
    except Exception as exc:
        book.add(
            check_id,
            "provenance",
            "FAIL",
            description,
            expected="readable and verifiable hash manifest",
            observed=type(exc).__name__,
            sources=[manifest],
            detail=str(exc),
        )
        return
    status = "FAIL" if malformed or mismatched else ("SKIP" if missing or not verified else "PASS")
    book.add(
        check_id,
        "provenance",
        status,
        description,
        expected={"malformed": 0, "missing": 0, "mismatched": 0},
        observed={
            "verified": len(verified),
            "missing": missing,
            "mismatched": mismatched,
            "malformed_lines": malformed,
        },
        sources=[manifest],
        detail="Missing staged files are incomplete; malformed lines or digest mismatches fail."
        if status != "PASS"
        else "",
    )


def validate_frozen_hashes(book: CheckBook, frozen: Path) -> None:
    manifest = frozen / "row_counts.json"
    if not manifest.exists():
        book.add(
            "provenance.frozen_product_hashes",
            "provenance",
            "SKIP",
            "Frozen Parquet files match hashes and sizes recorded by extraction.",
            expected="row_counts.json with product digests",
            observed=None,
            sources=[manifest],
            detail="The frozen-product extraction manifest has not been produced.",
        )
        return
    try:
        content = json.loads(manifest.read_text(encoding="utf-8"))
        products = content.get("products", {})
        if not isinstance(products, dict) or not products:
            raise ValueError("products mapping is empty")
        verified = []
        missing = []
        mismatched = []
        unrecorded = []
        for product, record in sorted(products.items()):
            path = frozen / f"{product}.parquet"
            expected_hash = record.get("parquet_sha256") if isinstance(record, dict) else None
            expected_size = record.get("parquet_size_bytes") if isinstance(record, dict) else None
            if not expected_hash or expected_size is None:
                unrecorded.append(product)
                continue
            if not path.exists():
                missing.append(path.name)
                continue
            actual_hash = sha256(path)
            actual_size = path.stat().st_size
            if actual_hash != expected_hash or actual_size != int(expected_size):
                mismatched.append(
                    {
                        "file": path.name,
                        "expected_hash": expected_hash,
                        "actual_hash": actual_hash,
                        "expected_size": int(expected_size),
                        "actual_size": actual_size,
                    }
                )
            else:
                verified.append(path.name)
        status = "FAIL" if mismatched else ("SKIP" if missing or unrecorded else "PASS")
        book.add(
            "provenance.frozen_product_hashes",
            "provenance",
            status,
            "Frozen Parquet files match hashes and sizes recorded by extraction.",
            expected={"mismatched": 0, "missing": 0, "unrecorded": 0},
            observed={
                "verified": verified,
                "missing": missing,
                "unrecorded": unrecorded,
                "mismatched": mismatched,
            },
            sources=[manifest] + [frozen / name for name in verified],
        )
    except Exception as exc:
        book.add(
            "provenance.frozen_product_hashes",
            "provenance",
            "FAIL",
            "Frozen Parquet files match hashes and sizes recorded by extraction.",
            expected="valid row_counts.json product digests",
            observed=type(exc).__name__,
            sources=[manifest],
            detail=str(exc),
        )


def validate_summary_input_hashes(
    book: CheckBook,
    check_id: str,
    summary_path: Path,
    search_dirs: Iterable[Path],
    description: str,
) -> None:
    if not summary_path.exists():
        book.add(
            check_id,
            "provenance",
            "SKIP",
            description,
            expected="summary JSON with input_hashes",
            observed=None,
            sources=[summary_path],
            detail="The derived summary has not been produced.",
        )
        return
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        hashes = summary.get("input_hashes")
        if not isinstance(hashes, dict) or not hashes:
            book.add(
                check_id,
                "provenance",
                "SKIP",
                description,
                expected="nonempty input_hashes mapping",
                observed=hashes,
                sources=[summary_path],
                detail="Input digests are not yet recorded.",
            )
            return
        verified = []
        missing = []
        mismatched = []
        for filename, expected in hashes.items():
            candidates = [directory / filename for directory in search_dirs]
            path = next((candidate for candidate in candidates if candidate.exists()), None)
            if path is None:
                missing.append(filename)
            elif sha256(path) != expected:
                mismatched.append(filename)
            else:
                verified.append(filename)
        status = "FAIL" if mismatched else ("SKIP" if missing else "PASS")
        book.add(
            check_id,
            "provenance",
            status,
            description,
            expected={"missing": 0, "mismatched": 0},
            observed={"verified": verified, "missing": missing, "mismatched": mismatched},
            sources=[summary_path],
        )
    except Exception as exc:
        book.add(
            check_id,
            "provenance",
            "FAIL",
            description,
            expected="readable summary and matching input digests",
            observed=type(exc).__name__,
            sources=[summary_path],
            detail=str(exc),
        )


def validate_provenance(
    book: CheckBook,
    inventory: Path,
    frozen: Path,
    review: Path,
    derived_unknown: Path,
    derived_known: Path | None,
    tables: Path | None,
) -> None:
    manifest_expectations = {
        "provenance.raw_manifest_presence": inventory / "raw_manifest.csv",
        "provenance.engineering_manifest_presence": inventory
        / "raw_engineering_manifest.csv",
        "provenance.l2_manifest_presence": inventory / "l2_manifest.csv",
        "provenance.frozen_file_ledger_presence": frozen / "file_status.csv",
        "provenance.review_ledger_presence": review / "review_and_mpc_status.csv",
    }
    for check_id, path in manifest_expectations.items():
        book.add(
            check_id,
            "provenance",
            "PASS" if path.exists() and path.stat().st_size > 0 else "SKIP",
            f"Manifest/ledger {path.name} exists and is nonempty.",
            expected="present and nonempty",
            observed={"exists": path.exists(), "size_bytes": path.stat().st_size if path.exists() else None},
            sources=[path],
            detail="Not-yet-collected manifests are incomplete, not zero-valued evidence."
            if not path.exists()
            else "",
        )

    validate_frozen_hashes(book, frozen)
    validate_hash_manifest(
        book,
        "provenance.review_hash_manifest",
        "Review-sample files match their standalone SHA-256 manifest.",
        review / "hashes.sha256",
    )
    provenance_dir = frozen.parent / "provenance"
    validate_hash_manifest(
        book,
        "provenance.production_hash_manifest",
        "Production-code, Gaia, environment, and footprint files match their SHA-256 manifest.",
        provenance_dir / "hashes.sha256",
    )
    validate_hash_manifest(
        book,
        "provenance.inventory_hash_manifest",
        "Inventory files match their standalone SHA-256 manifest.",
        inventory / "hashes.sha256",
    )
    if tables is not None:
        validate_hash_manifest(
            book,
            "provenance.table_hash_manifest",
            "Frozen table outputs match their standalone SHA-256 manifest.",
            tables / "hashes.sha256",
        )
    validate_summary_input_hashes(
        book,
        "provenance.derived_unknown_inputs",
        derived_unknown / "unknown_population_summary.json",
        [frozen, review, inventory],
        "The unknown-population summary input digests match frozen inputs.",
    )
    if derived_known is None:
        book.add(
            "provenance.derived_known_inputs",
            "provenance",
            "SKIP",
            "The known-population summary input digests match frozen inputs.",
            expected="known_population_summary.json with input_hashes",
            observed=None,
            detail="No derived-known directory was supplied.",
        )
    else:
        validate_summary_input_hashes(
            book,
            "provenance.derived_known_inputs",
            derived_known / "known_population_summary.json",
            [frozen, inventory],
            "The known-population summary input digests match frozen inputs.",
        )


def validate_derived_known(
    book: CheckBook,
    derived_known: Path | None,
    computed: dict[str, Any],
) -> None:
    if derived_known is None:
        book.add(
            "known.derived_summary_reconciliation",
            "derived known",
            "SKIP",
            "The derived-known summary reconciles to independently counted frozen rows.",
            expected="known_population_summary.json",
            observed=None,
            detail="No derived-known directory was supplied.",
        )
        return
    path = derived_known / "known_population_summary.json"
    if not path.exists():
        book.add(
            "known.derived_summary_reconciliation",
            "derived known",
            "SKIP",
            "The derived-known summary reconciles to independently counted frozen rows.",
            expected="known_population_summary.json",
            observed=None,
            sources=[path],
            detail="The derived-known analysis is not ready.",
        )
        return
    try:
        content = json.loads(path.read_text(encoding="utf-8"))
        pairs = {
            "source_prediction_rows_including_duplicates": computed.get(
                "known_prediction_rows"
            ),
            "matched_1arcsec_n": computed.get("known_matched_1arcsec_rows"),
        }
        unavailable = [key for key, value in pairs.items() if value is None]
        missing_keys = [key for key in pairs if key not in content]
        mismatches = {
            key: {"summary": content.get(key), "independent": value}
            for key, value in pairs.items()
            if value is not None and key in content and content.get(key) != value
        }
        status = "FAIL" if mismatches else (
            "SKIP" if unavailable or missing_keys else "PASS"
        )
        book.add(
            "known.derived_summary_reconciliation",
            "derived known",
            status,
            "The derived-known summary reconciles to independently counted frozen rows.",
            expected=pairs,
            observed={
                "summary_values": {key: content.get(key) for key in pairs},
                "unavailable_independent_values": unavailable,
                "missing_summary_keys": missing_keys,
                "mismatches": mismatches,
            },
            sources=[path],
        )
    except Exception as exc:
        book.add(
            "known.derived_summary_reconciliation",
            "derived known",
            "FAIL",
            "The derived-known summary reconciles to independently counted frozen rows.",
            expected="readable JSON summary",
            observed=type(exc).__name__,
            sources=[path],
            detail=str(exc),
        )


def validate_quality_mask_table(
    book: CheckBook,
    tables: Path | None,
    config: dict[str, Any] | None,
) -> None:
    if tables is None:
        path = Path("quality_mask.csv")
    else:
        path = tables / "quality_mask.csv"
    if tables is None or not path.exists():
        book.add(
            "quality.table_reconciliation",
            "quality mask",
            "SKIP",
            "The per-night quality-mask table exactly implements snapshot.json.",
            expected="quality_mask.csv",
            observed=None,
            sources=[path],
            detail="The frozen night tables have not been produced or supplied.",
        )
        return
    if config is None or pd is None:
        book.add(
            "quality.table_reconciliation",
            "quality mask",
            "SKIP",
            "The per-night quality-mask table exactly implements snapshot.json.",
            expected="readable table and configuration",
            observed=None,
            sources=[path],
            detail="Configuration or pandas is unavailable.",
        )
        return
    try:
        frame = load_csv(path, dtype={"night": "string"})
        required = {
            "night",
            "quality_code",
            "primary_science_included",
            "unknown_science_included",
        }
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"missing columns: {sorted(missing)}")
        frame["night"] = frame["night"].str.zfill(8)
        duplicate_nights = int(frame["night"].duplicated().sum())
        start = pd.Timestamp(config["observation_start"])
        end = pd.Timestamp(config["observation_end"])
        expected_nights = set(pd.date_range(start, end, freq="D").strftime("%Y%m%d"))
        observed_nights = set(frame["night"].dropna())
        rules = config["quality_rules"]
        excluded = set(rules["exclude_from_primary_science"])
        quarantine = set(rules["quarantine_unknown"])
        raw_without_l2 = set(rules["raw_without_l2"])
        primary = parse_bool_series(frame["primary_science_included"])
        unknown = parse_bool_series(frame["unknown_science_included"])
        expected_primary = ~frame["night"].isin(excluded)
        expected_unknown = expected_primary & ~frame["night"].isin(quarantine)
        wrong_primary = int((primary != expected_primary).sum())
        wrong_unknown = int((unknown != expected_unknown).sum())
        code = frame["quality_code"].fillna("").astype(str)
        missing_excluded_codes = int(
            (frame["night"].isin(excluded) & ~code.str.contains("excluded_primary")).sum()
        )
        missing_quarantine_codes = int(
            (frame["night"].isin(quarantine) & ~code.str.contains("quarantine_unknown")).sum()
        )
        missing_l2_codes = int(
            (frame["night"].isin(raw_without_l2) & ~code.str.contains("raw_without_l2")).sum()
        )
        absent_nights = sorted(expected_nights.difference(observed_nights))
        extra_nights = sorted(observed_nights.difference(expected_nights))
        contradictions = {
            "duplicate_nights": duplicate_nights,
            "wrong_primary_flags": wrong_primary,
            "wrong_unknown_flags": wrong_unknown,
            "missing_excluded_codes": missing_excluded_codes,
            "missing_quarantine_codes": missing_quarantine_codes,
            "missing_raw_without_l2_codes": missing_l2_codes,
            "extra_nights": extra_nights,
        }
        has_contradiction = any(
            value for key, value in contradictions.items() if key != "extra_nights"
        ) or bool(extra_nights)
        status = "FAIL" if has_contradiction else ("SKIP" if absent_nights else "PASS")
        book.add(
            "quality.table_reconciliation",
            "quality mask",
            status,
            "The per-night quality-mask table exactly implements snapshot.json.",
            expected={"night_rows": len(expected_nights), "contradictions": 0},
            observed={
                "night_rows": int(len(frame)),
                "absent_nights": absent_nights,
                **contradictions,
            },
            sources=[path],
            detail="Missing future rows are incomplete; contradictory existing flags fail.",
        )
    except Exception as exc:
        book.add(
            "quality.table_reconciliation",
            "quality mask",
            "FAIL",
            "The per-night quality-mask table exactly implements snapshot.json.",
            expected="readable quality_mask.csv",
            observed=type(exc).__name__,
            sources=[path],
            detail=str(exc),
        )


def validate_stage_units(book: CheckBook, tables: Path | None) -> None:
    path = (tables / "unknown_stage_definitions.csv") if tables is not None else Path(
        "unknown_stage_definitions.csv"
    )
    if tables is None or not path.exists():
        book.add(
            "units.stage_definitions",
            "units",
            "SKIP",
            "Every unknown-funnel stage has an explicit detection/tracklet/linkage unit.",
            expected=EXPECTED_STAGE_UNITS,
            observed=None,
            sources=[path],
            detail="The stage-definition table has not been produced or supplied.",
        )
        return
    if pd is None:
        book.add(
            "units.stage_definitions",
            "units",
            "SKIP",
            "Every unknown-funnel stage has an explicit detection/tracklet/linkage unit.",
            expected=EXPECTED_STAGE_UNITS,
            observed=None,
            sources=[path],
            detail="pandas is unavailable.",
        )
        return
    try:
        frame = load_csv(path, dtype="string")
        required = {"column", "unit"}
        missing_columns = required.difference(frame.columns)
        if missing_columns:
            raise ValueError(f"missing columns: {sorted(missing_columns)}")
        observed = {
            str(row.column): str(row.unit)
            for row in frame[["column", "unit"]].itertuples(index=False)
        }
        missing = sorted(set(EXPECTED_STAGE_UNITS).difference(observed))
        wrong = {
            key: {"expected": expected, "observed": observed.get(key)}
            for key, expected in EXPECTED_STAGE_UNITS.items()
            if key in observed and observed[key] != expected
        }
        status = "FAIL" if wrong else ("SKIP" if missing else "PASS")
        book.add(
            "units.stage_definitions",
            "units",
            status,
            "Every unknown-funnel stage has an explicit detection/tracklet/linkage unit.",
            expected=EXPECTED_STAGE_UNITS,
            observed={"definitions": observed, "missing": missing, "wrong": wrong},
            sources=[path],
            detail="Missing definitions are unfinished; contradictory units fail.",
        )
    except Exception as exc:
        book.add(
            "units.stage_definitions",
            "units",
            "FAIL",
            "Every unknown-funnel stage has an explicit detection/tracklet/linkage unit.",
            expected=EXPECTED_STAGE_UNITS,
            observed=type(exc).__name__,
            sources=[path],
            detail=str(exc),
        )


def validate_snapshot_summary(
    book: CheckBook,
    tables: Path | None,
    computed: dict[str, Any],
) -> None:
    path = (tables / "snapshot_summary.json") if tables is not None else Path(
        "snapshot_summary.json"
    )
    if tables is None or not path.exists():
        book.add(
            "tables.snapshot_summary_reconciliation",
            "tables",
            "SKIP",
            "The headline snapshot summary reconciles to independent computations.",
            expected="snapshot_summary.json",
            observed=None,
            sources=[path],
            detail="The frozen table summary has not been produced or supplied.",
        )
        return
    key_map = {
        "strict_raw_frames": "strict_raw_exposures",
        "strict_raw_nights": "strict_raw_nights",
        "raw_fields": "strict_raw_fields",
        "all_mp_fits": "engineering_raw_exposures",
        "l2_catalogs": "strict_l2_catalogs",
        "known_predictions": "known_prediction_rows",
        "known_matches_1arcsec": "known_matched_1arcsec_rows",
        "known_matches_1p5arcsec": "known_mask15_rows",
        "unknown_catalog_linkages": "unknown_catalog_linkages",
        "initial_review_selected_linkages": "initial_review_linkages",
        "posthoc_retained_linkages": "retained_review_linkages",
        "posthoc_retained_detection_memberships": "retained_detection_memberships",
    }
    try:
        summary = json.loads(path.read_text(encoding="utf-8"))
        missing_summary = []
        unavailable = []
        mismatches = {}
        compared = {}
        for summary_key, computed_key in key_map.items():
            if summary_key not in summary:
                missing_summary.append(summary_key)
            elif computed_key not in computed:
                unavailable.append(computed_key)
            else:
                compared[summary_key] = {
                    "summary": summary[summary_key],
                    "independent": computed[computed_key],
                }
                if summary[summary_key] != computed[computed_key]:
                    mismatches[summary_key] = compared[summary_key]
        status = "FAIL" if mismatches else (
            "SKIP" if missing_summary or unavailable else "PASS"
        )
        book.add(
            "tables.snapshot_summary_reconciliation",
            "tables",
            status,
            "The headline snapshot summary reconciles to independent computations.",
            expected="all available headline values equal independent counts",
            observed={
                "compared": compared,
                "mismatches": mismatches,
                "missing_summary_keys": missing_summary,
                "unavailable_independent_values": unavailable,
            },
            sources=[path],
        )
    except Exception as exc:
        book.add(
            "tables.snapshot_summary_reconciliation",
            "tables",
            "FAIL",
            "The headline snapshot summary reconciles to independent computations.",
            expected="readable JSON summary",
            observed=type(exc).__name__,
            sources=[path],
            detail=str(exc),
        )


def inspect_png(path: Path) -> dict[str, Any]:
    size = path.stat().st_size
    with path.open("rb") as handle:
        header = handle.read(24)
    if size <= 24 or len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("invalid or empty PNG signature/header")
    width, height = struct.unpack(">II", header[16:24])
    if width <= 0 or height <= 0:
        raise ValueError("non-positive PNG dimensions")
    return {"size_bytes": size, "width_px": width, "height_px": height}


def inspect_pdf(path: Path) -> dict[str, Any]:
    size = path.stat().st_size
    if size <= 8:
        raise ValueError("empty PDF")
    try:
        from pypdf import PdfReader
    except ImportError as exc:  # pragma: no cover - depends on runtime
        raise RuntimeError("pypdf is unavailable for PDF dimension validation") from exc
    reader = PdfReader(path)
    if not reader.pages:
        raise ValueError("PDF has no pages")
    dimensions = []
    for page in reader.pages:
        width = float(page.mediabox.width)
        height = float(page.mediabox.height)
        if width <= 0 or height <= 0:
            raise ValueError("non-positive PDF MediaBox dimensions")
        dimensions.append({"width_pt": width, "height_pt": height})
    return {"size_bytes": size, "pages": len(reader.pages), "page_dimensions": dimensions}


def validate_figures(book: CheckBook, figures: Path | None) -> None:
    ready = 0
    for index, stem in enumerate(FIGURE_STEMS, start=1):
        pdf = (figures / f"{stem}.pdf") if figures is not None else Path(f"{stem}.pdf")
        png = (figures / f"{stem}.png") if figures is not None else Path(f"{stem}.png")
        present = {"pdf": pdf.exists(), "png": png.exists()}
        if not any(present.values()):
            book.add(
                f"figures.fig{index:02d}",
                "figures",
                "SKIP",
                f"Figure {index} has nonempty PDF and PNG files with positive dimensions.",
                expected={"pdf": True, "png": True},
                observed=present,
                sources=[pdf, png],
                detail="This figure has not been produced.",
            )
            continue
        metadata: dict[str, Any] = {}
        errors: list[str] = []
        for kind, path, inspector in (("pdf", pdf, inspect_pdf), ("png", png, inspect_png)):
            if not path.exists():
                continue
            try:
                metadata[kind] = inspector(path)
            except RuntimeError as exc:
                metadata[kind] = {"inspection_unavailable": str(exc)}
            except Exception as exc:
                errors.append(f"{kind}: {type(exc).__name__}: {exc}")
        inspection_unavailable = any(
            isinstance(item, dict) and "inspection_unavailable" in item
            for item in metadata.values()
        )
        if errors:
            status = "FAIL"
        elif not all(present.values()) or inspection_unavailable:
            status = "SKIP"
        else:
            status = "PASS"
            ready += 1
        book.add(
            f"figures.fig{index:02d}",
            "figures",
            status,
            f"Figure {index} has nonempty PDF and PNG files with positive dimensions.",
            expected={"pdf": True, "png": True, "positive_dimensions": True},
            observed={"present": present, "metadata": metadata, "errors": errors},
            sources=[pdf, png],
            detail="A missing format or unavailable dimension backend is incomplete; a corrupt present file fails."
            if status != "PASS"
            else "",
        )
    book.add(
        "figures.complete_set",
        "figures",
        "PASS" if ready == len(FIGURE_STEMS) else "SKIP",
        "All 12 planned figures are available in validated PDF and PNG form.",
        expected=len(FIGURE_STEMS),
        observed=ready,
        sources=[figures] if figures is not None else [],
        detail=f"{len(FIGURE_STEMS) - ready} figure(s) remain incomplete."
        if ready != len(FIGURE_STEMS)
        else "",
    )


def validate(
    *,
    config_path: Path,
    inventory: Path,
    frozen: Path,
    review: Path,
    derived_unknown: Path,
    derived_known: Path | None,
    tables: Path | None,
    figures: Path | None,
) -> tuple[CheckBook, dict[str, Any]]:
    book = CheckBook()
    computed: dict[str, Any] = {}
    config = validate_config(book, config_path)
    validate_inventory(book, inventory, computed)
    validate_frozen_products(book, frozen, computed)
    memberships, _ = validate_unknown_memberships(book, derived_unknown, computed)
    validate_review_sample(
        book, review, derived_unknown, frozen, memberships, computed
    )
    validate_speed_units(book, derived_unknown)
    validate_derived_known(book, derived_known, computed)
    validate_quality_mask_table(book, tables, config)
    validate_stage_units(book, tables)
    validate_snapshot_summary(book, tables, computed)
    validate_provenance(
        book,
        inventory,
        frozen,
        review,
        derived_unknown,
        derived_known,
        tables,
    )
    validate_figures(book, figures)
    return book, computed


def compact_markdown(value: Any, limit: int = 320) -> str:
    if value is None:
        return "—"
    if isinstance(value, str):
        text = value
    else:
        text = json.dumps(json_value(value), sort_keys=True, ensure_ascii=False)
    text = text.replace("|", "\\|").replace("\n", " ")
    return text if len(text) <= limit else text[: limit - 1] + "…"


def build_report(payload: dict[str, Any]) -> str:
    counts = payload["status_counts"]
    lines = [
        "# Snapshot validation report",
        "",
        f"- Overall status: **{payload['overall_status']}**",
        f"- Generated (UTC): `{payload['generated_utc']}`",
        f"- Snapshot label: `{payload.get('snapshot_label') or 'unknown'}`",
        f"- Checks: {counts.get('PASS', 0)} PASS, {counts.get('FAIL', 0)} FAIL, {counts.get('SKIP', 0)} SKIP",
        "",
        "`SKIP` means the source artifact, author approval, or inspection backend is not yet available. It never means an observed value of zero. Expected headline values are comparison targets; observed values are independently computed.",
        "",
        "## Headline closure",
        "",
        "| Check | Status | Expected | Observed |",
        "|---|---:|---:|---:|",
    ]
    headline_ids = {
        "inventory.strict_raw_exposures",
        "inventory.strict_raw_nights",
        "inventory.strict_raw_fields",
        "inventory.engineering_raw_exposures",
        "inventory.engineering_raw_nights",
        "inventory.strict_l2_catalogs",
        "inventory.strict_l2_nights",
        "known.known_prediction_rows",
        "known.known_matched_1arcsec_rows",
        "known.known_mask15_rows",
        "unknown.catalog_linkages",
        "unknown.membership_rows",
        "unknown.globally_unique_detections",
        "unknown.nonempty_nights",
        "unknown.true_zero_nights",
        "review.initial_linkages",
        "review.retained_linkages",
        "review.rejected_linkages",
        "review.retained_detection_memberships",
    }
    for check in payload["checks"]:
        if check["id"] in headline_ids:
            lines.append(
                f"| `{check['id']}` | **{check['status']}** | {compact_markdown(check['expected'])} | {compact_markdown(check['observed'])} |"
            )

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for check in payload["checks"]:
        grouped[check["category"]].append(check)
    lines.extend(["", "## Detailed checks", ""])
    for category in sorted(grouped):
        lines.extend(
            [
                f"### {category.title()}",
                "",
                "| Check | Status | Expected | Observed | Note |",
                "|---|---:|---|---|---|",
            ]
        )
        for check in grouped[category]:
            detail = check["detail"]
            if check["sources"]:
                source_text = ", ".join(check["sources"])
                detail = f"{detail} Sources: {source_text}".strip()
            lines.append(
                f"| `{check['id']}` | **{check['status']}** | {compact_markdown(check['expected'])} | {compact_markdown(check['observed'])} | {compact_markdown(detail)} |"
            )
        lines.append("")

    failures = [check for check in payload["checks"] if check["status"] == "FAIL"]
    skips = [check for check in payload["checks"] if check["status"] == "SKIP"]
    lines.extend(["## Blocking failures", ""])
    if failures:
        for check in failures:
            lines.append(f"- `{check['id']}` — {check['description']}")
    else:
        lines.append("- None.")
    lines.extend(["", "## Incomplete / awaiting input", ""])
    if skips:
        for check in skips:
            lines.append(f"- `{check['id']}` — {check['description']}")
    else:
        lines.append("- None.")
    lines.append("")
    return "\n".join(lines)


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    project = Path(__file__).resolve().parents[1]
    final_inventory = project / "snapshot" / "inventory"
    partial_inventory = project / "snapshot" / "inventory_partial"
    default_inventory = final_inventory if final_inventory.exists() else partial_inventory
    parser = argparse.ArgumentParser(
        description="Recompute and validate P0/P1 closure for the frozen PASP snapshot."
    )
    parser.add_argument("--config", type=Path, default=project / "config" / "snapshot.json")
    parser.add_argument("--inventory", type=Path, default=default_inventory)
    parser.add_argument(
        "--frozen-products",
        type=Path,
        default=project / "snapshot" / "frozen_products",
    )
    parser.add_argument(
        "--review-sample",
        type=Path,
        default=project / "snapshot" / "review_sample",
    )
    parser.add_argument(
        "--derived-unknown",
        type=Path,
        default=project / "snapshot" / "derived_unknown",
    )
    parser.add_argument(
        "--derived-known",
        type=Path,
        default=project / "snapshot" / "derived_known",
        help="Optional derived-known directory; a missing directory is reported as SKIP.",
    )
    parser.add_argument(
        "--tables",
        type=Path,
        default=project / "tables",
        help="Optional frozen tables directory; a missing directory is reported as SKIP.",
    )
    parser.add_argument(
        "--figures",
        type=Path,
        default=project / "figures",
        help="Optional figure directory containing all 12 PDF/PNG pairs.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=project / "qa" / "validation_results.json",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=project / "reports" / "VALIDATION_REPORT.md",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = args.config.resolve()
    inventory = args.inventory.resolve()
    frozen = args.frozen_products.resolve()
    review = args.review_sample.resolve()
    derived_unknown = args.derived_unknown.resolve()
    derived_known = args.derived_known.resolve() if args.derived_known else None
    tables = args.tables.resolve() if args.tables else None
    figures = args.figures.resolve() if args.figures else None

    book, computed = validate(
        config_path=config_path,
        inventory=inventory,
        frozen=frozen,
        review=review,
        derived_unknown=derived_unknown,
        derived_known=derived_known,
        tables=tables,
        figures=figures,
    )
    status_counts = dict(Counter(check["status"] for check in book.checks))
    if status_counts.get("FAIL", 0):
        overall = "FAIL"
    elif status_counts.get("SKIP", 0):
        overall = "INCOMPLETE"
    else:
        overall = "PASS"
    snapshot_label = None
    try:
        snapshot_label = json.loads(config_path.read_text(encoding="utf-8")).get(
            "snapshot_label"
        )
    except Exception:
        pass
    payload = {
        "schema_version": "1.0",
        "generated_utc": utc_now(),
        "snapshot_label": snapshot_label,
        "overall_status": overall,
        "status_counts": {
            status: int(status_counts.get(status, 0))
            for status in ("PASS", "FAIL", "SKIP")
        },
        "inputs": {
            "config": str(config_path),
            "inventory": str(inventory),
            "frozen_products": str(frozen),
            "review_sample": str(review),
            "derived_unknown": str(derived_unknown),
            "derived_known": str(derived_known) if derived_known else None,
            "tables": str(tables) if tables else None,
            "figures": str(figures) if figures else None,
        },
        "acceptance_targets": EXPECTED_COUNTS,
        "independently_computed": json_value(computed),
        "checks": book.checks,
    }
    atomic_write_text(
        args.output_json.resolve(),
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
    )
    atomic_write_text(args.report.resolve(), build_report(payload))
    print(
        json.dumps(
            {
                "overall_status": overall,
                "status_counts": payload["status_counts"],
                "output_json": str(args.output_json.resolve()),
                "report": str(args.report.resolve()),
            },
            indent=2,
            sort_keys=True,
        )
    )
    if overall == "FAIL":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
