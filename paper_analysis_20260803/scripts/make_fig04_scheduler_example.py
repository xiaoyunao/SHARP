#!/usr/bin/env python3
"""Render Fig. 4 from frozen scheduler realization products.

The figure is deliberately evidence-limited.  It visualizes archived planned
and actual timestamps and only those Sun, Moon, altitude, or airmass columns
that are present in a frozen input.  It never reconstructs astronomical
conditions from coordinates, a site configuration, or a current ephemeris.

The historical near-Sun block is labelled as the RA-prefix implementation
that actually ran.  A planned exposure without a matched strict raw frame is
not attributed to weather, equipment, or observatory efficiency.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from figure_styles import (
    STATISTICS_COLORS,
    STATISTICS_FOUR_PANEL_FIGSIZE,
    apply_statistics_style,
    style_statistics_axis,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCHEDULER_DIR = PROJECT_ROOT / "snapshot" / "scheduler"
DEFAULT_EXAMPLE = DEFAULT_SCHEDULER_DIR / "scheduler_example_20260420.csv"
DEFAULT_NIGHT_SUMMARY = (
    DEFAULT_SCHEDULER_DIR / "scheduler_plan_realization_by_night.csv"
)
DEFAULT_OUTPUT_BASE = PROJECT_ROOT / "figures" / "fig04_scheduler_example"
DEFAULT_FIGURE_DATA_DIR = PROJECT_ROOT / "figure_data"
DEFAULT_NIGHT = "20260420"
SCHEMA_VERSION = "1.0"

NORMAL_COLOR = STATISTICS_COLORS["blue"]
NEAR_SUN_COLOR = STATISTICS_COLORS["orange"]
INK = STATISTICS_COLORS["ink"]
MID_GREY = STATISTICS_COLORS["mid_grey"]
LIGHT_GREY = STATISTICS_COLORS["light_grey"]

MODE_COLORS = {
    "normal": NORMAL_COLOR,
    "near_sun": NEAR_SUN_COLOR,
    "unplanned": MID_GREY,
    "unknown": MID_GREY,
}
MODE_MARKERS = {"normal": "o", "near_sun": "^", "unplanned": "D", "unknown": "s"}

EXAMPLE_REQUIRED = {
    "night",
    "field_id",
    "pair_id",
    "match_status",
    "mode",
    "cycle_id",
    "repeat",
    "planned_start_utc",
    "actual_start_utc",
}
SUMMARY_REQUIRED = {
    "night",
    "plan_exposure_n",
    "acquired_strict_frame_n",
    "multiset_matched_frame_n",
    "acquired_frame_plan_compliance",
    "planned_frame_realization",
}
VALID_MATCH_STATUS = {"matched", "planned_not_acquired", "acquired_not_planned"}

CONDITION_SPECS: dict[str, dict[str, Any]] = {
    "sun_separation": {
        "label": "Sun separation",
        "unit": "deg",
        "aliases": (
            "sun_sep_deg",
            "solar_elongation_deg",
            "sun_separation_deg",
            "solar_sep_deg",
            "sep_to_sun_deg",
        ),
        "prefer_time": "planned",
    },
    "moon_separation": {
        "label": "Moon separation",
        "unit": "deg",
        "aliases": ("moon_sep_deg", "moon_separation_deg", "sep_to_moon_deg"),
        "prefer_time": "planned",
    },
    "altitude": {
        "label": "Target altitude",
        "unit": "deg",
        "aliases": (
            "altitude_deg",
            "target_altitude_deg",
            "tel_elevation_deg",
            "elevation_deg",
            "alt_deg",
        ),
        "prefer_time": "actual",
    },
    "airmass": {
        "label": "Airmass",
        "unit": "",
        "aliases": ("airmass", "target_airmass", "secz"),
        "prefer_time": "actual",
    },
}

TIME_ALIASES = (
    "time_utc",
    "obs_time_utc",
    "start_utc",
    "planned_start_utc",
    "actual_start_utc",
)

FIGURE_DATA_NAMES = {
    "timeline": "fig04_scheduler_timeline.csv",
    "cadence": "fig04_scheduler_revisit_cadence.csv",
    "residuals": "fig04_scheduler_start_residuals.csv",
    "conditions": "fig04_scheduler_conditions.csv",
    "summary": "fig04_scheduler_night_summary.csv",
    "metadata": "fig04_scheduler_metadata.json",
}

NEAR_SUN_NOTE = (
    "Historical near-Sun rows are the executed RA-prefix implementation: an initial "
    "solar-separation ordering was overwritten by RA sorting before prefix selection."
)
ATTRIBUTION_NOTE = (
    "Planned–unacquired means only that no strict raw frame matched the planned "
    "(night, field) multiplicity; no weather or equipment cause is assigned."
)
ARCHIVE_NOTE = (
    "Missing condition tracks are marked ‘not archived’ and are not reconstructed "
    "from coordinates, site parameters, or a current ephemeris."
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def combined_hash(hashes: dict[str, str]) -> str:
    payload = "".join(f"{key}:{value}\n" for key, value in sorted(hashes.items()))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def normalize_text(value: Any) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace").strip()
    return str(value if value is not None else "").strip()


def normalize_night(value: Any) -> str:
    digits = re.sub(r"\D", "", normalize_text(value))
    return digits[:8] if len(digits) >= 8 else ""


def normalize_field_id(value: Any) -> str:
    text = normalize_text(value)
    if re.fullmatch(r"\d+(?:\.0+)?", text):
        return f"{int(float(text)):04d}"
    return text


def normalize_mode(value: Any) -> str:
    text = normalize_text(value).lower().replace("-", "_").replace(" ", "_")
    return text or "unknown"


def normalize_cycle(value: Any) -> str:
    text = normalize_text(value)
    if re.fullmatch(r"[-+]?\d+(?:\.0+)?", text):
        return str(int(float(text)))
    return text


def truthy(value: Any) -> bool:
    return normalize_text(value).lower() in {"1", "true", "yes", "y", "t"}


def natural_key(value: Any) -> tuple[Any, ...]:
    pieces = re.split(r"(\d+)", normalize_text(value))
    return tuple(int(piece) if piece.isdigit() else piece.lower() for piece in pieces)


def parse_time_series(frame: pd.DataFrame, column: str, source: Path) -> pd.Series:
    raw = frame[column].map(normalize_text)
    # Archived timestamps legitimately mix integral-second ISO strings with
    # microsecond-resolution camera times.  Pandas 2.x otherwise infers one
    # strict format from the first non-empty value and rejects the other form.
    parsed = pd.to_datetime(
        raw.mask(raw.eq("")), utc=True, errors="coerce", format="mixed"
    )
    invalid = raw.ne("") & parsed.isna()
    if invalid.any():
        examples = frame.loc[invalid, column].head().tolist()
        raise ValueError(f"{source}: invalid UTC timestamps in {column}: {examples}")
    return parsed


def numeric_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").replace([np.inf, -np.inf], np.nan)


def iso_utc(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    timestamp = pd.Timestamp(value)
    if timestamp.tzinfo is None:
        timestamp = timestamp.tz_localize("UTC")
    else:
        timestamp = timestamp.tz_convert("UTC")
    return timestamp.isoformat().replace("+00:00", "Z")


def read_example(path: Path, night: str) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False)
    missing = sorted(EXAMPLE_REQUIRED - set(frame.columns))
    if missing:
        raise ValueError(f"{path} lacks required Fig. 4 columns: {missing}")
    frame["night"] = frame["night"].map(normalize_night)
    frame = frame.loc[frame["night"].eq(night)].copy()
    if frame.empty:
        raise ValueError(f"{path} contains no rows for representative night {night}")
    if frame["pair_id"].map(normalize_text).duplicated().any():
        duplicates = frame.loc[
            frame["pair_id"].map(normalize_text).duplicated(False), "pair_id"
        ].head().tolist()
        raise ValueError(f"{path}: pair_id must be unique; duplicates={duplicates}")

    frame["match_status"] = frame["match_status"].map(normalize_text)
    unknown_status = sorted(set(frame["match_status"]) - VALID_MATCH_STATUS)
    if unknown_status:
        raise ValueError(f"{path}: unsupported match_status values: {unknown_status}")
    frame["field_id"] = frame["field_id"].map(normalize_field_id)
    frame["mode"] = frame["mode"].map(normalize_mode)
    frame["cycle_id"] = frame["cycle_id"].map(normalize_cycle)
    frame["repeat_number"] = numeric_series(frame, "repeat")
    planned_mask = frame["match_status"].isin({"matched", "planned_not_acquired"})
    invalid_repeat = planned_mask & (
        frame["repeat_number"].isna()
        | (frame["repeat_number"] < 1)
        | ~np.isclose(frame["repeat_number"], np.round(frame["repeat_number"]))
    )
    if invalid_repeat.any():
        examples = frame.loc[invalid_repeat, ["pair_id", "repeat"]].head().to_dict("records")
        raise ValueError(f"{path}: invalid planned repeat numbers: {examples}")
    frame.loc[planned_mask, "repeat_number"] = np.round(
        frame.loc[planned_mask, "repeat_number"]
    )

    frame["planned_dt"] = parse_time_series(frame, "planned_start_utc", path)
    frame["actual_dt"] = parse_time_series(frame, "actual_start_utc", path)
    if frame.loc[planned_mask, "planned_dt"].isna().any():
        examples = frame.loc[
            planned_mask & frame["planned_dt"].isna(), "pair_id"
        ].head().tolist()
        raise ValueError(f"{path}: planned rows lack planned_start_utc: {examples}")

    matched = frame["match_status"].eq("matched")
    both_timed = matched & frame["planned_dt"].notna() & frame["actual_dt"].notna()
    frame["computed_time_residual_s"] = np.nan
    frame.loc[both_timed, "computed_time_residual_s"] = (
        frame.loc[both_timed, "actual_dt"] - frame.loc[both_timed, "planned_dt"]
    ).dt.total_seconds()
    if "time_residual_s" in frame.columns:
        supplied = numeric_series(frame, "time_residual_s")
        check = both_timed & supplied.notna()
        mismatch = check & (
            (supplied - frame["computed_time_residual_s"]).abs() > 1.0e-3
        )
        if mismatch.any():
            examples = frame.loc[
                mismatch, ["pair_id", "time_residual_s", "computed_time_residual_s"]
            ].head().to_dict("records")
            raise ValueError(f"{path}: residual does not close to archived timestamps: {examples}")

    duplicate_plan = frame.loc[planned_mask].duplicated(
        ["field_id", "cycle_id", "repeat_number"], keep=False
    )
    if duplicate_plan.any():
        examples = frame.loc[planned_mask].loc[
            duplicate_plan, ["field_id", "cycle_id", "repeat"]
        ].head().to_dict("records")
        raise ValueError(f"{path}: duplicate planned field/cycle/repeat rows: {examples}")
    return frame.reset_index(drop=True)


def read_night_summary(path: Path, night: str, example: pd.DataFrame) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False)
    missing = sorted(SUMMARY_REQUIRED - set(frame.columns))
    if missing:
        raise ValueError(f"{path} lacks required Fig. 4 columns: {missing}")
    frame["night"] = frame["night"].map(normalize_night)
    selected = frame.loc[frame["night"].eq(night)].copy()
    if len(selected) != 1:
        raise ValueError(f"{path}: expected exactly one {night} row, found {len(selected)}")

    expected = {
        "plan_exposure_n": int(
            example["match_status"].isin({"matched", "planned_not_acquired"}).sum()
        ),
        "acquired_strict_frame_n": int(
            example["match_status"].isin({"matched", "acquired_not_planned"}).sum()
        ),
        "multiset_matched_frame_n": int(example["match_status"].eq("matched").sum()),
    }
    row = selected.iloc[0]
    for column, expected_value in expected.items():
        numeric = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
        if not np.isfinite(numeric) or int(round(numeric)) != expected_value:
            raise ValueError(
                f"{path}: {column}={row[column]!r}, but exposure accounting gives "
                f"{expected_value}"
            )
    acquired_ratio = (
        expected["multiset_matched_frame_n"] / expected["acquired_strict_frame_n"]
        if expected["acquired_strict_frame_n"]
        else np.nan
    )
    planned_ratio = (
        expected["multiset_matched_frame_n"] / expected["plan_exposure_n"]
        if expected["plan_exposure_n"]
        else np.nan
    )
    for ratio_column, expected_ratio in (
        ("acquired_frame_plan_compliance", acquired_ratio),
        ("planned_frame_realization", planned_ratio),
    ):
        supplied = pd.to_numeric(
            pd.Series([row[ratio_column]]), errors="coerce"
        ).iloc[0]
        if np.isfinite(expected_ratio) and (
            not np.isfinite(supplied)
            or not math.isclose(supplied, expected_ratio, abs_tol=1e-10)
        ):
            raise ValueError(
                f"{path}: {ratio_column}={row[ratio_column]!r}, "
                f"expected {expected_ratio:.12g}"
            )
    return selected.reset_index(drop=True)


def read_optional_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".fits", ".fit", ".fz"}:
        try:
            from astropy.table import Table
        except ImportError as exc:
            raise RuntimeError(f"astropy is required to read optional FITS table {path}") from exc
        return Table.read(path).to_pandas()
    return pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False)


def normalize_optional_table(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    if "night" in result.columns:
        result["_night"] = result["night"].map(normalize_night)
    if "field_id" in result.columns:
        result["_field_id"] = result["field_id"].map(normalize_field_id)
    if "sequence_id" in result.columns:
        result["_sequence_id"] = result["sequence_id"].map(normalize_field_id)
    return result


def first_numeric_alias(frame: pd.DataFrame, aliases: Sequence[str]) -> str | None:
    for alias in aliases:
        if alias in frame.columns and numeric_series(frame, alias).notna().any():
            return alias
    return None


def direct_condition_points(
    example: pd.DataFrame,
    condition: str,
    source_path: Path,
) -> list[dict[str, Any]]:
    spec = CONDITION_SPECS[condition]
    alias = first_numeric_alias(example, spec["aliases"])
    if alias is None:
        return []
    values = numeric_series(example, alias)
    preferred = example["planned_dt"] if spec["prefer_time"] == "planned" else example["actual_dt"]
    fallback = example["actual_dt"] if spec["prefer_time"] == "planned" else example["planned_dt"]
    timestamps = preferred.where(preferred.notna(), fallback)
    points: list[dict[str, Any]] = []
    for index in example.index[values.notna() & timestamps.notna()]:
        row = example.loc[index]
        points.append(
            condition_point(
                condition=condition,
                timestamp=timestamps.loc[index],
                value=values.loc[index],
                row=row,
                source_kind="scheduler_example",
                source_column=alias,
                source_path=source_path,
            )
        )
    return points


def l2_row_mapping(example: pd.DataFrame, l2: pd.DataFrame) -> dict[int, int]:
    if l2.empty:
        return {}
    working = normalize_optional_table(l2)
    if "strict_standard_catalog" in working.columns:
        strict = working["strict_standard_catalog"].map(truthy)
        working = working.loc[strict].copy()

    path_lookup: dict[str, list[int]] = defaultdict(list)
    name_lookup: dict[str, list[int]] = defaultdict(list)
    key_lookup: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    for index, row in working.iterrows():
        path_value = normalize_text(row.get("path", ""))
        file_value = normalize_text(row.get("file_name", ""))
        if path_value:
            path_lookup[path_value].append(index)
        if file_value:
            name_lookup[file_value].append(index)
        if {"_night", "_field_id", "_sequence_id"}.issubset(working.columns):
            key_lookup[
                (row["_night"], row["_field_id"], row["_sequence_id"])
            ].append(index)

    mapping: dict[int, int] = {}
    for example_index, row in example.iterrows():
        candidates: list[int] = []
        l2_path = normalize_text(row.get("l2_path", ""))
        l2_name = normalize_text(row.get("l2_file_name", ""))
        if l2_path:
            candidates = path_lookup.get(l2_path, [])
        if not candidates and l2_name:
            candidates = name_lookup.get(l2_name, [])
        if not candidates and key_lookup:
            key = (
                normalize_night(row.get("night", "")),
                normalize_field_id(row.get("field_id", "")),
                normalize_field_id(row.get("raw_sequence_id", "")),
            )
            candidates = key_lookup.get(key, [])
        if candidates:
            mapping[example_index] = sorted(candidates)[0]
    l2.attrs["_normalized_working"] = working
    return mapping


def l2_condition_points(
    example: pd.DataFrame,
    l2: pd.DataFrame | None,
    mapping: dict[int, int],
    condition: str,
    source_path: Path | None,
) -> list[dict[str, Any]]:
    if l2 is None or source_path is None or not mapping:
        return []
    working = l2.attrs.get("_normalized_working", l2)
    spec = CONDITION_SPECS[condition]
    alias = first_numeric_alias(working.loc[list(set(mapping.values()))], spec["aliases"])
    if alias is None:
        return []
    points: list[dict[str, Any]] = []
    for example_index, l2_index in mapping.items():
        value = pd.to_numeric(pd.Series([working.at[l2_index, alias]]), errors="coerce").iloc[0]
        if not np.isfinite(value):
            continue
        row = example.loc[example_index]
        timestamp = row["actual_dt"]
        if pd.isna(timestamp):
            for time_alias in ("obs_time_utc", "time_utc"):
                if time_alias in working.columns:
                    candidate = pd.to_datetime(
                        normalize_text(working.at[l2_index, time_alias]),
                        utc=True,
                        errors="coerce",
                    )
                    if not pd.isna(candidate):
                        timestamp = candidate
                        break
        if pd.isna(timestamp):
            continue
        points.append(
            condition_point(
                condition=condition,
                timestamp=timestamp,
                value=value,
                row=row,
                source_kind="l2_manifest",
                source_column=alias,
                source_path=source_path,
            )
        )
    return points


def footprint_condition_points(
    example: pd.DataFrame,
    footprint: pd.DataFrame | None,
    condition: str,
    source_path: Path | None,
    night: str,
) -> list[dict[str, Any]]:
    if footprint is None or source_path is None or footprint.empty:
        return []
    working = normalize_optional_table(footprint)
    spec = CONDITION_SPECS[condition]
    alias = first_numeric_alias(working, spec["aliases"])
    time_alias = next((name for name in TIME_ALIASES if name in working.columns), None)
    if alias is None or time_alias is None:
        return []
    if "_night" in working.columns:
        working = working.loc[working["_night"].eq(night)].copy()
    planned_fields = set(example.loc[example["planned_dt"].notna(), "field_id"])
    if "_field_id" in working.columns:
        working = working.loc[working["_field_id"].isin(planned_fields)].copy()
    if working.empty:
        return []
    timestamps = pd.to_datetime(
        working[time_alias].map(normalize_text).replace("", pd.NA),
        utc=True,
        errors="coerce",
    )
    values = numeric_series(working, alias)
    mode_by_field = (
        example.loc[example["planned_dt"].notna()]
        .groupby("field_id")["mode"]
        .agg(lambda items: items.iloc[0] if items.nunique() == 1 else "unknown")
        .to_dict()
    )
    points: list[dict[str, Any]] = []
    for index in working.index[timestamps.notna() & values.notna()]:
        field_id = (
            normalize_field_id(working.at[index, "field_id"])
            if "field_id" in working.columns
            else ""
        )
        points.append(
            {
                "condition": condition,
                "condition_label": spec["label"],
                "availability": "archived",
                "timestamp_utc": iso_utc(timestamps.loc[index]),
                "value": float(values.loc[index]),
                "unit": spec["unit"],
                "night": night,
                "field_id": field_id,
                "pair_id": "",
                "match_status": "",
                "mode": mode_by_field.get(field_id, "unknown"),
                "source_kind": "footprint_time_series",
                "source_column": alias,
                "source_path": str(source_path),
                "availability_reason": "",
            }
        )
    return points


def condition_point(
    *,
    condition: str,
    timestamp: Any,
    value: float,
    row: pd.Series,
    source_kind: str,
    source_column: str,
    source_path: Path,
) -> dict[str, Any]:
    spec = CONDITION_SPECS[condition]
    return {
        "condition": condition,
        "condition_label": spec["label"],
        "availability": "archived",
        "timestamp_utc": iso_utc(timestamp),
        "value": float(value),
        "unit": spec["unit"],
        "night": normalize_night(row.get("night", "")),
        "field_id": normalize_field_id(row.get("field_id", "")),
        "pair_id": normalize_text(row.get("pair_id", "")),
        "match_status": normalize_text(row.get("match_status", "")),
        "mode": normalize_mode(row.get("mode", "")),
        "source_kind": source_kind,
        "source_column": source_column,
        "source_path": str(source_path),
        "availability_reason": "",
    }


def build_conditions(
    example: pd.DataFrame,
    *,
    l2: pd.DataFrame | None,
    l2_mapping: dict[int, int],
    l2_path: Path | None,
    footprint: pd.DataFrame | None,
    footprint_path: Path | None,
    example_path: Path,
    night: str,
    input_hash: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for condition in CONDITION_SPECS:
        points = direct_condition_points(example, condition, example_path)
        if not points:
            points = l2_condition_points(example, l2, l2_mapping, condition, l2_path)
        if not points:
            points = footprint_condition_points(
                example, footprint, condition, footprint_path, night
            )
        if not points:
            spec = CONDITION_SPECS[condition]
            points = [
                {
                    "condition": condition,
                    "condition_label": spec["label"],
                    "availability": "not_archived",
                    "timestamp_utc": "",
                    "value": np.nan,
                    "unit": spec["unit"],
                    "night": night,
                    "field_id": "",
                    "pair_id": "",
                    "match_status": "",
                    "mode": "",
                    "source_kind": "",
                    "source_column": "",
                    "source_path": "",
                    "availability_reason": (
                        "no archived numeric column with an archived timestamp in "
                        "scheduler example, optional L2 manifest, or optional footprint"
                    ),
                }
            ]
        rows.extend(points)
    result = pd.DataFrame(rows)
    result["input_manifest_hash"] = input_hash
    return result


def build_timeline(example: pd.DataFrame, input_hash: str) -> pd.DataFrame:
    fields = [
        "night",
        "pair_id",
        "field_id",
        "cycle_id",
        "repeat",
        "mode",
        "match_status",
        "match_rule",
        "planned_start_utc",
        "actual_start_utc",
        "actual_time_source",
        "computed_time_residual_s",
        "raw_file_name",
        "raw_sequence_id",
        "l2_file_name",
    ]
    result = pd.DataFrame(index=example.index)
    for field in fields:
        result[field] = example[field] if field in example.columns else ""
    result["planned_start_utc"] = example["planned_dt"].map(iso_utc)
    result["actual_start_utc"] = example["actual_dt"].map(iso_utc)
    result["repeat"] = example["repeat_number"]
    result["input_manifest_hash"] = input_hash
    return result


def build_cadence(example: pd.DataFrame, input_hash: str) -> pd.DataFrame:
    planned = example.loc[
        example["match_status"].isin({"matched", "planned_not_acquired"})
    ].copy()
    rows: list[dict[str, Any]] = []
    for (field_id, cycle_id), group in planned.groupby(
        ["field_id", "cycle_id"], sort=False, dropna=False
    ):
        group = group.sort_values(["repeat_number", "planned_dt", "pair_id"])
        repeat_one = group.loc[group["repeat_number"].eq(1)]
        planned_anchor = (
            repeat_one["planned_dt"].iloc[0]
            if len(repeat_one) == 1
            else group["planned_dt"].iloc[0]
        )
        actual_anchor_rows = repeat_one.loc[
            repeat_one["match_status"].eq("matched") & repeat_one["actual_dt"].notna()
        ]
        actual_anchor = (
            actual_anchor_rows["actual_dt"].iloc[0]
            if len(actual_anchor_rows) == 1
            else pd.NaT
        )
        previous_planned = pd.NaT
        previous_actual = pd.NaT
        for _, item in group.iterrows():
            planned_elapsed = (item["planned_dt"] - planned_anchor).total_seconds() / 60.0
            actual_elapsed = (
                (item["actual_dt"] - actual_anchor).total_seconds() / 60.0
                if not pd.isna(item["actual_dt"]) and not pd.isna(actual_anchor)
                else np.nan
            )
            planned_interval = (
                (item["planned_dt"] - previous_planned).total_seconds() / 60.0
                if not pd.isna(previous_planned)
                else np.nan
            )
            actual_interval = (
                (item["actual_dt"] - previous_actual).total_seconds() / 60.0
                if not pd.isna(item["actual_dt"]) and not pd.isna(previous_actual)
                else np.nan
            )
            rows.append(
                {
                    "night": item["night"],
                    "field_id": field_id,
                    "cycle_id": cycle_id,
                    "repeat": int(item["repeat_number"]),
                    "mode": item["mode"],
                    "pair_id": item["pair_id"],
                    "match_status": item["match_status"],
                    "planned_start_utc": iso_utc(item["planned_dt"]),
                    "actual_start_utc": iso_utc(item["actual_dt"]),
                    "planned_elapsed_from_repeat1_min": planned_elapsed,
                    "actual_elapsed_from_matched_repeat1_min": actual_elapsed,
                    "planned_interval_from_previous_repeat_min": planned_interval,
                    "actual_interval_from_previous_available_repeat_min": actual_interval,
                    "actual_repeat1_anchor_available": not pd.isna(actual_anchor),
                    "input_manifest_hash": input_hash,
                }
            )
            previous_planned = item["planned_dt"]
            if not pd.isna(item["actual_dt"]):
                previous_actual = item["actual_dt"]
    return pd.DataFrame(rows)


def build_residuals(example: pd.DataFrame, input_hash: str) -> pd.DataFrame:
    selected = example.loc[
        example["match_status"].eq("matched")
        & example["planned_dt"].notna()
        & example["actual_dt"].notna()
    ].copy()
    result = pd.DataFrame(
        {
            "night": selected["night"],
            "pair_id": selected["pair_id"],
            "field_id": selected["field_id"],
            "cycle_id": selected["cycle_id"],
            "repeat": selected["repeat_number"],
            "mode": selected["mode"],
            "planned_start_utc": selected["planned_dt"].map(iso_utc),
            "actual_start_utc": selected["actual_dt"].map(iso_utc),
            "time_residual_s": selected["computed_time_residual_s"],
            "abs_time_residual_s": selected["computed_time_residual_s"].abs(),
            "match_rule": selected["match_rule"] if "match_rule" in selected else "",
            "input_manifest_hash": input_hash,
        }
    )
    return result.reset_index(drop=True)


def flatten_json(payload: Any, prefix: str = "") -> dict[str, Any]:
    flattened: dict[str, Any] = {}
    if isinstance(payload, dict):
        for key, value in payload.items():
            child = f"{prefix}.{key}" if prefix else str(key)
            flattened.update(flatten_json(value, child))
    else:
        flattened[prefix.lower()] = payload
    return flattened


def config_thresholds(config: dict[str, Any] | None) -> dict[str, float]:
    if config is None:
        return {}
    flat = flatten_json(config)
    aliases = {
        "altitude_min": (
            "minimum_altitude_deg",
            "min_altitude_deg",
            "altitude_min_deg",
        ),
        "airmass_max": ("maximum_airmass", "max_airmass", "airmass_max"),
    }
    result: dict[str, float] = {}
    for output, names in aliases.items():
        candidates = [
            value
            for key, value in flat.items()
            if any(key == name or key.endswith("." + name) for name in names)
        ]
        numeric = pd.to_numeric(pd.Series(candidates), errors="coerce").dropna()
        if len(numeric) == 1 and np.isfinite(numeric.iloc[0]):
            result[output] = float(numeric.iloc[0])
    return result


def mode_color(mode: Any) -> str:
    return MODE_COLORS.get(normalize_mode(mode), MID_GREY)


def mode_marker(mode: Any) -> str:
    return MODE_MARKERS.get(normalize_mode(mode), "s")


def time_limits(example: pd.DataFrame) -> tuple[pd.Timestamp, pd.Timestamp]:
    values = pd.concat([example["planned_dt"], example["actual_dt"]]).dropna()
    if values.empty:
        raise ValueError("Fig. 4 requires at least one archived planned or actual timestamp")
    padding = pd.Timedelta(minutes=10)
    return values.min() - padding, values.max() + padding


def plot_timeline(
    axis,
    example: pd.DataFrame,
    summary: pd.Series,
    x_limits: tuple[pd.Timestamp, pd.Timestamp],
    night: str,
) -> None:
    planned = example.loc[
        example["match_status"].isin({"matched", "planned_not_acquired"})
    ].copy()
    cycles = sorted(planned["cycle_id"].unique(), key=natural_key)
    actual_only = example["match_status"].eq("acquired_not_planned").any()
    labels = (["unplanned"] if actual_only else []) + cycles
    y_lookup = {label: index for index, label in enumerate(labels)}

    near = planned.loc[planned["mode"].eq("near_sun")]
    if not near.empty:
        axis.axvspan(
            near["planned_dt"].min(),
            near["planned_dt"].max(),
            color=NEAR_SUN_COLOR,
            alpha=0.10,
            linewidth=0,
            zorder=-3,
        )

    for _, row in planned.iterrows():
        y = y_lookup[row["cycle_id"]]
        color = mode_color(row["mode"])
        axis.scatter(
            row["planned_dt"],
            y - 0.10,
            marker="|",
            s=72,
            linewidth=1.35,
            color=color,
            zorder=3,
        )
        if row["match_status"] == "matched" and not pd.isna(row["actual_dt"]):
            axis.plot(
                [row["planned_dt"], row["actual_dt"]],
                [y - 0.10, y + 0.10],
                color=LIGHT_GREY,
                linewidth=0.55,
                alpha=0.65,
                zorder=1,
            )
            axis.scatter(
                row["actual_dt"],
                y + 0.10,
                marker=mode_marker(row["mode"]),
                s=17,
                facecolor="white",
                edgecolor=color,
                linewidth=0.8,
                zorder=4,
            )
        elif row["match_status"] == "planned_not_acquired":
            axis.scatter(
                row["planned_dt"],
                y + 0.10,
                marker="x",
                s=24,
                color=INK,
                linewidth=0.9,
                zorder=5,
            )

    if actual_only:
        for _, row in example.loc[
            example["match_status"].eq("acquired_not_planned") & example["actual_dt"].notna()
        ].iterrows():
            axis.scatter(
                row["actual_dt"],
                y_lookup["unplanned"],
                marker="D",
                s=18,
                facecolor="white",
                edgecolor=MID_GREY,
                linewidth=0.8,
                zorder=4,
            )

    tick_labels = (
        ["Unplanned"] + [f"Cycle {cycle}" for cycle in cycles]
        if actual_only
        else [f"Cycle {cycle}" for cycle in cycles]
    )
    axis.set_yticks(list(range(len(labels))), tick_labels, fontsize=14)
    axis.set_ylim(-0.6, len(labels) - 0.4)
    axis.set_xlim(*x_limits)
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone.utc))
    axis.set_xlabel(f"UTC on {night[:4]}-{night[4:6]}-{night[6:]}", fontsize=18)
    axis.set_ylabel("Scheduler cycle", fontsize=18)
    axis.set_title(
        "(a) Planned sequence and archived starts",
        loc="left",
        fontsize=20,
        fontweight="bold",
    )
    plan_n = int(float(summary["plan_exposure_n"]))
    actual_n = int(float(summary["acquired_strict_frame_n"]))
    matched_n = int(float(summary["multiset_matched_frame_n"]))
    axis.text(
        0.985,
        0.97,
        f"plan {plan_n:,}  |  strict raw {actual_n:,}  |  matched {matched_n:,}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=13,
        color=MID_GREY,
    )
    handles = [
        Line2D(
            [], [], color=NORMAL_COLOR, marker="|", linestyle="none", markersize=11,
            label="Normal plan",
        ),
        Line2D(
            [], [], color=NEAR_SUN_COLOR, marker="|", linestyle="none", markersize=11,
            label="Near-Sun plan",
        ),
        Line2D(
            [], [], color=INK, marker="o", markerfacecolor="white", linestyle="none",
            markersize=5, label="Archived actual start",
        ),
        Line2D(
            [], [], color=INK, marker="x", linestyle="none", markersize=5,
            label="Planned–unacquired",
        ),
    ]
    axis.legend(
        handles=handles,
        fontsize=12,
        loc="upper left",
        ncol=2,
        bbox_to_anchor=(0.0, 0.92),
    )
    style_statistics_axis(axis, tick_fontsize=14)


def quantile_triplet(values: Iterable[float]) -> tuple[float, float, float]:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    if not len(array):
        return np.nan, np.nan, np.nan
    percentiles = np.asarray(np.percentile(array, [16, 50, 84]), dtype=float)
    return float(percentiles[0]), float(percentiles[1]), float(percentiles[2])


def plot_cadence(axis, cadence: pd.DataFrame) -> None:
    repeat_values = sorted(
        int(value) for value in cadence["repeat"].dropna().unique() if int(value) >= 1
    )
    planned_group = cadence.groupby("repeat")["planned_elapsed_from_repeat1_min"]
    planned_x: list[int] = []
    planned_y: list[float] = []
    for repeat in repeat_values:
        values = pd.to_numeric(planned_group.get_group(repeat), errors="coerce").dropna()
        if len(values):
            planned_x.append(repeat)
            planned_y.append(float(values.median()))
    axis.plot(
        planned_x,
        planned_y,
        color=INK,
        linestyle="--",
        linewidth=1.7,
        marker="s",
        markersize=5.5,
        markerfacecolor="white",
        label="Planned median",
        zorder=5,
    )

    offsets = {"normal": -0.08, "near_sun": 0.08}
    for mode in ("normal", "near_sun"):
        mode_rows = cadence.loc[
            cadence["mode"].eq(mode)
            & cadence["actual_repeat1_anchor_available"].astype(bool)
        ]
        xs: list[float] = []
        medians: list[float] = []
        lower: list[float] = []
        upper: list[float] = []
        counts: list[int] = []
        for repeat in repeat_values:
            values = pd.to_numeric(
                mode_rows.loc[
                    mode_rows["repeat"].eq(repeat),
                    "actual_elapsed_from_matched_repeat1_min",
                ],
                errors="coerce",
            ).dropna()
            if values.empty:
                continue
            p16, median, p84 = quantile_triplet(values)
            x = repeat + offsets[mode]
            axis.scatter(
                np.full(len(values), x),
                values,
                s=7,
                alpha=0.13,
                color=mode_color(mode),
                linewidth=0,
                rasterized=True,
                zorder=1,
            )
            xs.append(x)
            medians.append(median)
            lower.append(median - p16)
            upper.append(p84 - median)
            counts.append(len(values))
        if not xs:
            continue
        axis.errorbar(
            xs,
            medians,
            yerr=np.vstack([lower, upper]),
            color=mode_color(mode),
            marker=mode_marker(mode),
            markerfacecolor="white" if mode == "near_sun" else mode_color(mode),
            markeredgecolor=mode_color(mode),
            markersize=6,
            linewidth=1.65,
            capsize=3,
            label=f"{mode.replace('_', '-')} actual p16–p84",
            zorder=6,
        )
        for x, median, count in zip(xs, medians, counts):
            axis.annotate(
                f"n={count}",
                (x, median),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=9.5,
                color=mode_color(mode),
            )

    axis.set_xticks(repeat_values)
    axis.set_xlim(min(repeat_values) - 0.35, max(repeat_values) + 0.35)
    axis.set_xlabel("Visit number within field/cycle", fontsize=18)
    axis.set_ylabel("Minutes since repeat 1", fontsize=18)
    axis.set_title(
        "(b) Three-visit cadence realization",
        loc="left",
        fontsize=20,
        fontweight="bold",
    )
    axis.legend(fontsize=11.5, loc="upper left")
    style_statistics_axis(axis, tick_fontsize=15)


def plot_residuals(
    axis,
    residuals: pd.DataFrame,
    x_limits: tuple[pd.Timestamp, pd.Timestamp],
    night: str,
) -> None:
    axis.axhline(0, color=INK, linewidth=1.1, alpha=0.72, zorder=1)
    if residuals.empty:
        axis.text(
            0.5,
            0.5,
            "Actual start timestamps\nnot archived",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=18,
            color=MID_GREY,
        )
    else:
        residuals = residuals.copy()
        residuals["planned_dt"] = pd.to_datetime(
            residuals["planned_start_utc"], utc=True, errors="coerce"
        )
        residuals["residual_min"] = pd.to_numeric(
            residuals["time_residual_s"], errors="coerce"
        ) / 60.0
        for mode in ("normal", "near_sun"):
            selected = residuals.loc[
                residuals["mode"].eq(mode)
                & residuals["planned_dt"].notna()
                & residuals["residual_min"].notna()
            ]
            if selected.empty:
                continue
            p16, median, p84 = quantile_triplet(selected["residual_min"])
            label = (
                f"{mode.replace('_', '-')}: median {median:+.2f} min "
                f"[p16, p84]=[{p16:+.2f}, {p84:+.2f}]"
            )
            axis.scatter(
                selected["planned_dt"],
                selected["residual_min"],
                s=17 if mode == "normal" else 23,
                marker=mode_marker(mode),
                facecolor=mode_color(mode) if mode == "normal" else "white",
                edgecolor=mode_color(mode),
                linewidth=0.7,
                alpha=0.64,
                rasterized=True,
                label=label,
                zorder=3,
            )
            axis.plot(
                [selected["planned_dt"].min(), selected["planned_dt"].max()],
                [median, median],
                color=mode_color(mode),
                linewidth=1.8,
                linestyle="--",
                zorder=4,
            )
        axis.legend(fontsize=10.5, loc="best")
    axis.set_xlim(*x_limits)
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone.utc))
    axis.set_xlabel(f"Planned UTC on {night[:4]}-{night[4:6]}-{night[6:]}", fontsize=18)
    axis.set_ylabel("Actual − planned start (min)", fontsize=18)
    axis.set_title(
        "(c) Plan–actual start residual",
        loc="left",
        fontsize=20,
        fontweight="bold",
    )
    style_statistics_axis(axis, tick_fontsize=14)


def plot_condition_track(
    axis,
    conditions: pd.DataFrame,
    condition: str,
    x_limits: tuple[pd.Timestamp, pd.Timestamp],
    *,
    threshold: float | None = None,
    threshold_label: str = "",
    show_x: bool,
) -> None:
    spec = CONDITION_SPECS[condition]
    selected = conditions.loc[
        conditions["condition"].eq(condition)
        & conditions["availability"].eq("archived")
    ].copy()
    selected["timestamp"] = pd.to_datetime(
        selected["timestamp_utc"], utc=True, errors="coerce"
    )
    selected["numeric_value"] = pd.to_numeric(selected["value"], errors="coerce")
    selected = selected.dropna(subset=["timestamp", "numeric_value"])
    axis.set_xlim(*x_limits)
    if selected.empty:
        axis.set_facecolor("#fafafa")
        axis.text(
            0.5,
            0.5,
            f"{spec['label']}: not archived",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=13.5,
            color=MID_GREY,
        )
        axis.set_yticks([])
    else:
        for mode in sorted(selected["mode"].map(normalize_mode).unique()):
            mode_rows = selected.loc[selected["mode"].map(normalize_mode).eq(mode)]
            axis.scatter(
                mode_rows["timestamp"],
                mode_rows["numeric_value"],
                s=10,
                marker=mode_marker(mode),
                facecolor=mode_color(mode) if mode != "near_sun" else "white",
                edgecolor=mode_color(mode),
                linewidth=0.55,
                alpha=0.58,
                rasterized=True,
                zorder=3,
            )
        if threshold is not None:
            axis.axhline(threshold, color=INK, linestyle="--", linewidth=1.0, alpha=0.75)
            axis.text(
                0.99,
                0.91,
                threshold_label,
                transform=axis.transAxes,
                ha="right",
                va="top",
                fontsize=9.5,
                color=INK,
            )
        unit = f"({spec['unit']})" if spec["unit"] else ""
        axis.set_ylabel(f"{spec['label']}\n{unit}".rstrip(), fontsize=11)
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone.utc))
    if show_x:
        axis.set_xlabel("UTC", fontsize=12)
    else:
        axis.tick_params(labelbottom=False)
    style_statistics_axis(axis, tick_fontsize=10.5, density=True)


def plot_figure(
    example: pd.DataFrame,
    summary: pd.DataFrame,
    cadence: pd.DataFrame,
    residuals: pd.DataFrame,
    conditions: pd.DataFrame,
    thresholds: dict[str, float],
    night: str,
) -> plt.Figure:
    apply_statistics_style()
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    figure = plt.figure(figsize=STATISTICS_FOUR_PANEL_FIGSIZE)
    grid = figure.add_gridspec(2, 2, hspace=0.34, wspace=0.25)
    timeline_axis = figure.add_subplot(grid[0, 0])
    cadence_axis = figure.add_subplot(grid[0, 1])
    residual_axis = figure.add_subplot(grid[1, 0])
    condition_grid = grid[1, 1].subgridspec(3, 1, hspace=0.10)
    condition_axes = [
        figure.add_subplot(condition_grid[index, 0]) for index in range(3)
    ]
    x_limits = time_limits(example)

    plot_timeline(timeline_axis, example, summary.iloc[0], x_limits, night)
    plot_cadence(cadence_axis, cadence)
    plot_residuals(residual_axis, residuals, x_limits, night)

    third_condition = (
        "altitude"
        if (
            conditions["condition"].eq("altitude")
            & conditions["availability"].eq("archived")
        ).any()
        else "airmass"
    )
    condition_axes[0].set_title(
        "(d) Archived observing conditions",
        loc="left",
        fontsize=20,
        fontweight="bold",
        pad=13,
    )
    plot_condition_track(
        condition_axes[0],
        conditions,
        "sun_separation",
        x_limits,
        show_x=False,
    )
    plot_condition_track(
        condition_axes[1],
        conditions,
        "moon_separation",
        x_limits,
        show_x=False,
    )
    threshold = thresholds.get(
        "altitude_min" if third_condition == "altitude" else "airmass_max"
    )
    threshold_label = (
        f"archived minimum = {threshold:g}°"
        if third_condition == "altitude" and threshold is not None
        else (
            f"archived maximum = {threshold:g}"
            if threshold is not None
            else ""
        )
    )
    plot_condition_track(
        condition_axes[2],
        conditions,
        third_condition,
        x_limits,
        threshold=threshold,
        threshold_label=threshold_label,
        show_x=True,
    )

    figure.text(
        0.5,
        0.034,
        NEAR_SUN_NOTE,
        ha="center",
        va="bottom",
        fontsize=13.2,
        color=INK,
    )
    figure.text(
        0.5,
        0.015,
        ATTRIBUTION_NOTE + "  " + ARCHIVE_NOTE,
        ha="center",
        va="bottom",
        fontsize=11.5,
        color=MID_GREY,
    )
    figure.subplots_adjust(left=0.07, right=0.985, top=0.975, bottom=0.09)
    return figure


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (pd.Timestamp, datetime)):
        return iso_utc(value)
    if pd.isna(value) if not isinstance(value, (str, bytes)) else False:
        return None
    return value


def stage_csv(path: Path, frame: pd.DataFrame) -> tuple[Path, Path]:
    temporary = path.with_name(path.name + ".inprogress")
    frame.to_csv(temporary, index=False)
    return temporary, path


def stage_json(path: Path, payload: dict[str, Any]) -> tuple[Path, Path]:
    temporary = path.with_name(path.name + ".inprogress")
    temporary.write_text(
        json.dumps(json_safe(payload), indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return temporary, path


def output_paths(output_base: Path, figure_data_dir: Path) -> dict[str, Path]:
    stem = output_base.with_suffix("") if output_base.suffix else output_base
    paths = {
        "png": stem.with_suffix(".png"),
        "pdf": stem.with_suffix(".pdf"),
    }
    paths.update(
        {key: figure_data_dir / name for key, name in FIGURE_DATA_NAMES.items()}
    )
    return paths


def validate_no_collisions(paths: dict[str, Path]) -> None:
    collisions: list[Path] = []
    for path in paths.values():
        if path.exists() or path.with_name(path.name + ".inprogress").exists():
            collisions.append(path)
    if collisions:
        raise FileExistsError(
            "refusing to overwrite existing Fig. 4 outputs: "
            + ", ".join(str(path) for path in collisions)
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scheduler-example", type=Path, default=DEFAULT_EXAMPLE)
    parser.add_argument("--night-summary", type=Path, default=DEFAULT_NIGHT_SUMMARY)
    parser.add_argument("--night", default=DEFAULT_NIGHT)
    parser.add_argument("--l2-manifest", type=Path, default=None)
    parser.add_argument("--footprint", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT_BASE)
    parser.add_argument("--figure-data-dir", type=Path, default=DEFAULT_FIGURE_DATA_DIR)
    return parser


def resolve_inputs(args: argparse.Namespace) -> tuple[dict[str, Path], str, dict[str, Path]]:
    night = normalize_night(args.night)
    if len(night) != 8:
        raise ValueError("--night must be YYYYMMDD")
    inputs = {
        "scheduler_example": args.scheduler_example.expanduser().resolve(strict=False),
        "night_summary": args.night_summary.expanduser().resolve(strict=False),
    }
    for name in ("l2_manifest", "footprint", "config"):
        value = getattr(args, name)
        if value is not None:
            inputs[name] = value.expanduser().resolve(strict=False)
    for name, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(f"{name}: {path}")
    output_base = args.output_base.expanduser().resolve(strict=False)
    figure_data_dir = args.figure_data_dir.expanduser().resolve(strict=False)
    outputs = output_paths(output_base, figure_data_dir)
    validate_no_collisions(outputs)
    return inputs, night, outputs


def run(args: argparse.Namespace) -> None:
    inputs, night, outputs = resolve_inputs(args)
    hashes = {name: sha256_file(path) for name, path in inputs.items()}
    input_hash = combined_hash(hashes)

    example = read_example(inputs["scheduler_example"], night)
    summary = read_night_summary(inputs["night_summary"], night, example)
    l2 = read_optional_table(inputs["l2_manifest"]) if "l2_manifest" in inputs else None
    footprint = read_optional_table(inputs["footprint"]) if "footprint" in inputs else None
    config = (
        json.loads(inputs["config"].read_text(encoding="utf-8"))
        if "config" in inputs
        else None
    )
    if config is not None and not isinstance(config, dict):
        raise ValueError(f"{inputs['config']}: config root must be a JSON object")

    l2_mapping = l2_row_mapping(example, l2) if l2 is not None else {}
    timeline = build_timeline(example, input_hash)
    cadence = build_cadence(example, input_hash)
    residuals = build_residuals(example, input_hash)
    conditions = build_conditions(
        example,
        l2=l2,
        l2_mapping=l2_mapping,
        l2_path=inputs.get("l2_manifest"),
        footprint=footprint,
        footprint_path=inputs.get("footprint"),
        example_path=inputs["scheduler_example"],
        night=night,
        input_hash=input_hash,
    )
    summary_output = summary.copy()
    summary_output["input_manifest_hash"] = input_hash
    thresholds = config_thresholds(config)

    condition_availability = {
        condition: (
            "archived"
            if (
                conditions["condition"].eq(condition)
                & conditions["availability"].eq("archived")
            ).any()
            else "not_archived"
        )
        for condition in CONDITION_SPECS
    }
    metadata = {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "representative_night": night,
        "input_sources": {
            name: {"path": str(path), "sha256": hashes[name]}
            for name, path in inputs.items()
        },
        "combined_input_manifest_hash": input_hash,
        "condition_availability": condition_availability,
        "condition_thresholds_read_directly_from_config": thresholds,
        "counts": {
            "timeline_rows": len(timeline),
            "cadence_rows": len(cadence),
            "timed_residual_rows": len(residuals),
            "condition_value_rows": int(conditions["availability"].eq("archived").sum()),
        },
        "chart_contract": {
            "analytical_question": (
                "What cadence did the archived plan specify, how closely do archived "
                "starts follow it, and which observing conditions were actually archived?"
            ),
            "panels": {
                "a": "planned/actual timeline by scheduler cycle and mode",
                "b": "repeat-relative planned median and actual p16/median/p84 cadence",
                "c": "actual minus planned start residual at matched-pair grain",
                "d": "archived exposure-level Sun, Moon, and altitude/airmass tracks",
            },
            "palette_policy": "hard two-root cap: GOTTA blue normal, orange near-Sun, neutrals",
            "non_color_distinction": "marker shape, open fill, line style, and direct labels",
            "surface": "static 300-dpi PNG plus vector PDF",
        },
        "definitions": {
            "cadence_anchor": (
                "actual elapsed time is reported only for field/cycle groups with a "
                "matched, timed repeat-1 exposure"
            ),
            "plan_actual_residual": "actual_start_utc minus planned_start_utc",
            "conditions": (
                "direct archived values only; no ephemeris, coordinate, altitude, or "
                "airmass reconstruction"
            ),
        },
        "historical_near_sun_implementation": NEAR_SUN_NOTE,
        "attribution_limit": ATTRIBUTION_NOTE,
        "archive_limit": ARCHIVE_NOTE,
    }

    for path in outputs.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    staged: list[tuple[Path, Path]] = []
    staged.append(stage_csv(outputs["timeline"], timeline))
    staged.append(stage_csv(outputs["cadence"], cadence))
    staged.append(stage_csv(outputs["residuals"], residuals))
    staged.append(stage_csv(outputs["conditions"], conditions))
    staged.append(stage_csv(outputs["summary"], summary_output))
    staged.append(stage_json(outputs["metadata"], metadata))

    figure = plot_figure(
        example,
        summary,
        cadence,
        residuals,
        conditions,
        thresholds,
        night,
    )
    png_temporary = outputs["png"].with_name(outputs["png"].name + ".inprogress")
    pdf_temporary = outputs["pdf"].with_name(outputs["pdf"].name + ".inprogress")
    figure.savefig(
        png_temporary,
        format="png",
        dpi=300,
        bbox_inches="tight",
    )
    figure.savefig(pdf_temporary, format="pdf", bbox_inches="tight")
    plt.close(figure)
    staged.extend(
        [(png_temporary, outputs["png"]), (pdf_temporary, outputs["pdf"])]
    )
    for temporary, final in staged:
        os.replace(temporary, final)

    print(
        f"[done] night={night} plan={summary.iloc[0]['plan_exposure_n']} "
        f"matched={summary.iloc[0]['multiset_matched_frame_n']} "
        f"conditions={condition_availability} output={outputs['pdf']}",
        flush=True,
    )


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        run(args)
    except (
        FileExistsError,
        FileNotFoundError,
        RuntimeError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
