#!/usr/bin/env python3
"""Analyze scheduler plan realization without modifying production products.

The acquired-frame plan-compliance denominator is the strict raw science-frame
cohort on nights for which at least one valid archived plan row is available
and at least one strict raw frame was acquired.  Strict raw frames elsewhere in
the requested interval remain visible as ``outside_plan_archive`` counts, but
are not silently treated as noncompliant.  Counts are matched as a multiset on
``(night, field_id)``:

``matched = sum(min(planned_count[key], acquired_count[key]))``.

For row-level timing diagnostics, rows within each ``(night, field_id)`` group
are matched one-to-one while preserving chronological order.  When all planned
and actual timestamps are available, dynamic programming chooses the monotonic
maximum-cardinality assignment with minimum total absolute time residual.  If
either side lacks timestamps, rows are paired by scheduler/sequence order.  The
timing assignment never changes the multiset compliance count.

Actual timestamps use raw ``obs_time_utc`` first, then the joined L2
``obs_time_utc`` or MJD.  File mtimes are deliberately not treated as exposure
times.  Actual inter-exposure ``overhead`` is reported only as a wall-clock gap
proxy (successive start-time difference minus the previous exposure duration),
not as weather loss, equipment efficiency, or a hardware timing measurement.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


SCHEMA_VERSION = "1.0"
DEFAULT_START = "2025-11-15"
DEFAULT_END = "2026-07-15"
DEFAULT_EXAMPLE_NIGHT = "20260420"
MJD_EPOCH = datetime(1858, 11, 17, tzinfo=timezone.utc)

KNOWN_LIMITATION_CODES = (
    "near_sun_solar_separation_sort_overwritten_by_ra_sort",
    "scheduler_history_ingestion_stalled_after_20260330_server_audit",
)
NEAR_SUN_STATUS = (
    "mode_label_is_end_of_night_block; current code does not retain minimum-"
    "solar-elongation ordering because the candidate list is re-sorted by RA"
)
HISTORY_STATUS = (
    "production history feedback was observed stalled after 20260330; recent-"
    "night cooldown and current coverage feedback must not be assumed active"
)
ATTRIBUTION_LIMIT = (
    "plan realization and wall-clock gap diagnostics only; do not attribute "
    "unmatched frames or gaps to weather or equipment efficiency"
)


CONDITION_FIELDS = (
    "sky_adu",
    "sky_rms_adu",
    "sky_mag",
    "limiting_mag",
    "seeing_arcsec",
    "fwhm_header",
    "ccd_zeropoint",
    "ccd_zeropoint_rms",
    "photometric_rms",
    "n_sources",
    "n_wcs",
    "n_stars",
    "n_match_calibration",
    "n_psf",
)


EXPOSURE_OUTPUT_FIELDS = (
    "night",
    "plan_archive_available",
    "plan_active_acquired_cohort",
    "acquired_frame_in_plan_active_denominator",
    "acquired_frame_outside_plan_archive",
    "field_id",
    "pair_id",
    "match_status",
    "match_rule",
    "mode",
    "cycle_id",
    "repeat",
    "plan_row",
    "plan_manifest_line",
    "planned_start_utc",
    "actual_start_utc",
    "actual_time_source",
    "time_residual_s",
    "abs_time_residual_s",
    "planned_exptime_s",
    "actual_exptime_s",
    "planned_interval_prev_same_cycle_field_s",
    "planned_baseline_from_first_same_cycle_field_s",
    "actual_interval_prev_matched_same_cycle_field_s",
    "actual_baseline_from_first_matched_same_cycle_field_s",
    "planned_inter_start_prev_global_s",
    "planned_gap_proxy_prev_global_s",
    "actual_inter_start_prev_global_s",
    "actual_gap_proxy_prev_global_s",
    "raw_file_name",
    "raw_path",
    "raw_sequence_id",
    "raw_manifest_line",
    "l2_join_status",
    "l2_file_name",
    "l2_path",
    "l2_manifest_line",
    *CONDITION_FIELDS,
    "moon_phase_fraction",
    "moon_sep_deg",
    "planned_ra_deg",
    "planned_dec_deg",
    "known_limitation_codes",
    "input_manifest_hash",
)


NIGHT_OUTPUT_FIELDS = (
    "night",
    "plan_archive_available",
    "plan_active_acquired_cohort",
    "plan_exposure_n",
    "acquired_strict_frame_n",
    "plan_active_acquired_frame_denominator_n",
    "plan_active_cohort_matched_frame_n",
    "plan_active_acquired_not_matched_n",
    "full_interval_acquired_without_plan_archive_n",
    "multiset_matched_frame_n",
    "planned_not_acquired_n",
    "acquired_not_planned_n",
    "acquired_frame_plan_compliance",
    "planned_frame_realization",
    "planned_unique_field_n",
    "acquired_unique_field_n",
    "actual_time_available_n",
    "matched_pair_with_two_times_n",
    "abs_time_residual_median_s",
    "abs_time_residual_p90_s",
    "planned_revisit_interval_n",
    "planned_revisit_interval_median_s",
    "actual_revisit_interval_n",
    "actual_revisit_interval_median_s",
    "planned_revisit_baseline_n",
    "planned_revisit_baseline_median_s",
    "actual_revisit_baseline_n",
    "actual_revisit_baseline_median_s",
    "actual_gap_proxy_n",
    "actual_gap_proxy_median_s",
    "actual_gap_proxy_p90_s",
    "matched_with_l2_n",
    "normal_planned_n",
    "normal_matched_n",
    "normal_planned_realization",
    "near_sun_planned_n",
    "near_sun_matched_n",
    "near_sun_planned_realization",
    "followup_planned_n",
    "followup_matched_n",
    "followup_planned_realization",
    "quality_code",
    "quality_reason",
    "primary_science_included",
    "near_sun_selection_status",
    "history_feedback_status",
    "attribution_limit",
    "input_manifest_hash",
)


MODE_BASE_FIELDS = (
    "mode",
    "planned_n",
    "matched_n",
    "planned_realization",
    "planned_realization_denominator_n",
    "fraction_of_full_interval_acquired_frames_matching_mode",
    "full_interval_acquired_frame_denominator_n",
    "fraction_of_plan_active_acquired_frames_matching_mode",
    "plan_active_acquired_frame_denominator_n",
    "fraction_of_all_matched_frames_in_mode",
    "all_matched_frame_denominator_n",
    "fraction_of_all_planned_frames_in_mode",
    "all_planned_frame_denominator_n",
    "night_n",
    "field_n",
    "timed_pair_n",
    "abs_time_residual_median_s",
    "abs_time_residual_p16_s",
    "abs_time_residual_p84_s",
    "abs_time_residual_p90_s",
    "planned_revisit_interval_n",
    "planned_revisit_interval_median_s",
    "actual_revisit_interval_n",
    "actual_revisit_interval_median_s",
    "planned_revisit_baseline_n",
    "planned_revisit_baseline_median_s",
    "actual_revisit_baseline_n",
    "actual_revisit_baseline_median_s",
    "actual_gap_proxy_n",
    "actual_gap_proxy_median_s",
    "actual_gap_proxy_p90_s",
)


def condition_mode_fields() -> tuple[str, ...]:
    fields: list[str] = []
    for name in CONDITION_FIELDS:
        fields.extend(
            (
                f"{name}_availability",
                f"{name}_n",
                f"{name}_median",
                f"{name}_p16",
                f"{name}_p84",
            )
        )
        if name == "n_sources":
            fields.append("n_sources_total")
    return tuple(fields)


MODE_OUTPUT_FIELDS = (
    *MODE_BASE_FIELDS,
    *condition_mode_fields(),
    "condition_interpretation",
    "near_sun_selection_status",
    "history_feedback_status",
    "attribution_limit",
    "input_manifest_hash",
)


FIGURE_NIGHT_FIELDS = (
    "night",
    "plan_archive_available",
    "plan_active_acquired_cohort",
    "plan_exposure_n",
    "acquired_strict_frame_n",
    "plan_active_acquired_frame_denominator_n",
    "plan_active_cohort_matched_frame_n",
    "full_interval_acquired_without_plan_archive_n",
    "multiset_matched_frame_n",
    "acquired_frame_plan_compliance",
    "planned_frame_realization",
    "normal_planned_n",
    "normal_matched_n",
    "near_sun_planned_n",
    "near_sun_matched_n",
    "planned_revisit_interval_median_s",
    "actual_revisit_interval_median_s",
    "actual_gap_proxy_median_s",
    "quality_code",
    "input_manifest_hash",
)


@dataclass(frozen=True)
class CsvRecord:
    line_number: int
    values: dict[str, str]


@dataclass(frozen=True)
class CsvTable:
    path: Path
    fieldnames: tuple[str, ...]
    records: tuple[CsvRecord, ...]


def parse_iso_date(value: str) -> date:
    try:
        return datetime.strptime(value, "%Y-%m-%d").date()
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected YYYY-MM-DD, got {value!r}") from exc


def iter_nights(start: date, end: date) -> Iterable[str]:
    current = start
    while current <= end:
        yield current.strftime("%Y%m%d")
        current += timedelta(days=1)


def normalize_night(value: Any) -> str:
    digits = re.sub(r"\D", "", str(value or ""))
    return digits[:8] if len(digits) >= 8 else ""


def normalize_field_id(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if re.fullmatch(r"\d+(?:\.0+)?", text):
        return f"{int(float(text)):04d}"
    if re.fullmatch(r"MP_?\d+", text, flags=re.IGNORECASE):
        digits = re.sub(r"\D", "", text)
        return f"{int(digits):04d}"
    return text


def truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y", "t"}


def parse_float(value: Any) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def parse_int(value: Any, default: int = 0) -> int:
    numeric = parse_float(value)
    return int(numeric) if numeric is not None else default


def parse_utc(value: Any) -> datetime | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def parse_mjd(value: Any) -> datetime | None:
    numeric = parse_float(value)
    if numeric is None:
        return None
    return MJD_EPOCH + timedelta(days=numeric)


def iso_utc(value: datetime | None) -> str:
    if value is None:
        return ""
    return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv_table(path: Path) -> CsvTable:
    if not path.is_file():
        raise FileNotFoundError(path)
    records: list[CsvRecord] = []
    with path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        fieldnames = tuple(str(name) for name in reader.fieldnames)
        for line_number, row in enumerate(reader, start=2):
            records.append(
                CsvRecord(
                    line_number=line_number,
                    values={str(key): str(value or "") for key, value in row.items()},
                )
            )
    return CsvTable(path=path, fieldnames=fieldnames, records=tuple(records))


def strict_rows(table: CsvTable, flag_name: str) -> list[CsvRecord]:
    if flag_name not in table.fieldnames:
        return list(table.records)
    return [record for record in table.records if truthy(record.values.get(flag_name))]


def base_exposure_name(name: str) -> str:
    text = Path(str(name or "")).name
    text = re.sub(r"\.fits\.gz$", "", text, flags=re.IGNORECASE)
    text = re.sub(r"\.fits$", "", text, flags=re.IGNORECASE)
    text = re.sub(r"_cat$", "", text, flags=re.IGNORECASE)
    return text


def prepare_l2(
    table: CsvTable,
    allowed_nights: set[str],
) -> tuple[
    list[dict[str, Any]],
    dict[tuple[str, str, str], list[dict[str, Any]]],
    dict[tuple[str, str], list[dict[str, Any]]],
]:
    rows: list[dict[str, Any]] = []
    by_key: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    by_name: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in strict_rows(table, "strict_standard_catalog"):
        night = normalize_night(record.values.get("night"))
        if night not in allowed_nights:
            continue
        field_id = normalize_field_id(record.values.get("field_id"))
        sequence_id = normalize_field_id(record.values.get("sequence_id"))
        row: dict[str, Any] = {
            **record.values,
            "_line": record.line_number,
            "_night": night,
            "_field_id": field_id,
            "_sequence_id": sequence_id,
        }
        rows.append(row)
        by_key[(night, field_id, sequence_id)].append(row)
        by_name[(night, base_exposure_name(record.values.get("file_name", "")))].append(row)
    return rows, by_key, by_name


def choose_l2(candidates: Sequence[dict[str, Any]]) -> dict[str, Any] | None:
    if not candidates:
        return None
    return max(
        candidates,
        key=lambda row: (
            bool(parse_utc(row.get("obs_time_utc")) or parse_mjd(row.get("mjd"))),
            sum(parse_float(row.get(name)) is not None for name in CONDITION_FIELDS),
            parse_utc(row.get("mtime_utc")) or datetime.min.replace(tzinfo=timezone.utc),
            -int(row["_line"]),
        ),
    )


def prepare_actuals(
    raw_table: CsvTable,
    *,
    allowed_nights: set[str],
    l2_by_key: dict[tuple[str, str, str], list[dict[str, Any]]],
    l2_by_name: dict[tuple[str, str], list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, record in enumerate(strict_rows(raw_table, "strict_standard_science")):
        night = normalize_night(record.values.get("night"))
        if night not in allowed_nights:
            continue
        field_id = normalize_field_id(record.values.get("field_id"))
        if not field_id:
            continue
        sequence_id = normalize_field_id(record.values.get("sequence_id"))
        candidates = l2_by_key.get((night, field_id, sequence_id), [])
        if not candidates:
            candidates = l2_by_name.get(
                (night, base_exposure_name(record.values.get("file_name", ""))),
                [],
            )
        l2 = choose_l2(candidates)
        actual_dt = parse_utc(record.values.get("obs_time_utc"))
        time_source = "raw_obs_time_utc" if actual_dt else ""
        if actual_dt is None and l2 is not None:
            actual_dt = parse_utc(l2.get("obs_time_utc"))
            if actual_dt is not None:
                time_source = "l2_obs_time_utc"
            else:
                actual_dt = parse_mjd(l2.get("mjd"))
                if actual_dt is not None:
                    time_source = "l2_mjd"
        exposure_s = parse_float(record.values.get("exposure_s"))
        if exposure_s is None and l2 is not None:
            exposure_s = parse_float(l2.get("exposure_s"))
        sequence_order = parse_int(sequence_id, default=10**9)
        row: dict[str, Any] = {
            "_night": night,
            "_field_id": field_id,
            "_sequence_id": sequence_id,
            "_actual_order": index,
            "_sequence_order": sequence_order,
            "_actual_dt": actual_dt,
            "_actual_time_source": time_source,
            "_actual_exptime_s": exposure_s,
            "_raw": record.values,
            "_raw_line": record.line_number,
            "_l2": l2,
            "_actual_inter_start_s": None,
            "_actual_gap_proxy_s": None,
        }
        rows.append(row)

    by_night: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_night[row["_night"]].append(row)
    for night_rows in by_night.values():
        timed = sorted(
            (row for row in night_rows if row["_actual_dt"] is not None),
            key=lambda row: (row["_actual_dt"], row["_sequence_order"], row["_actual_order"]),
        )
        for previous, current in zip(timed, timed[1:]):
            inter_start = (current["_actual_dt"] - previous["_actual_dt"]).total_seconds()
            current["_actual_inter_start_s"] = inter_start
            if previous["_actual_exptime_s"] is not None:
                current["_actual_gap_proxy_s"] = inter_start - previous["_actual_exptime_s"]
    return rows


def prepare_plans(table: CsvTable, allowed_nights: set[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, record in enumerate(table.records):
        night = normalize_night(record.values.get("night"))
        if night not in allowed_nights:
            continue
        field_id = normalize_field_id(record.values.get("field_id"))
        if not field_id:
            continue
        plan_order = parse_int(record.values.get("plan_row"), default=index)
        mode = str(record.values.get("mode", "")).strip().lower() or "unknown"
        row: dict[str, Any] = {
            "_night": night,
            "_field_id": field_id,
            "_plan_order": plan_order,
            "_planned_dt": parse_utc(record.values.get("start_utc")),
            "_planned_exptime_s": parse_float(record.values.get("exptime_s")),
            "_mode": mode,
            "_cycle_id": str(record.values.get("cycle_id", "")).strip(),
            "_repeat": parse_int(record.values.get("repeat"), default=0),
            "_plan": record.values,
            "_plan_line": record.line_number,
            "_planned_inter_start_s": None,
            "_planned_gap_proxy_s": None,
            "_planned_revisit_interval_s": None,
            "_planned_revisit_baseline_s": None,
        }
        rows.append(row)

    by_night: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_night[row["_night"]].append(row)
    for night_rows in by_night.values():
        timed = sorted(
            (row for row in night_rows if row["_planned_dt"] is not None),
            key=lambda row: (row["_planned_dt"], row["_plan_order"]),
        )
        for previous, current in zip(timed, timed[1:]):
            inter_start = (current["_planned_dt"] - previous["_planned_dt"]).total_seconds()
            current["_planned_inter_start_s"] = inter_start
            if previous["_planned_exptime_s"] is not None:
                current["_planned_gap_proxy_s"] = inter_start - previous["_planned_exptime_s"]

    revisit_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        revisit_groups[(row["_night"], row["_field_id"], row["_cycle_id"])].append(row)
    for group in revisit_groups.values():
        ordered = sorted(
            group,
            key=lambda row: (
                row["_planned_dt"] or datetime.max.replace(tzinfo=timezone.utc),
                row["_repeat"],
                row["_plan_order"],
            ),
        )
        first = ordered[0]["_planned_dt"] if ordered else None
        previous = None
        for row in ordered:
            current = row["_planned_dt"]
            if previous is not None and current is not None:
                row["_planned_revisit_interval_s"] = (current - previous).total_seconds()
            if first is not None and current is not None:
                row["_planned_revisit_baseline_s"] = (current - first).total_seconds()
            if current is not None:
                previous = current
    return rows


def sort_plan(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            row["_planned_dt"] or datetime.max.replace(tzinfo=timezone.utc),
            row["_plan_order"],
        ),
    )


def sort_actual(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            row["_actual_dt"] or datetime.max.replace(tzinfo=timezone.utc),
            row["_sequence_order"],
            row["_actual_order"],
        ),
    )


def sort_plan_by_scheduler_order(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            row["_plan_order"],
            row["_planned_dt"] or datetime.max.replace(tzinfo=timezone.utc),
        ),
    )


def sort_actual_by_sequence_order(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            row["_sequence_order"],
            row["_actual_order"],
        ),
    )


def minimum_cost_subsequence_pairs(
    short_times: Sequence[datetime],
    long_times: Sequence[datetime],
) -> list[tuple[int, int]]:
    """Match every short item to an ordered subsequence of the long items."""

    n_short = len(short_times)
    n_long = len(long_times)
    infinity = float("inf")
    dp = [[infinity] * (n_long + 1) for _ in range(n_short + 1)]
    take = [[False] * (n_long + 1) for _ in range(n_short + 1)]
    for j in range(n_long + 1):
        dp[0][j] = 0.0
    for i in range(1, n_short + 1):
        for j in range(1, n_long + 1):
            skip_cost = dp[i][j - 1]
            match_cost = dp[i - 1][j - 1] + abs(
                (short_times[i - 1] - long_times[j - 1]).total_seconds()
            )
            if match_cost <= skip_cost:
                dp[i][j] = match_cost
                take[i][j] = True
            else:
                dp[i][j] = skip_cost
    pairs: list[tuple[int, int]] = []
    i, j = n_short, n_long
    while i > 0 and j > 0:
        if take[i][j]:
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        else:
            j -= 1
    if i != 0:
        raise RuntimeError("failed to reconstruct monotonic timing assignment")
    pairs.reverse()
    return pairs


def match_group(
    plan_rows: Sequence[dict[str, Any]],
    actual_rows: Sequence[dict[str, Any]],
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[tuple[int, int]],
    str,
]:
    pair_count = min(len(plan_rows), len(actual_rows))
    if pair_count == 0:
        return list(plan_rows), list(actual_rows), [], "no_counterpart"
    all_timed = all(row["_planned_dt"] is not None for row in plan_rows) and all(
        row["_actual_dt"] is not None for row in actual_rows
    )
    if not all_timed:
        plans = sort_plan_by_scheduler_order(plan_rows)
        actuals = sort_actual_by_sequence_order(actual_rows)
        pairs = [(index, index) for index in range(pair_count)]
        return plans, actuals, pairs, "order_only_missing_timestamp"
    plans = sort_plan(plan_rows)
    actuals = sort_actual(actual_rows)
    if len(plans) <= len(actuals):
        raw_pairs = minimum_cost_subsequence_pairs(
            [row["_planned_dt"] for row in plans],
            [row["_actual_dt"] for row in actuals],
        )
        return (
            plans,
            actuals,
            raw_pairs,
            "monotonic_minimum_total_absolute_time_residual",
        )
    raw_pairs = minimum_cost_subsequence_pairs(
        [row["_actual_dt"] for row in actuals],
        [row["_planned_dt"] for row in plans],
    )
    pairs = [(plan_index, actual_index) for actual_index, plan_index in raw_pairs]
    return (
        plans,
        actuals,
        pairs,
        "monotonic_minimum_total_absolute_time_residual",
    )


def output_row(
    *,
    plan: dict[str, Any] | None,
    actual: dict[str, Any] | None,
    pair_id: str,
    status: str,
    rule: str,
    input_hash: str,
) -> dict[str, Any]:
    night = plan["_night"] if plan is not None else actual["_night"]
    field_id = plan["_field_id"] if plan is not None else actual["_field_id"]
    planned_dt = plan["_planned_dt"] if plan is not None else None
    actual_dt = actual["_actual_dt"] if actual is not None else None
    residual = (
        (actual_dt - planned_dt).total_seconds()
        if planned_dt is not None and actual_dt is not None
        else None
    )
    l2 = actual["_l2"] if actual is not None else None
    raw = actual["_raw"] if actual is not None else {}
    plan_values = plan["_plan"] if plan is not None else {}
    row: dict[str, Any] = {
        "night": night,
        "field_id": field_id,
        "pair_id": pair_id,
        "match_status": status,
        "match_rule": rule,
        "mode": plan["_mode"] if plan is not None else "unplanned",
        "cycle_id": plan["_cycle_id"] if plan is not None else "",
        "repeat": plan["_repeat"] if plan is not None else "",
        "plan_row": plan["_plan_order"] if plan is not None else "",
        "plan_manifest_line": plan["_plan_line"] if plan is not None else "",
        "planned_start_utc": iso_utc(planned_dt),
        "actual_start_utc": iso_utc(actual_dt),
        "actual_time_source": actual["_actual_time_source"] if actual is not None else "",
        "time_residual_s": residual,
        "abs_time_residual_s": abs(residual) if residual is not None else None,
        "planned_exptime_s": plan["_planned_exptime_s"] if plan is not None else None,
        "actual_exptime_s": actual["_actual_exptime_s"] if actual is not None else None,
        "planned_interval_prev_same_cycle_field_s": (
            plan["_planned_revisit_interval_s"] if plan is not None else None
        ),
        "planned_baseline_from_first_same_cycle_field_s": (
            plan["_planned_revisit_baseline_s"] if plan is not None else None
        ),
        "actual_interval_prev_matched_same_cycle_field_s": None,
        "actual_baseline_from_first_matched_same_cycle_field_s": None,
        "planned_inter_start_prev_global_s": (
            plan["_planned_inter_start_s"] if plan is not None else None
        ),
        "planned_gap_proxy_prev_global_s": (
            plan["_planned_gap_proxy_s"] if plan is not None else None
        ),
        "actual_inter_start_prev_global_s": (
            actual["_actual_inter_start_s"] if actual is not None else None
        ),
        "actual_gap_proxy_prev_global_s": (
            actual["_actual_gap_proxy_s"] if actual is not None else None
        ),
        "raw_file_name": raw.get("file_name", ""),
        "raw_path": raw.get("path", ""),
        "raw_sequence_id": actual["_sequence_id"] if actual is not None else "",
        "raw_manifest_line": actual["_raw_line"] if actual is not None else "",
        "l2_join_status": (
            "matched" if l2 is not None else ("missing" if actual else "not_applicable")
        ),
        "l2_file_name": l2.get("file_name", "") if l2 is not None else "",
        "l2_path": l2.get("path", "") if l2 is not None else "",
        "l2_manifest_line": l2.get("_line", "") if l2 is not None else "",
        "moon_phase_fraction": parse_float(plan_values.get("moon_phase_fraction")),
        "moon_sep_deg": parse_float(plan_values.get("moon_sep_deg")),
        "planned_ra_deg": parse_float(plan_values.get("ra")),
        "planned_dec_deg": parse_float(plan_values.get("dec")),
        "known_limitation_codes": ";".join(KNOWN_LIMITATION_CODES),
        "input_manifest_hash": input_hash,
    }
    for name in CONDITION_FIELDS:
        row[name] = parse_float(l2.get(name)) if l2 is not None else None
    return row


def build_matches(
    plans: Sequence[dict[str, Any]],
    actuals: Sequence[dict[str, Any]],
    *,
    input_hash: str,
) -> list[dict[str, Any]]:
    plan_archive_nights = {row["_night"] for row in plans}
    acquired_nights = {row["_night"] for row in actuals}
    plan_active_acquired_nights = plan_archive_nights & acquired_nights
    plans_by_key: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    actuals_by_key: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in plans:
        plans_by_key[(row["_night"], row["_field_id"])].append(row)
    for row in actuals:
        actuals_by_key[(row["_night"], row["_field_id"])].append(row)

    output: list[dict[str, Any]] = []
    for night, field_id in sorted(set(plans_by_key) | set(actuals_by_key)):
        plan_group, actual_group, pairs, rule = match_group(
            plans_by_key.get((night, field_id), []),
            actuals_by_key.get((night, field_id), []),
        )
        matched_plan = {plan_index for plan_index, _ in pairs}
        matched_actual = {actual_index for _, actual_index in pairs}
        for pair_number, (plan_index, actual_index) in enumerate(pairs, start=1):
            output.append(
                output_row(
                    plan=plan_group[plan_index],
                    actual=actual_group[actual_index],
                    pair_id=f"{night}:{field_id}:{pair_number:04d}",
                    status="matched",
                    rule=rule,
                    input_hash=input_hash,
                )
            )
        for plan_index, plan in enumerate(plan_group):
            if plan_index in matched_plan:
                continue
            output.append(
                output_row(
                    plan=plan,
                    actual=None,
                    pair_id=f"{night}:{field_id}:plan:{plan_index:04d}",
                    status="planned_not_acquired",
                    rule=rule,
                    input_hash=input_hash,
                )
            )
        for actual_index, actual in enumerate(actual_group):
            if actual_index in matched_actual:
                continue
            output.append(
                output_row(
                    plan=None,
                    actual=actual,
                    pair_id=f"{night}:{field_id}:actual:{actual_index:04d}",
                    status="acquired_not_planned",
                    rule=rule,
                    input_hash=input_hash,
                )
            )

    revisit_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in output:
        if row["match_status"] == "matched":
            revisit_groups[(row["night"], row["field_id"], str(row["cycle_id"]))].append(row)
    for group in revisit_groups.values():
        ordered = sorted(
            group,
            key=lambda row: (
                parse_utc(row["planned_start_utc"])
                or datetime.max.replace(tzinfo=timezone.utc),
                parse_int(row["repeat"]),
                parse_int(row["plan_row"]),
            ),
        )
        first_actual = parse_utc(ordered[0]["actual_start_utc"]) if ordered else None
        previous_actual = None
        for row in ordered:
            current_actual = parse_utc(row["actual_start_utc"])
            if previous_actual is not None and current_actual is not None:
                row["actual_interval_prev_matched_same_cycle_field_s"] = (
                    current_actual - previous_actual
                ).total_seconds()
            if first_actual is not None and current_actual is not None:
                row["actual_baseline_from_first_matched_same_cycle_field_s"] = (
                    current_actual - first_actual
                ).total_seconds()
            if current_actual is not None:
                previous_actual = current_actual

    for row in output:
        has_acquired_frame = row["match_status"] in {"matched", "acquired_not_planned"}
        archive_available = row["night"] in plan_archive_nights
        row["plan_archive_available"] = archive_available
        row["plan_active_acquired_cohort"] = (
            row["night"] in plan_active_acquired_nights
        )
        row["acquired_frame_in_plan_active_denominator"] = (
            has_acquired_frame and archive_available
        )
        row["acquired_frame_outside_plan_archive"] = (
            has_acquired_frame and not archive_available
        )

    return sorted(
        output,
        key=lambda row: (
            row["night"],
            parse_utc(row["planned_start_utc"])
            or parse_utc(row["actual_start_utc"])
            or datetime.max.replace(tzinfo=timezone.utc),
            row["field_id"],
            row["pair_id"],
        ),
    )


def percentile(values: Sequence[float], quantile: float) -> float | None:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return None
    if len(finite) == 1:
        return float(finite[0])
    position = (len(finite) - 1) * quantile
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return float(finite[low])
    fraction = position - low
    return float(finite[low] * (1.0 - fraction) + finite[high] * fraction)


def numeric_values(rows: Iterable[dict[str, Any]], column: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        numeric = parse_float(row.get(column))
        if numeric is not None:
            values.append(numeric)
    return values


def revisit_endpoint_baselines(
    rows: Iterable[dict[str, Any]],
    column: str,
) -> list[float]:
    """Return one end-to-end revisit baseline per night/field/cycle group."""

    maxima: dict[tuple[str, str, str], float] = {}
    for row in rows:
        value = parse_float(row.get(column))
        if value is None or value <= 0:
            continue
        key = (
            str(row.get("night", "")),
            str(row.get("field_id", "")),
            str(row.get("cycle_id", "")),
        )
        maxima[key] = max(value, maxima.get(key, value))
    return list(maxima.values())


def safe_ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def build_cohort_accounting(
    plans: Sequence[dict[str, Any]],
    actuals: Sequence[dict[str, Any]],
    matches: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    """Reconcile full-interval and plan-active acquired-frame cohorts."""

    plan_archive_nights = {row["_night"] for row in plans}
    acquired_nights = {row["_night"] for row in actuals}
    plan_active_nights = plan_archive_nights & acquired_nights
    matched_rows = [row for row in matches if row["match_status"] == "matched"]
    plan_active_actuals = [
        row for row in actuals if row["_night"] in plan_archive_nights
    ]
    outside_archive_actuals = [
        row for row in actuals if row["_night"] not in plan_archive_nights
    ]
    cohort_matched = [
        row for row in matched_rows if row["night"] in plan_active_nights
    ]
    matched_outside_archive = [
        row for row in matched_rows if row["night"] not in plan_archive_nights
    ]
    if matched_outside_archive:
        raise ValueError("matched rows exist outside the archived-plan night set")
    if len(cohort_matched) != len(matched_rows):
        raise ValueError("matched rows do not close within the plan-active acquired cohort")

    plan_active_denominator = len(plan_active_actuals)
    cohort_matched_n = len(cohort_matched)
    plan_active_not_matched = plan_active_denominator - cohort_matched_n
    if plan_active_not_matched < 0:
        raise ValueError("plan-active matched count exceeds its acquired-frame denominator")
    full_interval_unmatched = len(actuals) - cohort_matched_n
    if full_interval_unmatched != len(outside_archive_actuals) + plan_active_not_matched:
        raise ValueError("full-interval acquired-frame accounting does not close")

    return {
        "plan_archive_available_night_n": len(plan_archive_nights),
        "plan_archive_first_night": min(plan_archive_nights) if plan_archive_nights else None,
        "plan_archive_last_night": max(plan_archive_nights) if plan_archive_nights else None,
        "acquired_night_n_full_interval": len(acquired_nights),
        "plan_active_acquired_night_n": len(plan_active_nights),
        "plan_archive_without_acquired_night_n": len(plan_archive_nights - acquired_nights),
        "acquired_without_plan_archive_night_n": len(acquired_nights - plan_archive_nights),
        "strict_raw_frame_n_full_interval": len(actuals),
        "plan_active_acquired_frame_denominator_n": plan_active_denominator,
        "plan_active_cohort_matched_frame_n": cohort_matched_n,
        "plan_active_acquired_not_matched_n": plan_active_not_matched,
        "full_interval_acquired_without_plan_archive_n": len(outside_archive_actuals),
        "full_interval_acquired_not_matched_n": full_interval_unmatched,
        "acquired_frame_plan_compliance": safe_ratio(
            cohort_matched_n, plan_active_denominator
        ),
        "plan_row_n": len(plans),
        "planned_not_acquired_n": len(plans) - cohort_matched_n,
        "planned_frame_realization": safe_ratio(cohort_matched_n, len(plans)),
    }


def load_night_status(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    table = read_csv_table(path)
    return {
        normalize_night(record.values.get("night")): record.values
        for record in table.records
        if normalize_night(record.values.get("night"))
    }


def build_night_summary(
    nights: Sequence[str],
    plans: Sequence[dict[str, Any]],
    actuals: Sequence[dict[str, Any]],
    matches: Sequence[dict[str, Any]],
    *,
    night_status: dict[str, dict[str, str]],
    input_hash: str,
) -> list[dict[str, Any]]:
    plans_by_night: dict[str, list[dict[str, Any]]] = defaultdict(list)
    actuals_by_night: dict[str, list[dict[str, Any]]] = defaultdict(list)
    matches_by_night: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in plans:
        plans_by_night[row["_night"]].append(row)
    for row in actuals:
        actuals_by_night[row["_night"]].append(row)
    for row in matches:
        matches_by_night[row["night"]].append(row)

    plan_archive_nights = set(plans_by_night)
    acquired_nights = set(actuals_by_night)
    plan_active_acquired_nights = plan_archive_nights & acquired_nights

    output: list[dict[str, Any]] = []
    for night in nights:
        night_plans = plans_by_night.get(night, [])
        night_actuals = actuals_by_night.get(night, [])
        night_matches = matches_by_night.get(night, [])
        matched = [row for row in night_matches if row["match_status"] == "matched"]
        residuals = numeric_values(matched, "abs_time_residual_s")
        planned_intervals = numeric_values(
            night_matches, "planned_interval_prev_same_cycle_field_s"
        )
        actual_intervals = numeric_values(
            matched, "actual_interval_prev_matched_same_cycle_field_s"
        )
        planned_baselines = revisit_endpoint_baselines(
            night_matches, "planned_baseline_from_first_same_cycle_field_s"
        )
        actual_baselines = revisit_endpoint_baselines(
            matched, "actual_baseline_from_first_matched_same_cycle_field_s"
        )
        actual_gaps = numeric_values(night_matches, "actual_gap_proxy_prev_global_s")
        mode_plan = Counter(row["_mode"] for row in night_plans)
        mode_match = Counter(row["mode"] for row in matched)
        status = night_status.get(night, {})
        plan_archive_available = night in plan_archive_nights
        plan_active_acquired = night in plan_active_acquired_nights
        plan_active_denominator = len(night_actuals) if plan_archive_available else 0
        cohort_matched_n = len(matched) if plan_active_acquired else 0
        row: dict[str, Any] = {
            "night": night,
            "plan_archive_available": plan_archive_available,
            "plan_active_acquired_cohort": plan_active_acquired,
            "plan_exposure_n": len(night_plans),
            "acquired_strict_frame_n": len(night_actuals),
            "plan_active_acquired_frame_denominator_n": plan_active_denominator,
            "plan_active_cohort_matched_frame_n": cohort_matched_n,
            "plan_active_acquired_not_matched_n": (
                plan_active_denominator - cohort_matched_n
            ),
            "full_interval_acquired_without_plan_archive_n": (
                0 if plan_archive_available else len(night_actuals)
            ),
            "multiset_matched_frame_n": len(matched),
            "planned_not_acquired_n": len(night_plans) - len(matched),
            "acquired_not_planned_n": len(night_actuals) - len(matched),
            "acquired_frame_plan_compliance": (
                safe_ratio(cohort_matched_n, plan_active_denominator)
                if plan_active_acquired
                else None
            ),
            "planned_frame_realization": safe_ratio(len(matched), len(night_plans)),
            "planned_unique_field_n": len({row["_field_id"] for row in night_plans}),
            "acquired_unique_field_n": len({row["_field_id"] for row in night_actuals}),
            "actual_time_available_n": sum(row["_actual_dt"] is not None for row in night_actuals),
            "matched_pair_with_two_times_n": len(residuals),
            "abs_time_residual_median_s": percentile(residuals, 0.50),
            "abs_time_residual_p90_s": percentile(residuals, 0.90),
            "planned_revisit_interval_n": len(planned_intervals),
            "planned_revisit_interval_median_s": percentile(planned_intervals, 0.50),
            "actual_revisit_interval_n": len(actual_intervals),
            "actual_revisit_interval_median_s": percentile(actual_intervals, 0.50),
            "planned_revisit_baseline_n": len(planned_baselines),
            "planned_revisit_baseline_median_s": percentile(planned_baselines, 0.50),
            "actual_revisit_baseline_n": len(actual_baselines),
            "actual_revisit_baseline_median_s": percentile(actual_baselines, 0.50),
            "actual_gap_proxy_n": len(actual_gaps),
            "actual_gap_proxy_median_s": percentile(actual_gaps, 0.50),
            "actual_gap_proxy_p90_s": percentile(actual_gaps, 0.90),
            "matched_with_l2_n": sum(row["l2_join_status"] == "matched" for row in matched),
            "quality_code": status.get("quality_code", ""),
            "quality_reason": status.get("quality_reason", ""),
            "primary_science_included": status.get("primary_science_included", ""),
            "near_sun_selection_status": NEAR_SUN_STATUS,
            "history_feedback_status": HISTORY_STATUS,
            "attribution_limit": ATTRIBUTION_LIMIT,
            "input_manifest_hash": input_hash,
        }
        for mode in ("normal", "near_sun", "followup"):
            row[f"{mode}_planned_n"] = mode_plan[mode]
            row[f"{mode}_matched_n"] = mode_match[mode]
            row[f"{mode}_planned_realization"] = safe_ratio(mode_match[mode], mode_plan[mode])
        output.append(row)
    return output


def availability(column: str, l2_fieldnames: set[str], values: Sequence[float]) -> str:
    if column not in l2_fieldnames:
        return "NA_column_missing"
    if not values:
        return "NA_no_matched_l2_values"
    return "available"


def build_mode_summary(
    plans: Sequence[dict[str, Any]],
    actuals: Sequence[dict[str, Any]],
    matches: Sequence[dict[str, Any]],
    *,
    l2_fieldnames: set[str],
    input_hash: str,
) -> list[dict[str, Any]]:
    cohort = build_cohort_accounting(plans, actuals, matches)
    all_acquired = cohort["strict_raw_frame_n_full_interval"]
    plan_active_acquired = cohort["plan_active_acquired_frame_denominator_n"]
    all_matched = cohort["plan_active_cohort_matched_frame_n"]
    all_planned = cohort["plan_row_n"]
    modes = {"normal", "near_sun"}
    modes.update(row["_mode"] for row in plans)
    output: list[dict[str, Any]] = []
    for mode in sorted(modes, key=lambda value: (value not in {"normal", "near_sun"}, value)):
        mode_plans = [row for row in plans if row["_mode"] == mode]
        mode_rows = [row for row in matches if row["mode"] == mode]
        matched = [row for row in mode_rows if row["match_status"] == "matched"]
        residuals = numeric_values(matched, "abs_time_residual_s")
        planned_intervals = numeric_values(mode_rows, "planned_interval_prev_same_cycle_field_s")
        actual_intervals = numeric_values(
            matched, "actual_interval_prev_matched_same_cycle_field_s"
        )
        planned_baselines = revisit_endpoint_baselines(
            mode_rows, "planned_baseline_from_first_same_cycle_field_s"
        )
        actual_baselines = revisit_endpoint_baselines(
            matched, "actual_baseline_from_first_matched_same_cycle_field_s"
        )
        actual_gaps = numeric_values(matched, "actual_gap_proxy_prev_global_s")
        row: dict[str, Any] = {
            "mode": mode,
            "planned_n": len(mode_plans),
            "matched_n": len(matched),
            "planned_realization": safe_ratio(len(matched), len(mode_plans)),
            "planned_realization_denominator_n": len(mode_plans),
            "fraction_of_full_interval_acquired_frames_matching_mode": safe_ratio(
                len(matched), all_acquired
            ),
            "full_interval_acquired_frame_denominator_n": all_acquired,
            "fraction_of_plan_active_acquired_frames_matching_mode": safe_ratio(
                len(matched), plan_active_acquired
            ),
            "plan_active_acquired_frame_denominator_n": plan_active_acquired,
            "fraction_of_all_matched_frames_in_mode": safe_ratio(
                len(matched), all_matched
            ),
            "all_matched_frame_denominator_n": all_matched,
            "fraction_of_all_planned_frames_in_mode": safe_ratio(
                len(mode_plans), all_planned
            ),
            "all_planned_frame_denominator_n": all_planned,
            "night_n": len({item["_night"] for item in mode_plans}),
            "field_n": len({item["_field_id"] for item in mode_plans}),
            "timed_pair_n": len(residuals),
            "abs_time_residual_median_s": percentile(residuals, 0.50),
            "abs_time_residual_p16_s": percentile(residuals, 0.16),
            "abs_time_residual_p84_s": percentile(residuals, 0.84),
            "abs_time_residual_p90_s": percentile(residuals, 0.90),
            "planned_revisit_interval_n": len(planned_intervals),
            "planned_revisit_interval_median_s": percentile(planned_intervals, 0.50),
            "actual_revisit_interval_n": len(actual_intervals),
            "actual_revisit_interval_median_s": percentile(actual_intervals, 0.50),
            "planned_revisit_baseline_n": len(planned_baselines),
            "planned_revisit_baseline_median_s": percentile(planned_baselines, 0.50),
            "actual_revisit_baseline_n": len(actual_baselines),
            "actual_revisit_baseline_median_s": percentile(actual_baselines, 0.50),
            "actual_gap_proxy_n": len(actual_gaps),
            "actual_gap_proxy_median_s": percentile(actual_gaps, 0.50),
            "actual_gap_proxy_p90_s": percentile(actual_gaps, 0.90),
            "condition_interpretation": (
                "matched-exposure L2 metadata only; n_sources is L2 catalog row count, "
                "not asteroid yield or survey efficiency"
            ),
            "near_sun_selection_status": NEAR_SUN_STATUS,
            "history_feedback_status": HISTORY_STATUS,
            "attribution_limit": ATTRIBUTION_LIMIT,
            "input_manifest_hash": input_hash,
        }
        for column in CONDITION_FIELDS:
            values = numeric_values(matched, column)
            row[f"{column}_availability"] = availability(column, l2_fieldnames, values)
            row[f"{column}_n"] = len(values)
            row[f"{column}_median"] = percentile(values, 0.50) if values else "NA"
            row[f"{column}_p16"] = percentile(values, 0.16) if values else "NA"
            row[f"{column}_p84"] = percentile(values, 0.84) if values else "NA"
            if column == "n_sources":
                row["n_sources_total"] = sum(values) if values else "NA"
        output.append(row)
    return output


def json_safe_row(row: dict[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in row.items():
        if value == "NA":
            result[key] = None
        elif isinstance(value, float) and not math.isfinite(value):
            result[key] = None
        else:
            result[key] = value
    return result


def atomic_write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    temporary = path.with_name(path.name + ".inprogress")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {key: "" if row.get(key) is None else row.get(key, "") for key in fields}
            )
    os.replace(temporary, path)


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(path.name + ".inprogress")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory-dir", type=Path, required=True)
    parser.add_argument("--night-status", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-data-dir", type=Path, required=True)
    parser.add_argument("--start", type=parse_iso_date, default=parse_iso_date(DEFAULT_START))
    parser.add_argument("--end", type=parse_iso_date, default=parse_iso_date(DEFAULT_END))
    parser.add_argument("--example-night", default=DEFAULT_EXAMPLE_NIGHT)
    return parser


def validate_args(args: argparse.Namespace) -> tuple[Path, Path, Path, Path, Path, Path | None]:
    if args.start > args.end:
        raise ValueError("--start is after --end")
    example_night = normalize_night(args.example_night)
    if len(example_night) != 8:
        raise ValueError("--example-night must be YYYYMMDD")
    args.example_night = example_night
    inventory = args.inventory_dir.expanduser().resolve(strict=False)
    output = args.output_dir.expanduser().resolve(strict=False)
    figure_data = args.figure_data_dir.expanduser().resolve(strict=False)
    plan_path = inventory / "plan_rows.csv"
    raw_path = inventory / "raw_manifest.csv"
    l2_path = inventory / "l2_manifest.csv"
    for path in (plan_path, raw_path, l2_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    night_status = (
        args.night_status.expanduser().resolve(strict=False) if args.night_status else None
    )
    if night_status is not None and not night_status.is_file():
        raise FileNotFoundError(night_status)

    outputs = (
        output / "scheduler_plan_realization_by_night.csv",
        output / "scheduler_exposure_matches.csv",
        output / "scheduler_mode_summary.csv",
        output / "scheduler_mode_summary.json",
        output / f"scheduler_example_{args.example_night}.csv",
        figure_data / "fig04_scheduler_plan_realization_by_night.csv",
        figure_data / "fig04_scheduler_mode_summary.csv",
        figure_data / f"fig04_scheduler_example_{args.example_night}.csv",
    )
    collisions = [
        path
        for path in outputs
        if path.exists() or path.with_name(path.name + ".inprogress").exists()
    ]
    if collisions:
        raise FileExistsError(
            "refusing to overwrite existing scheduler-analysis outputs: "
            + ", ".join(str(path) for path in collisions)
        )
    output.mkdir(parents=True, exist_ok=True)
    figure_data.mkdir(parents=True, exist_ok=True)
    return plan_path, raw_path, l2_path, output, figure_data, night_status


def run(args: argparse.Namespace) -> None:
    plan_path, raw_path, l2_path, output, figure_data, night_status_path = validate_args(args)
    plan_table = read_csv_table(plan_path)
    raw_table = read_csv_table(raw_path)
    l2_table = read_csv_table(l2_path)
    for table, required in (
        (plan_table, {"night", "field_id"}),
        (raw_table, {"night", "field_id"}),
        (l2_table, {"night", "field_id"}),
    ):
        missing = sorted(required - set(table.fieldnames))
        if missing:
            raise ValueError(f"{table.path} missing required columns: {missing}")

    input_paths = [plan_path, raw_path, l2_path]
    if night_status_path is not None:
        input_paths.append(night_status_path)
    input_hashes = {path.name: sha256_file(path) for path in input_paths}
    combined_hash_text = "".join(
        f"{name}:{digest}\n" for name, digest in sorted(input_hashes.items())
    )
    combined_hash = hashlib.sha256(combined_hash_text.encode("ascii")).hexdigest()

    nights = list(iter_nights(args.start, args.end))
    allowed_nights = set(nights)
    _l2_rows, l2_by_key, l2_by_name = prepare_l2(l2_table, allowed_nights)
    plans = prepare_plans(plan_table, allowed_nights)
    actuals = prepare_actuals(
        raw_table,
        allowed_nights=allowed_nights,
        l2_by_key=l2_by_key,
        l2_by_name=l2_by_name,
    )
    matches = build_matches(plans, actuals, input_hash=combined_hash)
    cohort_accounting = build_cohort_accounting(plans, actuals, matches)
    night_status = load_night_status(night_status_path)
    night_rows = build_night_summary(
        nights,
        plans,
        actuals,
        matches,
        night_status=night_status,
        input_hash=combined_hash,
    )
    mode_rows = build_mode_summary(
        plans,
        actuals,
        matches,
        l2_fieldnames=set(l2_table.fieldnames),
        input_hash=combined_hash,
    )
    example_rows = [row for row in matches if row["night"] == args.example_night]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": iso_utc(datetime.now(timezone.utc)),
        "date_range_inclusive": {"start": args.start.isoformat(), "end": args.end.isoformat()},
        "example_night": args.example_night,
        "input_hashes": input_hashes,
        "combined_input_manifest_hash": combined_hash,
        "input_columns": {
            "plan_rows.csv": list(plan_table.fieldnames),
            "raw_manifest.csv": list(raw_table.fieldnames),
            "l2_manifest.csv": list(l2_table.fieldnames),
        },
        "counts": {
            "plan_rows_within_interval": len(plans),
            "strict_raw_rows_full_interval": len(actuals),
            "exposure_accounting_rows": len(matches),
            "matched_rows_plan_active_acquired_cohort": cohort_accounting[
                "plan_active_cohort_matched_frame_n"
            ],
            "plan_active_acquired_frame_denominator_n": cohort_accounting[
                "plan_active_acquired_frame_denominator_n"
            ],
            "full_interval_acquired_without_plan_archive_n": cohort_accounting[
                "full_interval_acquired_without_plan_archive_n"
            ],
            "nights_full_interval": len(nights),
        },
        "cohort_accounting": cohort_accounting,
        "definitions": {
            "plan_archive_available": (
                "true for a night with one or more valid rows in the frozen archived "
                "plan_rows.csv within the requested interval; absence of a row is not "
                "interpreted as weather or equipment state"
            ),
            "plan_active_acquired_cohort": (
                "nights with plan_archive_available=true and at least one strict raw "
                "science frame acquired"
            ),
            "acquired_frame_plan_compliance": (
                "sum over (night,field_id) of min(planned count, acquired strict raw count), "
                "divided only by strict raw frames on plan-active acquired-cohort nights"
            ),
            "full_interval_acquired_without_plan_archive": (
                "strict raw frames on nights without an archived plan row; retained as "
                "not-covered accounting and excluded from acquired-frame compliance"
            ),
            "planned_frame_realization": (
                "the same multiset-matched count divided by planned exposure count"
            ),
            "mode_fraction_denominators": {
                "planned_realization": "matched mode frames / planned frames in that mode",
                "fraction_of_full_interval_acquired_frames_matching_mode": (
                    "matched mode frames / all strict raw frames in the full interval"
                ),
                "fraction_of_plan_active_acquired_frames_matching_mode": (
                    "matched mode frames / strict raw frames in the plan-active acquired cohort"
                ),
                "fraction_of_all_matched_frames_in_mode": (
                    "matched mode frames / all matched frames in the plan-active acquired cohort"
                ),
                "fraction_of_all_planned_frames_in_mode": (
                    "planned frames in mode / all archived plan rows in the interval"
                ),
            },
            "row_matching": (
                "within each (night,field_id), maximum-cardinality one-to-one monotonic "
                "minimum-total-absolute-time-residual matching; scheduler/sequence order "
                "fallback when timestamps are incomplete"
            ),
            "actual_time": (
                "raw obs_time_utc, then joined L2 obs_time_utc, then joined L2 MJD; "
                "file mtime is never used as exposure time"
            ),
            "actual_gap_proxy": (
                "successive actual start-time difference minus previous exposure duration; "
                "a wall-clock gap proxy, not isolated hardware overhead"
            ),
            "revisit_interval": (
                "successive start-time difference within the same night, field_id, and "
                "scheduler cycle"
            ),
            "revisit_baseline": (
                "one end-to-end start-time span per night, field_id, and scheduler cycle; "
                "only groups with at least two timed visits contribute"
            ),
            "n_sources": "L2 catalog row count, not asteroid discovery yield or efficiency",
        },
        "known_as_executed_limitations": {
            "near_sun_sorting": {
                "status": "known_code_bug",
                "code": KNOWN_LIMITATION_CODES[0],
                "detail": NEAR_SUN_STATUS,
                "source": "survey/scheduler.py:186-203",
            },
            "history_feedback": {
                "status": "known_operational_stall_from_server_audit",
                "code": KNOWN_LIMITATION_CODES[1],
                "detail": HISTORY_STATUS,
                "source": "survey/run_daily.py:137-154 and server history audit",
            },
        },
        "attribution_limit": ATTRIBUTION_LIMIT,
        "mode_rows": [json_safe_row(row) for row in mode_rows],
    }

    atomic_write_csv(
        output / "scheduler_plan_realization_by_night.csv", night_rows, NIGHT_OUTPUT_FIELDS
    )
    atomic_write_csv(output / "scheduler_exposure_matches.csv", matches, EXPOSURE_OUTPUT_FIELDS)
    atomic_write_csv(output / "scheduler_mode_summary.csv", mode_rows, MODE_OUTPUT_FIELDS)
    atomic_write_json(output / "scheduler_mode_summary.json", payload)
    atomic_write_csv(
        output / f"scheduler_example_{args.example_night}.csv",
        example_rows,
        EXPOSURE_OUTPUT_FIELDS,
    )
    atomic_write_csv(
        figure_data / "fig04_scheduler_plan_realization_by_night.csv",
        night_rows,
        FIGURE_NIGHT_FIELDS,
    )
    atomic_write_csv(
        figure_data / "fig04_scheduler_mode_summary.csv",
        mode_rows,
        MODE_OUTPUT_FIELDS,
    )
    atomic_write_csv(
        figure_data / f"fig04_scheduler_example_{args.example_night}.csv",
        example_rows,
        EXPOSURE_OUTPUT_FIELDS,
    )
    compliance = cohort_accounting["acquired_frame_plan_compliance"]
    compliance_text = "NA" if compliance is None else f"{compliance:.8f}"
    print(
        f"[done] plan={len(plans)} actual_full_interval={len(actuals)} "
        f"plan_active_acquired_denominator="
        f"{cohort_accounting['plan_active_acquired_frame_denominator_n']} "
        f"cohort_matched={cohort_accounting['plan_active_cohort_matched_frame_n']} "
        f"acquired_compliance={compliance_text} "
        f"output={output}",
        flush=True,
    )


def synthetic_plan(
    night: str,
    field_id: str,
    mode: str,
    order: int,
) -> dict[str, Any]:
    """Return a minimal internal plan row for the synthetic cohort test."""

    return {
        "_night": night,
        "_field_id": field_id,
        "_plan_order": order,
        "_planned_dt": None,
        "_planned_exptime_s": 30.0,
        "_mode": mode,
        "_cycle_id": "1",
        "_repeat": order,
        "_plan": {},
        "_plan_line": order + 2,
        "_planned_inter_start_s": None,
        "_planned_gap_proxy_s": None,
        "_planned_revisit_interval_s": None,
        "_planned_revisit_baseline_s": None,
    }


def synthetic_actual(
    night: str,
    field_id: str,
    order: int,
) -> dict[str, Any]:
    """Return a minimal internal actual row for the synthetic cohort test."""

    return {
        "_night": night,
        "_field_id": field_id,
        "_sequence_id": str(order),
        "_actual_order": order,
        "_sequence_order": order,
        "_actual_dt": None,
        "_actual_time_source": "",
        "_actual_exptime_s": 30.0,
        "_raw": {"file_name": f"fixture_{night}_{order}.fits", "path": "fixture"},
        "_raw_line": order + 2,
        "_l2": None,
        "_actual_inter_start_s": None,
        "_actual_gap_proxy_s": None,
    }


def run_self_test() -> None:
    """Exercise cohort closure and every published fraction denominator."""

    raw_only = "20990101"
    active_partial = "20990102"
    plan_only = "20990103"
    active_complete = "20990104"
    plans = [
        synthetic_plan(active_partial, "0010", "normal", 0),
        synthetic_plan(active_partial, "0010", "normal", 1),
        synthetic_plan(active_partial, "0020", "normal", 2),
        synthetic_plan(active_partial, "0030", "normal", 3),
        synthetic_plan(plan_only, "0040", "normal", 4),
        synthetic_plan(plan_only, "0050", "normal", 5),
        synthetic_plan(active_complete, "0060", "near_sun", 6),
    ]
    actuals = [
        synthetic_actual(raw_only, "0090", 0),
        synthetic_actual(raw_only, "0091", 1),
        synthetic_actual(raw_only, "0092", 2),
        synthetic_actual(active_partial, "0010", 3),
        synthetic_actual(active_partial, "0020", 4),
        synthetic_actual(active_partial, "0099", 5),
        synthetic_actual(active_complete, "0060", 6),
    ]
    matches = build_matches(plans, actuals, input_hash="fixture")
    cohort = build_cohort_accounting(plans, actuals, matches)
    assert cohort["strict_raw_frame_n_full_interval"] == 7
    assert cohort["plan_archive_available_night_n"] == 3
    assert cohort["plan_active_acquired_night_n"] == 2
    assert cohort["plan_active_acquired_frame_denominator_n"] == 4
    assert cohort["plan_active_cohort_matched_frame_n"] == 3
    assert cohort["plan_active_acquired_not_matched_n"] == 1
    assert cohort["full_interval_acquired_without_plan_archive_n"] == 3
    assert math.isclose(cohort["acquired_frame_plan_compliance"], 0.75)

    night_rows = build_night_summary(
        [raw_only, active_partial, plan_only, active_complete],
        plans,
        actuals,
        matches,
        night_status={},
        input_hash="fixture",
    )
    by_night = {row["night"]: row for row in night_rows}
    assert by_night[raw_only]["plan_archive_available"] is False
    assert by_night[raw_only]["acquired_frame_plan_compliance"] is None
    assert by_night[raw_only]["full_interval_acquired_without_plan_archive_n"] == 3
    assert by_night[active_partial]["plan_active_acquired_frame_denominator_n"] == 3
    assert by_night[active_partial]["plan_active_cohort_matched_frame_n"] == 2
    assert math.isclose(
        by_night[active_partial]["acquired_frame_plan_compliance"], 2 / 3
    )
    assert by_night[plan_only]["plan_archive_available"] is True
    assert by_night[plan_only]["plan_active_acquired_cohort"] is False
    assert by_night[plan_only]["acquired_frame_plan_compliance"] is None

    mode_rows = build_mode_summary(
        plans,
        actuals,
        matches,
        l2_fieldnames=set(),
        input_hash="fixture",
    )
    by_mode = {row["mode"]: row for row in mode_rows}
    normal = by_mode["normal"]
    assert normal["planned_n"] == 6
    assert normal["matched_n"] == 2
    assert normal["planned_realization_denominator_n"] == 6
    assert normal["full_interval_acquired_frame_denominator_n"] == 7
    assert normal["plan_active_acquired_frame_denominator_n"] == 4
    assert normal["all_matched_frame_denominator_n"] == 3
    assert normal["all_planned_frame_denominator_n"] == 7
    assert math.isclose(
        normal["fraction_of_full_interval_acquired_frames_matching_mode"], 2 / 7
    )
    assert math.isclose(
        normal["fraction_of_plan_active_acquired_frames_matching_mode"], 2 / 4
    )
    assert math.isclose(normal["fraction_of_all_matched_frames_in_mode"], 2 / 3)
    assert math.isclose(normal["fraction_of_all_planned_frames_in_mode"], 6 / 7)

    acquired_output_rows = [
        row
        for row in matches
        if row["match_status"] in {"matched", "acquired_not_planned"}
    ]
    assert sum(row["acquired_frame_in_plan_active_denominator"] for row in acquired_output_rows) == 4
    assert sum(row["acquired_frame_outside_plan_archive"] for row in acquired_output_rows) == 3
    assert not any(
        row["plan_archive_available"]
        for row in acquired_output_rows
        if row["night"] == raw_only
    )
    print("[self-test] analyze_scheduler cohort denominators: PASS")


def main() -> None:
    if "--self-test" in sys.argv[1:]:
        if sys.argv[1:] != ["--self-test"]:
            raise SystemExit("--self-test must be used alone")
        run_self_test()
        return
    parser = build_parser()
    args = parser.parse_args()
    try:
        run(args)
    except (FileExistsError, FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
