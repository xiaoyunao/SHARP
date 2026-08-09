#!/usr/bin/env python3
"""Generate Fig. 12 from frozen, evidence-backed operations products.

The figure deliberately keeps four concepts separate:

* automated pipeline processing;
* wall-clock human-review proxies;
* MPC submission-to-reply intervals; and
* evidence quality / operational exceptions.

Eligibility is evaluated independently for every metric.  A value enters its
plotted p50/p90/p95 statistics only when it is finite, nonnegative, belongs to
the metric's ``normal_daily`` run class, and its night has neither a
restart/reboot anomaly nor an explicit primary-science exclusion.  A negative
human-review mtime proxy therefore excludes only that human metric, not valid
automated metrics from the same night.  Historical reruns, latest-mtime-only
values, negative metrics, and anomaly nights remain in the exported audit data.
CPU and peak-RAM values are shown only when the operations summary actually
contains them; otherwise the figure says N/A.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from figure_styles import (
    STATISTICS_COLORS,
    STATISTICS_FOUR_PANEL_FIGSIZE,
    apply_statistics_style,
    save_pdf_png,
    style_statistics_axis,
)


AUTO_COLOR = STATISTICS_COLORS["blue"]
HUMAN_COLOR = STATISTICS_COLORS["orange"]
MPC_COLOR = STATISTICS_COLORS["purple"]
EXCEPTION_COLOR = STATISTICS_COLORS["red"]
RAW_COLOR = "#9aa3ad"
RESOURCE_COLOR = "#eceff2"


@dataclass(frozen=True)
class MetricSpec:
    name: str
    label: str
    group: str
    source_unit: str
    class_field: str


METRICS: tuple[MetricSpec, ...] = (
    MetricSpec(
        "auto_l2_latest_to_known_start_hours",
        "L2 latest → known start",
        "automatic",
        "hours",
        "known_run_class",
    ),
    MetricSpec(
        "auto_known_start_to_ades_hours",
        "Known start → ADES file",
        "automatic",
        "hours",
        "known_run_class",
    ),
    MetricSpec(
        "auto_l2_latest_to_known_ades_hours",
        "L2 latest → known ADES",
        "automatic",
        "hours",
        "known_run_class",
    ),
    MetricSpec(
        "auto_known_reply_to_unknown_start_hours",
        "Known reply → unknown start",
        "automatic",
        "hours",
        "unknown_run_class",
    ),
    MetricSpec(
        "auto_known_mask15_to_unknown_start_hours",
        "Known mask → unknown start",
        "automatic",
        "hours",
        "unknown_run_class",
    ),
    MetricSpec(
        "auto_unknown_start_to_end_minutes",
        "Unknown wrapper runtime",
        "automatic",
        "minutes",
        "unknown_run_class",
    ),
    MetricSpec(
        "auto_unknown_start_to_catalog_hours",
        "Unknown start → catalog",
        "automatic",
        "hours",
        "unknown_run_class",
    ),
    MetricSpec(
        "auto_unknown_start_to_review_package_hours",
        "Unknown start → review package",
        "automatic",
        "hours",
        "unknown_run_class",
    ),
    MetricSpec(
        "human_review_package_to_review_csv_hours",
        "Review package → review CSV",
        "human",
        "hours",
        "unknown_run_class",
    ),
    MetricSpec(
        "human_review_package_to_submit_csv_hours",
        "Review package → submit CSV",
        "human",
        "hours",
        "unknown_run_class",
    ),
    MetricSpec(
        "human_review_csv_to_submit_csv_hours",
        "Review CSV → submit CSV",
        "human",
        "hours",
        "unknown_run_class",
    ),
    MetricSpec(
        "mpc_known_ades_to_reply_hours",
        "Known ADES → MPC reply",
        "mpc",
        "hours",
        "known_run_class",
    ),
    MetricSpec(
        "mpc_unknown_ades_to_reply_hours",
        "Unknown ADES → MPC reply",
        "mpc",
        "hours",
        "unknown_run_class",
    ),
)


TIMELINE_FIELDS: tuple[tuple[str, str, str, str], ...] = (
    ("l2_latest_mtime_utc", "L2 latest", "automatic", "l2_latest_mtime"),
    ("known_start_utc", "Known start", "automatic", "known_start"),
    ("known_ades_mtime_utc", "Known ADES", "automatic", "known_ades_mtime"),
    ("known_reply_mtime_utc", "Known reply", "mpc", "known_reply_mtime"),
    ("unknown_start_utc", "Unknown start", "automatic", "unknown_start"),
    ("unknown_end_utc", "Unknown end", "automatic", "unknown_end"),
    (
        "unknown_catalog_mtime_utc",
        "Unknown catalog",
        "automatic",
        "unknown_catalog_mtime",
    ),
    (
        "review_manifest_mtime_utc",
        "Review pkg",
        "automatic",
        "review_manifest_mtime",
    ),
    ("review_csv_mtime_utc", "Review CSV", "human", "review_csv_mtime"),
    ("submit_csv_mtime_utc", "Submit selection", "human", "submit_csv_mtime"),
    ("unknown_ades_mtime_utc", "Unknown ADES", "human", "unknown_ades_mtime"),
    ("unknown_reply_mtime_utc", "Unknown reply", "mpc", "unknown_reply_mtime"),
)


GROUP_COLOR = {
    "automatic": AUTO_COLOR,
    "human": HUMAN_COLOR,
    "mpc": MPC_COLOR,
}
GROUP_LABEL = {
    "automatic": "Automated pipeline",
    "human": "Human review proxy",
    "mpc": "MPC exchange",
}


RESTART_PATTERN = re.compile(
    r"\b(?:reboot(?:ed|ing)?|restart(?:ed|ing)?|maintenance|power\s+cycle|"
    r"watcher\s+(?:down|failed|stopped)|daemon\s+(?:down|failed|stopped))\b",
    re.IGNORECASE,
)


def normalize_night(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.replace(r"\D", "", regex=True)
        .str.slice(0, 8)
        .str.zfill(8)
    )


def optional_bool(value: Any) -> bool | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "t"}:
        return True
    if text in {"0", "false", "no", "n", "f"}:
        return False
    return None


def finite_number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def quantile(values: Iterable[float], probability: float) -> float | None:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return None
    return float(np.percentile(array, probability * 100.0))


def distribution(values: Iterable[Any]) -> dict[str, Any]:
    nonnegative: list[float] = []
    negative = 0
    for value in values:
        number = finite_number(value)
        if number is None:
            continue
        if number < 0:
            negative += 1
        else:
            nonnegative.append(number)
    nonnegative.sort()
    return {
        "n": len(nonnegative),
        "n_negative_excluded": negative,
        "p50": quantile(nonnegative, 0.50),
        "p90": quantile(nonnegative, 0.90),
        "p95": quantile(nonnegative, 0.95),
    }


def close_or_both_none(left: Any, right: Any, tolerance: float = 1e-9) -> bool:
    if left is None and right is None:
        return True
    if left is None or right is None:
        return False
    try:
        return math.isclose(float(left), float(right), rel_tol=tolerance, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def validate_summary_against_nightly(
    nightly: pd.DataFrame,
    summary: dict[str, Any],
) -> None:
    """Ensure the supplied summary is derived from the supplied nightly table."""

    categories = summary.get("latency_categories")
    if not isinstance(categories, dict):
        raise ValueError("operations_latency_summary.json lacks latency_categories")
    mismatches: list[str] = []
    for spec in METRICS:
        if spec.name not in nightly.columns:
            raise ValueError(f"pipeline_latency_by_night.csv lacks {spec.name}")
        if spec.class_field not in nightly.columns:
            raise ValueError(f"pipeline_latency_by_night.csv lacks {spec.class_field}")
        raw = nightly.loc[
            nightly[spec.class_field].astype("string").eq("normal_daily"), spec.name
        ]
        recomputed = distribution(raw)
        supplied = (
            categories.get(spec.group, {})
            .get(spec.name, {})
            .get("normal_daily")
        )
        if not isinstance(supplied, dict):
            mismatches.append(f"{spec.name}: missing normal_daily summary")
            continue
        comparisons = {
            "n": (recomputed["n"], supplied.get("n")),
            "n_negative_excluded": (
                recomputed["n_negative_excluded"],
                supplied.get("n_negative_excluded"),
            ),
            "p50": (recomputed["p50"], supplied.get("median")),
            "p90": (recomputed["p90"], supplied.get("p90")),
            "p95": (recomputed["p95"], supplied.get("p95")),
        }
        for name, (observed, expected) in comparisons.items():
            if name in {"n", "n_negative_excluded"}:
                if observed != expected:
                    mismatches.append(
                        f"{spec.name}.{name}: nightly={observed!r}, summary={expected!r}"
                    )
            elif not close_or_both_none(observed, expected):
                mismatches.append(
                    f"{spec.name}.{name}: nightly={observed!r}, summary={expected!r}"
                )
    if mismatches:
        preview = "; ".join(mismatches[:8])
        if len(mismatches) > 8:
            preview += f"; ... ({len(mismatches)} mismatches total)"
        raise ValueError(f"nightly/summary latency mismatch: {preview}")


def load_inputs(
    latency_path: Path,
    summary_path: Path,
    evidence_path: Path,
    night_status_path: Path | None,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    for path in (latency_path, summary_path, evidence_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    nightly = pd.read_csv(latency_path, dtype={"night": "string"})
    if "night" not in nightly:
        raise ValueError("pipeline latency table lacks night")
    nightly["night"] = normalize_night(nightly["night"])
    if nightly["night"].isna().any() or nightly["night"].duplicated().any():
        raise ValueError("pipeline latency table requires one complete row per night")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    evidence = pd.read_csv(evidence_path, dtype={"night": "string"})
    required_evidence = {"night", "event_type", "timestamp_utc", "source_kind"}
    missing_evidence = required_evidence.difference(evidence.columns)
    if missing_evidence:
        raise ValueError(
            f"operations event evidence lacks columns: {sorted(missing_evidence)}"
        )
    evidence["night"] = normalize_night(evidence["night"])
    evidence["timestamp_parsed_utc"] = pd.to_datetime(
        evidence["timestamp_utc"], errors="coerce", utc=True
    )

    if night_status_path is not None:
        if not night_status_path.is_file():
            raise FileNotFoundError(night_status_path)
        status = pd.read_csv(night_status_path, dtype={"night": "string"})
        if "night" not in status:
            raise ValueError("night_status table lacks night")
        status["night"] = normalize_night(status["night"])
        if status["night"].duplicated().any():
            raise ValueError("night_status table contains duplicate nights")
        fields = [
            name
            for name in (
                "quality_code",
                "quality_reason",
                "primary_science_included",
                "unknown_science_included",
                "known_mpc_state",
                "unknown_mpc_state",
            )
            if name in status.columns
        ]
        merged = nightly.merge(
            status[["night"] + fields],
            on="night",
            how="left",
            validate="one_to_one",
            suffixes=("", "_night_status"),
        )
        for field in fields:
            incoming = f"{field}_night_status"
            if incoming not in merged:
                continue
            if field in merged:
                present = merged[incoming].notna() & merged[incoming].astype(str).str.strip().ne("")
                merged.loc[present, field] = merged.loc[present, incoming]
                merged = merged.drop(columns=[incoming])
            else:
                merged = merged.rename(columns={incoming: field})
        nightly = merged

    required_nightly = {
        "known_run_class",
        "unknown_run_class",
        "operations_run_class",
        "timing_consistency_flags",
        "cpu_accounting_status",
        "ram_accounting_status",
    }
    required_nightly.update(spec.name for spec in METRICS)
    missing_nightly = required_nightly.difference(nightly.columns)
    if missing_nightly:
        raise ValueError(
            f"pipeline latency table lacks columns: {sorted(missing_nightly)}"
        )
    validate_summary_against_nightly(nightly, summary)
    return nightly, summary, evidence


def find_restart_nights(
    nightly: pd.DataFrame,
    evidence: pd.DataFrame,
) -> tuple[set[str], pd.DataFrame]:
    searchable_columns = [
        name
        for name in (
            "event_type",
            "line_text",
            "outcome",
            "outcome_tags",
            "selection_role",
            "source_relative_path",
        )
        if name in evidence.columns
    ]
    if searchable_columns:
        search_text = evidence[searchable_columns].fillna("").astype(str).agg(" | ".join, axis=1)
        event_mask = search_text.str.contains(RESTART_PATTERN, na=False)
    else:
        event_mask = pd.Series(False, index=evidence.index)
    restart_events = evidence.loc[event_mask].copy()
    nights = set(restart_events["night"].dropna().astype(str))
    status_columns = [
        name for name in ("quality_code", "quality_reason") if name in nightly.columns
    ]
    if status_columns:
        status_text = nightly[status_columns].fillna("").astype(str).agg(" | ".join, axis=1)
        status_mask = status_text.str.contains(RESTART_PATTERN, na=False)
        nights.update(nightly.loc[status_mask, "night"].dropna().astype(str))
    return nights, restart_events


def build_latency_values(
    nightly: pd.DataFrame,
    restart_nights: set[str],
) -> pd.DataFrame:
    numeric = nightly.copy()
    for spec in METRICS:
        numeric[spec.name] = pd.to_numeric(numeric[spec.name], errors="coerce")
    any_negative = numeric[[spec.name for spec in METRICS]].lt(0).any(axis=1)
    consistency_negative = (
        numeric["timing_consistency_flags"]
        .fillna("")
        .astype(str)
        .str.contains("_negative", regex=False)
    )
    negative_nights = set(numeric.loc[any_negative | consistency_negative, "night"])

    rows: list[dict[str, Any]] = []
    for row in numeric.to_dict(orient="records"):
        night = str(row["night"])
        primary = optional_bool(row.get("primary_science_included"))
        for spec in METRICS:
            original = finite_number(row.get(spec.name))
            run_class = str(row.get(spec.class_field, "") or "")
            value_hours = (
                original / 60.0
                if original is not None and spec.source_unit == "minutes"
                else original
            )
            metric_negative = original is not None and original < 0
            reasons: list[str] = []
            if original is None:
                reasons.append("missing_value")
            if metric_negative:
                reasons.append("negative_latency")
            if run_class != "normal_daily":
                reasons.append(run_class or "run_class_missing")
            if night in restart_nights:
                reasons.append("restart_or_reboot_anomaly")
            if primary is False:
                reasons.append("primary_science_excluded")
            usable = (
                original is not None
                and original >= 0
                and run_class == "normal_daily"
                and night not in restart_nights
                and primary is not False
            )
            rows.append(
                {
                    "night": night,
                    "group": spec.group,
                    "metric": spec.name,
                    "label": spec.label,
                    "source_unit": spec.source_unit,
                    "display_unit": "hours",
                    "value_original": original,
                    "value_hours": value_hours,
                    "class_field": spec.class_field,
                    "run_class": run_class,
                    "usable_normal": usable,
                    "exclusion_reason": ";".join(reasons),
                    "metric_negative": metric_negative,
                    "night_has_any_negative_metric": night in negative_nights,
                    "restart_anomaly_night": night in restart_nights,
                    "primary_science_included": primary,
                    "quality_code": row.get("quality_code", ""),
                    "timing_consistency_flags": row.get(
                        "timing_consistency_flags", ""
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_latency_statistics(values: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for spec in METRICS:
        metric = values.loc[values["metric"].eq(spec.name)]
        eligible = pd.to_numeric(
            metric.loc[metric["usable_normal"], "value_hours"], errors="coerce"
        ).dropna()
        rows.append(
            {
                "group": spec.group,
                "metric": spec.name,
                "label": spec.label,
                "display_unit": "hours",
                "usable_normal_n": int(len(eligible)),
                "p50_hours": quantile(eligible, 0.50),
                "p90_hours": quantile(eligible, 0.90),
                "p95_hours": quantile(eligible, 0.95),
                "available_value_n": int(metric["value_hours"].notna().sum()),
                "negative_value_n": int(
                    pd.to_numeric(metric["value_original"], errors="coerce").lt(0).sum()
                ),
                "historical_rerun_value_n": int(
                    (metric["run_class"].eq("historical_rerun") & metric["value_hours"].notna()).sum()
                ),
                "latest_mtime_only_value_n": int(
                    (metric["run_class"].eq("latest_mtime_only") & metric["value_hours"].notna()).sum()
                ),
                "restart_anomaly_value_n": int(
                    (metric["restart_anomaly_night"] & metric["value_hours"].notna()).sum()
                ),
                "percentile_filter": (
                    "per-metric finite_nonnegative; metric run_class=normal_daily; "
                    "exclude restart/reboot anomaly or explicit primary-science exclusion; "
                    "a negative different metric does not exclude this metric"
                ),
            }
        )
    return pd.DataFrame(rows)


def lookup_evidence(
    evidence: pd.DataFrame,
    night: str,
    event_type: str,
    timestamp: pd.Timestamp,
) -> dict[str, Any]:
    candidates = evidence.loc[
        evidence["night"].eq(night) & evidence["event_type"].astype(str).eq(event_type)
    ].copy()
    if candidates.empty:
        return {}
    valid = candidates["timestamp_parsed_utc"].notna()
    if valid.any():
        delta = (
            candidates.loc[valid, "timestamp_parsed_utc"] - timestamp
        ).abs()
        index = delta.idxmin()
    else:
        index = candidates.index[0]
    row = candidates.loc[index]
    return {
        key: row.get(key, "")
        for key in (
            "event_id",
            "source_kind",
            "source_group",
            "source_relative_path",
            "selected_for_night",
            "selection_role",
            "run_class",
            "outcome",
        )
    }


def choose_example_night(
    nightly: pd.DataFrame,
    values: pd.DataFrame,
    restart_nights: set[str],
) -> str | None:
    timeline_available = []
    for row in nightly.to_dict(orient="records"):
        night = str(row["night"])
        night_values = values.loc[values["night"].eq(night)]
        timestamp_n = sum(
            pd.notna(pd.to_datetime(row.get(field), errors="coerce", utc=True))
            for field, _label, _group, _event in TIMELINE_FIELDS
        )
        group_n = int(
            night_values.loc[night_values["usable_normal"], "group"].nunique()
        )
        usable_metric_n = int(night_values["usable_normal"].sum())
        negative_metric_n = int(night_values["metric_negative"].sum())
        normal_n = sum(
            str(row.get(field, "")) == "normal_daily"
            for field in ("known_run_class", "unknown_run_class")
        )
        primary_excluded = optional_bool(row.get("primary_science_included")) is False
        score = timestamp_n + 5 * group_n + 3 * normal_n + usable_metric_n
        usable_candidate = (
            normal_n > 0
            and usable_metric_n > 0
            and night not in restart_nights
            and not primary_excluded
            and timestamp_n >= 3
        )
        timeline_available.append(
            (usable_candidate, score, -negative_metric_n, night)
        )
    eligible = [item for item in timeline_available if item[0]]
    if not eligible:
        return None
    return max(eligible, key=lambda item: (item[1], item[2], item[3]))[3]


def build_example_timeline(
    nightly: pd.DataFrame,
    evidence: pd.DataFrame,
    values: pd.DataFrame,
    restart_nights: set[str],
) -> pd.DataFrame:
    night = choose_example_night(nightly, values, restart_nights)
    columns = [
        "night",
        "event_order",
        "field",
        "event_type",
        "label",
        "group",
        "timestamp_utc",
        "elapsed_hours",
        "elapsed_log1p",
        "event_id",
        "source_kind",
        "source_group",
        "source_relative_path",
        "selected_for_night",
        "selection_role",
        "run_class",
        "outcome",
    ]
    if night is None:
        return pd.DataFrame(columns=columns)
    row = nightly.loc[nightly["night"].eq(night)].iloc[0]
    events: list[dict[str, Any]] = []
    for field, label, group, event_type in TIMELINE_FIELDS:
        timestamp = pd.to_datetime(row.get(field), errors="coerce", utc=True)
        if pd.isna(timestamp):
            continue
        event = {
            "night": night,
            "field": field,
            "event_type": event_type,
            "label": label,
            "group": group,
            "timestamp_utc": timestamp,
        }
        event.update(lookup_evidence(evidence, night, event_type, timestamp))
        events.append(event)
    if not events:
        return pd.DataFrame(columns=columns)
    events.sort(key=lambda item: (item["timestamp_utc"], item["label"]))
    first = events[0]["timestamp_utc"]
    for order, event in enumerate(events, start=1):
        elapsed = (event["timestamp_utc"] - first).total_seconds() / 3600.0
        event["event_order"] = order
        event["elapsed_hours"] = elapsed
        event["elapsed_log1p"] = math.log1p(max(elapsed, 0.0))
        event["timestamp_utc"] = event["timestamp_utc"].isoformat()
    return pd.DataFrame(events).reindex(columns=columns)


def event_stage_nights(evidence: pd.DataFrame) -> dict[str, set[str]]:
    stages: dict[str, set[str]] = {group: set() for group in GROUP_LABEL}
    automatic_pattern = re.compile(
        r"^(?:l2_|known_start|unknown_start|unknown_end|known_(?:all|match1|mask15|ades)_mtime|"
        r"unknown_(?:catalog|summary)_mtime|review_manifest_mtime)"
    )
    human_pattern = re.compile(
        r"^(?:review_csv_mtime|submit_csv_mtime|unknown_ades_mtime)"
    )
    mpc_pattern = re.compile(r"^(?:known_reply_mtime|unknown_reply_mtime)")
    for row in evidence[["night", "event_type"]].dropna().itertuples(index=False):
        event_type = str(row.event_type)
        night = str(row.night)
        if automatic_pattern.search(event_type):
            stages["automatic"].add(night)
        if human_pattern.search(event_type):
            stages["human"].add(night)
        if mpc_pattern.search(event_type):
            stages["mpc"].add(night)
    return stages


def resource_status(summary: dict[str, Any], key: str) -> tuple[str, Any]:
    record = summary.get("resource_accounting", {}).get(key, {})
    if not isinstance(record, dict):
        return "unavailable_not_reported", None
    return str(record.get("status") or "unavailable_not_reported"), record.get("value")


def build_evidence_matrix(
    nightly: pd.DataFrame,
    values: pd.DataFrame,
    evidence: pd.DataFrame,
    restart_nights: set[str],
    summary: dict[str, Any],
) -> pd.DataFrame:
    stage_events = event_stage_nights(evidence)
    row_specs = (
        ("Usable normal", "usable_normal"),
        ("Historical rerun", "historical_rerun"),
        ("Latest mtime only", "latest_mtime_only"),
        ("Negative latency", "negative_latency"),
        ("Restart anomaly", "restart_anomaly"),
        ("No evidence", "no_evidence"),
    )
    records: list[dict[str, Any]] = []
    for row_label, status in row_specs:
        for group in ("automatic", "human", "mpc"):
            subset = values.loc[values["group"].eq(group)]
            if status == "usable_normal":
                nights = set(subset.loc[subset["usable_normal"], "night"])
            elif status in {"historical_rerun", "latest_mtime_only", "no_evidence"}:
                nights = set(subset.loc[subset["run_class"].eq(status), "night"])
            elif status == "negative_latency":
                nights = set(
                    subset.loc[
                        subset["metric_negative"],
                        "night",
                    ]
                )
            else:
                stage_nights = set(
                    subset.loc[subset["value_hours"].notna(), "night"]
                ).union(stage_events[group])
                nights = restart_nights.intersection(stage_nights)
            records.append(
                {
                    "status_row": row_label,
                    "status_code": status,
                    "stage": group,
                    "stage_label": GROUP_LABEL[group],
                    "night_count": len(nights),
                    "night_keys": ";".join(sorted(nights)),
                    "resource_status": "",
                    "resource_value": np.nan,
                }
            )
        for resource in ("cpu", "ram"):
            status_text, value = resource_status(summary, resource)
            records.append(
                {
                    "status_row": row_label,
                    "status_code": status,
                    "stage": resource,
                    "stage_label": resource.upper(),
                    "night_count": np.nan,
                    "night_keys": "",
                    "resource_status": status_text,
                    "resource_value": value,
                }
            )
    return pd.DataFrame(records)


def build_exception_nights(
    nightly: pd.DataFrame,
    values: pd.DataFrame,
    restart_nights: set[str],
    restart_events: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    restart_counts = (
        restart_events.groupby("night").size().to_dict()
        if not restart_events.empty
        else {}
    )
    for row in nightly.to_dict(orient="records"):
        night = str(row["night"])
        subset = values.loc[values["night"].eq(night)]
        negative_metrics = sorted(
            subset.loc[
                pd.to_numeric(subset["value_original"], errors="coerce").lt(0),
                "metric",
            ].unique()
        )
        flags: list[str] = []
        for component, field in (
            ("known", "known_run_class"),
            ("unknown", "unknown_run_class"),
        ):
            run_class = str(row.get(field, "") or "")
            if run_class in {
                "historical_rerun",
                "latest_mtime_only",
                "no_evidence",
                "nonstandard_date_relation",
            }:
                flags.append(f"{component}_{run_class}")
        if negative_metrics:
            flags.append("negative_timing")
        if night in restart_nights:
            flags.append("restart_or_reboot_anomaly")
        if optional_bool(row.get("primary_science_included")) is False:
            flags.append("primary_science_excluded")
        if not flags:
            continue
        rows.append(
            {
                "night": night,
                "exception_flags": ";".join(flags),
                "known_run_class": row.get("known_run_class", ""),
                "unknown_run_class": row.get("unknown_run_class", ""),
                "operations_run_class": row.get("operations_run_class", ""),
                "negative_metrics": ";".join(negative_metrics),
                "timing_consistency_flags": row.get("timing_consistency_flags", ""),
                "restart_event_count": int(restart_counts.get(night, 0)),
                "quality_code": row.get("quality_code", ""),
                "quality_reason": row.get("quality_reason", ""),
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "night",
            "exception_flags",
            "known_run_class",
            "unknown_run_class",
            "operations_run_class",
            "negative_metrics",
            "timing_consistency_flags",
            "restart_event_count",
            "quality_code",
            "quality_reason",
        ],
    )


def elapsed_ticks(maximum: float) -> tuple[list[float], list[str]]:
    candidates = [0, 0.25, 1, 3, 6, 12, 24, 72, 168, 336, 720, 1440]
    chosen = [value for value in candidates if value <= maximum * 1.02]
    if not chosen or chosen[-1] < maximum * 0.65:
        chosen.append(maximum)
    positions = [math.log1p(max(value, 0.0)) for value in chosen]
    labels = [
        f"{value:g}"
        if value < 24
        else (f"{value / 24:g} d" if value % 24 == 0 else f"{value:g}")
        for value in chosen
    ]
    return positions, labels


def plot_timeline(axis, timeline: pd.DataFrame) -> None:
    if timeline.empty:
        axis.set_title(
            "(a) Actual event sequence example",
            loc="left",
            fontsize=20,
            weight="bold",
        )
        axis.text(
            0.5,
            0.54,
            "No normal-daily night has enough\ntraceable timestamps for an event chain.",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=18,
            color=STATISTICS_COLORS["mid_grey"],
        )
        axis.text(
            0.5,
            0.35,
            "Historical reruns and mtime-only evidence are not used as examples.",
            transform=axis.transAxes,
            ha="center",
            fontsize=13,
            color=EXCEPTION_COLOR,
        )
        axis.set_axis_off()
        return
    example_night = str(timeline.iloc[0]["night"])
    axis.set_title(
        f"(a) Actual event sequence example — {example_night}",
        loc="left",
        fontsize=20,
        weight="bold",
    )
    lanes = {"automatic": 2, "human": 1, "mpc": 0}
    maximum = max(float(timeline["elapsed_hours"].max()), 0.1)
    for group, y in lanes.items():
        axis.hlines(
            y,
            0,
            math.log1p(maximum),
            color="#d4d9de",
            linewidth=1.0,
            zorder=0,
        )
    ordered = timeline.sort_values("event_order")
    sequence_x = ordered["elapsed_log1p"].to_numpy(dtype=float)
    sequence_y = ordered["group"].map(lanes).to_numpy(dtype=float)
    axis.plot(
        sequence_x,
        sequence_y,
        color=STATISTICS_COLORS["mid_grey"],
        linewidth=1.0,
        alpha=0.55,
        zorder=1,
    )
    marker = {"automatic": "o", "human": "s", "mpc": "D"}
    for group, frame in ordered.groupby("group", sort=False):
        axis.scatter(
            frame["elapsed_log1p"],
            frame["group"].map(lanes),
            s=68,
            marker=marker[group],
            color=GROUP_COLOR[group],
            edgecolor="white",
            linewidth=0.7,
            zorder=4,
            label=GROUP_LABEL[group],
        )
    # Dense wrapper events can be only seconds apart and share the same log1p
    # x-position.  Keep their points/evidence in the timeline but label only
    # representative milestones so the exported panel remains legible.
    omitted_labels = {
        "Known start",
        "Known ADES",
        "Known reply",
        "Unknown start",
        "Unknown end",
        "Unknown catalog",
        "Review pkg",
        "Submit selection",
        "Unknown ADES",
    }
    labelled_by_group = {group: 0 for group in lanes}
    maximum_x = max(float(ordered["elapsed_log1p"].max()), 0.1)
    for _index, row in ordered.reset_index(drop=True).iterrows():
        if row["label"] in omitted_labels:
            continue
        y = lanes[row["group"]]
        group_index = labelled_by_group[row["group"]]
        labelled_by_group[row["group"]] += 1
        if row["group"] == "automatic":
            vertical = (14, 31)[group_index % 2]
        elif row["group"] == "human":
            vertical = (-18, 14)[group_index % 2]
        else:
            vertical = (-18, -35)[group_index % 2]
        horizontal_alignment = (
            "right"
            if float(row["elapsed_log1p"]) > 0.78 * maximum_x
            else "left"
        )
        axis.annotate(
            row["label"],
            (row["elapsed_log1p"], y),
            xytext=(0, vertical),
            textcoords="offset points",
            ha=horizontal_alignment,
            va="bottom" if vertical > 0 else "top",
            rotation=48,
            fontsize=9.5,
            color=STATISTICS_COLORS["ink"],
        )
    tick_positions, tick_labels = elapsed_ticks(maximum)
    axis.set_xticks(tick_positions, tick_labels)
    axis.set_xlim(-0.03, math.log1p(maximum) * 1.04 + 0.02)
    axis.set_ylim(-0.55, 2.55)
    axis.set_yticks(
        [2, 1, 0],
        [GROUP_LABEL["automatic"], GROUP_LABEL["human"], GROUP_LABEL["mpc"]],
        fontsize=13,
    )
    axis.set_xlabel(
        "Elapsed time since first event (hours; log1p spacing)",
        fontsize=15,
    )
    axis.text(
        0.99,
        0.03,
        "UTC evidence; unlabeled points are intermediate persisted events",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=10.5,
        color=STATISTICS_COLORS["mid_grey"],
    )
    axis.grid(axis="x", alpha=0.16, linewidth=0.6)
    axis.grid(axis="y", visible=False)
    axis.tick_params(axis="x", labelsize=12)
    axis.spines[["top", "right"]].set_visible(False)


def percentile_handles() -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=STATISTICS_COLORS["ink"],
            markeredgecolor=STATISTICS_COLORS["ink"],
            markersize=7,
            label="p50",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="none",
            markerfacecolor="white",
            markeredgecolor=STATISTICS_COLORS["ink"],
            markersize=6,
            label="p90",
        ),
        Line2D(
            [0],
            [0],
            marker="^",
            color="none",
            markerfacecolor="white",
            markeredgecolor=STATISTICS_COLORS["ink"],
            markersize=7,
            label="p95",
        ),
        Line2D(
            [0],
            [0],
            marker=".",
            color="none",
            markerfacecolor=RAW_COLOR,
            markeredgecolor=RAW_COLOR,
            markersize=6,
            alpha=0.55,
            label="eligible nights",
        ),
    ]


def plot_latency_panel(
    axis,
    values: pd.DataFrame,
    statistics: pd.DataFrame,
    groups: tuple[str, ...],
    title: str,
) -> None:
    specs = [spec for spec in METRICS if spec.group in groups]
    y_positions = np.arange(len(specs))[::-1]
    maximum = 0.0
    for y, spec in zip(y_positions, specs):
        metric_values = values.loc[
            values["metric"].eq(spec.name) & values["usable_normal"],
            "value_hours",
        ]
        numeric = pd.to_numeric(metric_values, errors="coerce").dropna().to_numpy()
        row = statistics.loc[statistics["metric"].eq(spec.name)].iloc[0]
        color = GROUP_COLOR[spec.group]
        if numeric.size:
            offsets = (
                np.linspace(-0.16, 0.16, numeric.size)
                if numeric.size > 1
                else np.array([0.0])
            )
            axis.scatter(
                numeric,
                y + offsets,
                s=10,
                color=RAW_COLOR,
                alpha=0.34,
                linewidths=0,
                rasterized=True,
                zorder=1,
            )
            p50 = float(row["p50_hours"])
            p90 = float(row["p90_hours"])
            p95 = float(row["p95_hours"])
            axis.hlines(y, p50, p95, color=color, linewidth=3.0, alpha=0.85, zorder=3)
            axis.scatter([p50], [y], s=64, color=color, edgecolor="white", linewidth=0.7, zorder=5)
            axis.scatter([p90], [y], s=48, marker="D", facecolor="white", edgecolor=color, linewidth=1.8, zorder=5)
            axis.scatter([p95], [y], s=62, marker="^", facecolor="white", edgecolor=color, linewidth=1.8, zorder=5)
            maximum = max(maximum, float(np.nanmax(numeric)), p95)
        else:
            axis.text(
                0.98,
                y + 0.10,
                "no usable normal nights",
                transform=axis.get_yaxis_transform(),
                ha="right",
                va="bottom",
                fontsize=10,
                color=EXCEPTION_COLOR,
            )
    labels = []
    for spec in specs:
        n = int(
            statistics.loc[
                statistics["metric"].eq(spec.name), "usable_normal_n"
            ].iloc[0]
        )
        labels.append(f"{spec.label}  (N={n})")
    axis.set_yticks(y_positions, labels, fontsize=11.5)
    axis.set_xscale("symlog", linthresh=0.1, linscale=1.0)
    axis.set_xlim(left=0, right=max(1.0, maximum * 1.18))
    axis.set_xlabel(
        "Latency (hours; symlog scale below 0.1 h)",
        fontsize=16,
    )
    axis.set_title(title, loc="left", fontsize=20, weight="bold")
    style_statistics_axis(axis, tick_fontsize=13)
    axis.grid(axis="y", visible=False)
    axis.axvline(0, color=STATISTICS_COLORS["ink"], linewidth=0.8, alpha=0.45)
    handles = percentile_handles()
    if len(groups) > 1:
        handles.extend(
            [
                Patch(facecolor=GROUP_COLOR[group], edgecolor="none", label=GROUP_LABEL[group])
                for group in groups
            ]
        )
    axis.legend(handles=handles, fontsize=11, loc="upper left", ncol=2)


def concise_resource(status: str, value: Any) -> str:
    if value is not None and str(value).strip() not in {"", "None", "nan"}:
        return f"{status.replace('_', ' ')}: {value}"
    return status.replace("_", " ")


def plot_evidence_matrix(
    axis,
    matrix: pd.DataFrame,
    summary: dict[str, Any],
) -> None:
    row_order = [
        "Usable normal",
        "Historical rerun",
        "Latest mtime only",
        "Negative latency",
        "Restart anomaly",
        "No evidence",
    ]
    column_order = ["automatic", "human", "mpc", "cpu", "ram"]
    array = np.full((len(row_order), len(column_order)), np.nan)
    for row_index, row_label in enumerate(row_order):
        for column_index, stage in enumerate(column_order):
            record = matrix.loc[
                matrix["status_row"].eq(row_label) & matrix["stage"].eq(stage)
            ]
            if not record.empty and pd.notna(record.iloc[0]["night_count"]):
                array[row_index, column_index] = float(record.iloc[0]["night_count"])
    masked = np.ma.masked_invalid(array)
    maximum = max(1.0, float(np.nanmax(array[:, :3])) if np.isfinite(array[:, :3]).any() else 1.0)
    axis.imshow(masked, cmap="Blues", vmin=0, vmax=maximum, aspect="auto")
    for row_index in range(array.shape[0]):
        for column_index in range(array.shape[1]):
            value = array[row_index, column_index]
            if np.isfinite(value):
                color = "white" if value > maximum * 0.58 else STATISTICS_COLORS["ink"]
                text = f"{int(value)}"
            else:
                color = STATISTICS_COLORS["mid_grey"]
                text = "N/A"
                axis.add_patch(
                    plt.Rectangle(
                        (column_index - 0.5, row_index - 0.5),
                        1,
                        1,
                        facecolor=RESOURCE_COLOR,
                        edgecolor="white",
                        linewidth=1.0,
                    )
                )
            axis.text(
                column_index,
                row_index,
                text,
                ha="center",
                va="center",
                fontsize=15,
                color=color,
                weight="bold",
            )
    axis.set_xticks(
        np.arange(len(column_order)),
        ["Automated", "Human", "MPC", "CPU", "RAM"],
        fontsize=12,
    )
    axis.set_yticks(np.arange(len(row_order)), row_order, fontsize=12)
    axis.set_xticks(np.arange(-0.5, len(column_order), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(row_order), 1), minor=True)
    axis.grid(which="minor", color="white", linestyle="-", linewidth=1.2)
    axis.tick_params(which="minor", bottom=False, left=False)
    axis.set_title(
        "(d) Evidence and exception status (distinct nights)",
        loc="left",
        fontsize=20,
        weight="bold",
    )
    cpu_status, cpu_value = resource_status(summary, "cpu")
    ram_status, ram_value = resource_status(summary, "ram")
    axis.text(
        0.0,
        -0.13,
        f"CPU: {concise_resource(cpu_status, cpu_value)}  |  "
        f"RAM: {concise_resource(ram_status, ram_value)}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        color=STATISTICS_COLORS["mid_grey"],
    )
    axis.text(
        0.0,
        -0.21,
        "Exception rows may overlap; CPU/RAM cells are N/A, never zero-filled.",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        color=EXCEPTION_COLOR,
    )


def write_figure_data(
    directory: Path,
    values: pd.DataFrame,
    statistics: pd.DataFrame,
    timeline: pd.DataFrame,
    matrix: pd.DataFrame,
    exceptions: pd.DataFrame,
) -> list[Path]:
    directory.mkdir(parents=True, exist_ok=True)
    outputs = [
        directory / "fig12_latency_values.csv",
        directory / "fig12_latency_statistics.csv",
        directory / "fig12_example_timeline.csv",
        directory / "fig12_evidence_status_matrix.csv",
        directory / "fig12_exception_nights.csv",
    ]
    for frame, path in zip(
        (values, statistics, timeline, matrix, exceptions), outputs
    ):
        frame.to_csv(path, index=False)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate the evidence-separated PASP operations timeline figure."
    )
    parser.add_argument("--latency-by-night", type=Path, required=True)
    parser.add_argument("--latency-summary", type=Path, required=True)
    parser.add_argument("--event-evidence", type=Path, required=True)
    parser.add_argument("--night-status", type=Path)
    parser.add_argument("--output-base", type=Path, required=True)
    parser.add_argument("--figure-data-dir", type=Path, required=True)
    args = parser.parse_args()

    nightly, summary, evidence = load_inputs(
        args.latency_by_night.resolve(),
        args.latency_summary.resolve(),
        args.event_evidence.resolve(),
        args.night_status.resolve() if args.night_status else None,
    )
    restart_nights, restart_events = find_restart_nights(nightly, evidence)
    values = build_latency_values(nightly, restart_nights)
    statistics = build_latency_statistics(values)
    timeline = build_example_timeline(
        nightly, evidence, values, restart_nights
    )
    matrix = build_evidence_matrix(
        nightly, values, evidence, restart_nights, summary
    )
    exceptions = build_exception_nights(
        nightly, values, restart_nights, restart_events
    )
    figure_data_paths = write_figure_data(
        args.figure_data_dir.resolve(),
        values,
        statistics,
        timeline,
        matrix,
        exceptions,
    )

    apply_statistics_style()
    figure, axes = plt.subplots(
        2,
        2,
        figsize=STATISTICS_FOUR_PANEL_FIGSIZE,
        gridspec_kw={"height_ratios": [0.93, 1.07]},
    )
    plot_timeline(axes[0, 0], timeline)
    plot_latency_panel(
        axes[0, 1],
        values,
        statistics,
        ("automatic",),
        "(b) Automated latency: usable normal nights",
    )
    plot_latency_panel(
        axes[1, 0],
        values,
        statistics,
        ("human", "mpc"),
        "(c) Human-review proxies and MPC round trips",
    )
    plot_evidence_matrix(axes[1, 1], matrix, summary)
    figure.text(
        0.5,
        0.014,
        (
            "p50/p90/p95 eligibility is per metric: finite nonnegative, normal-daily, no "
            "restart/reboot or explicit primary-quality exclusion. A negative human proxy "
            "does not remove valid automated metrics from the same night. Historical reruns "
            "and latest-mtime-only values remain audit context. Human intervals are wall-clock mtime proxies; "
            "inventory product timestamps are latest mtimes, not creation times; MPC intervals "
            "include transport and external response time."
        ),
        ha="center",
        va="bottom",
        fontsize=13.5,
        color=STATISTICS_COLORS["mid_grey"],
        wrap=True,
    )
    figure.tight_layout(rect=[0.02, 0.065, 0.99, 0.99], h_pad=2.8, w_pad=2.0)
    png_path, pdf_path = save_pdf_png(figure, args.output_base.resolve())

    cpu_status, cpu_value = resource_status(summary, "cpu")
    ram_status, ram_value = resource_status(summary, "ram")
    print(
        json.dumps(
            {
                "pdf": str(pdf_path),
                "png": str(png_path),
                "figure_data": [str(path) for path in figure_data_paths],
                "example_night": (
                    str(timeline.iloc[0]["night"]) if not timeline.empty else None
                ),
                "usable_normal_counts": {
                    row.metric: int(row.usable_normal_n)
                    for row in statistics.itertuples(index=False)
                },
                "exception_nights": int(len(exceptions)),
                "restart_anomaly_nights": sorted(restart_nights),
                "resource_accounting": {
                    "cpu": {"status": cpu_status, "value": cpu_value},
                    "ram": {"status": ram_status, "value": ram_value},
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
