#!/usr/bin/env python3
"""Build an evidence-backed operations and latency accounting table.

Inputs are the frozen server inventory (``l2_manifest.csv`` and
``nightly_products.csv``), locally copied production logs beneath
``operations_logs/{known_daily,unknown_daily}``, and an optional per-night
quality/status table.  No production path is accessed by this script.

The wrapper logs label wall-clock timestamps as CST.  In this project CST means
China Standard Time (UTC+08:00), not North American Central Standard Time.  All
such timestamps are converted explicitly to UTC before durations are computed.

File mtimes are *latest mtimes*, not creation times.  A log start is paired to a
product mtime only when the segment explicitly submitted known jobs or ran the
unknown extraction and lies within ``--max-event-product-gap-hours``.  Other
nights remain visible as ``historical_rerun``, ``latest_mtime_only``, or
``no_evidence`` and are not pooled with the normal-daily latency cohort.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import socket
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date, datetime, time, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


SCHEMA_VERSION = "1.0"
DEFAULT_START = "2025-11-15"
DEFAULT_END = "2026-07-15"
CHINA_STANDARD_TIME = timezone(timedelta(hours=8), name="CST+08:00")

LOG_TIMESTAMP_RE = re.compile(
    r"^\[(?P<stamp>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})"
    r"(?:\s+(?P<zone>CST|UTC))?\]\s*(?P<payload>.*)$"
)
KNOWN_START_RE = re.compile(r"\bknown asteroid for\s+(?P<night>\d{8})\b", re.IGNORECASE)
TARGET_NIGHT_RE = re.compile(r"\btarget_night=(?P<night>\d{8})\b")
RUN_DATE_RE = re.compile(r"\brun_date=(?P<run_date>\d{4}-\d{2}-\d{2})\b")


INVENTORY_EVENTS = (
    ("known_all_mtime", "known_all_mtime_utc"),
    ("known_match1_mtime", "known_match1_mtime_utc"),
    ("known_mask15_mtime", "known_mask15_mtime_utc"),
    ("known_ades_mtime", "known_ades_mtime_utc"),
    ("known_reply_mtime", "known_reply_mtime_utc"),
    ("unknown_catalog_mtime", "unknown_links_mtime_utc"),
    ("unknown_summary_mtime", "summary_mtime_utc"),
    ("review_manifest_mtime", "review_manifest_mtime_utc"),
    ("review_csv_mtime", "review_mtime_utc"),
    ("submit_csv_mtime", "submit_mtime_utc"),
    ("unknown_ades_mtime", "unknown_ades_mtime_utc"),
    ("unknown_reply_mtime", "unknown_reply_mtime_utc"),
)


EVENT_COLUMNS = (
    "event_id",
    "night",
    "event_type",
    "timestamp_utc",
    "timestamp_as_logged",
    "timezone_basis",
    "source_kind",
    "source_group",
    "source_path",
    "source_relative_path",
    "source_sha256",
    "line_number",
    "inventory_column",
    "line_text",
    "target_resolution",
    "segment_id",
    "run_date",
    "target_lag_days",
    "run_class",
    "trigger_class",
    "outcome",
    "outcome_tags",
    "selected_for_night",
    "selection_role",
)


NIGHT_OUTPUT_COLUMNS = (
    "night",
    "quality_code",
    "quality_reason",
    "primary_science_included",
    "unknown_science_included",
    "known_mpc_state",
    "unknown_mpc_state",
    "operations_run_class",
    "known_run_class",
    "unknown_run_class",
    "known_start_utc",
    "known_start_event_id",
    "known_start_outcome",
    "known_trigger_class",
    "known_start_candidate_n",
    "known_pairing_note",
    "unknown_start_utc",
    "unknown_start_event_id",
    "unknown_start_outcome",
    "unknown_trigger_class",
    "unknown_start_candidate_n",
    "unknown_pairing_note",
    "unknown_end_utc",
    "unknown_end_event_id",
    "l2_latest_mtime_utc",
    "l2_latest_evidence_id",
    "known_all_mtime_utc",
    "known_match1_mtime_utc",
    "known_mask15_mtime_utc",
    "known_ades_mtime_utc",
    "known_reply_mtime_utc",
    "unknown_catalog_mtime_utc",
    "unknown_summary_mtime_utc",
    "review_manifest_mtime_utc",
    "review_csv_mtime_utc",
    "submit_csv_mtime_utc",
    "unknown_ades_mtime_utc",
    "unknown_reply_mtime_utc",
    "auto_l2_latest_to_known_start_hours",
    "auto_known_start_to_ades_hours",
    "auto_l2_latest_to_known_ades_hours",
    "auto_known_reply_to_unknown_start_hours",
    "auto_known_mask15_to_unknown_start_hours",
    "auto_unknown_start_to_end_minutes",
    "auto_unknown_start_to_catalog_hours",
    "auto_unknown_start_to_review_package_hours",
    "human_review_package_to_review_csv_hours",
    "human_review_package_to_submit_csv_hours",
    "human_review_csv_to_submit_csv_hours",
    "mpc_known_ades_to_reply_hours",
    "mpc_unknown_ades_to_reply_hours",
    "end_to_end_l2_latest_to_known_reply_hours",
    "end_to_end_l2_latest_to_unknown_catalog_hours",
    "end_to_end_l2_latest_to_review_package_hours",
    "end_to_end_l2_latest_to_unknown_reply_hours",
    "timing_consistency_flags",
    "latest_mtime_caveat",
    "cpu_accounting_status",
    "ram_accounting_status",
)


LATENCY_GROUPS: dict[str, tuple[tuple[str, str, str], ...]] = {
    "automatic": (
        ("auto_l2_latest_to_known_start_hours", "hours", "known_run_class"),
        ("auto_known_start_to_ades_hours", "hours", "known_run_class"),
        ("auto_l2_latest_to_known_ades_hours", "hours", "known_run_class"),
        ("auto_known_reply_to_unknown_start_hours", "hours", "unknown_run_class"),
        ("auto_known_mask15_to_unknown_start_hours", "hours", "unknown_run_class"),
        ("auto_unknown_start_to_end_minutes", "minutes", "unknown_run_class"),
        ("auto_unknown_start_to_catalog_hours", "hours", "unknown_run_class"),
        ("auto_unknown_start_to_review_package_hours", "hours", "unknown_run_class"),
    ),
    "human": (
        ("human_review_package_to_review_csv_hours", "hours", "unknown_run_class"),
        ("human_review_package_to_submit_csv_hours", "hours", "unknown_run_class"),
        ("human_review_csv_to_submit_csv_hours", "hours", "unknown_run_class"),
    ),
    "mpc": (
        ("mpc_known_ades_to_reply_hours", "hours", "known_run_class"),
        ("mpc_unknown_ades_to_reply_hours", "hours", "unknown_run_class"),
    ),
    "end_to_end": (
        ("end_to_end_l2_latest_to_known_reply_hours", "hours", "known_run_class"),
        ("end_to_end_l2_latest_to_unknown_catalog_hours", "hours", "unknown_run_class"),
        ("end_to_end_l2_latest_to_review_package_hours", "hours", "unknown_run_class"),
        ("end_to_end_l2_latest_to_unknown_reply_hours", "hours", "unknown_run_class"),
    ),
}


@dataclass(frozen=True)
class CsvRecord:
    line_number: int
    values: dict[str, str]


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def iso_utc(value: datetime | None) -> str:
    if value is None:
        return ""
    return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def parse_inventory_timestamp(value: Any) -> datetime | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def parse_log_timestamp(match: re.Match[str]) -> tuple[datetime, str, str]:
    naive = datetime.strptime(match.group("stamp"), "%Y-%m-%d %H:%M:%S")
    zone = match.group("zone") or ""
    if zone == "UTC":
        aware = naive.replace(tzinfo=timezone.utc)
        basis = "explicit_UTC"
    elif zone == "CST":
        aware = naive.replace(tzinfo=CHINA_STANDARD_TIME)
        basis = "explicit_CST_interpreted_as_UTC+08:00"
    else:
        aware = naive.replace(tzinfo=CHINA_STANDARD_TIME)
        basis = "unlabelled_log_timestamp_assumed_CST_UTC+08:00"
    return aware.astimezone(timezone.utc), naive.isoformat(sep=" "), basis


def read_csv_records(path: Path, required: bool = True) -> list[CsvRecord]:
    if not path.is_file():
        if required:
            raise FileNotFoundError(path)
        return []
    records: list[CsvRecord] = []
    with path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        for line_number, row in enumerate(reader, start=2):
            records.append(
                CsvRecord(
                    line_number=line_number,
                    values={str(key): str(value or "") for key, value in row.items()},
                )
            )
    return records


def event_id(event: dict[str, Any]) -> str:
    key = "|".join(
        str(event.get(name, ""))
        for name in (
            "source_relative_path",
            "line_number",
            "night",
            "event_type",
            "timestamp_utc",
            "inventory_column",
        )
    )
    return "evt_" + hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]


def classify_run(night: str, timestamp_utc: datetime) -> tuple[int, str]:
    target = datetime.strptime(night, "%Y%m%d").date()
    local_date = timestamp_utc.astimezone(CHINA_STANDARD_TIME).date()
    lag = (local_date - target).days
    if lag == 1:
        return lag, "normal_daily"
    if lag > 1:
        return lag, "historical_rerun"
    return lag, "nonstandard_date_relation"


def trigger_class(path: Path, timestamp_utc: datetime) -> str:
    local = timestamp_utc.astimezone(CHINA_STANDARD_TIME)
    dated_name = path.name.startswith(local.strftime("%Y%m%d"))
    in_cron_window = time(8, 50) <= local.time() <= time(9, 15)
    if dated_name and in_cron_window and "daily_pipeline" in path.name:
        return "scheduled_0900_candidate"
    return "trigger_not_established"


def source_priority(event: dict[str, Any]) -> int:
    name = Path(str(event.get("source_path", ""))).name
    if name.endswith("_daily_pipeline.log"):
        return 3
    if name.endswith("_unknown_daily.log"):
        return 2
    if name.endswith("_daily.log"):
        return 1
    return 0


def new_log_event(
    *,
    night: str,
    event_type: str,
    timestamp_utc: datetime,
    timestamp_as_logged: str,
    timezone_basis: str,
    path: Path,
    relative_path: str,
    source_group: str,
    source_hash: str,
    line_number: int,
    line_text: str,
    target_resolution: str,
    segment_id: str,
    run_date: str = "",
) -> dict[str, Any]:
    lag, run_class = classify_run(night, timestamp_utc)
    event: dict[str, Any] = {
        "night": night,
        "event_type": event_type,
        "timestamp_utc": iso_utc(timestamp_utc),
        "timestamp_as_logged": timestamp_as_logged,
        "timezone_basis": timezone_basis,
        "source_kind": "log_line",
        "source_group": source_group,
        "source_path": str(path.resolve()),
        "source_relative_path": relative_path,
        "source_sha256": source_hash,
        "line_number": line_number,
        "inventory_column": "",
        "line_text": line_text.rstrip("\n"),
        "target_resolution": target_resolution,
        "segment_id": segment_id,
        "run_date": run_date,
        "target_lag_days": lag,
        "run_class": run_class,
        "trigger_class": trigger_class(path, timestamp_utc),
        "outcome": "observed",
        "outcome_tags": "",
        "selected_for_night": False,
        "selection_role": "",
        "_timestamp_dt": timestamp_utc,
        "_tags": set(),
    }
    event["event_id"] = event_id(event)
    return event


def apply_outcome_tag(event: dict[str, Any] | None, line: str) -> None:
    if event is None:
        return
    lower = line.lower()
    tags: set[str] = event["_tags"]
    if "[run] submit night" in lower or "[submit]" in lower:
        tags.add("submitted_jobs")
    if "run unknown extraction" in lower or "[run] mask_gaia" in lower:
        tags.add("extraction_run")
    if "[skip]" in lower:
        tags.add("skipped")
    if "missing night dir" in lower:
        tags.add("missing_night")
    if "known report not ready" in lower or "known mask15 not ready" in lower:
        tags.add("known_not_ready")
    if "already exists" in lower or "already complete" in lower:
        tags.add("existing_product")
    if "unknown review package already exists" in lower:
        tags.add("existing_review_package")
    if "unknown_count=0" in lower or "zero unknown" in lower:
        tags.add("zero_unknown")
    if "exit_code=0" in lower:
        tags.add("exit_zero")


def finalize_event_outcome(event: dict[str, Any]) -> None:
    tags: set[str] = event["_tags"]
    event["outcome_tags"] = ";".join(sorted(tags))
    if event["event_type"] == "known_start":
        if "submitted_jobs" in tags:
            event["outcome"] = "submitted_jobs"
        elif "missing_night" in tags:
            event["outcome"] = "skipped_missing_night"
        elif "skipped" in tags:
            event["outcome"] = "skipped_or_already_complete"
    elif event["event_type"] == "unknown_start":
        if "extraction_run" in tags and "completed" in tags:
            event["outcome"] = "completed_extraction"
        elif "extraction_run" in tags:
            event["outcome"] = "extraction_started_no_explicit_end"
        elif "existing_review_package" in tags:
            event["outcome"] = "existing_review_package"
        elif "known_not_ready" in tags:
            event["outcome"] = "known_not_ready"
        elif "skipped" in tags:
            event["outcome"] = "skipped"
    elif event["event_type"] == "unknown_end":
        event["outcome"] = "completed"


def parse_log_file(
    path: Path,
    *,
    logs_root: Path,
    source_hash: str,
) -> list[dict[str, Any]]:
    relative_path = str(path.relative_to(logs_root))
    source_group = relative_path.split("/", 1)[0]
    events: list[dict[str, Any]] = []
    pending_unknown: dict[str, Any] | None = None
    active_unknown: dict[str, Any] | None = None
    active_known: dict[str, Any] | None = None

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            stamp_match = LOG_TIMESTAMP_RE.match(line.rstrip("\n"))
            payload = stamp_match.group("payload") if stamp_match else line.strip()

            known_match = KNOWN_START_RE.search(payload)
            if stamp_match and known_match:
                timestamp_utc, logged, basis = parse_log_timestamp(stamp_match)
                night = known_match.group("night")
                segment_id = f"{relative_path}:{line_number}:known:{night}"
                active_known = new_log_event(
                    night=night,
                    event_type="known_start",
                    timestamp_utc=timestamp_utc,
                    timestamp_as_logged=logged,
                    timezone_basis=basis,
                    path=path,
                    relative_path=relative_path,
                    source_group=source_group,
                    source_hash=source_hash,
                    line_number=line_number,
                    line_text=line,
                    target_resolution="known_marker_contains_target_night",
                    segment_id=segment_id,
                )
                events.append(active_known)
                active_unknown = None
                pending_unknown = None
                continue

            if stamp_match and "unknown daily start" in payload.lower():
                timestamp_utc, logged, basis = parse_log_timestamp(stamp_match)
                pending_unknown = {
                    "timestamp_utc": timestamp_utc,
                    "timestamp_as_logged": logged,
                    "timezone_basis": basis,
                    "line_number": line_number,
                    "line_text": line,
                }
                active_unknown = None
                active_known = None
                continue

            target_match = TARGET_NIGHT_RE.search(line)
            if target_match and pending_unknown is not None:
                night = target_match.group("night")
                run_date_match = RUN_DATE_RE.search(line)
                segment_id = (
                    f"{relative_path}:{pending_unknown['line_number']}:unknown:{night}"
                )
                active_unknown = new_log_event(
                    night=night,
                    event_type="unknown_start",
                    timestamp_utc=pending_unknown["timestamp_utc"],
                    timestamp_as_logged=pending_unknown["timestamp_as_logged"],
                    timezone_basis=pending_unknown["timezone_basis"],
                    path=path,
                    relative_path=relative_path,
                    source_group=source_group,
                    source_hash=source_hash,
                    line_number=pending_unknown["line_number"],
                    line_text=pending_unknown["line_text"],
                    target_resolution=f"following_target_night_line_{line_number}",
                    segment_id=segment_id,
                    run_date=run_date_match.group("run_date") if run_date_match else "",
                )
                events.append(active_unknown)
                pending_unknown = None
                continue

            if target_match and active_known is not None:
                run_date_match = RUN_DATE_RE.search(line)
                if target_match.group("night") == active_known["night"]:
                    if run_date_match:
                        active_known["run_date"] = run_date_match.group("run_date")
                else:
                    active_known["_tags"].add("target_night_mismatch")
                continue

            if stamp_match and "unknown daily done" in payload.lower():
                if active_unknown is not None:
                    timestamp_utc, logged, basis = parse_log_timestamp(stamp_match)
                    active_unknown["_tags"].add("completed")
                    end_event = new_log_event(
                        night=active_unknown["night"],
                        event_type="unknown_end",
                        timestamp_utc=timestamp_utc,
                        timestamp_as_logged=logged,
                        timezone_basis=basis,
                        path=path,
                        relative_path=relative_path,
                        source_group=source_group,
                        source_hash=source_hash,
                        line_number=line_number,
                        line_text=line,
                        target_resolution="active_unknown_segment",
                        segment_id=active_unknown["segment_id"],
                        run_date=active_unknown["run_date"],
                    )
                    end_event["_tags"].add("completed")
                    events.append(end_event)
                    active_unknown = None
                continue

            if stamp_match and "daily pipeline done" in payload.lower():
                active_known = None
                active_unknown = None
                pending_unknown = None
                continue

            apply_outcome_tag(active_known, line)
            apply_outcome_tag(active_unknown, line)

    for event in events:
        finalize_event_outcome(event)
    return events


def collect_log_hashes(logs_root: Path) -> tuple[list[dict[str, Any]], dict[Path, str]]:
    rows: list[dict[str, Any]] = []
    hashes: dict[Path, str] = {}
    for group in ("known_daily", "unknown_daily"):
        group_dir = logs_root / group
        if not group_dir.is_dir():
            continue
        for path in sorted(item for item in group_dir.rglob("*") if item.is_file()):
            digest = sha256_file(path)
            hashes[path] = digest
            stat_result = path.stat()
            rows.append(
                {
                    "source_group": group,
                    "relative_path": str(path.relative_to(logs_root)),
                    "size_bytes": stat_result.st_size,
                    "mtime_utc": iso_utc(
                        datetime.fromtimestamp(stat_result.st_mtime, tz=timezone.utc)
                    ),
                    "sha256": digest,
                    "parsed_as_log": path.suffix.lower() == ".log",
                }
            )
    return rows, hashes


def inventory_event(
    *,
    night: str,
    event_type: str,
    timestamp: datetime,
    manifest_path: Path,
    manifest_hash: str,
    line_number: int,
    column: str,
    raw_value: str,
    extra_text: str = "",
) -> dict[str, Any]:
    event: dict[str, Any] = {
        "night": night,
        "event_type": event_type,
        "timestamp_utc": iso_utc(timestamp),
        "timestamp_as_logged": raw_value,
        "timezone_basis": "inventory_field_explicit_or_assumed_UTC",
        "source_kind": "inventory_mtime",
        "source_group": "inventory",
        "source_path": str(manifest_path.resolve()),
        "source_relative_path": manifest_path.name,
        "source_sha256": manifest_hash,
        "line_number": line_number,
        "inventory_column": column,
        "line_text": f"{column}={raw_value}" + (f"; {extra_text}" if extra_text else ""),
        "target_resolution": "inventory_night_column",
        "segment_id": "",
        "run_date": "",
        "target_lag_days": "",
        "run_class": "latest_mtime_evidence",
        "trigger_class": "not_applicable",
        "outcome": "latest_mtime",
        "outcome_tags": "",
        "selected_for_night": True,
        "selection_role": event_type,
        "_timestamp_dt": timestamp,
        "_tags": set(),
    }
    event["event_id"] = event_id(event)
    return event


def build_inventory_events(
    l2_records: Sequence[CsvRecord],
    product_records: Sequence[CsvRecord],
    *,
    l2_path: Path,
    products_path: Path,
    input_hashes: dict[str, str],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, dict[str, Any]]]]:
    events: list[dict[str, Any]] = []
    by_night: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)

    latest_l2: dict[str, tuple[datetime, CsvRecord]] = {}
    for record in l2_records:
        night = normalize_night(record.values.get("night"))
        timestamp = parse_inventory_timestamp(record.values.get("mtime_utc"))
        if not night or timestamp is None:
            continue
        previous = latest_l2.get(night)
        if previous is None or timestamp > previous[0]:
            latest_l2[night] = (timestamp, record)

    for night, (timestamp, record) in latest_l2.items():
        raw = record.values.get("mtime_utc", "")
        event = inventory_event(
            night=night,
            event_type="l2_latest_mtime",
            timestamp=timestamp,
            manifest_path=l2_path,
            manifest_hash=input_hashes[l2_path.name],
            line_number=record.line_number,
            column="mtime_utc",
            raw_value=raw,
            extra_text=f"path={record.values.get('path', '')}",
        )
        events.append(event)
        by_night[night][event["event_type"]] = event

    for record in product_records:
        night = normalize_night(record.values.get("night"))
        if not night:
            continue
        for event_type, column in INVENTORY_EVENTS:
            raw = record.values.get(column, "")
            timestamp = parse_inventory_timestamp(raw)
            if timestamp is None:
                continue
            event = inventory_event(
                night=night,
                event_type=event_type,
                timestamp=timestamp,
                manifest_path=products_path,
                manifest_hash=input_hashes[products_path.name],
                line_number=record.line_number,
                column=column,
                raw_value=raw,
            )
            events.append(event)
            previous = by_night[night].get(event_type)
            if previous is None or event["_timestamp_dt"] > previous["_timestamp_dt"]:
                by_night[night][event_type] = event
    return events, by_night


def select_start_event(
    events: Sequence[dict[str, Any]],
    *,
    event_type: str,
    endpoint: datetime | None,
    max_gap_hours: float,
) -> tuple[dict[str, Any] | None, str, int]:
    if event_type == "known_start":
        eligible = [event for event in events if "submitted_jobs" in event["_tags"]]
    else:
        eligible = [event for event in events if "extraction_run" in event["_tags"]]
    candidate_n = len(eligible)
    if not eligible:
        return None, "no_explicit_processing_segment", candidate_n

    if endpoint is None:
        selected = max(
            eligible,
            key=lambda event: (event["_timestamp_dt"], source_priority(event)),
        )
        return selected, "explicit_segment_without_product_endpoint", candidate_n

    paired = []
    for event in eligible:
        delta_hours = (endpoint - event["_timestamp_dt"]).total_seconds() / 3600.0
        if 0.0 <= delta_hours <= max_gap_hours:
            paired.append(event)
    if not paired:
        return None, "explicit_segments_not_pairable_to_latest_product_mtime", candidate_n
    selected = max(
        paired,
        key=lambda event: (event["_timestamp_dt"], source_priority(event)),
    )
    return selected, "explicit_segment_paired_to_latest_product_mtime", candidate_n


def duration(start: datetime | None, end: datetime | None, divisor: float = 3600.0) -> float | None:
    if start is None or end is None:
        return None
    return (end - start).total_seconds() / divisor


def combine_run_class(known_class: str, unknown_class: str) -> str:
    active = sorted({value for value in (known_class, unknown_class) if value != "no_evidence"})
    if not active:
        return "no_evidence"
    if len(active) == 1:
        return active[0]
    return "mixed_" + "__".join(active)


def truthy(value: Any) -> bool | None:
    text = str(value or "").strip().lower()
    if text in {"1", "true", "yes", "y", "t"}:
        return True
    if text in {"0", "false", "no", "n", "f"}:
        return False
    return None


def mark_selected(event: dict[str, Any] | None, role: str) -> None:
    if event is None:
        return
    event["selected_for_night"] = True
    event["selection_role"] = role


def build_night_rows(
    nights: Sequence[str],
    *,
    log_events: Sequence[dict[str, Any]],
    inventory_by_night: dict[str, dict[str, dict[str, Any]]],
    night_status: dict[str, dict[str, str]],
    max_gap_hours: float,
) -> list[dict[str, Any]]:
    log_by_night: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for event in log_events:
        log_by_night[event["night"]].append(event)
    end_by_segment = {
        event["segment_id"]: event
        for event in log_events
        if event["event_type"] == "unknown_end"
    }

    rows: list[dict[str, Any]] = []
    for night in nights:
        inventory = inventory_by_night.get(night, {})

        def endpoint(name: str) -> datetime | None:
            event = inventory.get(name)
            return event["_timestamp_dt"] if event else None

        known_completion = endpoint("known_ades_mtime") or endpoint("known_all_mtime")
        unknown_completion = (
            endpoint("unknown_catalog_mtime")
            or endpoint("review_manifest_mtime")
            or endpoint("unknown_summary_mtime")
        )
        night_logs = log_by_night.get(night, [])
        known_candidates = [event for event in night_logs if event["event_type"] == "known_start"]
        unknown_candidates = [event for event in night_logs if event["event_type"] == "unknown_start"]
        known_start, known_note, known_candidate_n = select_start_event(
            known_candidates,
            event_type="known_start",
            endpoint=known_completion,
            max_gap_hours=max_gap_hours,
        )
        unknown_start, unknown_note, unknown_candidate_n = select_start_event(
            unknown_candidates,
            event_type="unknown_start",
            endpoint=unknown_completion,
            max_gap_hours=max_gap_hours,
        )
        unknown_end = end_by_segment.get(unknown_start["segment_id"]) if unknown_start else None
        mark_selected(known_start, "known_processing_start")
        mark_selected(unknown_start, "unknown_processing_start")
        mark_selected(unknown_end, "unknown_processing_end")

        known_has_mtime = any(
            name in inventory
            for name in (
                "known_all_mtime",
                "known_match1_mtime",
                "known_mask15_mtime",
                "known_ades_mtime",
                "known_reply_mtime",
            )
        )
        unknown_has_mtime = any(
            name in inventory
            for name in (
                "unknown_catalog_mtime",
                "unknown_summary_mtime",
                "review_manifest_mtime",
                "review_csv_mtime",
                "submit_csv_mtime",
                "unknown_ades_mtime",
                "unknown_reply_mtime",
            )
        )
        known_class = (
            known_start["run_class"]
            if known_start
            else ("latest_mtime_only" if known_has_mtime else "no_evidence")
        )
        unknown_class = (
            unknown_start["run_class"]
            if unknown_start
            else ("latest_mtime_only" if unknown_has_mtime else "no_evidence")
        )

        l2 = endpoint("l2_latest_mtime")
        known_start_dt = known_start["_timestamp_dt"] if known_start else None
        unknown_start_dt = unknown_start["_timestamp_dt"] if unknown_start else None
        unknown_end_dt = unknown_end["_timestamp_dt"] if unknown_end else None

        metric_pairs = {
            "auto_l2_latest_to_known_start_hours": (l2, known_start_dt, 3600.0),
            "auto_known_start_to_ades_hours": (
                known_start_dt,
                endpoint("known_ades_mtime"),
                3600.0,
            ),
            "auto_l2_latest_to_known_ades_hours": (
                l2,
                endpoint("known_ades_mtime"),
                3600.0,
            ),
            "auto_known_reply_to_unknown_start_hours": (
                endpoint("known_reply_mtime"),
                unknown_start_dt,
                3600.0,
            ),
            "auto_known_mask15_to_unknown_start_hours": (
                endpoint("known_mask15_mtime"),
                unknown_start_dt,
                3600.0,
            ),
            "auto_unknown_start_to_end_minutes": (
                unknown_start_dt,
                unknown_end_dt,
                60.0,
            ),
            "auto_unknown_start_to_catalog_hours": (
                unknown_start_dt,
                endpoint("unknown_catalog_mtime"),
                3600.0,
            ),
            "auto_unknown_start_to_review_package_hours": (
                unknown_start_dt,
                endpoint("review_manifest_mtime"),
                3600.0,
            ),
            "human_review_package_to_review_csv_hours": (
                endpoint("review_manifest_mtime"),
                endpoint("review_csv_mtime"),
                3600.0,
            ),
            "human_review_package_to_submit_csv_hours": (
                endpoint("review_manifest_mtime"),
                endpoint("submit_csv_mtime"),
                3600.0,
            ),
            "human_review_csv_to_submit_csv_hours": (
                endpoint("review_csv_mtime"),
                endpoint("submit_csv_mtime"),
                3600.0,
            ),
            "mpc_known_ades_to_reply_hours": (
                endpoint("known_ades_mtime"),
                endpoint("known_reply_mtime"),
                3600.0,
            ),
            "mpc_unknown_ades_to_reply_hours": (
                endpoint("unknown_ades_mtime"),
                endpoint("unknown_reply_mtime"),
                3600.0,
            ),
            "end_to_end_l2_latest_to_known_reply_hours": (
                l2,
                endpoint("known_reply_mtime"),
                3600.0,
            ),
            "end_to_end_l2_latest_to_unknown_catalog_hours": (
                l2,
                endpoint("unknown_catalog_mtime"),
                3600.0,
            ),
            "end_to_end_l2_latest_to_review_package_hours": (
                l2,
                endpoint("review_manifest_mtime"),
                3600.0,
            ),
            "end_to_end_l2_latest_to_unknown_reply_hours": (
                l2,
                endpoint("unknown_reply_mtime"),
                3600.0,
            ),
        }
        metrics = {
            name: duration(start, finish, divisor)
            for name, (start, finish, divisor) in metric_pairs.items()
        }
        flags = sorted(name + "_negative" for name, value in metrics.items() if value is not None and value < 0)

        status = night_status.get(night, {})
        row: dict[str, Any] = {
            "night": night,
            "quality_code": status.get("quality_code", ""),
            "quality_reason": status.get("quality_reason", ""),
            "primary_science_included": status.get("primary_science_included", ""),
            "unknown_science_included": status.get("unknown_science_included", ""),
            "known_mpc_state": status.get("known_mpc_state", ""),
            "unknown_mpc_state": status.get("unknown_mpc_state", ""),
            "operations_run_class": combine_run_class(known_class, unknown_class),
            "known_run_class": known_class,
            "unknown_run_class": unknown_class,
            "known_start_utc": known_start["timestamp_utc"] if known_start else "",
            "known_start_event_id": known_start["event_id"] if known_start else "",
            "known_start_outcome": known_start["outcome"] if known_start else "",
            "known_trigger_class": known_start["trigger_class"] if known_start else "",
            "known_start_candidate_n": known_candidate_n,
            "known_pairing_note": known_note,
            "unknown_start_utc": unknown_start["timestamp_utc"] if unknown_start else "",
            "unknown_start_event_id": unknown_start["event_id"] if unknown_start else "",
            "unknown_start_outcome": unknown_start["outcome"] if unknown_start else "",
            "unknown_trigger_class": unknown_start["trigger_class"] if unknown_start else "",
            "unknown_start_candidate_n": unknown_candidate_n,
            "unknown_pairing_note": unknown_note,
            "unknown_end_utc": unknown_end["timestamp_utc"] if unknown_end else "",
            "unknown_end_event_id": unknown_end["event_id"] if unknown_end else "",
            "l2_latest_mtime_utc": iso_utc(l2),
            "l2_latest_evidence_id": inventory.get("l2_latest_mtime", {}).get("event_id", ""),
            "timing_consistency_flags": ";".join(flags),
            "latest_mtime_caveat": (
                "inventory timestamps are latest mtimes and may reflect reruns or copies"
                if inventory
                else ""
            ),
            "cpu_accounting_status": "unavailable_slurm_accounting_disabled",
            "ram_accounting_status": "unavailable_slurm_accounting_disabled",
        }
        for event_type, _column in INVENTORY_EVENTS:
            row[f"{event_type}_utc"] = iso_utc(endpoint(event_type))
        row.update(metrics)
        rows.append(row)
    return rows


def percentile(sorted_values: Sequence[float], quantile: float) -> float | None:
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    position = (len(sorted_values) - 1) * quantile
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(sorted_values[lower])
    weight = position - lower
    return float(sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight)


def distribution(values: Iterable[Any]) -> dict[str, Any]:
    finite: list[float] = []
    negative = 0
    for value in values:
        if value in (None, ""):
            continue
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(numeric):
            continue
        if numeric < 0:
            negative += 1
            continue
        finite.append(numeric)
    finite.sort()
    return {
        "n": len(finite),
        "n_negative_excluded": negative,
        "min": finite[0] if finite else None,
        "p16": percentile(finite, 0.16),
        "median": percentile(finite, 0.50),
        "p84": percentile(finite, 0.84),
        "p90": percentile(finite, 0.90),
        "p95": percentile(finite, 0.95),
        "max": finite[-1] if finite else None,
    }


def build_summary(
    rows: Sequence[dict[str, Any]],
    events: Sequence[dict[str, Any]],
    *,
    log_hash_rows: Sequence[dict[str, Any]],
    input_hashes: dict[str, str],
    start: date,
    end: date,
    max_gap_hours: float,
) -> dict[str, Any]:
    latency: dict[str, Any] = {}
    cohorts = ("normal_daily", "historical_rerun", "latest_mtime_only", "no_evidence")
    for category, metric_specs in LATENCY_GROUPS.items():
        category_result: dict[str, Any] = {}
        for metric, unit, class_field in metric_specs:
            metric_result: dict[str, Any] = {
                "unit": unit,
                "all_nonnegative": distribution(row.get(metric) for row in rows),
            }
            for cohort in cohorts:
                metric_result[cohort] = distribution(
                    row.get(metric) for row in rows if row.get(class_field) == cohort
                )
            metric_result["primary_science_included"] = distribution(
                row.get(metric)
                for row in rows
                if truthy(row.get("primary_science_included")) is True
            )
            category_result[metric] = metric_result
        latency[category] = category_result

    hash_lines = "".join(
        f"{row['sha256']}  {row['relative_path']}\n"
        for row in sorted(log_hash_rows, key=lambda item: item["relative_path"])
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": iso_utc(datetime.now(timezone.utc)),
        "host": socket.gethostname(),
        "python": sys.version.replace("\n", " "),
        "date_range_inclusive": {"start": start.isoformat(), "end": end.isoformat()},
        "input_hashes": input_hashes,
        "logs": {
            "file_count": len(log_hash_rows),
            "parsed_log_count": sum(bool(row["parsed_as_log"]) for row in log_hash_rows),
            "combined_sorted_sha256_manifest_hash": hashlib.sha256(
                hash_lines.encode("utf-8")
            ).hexdigest(),
        },
        "event_counts": {
            "total": len(events),
            "by_type": dict(sorted(Counter(event["event_type"] for event in events).items())),
            "selected_log_events": sum(
                event["source_kind"] == "log_line" and event["selected_for_night"]
                for event in events
            ),
        },
        "night_counts": {
            "total": len(rows),
            "operations_run_class": dict(
                sorted(Counter(row["operations_run_class"] for row in rows).items())
            ),
            "known_run_class": dict(sorted(Counter(row["known_run_class"] for row in rows).items())),
            "unknown_run_class": dict(
                sorted(Counter(row["unknown_run_class"] for row in rows).items())
            ),
        },
        "pairing": {
            "max_event_product_gap_hours": max_gap_hours,
            "rule": (
                "latest explicit successful segment at or before the relevant latest product "
                "mtime, within the configured gap"
            ),
        },
        "time_handling": {
            "log_timezone": "CST explicitly interpreted as China Standard Time UTC+08:00",
            "output_timezone": "UTC",
            "inventory_timezone": "ISO offset honored; timezone-naive inventory values assumed UTC",
        },
        "latency_categories": latency,
        "resource_accounting": {
            "cpu": {"status": "unavailable_slurm_accounting_disabled", "value": None},
            "ram": {"status": "unavailable_slurm_accounting_disabled", "value": None},
        },
        "caveats": [
            "Inventory timestamps are latest mtimes, not immutable creation timestamps.",
            "Historical reruns and latest-mtime-only nights are reported separately from normal daily runs.",
            "Human-review intervals are wall-clock mtime proxies, not active-person labor time.",
            "MPC intervals include transport and external queue/response time.",
            "CPU and peak-RAM usage cannot be recovered because Slurm accounting is disabled.",
            "Negative durations remain in the per-night CSV as audit flags but are excluded from summaries.",
        ],
    }


def atomic_write_text(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".inprogress")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def atomic_write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    temporary = path.with_name(path.name + ".inprogress")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if row.get(key) is None else row.get(key, "") for key in fields})
    os.replace(temporary, path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory-dir", type=Path, required=True)
    parser.add_argument(
        "--logs-dir",
        type=Path,
        required=True,
        help="Local snapshot directory containing known_daily/ and unknown_daily/.",
    )
    parser.add_argument("--night-status", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start", type=parse_iso_date, default=parse_iso_date(DEFAULT_START))
    parser.add_argument("--end", type=parse_iso_date, default=parse_iso_date(DEFAULT_END))
    parser.add_argument(
        "--max-event-product-gap-hours",
        type=float,
        default=72.0,
        help="Largest allowed explicit-start to latest-product-mtime pairing gap.",
    )
    return parser


def validate_paths(args: argparse.Namespace) -> tuple[Path, Path, Path, Path, Path]:
    if args.start > args.end:
        raise ValueError("--start is after --end")
    if args.max_event_product_gap_hours <= 0:
        raise ValueError("--max-event-product-gap-hours must be positive")
    inventory_dir = args.inventory_dir.expanduser().resolve(strict=False)
    logs_dir = args.logs_dir.expanduser().resolve(strict=False)
    output_dir = args.output_dir.expanduser().resolve(strict=False)
    l2_path = inventory_dir / "l2_manifest.csv"
    products_path = inventory_dir / "nightly_products.csv"
    for required in (l2_path, products_path):
        if not required.is_file():
            raise FileNotFoundError(required)
    if not logs_dir.is_dir():
        raise FileNotFoundError(logs_dir)

    output_names = (
        "pipeline_latency_by_night.csv",
        "operations_latency_summary.json",
        "operations_event_evidence.csv",
        "operations_log_hashes.csv",
        "operations_logs.sha256",
    )
    collisions = [
        output_dir / name
        for name in output_names
        if (output_dir / name).exists()
        or (output_dir / f"{name}.inprogress").exists()
    ]
    if collisions:
        raise FileExistsError(
            "refusing to overwrite existing outputs: "
            + ", ".join(str(path) for path in collisions)
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    return inventory_dir, logs_dir, output_dir, l2_path, products_path


def run(args: argparse.Namespace) -> None:
    _inventory_dir, logs_dir, output_dir, l2_path, products_path = validate_paths(args)
    nights = list(iter_nights(args.start, args.end))
    night_set = set(nights)
    night_status_path = (
        args.night_status.expanduser().resolve(strict=False) if args.night_status else None
    )
    if night_status_path is not None and not night_status_path.is_file():
        raise FileNotFoundError(night_status_path)

    l2_records = read_csv_records(l2_path)
    product_records = read_csv_records(products_path)
    status_records = read_csv_records(night_status_path, required=False) if night_status_path else []
    night_status = {
        normalize_night(record.values.get("night")): record.values
        for record in status_records
        if normalize_night(record.values.get("night"))
    }
    input_paths = [l2_path, products_path]
    if night_status_path:
        input_paths.append(night_status_path)
    input_hashes = {path.name: sha256_file(path) for path in input_paths}

    log_hash_rows, log_hashes = collect_log_hashes(logs_dir)
    log_events: list[dict[str, Any]] = []
    for path, digest in sorted(log_hashes.items(), key=lambda item: str(item[0])):
        if path.suffix.lower() != ".log":
            continue
        log_events.extend(parse_log_file(path, logs_root=logs_dir, source_hash=digest))
    log_events = [event for event in log_events if event["night"] in night_set]

    inventory_events, inventory_by_night = build_inventory_events(
        l2_records,
        product_records,
        l2_path=l2_path,
        products_path=products_path,
        input_hashes=input_hashes,
    )
    inventory_events = [event for event in inventory_events if event["night"] in night_set]
    all_events = log_events + inventory_events
    rows = build_night_rows(
        nights,
        log_events=log_events,
        inventory_by_night=inventory_by_night,
        night_status=night_status,
        max_gap_hours=args.max_event_product_gap_hours,
    )
    summary = build_summary(
        rows,
        all_events,
        log_hash_rows=log_hash_rows,
        input_hashes=input_hashes,
        start=args.start,
        end=args.end,
        max_gap_hours=args.max_event_product_gap_hours,
    )

    event_rows = sorted(
        all_events,
        key=lambda event: (
            event["night"],
            event["timestamp_utc"],
            event["event_type"],
            event["source_relative_path"],
            int(event["line_number"]),
        ),
    )
    atomic_write_csv(output_dir / "pipeline_latency_by_night.csv", rows, NIGHT_OUTPUT_COLUMNS)
    atomic_write_csv(output_dir / "operations_event_evidence.csv", event_rows, EVENT_COLUMNS)
    atomic_write_csv(
        output_dir / "operations_log_hashes.csv",
        log_hash_rows,
        ("source_group", "relative_path", "size_bytes", "mtime_utc", "sha256", "parsed_as_log"),
    )
    hash_manifest = "".join(
        f"{row['sha256']}  {row['relative_path']}\n"
        for row in sorted(log_hash_rows, key=lambda item: item["relative_path"])
    )
    atomic_write_text(output_dir / "operations_logs.sha256", hash_manifest)
    atomic_write_text(
        output_dir / "operations_latency_summary.json",
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
    )
    print(
        f"[done] nights={len(rows)} events={len(all_events)} logs={len(log_hash_rows)} "
        f"output={output_dir}",
        flush=True,
    )


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        run(args)
    except (FileExistsError, FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
