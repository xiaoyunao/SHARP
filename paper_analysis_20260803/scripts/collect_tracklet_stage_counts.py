#!/usr/bin/env python3
"""Collect per-night edge/static/tracklet counts from frozen production logs.

This collector is intentionally small and read-only.  It parses the explicit
per-group count lines written by ``make_tracklet_and_idx_v3.py`` and writes one
auditable row per night.  No production file is modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timedelta, timezone
from pathlib import Path


GROUP_RE = re.compile(r"^\[group\s+(\d+)\]\s+start\s+n_exp=(\d+)\s*$")
LIST_PATTERNS = {
    "grouped_gaia_input_detection_n": re.compile(r"^\[group\s+(\d+)\]\s+raw Ns=\[([^]]*)\]"),
    "common_area_survivor_detection_n": re.compile(
        r"^\[group\s+(\d+)\]\s+after common-area Ns=\[([^]]*)\]"
    ),
    "edge_shell_survivor_detection_n": re.compile(
        r"^\[group\s+(\d+)\]\s+after edge-shell filter:\s*\[([^]]*)\]"
    ),
}
STATIC_RE = re.compile(r"^\[group\s+(\d+)\]\s+after static removal:\s*(.*)$")
STATIC_VALUE_RE = re.compile(r"exp\d+=(\d+)")
DONE_RE = re.compile(r"^\[group\s+(\d+)\]\s+DONE tracklets=(\d+)\b")
FAILED_RE = re.compile(r"^\[group\s+(\d+)\].*(?:FAILED|ERROR|Traceback)", re.IGNORECASE)
PARAM_RE = re.compile(r"^params:\s*(.*)$")
DECLARED_GROUPS_RE = re.compile(r"\bn_groups=(\d+)\b")
SKIP_RE = re.compile(r"^\[group\s+(\d+)\]\s+skip:\s*(.*)$")


def nights(start: str, end: str):
    current = datetime.strptime(start, "%Y-%m-%d").date()
    finish = datetime.strptime(end, "%Y-%m-%d").date()
    while current <= finish:
        yield current.strftime("%Y%m%d")
        current += timedelta(days=1)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_int_list(text: str) -> list[int]:
    if not text.strip():
        return []
    return [int(part.strip()) for part in text.split(",") if part.strip()]


def parse_log(path: Path, night: str) -> dict[str, object]:
    record: dict[str, object] = {
        "night": night,
        "path": str(path),
        "exists": path.exists(),
        "mtime_utc": "",
        "sha256": "",
        "parse_status": "missing",
        "parameter_line": "",
        "declared_group_n": "",
        "skip_reason_counts": "{}",
        "group_started_n": 0,
        "group_complete_n": 0,
        "group_failed_n": 0,
        "exposure_group_membership_n": 0,
        "grouped_gaia_input_detection_n": 0,
        "common_area_survivor_detection_n": 0,
        "edge_shell_survivor_detection_n": 0,
        "static_survivor_detection_n": 0,
        "tracklet_n": 0,
        "missing_stage_group_n": 0,
        "duplicate_stage_line_n": 0,
    }
    if not path.exists():
        return record

    record["mtime_utc"] = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()
    record["sha256"] = sha256(path)
    per_group: dict[str, dict[str, object]] = {}
    duplicate_lines = 0
    skip_reasons: Counter[str] = Counter()

    def group(group_id: str) -> dict[str, object]:
        return per_group.setdefault(group_id, {})

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        declared = DECLARED_GROUPS_RE.search(line)
        if declared and record["declared_group_n"] == "":
            record["declared_group_n"] = int(declared.group(1))
        parameter = PARAM_RE.match(line)
        if parameter:
            record["parameter_line"] = parameter.group(1).strip()
        started = GROUP_RE.match(line)
        if started:
            target = group(started.group(1))
            if "n_exp" in target:
                duplicate_lines += 1
            target["n_exp"] = int(started.group(2))
            continue
        matched_stage = False
        for key, pattern in LIST_PATTERNS.items():
            found = pattern.match(line)
            if not found:
                continue
            target = group(found.group(1))
            if key in target:
                duplicate_lines += 1
            target[key] = sum(parse_int_list(found.group(2)))
            matched_stage = True
            break
        if matched_stage:
            continue
        static = STATIC_RE.match(line)
        if static:
            target = group(static.group(1))
            if "static_survivor_detection_n" in target:
                duplicate_lines += 1
            target["static_survivor_detection_n"] = sum(
                int(value) for value in STATIC_VALUE_RE.findall(static.group(2))
            )
            continue
        done = DONE_RE.match(line)
        if done:
            target = group(done.group(1))
            if "tracklet_n" in target:
                duplicate_lines += 1
            target["tracklet_n"] = int(done.group(2))
            target["complete"] = True
            continue
        skipped = SKIP_RE.match(line)
        if skipped:
            target = group(skipped.group(1))
            reason = skipped.group(2).strip()
            skip_reasons[reason] += 1
            target["skip_reason"] = reason
            target.setdefault("tracklet_n", 0)
            if reason.startswith("fewer than 2 exposures"):
                target.setdefault("edge_shell_survivor_detection_n", 0)
                target.setdefault("static_survivor_detection_n", 0)
            target["complete"] = True
            continue
        failed = FAILED_RE.match(line)
        if failed:
            group(failed.group(1))["failed"] = True

    required = tuple(LIST_PATTERNS) + ("static_survivor_detection_n", "tracklet_n")
    missing = Counter()
    for values in per_group.values():
        record["exposure_group_membership_n"] += int(values.get("n_exp", 0))
        for key in required:
            if key in values:
                record[key] += int(values[key])
            else:
                missing[key] += 1
        record["group_complete_n"] += int(bool(values.get("complete", False)))
        record["group_failed_n"] += int(bool(values.get("failed", False)))
    record["group_started_n"] = len(per_group)
    record["skip_reason_counts"] = json.dumps(dict(sorted(skip_reasons.items())), sort_keys=True)
    record["missing_stage_group_n"] = max(missing.values(), default=0)
    record["duplicate_stage_line_n"] = duplicate_lines
    if not per_group and record["declared_group_n"] == 0:
        record["parse_status"] = "zero_groups"
    elif not per_group:
        record["parse_status"] = "no_group_lines"
    elif missing or duplicate_lines or record["group_failed_n"]:
        record["parse_status"] = "incomplete_or_ambiguous"
    else:
        record["parse_status"] = "ok"
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default="/pipeline/xiaoyunao/data/heliolincrr")
    parser.add_argument("--start", default="2025-11-15")
    parser.add_argument("--end", default="2026-07-15")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output}")
    output.mkdir(parents=True)
    data_root = Path(args.data_root)
    rows = [
        parse_log(data_root / night / "tracklets_linreproj" / "make_tracklet.log", night)
        for night in nights(args.start, args.end)
    ]
    fieldnames = list(rows[0]) if rows else []
    csv_path = output / "tracklet_stage_counts_by_night.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    status_counts = Counter(str(row["parse_status"]) for row in rows)
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_root": str(data_root),
        "start": args.start,
        "end": args.end,
        "night_rows": len(rows),
        "status_counts": dict(sorted(status_counts.items())),
        "ok_nights": status_counts.get("ok", 0),
        "true_missing_log_nights": status_counts.get("missing", 0),
        "totals_for_complete_nights": {
            key: sum(
                int(row[key])
                for row in rows
                if row["parse_status"] in {"ok", "zero_groups"}
            )
            for key in tuple(LIST_PATTERNS) + ("static_survivor_detection_n", "tracklet_n")
        },
        "output_sha256": sha256(csv_path),
    }
    summary_path = output / "tracklet_stage_counts_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "hashes.sha256").write_text(
        f"{sha256(csv_path)}  {csv_path.name}\n{sha256(summary_path)}  {summary_path.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
