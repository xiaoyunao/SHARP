#!/usr/bin/env python3
"""Export compact, source-traceable SHARP annual plotting data (stdlib only).

Examples:
  python3 export_plot_data.py /path/to/SHARP --output /tmp/sharp-plot-bundle
  python3 export_plot_data.py /path/to/selected_nominal --output plot_bundle \
      --source-commit a5cb9454e56771ebe4211b7a544412056cce4c40
  python3 export_plot_data.py /path/to/one-night-fixture --output smoke_bundle \
      --expected-nights 1

The input is either a repository root, its targeted_observing directory, or a
directory containing YYYY-MM-DD.json plans. The default requires 365 consecutive
nights, including explicit zero-plan nights. This script does not schedule,
download ephemerides, calculate positions, certify photometry, or run simulations.

Outputs: fields.csv, exposures.csv.gz, night_index.json, manifest.json.
CSV list/dict cells use compact JSON; empty optional cells mean unavailable, not
zero. All required input fields are checked; an error names the file/JSON path.
The output directory may be reused only with --overwrite (four bundle files only).
"""
from __future__ import annotations

import argparse
import csv
from datetime import date, datetime, timedelta, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict


VERSION = "1.0"
DATE_FILE = re.compile(r"\d{4}-\d{2}-\d{2}\.json\Z")
OUTPUT_FILES = ("fields.csv", "exposures.csv.gz", "night_index.json", "manifest.json")
FIELD_COLUMNS = (
    "date", "block_id", "mode", "field_id", "field_version", "ra_deg", "dec_deg",
    "pa_deg", "core_targets_json", "secondary_targets_json", "anchor_object",
    "footprint_ra_deg_json", "footprint_dec_deg_json", "visits", "filters_json",
    "first_start_utc", "first_mid_utc", "last_mid_utc", "last_end_utc",
    "last_release_utc", "exposure_sum_s",
)
EXPOSURE_REQUIRED = (
    "block_id", "mode", "field_id", "field_version", "ra_deg", "dec_deg", "filter",
    "visit", "overhead_start_utc", "start_utc", "mid_utc", "end_utc", "release_utc",
    "exposure_s", "readout_s", "slew_s", "filter_change_s", "target_states",
)
EXPOSURE_COLUMNS = (
    "date", "row_in_block", "source_exposure_id",
    *[key for key in EXPOSURE_REQUIRED if key != "target_states"], "target_states_json",
)


class InputError(ValueError):
    pass


def fail(where, message):
    raise InputError(f"{where}: {message}")


def need(obj, key, where, kind=None):
    if not isinstance(obj, dict):
        fail(where, "expected a JSON object")
    if key not in obj:
        fail(f"{where}.{key}", "required field is missing; no fallback value was generated")
    value = obj[key]
    if kind is not None and not isinstance(value, kind):
        fail(f"{where}.{key}", f"expected {kind.__name__}, got {type(value).__name__}")
    return value


def text_value(value, where):
    if not isinstance(value, str) or not value.strip():
        fail(where, "expected a non-empty string")
    return value


def number(value, where):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        fail(where, "expected a finite number")
    return float(value)


def utc_seconds(value, where):
    text_value(value, where)
    try:
        t = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        fail(where, f"invalid ISO timestamp: {exc}")
    if t.tzinfo is None or t.utcoffset() != timedelta(0):
        fail(where, "expected an explicit UTC timestamp (Z or +00:00)")
    return t.timestamp()


def json_text(value):
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def json_file(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def id_list(value, where):
    if not isinstance(value, list):
        fail(where, "expected an explicit target-ID list (empty is allowed)")
    for i, item in enumerate(value):
        if isinstance(item, bool) or not isinstance(item, (str, int)) or not str(item).strip():
            fail(f"{where}[{i}]", "expected a non-empty target ID")
    return value


def coordinates(obj, where):
    ra = number(need(obj, "ra_deg", where), where + ".ra_deg")
    dec = number(need(obj, "dec_deg", where), where + ".dec_deg")
    if not (0 <= ra < 360 and -90 <= dec <= 90):
        fail(where, "coordinate outside RA [0,360), Dec [-90,90] degrees")
    return ra, dec


def resolve_input(path):
    path = path.expanduser().resolve()
    if not path.is_dir():
        fail(str(path), "input directory does not exist")
    for candidate in (path / "targeted_observing/products/annual/selected_nominal",
                      path / "products/annual/selected_nominal", path):
        if candidate.is_dir() and any(DATE_FILE.fullmatch(p.name) for p in candidate.iterdir()):
            return candidate
    fail(str(path), "no YYYY-MM-DD.json plans found; expected products/annual/selected_nominal")


def commit_metadata(plan_dir, declared):
    observed = None
    for parent in (plan_dir, *plan_dir.parents):
        if (parent / ".git").exists():
            try:
                observed = subprocess.run(
                    ["git", "-C", str(parent), "rev-parse", "HEAD"], check=True,
                    capture_output=True, text=True, timeout=10,
                ).stdout.strip()
            except (OSError, subprocess.SubprocessError):
                pass
            break
    if declared and observed and declared.lower() != observed.lower():
        fail("--source-commit", f"declared commit differs from local repository HEAD {observed}")
    if not declared and not observed:
        fail("--source-commit", "source commit is required outside a Git worktree; no provenance was inferred")
    return {"commit": declared or observed,
            "commit_basis": "caller declaration" if declared else "local git HEAD" if observed else "unavailable",
            "local_git_head": observed,
            "caveat": "A code commit alone does not prove which code generated plans; input file hashes and source_policy are retained."}


def validate_plan(plan, filename):
    """Validate the export contract, not astronomical or scientific correctness."""
    day = need(plan, "date", filename, str)
    if day != filename[:-5]:
        fail(filename + ".date", "date must match the input filename")
    config = need(plan, "config", filename, dict)
    site = need(config, "site", filename + ".config", dict)
    for key in ("lat_deg", "lon_deg", "height_m"):
        number(need(site, key, filename + ".config.site"), filename + ".config.site." + key)
    number(need(config, "fov_deg", filename + ".config"), filename + ".config.fov_deg")
    night = need(plan, "night", filename, dict)
    for key in ("evening", "dark", "morning"):
        values = need(night, key, filename + ".night", list)
        if len(values) != 2:
            fail(filename + ".night." + key, "expected [start,end] in UTC Unix seconds")
        lo, hi = [number(v, filename + ".night." + key) for v in values]
        if hi < lo:
            fail(filename + ".night." + key, "window end precedes start")
    text_value(need(plan, "status", filename), filename + ".status")
    need(plan, "observations_completed", filename, bool)
    need(plan, "skips", filename, list)
    free = need(plan, "free_intervals", filename, list)
    for j, item in enumerate(free):
        where = f"{filename}.free_intervals[{j}]"
        a = utc_seconds(need(item, "start_utc", where), where + ".start_utc")
        b = utc_seconds(need(item, "end_utc", where), where + ".end_utc")
        duration = number(need(item, "duration_s", where), where + ".duration_s")
        if b < a or abs((b - a) - duration) > .01:
            fail(where, "free-interval times and duration disagree")
    blocks = need(plan, "blocks", filename, list)
    ids = set()
    for j, block in enumerate(blocks):
        where = f"{filename}.blocks[{j}]"
        blockid = text_value(need(block, "block_id", where), where + ".block_id")
        if blockid in ids:
            fail(where + ".block_id", "duplicate block ID within this night")
        ids.add(blockid)
        text_value(need(block, "mode", where), where + ".mode")
        for key in ("not_before_utc", "hard_deadline_utc", "release_utc"):
            utc_seconds(need(block, key, where), where + "." + key)
        need(block, "fields", where, list)
        need(block, "rows", where, list)
    return config


def export_night(plan, filename, field_writer, exposure_writer, optional_missing):
    field_count = exposure_count = 0
    for bi, block in enumerate(plan["blocks"]):
        where = f"{filename}.blocks[{bi}]"
        fields = {}
        for fi, f in enumerate(block["fields"]):
            fw = f"{where}.fields[{fi}]"
            fid = text_value(need(f, "field_id", fw), fw + ".field_id")
            if fid in fields:
                fail(fw, "duplicate field ID within this block")
            text_value(need(f, "version", fw), fw + ".version")
            coordinates(f, fw)
            number(need(f, "pa_deg", fw), fw + ".pa_deg")
            id_list(need(f, "core_targets", fw), fw + ".core_targets")
            if "secondary_targets" not in f:
                if block["mode"] == "trojan":
                    fail(fw + ".secondary_targets", "required for final Trojan plans; explicitly use [] if none")
                optional_missing["twilight_secondary_targets_absent_means_empty"] += 1
            else:
                id_list(f["secondary_targets"], fw + ".secondary_targets")
            fields[fid] = f
        selected = defaultdict(list)
        for ri, row in enumerate(block["rows"], 1):
            rw = f"{where}.rows[{ri - 1}]"
            for key in EXPOSURE_REQUIRED:
                need(row, key, rw)
            fid = row["field_id"]
            if fid not in fields:
                fail(rw + ".field_id", "no matching field definition in this block")
            f = fields[fid]
            if row["block_id"] != block["block_id"] or row["mode"] != block["mode"] or row["field_version"] != f["version"]:
                fail(rw, "row block/mode/version disagrees with its containing block or field")
            ra, dec = coordinates(row, rw)
            if abs(ra - f["ra_deg"]) > 1e-8 or abs(dec - f["dec_deg"]) > 1e-8:
                fail(rw, "exposure coordinates differ from its fixed field; export needs an explicit moving-field schema")
            text_value(row["filter"], rw + ".filter")
            if isinstance(row["visit"], bool) or not isinstance(row["visit"], int) or row["visit"] < 1:
                fail(rw + ".visit", "expected a positive integer")
            ts = [utc_seconds(row[k], rw + "." + k) for k in
                  ("overhead_start_utc", "start_utc", "mid_utc", "end_utc", "release_utc")]
            if any(a > b for a, b in zip(ts, ts[1:])):
                fail(rw, "exposure timestamps are not ordered")
            for key in ("exposure_s", "readout_s", "slew_s", "filter_change_s"):
                if number(row[key], rw + "." + key) < 0:
                    fail(rw + "." + key, "duration must be nonnegative")
            if abs((ts[3] - ts[1]) - row["exposure_s"]) > .003 or abs((ts[4] - ts[3]) - row["readout_s"]) > .003:
                fail(rw, "exposure/readout duration disagrees with timestamps")
            if not isinstance(row["target_states"], dict):
                fail(rw + ".target_states", "expected a JSON object")
            out = {key: row[key] for key in EXPOSURE_REQUIRED if key != "target_states"}
            out.update(date=plan["date"], row_in_block=ri, source_exposure_id=row.get("exposure_id", ""),
                       target_states_json=json_text(row["target_states"]))
            exposure_writer.writerow(out)
            selected[fid].append((ts[1], row))
            exposure_count += 1
        for fid, f in fields.items():
            if not selected[fid]:
                fail(where + ".fields." + fid, "field has no exposures; no first/last time can be exported")
            rows = [r for _, r in sorted(selected[fid], key=lambda x: x[0])]
            first = rows[0]
            last = max(rows, key=lambda r: utc_seconds(r["end_utc"], where))
            out = {"date": plan["date"], "block_id": block["block_id"], "mode": block["mode"],
                   "field_id": fid, "field_version": f["version"], "ra_deg": f["ra_deg"], "dec_deg": f["dec_deg"],
                   "pa_deg": f["pa_deg"], "core_targets_json": json_text(f["core_targets"]),
                   "secondary_targets_json": json_text(f.get("secondary_targets", [])), "anchor_object": f.get("anchor_object", ""),
                   "visits": len(rows), "filters_json": json_text(sorted({r["filter"] for r in rows})),
                   "first_start_utc": first["start_utc"], "first_mid_utc": first["mid_utc"],
                   "last_mid_utc": last["mid_utc"], "last_end_utc": last["end_utc"],
                   "last_release_utc": max(rows, key=lambda r: utc_seconds(r["release_utc"], where))["release_utc"],
                   "exposure_sum_s": sum(r["exposure_s"] for r in rows)}
            for key in ("footprint_ra_deg", "footprint_dec_deg"):
                out[key + "_json"] = json_text(f[key]) if key in f else ""
                if key not in f:
                    optional_missing[key] += 1
            field_writer.writerow(out)
            field_count += 1
    return field_count, exposure_count


def run(args):
    plan_dir = resolve_input(args.input)
    paths = sorted(p for p in plan_dir.iterdir() if p.is_file() and DATE_FILE.fullmatch(p.name))
    days = [date.fromisoformat(p.stem) for p in paths]
    if len(paths) != args.expected_nights:
        fail(str(plan_dir), f"found {len(paths)} dated plans; expected {args.expected_nights}. Use --expected-nights only for an explicitly smaller/larger dataset")
    if any(b - a != timedelta(days=1) for a, b in zip(days, days[1:])):
        fail(str(plan_dir), "calendar has missing nights; zero-observation nights still need explicit plan JSON")
    if args.start_date and days[0].isoformat() != args.start_date:
        fail("--start-date", f"first input date is {days[0]}")
    stop = (days[-1] + timedelta(days=1)).isoformat()
    if args.stop_date_exclusive and stop != args.stop_date_exclusive:
        fail("--stop-date-exclusive", f"input exclusive stop is {stop}")
    provenance = commit_metadata(plan_dir, args.source_commit)
    output = args.output.expanduser().resolve()
    if output == plan_dir:
        fail("--output", "must differ from the source plan directory")
    if output.exists() and not output.is_dir():
        fail("--output", "output exists and is not a directory")
    existing = [name for name in OUTPUT_FILES if (output / name).exists()]
    if existing and not args.overwrite:
        fail("--output", f"bundle files already exist: {', '.join(existing)}; pass --overwrite to replace these files only")
    output.parent.mkdir(parents=True, exist_ok=True)
    configs = {}; nights = []; inputs = []; optional_missing = Counter()
    counts = {"nights": len(paths), "blocks": 0, "field_visits": 0, "exposures": 0, "zero_plan_nights": 0}
    with tempfile.TemporaryDirectory(prefix=".plot-export-", dir=output.parent) as tmp:
        stage = Path(tmp)
        with (stage / "fields.csv").open("w", encoding="utf-8", newline="") as fh, (stage / "exposures.csv.gz").open("wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz, io.TextIOWrapper(gz, encoding="utf-8", newline="") as eh:
                fw = csv.DictWriter(fh, fieldnames=FIELD_COLUMNS); fw.writeheader()
                ew = csv.DictWriter(eh, fieldnames=EXPOSURE_COLUMNS); ew.writeheader()
                for path in paths:
                    data = path.read_bytes()
                    try:
                        plan = json.loads(data)
                    except (ValueError, UnicodeError) as exc:
                        fail(path.name, f"cannot parse JSON: {exc}")
                    config = validate_plan(plan, path.name)
                    config_id = sha256(json_text(config).encode("utf-8"))
                    # Small configurations aid altitude/footprint reconstruction; no caches are read.
                    configs[config_id] = config
                    nf, ne = export_night(plan, path.name, fw, ew, optional_missing)
                    counts["field_visits"] += nf; counts["exposures"] += ne
                    counts["blocks"] += len(plan["blocks"]); counts["zero_plan_nights"] += not plan["blocks"]
                    digest = sha256(data)
                    inputs.append({"file": path.name, "sha256": digest, "bytes": len(data)})
                    night = {key: plan[key] for key in ("date", "night", "free_intervals", "skips", "status", "observations_completed")}
                    night.update(config_id=config_id, source_file=path.name, source_sha256=digest,
                                 blocks=[{k: b[k] for k in ("block_id", "mode", "not_before_utc", "hard_deadline_utc", "release_utc")} for b in plan["blocks"]])
                    for key in ("source_policy", "schema_version", "free_intervals_scope", "nominal_feedback", "validation",
                                "field_version_start", "field_version_stop", "iers", "reviewed_head", "reconstruction", "smoke_fixture"):
                        if key in plan:
                            night[key] = plan[key]
                    nights.append(night)
        json_file(stage / "night_index.json", {"schema_version": 1, "nights": nights})
        manifest = {
            "schema_version": 1, "exporter_version": VERSION,
            "exported_utc": datetime.now(timezone.utc).isoformat(),
            "exporter_sha256": sha256(Path(__file__).read_bytes()),
            "source": {**provenance, "directory_name": plan_dir.name, "plan_files": inputs,
                       "source_policies": sorted({p["source_policy"] for p in nights if p.get("source_policy")}),
                       "statuses": sorted({p["status"] for p in nights}),
                       "observations_completed_values": sorted({p["observations_completed"] for p in nights})},
            "calendar": {"start_date": str(days[0]), "stop_date_exclusive": stop,
                         "expected_nights": args.expected_nights, "consecutive": True},
            "counts": counts, "configurations": configs,
            "conventions": {
                "dates": "Original observing-night labels, not necessarily exposure UTC dates.",
                "timestamps": "Copied UTC ISO strings; night.evening/dark/morning retain UTC Unix seconds.",
                "coordinates": "Copied equatorial field centers in degrees; no positions calculated or propagated.",
                "csv_json_cells": "Columns ending _json contain JSON; empty optional cells mean unavailable.",
                "fields_grain": "One row per night / block / field; visits is the actual number of copied exposure rows.",
                "exposure_order": "Original block row order; row_in_block is 1-based. source_exposure_id is empty when not supplied.",
                "target_ids": "core_targets and secondary_targets retain source membership; no ephemerides or scientific classifications inferred.",
                "missing_secondary_targets": "Allowed only for non-Trojan fields and represented as []; counted below.",
                "footprints": "Copied when supplied; missing vertices are not fabricated. Use config FOV and source position angle if computing a separate schematic.",
            },
            "optional_missing_counts": dict(optional_missing),
            "scope": "Plotting-only export. Counts describe source plans, not executed observations, detections, lightcurve yield, or an independent geometry validation.",
            "not_included": ["Ephemeris/stellar caches", "200-weather realization products", "target sky positions absent from plans", "new altitude or Sun/Moon calculations"],
            "outputs": {name: {"bytes": (stage / name).stat().st_size, "sha256": sha256((stage / name).read_bytes())} for name in OUTPUT_FILES if name != "manifest.json"},
        }
        json_file(stage / "manifest.json", manifest)
        output.mkdir(parents=True, exist_ok=True)
        for name in OUTPUT_FILES:  # Manifest is last: only a complete, validated bundle is published.
            shutil.move(str(stage / name), str(output / name))
    print(json.dumps({"output": str(output), "counts": counts, "calendar": manifest["calendar"],
                      "source_commit": provenance["commit"]}, ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=Path, help="SHARP repo root, targeted_observing, or dated-plan directory")
    parser.add_argument("--output", required=True, type=Path, help="directory for the four plot-bundle files")
    parser.add_argument("--expected-nights", type=int, default=365, help="required consecutive plan count (default: 365; use 1 for a smoke fixture)")
    parser.add_argument("--start-date", help="optional required first night, YYYY-MM-DD")
    parser.add_argument("--stop-date-exclusive", help="optional required exclusive last date, YYYY-MM-DD")
    parser.add_argument("--source-commit", help="full 40-character source commit, for copied products without .git")
    parser.add_argument("--overwrite", action="store_true", help="replace the four existing bundle files after validating all inputs")
    args = parser.parse_args()
    try:
        if args.expected_nights < 1:
            fail("--expected-nights", "must be positive")
        if args.source_commit and not re.fullmatch(r"[0-9a-fA-F]{40}", args.source_commit):
            fail("--source-commit", "expected a full 40-character Git SHA")
        for key in ("start_date", "stop_date_exclusive"):
            value = getattr(args, key)
            if value:
                date.fromisoformat(value)
        run(args)
    except (InputError, OSError, ValueError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
