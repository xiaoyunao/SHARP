#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DONE_STATUSES = {"exported", "validated", "submitted", "no_observations", "dry_run"}
PENDING_STATUSES = {"invalid_submit_removed"}
NON_RETRYABLE_FAILURE_PATTERNS = (
    "blank is_real values",
    "missing is_real values",
    "unknown is_real values",
    "duplicate review keys",
)


def valid_night(value: str) -> bool:
    return len(value) == 8 and value.isdigit()


def load_state(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"nights": {}}
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        return {"nights": {}}
    data.setdefault("nights", {})
    return data


def save_state(path: Path, state: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)


def submit_csv_signature(path: Path) -> dict[str, Any]:
    st = path.stat()
    return {"path": str(path), "mtime_ns": st.st_mtime_ns, "size": st.st_size}


def file_signature(path: Path) -> dict[str, Any]:
    st = path.stat()
    return {"path": str(path), "mtime_ns": st.st_mtime_ns, "size": st.st_size}


def discover_review_packages(review_root: Path, start: str, end: str, only_new: set[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for night_dir in sorted(review_root.iterdir() if review_root.exists() else []):
        if not night_dir.is_dir() or not valid_night(night_dir.name):
            continue
        night = night_dir.name
        if start and night < start:
            continue
        if end and night > end:
            continue
        if only_new and night not in only_new:
            continue
        manifest = night_dir / f"{night}_unknown_review_manifest.json"
        if manifest.exists():
            out[night] = manifest
    return out


def load_manifest(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"manifest_error": str(exc)}
    return data if isinstance(data, dict) else {"manifest_error": "manifest is not a JSON object"}


def manifest_unknown_count(manifest: dict[str, Any]) -> int:
    for key in ("n_catalog_rows", "review_full_rows"):
        value = manifest.get(key)
        if value is not None:
            try:
                return int(value)
            except (TypeError, ValueError):
                pass
    rows = manifest.get("rows")
    return len(rows) if isinstance(rows, list) else 0


def discover_submit_csvs(
    review_root: Path,
    start: str,
    end: str,
    only_new: set[str],
    packages: dict[str, Path],
    allow_missing_package: bool,
) -> list[tuple[str, Path]]:
    out: list[tuple[str, Path]] = []
    for night_dir in sorted(review_root.iterdir() if review_root.exists() else []):
        if not night_dir.is_dir() or not valid_night(night_dir.name):
            continue
        night = night_dir.name
        if start and night < start:
            continue
        if end and night > end:
            continue
        if only_new and night not in only_new:
            continue
        if not allow_missing_package and night not in packages:
            continue
        submit_csv = night_dir / f"{night}_submit.csv"
        if submit_csv.exists():
            out.append((night, submit_csv))
    return out


def already_done(
    record: dict[str, Any] | None,
    submit_signature: dict[str, Any] | None,
    manifest_signature: dict[str, Any] | None,
    allow_retry_failed: bool,
) -> bool:
    if not record:
        return False
    if submit_signature is not None and record.get("submit_csv") != submit_signature:
        return False
    if manifest_signature is not None and record.get("review_package") != manifest_signature:
        return False
    status = str(record.get("status") or "")
    if status in DONE_STATUSES:
        return True
    if status in PENDING_STATUSES:
        return False
    if status == "failed" and not failure_is_retryable(record):
        return True
    if status == "failed" and not allow_retry_failed:
        return True
    return False


def failure_is_retryable(record: dict[str, Any]) -> bool:
    text = "\n".join(str(record.get(key) or "") for key in ("stderr_tail", "stdout_tail"))
    return not any(pattern in text for pattern in NON_RETRYABLE_FAILURE_PATTERNS)


def run_one(args: argparse.Namespace, night: str, submit_csv: Path) -> dict[str, Any]:
    script = Path(__file__).with_name("submit_reviewed_unknown.py")
    cmd = [
        sys.executable,
        str(script),
        night,
        "--processed-root",
        args.processed_root,
        "--review-root",
        args.review_root,
        "--submit-csv",
        str(submit_csv),
    ]
    if args.validate or args.submit:
        cmd.append("--validate")
    if args.submit:
        cmd.append("--submit")
    if args.no_logsnr:
        cmd.append("--no-logsnr")
    if args.dry_run:
        return {"status": "dry_run", "cmd": cmd, "returncode": 0}

    proc = subprocess.run(cmd, text=True, capture_output=True)
    status = "failed"
    if proc.returncode == 0:
        if "[SKIP] No unknown ADES obsData rows exported" in proc.stdout:
            status = "no_observations"
        elif args.submit:
            status = "submitted"
        elif args.validate:
            status = "validated"
        else:
            status = "exported"
    return {
        "status": status,
        "cmd": cmd,
        "returncode": proc.returncode,
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
    }


def run_followup_update(args: argparse.Namespace, night: str) -> dict[str, Any]:
    if not args.enable_followup or args.dry_run:
        return {"followup_status": "disabled"}
    repo_root = Path(__file__).resolve().parents[1]
    cmd = [
        sys.executable,
        "-m",
        "survey.apply_followup",
        "--workspace",
        args.followup_workspace,
        "--footprints",
        args.followup_footprints,
        "--processed-root",
        args.processed_root,
        "--review-root",
        args.review_root,
        "--publish-dir",
        args.followup_publish_dir,
        "--state",
        args.followup_state,
        "--start-night",
        args.followup_start_night,
        "--only-ingest-night",
        night,
        "--max-age-days",
        str(args.followup_max_age_days),
        "--obs-per-night",
        str(args.followup_obs_per_night),
    ]
    if args.followup_plan_date:
        cmd.extend(["--date", args.followup_plan_date])
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    proc = subprocess.run(cmd, text=True, capture_output=True, cwd=str(repo_root), env=env)
    return {
        "followup_status": "ok" if proc.returncode == 0 else "failed",
        "followup_cmd": cmd,
        "followup_returncode": proc.returncode,
        "followup_stdout_tail": proc.stdout[-4000:],
        "followup_stderr_tail": proc.stderr[-4000:],
    }


def remove_invalid_submit_if_needed(result: dict[str, Any], submit_csv: Path) -> None:
    if result.get("status") != "failed" or failure_is_retryable(result):
        return
    try:
        submit_csv.unlink()
    except FileNotFoundError:
        pass
    result["status"] = "invalid_submit_removed"
    result["removed_submit_csv"] = str(submit_csv)
    result["reason"] = "invalid submit CSV removed; waiting for regenerated complete submit CSV"


def mark_empty_package(
    args: argparse.Namespace,
    state: dict[str, Any],
    night: str,
    manifest_path: Path,
    manifest: dict[str, Any],
) -> bool:
    signature = file_signature(manifest_path)
    nights_state = state.setdefault("nights", {})
    record = nights_state.get(night)
    if already_done(record, None, signature, args.retry_failed):
        return False
    result = {
        "status": "no_observations",
        "night": night,
        "review_package": signature,
        "manifest_counts": {
            "n_catalog_rows": manifest.get("n_catalog_rows"),
            "review_full_rows": manifest.get("review_full_rows"),
            "review_ades_rows": manifest.get("review_ades_rows"),
        },
        "reason": "review package contains zero unknown links",
        "updated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    nights_state[night] = result
    save_state(Path(args.state), state)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True), flush=True)
    return True


def refresh_summary(
    args: argparse.Namespace,
    state: dict[str, Any],
    packages: dict[str, Path],
) -> dict[str, Any]:
    nights_state = state.setdefault("nights", {})
    pending_submit: list[str] = []
    failed: list[str] = []
    complete: list[str] = []
    manifest_errors: list[str] = []
    for night, manifest_path in sorted(packages.items()):
        manifest = load_manifest(manifest_path)
        if manifest.get("manifest_error"):
            manifest_errors.append(night)
            continue
        unknown_count = manifest_unknown_count(manifest)
        record = nights_state.get(night) or {}
        status = str(record.get("status") or "")
        manifest_sig = file_signature(manifest_path)
        submit_csv = manifest_path.parent / f"{night}_submit.csv"
        submit_sig = submit_csv_signature(submit_csv) if submit_csv.exists() else None
        current_record = bool(record.get("review_package") == manifest_sig)
        if status != "no_observations" and submit_sig is not None:
            current_record = current_record and bool(record.get("submit_csv") == submit_sig)
        if status == "failed" and current_record:
            failed.append(night)
            continue
        if status in DONE_STATUSES and current_record:
            complete.append(night)
            continue
        if unknown_count == 0:
            pending_submit.append(night)
            continue
        if not submit_csv.exists():
            pending_submit.append(night)
        else:
            pending_submit.append(night)
    summary = {
        "updated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "review_packages": len(packages),
        "complete": len(complete),
        "pending": len(pending_submit),
        "failed": len(failed),
        "manifest_errors": len(manifest_errors),
        "pending_nights": pending_submit[:50],
        "failed_nights": failed[:50],
        "manifest_error_nights": manifest_errors[:50],
    }
    state["summary"] = summary
    save_state(Path(args.state), state)
    return summary


def scan_once(args: argparse.Namespace, state: dict[str, Any]) -> int:
    only_new = set(args.only_night or [])
    packages = discover_review_packages(Path(args.review_root), args.start, args.end, only_new)
    for night, manifest_path in sorted(packages.items()):
        manifest = load_manifest(manifest_path)
        if manifest.get("manifest_error"):
            continue
        if manifest_unknown_count(manifest) == 0:
            mark_empty_package(args, state, night, manifest_path, manifest)

    found = discover_submit_csvs(
        Path(args.review_root),
        args.start,
        args.end,
        only_new,
        packages,
        args.allow_missing_review_package,
    )
    nights_state = state.setdefault("nights", {})
    processed = 0
    for night, submit_csv in found:
        signature = submit_csv_signature(submit_csv)
        manifest_signature = file_signature(packages[night]) if night in packages else None
        record = nights_state.get(night)
        if already_done(record, signature, manifest_signature, args.retry_failed):
            continue
        result = run_one(args, night, submit_csv)
        remove_invalid_submit_if_needed(result, submit_csv)
        if result.get("status") in DONE_STATUSES:
            result.update(run_followup_update(args, night))
        result.update(
            {
                "night": night,
                "submit_csv": signature,
                "review_package": manifest_signature,
                "updated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            }
        )
        nights_state[night] = result
        save_state(Path(args.state), state)
        processed += 1
        print(json.dumps(result, ensure_ascii=False, sort_keys=True), flush=True)
    refresh_summary(args, state, packages)
    return processed


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Watch reviewed unknown submit CSVs and export/validate/submit ADES.")
    ap.add_argument("--review-root", default="/pipeline/xiaoyunao/heliolincrr/review_packages")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--state", default="/pipeline/xiaoyunao/data/heliolincrr/review_submit_state.json")
    ap.add_argument("--start", default="")
    ap.add_argument("--end", default="")
    ap.add_argument("--only-night", action="append", default=[])
    ap.add_argument("--validate", action="store_true", help="Run MPC test validation before optional submit")
    ap.add_argument("--submit", action="store_true", help="Submit to MPC after validation/export")
    ap.add_argument("--no-logsnr", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument(
        "--allow-missing-review-package",
        action="store_true",
        help="Allow processing submit CSVs even if the corresponding review manifest is missing.",
    )
    ap.add_argument("--follow", action="store_true", help="Keep scanning until interrupted")
    ap.add_argument("--exit-when-complete", action="store_true", help="In --follow mode, exit when all discovered review packages are done")
    ap.add_argument("--interval-sec", type=int, default=300)
    ap.add_argument("--enable-followup", action="store_true", help="Update the follow-up state and current survey plan after a submit CSV is processed.")
    ap.add_argument("--followup-state", default="/pipeline/xiaoyunao/survey/runtime/followup/followup_state.json")
    ap.add_argument("--followup-start-night", default="20260624")
    ap.add_argument("--followup-workspace", default="/pipeline/xiaoyunao/survey/runtime")
    ap.add_argument("--followup-footprints", default="/pipeline/xiaoyunao/survey/footprints/survey_fov_footprints_with_visibility.fits")
    ap.add_argument("--followup-publish-dir", default="/pipeline/xiaoyunao/script")
    ap.add_argument("--followup-plan-date", default="", help="Optional YYYY-MM-DD plan date; default is current Asia/Shanghai date.")
    ap.add_argument("--followup-max-age-days", type=int, default=10)
    ap.add_argument("--followup-obs-per-night", type=int, default=5)
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    if args.start and not valid_night(args.start):
        raise SystemExit("--start must be YYYYMMDD")
    if args.end and not valid_night(args.end):
        raise SystemExit("--end must be YYYYMMDD")
    if args.start and args.end and args.start > args.end:
        raise SystemExit("--start must be <= --end")

    state_path = Path(args.state)
    state = load_state(state_path)
    while True:
        processed = scan_once(args, state)
        summary = state.get("summary") or {}
        if args.follow and args.exit_when_complete and summary.get("review_packages", 0) > 0:
            if summary.get("pending", 0) == 0 and summary.get("failed", 0) == 0 and summary.get("manifest_errors", 0) == 0:
                print(
                    json.dumps(
                        {
                            "status": "complete",
                            "summary": summary,
                            "updated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                        },
                        sort_keys=True,
                    ),
                    flush=True,
                )
                break
        if not args.follow:
            break
        if processed == 0:
            print(
                json.dumps(
                    {
                        "status": "idle",
                        "summary": summary,
                        "updated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
        time.sleep(max(1, int(args.interval_sec)))


if __name__ == "__main__":
    main()
