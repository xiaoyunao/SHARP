#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
import shutil
from zoneinfo import ZoneInfo

from astropy.table import Table

from .followup import (
    DEFAULT_MAX_AGE_DAYS,
    DEFAULT_OBS_PER_NIGHT,
    DEFAULT_START_NIGHT,
    expire_old_sources,
    ingest_reviewed_sources,
    load_state,
    reconcile_observed_followups,
    save_state,
    schedule_followups,
)
from .io_utils import write_plan_json, write_plan_txt
from .run_daily import prepare_footprints


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Inject follow-up observations into an already generated survey plan.")
    ap.add_argument("--date", default="", help="Plan date in YYYY-MM-DD. Defaults to today in Asia/Shanghai.")
    ap.add_argument("--night", default="", help="Plan night in YYYYMMDD. Overrides --date.")
    ap.add_argument("--workspace", default="/pipeline/xiaoyunao/survey/runtime")
    ap.add_argument("--footprints", default="/pipeline/xiaoyunao/survey/footprints/survey_fov_footprints_with_visibility.fits")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--review-root", default="/pipeline/xiaoyunao/heliolincrr/review_packages")
    ap.add_argument("--publish-dir", default="/pipeline/xiaoyunao/script")
    ap.add_argument("--state", default="/pipeline/xiaoyunao/survey/runtime/followup/followup_state.json")
    ap.add_argument("--start-night", default=DEFAULT_START_NIGHT)
    ap.add_argument("--only-ingest-night", default="", help="Only ingest reviewed true sources from this data night.")
    ap.add_argument("--max-age-days", type=int, default=DEFAULT_MAX_AGE_DAYS)
    ap.add_argument("--obs-per-night", type=int, default=DEFAULT_OBS_PER_NIGHT)
    ap.add_argument("--dry-run", action="store_true", help="Do not write plan/state outputs.")
    return ap.parse_args()


def resolve_night(args: argparse.Namespace) -> str:
    if args.night:
        return args.night
    if args.date:
        return datetime.strptime(args.date, "%Y-%m-%d").strftime("%Y%m%d")
    return datetime.now(ZoneInfo("Asia/Shanghai")).strftime("%Y%m%d")


def main() -> None:
    args = parse_args()
    plan_night = resolve_night(args)
    workspace = Path(args.workspace)
    plan_json = workspace / "plans" / f"{plan_night}_plan.json"
    plan_txt = workspace / "plans" / f"{plan_night}_plan.txt"
    if not plan_json.exists():
        raise FileNotFoundError(f"Missing plan JSON: {plan_json}")

    state_path = Path(args.state)
    state = load_state(state_path, args.start_night)
    ingest_stats = ingest_reviewed_sources(
        state=state,
        processed_root=Path(args.processed_root),
        review_root=Path(args.review_root),
        start_night=args.start_night,
        through_night=plan_night,
        only_night=args.only_ingest_night or None,
    )
    reconcile_stats = reconcile_observed_followups(
        state=state,
        processed_root=Path(args.processed_root),
        plan_night=plan_night,
        obs_per_night=args.obs_per_night,
    )
    expire_stats = expire_old_sources(state, plan_night, max_age_days=args.max_age_days)

    plan = json.loads(plan_json.read_text(encoding="utf-8"))
    if not isinstance(plan, list):
        raise ValueError(f"Plan JSON must contain a list: {plan_json}")
    footprints = prepare_footprints(Table.read(args.footprints))
    new_plan, schedule_stats = schedule_followups(
        plan=plan,
        state=state,
        footprints=footprints,
        plan_night=plan_night,
        obs_per_night=args.obs_per_night,
    )

    summary = {
        "plan_night": plan_night,
        "plan_json": str(plan_json),
        "state": str(state_path),
        "ingest": ingest_stats,
        "reconcile": reconcile_stats,
        "expire": expire_stats,
        "schedule": schedule_stats,
        "dry_run": bool(args.dry_run),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))

    if args.dry_run:
        return
    save_state(state_path, state)
    write_plan_json(new_plan, plan_json)
    write_plan_txt(new_plan, plan_txt)
    if args.publish_dir:
        publish_dir = Path(args.publish_dir)
        publish_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(plan_txt, publish_dir / plan_txt.name)
        shutil.copy2(plan_txt, publish_dir / "current_plan.txt")


if __name__ == "__main__":
    main()
