#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from survey.followup import associate_unknown_links, load_state, save_state


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Associate newly generated unknown links with active follow-up sources.")
    ap.add_argument("night", help="Unknown-link data night in YYYYMMDD format")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--catalog", default="", help="Default: /processed1/<night>/L4/<night>_unknown_links.json")
    ap.add_argument("--state", default="/pipeline/xiaoyunao/survey/runtime/followup/followup_state.json")
    ap.add_argument("--max-sep-arcsec", type=float, default=30.0)
    ap.add_argument("--dry-run", action="store_true")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    catalog = Path(args.catalog) if args.catalog else Path(args.processed_root) / args.night / "L4" / f"{args.night}_unknown_links.json"
    state_path = Path(args.state)
    state = load_state(state_path)
    stats = associate_unknown_links(
        state=state,
        night=args.night,
        catalog_path=catalog,
        max_sep_arcsec=args.max_sep_arcsec,
    )
    stats.update({"night": args.night, "catalog": str(catalog), "state": str(state_path), "dry_run": bool(args.dry_run)})
    print(json.dumps(stats, indent=2, sort_keys=True))
    if not args.dry_run:
        save_state(state_path, state)


if __name__ == "__main__":
    main()
