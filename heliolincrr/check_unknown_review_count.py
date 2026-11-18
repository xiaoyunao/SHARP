#!/usr/bin/env python3
"""Limit the human-review queue after validating current morphology results."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from unknown_morphology import paths_for
from unknown_review_gate import reviewable_rows

OVERSIZED_NIGHT = 20


def check_count(night: str, processed_root: Path, maximum: int) -> int:
    raw = json.loads(paths_for(night, processed_root)['raw'].read_text())
    selected, report = reviewable_rows(night, processed_root, raw)
    print(f"[info] unknown links after screening={len(selected)} "
          f"raw={len(raw)} counts={report['counts']} threshold={maximum}", flush=True)
    if maximum > 0 and len(selected) > maximum:
        print(f"[skip] {night} reviewable unknown links {len(selected)} exceeds "
              f"MAX_UNKNOWN_LINKS_AFTER_KNOWN={maximum} after screening; "
              "raw and screening products preserved", flush=True)
        return OVERSIZED_NIGHT
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('night')
    parser.add_argument('--processed-root', type=Path, default=Path('/processed1'))
    parser.add_argument('--maximum', type=int, default=200,
                        help='Maximum pass + review candidates; <=0 disables the count limit.')
    args = parser.parse_args()
    return check_count(args.night, args.processed_root, args.maximum)


if __name__ == '__main__':
    raise SystemExit(main())
