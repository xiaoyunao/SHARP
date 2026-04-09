#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from datetime import timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

from astropy.io import fits
from astropy.time import Time


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a known-asteroid manifest by filtering L2 files to the intended observing night."
    )
    parser.add_argument("night", help="Target observing night in YYYYMMDD")
    parser.add_argument("--l2-dir", required=True, help="Directory containing L2 catalog files")
    parser.add_argument("--out", required=True, help="Output manifest path")
    parser.add_argument("--file-regex", default=r".*MP.*\.(fits|fits\.gz)$")
    parser.add_argument("--hdu", type=int, default=1, help="HDU holding the catalog/header")
    parser.add_argument("--obs-date-key", default="OBS_DATE", help="Preferred observation-time keyword")
    parser.add_argument("--timezone", default="Asia/Shanghai", help="Local timezone used to assign observing nights")
    parser.add_argument(
        "--night-rollover-hour",
        type=int,
        default=12,
        help="Local hour before which exposures are assigned to the previous observing night",
    )
    return parser.parse_args()


def read_obs_time(path: Path, hdu: int, obs_date_key: str) -> Time | None:
    try:
        with fits.open(path, memmap=False) as hdul:
            header = hdul[hdu].header
    except Exception:
        return None

    for key in (obs_date_key, "DATE-OBS", "DATEOBS"):
        if not key or key not in header:
            continue
        try:
            return Time(header[key], format="isot", scale="utc")
        except Exception:
            continue
    return None


def assign_observing_night(obs_time: Time, timezone_name: str, rollover_hour: int) -> str:
    tz = ZoneInfo(timezone_name)
    local_dt = obs_time.to_datetime(timezone=tz)
    if local_dt.hour < rollover_hour:
        local_dt = local_dt - timedelta(days=1)
    return local_dt.strftime("%Y%m%d")


def main() -> int:
    args = parse_args()
    l2_dir = Path(args.l2_dir)
    out_path = Path(args.out)
    file_re = re.compile(args.file_regex)

    if not l2_dir.is_dir():
        print(f"[FATAL] missing L2 dir: {l2_dir}", file=sys.stderr)
        return 2

    files = sorted(path for path in l2_dir.iterdir() if path.is_file() and file_re.fullmatch(path.name))

    kept: list[str] = []
    skipped_wrong_night = 0
    skipped_no_time = 0

    for path in files:
        obs_time = read_obs_time(path, args.hdu, args.obs_date_key)
        if obs_time is None:
            skipped_no_time += 1
            continue
        assigned_night = assign_observing_night(obs_time, args.timezone, args.night_rollover_hour)
        if assigned_night != args.night:
            skipped_wrong_night += 1
            continue
        kept.append(path.name)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(f"{name}\n" for name in kept), encoding="utf-8")

    print(
        f"[MANIFEST] night={args.night} total={len(files)} kept={len(kept)} "
        f"skipped_wrong_night={skipped_wrong_night} skipped_no_time={skipped_no_time}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
