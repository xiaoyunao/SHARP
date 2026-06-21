#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Export, validate, and optionally submit web-reviewed unknown candidates for one night."
    )
    ap.add_argument("night", help="Target night in YYYYMMDD format")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--review-root", default="/pipeline/xiaoyunao/heliolincrr/review_packages")
    ap.add_argument("--submit-csv", default="", help="Default: <review-root>/<night>/<night>_submit.csv")
    ap.add_argument("--out", default="", help="Optional ADES PSV output path")
    ap.add_argument("--validate", action="store_true", help="Send the generated PSV to MPC test validation")
    ap.add_argument("--submit", action="store_true", help="Submit the generated PSV to MPC")
    ap.add_argument("--no-stats", action="store_true", help="Do not write the default stats JSON")
    ap.add_argument("--no-logsnr", action="store_true", help="Do not include ADES logSNR in the generated PSV")
    ap.add_argument("--response-out", default="", help="Optional MPC response output path")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    script = Path(__file__).with_name("export_unknown_ades.py")
    cmd = [
        sys.executable,
        str(script),
        args.night,
        "--processed-root",
        args.processed_root,
        "--review-root",
        args.review_root,
    ]
    if args.submit_csv:
        cmd.extend(["--submit-csv", args.submit_csv])
    else:
        cmd.append("--submit-csv")
    if args.out:
        cmd.extend(["--out", args.out])
    if args.no_stats:
        cmd.extend(["--stats-out", ""])
    if not args.no_logsnr:
        cmd.append("--include-logsnr")
    if args.response_out:
        cmd.extend(["--response-out", args.response_out])
    elif args.submit:
        cmd.extend([
            "--response-out",
            str(Path(args.processed_root) / args.night / "L4" / f"{args.night}_unknown_mpc_reply.txt"),
        ])
    elif args.validate:
        cmd.extend([
            "--response-out",
            str(Path(args.processed_root) / args.night / "L4" / f"{args.night}_unknown_validate_reply.txt"),
        ])
    if args.validate:
        cmd.append("--validate")
    if args.submit:
        cmd.append("--submit")

    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
