#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path


def valid_night(value: str) -> bool:
    return len(value) == 8 and value.isdigit()


def discover_nights(processed_root: Path, start: str, end: str) -> list[str]:
    nights: list[str] = []
    for path in sorted(processed_root.iterdir()):
        if path.is_dir() and valid_night(path.name) and start <= path.name <= end:
            nights.append(path.name)
    return nights


def run_cmd(cmd: list[str], dry_run: bool) -> None:
    print("[run] " + " ".join(cmd), flush=True)
    if not dry_run:
        subprocess.run(cmd, check=True)


def count_json_rows(path: Path) -> int:
    if not path.exists():
        return 0
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"Expected JSON list: {path}")
    return len(data)


def remove_unknown_outputs(processed_root: Path, night: str) -> None:
    l4 = processed_root / night / "L4"
    for path in (
        l4 / f"{night}_unknown_links.json",
        l4 / f"{night}_unknown_links.fits",
        l4 / f"{night}_unknown_links_ades.psv",
        l4 / f"{night}_unknown_mpc_reply.txt",
        l4 / f"{night}_unknown_mpc_validate_reply.txt",
        l4 / f"{night}_unknown_validate_reply.txt",
    ):
        if path.exists():
            path.unlink()
            print(f"[remove] {path}", flush=True)


def clean_review_artifacts(args: argparse.Namespace, night: str) -> None:
    plot_dir = Path(args.plot_root) / night
    if plot_dir.exists():
        for path in plot_dir.glob(f"unknown_link_*_{night}.gif"):
            path.unlink()
            print(f"[remove] {path}", flush=True)

    review_dir = Path(args.review_root) / night
    if review_dir.exists():
        shutil.rmtree(review_dir)
        print(f"[remove] {review_dir}", flush=True)


def remask_one(args: argparse.Namespace, night: str) -> dict[str, object]:
    processed_root = Path(args.processed_root)
    root_out = Path(args.root_out)
    script_dir = Path(__file__).resolve().parent
    l4 = processed_root / night / "L4"
    analysis_dir = root_out / night / "analysis"
    rr_dir = root_out / night / args.rr_subdir
    tracklets = root_out / night / "tracklets_linreproj" / f"tracklets_{night}_ALL.fits"
    matched = l4 / f"{night}_matched_asteroids{args.mask_matched_suffix}.fits"
    summary_json = analysis_dir / f"{night}_single_night_summary.json"
    unknown_json = l4 / f"{night}_unknown_links.json"
    unknown_fits = l4 / f"{night}_unknown_links.fits"

    required = [
        matched,
        tracklets,
        rr_dir / "links_tracklets.fits",
        rr_dir / "linkage_members.fits",
        rr_dir / "orbit_confirm" / "orbit_links.fits",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        return {"night": night, "status": "skip", "reason": "missing_inputs", "missing": ";".join(missing)}

    if args.dry_run:
        return {"night": night, "status": "dry_run", "reason": "", "matched": str(matched)}

    analysis_dir.mkdir(parents=True, exist_ok=True)
    l4.mkdir(parents=True, exist_ok=True)

    run_cmd(
        [
            sys.executable,
            str(script_dir / "summarize_single_night.py"),
            night,
            "--processed-root",
            str(processed_root),
            "--root-out",
            str(root_out),
            "--rr-subdir",
            args.rr_subdir,
            "--matched",
            str(matched),
            "--summary-json",
            str(summary_json),
            "--unknown-json",
            str(unknown_json),
            "--unknown-fits",
            str(unknown_fits),
        ],
        dry_run=False,
    )

    unknown_count = count_json_rows(unknown_json)
    if args.max_unknown_links_after_known > 0 and unknown_count > args.max_unknown_links_after_known:
        remove_unknown_outputs(processed_root, night)
        return {
            "night": night,
            "status": "skip",
            "reason": "unknown_links_after_known_gt_limit",
            "unknown_count": unknown_count,
            "limit": args.max_unknown_links_after_known,
            "matched": str(matched),
        }

    assign_stats: dict[str, object] = {}
    if not args.no_assign_trksub:
        cmd = [
            sys.executable,
            str(script_dir / "assign_unknown_trksub.py"),
            night,
            "--mode",
            args.mode,
            "--catalog",
            str(unknown_json),
            "--fits",
            str(unknown_fits),
            "--summary-json",
            str(summary_json),
            "--history",
            args.trksub_history,
        ]
        print("[run] " + " ".join(cmd), flush=True)
        proc = subprocess.run(cmd, check=True, text=True, capture_output=True)
        print(proc.stdout, end="", flush=True)
        assign_stats = json.loads(proc.stdout)

    package_stats: dict[str, object] = {}
    if not args.skip_package:
        if args.clean_review_artifacts:
            clean_review_artifacts(args, night)

        if not args.skip_plots:
            cmd = [
                sys.executable,
                str(script_dir / "plot_unknown_links.py"),
                night,
                "--processed-root",
                str(processed_root),
                "--root-out",
                str(root_out),
                "--plot-root",
                args.plot_root,
                "--catalog",
                str(unknown_json),
                "--gif-size",
                str(args.gif_size),
                "--gif-duration",
                str(args.gif_duration),
            ]
            if args.plot_limit_links > 0:
                cmd.extend(["--limit-links", str(args.plot_limit_links)])
            print("[run] " + " ".join(cmd), flush=True)
            subprocess.run(cmd, check=True)

        cmd = [
            sys.executable,
            str(script_dir / "package_unknown_review.py"),
            night,
            "--processed-root",
            str(processed_root),
            "--root-out",
            str(root_out),
            "--plot-root",
            args.plot_root,
            "--out-root",
            args.review_root,
            "--make-tar",
        ]
        print("[run] " + " ".join(cmd), flush=True)
        proc = subprocess.run(cmd, check=True, text=True, capture_output=True)
        print(proc.stdout, end="", flush=True)
        package_stats = json.loads(proc.stdout)

    return {
        "night": night,
        "status": "done",
        "reason": "",
        "unknown_count": unknown_count,
        "matched": str(matched),
        "assigned_new": assign_stats.get("assigned_new", ""),
        "reused_existing": assign_stats.get("reused_existing", ""),
        "review_full_rows": package_stats.get("review_full_rows", ""),
        "review_ades_rows": package_stats.get("review_ades_rows", ""),
        "n_gifs_missing": package_stats.get("n_gifs_missing", ""),
    }


def write_status(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "night",
        "status",
        "reason",
        "unknown_count",
        "limit",
        "assigned_new",
        "reused_existing",
        "review_full_rows",
        "review_ades_rows",
        "n_gifs_missing",
        "matched",
        "missing",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Rebuild unknown catalogs using the wider known-asteroid mask match table.")
    ap.add_argument("start", help="Start night YYYYMMDD")
    ap.add_argument("end", help="End night YYYYMMDD")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--root-out", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--rr-subdir", default="rr_links")
    ap.add_argument("--mask-matched-suffix", default="_mask15")
    ap.add_argument("--mode", default="single-night")
    ap.add_argument("--trksub-history", default="/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl")
    ap.add_argument("--plot-root", default="/pipeline/xiaoyunao/heliolincrr/plots")
    ap.add_argument("--review-root", default="/pipeline/xiaoyunao/heliolincrr/review_packages")
    ap.add_argument("--status-out", default="")
    ap.add_argument("--max-unknown-links-after-known", type=int, default=200)
    ap.add_argument("--skip-plots", action="store_true", help="Do not regenerate unknown GIFs before packaging")
    ap.add_argument(
        "--no-clean-review-artifacts",
        dest="clean_review_artifacts",
        action="store_false",
        help="Keep existing unknown GIFs and review package files before rebuilding",
    )
    ap.add_argument("--gif-size", type=int, default=280)
    ap.add_argument("--gif-duration", type=float, default=0.6)
    ap.add_argument("--plot-limit-links", type=int, default=0)
    ap.add_argument("--exclude-night", action="append", default=[])
    ap.add_argument("--no-assign-trksub", action="store_true")
    ap.add_argument("--skip-package", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.set_defaults(clean_review_artifacts=True)
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    if not (valid_night(args.start) and valid_night(args.end)):
        raise SystemExit("start/end must be YYYYMMDD")
    if args.start > args.end:
        raise SystemExit("start must be <= end")

    processed_root = Path(args.processed_root)
    nights = [night for night in discover_nights(processed_root, args.start, args.end) if night not in set(args.exclude_night)]
    rows: list[dict[str, object]] = []
    for night in nights:
        print(f"[night] {night}", flush=True)
        try:
            row = remask_one(args, night)
        except Exception as exc:
            row = {"night": night, "status": "error", "reason": type(exc).__name__, "missing": str(exc)}
        print(json.dumps(row, ensure_ascii=False, sort_keys=True), flush=True)
        rows.append(row)

    status_out = Path(args.status_out) if args.status_out else (
        Path(args.root_out) / "batch_logs" / f"remask_unknown_with_known_{args.start}_{args.end}_status.tsv"
    )
    write_status(status_out, rows)
    print(f"[write] {status_out}", flush=True)


if __name__ == "__main__":
    main()
