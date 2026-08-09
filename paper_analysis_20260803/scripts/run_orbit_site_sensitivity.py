#!/usr/bin/env python3
"""Re-run frozen short-arc fits with an alternate observer location.

The production module is imported read-only and its three site globals are
overridden in the analysis process.  All outputs are written below a new,
non-production ``--output-root``; existing per-night results are never
overwritten unless ``--resume`` finds that they are already complete and skips
them.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

from astropy.io import fits


def iter_nights(start: str, end: str):
    value = datetime.strptime(start, "%Y-%m-%d").date()
    finish = datetime.strptime(end, "%Y-%m-%d").date()
    while value <= finish:
        yield value.strftime("%Y%m%d")
        value += timedelta(days=1)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def nrows(path: Path) -> int:
    return int(fits.getheader(path, 1).get("NAXIS2", 0)) if path.exists() else -1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-module", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--start", default="2025-11-15")
    parser.add_argument("--end", default="2026-07-15")
    parser.add_argument("--longitude-deg-east", type=float, required=True)
    parser.add_argument("--latitude-deg-north", type=float, required=True)
    parser.add_argument("--height-m", type=float, required=True)
    parser.add_argument("--cores", type=int, default=8)
    parser.add_argument("--limit-nights", type=int)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    source = args.source_module.resolve()
    data_root = args.data_root.resolve()
    output = args.output_root.resolve()
    for production in (source.parent.resolve(), data_root):
        if output == production or production in output.parents or output in production.parents:
            raise ValueError("output root must not overlap a production input")
    output.mkdir(parents=True, exist_ok=True)

    spec = importlib.util.spec_from_file_location("orbit_site_sensitivity_module", source)
    if spec is None or spec.loader is None:
        raise ImportError(source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.MPC327_LON_DEG = float(args.longitude_deg_east)
    module.MPC327_LAT_DEG = float(args.latitude_deg_north)
    module.MPC327_ALT_M = float(args.height_m)

    statuses = []
    candidates = []
    for night in iter_nights(args.start, args.end):
        rr_dir = data_root / night / "rr_links"
        if (rr_dir / "orbit_confirm" / "orbit_links.fits").exists():
            candidates.append((night, rr_dir))
    if args.limit_nights is not None:
        candidates = candidates[: args.limit_nights]

    original_argv = sys.argv[:]
    try:
        for index, (night, rr_dir) in enumerate(candidates, 1):
            night_out = output / night
            links_out = night_out / "orbit_links.fits"
            residual_out = night_out / "orbit_obs_residuals.fits"
            if args.resume and links_out.exists() and residual_out.exists():
                statuses.append(
                    {
                        "night": night,
                        "status": "skipped_complete",
                        "elapsed_s": 0.0,
                        "link_rows": nrows(links_out),
                        "residual_rows": nrows(residual_out),
                    }
                )
                continue
            if night_out.exists():
                raise FileExistsError(f"refusing to overwrite incomplete output: {night_out}")
            started = time.monotonic()
            status = "ok"
            error = ""
            try:
                sys.argv = [
                    str(source),
                    "--rr-dir",
                    str(rr_dir),
                    "--profile",
                    "single-night",
                    "--outdir",
                    str(night_out),
                    "--cores",
                    str(args.cores),
                    "--quiet",
                    "--set-threads-1",
                ]
                module.main()
            except Exception as exc:
                status = "error"
                error = f"{type(exc).__name__}: {exc}"
            elapsed = time.monotonic() - started
            statuses.append(
                {
                    "night": night,
                    "status": status,
                    "error": error,
                    "elapsed_s": elapsed,
                    "link_rows": nrows(links_out),
                    "residual_rows": nrows(residual_out),
                }
            )
            print(
                f"[site] nights={index}/{len(candidates)} night={night} status={status} elapsed_s={elapsed:.2f}",
                flush=True,
            )
    finally:
        sys.argv = original_argv

    fields = sorted({key for row in statuses for key in row})
    with (output / "run_status.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(statuses)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_module": str(source),
        "source_module_sha256": sha256(source),
        "observer": {
            "longitude_deg_east": args.longitude_deg_east,
            "latitude_deg_north": args.latitude_deg_north,
            "height_m": args.height_m,
        },
        "cores": args.cores,
        "candidate_nights": len(candidates),
        "completed_nights": sum(row["status"] in {"ok", "skipped_complete"} for row in statuses),
        "error_nights": sum(row["status"] == "error" for row in statuses),
        "elapsed_s": sum(float(row["elapsed_s"]) for row in statuses),
    }
    (output / "run_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
