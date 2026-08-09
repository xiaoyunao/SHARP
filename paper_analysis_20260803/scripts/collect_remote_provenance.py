#!/usr/bin/env python3
"""Collect versioned, read-only provenance from the SHARP production host."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from astropy.io import fits


CODE_FILES = [
    "survey/config.py",
    "survey/scheduler.py",
    "survey/run_daily.py",
    "survey/run_daily.sh",
    "known_asteroid/match_single_night.py",
    "known_asteroid/slurm_match_one_file.sh",
    "known_asteroid/export_ades.py",
    "known_asteroid/slurm_merge_submit.sh",
    "heliolincrr/run_daily_pipeline.sh",
    "heliolincrr/run_daily_unknown.sh",
    "heliolincrr/run_single_night.sh",
    "heliolincrr/mask_gaia.py",
    "heliolincrr/make_tracklet_linreproj.py",
    "heliolincrr/run_linear_links_from_tracklets.py",
    "heliolincrr/orbit_confirm_links.py",
    "heliolincrr/summarize_single_night.py",
    "heliolincrr/export_unknown_ades.py",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def iso_mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()


def command(args: list[str], timeout: int = 60) -> dict[str, Any]:
    try:
        result = subprocess.run(
            args,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        return {
            "command": args,
            "returncode": result.returncode,
            "stdout": result.stdout.strip(),
            "stderr": result.stderr.strip(),
        }
    except Exception as exc:
        return {
            "command": args,
            "returncode": None,
            "stdout": "",
            "stderr": f"{type(exc).__name__}: {exc}",
        }


def file_row(path: Path, category: str, hash_contents: bool = True) -> dict[str, Any]:
    row: dict[str, Any] = {
        "category": category,
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": "",
        "mtime_utc": "",
        "sha256": "",
        "hash_status": "not_found",
    }
    if not path.exists():
        return row
    stat = path.stat()
    row.update({"size_bytes": stat.st_size, "mtime_utc": iso_mtime(path)})
    if hash_contents:
        try:
            row.update({"sha256": sha256(path), "hash_status": "ok"})
        except Exception as exc:
            row["hash_status"] = f"error:{type(exc).__name__}:{exc}"
    else:
        row["hash_status"] = "manifest_metadata_only"
    return row


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pipeline-root", type=Path, default=Path("/pipeline/xiaoyunao"))
    parser.add_argument("--gaia-root", type=Path, default=Path("/pipeline/ref/healpix"))
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    pipeline = args.pipeline_root

    file_rows = [file_row(pipeline / rel, "production_code") for rel in CODE_FILES]
    dependency_paths = [
        pipeline / "known_asteroid" / "astorb.dat",
        pipeline / "known_asteroid" / "de432s.bsp",
        pipeline / "survey" / "footprints" / "survey_fov_footprints_with_visibility.fits",
    ]
    file_rows.extend(file_row(path, "dependency") for path in dependency_paths)
    write_csv(out / "production_file_hashes.csv", file_rows)

    footprint = dependency_paths[-1]
    footprint_meta: dict[str, Any] = {"path": str(footprint), "exists": footprint.exists()}
    if footprint.exists():
        with fits.open(footprint, memmap=True) as hdul:
            footprint_meta["n_rows"] = int(len(hdul[1].data))
            footprint_meta["columns"] = list(hdul[1].columns.names)
        shutil.copy2(footprint, out / "survey_fov_footprints_with_visibility.fits")

    gaia_rows: list[dict[str, Any]] = []
    for path in sorted(args.gaia_root.glob("*.fits*")):
        gaia_rows.append(file_row(path, "gaia_healpix_tile", hash_contents=False))
    write_csv(out / "gaia_tile_manifest.csv", gaia_rows)
    gaia_manifest_hash = sha256(out / "gaia_tile_manifest.csv")

    packages = {}
    for name in (
        "numpy",
        "astropy",
        "astro-aleph",
        "pandas",
        "scipy",
        "healpy",
        "pyarrow",
    ):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None

    environment = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python_executable": sys.executable,
        "python_version": sys.version,
        "packages": packages,
        "environment_name": os.environ.get("CONDA_DEFAULT_ENV", ""),
        "conda_list": command([sys.executable, "-m", "pip", "freeze"], timeout=120),
        "short_arc_environment": command(
            [
                "/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python",
                "-c",
                (
                    "import importlib.metadata as m,json,sys;"
                    "names=['numpy','astropy','poliastro','pandas','scipy','pyarrow'];"
                    "print(json.dumps({'python':sys.version,'executable':sys.executable,"
                    "'packages':{n:(m.version(n) if n in {d.metadata['Name'] for d in m.distributions()} else None) for n in names}}))"
                ),
            ],
            timeout=120,
        ),
        "git_pipeline_root": command(["git", "-C", str(pipeline), "rev-parse", "HEAD"]),
        "crontab": command(["crontab", "-l"]),
        "slurm_accounting_probe": command(
            ["sacct", "-S", "2026-07-15", "-E", "2026-07-16", "-n", "-X"],
            timeout=30,
        ),
        "footprint": footprint_meta,
        "gaia_tiles": {
            "root": str(args.gaia_root),
            "count": len(gaia_rows),
            "total_size_bytes": int(sum(int(row["size_bytes"]) for row in gaia_rows)),
            "manifest_sha256": gaia_manifest_hash,
            "release_provenance_in_files": "not_established",
        },
        "model_notes": {
            "known_propagation": "astro_aleph with local de432s.bsp modification",
            "short_arc_orbit": "poliastro environment; astropy solar_system_ephemeris=builtin",
        },
    }
    (out / "remote_environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    hash_lines = []
    for path in sorted(out.iterdir()):
        if path.is_file() and path.name != "hashes.sha256":
            hash_lines.append(f"{sha256(path)}  {path.name}")
    (out / "hashes.sha256").write_text("\n".join(hash_lines) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "production_code_files": len(CODE_FILES),
                "gaia_tiles": len(gaia_rows),
                "footprint_rows": footprint_meta.get("n_rows"),
                "output_dir": str(out),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
