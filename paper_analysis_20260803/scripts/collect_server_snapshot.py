#!/usr/bin/env python3
"""Collect a compact, read-only inventory of the SHARP production snapshot.

The script is designed to run on the production server.  It reads FITS headers,
small JSON/CSV/text products, file metadata, and archived scheduler plans.  It
never writes beneath any production root; every output goes to ``--output-dir``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import socket
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable

from astropy.io import fits


STRICT_RAW_RE = re.compile(r"^OBJ_MP_(\d{4})_(\d{4})\.fits$")
STRICT_L2_RE = re.compile(r"^OBJ_MP_(\d{4})_(\d{4})_cat\.fits(?:\.gz)?$")
ANY_MP_L2_RE = re.compile(r".*MP.*\.(?:fits|fits\.gz)$")
SUBMISSION_ID_RE = re.compile(r"Submission ID is\s+([^\s]+)", re.IGNORECASE)


def iso_mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()


def iter_nights(start: str, end: str) -> Iterable[str]:
    current = datetime.strptime(start, "%Y-%m-%d").date()
    finish = datetime.strptime(end, "%Y-%m-%d").date()
    while current <= finish:
        yield current.strftime("%Y%m%d")
        current += timedelta(days=1)


def safe_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, (str, int, float, bool)):
        return value
    try:
        return value.item()
    except Exception:
        return str(value)


def header_value(header: fits.Header, *keys: str) -> Any:
    for key in keys:
        if key in header and header[key] not in (None, ""):
            return safe_value(header[key])
    return ""


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fits_nrows(path: Path, hdu: int = 1) -> int:
    if not path.exists():
        return 0
    try:
        return int(fits.getheader(path, hdu).get("NAXIS2", 0))
    except Exception:
        return -1


def collect_raw(
    raw_root: Path,
    nights: list[str],
    read_headers: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    strict_rows: list[dict[str, Any]] = []
    engineering_rows: list[dict[str, Any]] = []
    for index, night in enumerate(nights, 1):
        obj_dir = raw_root / night / "OBJ"
        if not obj_dir.is_dir():
            continue
        for path in sorted(obj_dir.iterdir()):
            if not path.is_file() or "MP" not in path.name or path.suffix.lower() != ".fits":
                continue
            match = STRICT_RAW_RE.fullmatch(path.name)
            row: dict[str, Any] = {
                "night": night,
                "file_name": path.name,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "mtime_utc": iso_mtime(path),
                "strict_standard_science": bool(match),
                "field_id": match.group(1) if match else "",
                "sequence_id": match.group(2) if match else "",
                "header_read_status": "not_requested" if not read_headers else "pending",
                "header_ok": "",
                "header_error": "",
                "exposure_s": 30.0 if match else "",
            }
            if read_headers:
                try:
                    header = fits.getheader(path, 0)
                    row.update(
                        {
                            "header_read_status": "ok",
                            "header_ok": True,
                            "obs_time_utc": header_value(header, "OBS_DATE", "DATE-OBS", "DATEOBS", "EXPSTA"),
                            "exposure_s": header_value(header, "EXPOSURE", "EXPREQ", "EXPTIME"),
                            "object": header_value(header, "OBJECT"),
                            "obj_ra_header": header_value(header, "OBJ_RA"),
                            "obj_dec_deg": header_value(header, "OBJ_DEC"),
                            "tel_ra_header": header_value(header, "TEL_RA"),
                            "tel_dec_deg": header_value(header, "TEL_DEC"),
                            "tel_elevation_deg": header_value(header, "TEL_EL"),
                            "filter": header_value(header, "FILTER", "BAND"),
                            "telescope": header_value(header, "TELESCOP"),
                            "observatory_header": header_value(header, "OBSERVAT"),
                            "obslon_header_deg": header_value(header, "OBSLON"),
                            "obslat_header_deg": header_value(header, "OBSLAT"),
                            "obsalt_header_m": header_value(header, "OBSALT"),
                        }
                    )
                except Exception as exc:
                    row["header_read_status"] = "error"
                    row["header_ok"] = False
                    row["header_error"] = f"{type(exc).__name__}: {exc}"
            engineering_rows.append(row)
            if match:
                strict_rows.append(dict(row))
        if index % 10 == 0:
            print(f"[raw] nights={index}/{len(nights)} strict_rows={len(strict_rows)}", flush=True)
    return strict_rows, engineering_rows


L2_HEADER_FIELDS = {
    "obs_time_utc": ("OBS_DATE", "DATE-OBS", "DATEOBS", "EXPSTA"),
    "mjd": ("MJD",),
    "exposure_s": ("EXPOSURE", "EXPREQ", "EXPTIME"),
    "filter": ("FILTER", "BAND"),
    "center_ra_deg": ("CEN_RA", "CRVAL1"),
    "center_dec_deg": ("CEN_DEC", "CRVAL2"),
    "crval1_deg": ("CRVAL1",),
    "crval2_deg": ("CRVAL2",),
    "crpix1_px": ("CRPIX1",),
    "crpix2_px": ("CRPIX2",),
    "ctype1": ("CTYPE1",),
    "ctype2": ("CTYPE2",),
    "cunit1": ("CUNIT1",),
    "cunit2": ("CUNIT2",),
    "cd1_1": ("CD1_1",),
    "cd1_2": ("CD1_2",),
    "cd2_1": ("CD2_1",),
    "cd2_2": ("CD2_2",),
    "pc1_1": ("PC1_1",),
    "pc1_2": ("PC1_2",),
    "pc2_1": ("PC2_1",),
    "pc2_2": ("PC2_2",),
    "cdelt1": ("CDELT1",),
    "cdelt2": ("CDELT2",),
    "seeing_arcsec": ("SEEING",),
    "fwhm_header": ("FWHM",),
    "sky_adu": ("SKY",),
    "sky_rms_adu": ("SKYRMS",),
    "sky_mag": ("SKYMAG",),
    "limiting_mag": ("MLIM",),
    "ccd_zeropoint": ("CCDZP",),
    "ccd_zeropoint_rms": ("CCDZPRMS",),
    "ra_offset_header": ("RA_OFF",),
    "dec_offset_header": ("DEC_OFF",),
    "ra_rms_header": ("RA_RMS",),
    "dec_rms_header": ("DEC_RMS",),
    "n_wcs": ("N_WCS",),
    "n_stars": ("NSTARS",),
    "n_match_calibration": ("NMATCH",),
    "n_psf": ("NS_PSF",),
    "photometric_rms": ("PHOTIRMS",),
    "photometric_link": ("PHOTLINK",),
    "astrometric_catalog_header": ("CALI_REF",),
    "telescope": ("TELESCOP",),
    "observatory_header": ("OBSERVAT",),
    "obslon_header_deg": ("OBSLON",),
    "obslat_header_deg": ("OBSLAT",),
    "obsalt_header_m": ("OBSALT",),
}


def collect_l2_path(path: Path, night: str) -> dict[str, Any]:
    match = STRICT_L2_RE.fullmatch(path.name)
    row: dict[str, Any] = {
        "night": night,
        "file_name": path.name,
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "mtime_utc": iso_mtime(path),
        "strict_standard_catalog": bool(match),
        "field_id": match.group(1) if match else "",
        "sequence_id": match.group(2) if match else "",
        "header_ok": False,
        "header_error": "",
    }
    try:
        header = fits.getheader(path, 1)
        row["header_ok"] = True
        row["n_sources"] = int(header.get("NAXIS2", 0))
        row["row_width_bytes"] = int(header.get("NAXIS1", 0))
        for out_name, keys in L2_HEADER_FIELDS.items():
            row[out_name] = header_value(header, *keys)
    except Exception as exc:
        row["header_error"] = f"{type(exc).__name__}: {exc}"
    return row


def collect_l2(processed_root: Path, nights: list[str], workers: int = 16) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    inputs: list[tuple[Path, str]] = []
    for night in nights:
        l2_dir = processed_root / night / "L2"
        if l2_dir.is_dir():
            inputs.extend(
                (path, night)
                for path in sorted(l2_dir.iterdir())
                if path.is_file() and ANY_MP_L2_RE.fullmatch(path.name)
            )
    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        for index, row in enumerate(
            executor.map(lambda item: collect_l2_path(*item), inputs), 1
        ):
            rows.append(row)
            if index % 1000 == 0:
                print(f"[l2] files={index}/{len(inputs)}", flush=True)
    return rows


def collect_plans(plan_root: Path, start: str, end: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    files: list[dict[str, Any]] = []
    start_key = start.replace("-", "")
    end_key = end.replace("-", "")
    for path in sorted(plan_root.glob("*_plan.json")):
        night = path.name[:8]
        if not (start_key <= night <= end_key):
            continue
        record = {
            "night": night,
            "path": str(path),
            "size_bytes": path.stat().st_size,
            "mtime_utc": iso_mtime(path),
            "json_ok": False,
            "n_rows": 0,
            "error": "",
        }
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(payload, list):
                raise TypeError(f"expected list, got {type(payload).__name__}")
            record["json_ok"] = True
            record["n_rows"] = len(payload)
            for row_index, item in enumerate(payload):
                out = {"night": night, "plan_row": row_index, "plan_file_mtime_utc": record["mtime_utc"]}
                out.update({str(key): safe_value(value) for key, value in item.items()})
                rows.append(out)
        except Exception as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
        files.append(record)
    return rows, files


def count_ades_rows(path: Path) -> int:
    if not path.exists():
        return 0
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith(("#", "!")) or text.startswith("permID") or text.startswith("trkSub"):
                continue
            if "|" in text:
                count += 1
    return count


def reply_info(path: Path) -> tuple[str, str]:
    if not path.exists():
        return "", ""
    text = path.read_text(encoding="utf-8", errors="replace").strip()
    match = SUBMISSION_ID_RE.search(text)
    return (match.group(1) if match else "", text.replace("\n", " ")[:1000])


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        return payload if isinstance(payload, dict) else None
    except Exception:
        return None


def parse_review_csv(path: Path) -> tuple[int, int, int, int]:
    if not path.exists():
        return 0, 0, 0, 0
    total = real = false = blank = 0
    with path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        for row in csv.DictReader(handle):
            total += 1
            value = str(row.get("is_real", "")).strip().lower()
            if value in {"1", "true", "t", "yes", "y"}:
                real += 1
            elif value in {"0", "false", "f", "no", "n"}:
                false += 1
            else:
                blank += 1
    return total, real, false, blank


def collect_nightly_products(
    processed_root: Path,
    heliolinc_data_root: Path,
    review_root: Path,
    nights: list[str],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for night in nights:
        l4 = processed_root / night / "L4"
        analysis = heliolinc_data_root / night / "analysis"
        review_dir = review_root / night
        paths = {
            "known_all": l4 / f"{night}_all_asteroids.fits",
            "known_match1": l4 / f"{night}_matched_asteroids.fits",
            "known_mask15": l4 / f"{night}_matched_asteroids_mask15.fits",
            "known_ades": l4 / f"{night}_matched_asteroids_ades.psv",
            "known_reply": l4 / f"{night}_mpc_reply.txt",
            "known_status": l4 / f"{night}_known_asteroid_status.json",
            "unknown_links": l4 / f"{night}_unknown_links.fits",
            "unknown_json": l4 / f"{night}_unknown_links.json",
            "unknown_ades": l4 / f"{night}_unknown_links_submit_ades.psv",
            "unknown_reply": l4 / f"{night}_unknown_mpc_reply.txt",
            "unknown_validate_reply": l4 / f"{night}_unknown_validate_reply.txt",
            "summary": analysis / f"{night}_single_night_summary.json",
            "review_manifest": review_dir / f"{night}_unknown_review_manifest.json",
            "review_csv": review_dir / f"{night}_unknown_review.csv",
            "submit_csv": review_dir / f"{night}_submit.csv",
        }
        if not any(path.exists() for path in paths.values()):
            continue
        row: dict[str, Any] = {"night": night}
        for name in ("known_all", "known_match1", "known_mask15", "unknown_links"):
            path = paths[name]
            row[f"{name}_exists"] = path.exists()
            row[f"{name}_rows"] = fits_nrows(path) if path.exists() else 0
            row[f"{name}_mtime_utc"] = iso_mtime(path) if path.exists() else ""
        for name in ("known_ades", "unknown_ades"):
            path = paths[name]
            row[f"{name}_exists"] = path.exists()
            row[f"{name}_rows"] = count_ades_rows(path)
            row[f"{name}_mtime_utc"] = iso_mtime(path) if path.exists() else ""
        for prefix, name in (("known", "known_reply"), ("unknown", "unknown_reply"), ("unknown_validate", "unknown_validate_reply")):
            path = paths[name]
            submission_id, text = reply_info(path)
            row[f"{prefix}_reply_exists"] = path.exists()
            row[f"{prefix}_submission_id"] = submission_id
            row[f"{prefix}_reply_text"] = text
            row[f"{prefix}_reply_mtime_utc"] = iso_mtime(path) if path.exists() else ""
        for prefix, name in (("review", "review_csv"), ("submit", "submit_csv")):
            total, real, false, blank = parse_review_csv(paths[name])
            row[f"{prefix}_csv_exists"] = paths[name].exists()
            row[f"{prefix}_rows"] = total
            row[f"{prefix}_real"] = real
            row[f"{prefix}_false"] = false
            row[f"{prefix}_blank"] = blank
            row[f"{prefix}_mtime_utc"] = iso_mtime(paths[name]) if paths[name].exists() else ""
        summary = read_json(paths["summary"])
        row["summary_exists"] = summary is not None
        row["summary_mtime_utc"] = iso_mtime(paths["summary"]) if paths["summary"].exists() else ""
        if summary:
            counts = summary.get("counts", {})
            for key, value in counts.items():
                if not isinstance(value, (dict, list)):
                    row[f"summary_{key}"] = safe_value(value)
            for group in ("link_class_counts", "fit_ok_link_class_counts", "is_good_link_class_counts"):
                for key, value in summary.get(group, {}).items():
                    row[f"{group}_{key}"] = safe_value(value)
            requested = summary.get("requested_metrics", {})
            for stage, metrics in requested.items():
                if isinstance(metrics, dict):
                    for key, value in metrics.items():
                        if not isinstance(value, (dict, list)):
                            row[f"requested_{stage}_{key}"] = safe_value(value)
        manifest = read_json(paths["review_manifest"])
        row["review_manifest_exists"] = manifest is not None
        row["review_manifest_mtime_utc"] = iso_mtime(paths["review_manifest"]) if paths["review_manifest"].exists() else ""
        if manifest:
            for key in ("n_catalog_rows", "n_gifs_copied", "n_gifs_missing", "review_full_rows", "review_ades_rows"):
                row[f"review_manifest_{key}"] = safe_value(manifest.get(key, ""))
        status = read_json(paths["known_status"])
        row["known_status_exists"] = status is not None
        if status:
            row["known_status_complete"] = safe_value(status.get("known_complete", ""))
        rows.append(row)
    return rows


def collect_mask_gaia_logs(heliolinc_data_root: Path, nights: list[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for night in nights:
        path = heliolinc_data_root / night / "mask_gaia" / "mask_gaia.log"
        if not path.exists():
            continue
        record: dict[str, Any] = {
            "night": night,
            "path": str(path),
            "mtime_utc": iso_mtime(path),
            "parsed_rows": 0,
            "prefilter_rows": 0,
            "gaia_survivor_rows": 0,
            "gaia_matched_rows": 0,
            "status_ok_files": 0,
            "status_counts": "",
        }
        statuses: Counter[str] = Counter()
        in_table = False
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("status\t"):
                in_table = True
                continue
            if not in_table or not line or line.startswith("-"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            status = parts[0]
            try:
                n_in, n_out, n_match = int(parts[1]), int(parts[2]), int(parts[4])
            except ValueError:
                continue
            statuses[status] += 1
            record["parsed_rows"] += 1
            record["prefilter_rows"] += n_in
            record["gaia_survivor_rows"] += n_out
            record["gaia_matched_rows"] += n_match
        record["status_ok_files"] = statuses.get("ok", 0)
        record["status_counts"] = json.dumps(dict(sorted(statuses.items())), sort_keys=True)
        rows.append(record)
    return rows


def file_metadata_rows(paths: Iterable[Path], category: str) -> list[dict[str, Any]]:
    rows = []
    for path in sorted({path for path in paths if path.exists()}):
        rows.append(
            {
                "category": category,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "mtime_utc": iso_mtime(path),
            }
        )
    return rows


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_output_hashes(output_dir: Path) -> None:
    rows = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file() or path.name == "hashes.sha256":
            continue
        rows.append((sha256_file(path), path.relative_to(output_dir).as_posix()))
    with (output_dir / "hashes.sha256").open("w", encoding="utf-8") as handle:
        for digest, relative in rows:
            handle.write(f"{digest}  {relative}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--start", default="2025-11-15")
    parser.add_argument("--end", default="2026-07-15")
    parser.add_argument("--plan-end", default="2026-08-03")
    parser.add_argument("--raw-root", default="/raw1")
    parser.add_argument("--processed-root", default="/processed1")
    parser.add_argument("--pipeline-root", default="/pipeline/xiaoyunao")
    parser.add_argument(
        "--read-raw-headers",
        action="store_true",
        help="Read every raw FITS header (slow on the production NFS mount).",
    )
    parser.add_argument("--l2-workers", type=int, default=16)
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    nights = list(iter_nights(args.start, args.end))
    pipeline_root = Path(args.pipeline_root)
    processed_root = Path(args.processed_root)
    heliolinc_data_root = pipeline_root / "data" / "heliolincrr"
    review_root = pipeline_root / "heliolincrr" / "review_packages"

    run_meta = {
        "collector_schema_version": "1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version,
        "command": " ".join(sys.argv),
        "start": args.start,
        "end": args.end,
        "plan_end": args.plan_end,
        "roots": {
            "raw": str(Path(args.raw_root)),
            "processed": str(processed_root),
            "pipeline": str(pipeline_root),
        },
    }
    (output_dir / "collector_run.json").write_text(json.dumps(run_meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    strict_raw, engineering_raw = collect_raw(
        Path(args.raw_root), nights, read_headers=args.read_raw_headers
    )
    write_csv(output_dir / "raw_manifest.csv", strict_raw)
    write_csv(output_dir / "raw_engineering_manifest.csv", engineering_raw)

    l2_rows = collect_l2(processed_root, nights, workers=args.l2_workers)
    write_csv(output_dir / "l2_manifest.csv", l2_rows)

    plan_rows, plan_files = collect_plans(
        pipeline_root / "survey" / "runtime" / "plans",
        args.start,
        args.plan_end,
    )
    write_csv(output_dir / "plan_rows.csv", plan_rows)
    write_csv(output_dir / "plan_files.csv", plan_files)

    nightly = collect_nightly_products(processed_root, heliolinc_data_root, review_root, nights)
    write_csv(output_dir / "nightly_products.csv", nightly)

    mask_logs = collect_mask_gaia_logs(heliolinc_data_root, nights)
    write_csv(output_dir / "mask_gaia_stage_counts.csv", mask_logs)

    source_paths = []
    source_paths.extend((processed_root / row["night"] / "L4" / f"{row['night']}_known_asteroid_status.json") for row in nightly)
    source_paths.extend((heliolinc_data_root / row["night"] / "analysis" / f"{row['night']}_single_night_summary.json") for row in nightly)
    source_paths.extend((review_root / row["night"] / f"{row['night']}_unknown_review_manifest.json") for row in nightly)
    source_metadata = file_metadata_rows(source_paths, "small_production_product")
    write_csv(output_dir / "source_product_metadata.csv", source_metadata)

    summary = {
        "strict_raw_rows": len(strict_raw),
        "strict_raw_nights": len({row["night"] for row in strict_raw}),
        "all_mp_fits_rows": len(engineering_raw),
        "all_mp_fits_nights": len({row["night"] for row in engineering_raw}),
        "l2_mp_rows": len(l2_rows),
        "l2_mp_nights": len({row["night"] for row in l2_rows}),
        "plan_rows": len(plan_rows),
        "plan_nights": len({row["night"] for row in plan_rows}),
        "nightly_product_rows": len(nightly),
        "mask_gaia_log_nights": len(mask_logs),
    }
    (output_dir / "collector_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_output_hashes(output_dir)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
