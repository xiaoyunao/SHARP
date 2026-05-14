#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from astropy.table import Table


SUBMIT_URL = "https://minorplanetcenter.net/submit_psv"
VALIDATE_URL = "https://www.minorplanetcenter.net/submit_psv_test"
TRK_SUB_RE = re.compile(r"^[-A-Za-z0-9_]{1,8}$")


def split_semicolon(value: object) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def mjd_to_obs_time_z(mjd: float) -> str:
    if not np.isfinite(float(mjd)):
        return ""
    mjd0 = datetime(1858, 11, 17, tzinfo=timezone.utc)
    return (mjd0 + timedelta(days=float(mjd))).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def validate_trk_sub(value: object) -> str:
    trk_sub = str(value or "").strip()
    if not trk_sub:
        return ""
    if not TRK_SUB_RE.fullmatch(trk_sub):
        raise ValueError(
            f"Invalid trkSub {trk_sub!r}: must be 1-8 characters using A-Z, a-z, 0-9, '_' or '-'."
        )
    return trk_sub


def read_review_csv(path: str | None, require_review: bool) -> tuple[dict[str, int], bool]:
    if not path:
        return {}, require_review
    review_path = Path(path)
    text = review_path.read_text(encoding="utf-8")
    sample = text[:2048]
    dialect = csv.Sniffer().sniff(sample, delimiters=",\t|") if sample.strip() else csv.excel
    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    out: dict[str, int] = {}
    for row in rows:
        key = (
            row.get("trk_sub")
            or row.get("trkSub")
            or row.get("tracklet_id")
            or row.get("linkage_id")
            or ""
        )
        key = str(key).strip()
        if not key:
            continue
        value = row.get("is_real") or row.get("real") or row.get("truth") or row.get("valid") or ""
        try:
            out[key] = 1 if int(str(value).strip()) == 1 else 0
        except ValueError as exc:
            raise ValueError(f"Invalid review value for {key!r}: {value!r}; expected 0 or 1") from exc
    return out, require_review


def build_header(args: argparse.Namespace) -> str:
    lines = ["# version=2017"]
    lines.append("# observatory")
    lines.append(f"! mpcCode {args.obs_code}")
    if args.obs_name:
        lines.append(f"! name {args.obs_name}")

    lines.append("# submitter")
    lines.append(f"! name {args.submitter_name}")
    if args.institution:
        lines.append(f"! institution {args.institution}")

    for block_name, values in (
        ("observers", args.observers),
        ("measurers", args.measurers),
        ("coinvestigators", args.coinvestigators),
        ("collaborators", args.collaborators),
    ):
        names = split_semicolon(values)
        if names:
            lines.append(f"# {block_name}")
            for name in names:
                lines.append(f"! name {name}")

    telescope_lines = []
    for key, value in (
        ("name", args.telescope_name),
        ("aperture", args.aperture),
        ("design", args.design),
        ("detector", args.detector),
        ("fRatio", args.f_ratio),
        ("filter", args.filter_name),
        ("arraySize", args.array_size),
        ("pixelScale", args.pixel_scale),
    ):
        if value:
            telescope_lines.append(f"! {key} {value}")
    if telescope_lines:
        lines.append("# telescope")
        lines.extend(telescope_lines)

    if args.funding_source:
        lines.append(f"# fundingSource {args.funding_source}")
    comments = split_semicolon(args.comment)
    if comments:
        lines.append("# comment")
        for line in comments:
            lines.append(f"! line {line}")
    lines.append("")
    return "\n".join(lines)


def load_unknown_rows(path: Path) -> list[dict[str, object]]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError(f"Unknown catalog must contain a JSON list: {path}")
    return rows


def should_keep(row: dict[str, object], review: dict[str, int], require_review: bool) -> bool:
    if not review and not require_review:
        return True
    keys = [
        str(row.get("trk_sub", "")).strip(),
        str(row.get("trkSub", "")).strip(),
        str(row.get("linkage_id", "")).strip(),
    ]
    keys.extend(split_semicolon(row.get("tracklet_ids", "")))
    keys = [key for key in keys if key]
    for key in keys:
        if key in review:
            return review[key] == 1
    return not require_review


def load_detection_row(processed_root: Path, night: str, image_name: str, obj_id: int) -> dict[str, float] | None:
    cat_name = image_name.replace(".fits.gz", "_cat.fits.gz")
    path = processed_root / night / "L2" / cat_name
    if not path.exists():
        return None
    table = Table.read(path)
    if "objID" not in table.colnames:
        return None
    match = table[np.asarray(table["objID"], dtype=np.int64) == int(obj_id)]
    if len(match) == 0:
        return None
    row = match[0]
    required = ["MJD", "RA_Win", "DEC_Win", "RAErr_Win", "DECErr_Win", "Mag_Aper4", "MagErr_Aper4"]
    if any(col not in table.colnames for col in required):
        return None
    return {col: float(row[col]) for col in required}


def build_obs_rows(args: argparse.Namespace) -> tuple[list[dict[str, object]], dict[str, int]]:
    catalog_path = Path(args.catalog)
    processed_root = Path(args.processed_root)
    review, require_review = read_review_csv(args.review_csv, bool(args.require_review))
    rows = load_unknown_rows(catalog_path)

    obs_rows: list[dict[str, object]] = []
    stats = {
        "catalog_links": len(rows),
        "review_rejected_links": 0,
        "missing_trk_sub_links": 0,
        "missing_detection_rows": 0,
        "exported_observations": 0,
    }

    for row in rows:
        if not should_keep(row, review, require_review):
            stats["review_rejected_links"] += 1
            continue
        trk_sub = validate_trk_sub(row.get("trk_sub") or row.get("trkSub"))
        if not trk_sub:
            stats["missing_trk_sub_links"] += 1
            continue
        image_names = split_semicolon(row.get("image_names"))
        objids = [int(x) for x in split_semicolon(row.get("objids"))]
        fallback_mjds = [float(x) for x in split_semicolon(row.get("mjds"))]
        fallback_ras = [float(x) for x in split_semicolon(row.get("ras_deg"))]
        fallback_decs = [float(x) for x in split_semicolon(row.get("decs_deg"))]
        for idx, image_name in enumerate(image_names):
            if idx >= len(objids):
                stats["missing_detection_rows"] += 1
                continue
            det = load_detection_row(processed_root, args.night, image_name, objids[idx])
            if det is None:
                stats["missing_detection_rows"] += 1
                continue
            mjd = det.get("MJD", fallback_mjds[idx] if idx < len(fallback_mjds) else np.nan)
            ra = det.get("RA_Win", fallback_ras[idx] if idx < len(fallback_ras) else np.nan)
            dec = det.get("DEC_Win", fallback_decs[idx] if idx < len(fallback_decs) else np.nan)
            rms_ra = det["RAErr_Win"]
            rms_dec = det["DECErr_Win"]
            if args.err_unit == "deg":
                rms_ra *= 3600.0
                rms_dec *= 3600.0
            elif args.err_unit == "mas":
                rms_ra /= 1000.0
                rms_dec /= 1000.0
            elif args.err_unit != "arcsec":
                raise ValueError("err-unit must be one of: arcsec, deg, mas")
            mag = det["Mag_Aper4"]
            rms_mag = det["MagErr_Aper4"]
            finite = np.all(np.isfinite([mjd, ra, dec, rms_ra, rms_dec, mag, rms_mag]))
            if not finite or rms_mag <= 0:
                stats["missing_detection_rows"] += 1
                continue
            obs_rows.append(
                {
                    "trkSub": trk_sub,
                    "mode": args.mode,
                    "stn": args.obs_code,
                    "obsTime": mjd_to_obs_time_z(mjd),
                    "ra": f"{ra:.10f}",
                    "dec": f"{dec:.10f}",
                    "rmsRA": f"{rms_ra:.4f}",
                    "rmsDec": f"{rms_dec:.4f}",
                    "astCat": args.astcat,
                    "mag": f"{mag:.3f}",
                    "rmsMag": f"{rms_mag:.3f}",
                    "band": args.band,
                    "photCat": args.photcat,
                }
            )
    stats["exported_observations"] = len(obs_rows)
    return obs_rows, stats


def write_psv(args: argparse.Namespace, obs_rows: list[dict[str, object]]) -> None:
    cols = ["trkSub", "mode", "stn", "obsTime", "ra", "dec", "rmsRA", "rmsDec", "astCat", "mag", "rmsMag", "band", "photCat"]
    lines = [build_header(args).rstrip("\n"), " | ".join(cols)]
    for row in obs_rows:
        lines.append("|".join(str(row[col]) for col in cols))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_curl(url: str, args: argparse.Namespace) -> None:
    cmd = [
        "curl",
        "-sS",
        "-L",
        url,
        "-F",
        f"ack={args.ack}",
        "-F",
        f"ac2={args.ac2_email}",
        "-F",
        f"obj_type={args.obj_type}",
        "-F",
        f"source=@{args.out}",
    ]
    print("[RUN]", " ".join(cmd))
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        if proc.stderr.strip():
            print("[CURL-ERR]", proc.stderr.strip())
        raise RuntimeError(f"curl failed with return code {proc.returncode}")
    reply = proc.stdout.strip()
    print("[MPC-REPLY]")
    print(reply)
    if args.response_out:
        Path(args.response_out).write_text(reply + "\n", encoding="utf-8")
        print(f"[OK] Wrote reply: {args.response_out}")


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Export fit_ok unknown single-night links to ADES PSV using trkSub.")
    ap.add_argument("night", help="Target night in YYYYMMDD format")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--catalog", default=None, help="Default: /processed1/<night>/L4/<night>_unknown_links.json")
    ap.add_argument("--out", default=None, help="Default: /processed1/<night>/L4/<night>_unknown_links_ades.psv")
    ap.add_argument("--review-csv", default="", help="Optional CSV with tracklet_id/trk_sub/linkage_id and is_real columns")
    ap.add_argument("--require-review", action="store_true", help="Only export links explicitly marked is_real=1")

    ap.add_argument("--mode", default="CCD")
    ap.add_argument("--obs-code", default="327")
    ap.add_argument("--obs-name", default="Xinglong Station")
    ap.add_argument("--submitter-name", default="Y.-A. Xiao")
    ap.add_argument("--institution", default="NAOC")
    ap.add_argument("--observers", default="Xiangnan Guan;Pengfei Liu;Xiaoming Teng;WSGP Team")
    ap.add_argument("--measurers", default="Niu Li;Yun-Ao Xiao;Jingyi Zhang;Hu Zou;WSGP Team")
    ap.add_argument("--comment", default="")
    ap.add_argument("--telescope-name", default="60/90cm Schmidt telescope")
    ap.add_argument("--aperture", default="0.6")
    ap.add_argument("--design", default="Schmidt")
    ap.add_argument("--detector", default="CCD")
    ap.add_argument("--f-ratio", dest="f_ratio", default="3")
    ap.add_argument("--filter-name", default="unfiltered")
    ap.add_argument("--array-size", default="9216 x 9232")
    ap.add_argument("--pixel-scale", default="1.15")
    ap.add_argument("--coinvestigators", default="Hu Zou")
    ap.add_argument("--collaborators", default="Yun-Ao Xiao;Niu Li;Jingyi Zhang;Wenxiong Li;Zhaobin Chen;Shufei Liu;WSGP Team")
    ap.add_argument(
        "--funding-source",
        default="The 60/90-cm Schmidt telescope situated at the Xinglong station is overseen by the research team of the Wide-field Survey and Galaxy Physics (WSGP) at NAOC.",
    )
    ap.add_argument("--astcat", default="Gaia3E")
    ap.add_argument("--band", default="G")
    ap.add_argument("--photcat", default="Gaia3E")
    ap.add_argument("--err-unit", default="deg", choices=["arcsec", "deg", "mas"])
    ap.add_argument("--ac2-email", default="wsgp2024@163.com")
    ap.add_argument("--ack", default="Unknown-object ADES submission")
    ap.add_argument("--obj-type", default="NEO")
    ap.add_argument("--response-out", default="")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--submit", action="store_true")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    if args.catalog is None:
        args.catalog = str(Path(args.processed_root) / args.night / "L4" / f"{args.night}_unknown_links.json")
    if args.out is None:
        args.out = str(Path(args.processed_root) / args.night / "L4" / f"{args.night}_unknown_links_ades.psv")

    obs_rows, stats = build_obs_rows(args)
    write_psv(args, obs_rows)
    print(f"[OK] Wrote PSV: {args.out}")
    print(json.dumps(stats, indent=2, sort_keys=True))
    if not obs_rows:
        print("[SKIP] No unknown ADES obsData rows exported; skip validate/submit")
        return
    if args.validate:
        run_curl(VALIDATE_URL, args)
    if args.submit:
        run_curl(SUBMIT_URL, args)


if __name__ == "__main__":
    main()
