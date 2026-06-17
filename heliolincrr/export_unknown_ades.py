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


def read_review_csv(path: str | None, require_values: bool = False) -> tuple[dict[str, int], dict[str, object]]:
    meta: dict[str, object] = {
        "path": str(path or ""),
        "rows": 0,
        "decisions": 0,
        "blank_values": 0,
        "duplicate_keys": 0,
        "blank_samples": [],
    }
    if not path:
        return {}, meta
    review_path = Path(path)
    text = review_path.read_text(encoding="utf-8")
    sample = text[:2048]
    dialect = csv.Sniffer().sniff(sample, delimiters=",\t|") if sample.strip() else csv.excel
    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    out: dict[str, int] = {}
    for row in rows:
        meta["rows"] = int(meta["rows"]) + 1
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
        value_text = str(value).strip()
        if value_text == "":
            meta["blank_values"] = int(meta["blank_values"]) + 1
            samples = meta["blank_samples"]
            if isinstance(samples, list) and len(samples) < 10:
                samples.append(key)
            continue
        try:
            decision = int(value_text)
        except ValueError as exc:
            raise ValueError(f"Invalid review value for {key!r}: {value!r}; expected 0 or 1") from exc
        if decision not in (0, 1):
            raise ValueError(f"Invalid review value for {key!r}: {value!r}; expected 0 or 1")
        if key in out:
            meta["duplicate_keys"] = int(meta["duplicate_keys"]) + 1
            if out[key] != decision:
                raise ValueError(f"Conflicting duplicate review decision for {key!r}")
        out[key] = decision
    meta["decisions"] = len(out)
    if require_values and int(meta["blank_values"]) > 0:
        raise ValueError(
            f"Review CSV {review_path} has {meta['blank_values']} blank is_real values; "
            f"sample keys: {meta['blank_samples']}"
        )
    return out, meta


def build_header(args: argparse.Namespace) -> str:
    lines = ["# version=2017"]

    lines.append("# observatory")
    lines.append(f"! mpcCode {args.obs_code}")
    if args.obs_name:
        lines.append(f"! name {args.obs_name}")

    lines.append("# submitter")
    if args.submitter_name:
        lines.append(f"! name {args.submitter_name}")
    if args.institution:
        lines.append(f"! institution {args.institution}")

    obs_list = split_semicolon(args.observers)
    if obs_list:
        lines.append("# observers")
        for name in obs_list:
            lines.append(f"! name {name}")

    meas_list = split_semicolon(args.measurers)
    if meas_list:
        lines.append("# measurers")
        for name in meas_list:
            lines.append(f"! name {name}")

    telescope_lines: list[str] = []
    if args.telescope_name:
        telescope_lines.append(f"! name {args.telescope_name}")
    if args.aperture:
        telescope_lines.append(f"! aperture {args.aperture}")
    if args.design:
        telescope_lines.append(f"! design {args.design}")
    if args.detector:
        telescope_lines.append(f"! detector {args.detector}")
    if args.f_ratio:
        telescope_lines.append(f"! fRatio {args.f_ratio}")
    if args.filter_name:
        telescope_lines.append(f"! filter {args.filter_name}")
    if args.array_size:
        telescope_lines.append(f"! arraySize {args.array_size}")
    if args.pixel_scale:
        telescope_lines.append(f"! pixelScale {args.pixel_scale}")
    if telescope_lines:
        lines.append("# telescope")
        lines.extend(telescope_lines)

    coinvestigator_list = split_semicolon(args.coinvestigators)
    if coinvestigator_list:
        lines.append("# coinvestigators")
        for name in coinvestigator_list:
            lines.append(f"! name {name}")

    collaborator_list = split_semicolon(args.collaborators)
    if collaborator_list:
        lines.append("# collaborators")
        for name in collaborator_list:
            lines.append(f"! name {name}")

    if args.funding_source:
        lines.append(f"# fundingSource {args.funding_source}")

    comment_list = split_semicolon(args.comment)
    if comment_list:
        lines.append("# comment")
        for line in comment_list:
            lines.append(f"! line {line}")

    lines.append("")
    return "\n".join(lines)


def compute_logsnr(flux: float | None, flux_err: float | None) -> str:
    if flux is None or flux_err is None:
        return ""
    if not np.isfinite(flux) or not np.isfinite(flux_err) or flux_err == 0:
        return ""
    snr = flux / flux_err
    if not np.isfinite(snr) or snr <= 0:
        return ""
    return f"{np.log10(snr):.2f}"


def load_unknown_rows(path: Path) -> list[dict[str, object]]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError(f"Unknown catalog must contain a JSON list: {path}")
    return rows


def review_keys_for_row(row: dict[str, object]) -> list[str]:
    keys = [
        str(row.get("trk_sub", "")).strip(),
        str(row.get("trkSub", "")).strip(),
        str(row.get("linkage_id", "")).strip(),
    ]
    keys.extend(split_semicolon(row.get("tracklet_ids", "")))
    return [key for key in keys if key]


def review_decision_for_row(row: dict[str, object], review: dict[str, int]) -> int | None:
    for key in review_keys_for_row(row):
        if key in review:
            return review[key]
    return None


def should_keep(row: dict[str, object], review: dict[str, int], require_review: bool) -> bool:
    if not review and not require_review:
        return True
    decision = review_decision_for_row(row, review)
    if decision is None:
        return not require_review
    return decision == 1


def filter_catalog_rows(
    rows: list[dict[str, object]],
    review: dict[str, int],
    require_review: bool,
    require_complete_review: bool,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    kept: list[dict[str, object]] = []
    missing_samples: list[str] = []
    stats: dict[str, object] = {
        "catalog_links": len(rows),
        "review_accepted_links": 0,
        "review_rejected_links": 0,
        "review_missing_links": 0,
        "review_missing_samples": missing_samples,
    }
    for row in rows:
        decision = review_decision_for_row(row, review)
        if decision is None:
            if review or require_review or require_complete_review:
                stats["review_missing_links"] = int(stats["review_missing_links"]) + 1
                if len(missing_samples) < 10:
                    keys = review_keys_for_row(row)
                    missing_samples.append(keys[0] if keys else str(row.get("linkage_id", "")))
            if require_review or require_complete_review:
                continue
            kept.append(row)
            continue
        if decision == 1:
            stats["review_accepted_links"] = int(stats["review_accepted_links"]) + 1
            kept.append(row)
        else:
            stats["review_rejected_links"] = int(stats["review_rejected_links"]) + 1
    if require_complete_review and int(stats["review_missing_links"]) > 0:
        raise ValueError(
            f"Review CSV is incomplete: {stats['review_missing_links']} catalog links have no 0/1 decision; "
            f"sample keys: {missing_samples}"
        )
    return kept, stats


def write_filtered_catalog(path: str, rows: list[dict[str, object]]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_filtered_fits(path: str, source_fits: str, kept_rows: list[dict[str, object]]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    keep_keys = {
        key
        for row in kept_rows
        for key in (
            str(row.get("trk_sub", "")).strip(),
            str(row.get("trkSub", "")).strip(),
            str(row.get("linkage_id", "")).strip(),
        )
        if key
    }
    src = Path(source_fits)
    if src.exists():
        table = Table.read(src)
        key_col = next((col for col in ("trk_sub", "trkSub", "linkage_id") if col in table.colnames), "")
        if key_col:
            mask = np.asarray([str(value).strip() in keep_keys for value in table[key_col]], dtype=bool)
            table[mask].write(out, overwrite=True)
            return
    Table(rows=kept_rows).write(out, overwrite=True)
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
    out = {col: float(row[col]) for col in required}
    for optional in ("Flux_Aper4", "FluxErr_Aper4"):
        if optional in table.colnames:
            out[optional] = float(row[optional])
    return out


def build_obs_rows(
    args: argparse.Namespace,
    return_catalog: bool = False,
) -> tuple[list[dict[str, object]], dict[str, object]] | tuple[list[dict[str, object]], dict[str, object], list[dict[str, object]]]:
    catalog_path = Path(args.catalog)
    processed_root = Path(args.processed_root)
    require_complete_review = bool(getattr(args, "require_complete_review", False))
    review, review_meta = read_review_csv(args.review_csv, require_values=require_complete_review)
    require_review = bool(args.require_review)
    rows = load_unknown_rows(catalog_path)
    export_rows, review_stats = filter_catalog_rows(rows, review, require_review, require_complete_review)

    obs_rows: list[dict[str, object]] = []
    stats: dict[str, object] = {
        **review_stats,
        "review_csv": str(args.review_csv or ""),
        "review_csv_rows": review_meta["rows"],
        "review_csv_decisions": review_meta["decisions"],
        "review_csv_blank_values": review_meta["blank_values"],
        "missing_trk_sub_links": 0,
        "missing_detection_rows": 0,
        "exported_observations": 0,
    }

    for row in export_rows:
        trk_sub = validate_trk_sub(row.get("trk_sub") or row.get("trkSub"))
        if not trk_sub:
            stats["missing_trk_sub_links"] = int(stats["missing_trk_sub_links"]) + 1
            continue
        image_names = split_semicolon(row.get("image_names"))
        objids = [int(x) for x in split_semicolon(row.get("objids"))]
        fallback_mjds = [float(x) for x in split_semicolon(row.get("mjds"))]
        fallback_ras = [float(x) for x in split_semicolon(row.get("ras_deg"))]
        fallback_decs = [float(x) for x in split_semicolon(row.get("decs_deg"))]
        for idx, image_name in enumerate(image_names):
            if idx >= len(objids):
                stats["missing_detection_rows"] = int(stats["missing_detection_rows"]) + 1
                continue
            det = load_detection_row(processed_root, args.night, image_name, objids[idx])
            if det is None:
                stats["missing_detection_rows"] = int(stats["missing_detection_rows"]) + 1
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
                stats["missing_detection_rows"] = int(stats["missing_detection_rows"]) + 1
                continue
            obs_rows.append(
                {
                    "trkSub": trk_sub,
                    "mode": args.mode,
                    "stn": args.obs_code,
                    "prog": args.prog,
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
                    "logSNR": compute_logsnr(det.get("Flux_Aper4"), det.get("FluxErr_Aper4")),
                }
            )
    stats["exported_observations"] = len(obs_rows)
    stats["exported_links"] = len(export_rows)
    if return_catalog:
        return obs_rows, stats, export_rows
    return obs_rows, stats


def write_psv(args: argparse.Namespace, obs_rows: list[dict[str, object]]) -> None:
    cols = ["trkSub", "mode", "stn"]
    if args.prog:
        cols.append("prog")
    cols.extend(["obsTime", "ra", "dec", "rmsRA", "rmsDec", "astCat", "mag", "rmsMag", "band", "photCat"])
    if args.include_logsnr:
        cols.append("logSNR")
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
    ap.add_argument("--review-root", default="/pipeline/xiaoyunao/heliolincrr/review_packages")
    ap.add_argument("--catalog", default=None, help="Default: /processed1/<night>/L4/<night>_unknown_links.json")
    ap.add_argument("--unknown-fits", default=None, help="Default: /processed1/<night>/L4/<night>_unknown_links.fits")
    ap.add_argument("--out", default=None, help="Default: /processed1/<night>/L4/<night>_unknown_links_ades.psv")
    ap.add_argument("--review-csv", default="", help="Optional CSV with tracklet_id/trk_sub/linkage_id and is_real columns")
    ap.add_argument("--require-review", action="store_true", help="Only export links explicitly marked is_real=1")
    ap.add_argument(
        "--submit-csv",
        nargs="?",
        const="auto",
        default="",
        help="Use final web-submit CSV. With no value, defaults to <review-root>/<night>/<night>_submit.csv. Implies --require-review and --require-complete-review.",
    )
    ap.add_argument("--require-complete-review", action="store_true", help="Fail unless every catalog link has an explicit 0/1 decision")
    ap.add_argument("--filtered-catalog", default=None, help="Optional reviewed/masked unknown JSON output")
    ap.add_argument("--filtered-fits", default=None, help="Optional reviewed/masked unknown FITS output")
    ap.add_argument("--stats-out", default=None, help="Optional JSON stats output")

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
    ap.add_argument("--prog", default="", help="Optional MPC program code value")
    ap.add_argument("--include-logsnr", action="store_true", help="Include logSNR when flux columns exist")
    ap.add_argument("--ac2-email", default="wsgp2024@163.com")
    ap.add_argument("--ack", default="Unknown-object ADES submission")
    ap.add_argument("--obj-type", default="NEO")
    ap.add_argument("--response-out", default="")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--submit", action="store_true")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    l4_dir = Path(args.processed_root) / args.night / "L4"
    if args.catalog is None:
        args.catalog = str(l4_dir / f"{args.night}_unknown_links.json")
    if args.unknown_fits is None:
        args.unknown_fits = str(l4_dir / f"{args.night}_unknown_links.fits")
    using_submit_csv = bool(args.submit_csv)
    if using_submit_csv:
        if args.review_csv:
            raise ValueError("Use either --submit-csv or --review-csv, not both.")
        if args.submit_csv == "auto":
            args.submit_csv = str(Path(args.review_root) / args.night / f"{args.night}_submit.csv")
        args.review_csv = args.submit_csv
        args.require_review = True
        args.require_complete_review = True
    if args.out is None:
        suffix = "unknown_links_submit_ades.psv" if using_submit_csv else "unknown_links_ades.psv"
        args.out = str(l4_dir / f"{args.night}_{suffix}")
    if using_submit_csv:
        if args.filtered_catalog is None:
            args.filtered_catalog = str(l4_dir / f"{args.night}_unknown_links_submit_masked.json")
        if args.filtered_fits is None:
            args.filtered_fits = str(l4_dir / f"{args.night}_unknown_links_submit_masked.fits")
        if args.stats_out is None:
            args.stats_out = str(l4_dir / f"{args.night}_unknown_links_submit_stats.json")

    obs_rows, stats, kept_rows = build_obs_rows(args, return_catalog=True)
    if args.filtered_catalog:
        write_filtered_catalog(args.filtered_catalog, kept_rows)
        stats["filtered_catalog"] = args.filtered_catalog
    if args.filtered_fits:
        write_filtered_fits(args.filtered_fits, args.unknown_fits, kept_rows)
        stats["filtered_fits"] = args.filtered_fits
    write_psv(args, obs_rows)
    print(f"[OK] Wrote PSV: {args.out}")
    if args.stats_out:
        Path(args.stats_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.stats_out).write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"[OK] Wrote stats: {args.stats_out}")
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
