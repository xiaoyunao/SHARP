#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path
import re

import numpy as np
from astropy.table import Table


SUBMIT_URL = "https://minorplanetcenter.net/submit_psv"
VALIDATE_URL = "https://www.minorplanetcenter.net/submit_psv_test"
PROVID_RE = re.compile(r"^\d{4}\s+[A-Z]{1,2}\d{0,3}[A-Z]?$")


def safe_str(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (bytes, bytearray)):
        value = value.decode(errors="ignore")
    return str(value).strip()


def clean_permid(value) -> str:
    s = safe_str(value)
    if s == "" or s.lower() in {"nan", "none", "null", "--", "-", "—"}:
        return ""
    try:
        return str(int(float(s)))
    except Exception:
        return ""


def clean_provid(value) -> str:
    s = safe_str(value)
    if s == "" or s.lower() in {"nan", "none", "null", "--", "-", "—"}:
        return ""
    return s if PROVID_RE.match(s) else ""


def split_semicolon_list(value: str) -> list[str]:
    if not value:
        return []
    out: list[str] = []
    for item in value.split(";"):
        item = item.strip()
        if item:
            out.append(item)
    return out


def mjd_to_obs_time_z(mjd: np.ndarray) -> np.ndarray:
    mjd = np.asarray(mjd, dtype=float)
    out = np.empty(mjd.shape, dtype=object)
    mjd0 = datetime(1858, 11, 17, tzinfo=timezone.utc)
    for idx, value in np.ndenumerate(mjd):
        if not np.isfinite(value):
            out[idx] = ""
            continue
        dt = mjd0 + timedelta(days=float(value))
        out[idx] = dt.isoformat(timespec="milliseconds").replace("+00:00", "Z")
    return out


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

    obs_list = split_semicolon_list(args.observers)
    if obs_list:
        lines.append("# observers")
        for name in obs_list:
            lines.append(f"! name {name}")

    meas_list = split_semicolon_list(args.measurers)
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

    coinvestigator_list = split_semicolon_list(args.coinvestigators)
    if coinvestigator_list:
        lines.append("# coinvestigators")
        for name in coinvestigator_list:
            lines.append(f"! name {name}")

    collaborator_list = split_semicolon_list(args.collaborators)
    if collaborator_list:
        lines.append("# collaborators")
        for name in collaborator_list:
            lines.append(f"! name {name}")

    if args.funding_source:
        lines.append(f"# fundingSource {args.funding_source}")

    comment_list = split_semicolon_list(args.comment)
    if comment_list:
        lines.append("# comment")
        for line in comment_list:
            lines.append(f"! line {line}")

    lines.append("")
    return "\n".join(lines)


def compute_logsnr(table: Table) -> np.ndarray | None:
    if "Flux_Aper4" not in table.colnames or "FluxErr_Aper4" not in table.colnames:
        return None
    flux = np.array(table["Flux_Aper4"], dtype=float)
    flux_err = np.array(table["FluxErr_Aper4"], dtype=float)
    snr = np.divide(
        flux,
        flux_err,
        out=np.full_like(flux, np.nan, dtype=float),
        where=(flux_err != 0),
    )
    out = np.full(len(table), "", dtype=object)
    valid = np.isfinite(snr) & (snr > 0)
    out[valid] = np.array([f"{v:.2f}" for v in np.log10(snr[valid])], dtype=object)
    return out


def export_psv(args: argparse.Namespace) -> int:
    table = Table.read(args.fits)

    required = [
        "number",
        "name",
        "epoch",
        "RA_Win",
        "DEC_Win",
        "RAErr_Win",
        "DECErr_Win",
        "Mag_Aper4",
        "MagErr_Aper4",
    ]
    missing = [col for col in required if col not in table.colnames]
    if missing:
        raise KeyError(f"Missing required columns in FITS: {missing}")

    perm_id = np.array([clean_permid(v) for v in table["number"]], dtype=object)
    prov_id = np.array([clean_provid(v) for v in table["name"]], dtype=object)
    obj_key = np.where(perm_id != "", perm_id, prov_id).astype(object)

    obs_time = mjd_to_obs_time_z(np.array(table["epoch"], dtype=float))
    ra = np.array(table["RA_Win"], dtype=float)
    dec = np.array(table["DEC_Win"], dtype=float)
    rms_ra = np.array(table["RAErr_Win"], dtype=float)
    rms_dec = np.array(table["DECErr_Win"], dtype=float)
    if args.err_unit == "deg":
        rms_ra *= 3600.0
        rms_dec *= 3600.0
    elif args.err_unit == "mas":
        rms_ra /= 1000.0
        rms_dec /= 1000.0
    elif args.err_unit != "arcsec":
        raise ValueError("err-unit must be one of: arcsec, deg, mas")

    mag = np.array(table["Mag_Aper4"], dtype=float)
    rms_mag = np.array(table["MagErr_Aper4"], dtype=float)

    has_id = obj_key != ""
    finite_core = (
        np.isfinite(ra)
        & np.isfinite(dec)
        & np.isfinite(rms_ra)
        & np.isfinite(rms_dec)
        & np.isfinite(mag)
        & np.isfinite(rms_mag)
        & (rms_mag > 0.0)
    )
    ok = has_id & finite_core

    counts = {}
    for key in obj_key[ok]:
        counts[key] = counts.get(key, 0) + 1
    per_object_count = np.array([counts.get(key, 0) for key in obj_key], dtype=int)
    ok &= per_object_count >= int(args.min_observations)

    idx = np.where(ok)[0]
    header = build_header(args)

    include_prog = bool(args.prog)
    include_logsnr = bool(args.include_logsnr)
    logsnr = compute_logsnr(table) if include_logsnr else None
    if include_logsnr and logsnr is None:
        include_logsnr = False

    cols = [
        "permID",
        "provID",
        "mode",
        "stn",
    ]
    if include_prog:
        cols.append("prog")
    cols.extend(
        [
            "obsTime",
            "ra",
            "dec",
            "rmsRA",
            "rmsDec",
            "astCat",
            "mag",
            "rmsMag",
            "band",
            "photCat",
        ]
    )
    if include_logsnr:
        cols.append("logSNR")

    lines = [header.rstrip("\n"), " | ".join(cols)]
    for i in idx:
        row = [
            perm_id[i],
            prov_id[i],
            args.mode,
            args.obs_code,
        ]
        if include_prog:
            row.append(args.prog)
        row.extend(
            [
                obs_time[i],
                f"{ra[i]:.10f}",
                f"{dec[i]:.10f}",
                f"{rms_ra[i]:.4f}",
                f"{rms_dec[i]:.4f}",
                args.astcat,
                f"{mag[i]:.3f}",
                f"{rms_mag[i]:.3f}",
                args.band,
                args.photcat,
            ]
        )
        if include_logsnr:
            row.append(logsnr[i])
        lines.append("|".join(row))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    exported_objects = len(set(obj_key[idx]))
    dropped_singletons = len(set(obj_key[has_id & (per_object_count < int(args.min_observations))]))
    print(f"[OK] Read FITS: {args.fits}")
    print(f"[OK] Wrote PSV: {args.out}")
    print(
        f"[INFO] Rows: input={len(table)}, exported={len(idx)}, dropped={len(table) - len(idx)}"
    )
    print(
        f"[INFO] Objects: exported={exported_objects}, dropped_below_min_obs={dropped_singletons}, min_obs={args.min_observations}"
    )
    return int(len(idx))


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
    ap = argparse.ArgumentParser(
        description="Export a conservative ADES PSV file for known-object submissions."
    )
    ap.add_argument("--fits", required=True, help="Input matched FITS file")
    ap.add_argument("--out", required=True, help="Output PSV file")

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
    ap.add_argument(
        "--collaborators",
        default="Yun-Ao Xiao;Niu Li;Jingyi Zhang;Wenxiong Li;Zhaobin Chen;Shufei Liu;WSGP Team",
    )
    ap.add_argument(
        "--funding-source",
        default="The 60/90-cm Schmidt telescope situated at the Xinglong station is overseen by the research team of the Wide-field Survey and Galaxy Physics (WSGP) at NAOC.",
    )

    ap.add_argument("--astcat", default="Gaia3E")
    ap.add_argument("--band", default="G")
    ap.add_argument("--photcat", default="Gaia3E")
    ap.add_argument("--err-unit", default="deg", choices=["arcsec", "deg", "mas"])
    ap.add_argument("--min-observations", type=int, default=3)
    ap.add_argument("--prog", default="", help="Optional MPC program code value")
    ap.add_argument("--include-logsnr", action="store_true", help="Include logSNR when flux columns exist")

    ap.add_argument("--ac2-email", default="wsgp2024@163.com")
    ap.add_argument("--ack", default="Known-object ADES submission")
    ap.add_argument("--obj-type", default="MBA")
    ap.add_argument("--response-out", default="")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--submit", action="store_true")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    exported_rows = export_psv(args)
    if exported_rows == 0:
        print("[SKIP] No ADES obsData rows exported; skip validate/submit")
        return
    if args.validate:
        run_curl(VALIDATE_URL, args)
    if args.submit:
        run_curl(SUBMIT_URL, args)


if __name__ == "__main__":
    main()
