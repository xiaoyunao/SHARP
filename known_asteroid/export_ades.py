#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import tempfile
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


def row_quality(rms_ra: np.ndarray, rms_dec: np.ndarray, rms_mag: np.ndarray) -> np.ndarray:
    return (
        np.square(rms_ra, dtype=float)
        + np.square(rms_dec, dtype=float)
        + np.square(rms_mag, dtype=float)
    )


def dedup_current_rows(
    idx: np.ndarray,
    obj_key: np.ndarray,
    obs_time: np.ndarray,
    quality: np.ndarray,
) -> tuple[np.ndarray, int]:
    best_for_key: dict[tuple[str, str], tuple[float, int]] = {}
    duplicate_rows = 0
    for i in idx:
        key = (str(obj_key[i]), str(obs_time[i]))
        score = float(quality[i])
        prev = best_for_key.get(key)
        if prev is None or score < prev[0]:
            if prev is not None:
                duplicate_rows += 1
            best_for_key[key] = (score, int(i))
        else:
            duplicate_rows += 1
    dedup_idx = np.array(sorted(item[1] for item in best_for_key.values()), dtype=int)
    return dedup_idx, duplicate_rows


def load_submitted_history_keys(processed_root: str, current_night: str, current_out: str) -> set[tuple[str, str]]:
    if not processed_root:
        return set()
    root = Path(processed_root)
    current_out_path = Path(current_out).resolve() if current_out else None
    submitted: set[tuple[str, str]] = set()
    for path in sorted(root.glob("20*/L4/*_matched_asteroids_ades.psv")):
        if current_out_path is not None and path.resolve() == current_out_path:
            continue
        night = path.name.split("_", 1)[0]
        if current_night and night == current_night:
            continue
        try:
            with path.open("r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#") or line.startswith("!"):
                        continue
                    if "permID" in line and "obsTime" in line:
                        continue
                    parts = [p.strip() for p in line.rstrip("\n").split("|")]
                    if len(parts) < 6:
                        continue
                    obj = parts[0] or parts[1]
                    obs_time = parts[4]
                    if obj and obs_time:
                        submitted.add((obj, obs_time))
        except Exception:
            continue
    return submitted


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
        "Flag",
        "Flux_Aper4",
        "FluxErr_Aper4",
    ]
    missing = [col for col in required if col not in table.colnames]
    if missing:
        raise KeyError(f"Missing required columns in FITS: {missing}")

    perm_id = np.array([clean_permid(v) for v in table["number"]], dtype=object)
    prov_id = np.array([clean_provid(v) for v in table["name"]], dtype=object)
    obj_key = np.where(perm_id != "", perm_id, prov_id).astype(object)

    def numeric(name: str) -> np.ndarray:
        return np.asarray(np.ma.filled(np.ma.asarray(table[name], dtype=float), np.nan))

    obs_time = mjd_to_obs_time_z(numeric("epoch"))
    ra = numeric("RA_Win")
    dec = numeric("DEC_Win")
    rms_ra = numeric("RAErr_Win")
    rms_dec = numeric("DECErr_Win")
    if args.err_unit == "deg":
        rms_ra *= 3600.0
        rms_dec *= 3600.0
    elif args.err_unit == "mas":
        rms_ra /= 1000.0
        rms_dec /= 1000.0
    elif args.err_unit != "arcsec":
        raise ValueError("err-unit must be one of: arcsec, deg, mas")

    mag = numeric("Mag_Aper4")
    rms_mag = numeric("MagErr_Aper4")
    flux = numeric("Flux_Aper4")
    flux_err = numeric("FluxErr_Aper4")
    # Gate the actual serialized values as well as the input uncertainties.
    # Do not invent a precision floor for invalid/saturated detections.
    rms_ra_text = np.array([f"{v:.4f}" for v in rms_ra])
    rms_dec_text = np.array([f"{v:.4f}" for v in rms_dec])
    rms_mag_text = np.array([f"{v:.3f}" for v in rms_mag])

    has_id = obj_key != ""
    finite_core = (
        np.isfinite(ra)
        & np.isfinite(dec)
        & np.isfinite(rms_ra)
        & np.isfinite(rms_dec)
        & np.isfinite(mag)
        & np.isfinite(rms_mag)
        & (rms_ra > 0.0)
        & (rms_dec > 0.0)
        & (rms_mag > 0.0)
    )
    gates = {
        'missing_identity': has_id,
        'invalid_core_or_uncertainty': finite_core,
        'invalid_time_or_coordinates': (obs_time != '') & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90),
        'nonzero_or_missing_flag': numeric('Flag') == 0,
        'invalid_photometry': np.isfinite(flux) & (flux > 0) & np.isfinite(flux_err) & (flux_err > 0)
                              & (mag != 99.0) & (rms_mag < 99.0),
        'uncertainty_rounds_to_zero': (rms_ra_text.astype(float) > 0) & (rms_dec_text.astype(float) > 0)
                                     & (rms_mag_text.astype(float) > 0),
    }
    ok = np.logical_and.reduce(list(gates.values()))
    reasons = [';'.join(name for name, valid in gates.items() if not valid[i]) for i in range(len(table))]

    idx = np.where(ok)[0]
    quality = row_quality(rms_ra, rms_dec, rms_mag)
    before = set(idx)
    idx, duplicate_rows = dedup_current_rows(idx, obj_key, obs_time, quality)
    for i in before - set(idx):
        reasons[i] = 'duplicate_object_epoch'

    history_removed = 0
    if args.dedup_history_root:
        history_keys = load_submitted_history_keys(
            processed_root=args.dedup_history_root,
            current_night=args.current_night,
            current_out=args.out,
        )
        if history_keys:
            keep_mask = np.array(
                [(str(obj_key[i]), str(obs_time[i])) not in history_keys for i in idx],
                dtype=bool,
            )
            history_removed = int(np.count_nonzero(~keep_mask))
            for i in idx[~keep_mask]:
                reasons[i] = 'present_in_history'
            idx = idx[keep_mask]

    # Minimum arc length applies after quality cuts and both deduplication steps.
    counts: dict[str, int] = {}
    for key in obj_key[idx]:
        counts[key] = counts.get(key, 0) + 1
    keep = np.array([counts[obj_key[i]] >= int(args.min_observations) for i in idx], dtype=bool)
    dropped_singletons = len(set(obj_key[idx[~keep]]))
    for i in idx[~keep]:
        reasons[i] = 'below_minimum_observations_after_filtering'
    idx = idx[keep]

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
                rms_ra_text[i],
                rms_dec_text[i],
                args.astcat,
                f"{mag[i]:.3f}",
                rms_mag_text[i],
                args.band,
                args.photcat,
            ]
        )
        if include_logsnr:
            row.append(logsnr[i])
        lines.append("|".join(row))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = "\n".join(lines) + "\n"
    proof = validate_payload(payload) if len(idx) else {'status': 'no_observations', 'network_requests': 0}
    # Never replace an existing report with an invalid or partially written PSV.
    with tempfile.NamedTemporaryFile('w', encoding='utf-8', dir=out_path.parent, delete=False) as stream:
        temp = Path(stream.name)
        stream.write(payload)
    temp.chmod(out_path.stat().st_mode & 0o777 if out_path.exists() else 0o644)
    temp.replace(out_path)
    stats = {'input_fits': str(args.fits), 'input_rows': len(table), 'exported_rows': len(idx),
             'exported_objects': len(set(obj_key[idx])), 'dropped_rows': len(table) - len(idx),
             'quality_failures': {name: int(np.count_nonzero(~v)) for name, v in gates.items()},
             'duplicate_rows': duplicate_rows, 'history_removed': history_removed,
             'objects_below_minimum_after_filtering': dropped_singletons, 'validation': proof}
    Path(str(out_path) + '.audit.json').write_text(json.dumps(stats, indent=2) + '\n')
    with Path(str(out_path) + '.excluded.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['fits_row', 'object', 'obsTime', 'source_file', 'objID', 'reason'])
        for i, reason in enumerate(reasons):
            if reason:
                writer.writerow([i, obj_key[i], obs_time[i],
                                 safe_str(table['source_file'][i]) if 'source_file' in table.colnames else '',
                                 safe_str(table['objID'][i]) if 'objID' in table.colnames else '', reason])

    exported_objects = len(set(obj_key[idx]))
    print(f"[OK] Read FITS: {args.fits}")
    print(f"[OK] Wrote PSV: {args.out}")
    print(
        f"[INFO] Rows: input={len(table)}, exported={len(idx)}, dropped={len(table) - len(idx)}"
    )
    if duplicate_rows:
        print(f"[INFO] Dedup within current export: dropped_duplicate_rows={duplicate_rows}")
    if history_removed:
        print(f"[INFO] Dedup against submitted history: dropped_already_submitted_rows={history_removed}")
    print(
        f"[INFO] Objects: exported={exported_objects}, dropped_below_min_obs={dropped_singletons}, min_obs={args.min_observations}"
    )
    return int(len(idx))


def validate_payload(payload: str) -> dict:
    python = os.environ.get('MPC_ADES_VALIDATOR_PYTHON',
                            '/pipeline/xiaoyunao/data/heliolincrr/digest2_runtime/env/bin/python')
    proc = subprocess.run([python, str(Path(__file__).with_name('validate_ades.py')), '--stdin'],
                          input=payload, text=True, capture_output=True, timeout=120)
    if proc.returncode:
        raise ValueError('Local ADES validation failed; no submission allowed:\n' + proc.stderr)
    proof = json.loads(proc.stdout)
    if proof.get('status') != 'local_validated' or proof.get('psv_sha256') != hashlib.sha256(payload.encode()).hexdigest():
        raise ValueError('Local ADES validation did not confirm the exact PSV payload')
    return proof


def run_curl(url: str, args: argparse.Namespace) -> None:
    hold = Path(os.environ.get('MPC_SUBMISSION_HOLD_FILE',
                               '/pipeline/xiaoyunao/data/MPC_SUBMISSION_HOLD'))
    if url == SUBMIT_URL and hold.exists():
        raise RuntimeError(f'MPC submission paused: {hold}; no request sent')
    payload = Path(args.out).read_text()
    validate_payload(payload)
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
        "source=<-",
    ]
    print("[RUN]", " ".join(cmd))
    proc = subprocess.run(cmd, input=payload, text=True, capture_output=True)
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
    ap.add_argument("--validate", action="store_true", help="Local ADES schema validation only; no network request (also applied to every nonempty export)")
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--dedup-history-root", default="", help="Optional processed root used to filter rows already present in prior ADES files")
    ap.add_argument("--current-night", default="", help="Optional YYYYMMDD used to exclude the current night from history dedup")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    exported_rows = export_psv(args)
    if exported_rows == 0:
        print("[SKIP] No ADES obsData rows exported; skip validate/submit")
        return
    if args.validate:
        print('[LOCAL-VALIDATED] ADES schema and measurement quality; no network request')
    if args.submit:
        run_curl(SUBMIT_URL, args)


if __name__ == "__main__":
    main()
