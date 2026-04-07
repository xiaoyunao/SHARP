#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import os
import re
import traceback
import warnings

import astropy.units as u
import numpy as np
from aleph.Query import Query
from astropy.coordinates import Angle, EarthLocation, SkyCoord, match_coordinates_sky
from astropy.io import fits
from astropy.table import Table, hstack, vstack
from astropy.time import Time
from astropy.wcs import WCS


warnings.filterwarnings("ignore")


def query_asteroids(
    q: Query,
    field_center: SkyCoord,
    field_radius: u.Quantity,
    epoch: Time,
    observer: EarthLocation,
    mag_limit: float,
    njobs: int,
    confidence_radius: u.Quantity,
) -> Table | None:
    ephs = q.query_mixed_cat(
        field_center=field_center,
        radius=field_radius,
        epoch=epoch,
        observer=observer,
        njobs=njobs,
        confidence_radius=confidence_radius,
    )
    if ephs is None or len(ephs) == 0:
        return None

    ephs = ephs[ephs["V"] < mag_limit]
    if len(ephs) == 0:
        return None

    out = Table()
    out["name"] = np.array(ephs["name"], dtype="U32")
    out["number"] = np.array(ephs["number"], dtype="f8")
    out["ra"] = np.array(ephs["ra"], dtype="f8")
    out["dec"] = np.array(ephs["dec"], dtype="f8")
    out["mag"] = np.array(ephs["V"], dtype="f4")
    return out


def log(msg: str, log_path: str | None) -> None:
    print(msg, flush=True)
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(msg.rstrip() + "\n")


def parse_epoch(
    header,
    obs_date_key: str,
    exptime_key: str | None,
    exptime_default: float,
) -> tuple[Time | None, float | None]:
    t0 = None
    for key in (obs_date_key, "DATE-OBS", "DATEOBS"):
        if not key or key not in header:
            continue
        try:
            t0 = Time(header[key], format="isot", scale="utc")
            break
        except Exception:
            continue
    if t0 is None:
        return None, None

    exptime = None
    if exptime_key:
        try:
            exptime = float(header[exptime_key])
        except Exception:
            exptime = None
    if exptime is None:
        exptime = float(exptime_default)

    epoch = t0 + (exptime / 86400.0) / 2.0
    return epoch, float(epoch.mjd)


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Match known asteroids for one night.")
    p.add_argument("night", help="Night folder name, e.g. 20250304")

    p.add_argument("--root", default="/processed1", help="Root directory containing <NIGHT>/L2")
    p.add_argument("--outdir", default="./products", help="Output directory")
    p.add_argument("--astorb", default="../astorb.dat", help="Path to astorb.dat")
    p.add_argument("--file", default=None, help="Process only one catalog filename in <NIGHT>/L2")
    p.add_argument(
        "--parts-dir",
        default=None,
        help="When --file is used, write per-file outputs into this directory instead of night-level outputs",
    )
    p.add_argument(
        "--file-regex",
        default=r".*MP.*\.(fits|fits\.gz)$",
        help=r"Regex used to select catalog files",
    )
    p.add_argument("--hdu", type=int, default=1, help="HDU containing catalog table")

    p.add_argument("--sep-arcsec", type=float, default=1.0, help="Sky match radius in arcsec")
    p.add_argument("--mag-limit", type=float, default=22.5, help="Predicted V magnitude limit")
    p.add_argument(
        "--magdiff",
        type=float,
        default=99.0,
        help="Require |predicted mag - measured mag| < magdiff",
    )
    p.add_argument("--cat-mag-col", default="Mag_Kron", help="Measured magnitude column")
    p.add_argument("--njobs", type=int, default=1, help="Aleph query parallel jobs (keep at 1 for file-level parallelism)")
    p.add_argument("--confidence-pad-deg", type=float, default=0.5)
    p.add_argument("--corner-pad-deg", type=float, default=0.05)

    p.add_argument("--nx", type=int, default=9216, help="Image width in pixels")
    p.add_argument("--ny", type=int, default=9232, help="Image height in pixels")
    p.add_argument("--obs-date-key", default="OBS_DATE", help="Observation time keyword")
    p.add_argument("--exptime-key", default=None, help="Exposure time keyword")
    p.add_argument("--exptime-sec", type=float, default=30.0, help="Fallback exposure time")

    p.add_argument("--lon", type=float, default=117.575, help="Observatory longitude")
    p.add_argument("--lat", type=float, default=40.393, help="Observatory latitude")
    p.add_argument("--height", type=float, default=960.0, help="Observatory height in m")

    p.add_argument("--log", default=None, help="Optional append log path")
    p.add_argument("--verbose", action="store_true", help="Print per-file progress")
    return p


def file_tag(filename: str) -> str:
    if filename.endswith(".fits.gz"):
        return filename[:-8]
    if filename.endswith(".fits"):
        return filename[:-5]
    return filename


def run_match(args: argparse.Namespace) -> tuple[str | None, str | None]:
    outdir_real = os.path.realpath(args.outdir)
    known_real = os.path.realpath(os.path.dirname(args.astorb))
    try:
        if os.path.commonpath([outdir_real, known_real]) == known_real:
            raise SystemExit(
                f"[FATAL] outdir points inside the astorb/config directory: {outdir_real}"
            )
    except Exception:
        pass

    night = args.night
    if not re.fullmatch(r"\d{8}", night):
        raise SystemExit(f"[FATAL] night must be YYYYMMDD, got: {night}")

    l2_dir = os.path.join(args.root, night, "L2")
    if not os.path.isdir(l2_dir):
        raise SystemExit(f"[FATAL] L2 dir not found: {l2_dir}")

    os.makedirs(args.outdir, exist_ok=True)

    q = Query(service="Lowell", filename=args.astorb)
    observer = EarthLocation(
        lon=args.lon * u.deg,
        lat=args.lat * u.deg,
        height=args.height * u.m,
    )

    file_re = re.compile(args.file_regex)
    files = sorted([f for f in os.listdir(l2_dir) if file_re.fullmatch(f)])
    if args.file:
        if args.file not in files:
            raise SystemExit(f"[FATAL] requested file not found or not matched by regex: {args.file}")
        files = [args.file]
    if not files:
        log(f"[WARN] No files matched in {l2_dir} with regex {args.file_regex}", args.log)
        return None, None

    sep_limit = args.sep_arcsec * u.arcsec
    all_asteroids: list[Table] = []
    all_matches: list[Table] = []

    log(f"[START] night={night} n_files={len(files)} l2={l2_dir}", args.log)

    for i, fname in enumerate(files, 1):
        fpath = os.path.join(l2_dir, fname)
        if args.verbose:
            log(f"[INFO] ({i}/{len(files)}) {fname}", args.log)

        try:
            with fits.open(fpath, memmap=False) as hdul:
                header = hdul[args.hdu].header
                cat = Table(hdul[args.hdu].data)
                w = WCS(header)

                nx = int(header.get("NAXIS1", args.nx))
                ny = int(header.get("NAXIS2", args.ny))
                pix_corners = np.array([[1, 1], [1, ny], [nx, 1], [nx, ny]], dtype=float)
                world = w.all_pix2world(pix_corners, 1)
                ra_vals, dec_vals = world[:, 0], world[:, 1]

                corner_coords = SkyCoord(ra=ra_vals * u.deg, dec=dec_vals * u.deg)
                center_world = w.pixel_to_world((nx + 1) / 2.0, (ny + 1) / 2.0)
                field_center = SkyCoord(ra=center_world.ra, dec=center_world.dec)
                max_diff = field_center.separation(corner_coords).max() + args.corner_pad_deg * u.deg
                confidence_radius = (max_diff.to(u.deg).value + args.confidence_pad_deg) * u.deg

                epoch, mjd = parse_epoch(
                    header,
                    args.obs_date_key,
                    args.exptime_key,
                    args.exptime_sec,
                )
                if epoch is None:
                    log(f"[WARN] failed to parse {args.obs_date_key} for {fname}; skip", args.log)
                    continue
        except Exception as exc:
            log(f"[ERROR] failed to read {fname}: {exc}", args.log)
            if args.verbose:
                traceback.print_exc()
            continue

        try:
            ephs = query_asteroids(
                q=q,
                field_center=field_center,
                field_radius=max_diff,
                epoch=epoch,
                observer=observer,
                mag_limit=args.mag_limit,
                njobs=args.njobs,
                confidence_radius=confidence_radius,
            )
            if ephs is None or len(ephs) == 0:
                continue

            x_pix, y_pix = w.all_world2pix(ephs["ra"], ephs["dec"], 1)
            inside = (x_pix >= 1) & (x_pix <= nx) & (y_pix >= 1) & (y_pix <= ny)
            if not np.any(inside):
                continue

            ephs = ephs[inside]
            ephs["source_file"] = np.array([fname] * len(ephs), dtype="U128")
            ephs["epoch"] = np.array([mjd] * len(ephs), dtype="f8")
            all_asteroids.append(ephs)

            asteroid_coords = SkyCoord(ephs["ra"] * u.deg, ephs["dec"] * u.deg)
            cat_coords = SkyCoord(cat["RA_Win"] * u.deg, cat["DEC_Win"] * u.deg)
            idx, sep2d, _ = match_coordinates_sky(asteroid_coords, cat_coords)

            matched = sep2d < sep_limit
            if args.cat_mag_col not in cat.colnames:
                log(
                    f"[WARN] catalog missing {args.cat_mag_col} in {fname}; skip matches for this frame",
                    args.log,
                )
                matched &= False
            else:
                det_mag = np.array(cat[args.cat_mag_col][idx], dtype="f4")
                pred_mag = np.array(ephs["mag"], dtype="f4")
                matched &= (
                    np.isfinite(det_mag)
                    & np.isfinite(pred_mag)
                    & (np.abs(pred_mag - det_mag) < float(args.magdiff))
                )

            if np.any(matched):
                combined = hstack([ephs[matched], cat[idx[matched]]], join_type="exact")
                all_matches.append(combined)
        except Exception as exc:
            log(f"[ERROR] query/match failed on {fname}: {exc}", args.log)
            if args.verbose:
                traceback.print_exc()
            continue

    if args.file and args.parts_dir:
        os.makedirs(args.parts_dir, exist_ok=True)
        stem = file_tag(args.file)
        out_all = os.path.join(args.parts_dir, f"{stem}_all_asteroids.fits")
        out_matched = os.path.join(args.parts_dir, f"{stem}_matched_asteroids.fits")
    else:
        out_all = os.path.join(args.outdir, f"{night}_all_asteroids.fits")
        out_matched = os.path.join(args.outdir, f"{night}_matched_asteroids.fits")

    if all_asteroids:
        total_ast = vstack(all_asteroids)
        total_ast.write(out_all, overwrite=True)
        log(f"[WRITE] {out_all} (n={len(total_ast)})", args.log)
    else:
        out_all = None
        log("[WRITE] skipped all_asteroids (no results)", args.log)

    if all_matches:
        total_match = vstack(all_matches)
        total_match.write(out_matched, overwrite=True)
        log(f"[WRITE] {out_matched} (n={len(total_match)})", args.log)
    else:
        out_matched = None
        log("[WRITE] skipped matched_asteroids (no matches)", args.log)

    log(
        f"[DONE] night={night} ast_frames={len(all_asteroids)} match_frames={len(all_matches)}",
        args.log,
    )
    return out_all, out_matched


def main() -> None:
    args = build_argparser().parse_args()
    run_match(args)


if __name__ == "__main__":
    main()
