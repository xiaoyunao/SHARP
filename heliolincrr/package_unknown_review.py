#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import tarfile
from argparse import Namespace
from pathlib import Path
from typing import Any

import numpy as np
from astropy.table import Table

from export_unknown_ades import build_obs_rows

DETAIL_COLUMNS = [
    ("trk_sub", "U7"),
    ("linkage_id", np.int64),
    ("detection_index", np.int64),
    ("image_name", "U128"),
    ("catalog_name", "U128"),
    ("objID", np.int64),
    ("source_tracklet_ids", "U2048"),
    ("link_tracklet_ids", "U2048"),
    ("groups", "U512"),
    ("exp_pairs", "U512"),
    ("mjd", np.float64),
    ("RA_Win", np.float64),
    ("DEC_Win", np.float64),
    ("RA_PSF", np.float64),
    ("DEC_PSF", np.float64),
    ("RAErr_Win", np.float64),
    ("DECErr_Win", np.float64),
    ("RAErr_PSF", np.float64),
    ("DECErr_PSF", np.float64),
    ("X_Win", np.float64),
    ("Y_Win", np.float64),
    ("Mag_Aper4", np.float64),
    ("MagErr_Aper4", np.float64),
    ("Mag_PSF", np.float64),
    ("MagErr_PSF", np.float64),
    ("Flag", np.int64),
    ("n_obs", np.int64),
    ("n_tracklets", np.int64),
    ("rms_arcsec", np.float64),
    ("med_arcsec", np.float64),
    ("max_arcsec", np.float64),
    ("a_au", np.float64),
    ("ecc", np.float64),
    ("inc_deg", np.float64),
    ("raan_deg", np.float64),
    ("argp_deg", np.float64),
    ("nu_deg", np.float64),
    ("best_v1_kms", np.float64),
    ("lin_rms_arcsec", np.float64),
    ("lin_speed_arcsec_per_day", np.float64),
    ("lin_dir_deg", np.float64),
    ("fallback_ra_deg", np.float64),
    ("fallback_dec_deg", np.float64),
]


def split_semicolon(value: Any) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def load_json_list(path: Path) -> list[dict[str, Any]]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError(f"Expected JSON list: {path}")
    return rows


def gif_name_for_link(night: str, linkage_id: int) -> str:
    return f"unknown_link_{linkage_id:04d}_{night}.gif"


def safe_trk_sub(row: dict[str, Any]) -> str:
    trk_sub = str(row.get("trk_sub") or row.get("trkSub") or "").strip()
    if trk_sub:
        return trk_sub
    linkage_id = str(row.get("linkage_id", "")).strip()
    return f"link{int(linkage_id):04d}" if linkage_id else "unknown"


def finite_float(value: Any) -> float:
    try:
        if value is np.ma.masked:
            return np.nan
        out = float(value)
    except Exception:
        return np.nan
    return out if np.isfinite(out) else np.nan


def int_value(value: Any, default: int = -1) -> int:
    try:
        if value is np.ma.masked:
            return default
        return int(value)
    except Exception:
        return default


def catalog_name_for_image(image_name: str) -> str:
    return Path(image_name).name.replace(".fits.gz", "_cat.fits.gz")


def image_name_for_catalog(value: Any) -> str:
    name = Path(str(value)).name
    return name.replace("_cat.fits.gz", ".fits.gz")


def table_from_detail_rows(rows: list[dict[str, Any]]) -> Table:
    names = [name for name, _ in DETAIL_COLUMNS]
    if rows:
        return Table(rows=[{name: row.get(name) for name in names} for row in rows], names=names)
    return Table({name: np.asarray([], dtype=dtype) for name, dtype in DETAIL_COLUMNS})


def load_tracklet_detection_index(night: str, root_out: Path) -> dict[tuple[str, int], dict[str, list[str]]]:
    path = root_out / night / "tracklets_linreproj" / f"tracklets_{night}_ALL.fits"
    if not path.exists():
        return {}
    table = Table.read(path)
    out: dict[tuple[str, int], dict[str, list[str]]] = {}
    for row in table:
        tid = str(row["tracklet_id"]) if "tracklet_id" in table.colnames else ""
        group = str(int_value(row["group"])) if "group" in table.colnames else ""
        exp_pair = ""
        if "exp_i" in table.colnames and "exp_j" in table.colnames:
            exp_pair = f"{int_value(row['exp_i'])}->{int_value(row['exp_j'])}"
        for endpoint in (1, 2):
            file_col = f"file{endpoint}"
            obj_col = f"objID{endpoint}"
            if file_col not in table.colnames or obj_col not in table.colnames:
                continue
            key = (image_name_for_catalog(row[file_col]), int_value(row[obj_col]))
            item = out.setdefault(key, {"tracklet_ids": [], "groups": [], "exp_pairs": []})
            if tid and tid not in item["tracklet_ids"]:
                item["tracklet_ids"].append(tid)
            if group and group not in item["groups"]:
                item["groups"].append(group)
            if exp_pair and exp_pair not in item["exp_pairs"]:
                item["exp_pairs"].append(exp_pair)
    return out


def load_l2_detection(processed_root: Path, night: str, image_name: str, obj_id: int) -> dict[str, Any]:
    path = processed_root / night / "L2" / catalog_name_for_image(image_name)
    if not path.exists():
        return {}
    table = Table.read(path)
    if "objID" not in table.colnames:
        return {}
    matches = table[np.asarray(table["objID"], dtype=np.int64) == int(obj_id)]
    if len(matches) == 0:
        return {}
    row = matches[0]
    return {col: row[col] for col in table.colnames}


def write_review_detail_fits(
    args: argparse.Namespace,
    rows: list[dict[str, Any]],
    out_dir: Path,
) -> tuple[str, int]:
    night = str(args.night)
    processed_root = Path(args.processed_root)
    root_out = Path(args.root_out)
    tracklet_index = load_tracklet_detection_index(night, root_out)
    fits_path = out_dir / f"{night}_unknown_review_full.fits"
    detail_rows: list[dict[str, Any]] = []

    for row in rows:
        trk_sub = safe_trk_sub(row)
        linkage_id = int_value(row.get("linkage_id"))
        image_names = split_semicolon(row.get("image_names"))
        objids = [int_value(x) for x in split_semicolon(row.get("objids"))]
        fallback_ras = [finite_float(x) for x in split_semicolon(row.get("ras_deg"))]
        fallback_decs = [finite_float(x) for x in split_semicolon(row.get("decs_deg"))]
        link_tracklet_ids = split_semicolon(row.get("tracklet_ids"))

        for idx, image_name in enumerate(image_names):
            obj_id = objids[idx] if idx < len(objids) else -1
            l2 = load_l2_detection(processed_root, night, image_name, obj_id)
            tracklet_item = tracklet_index.get((Path(image_name).name, obj_id), {})
            detail_rows.append(
                {
                    "trk_sub": trk_sub,
                    "linkage_id": linkage_id,
                    "detection_index": idx,
                    "image_name": Path(image_name).name,
                    "catalog_name": catalog_name_for_image(image_name),
                    "objID": obj_id,
                    "source_tracklet_ids": ";".join(tracklet_item.get("tracklet_ids", [])),
                    "link_tracklet_ids": ";".join(link_tracklet_ids),
                    "groups": ";".join(tracklet_item.get("groups", [])),
                    "exp_pairs": ";".join(tracklet_item.get("exp_pairs", [])),
                    "mjd": finite_float(l2.get("MJD", np.nan)),
                    "RA_Win": finite_float(l2.get("RA_Win", np.nan)),
                    "DEC_Win": finite_float(l2.get("DEC_Win", np.nan)),
                    "RA_PSF": finite_float(l2.get("RA_PSF", np.nan)),
                    "DEC_PSF": finite_float(l2.get("DEC_PSF", np.nan)),
                    "RAErr_Win": finite_float(l2.get("RAErr_Win", np.nan)),
                    "DECErr_Win": finite_float(l2.get("DECErr_Win", np.nan)),
                    "RAErr_PSF": finite_float(l2.get("RAErr_PSF", np.nan)),
                    "DECErr_PSF": finite_float(l2.get("DECErr_PSF", np.nan)),
                    "X_Win": finite_float(l2.get("X_Win", l2.get("XWIN_IMAGE", np.nan))),
                    "Y_Win": finite_float(l2.get("Y_Win", l2.get("YWIN_IMAGE", np.nan))),
                    "Mag_Aper4": finite_float(l2.get("Mag_Aper4", np.nan)),
                    "MagErr_Aper4": finite_float(l2.get("MagErr_Aper4", np.nan)),
                    "Mag_PSF": finite_float(l2.get("Mag_PSF", np.nan)),
                    "MagErr_PSF": finite_float(l2.get("MagErr_PSF", np.nan)),
                    "Flag": int_value(l2.get("Flag", -1)),
                    "n_obs": int_value(row.get("n_obs")),
                    "n_tracklets": int_value(row.get("n_tracklets")),
                    "rms_arcsec": finite_float(row.get("rms_arcsec")),
                    "med_arcsec": finite_float(row.get("med_arcsec")),
                    "max_arcsec": finite_float(row.get("max_arcsec")),
                    "a_au": finite_float(row.get("a_au")),
                    "ecc": finite_float(row.get("ecc")),
                    "inc_deg": finite_float(row.get("inc_deg")),
                    "raan_deg": finite_float(row.get("raan_deg")),
                    "argp_deg": finite_float(row.get("argp_deg")),
                    "nu_deg": finite_float(row.get("nu_deg")),
                    "best_v1_kms": finite_float(row.get("best_v1_kms")),
                    "lin_rms_arcsec": finite_float(row.get("lin_rms_arcsec")),
                    "lin_speed_arcsec_per_day": finite_float(row.get("lin_speed_arcsec_per_day")),
                    "lin_dir_deg": finite_float(row.get("lin_dir_deg")),
                    "fallback_ra_deg": fallback_ras[idx] if idx < len(fallback_ras) else np.nan,
                    "fallback_dec_deg": fallback_decs[idx] if idx < len(fallback_decs) else np.nan,
                }
            )

    table_from_detail_rows(detail_rows).write(fits_path, overwrite=True)
    return str(fits_path), len(detail_rows)


def write_review_ades_fits(args: argparse.Namespace, catalog: Path, out_dir: Path) -> tuple[str, int]:
    fits_path = out_dir / f"{args.night}_unknown_review_ades.fits"
    export_args = Namespace(
        night=str(args.night),
        processed_root=str(args.processed_root),
        catalog=str(catalog),
        review_csv="",
        require_review=False,
        mode="CCD",
        obs_code="327",
        prog="",
        astcat="Gaia3E",
        band="G",
        photcat="Gaia3E",
        err_unit="deg",
    )
    obs_rows, _ = build_obs_rows(export_args)
    cols = ["trkSub", "mode", "stn", "obsTime", "ra", "dec", "rmsRA", "rmsDec", "astCat", "mag", "rmsMag", "band", "photCat"]
    if obs_rows:
        Table(rows=[{col: row.get(col, "") for col in cols} for row in obs_rows]).write(fits_path, overwrite=True)
    else:
        Table({col: np.asarray([], dtype="U64") for col in cols}).write(fits_path, overwrite=True)
    return str(fits_path), len(obs_rows)


def build_package(args: argparse.Namespace) -> dict[str, Any]:
    night = str(args.night)
    processed_root = Path(args.processed_root)
    plot_root = Path(args.plot_root)
    out_root = Path(args.out_root)
    catalog = Path(args.catalog) if args.catalog else processed_root / night / "L4" / f"{night}_unknown_links.json"
    plot_dir = plot_root / night
    out_dir = Path(args.out_dir) if args.out_dir else out_root / night
    gif_out_dir = out_dir / "gifs"
    out_dir.mkdir(parents=True, exist_ok=True)
    gif_out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_json_list(catalog)
    review_csv = out_dir / f"{night}_unknown_review.csv"
    manifest_json = out_dir / f"{night}_unknown_review_manifest.json"
    review_ades_fits, review_ades_rows = write_review_ades_fits(args, catalog, out_dir)
    review_full_fits, review_full_rows = write_review_detail_fits(args, rows, out_dir)

    copied = 0
    missing = 0
    review_rows: list[dict[str, Any]] = []
    manifest_rows: list[dict[str, Any]] = []

    for row in rows:
        linkage_id = int(row["linkage_id"])
        trk_sub = safe_trk_sub(row)
        src_gif = plot_dir / gif_name_for_link(night, linkage_id)
        dst_gif_name = f"{trk_sub}_link{linkage_id:04d}_{night}.gif"
        dst_gif = gif_out_dir / dst_gif_name
        if src_gif.exists():
            shutil.copy2(src_gif, dst_gif)
            copied += 1
            gif_status = "copied"
        else:
            missing += 1
            gif_status = "missing"

        review_rows.append({"tracklet_id": trk_sub, "is_real": ""})
        manifest_rows.append(
            {
                "tracklet_id": trk_sub,
                "linkage_id": linkage_id,
                "gif": str(dst_gif.relative_to(out_dir)) if src_gif.exists() else "",
                "source_gif": str(src_gif),
                "gif_status": gif_status,
                "n_obs": row.get("n_obs"),
                "n_tracklets": row.get("n_tracklets"),
                "rms_arcsec": row.get("rms_arcsec"),
                "a_au": row.get("a_au"),
                "ecc": row.get("ecc"),
                "inc_deg": row.get("inc_deg"),
                "tracklet_ids": split_semicolon(row.get("tracklet_ids")),
                "image_names": split_semicolon(row.get("image_names")),
                "objids": split_semicolon(row.get("objids")),
                "mjds": split_semicolon(row.get("mjds")),
                "ras_deg": split_semicolon(row.get("ras_deg")),
                "decs_deg": split_semicolon(row.get("decs_deg")),
            }
        )

    with review_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["tracklet_id", "is_real"])
        writer.writeheader()
        writer.writerows(review_rows)

    manifest = {
        "night": night,
        "catalog": str(catalog),
        "plot_dir": str(plot_dir),
        "out_dir": str(out_dir),
        "review_csv": str(review_csv),
        "review_ades_fits": review_ades_fits,
        "review_ades_rows": review_ades_rows,
        "review_full_fits": review_full_fits,
        "review_full_rows": review_full_rows,
        "n_catalog_rows": len(rows),
        "n_gifs_copied": copied,
        "n_gifs_missing": missing,
        "rows": manifest_rows,
    }
    manifest_json.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    tar_path = ""
    if args.make_tar:
        tar_file = out_root / f"{night}_unknown_review.tar.gz"
        tar_file.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_file, "w:gz") as tf:
            tf.add(out_dir, arcname=out_dir.name)
        tar_path = str(tar_file)

    return {
        "night": night,
        "out_dir": str(out_dir),
        "review_csv": str(review_csv),
        "manifest_json": str(manifest_json),
        "review_ades_fits": review_ades_fits,
        "review_ades_rows": review_ades_rows,
        "review_full_fits": review_full_fits,
        "review_full_rows": review_full_rows,
        "tar": tar_path,
        "n_catalog_rows": len(rows),
        "n_gifs_copied": copied,
        "n_gifs_missing": missing,
    }


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Package unknown-link GIFs and create a two-column human-review CSV.")
    ap.add_argument("night", help="Target night in YYYYMMDD format")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--root-out", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--plot-root", default="/pipeline/xiaoyunao/heliolincrr/plots")
    ap.add_argument("--out-root", default="/pipeline/xiaoyunao/heliolincrr/review_packages")
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--catalog", default="", help="Default: /processed1/<night>/L4/<night>_unknown_links.json")
    ap.add_argument("--make-tar", action="store_true", help="Also write <night>_unknown_review.tar.gz under out-root")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    print(json.dumps(build_package(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
