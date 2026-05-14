#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import tarfile
from pathlib import Path
from typing import Any


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
        "tar": tar_path,
        "n_catalog_rows": len(rows),
        "n_gifs_copied": copied,
        "n_gifs_missing": missing,
    }


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Package unknown-link GIFs and create a two-column human-review CSV.")
    ap.add_argument("night", help="Target night in YYYYMMDD format")
    ap.add_argument("--processed-root", default="/processed1")
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
