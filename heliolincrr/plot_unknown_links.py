#!/usr/bin/env python3
from __future__ import annotations

import argparse
import io
import json
import math
import re
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval
from astropy.wcs import FITSFixedWarning, WCS
from PIL import Image


plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.size"] = 12
warnings.filterwarnings("ignore", category=FITSFixedWarning)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Visualize fit_ok unknown links as cutout GIFs.")
    ap.add_argument("night", help="Target night in YYYYMMDD format")
    ap.add_argument("--processed-root", default="/processed1")
    ap.add_argument("--root-out", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--plot-root", default="/pipeline/xiaoyunao/heliolincrr/plots")
    ap.add_argument("--rr-subdir", default="rr_links")
    ap.add_argument("--catalog", default=None, help="Optional unknown-link catalog JSON (default: /processed1/<night>/L4/<night>_unknown_links.json)")
    ap.add_argument("--gif-size", type=int, default=280)
    ap.add_argument("--gif-duration", type=float, default=0.6)
    ap.add_argument("--limit-links", type=int, default=0)
    return ap.parse_args()


def extract_cutout(data: np.ndarray, x_center: float, y_center: float, size: int) -> tuple[np.ndarray, int, int]:
    half = size // 2
    xc = int(round(x_center))
    yc = int(round(y_center))
    x1 = xc - half
    x2 = xc + half
    y1 = yc - half
    y2 = yc + half
    xs1 = max(0, x1)
    xs2 = min(data.shape[1], x2)
    ys1 = max(0, y1)
    ys2 = min(data.shape[0], y2)
    out = np.full((size, size), np.nan, dtype=float)
    cx1 = xs1 - x1
    cy1 = ys1 - y1
    cx2 = cx1 + max(0, xs2 - xs1)
    cy2 = cy1 + max(0, ys2 - ys1)
    if xs2 > xs1 and ys2 > ys1:
        out[cy1:cy2, cx1:cx2] = data[ys1:ys2, xs1:xs2]
    return out, x1, y1


def draw_hollow_cross(ax: plt.Axes, x: float, y: float, color: str = "#00ff66") -> None:
    arm = 22
    gap = 8
    lw = 1.8
    ax.plot([x - arm, x - gap], [y, y], color=color, lw=lw)
    ax.plot([x + gap, x + arm], [y, y], color=color, lw=lw)
    ax.plot([x, x], [y - arm, y - gap], color=color, lw=lw)
    ax.plot([x, x], [y + gap, y + arm], color=color, lw=lw)


def buffer_to_array(fig: plt.Figure) -> np.ndarray:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight", pad_inches=0.04, facecolor="black")
    buf.seek(0)
    with Image.open(buf) as img:
        return np.asarray(img.convert("RGB"))


def save_gif(path: Path, frames: list[np.ndarray], duration: float) -> None:
    pil_frames = [Image.fromarray(frame) for frame in frames]
    pil_frames[0].save(
        path,
        save_all=True,
        append_images=pil_frames[1:],
        duration=max(1, int(round(duration * 1000))),
        loop=1,
    )


def parse_unknown_catalog(catalog_path: Path, processed_root: Path, night: str, limit_links: int) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    rows = json.loads(catalog_path.read_text(encoding="utf-8"))
    if limit_links > 0:
        rows = rows[:limit_links]
    det_rows: list[dict[str, object]] = []
    link_meta: dict[str, dict[str, object]] = {}
    for row in rows:
        lid = int(row["linkage_id"])
        image_names = [x for x in str(row.get("image_names", "")).split(";") if x]
        objids = [int(x) for x in str(row.get("objids", "")).split(";") if x]
        mjds = [float(x) for x in str(row.get("mjds", "")).split(";") if x]
        ras = [float(x) for x in str(row.get("ras_deg", "")).split(";") if x]
        decs = [float(x) for x in str(row.get("decs_deg", "")).split(";") if x]
        groups = [int(x) for x in str(row.get("groups", "")).split(";") if x]
        exp_pairs = [x for x in str(row.get("exp_pairs", "")).split(";") if x]
        for idx, image_name in enumerate(image_names):
            exp_i, exp_j = (-1, -1)
            if idx < len(exp_pairs) and "->" in exp_pairs[idx]:
                a, b = exp_pairs[idx].split("->", 1)
                exp_i, exp_j = int(a), int(b)
            det_rows.append(
                {
                    "linkage_id": lid,
                    "tracklet_id": "",
                    "endpoint": idx + 1,
                    "image_name": image_name,
                    "image_path": str(processed_root / night / "L1" / image_name),
                    "objID": objids[idx] if idx < len(objids) else -1,
                    "mjd": mjds[idx] if idx < len(mjds) else np.nan,
                    "ra_deg": ras[idx] if idx < len(ras) else np.nan,
                    "dec_deg": decs[idx] if idx < len(decs) else np.nan,
                    "group": groups[idx] if idx < len(groups) else -1,
                    "exp_i": exp_i,
                    "exp_j": exp_j,
                }
            )
        link_meta[str(lid)] = dict(row)
    return det_rows, link_meta


def build_gifs(
    rows: list[dict[str, object]],
    out_dir: Path,
    night: str,
    gif_size: int,
    duration: float,
    link_meta: dict[str, dict[str, object]],
) -> int:
    grouped: dict[int, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[int(row["linkage_id"])].append(row)

    n_written = 0
    summary_rows: list[dict[str, object]] = []
    for lid in sorted(grouped):
        dets = sorted(grouped[lid], key=lambda item: (item["mjd"], item["tracklet_id"], item["endpoint"]))
        mid = dets[len(dets) // 2]
        center_world = np.array([[mid["ra_deg"], mid["dec_deg"]]], dtype=float)
        frames = []
        missing_files = 0
        for row in dets:
            image_path = Path(str(row["image_path"]))
            if not image_path.exists():
                missing_files += 1
                continue
            with fits.open(image_path, memmap=False) as hdul:
                data = hdul[1].data.astype(float)
                wcs = WCS(hdul[1].header)
                center_pix = wcs.all_world2pix(center_world, 1)[0]
                marker_pix = wcs.all_world2pix(np.array([[row["ra_deg"], row["dec_deg"]]], dtype=float), 1)[0]
            cutout, x_left, y_bottom = extract_cutout(data, float(center_pix[0]), float(center_pix[1]), gif_size)
            valid = np.isfinite(cutout)
            vmin, vmax = ZScaleInterval().get_limits(cutout[valid]) if np.any(valid) else (0.0, 1.0)
            fig, ax = plt.subplots(figsize=(4.4, 4.4), facecolor="black")
            ax.imshow(cutout, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
            ax.axis("off")
            mx = float(marker_pix[0]) - x_left
            my = float(marker_pix[1]) - y_bottom
            if 0.0 <= mx < gif_size and 0.0 <= my < gif_size:
                draw_hollow_cross(ax, mx, my)
            ax.text(
                0.03,
                0.96,
                f"Unknown Link {lid}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=11,
                color="white",
                bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none"},
            )
            ax.text(
                0.03,
                0.10,
                f"{row['image_name']} | objID={row['objID']}",
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=9,
                color="#cdeffd",
                bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none"},
            )
            ax.text(
                0.03,
                0.04,
                f"group={row['group']} exp={row['exp_i']}->{row['exp_j']} mjd={row['mjd']:.6f}",
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=9,
                color="#cdeffd",
                bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none"},
            )
            frames.append(buffer_to_array(fig))
            plt.close(fig)
        if not frames:
            continue
        frames = frames + [frames[-1], frames[-1]]
        gif_name = f"unknown_link_{lid:04d}_{night}.gif"
        save_gif(out_dir / gif_name, frames, duration=duration)
        summary = dict(link_meta.get(str(lid), {}))
        summary.update(
            {
                "linkage_id": lid,
                "gif": gif_name,
                "n_frames": len(frames) - 2,
                "missing_files": missing_files,
                "image_names": [row["image_name"] for row in dets],
                "tracklet_ids": sorted({row["tracklet_id"] for row in dets}),
            }
        )
        summary_rows.append(summary)
        n_written += 1

    summary_path = out_dir / f"{night}_unknown_link_summary.json"
    summary_path.write_text(json.dumps(summary_rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return n_written


def main() -> None:
    args = parse_args()
    night = args.night
    processed_root = Path(args.processed_root)
    root_out = Path(args.root_out)
    out_dir = Path(args.plot_root) / night
    out_dir.mkdir(parents=True, exist_ok=True)
    catalog_path = Path(args.catalog) if args.catalog else (processed_root / night / "L4" / f"{night}_unknown_links.json")

    if not (processed_root / night / "L1").exists():
        raise SystemExit(f"Missing L1 dir: {processed_root / night / 'L1'}")
    if not catalog_path.exists():
        raise SystemExit(f"Missing unknown-link catalog: {catalog_path}")

    rows, link_meta = parse_unknown_catalog(
        catalog_path=catalog_path,
        processed_root=processed_root,
        night=night,
        limit_links=int(args.limit_links),
    )
    if not rows:
        raise SystemExit("No fit_ok all_non_asteroid links found.")

    n_written = build_gifs(
        rows=rows,
        out_dir=out_dir,
        night=night,
        gif_size=int(args.gif_size),
        duration=float(args.gif_duration),
        link_meta=link_meta,
    )
    print(f"[write] {out_dir}  n_gifs={n_written}")


if __name__ == "__main__":
    main()
