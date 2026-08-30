#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import math

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter


HERE = Path(__file__).resolve().parent
SOURCE_DIR = HERE / "source_data"
IMAGE_PATH = SOURCE_DIR / "stpxl-0592_20250204_0001_3.fits.gz"
MATCH_PATH = SOURCE_DIR / "20250204_matched_asteroids.fits"
LOCAL_TOTAL = HERE.parent / "gotta_asteroids.fits"
SOURCE_FILE = "stpxl-0592_20250204_0001_3_cat.fits.gz"
IMAGE_LABEL = "stpxl-0592_20250204_0001_3"
ADOPTED_MAG = "Mag_Aper4"

DRAFT_TARGET_INDICES = [2, 9, 12, 21]


def as_str(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore")
    return str(value)


def robust_limits(data: np.ndarray, lo: float = 0.5, hi: float = 99.7) -> tuple[float, float]:
    finite = np.asarray(data, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 0.0, 1.0
    vmin, vmax = np.nanpercentile(finite, [lo, hi])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        med = float(np.nanmedian(finite))
        std = float(np.nanstd(finite)) or 1.0
        return med - std, med + 5.0 * std
    return float(vmin), float(vmax)


def log_stretch(data: np.ndarray, vmin: float, vmax: float, scale: float = 1200.0) -> np.ndarray:
    arr = np.asarray(data, dtype=float)
    norm = np.clip((arr - vmin) / (vmax - vmin), 0.0, 1.0)
    norm[~np.isfinite(norm)] = 0.0
    return np.log1p(scale * norm) / np.log1p(scale)


def linear_stretch(data: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    arr = np.asarray(data, dtype=float)
    norm = np.clip((arr - vmin) / (vmax - vmin), 0.0, 1.0)
    norm[~np.isfinite(norm)] = 0.0
    return norm


def zscale_limits(data: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(data, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 0.0, 1.0
    interval = ZScaleInterval(contrast=0.25, krej=2.5)
    return tuple(float(x) for x in interval.get_limits(finite))


def oc_separation_arcsec(row) -> float:
    dra = (float(row["RA_Win"]) - float(row["ra"])) * np.cos(np.deg2rad(float(row["dec"]))) * 3600.0
    ddec = (float(row["DEC_Win"]) - float(row["dec"])) * 3600.0
    return float(np.hypot(dra, ddec))


def object_id(row) -> str:
    number = row["number"]
    if np.ma.is_masked(number):
        return as_str(row["name"])
    try:
        value = float(number)
    except Exception:
        return as_str(row["name"])
    if not np.isfinite(value):
        return as_str(row["name"])
    return f"({int(round(value))}) {as_str(row['name'])}"


def source_rows() -> list[dict]:
    matched = Table.read(MATCH_PATH)
    total = Table.read(LOCAL_TOTAL, memmap=True)
    sub = matched[matched["source_file"] == SOURCE_FILE]
    sub.sort(ADOPTED_MAG)

    local_names = np.array([as_str(x) for x in total["name"]])
    local_epoch = np.asarray(total["epoch"], dtype=float)
    rows = []
    for idx, row in enumerate(sub, 1):
        name = as_str(row["name"])
        local_mask = (local_names == name) & (np.abs(local_epoch - float(row["epoch"])) < 0.01)
        local = total[local_mask]
        rate = float(local[0]["ang_rate_arcsec_hour"]) if len(local) else np.nan
        rows.append(
            {
                "idx": idx,
                "name": name,
                "object_id": object_id(row),
                "x": float(row["X_Win"]),
                "y": float(row["Y_Win"]),
                "g_aper": float(row[ADOPTED_MAG]),
                "rate": rate,
                "sep": oc_separation_arcsec(row),
            }
        )
    return rows


def cutout(data: np.ndarray, x: float, y: float, size: int = 180) -> np.ndarray:
    ny, nx = data.shape
    half = size // 2
    cx = int(round(x)) - 1
    cy = int(round(y)) - 1
    x0, x1 = cx - half, cx + half
    y0, y1 = cy - half, cy + half
    out = np.full((size, size), np.nan, dtype=float)
    sx0, sx1 = max(0, x0), min(nx, x1)
    sy0, sy1 = max(0, y0), min(ny, y1)
    dx0, dy0 = sx0 - x0, sy0 - y0
    if sx1 > sx0 and sy1 > sy0:
        out[dy0 : dy0 + sy1 - sy0, dx0 : dx0 + sx1 - sx0] = data[sy0:sy1, sx0:sx1]
    return out


def draw_label(ax, x: float, y: float, text: str, size: float = 10.0, ha: str = "left", va: str = "top") -> None:
    ax.text(
        x,
        y,
        text,
        ha=ha,
        va=va,
        color="white",
        fontsize=size,
        fontweight="bold",
        linespacing=1.15,
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "black", "edgecolor": "none", "alpha": 0.64},
    )


def clear_frame(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_solid_crosshair(
    ax,
    x: float,
    y: float,
    color: str,
    half_length: float = 17.0,
    linewidth: float = 2.0,
) -> None:
    ax.plot([x - half_length, x + half_length], [y, y], color=color, linewidth=linewidth, solid_capstyle="butt")
    ax.plot([x, x], [y - half_length, y + half_length], color=color, linewidth=linewidth, solid_capstyle="butt")


def draw_gapped_crosshair(
    ax,
    x: float,
    y: float,
    color: str,
    gap: float = 11.0,
    length: float = 17.0,
    linewidth: float = 2.4,
) -> None:
    for x0, x1, y0, y1 in (
        (x - gap - length, x - gap, y, y),
        (x + gap, x + gap + length, y, y),
        (x, x, y - gap - length, y - gap),
        (x, x, y + gap, y + gap + length),
    ):
        ax.plot([x0, x1], [y0, y1], color=color, linewidth=linewidth, solid_capstyle="butt")


def marker_inside_frame(row: dict, nx: int, ny: int, half_length: float) -> bool:
    return half_length <= row["x"] <= nx - half_length and half_length <= row["y"] <= ny - half_length


def rate_text(rate: float) -> str:
    return "--" if not np.isfinite(rate) else f"{rate:.1f}"


def load_image() -> tuple[np.ndarray, float]:
    with fits.open(IMAGE_PATH, memmap=True) as hdul:
        data = np.asarray(hdul["IMG"].data, dtype=float)
        mjd = float(hdul["IMG"].header["MJD"])
    return data, mjd


def full_frame_view(data: np.ndarray, factor: int = 6) -> tuple[np.ndarray, tuple[int, int]]:
    ny, nx = data.shape
    view = data[: ny - ny % factor, : nx - nx % factor].reshape(ny // factor, factor, nx // factor, factor).mean(axis=(1, 3))
    return view, (nx, ny)


def flatten_cutout_background(cut: np.ndarray, filter_size: int = 61) -> np.ndarray:
    arr = np.asarray(cut, dtype=float)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return arr
    fill = float(np.nanmedian(arr[finite]))
    filled = np.where(finite, arr, fill)
    background = median_filter(filled, size=filter_size, mode="nearest")
    residual = arr - background
    residual[~finite] = np.nan
    return residual


def display_cutouts(cuts: list[np.ndarray]) -> list[np.ndarray]:
    return [flatten_cutout_background(cut) for cut in cuts]


def shared_cutout_limits(cuts: list[np.ndarray]) -> tuple[float, float]:
    vals = np.concatenate([np.asarray(cut, dtype=float).ravel() for cut in cuts])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 0.0, 1.0
    center = float(np.nanmedian(vals))
    mad = float(np.nanmedian(np.abs(vals - center)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.nanstd(vals))
    if not np.isfinite(sigma) or sigma <= 0:
        return robust_limits(vals)
    return center - 1.0 * sigma, center + 7.0 * sigma


def plot_four_target_draft(data: np.ndarray, mjd: float, rows: list[dict]) -> None:
    by_idx = {row["idx"]: row for row in rows}
    selected = [by_idx[idx] for idx in DRAFT_TARGET_INDICES]
    if len(selected) != 4:
        raise RuntimeError(f"Draft target selection returned {len(selected)} rows")

    view, (nx, ny) = full_frame_view(data)
    full_vmin, full_vmax = zscale_limits(view)
    selected_cuts = [cutout(data, row["x"], row["y"]) for row in selected]
    selected_display_cuts = display_cutouts(selected_cuts)
    cut_vmin, cut_vmax = shared_cutout_limits(selected_display_cuts)
    colors = ["#ffcc33", "#00b4d8", "#fb5607", "#80ed99"]

    fig = plt.figure(figsize=(15.2, 8.0), dpi=220)
    left = 0.015
    bottom = 0.035
    height = 0.93
    gap = 0.018
    cut_gap = 0.018
    left_width = 0.475
    cut_size = (1.0 - left - left_width - gap - 0.015 - cut_gap) / 2.0
    row_gap = cut_gap
    cut_height = (height - row_gap) / 2.0
    right_left = left + left_width + gap

    ax_full = fig.add_axes([left, bottom, left_width, height])
    ax_full.imshow(
        linear_stretch(view, full_vmin, full_vmax),
        cmap="gray",
        origin="lower",
        extent=[0, nx, 0, ny],
        interpolation="nearest",
    )
    clear_frame(ax_full)
    draw_label(ax_full, 130, ny - 150, f"{IMAGE_LABEL}\nMJD = {mjd:.8f}", size=11.5, va="top")

    selected_indices = {row["idx"] for row in selected}
    other_marker_half = 58.0
    for row in rows:
        if row["idx"] not in selected_indices and marker_inside_frame(row, nx, ny, other_marker_half):
            draw_solid_crosshair(ax_full, row["x"], row["y"], color="#2dd4bf", half_length=other_marker_half, linewidth=1.2)

    for row, color in zip(selected, colors):
        draw_gapped_crosshair(ax_full, row["x"], row["y"], color=color, gap=72, length=160, linewidth=1.5)
        label_x = min(max(row["x"] + 155, 90), nx - 1450)
        label_y = min(max(row["y"] + 95, 170), ny - 160)
        draw_label(ax_full, label_x, label_y, row["name"], size=9.8, va="center")

    for idx, (row, color) in enumerate(zip(selected, colors), 1):
        col = (idx - 1) % 2
        r = (idx - 1) // 2
        x0 = right_left + col * (cut_size + cut_gap)
        y0 = bottom + (1 - r) * (cut_height + row_gap)
        ax = fig.add_axes([x0, y0, cut_size, cut_height])
        cut = selected_display_cuts[idx - 1]
        ax.imshow(linear_stretch(cut, cut_vmin, cut_vmax), cmap="gray", origin="lower", interpolation="nearest")
        clear_frame(ax)
        cy = cx = cut.shape[0] / 2.0 - 0.5
        draw_gapped_crosshair(ax, cx, cy, color=color, gap=10, length=16, linewidth=2.4)
        label = (
            f"{row['object_id']}\n"
            f"$\\mathbf{{g}}_{{\\mathbf{{aper}}}}$ = {row['g_aper']:.2f} mag\n"
            f"$\\boldsymbol{{\\mu}}$ = {rate_text(row['rate'])} arcsec $\\mathbf{{hr}}^{{\\mathbf{{-1}}}}$\n"
            f"O-C = {row['sep']:.2f} arcsec"
        )
        draw_label(ax, 7, cut.shape[0] - 7, label, size=9.0, va="top")

    fig.savefig(HERE / "fig1_review_cutouts.png", bbox_inches="tight", pad_inches=0.0)
    fig.savefig(HERE / "fig1_review_cutouts.pdf", bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)


def plot_all_cutouts(data: np.ndarray, rows: list[dict]) -> None:
    ncols = 4
    nrows = math.ceil(len(rows) / ncols)
    fig = plt.figure(figsize=(14.6, 3.25 * nrows), dpi=220)
    margin_x = 0.012
    margin_y = 0.012
    gap = 0.010
    cell_w = (1.0 - 2 * margin_x - (ncols - 1) * gap) / ncols
    cell_h = (1.0 - 2 * margin_y - (nrows - 1) * gap) / nrows
    colors = ["#ffcc33", "#00b4d8", "#fb5607", "#80ed99", "#f15bb5", "#9b5de5", "#00f5d4", "#fee440"]
    all_cuts = [cutout(data, row["x"], row["y"]) for row in rows]
    all_display_cuts = display_cutouts(all_cuts)
    cut_vmin, cut_vmax = shared_cutout_limits(all_display_cuts)

    for i, row in enumerate(rows):
        col = i % ncols
        r = i // ncols
        x0 = margin_x + col * (cell_w + gap)
        y0 = 1.0 - margin_y - (r + 1) * cell_h - r * gap
        ax = fig.add_axes([x0, y0, cell_w, cell_h])
        cut = all_display_cuts[i]
        ax.imshow(linear_stretch(cut, cut_vmin, cut_vmax), cmap="gray", origin="lower", interpolation="nearest")
        clear_frame(ax)
        color = colors[i % len(colors)]
        cy = cx = cut.shape[0] / 2.0 - 0.5
        draw_solid_crosshair(ax, cx, cy, color=color, half_length=12.5, linewidth=1.8)
        label = (
            f"{row['idx']}. {row['object_id']}\n"
            f"$g_{{\\rm aper}}$={row['g_aper']:.2f}  "
            f"$\\mu$={rate_text(row['rate'])} arcsec hr$^{{-1}}$\n"
            f"O-C={row['sep']:.2f} arcsec"
        )
        draw_label(ax, 6, cut.shape[0] - 6, label, size=7.8, va="top")

    fig.savefig(HERE / "fig1_all_cutouts_selection.png", bbox_inches="tight", pad_inches=0.0)
    fig.savefig(HERE / "fig1_all_cutouts_selection.pdf", bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)


def main() -> None:
    rows = source_rows()
    data, mjd = load_image()
    plot_four_target_draft(data, mjd, rows)
    plot_all_cutouts(data, rows)


if __name__ == "__main__":
    main()
