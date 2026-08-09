#!/usr/bin/env python3
"""Generate Fig. 6 from two frozen unknown-link review GIFs.

Only two frozen inputs are read:

* ``snapshot/review_sample/review_gifs/``;
* ``snapshot/derived_unknown/unknown_all_links.csv``.

The three frames from each selected GIF are converted to grayscale, the baked
green review marker is removed, a 61-pixel median background is subtracted,
and one shared median/MAD stretch is applied to all six panels.  This script
does not use a per-frame ZScale.  The compact flow at right is explicitly
pipeline context; it does not claim that the displayed GIF frames are archived
intermediate images from those stages.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
from PIL import Image, ImageSequence
from scipy.ndimage import binary_dilation, median_filter

from figure_styles import (
    CUTOUT_CROSSHAIR_COLORS,
    CUTOUT_GAPPED_CROSSHAIR_STYLE,
    STATISTICS_COLORS,
    WORKFLOW_COLORS,
    add_cutout_label,
    apply_statistics_style,
    clear_cutout_axis,
    cutout_imshow_kwargs,
    draw_gapped_crosshair,
    flatten_cutout_background,
    save_pdf_png,
    shared_cutout_limits,
)


SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parent
DEFAULT_GIF_DIR = ANALYSIS_ROOT / "snapshot" / "review_sample" / "review_gifs"
DEFAULT_LINK_TABLE = ANALYSIS_ROOT / "snapshot" / "derived_unknown" / "unknown_all_links.csv"
DEFAULT_OUTPUT_STEM = ANALYSIS_ROOT / "figures" / "fig06_unknown_pipeline_examples"

RETAINED_KEY = ("20251121", "000004r", "598")
REJECTED_KEY = ("20260529", "00001gB", "128")
RETAINED_GIF = "20251121_000004r_link0598.gif"
REJECTED_GIF = "20260529_00001gB_link0128.gif"

FRAME_COUNT = 3
CROP_SIZE = 320


@dataclass(frozen=True)
class Example:
    key: tuple[str, str, str]
    gif_name: str
    status_title: str
    marker_color: str
    row: dict[str, str]
    frames: tuple[np.ndarray, ...]
    marker_centers: tuple[tuple[float, float], ...]
    selected_indices: tuple[int, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gif-dir", type=Path, default=DEFAULT_GIF_DIR)
    parser.add_argument("--link-table", type=Path, default=DEFAULT_LINK_TABLE)
    parser.add_argument("--output-stem", type=Path, default=DEFAULT_OUTPUT_STEM)
    return parser.parse_args()


def load_selected_rows(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    """Read the two required linkage rows from the frozen aggregate table."""

    if not path.is_file():
        raise FileNotFoundError(f"frozen unknown-link table not found: {path}")
    wanted = {RETAINED_KEY, REJECTED_KEY}
    selected: dict[tuple[str, str, str], dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required_columns = {
            "night",
            "trk_sub",
            "linkage_id",
            "n_obs",
            "mjds",
            "rms_arcsec",
            "speed_arcsec_per_hour",
            "median_mag_aper4",
            "final_paper_status",
            "posthoc_retained",
        }
        missing = required_columns.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"unknown-link table lacks required columns: {sorted(missing)}")
        for row in reader:
            key = (row["night"], row["trk_sub"], row["linkage_id"])
            if key in wanted:
                if key in selected:
                    raise ValueError(f"duplicate selected linkage row: {key}")
                selected[key] = row
    missing_rows = wanted.difference(selected)
    if missing_rows:
        raise ValueError(f"selected linkages absent from frozen table: {sorted(missing_rows)}")

    retained = selected[RETAINED_KEY]
    rejected = selected[REJECTED_KEY]
    if retained["final_paper_status"] != "retained_after_posthoc_audit" or retained["posthoc_retained"] != "True":
        raise ValueError("retained example does not have the expected frozen audit status")
    if rejected["final_paper_status"] != "rejected_posthoc" or rejected["posthoc_retained"] != "False":
        raise ValueError("rejected example does not have the expected frozen audit status")
    return selected


def green_marker_mask(rgb: np.ndarray) -> np.ndarray:
    """Identify the baked green gapped crosshair in a frozen review GIF frame."""

    red = rgb[..., 0].astype(np.int16)
    green = rgb[..., 1].astype(np.int16)
    blue = rgb[..., 2].astype(np.int16)
    mask = (green > 120) & ((green - red) > 45) & ((green - blue) > 20)
    if int(np.count_nonzero(mask)) < 40:
        raise ValueError("could not identify the frozen GIF review crosshair")
    return mask


def marker_center(mask: np.ndarray) -> tuple[float, float]:
    """Estimate the crosshair center from the symmetric marker extent."""

    yy, xx = np.where(mask)
    return (0.5 * (float(xx.min()) + float(xx.max())), 0.5 * (float(yy.min()) + float(yy.max())))


def remove_marker_and_crop(rgb: np.ndarray) -> tuple[np.ndarray, tuple[float, float]]:
    """Remove the review overlay and return a target-centered square grayscale crop."""

    mask = green_marker_mask(rgb)
    center_x, center_y = marker_center(mask)

    # Rec. 709 luminance; the source is an 8-bit review rendering, not raw CCD data.
    gray = (
        0.2126 * rgb[..., 0].astype(float)
        + 0.7152 * rgb[..., 1].astype(float)
        + 0.0722 * rgb[..., 2].astype(float)
    )
    expanded_mask = binary_dilation(mask, iterations=2)
    local_background = median_filter(gray, size=11, mode="nearest")
    gray[expanded_mask] = local_background[expanded_mask]

    half = CROP_SIZE // 2
    center_x_int = int(round(center_x))
    center_y_int = int(round(center_y))
    x0 = center_x_int - half
    y0 = center_y_int - half
    x1 = x0 + CROP_SIZE
    y1 = y0 + CROP_SIZE
    if x0 < 0 or y0 < 0 or x1 > gray.shape[1] or y1 > gray.shape[0]:
        raise ValueError("target-centered GIF crop would extend outside the frozen frame")
    crop = gray[y0:y1, x0:x1]
    center_in_crop = (center_x - x0, center_y - y0)
    return crop, center_in_crop


def select_frame_indices(n_frames: int, count: int = FRAME_COUNT) -> tuple[int, ...]:
    if n_frames < count:
        raise ValueError(f"review GIF has {n_frames} frames; {count} are required")
    if n_frames == count:
        return tuple(range(count))
    return tuple(int(value) for value in np.linspace(0, n_frames - 1, count).round())


def load_gif_frames(path: Path) -> tuple[tuple[np.ndarray, ...], tuple[tuple[float, float], ...], tuple[int, ...]]:
    """Load exactly three target-centered frames from one frozen review GIF."""

    if not path.is_file():
        raise FileNotFoundError(f"frozen review GIF not found: {path}")
    with Image.open(path) as image:
        rgb_frames = [np.asarray(frame.convert("RGB"), dtype=np.uint8) for frame in ImageSequence.Iterator(image)]
    selected_indices = select_frame_indices(len(rgb_frames))
    crops: list[np.ndarray] = []
    centers: list[tuple[float, float]] = []
    for index in selected_indices:
        crop, center = remove_marker_and_crop(rgb_frames[index])
        crops.append(flatten_cutout_background(crop))
        centers.append(center)
    return tuple(crops), tuple(centers), selected_indices


def make_example(
    gif_dir: Path,
    row: dict[str, str],
    *,
    key: tuple[str, str, str],
    gif_name: str,
    status_title: str,
    marker_color: str,
) -> Example:
    frames, centers, indices = load_gif_frames(gif_dir / gif_name)
    return Example(key, gif_name, status_title, marker_color, row, frames, centers, indices)


def workflow_rgba(name: str) -> tuple[float, float, float, float]:
    return tuple(component / 255.0 for component in WORKFLOW_COLORS[name])


def draw_context_flow(ax) -> None:
    """Draw the contextual pipeline labels, clearly separated from the images."""

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.5,
        0.965,
        "PIPELINE CONTEXT",
        ha="center",
        va="top",
        fontsize=17,
        fontweight="bold",
        color=workflow_rgba("ink"),
    )

    stages = [
        (0.80, "L2 source catalogs", "purple"),
        (0.61, "Gaia + static-source\nfiltering", "green"),
        (0.42, "Tracklet + shared-endpoint\nlink construction", "green"),
        (0.23, "Review package\n+ human audit", "orange"),
    ]
    boxes: list[tuple[float, float]] = []
    for y, label, color_kind in stages:
        x0, width, height = 0.10, 0.80, 0.115
        box = FancyBboxPatch(
            (x0, y - height / 2),
            width,
            height,
            boxstyle="round,pad=0.012,rounding_size=0.025",
            facecolor=workflow_rgba(f"{color_kind}_fill"),
            edgecolor=workflow_rgba(f"{color_kind}_line"),
            linewidth=1.8,
        )
        ax.add_patch(box)
        ax.text(
            0.5,
            y,
            label,
            ha="center",
            va="center",
            fontsize=12.0,
            fontweight="bold",
            color=workflow_rgba("ink"),
            linespacing=1.1,
        )
        boxes.append((y - height / 2, y + height / 2))

    for (source_bottom, _), (_, target_top) in zip(boxes[:-1], boxes[1:]):
        arrow = FancyArrowPatch(
            (0.5, source_bottom),
            (0.5, target_top),
            arrowstyle="-|>",
            mutation_scale=15,
            linewidth=1.6,
            color=workflow_rgba("arrow"),
            shrinkA=2,
            shrinkB=2,
        )
        ax.add_patch(arrow)

    note = (
        "Context only: the six panels are three epochs\n"
        "from each frozen review GIF. They are not\n"
        "independently archived images from the four\n"
        "processing stages."
    )
    ax.text(
        0.5,
        0.055,
        note,
        ha="center",
        va="bottom",
        fontsize=9.4,
        color=STATISTICS_COLORS["mid_grey"],
        linespacing=1.25,
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": "#f7f7f7",
            "edgecolor": "#b5b5b5",
            "linewidth": 0.8,
        },
    )


def metric_line(example: Example) -> str:
    row = example.row
    return (
        fr"$N_{{\mathrm{{obs}}}}={int(float(row['n_obs']))}$   |   "
        f"speed = {float(row['speed_arcsec_per_hour']):.1f} arcsec h$^{{-1}}$   |   "
        fr"median $g_{{\mathrm{{aper}}}}$ = {float(row['median_mag_aper4']):.2f}   |   "
        f"fit RMS = {float(row['rms_arcsec']):.2f} arcsec"
    )


def relative_hours(example: Example) -> tuple[float, ...]:
    mjds = tuple(float(value) for value in example.row["mjds"].split(";") if value)
    if len(mjds) < FRAME_COUNT:
        raise ValueError(f"selected linkage has only {len(mjds)} stored epochs")
    chosen = tuple(mjds[index] for index in example.selected_indices)
    first = chosen[0]
    return tuple((value - first) * 24.0 for value in chosen)


def make_figure(examples: Sequence[Example]):
    """Build the two-row cutout montage and contextual flow in memory."""

    if len(examples) != 2:
        raise ValueError("Fig. 6 requires exactly two linkage examples")
    all_frames = [frame for example in examples for frame in example.frames]
    vmin, vmax = shared_cutout_limits(all_frames)

    apply_statistics_style()
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(18, 10.5), facecolor="white")
    grid = fig.add_gridspec(
        2,
        4,
        width_ratios=(1.0, 1.0, 1.0, 0.90),
        left=0.045,
        right=0.985,
        bottom=0.075,
        top=0.775,
        wspace=0.055,
        hspace=0.35,
    )

    fig.text(
        0.5,
        0.955,
        "Representative Unknown-Link Review Sequences",
        ha="center",
        va="center",
        fontsize=27,
    )
    fig.text(
        0.5,
        0.918,
        "Median-filtered backgrounds and one shared grayscale stretch across all six displayed epochs",
        ha="center",
        va="center",
        fontsize=15,
        color=STATISTICS_COLORS["mid_grey"],
    )

    header_y = (0.845, 0.425)
    metric_y = (0.812, 0.392)
    panel_letters = ("a", "b")
    for row_index, example in enumerate(examples):
        date_label = f"{example.key[0][:4]}-{example.key[0][4:6]}-{example.key[0][6:]}"
        fig.text(
            0.052,
            header_y[row_index],
            f"({panel_letters[row_index]}) {example.status_title} — {date_label}, "
            f"trkSub {example.key[1]}, link {int(example.key[2])}",
            ha="left",
            va="center",
            fontsize=16.5,
            fontweight="bold",
            color=example.marker_color,
        )
        fig.text(
            0.052,
            metric_y[row_index],
            metric_line(example),
            ha="left",
            va="center",
            fontsize=13.0,
            color=STATISTICS_COLORS["ink"],
        )

        delta_hours = relative_hours(example)
        for frame_index, (frame, center, delta_hour) in enumerate(
            zip(example.frames, example.marker_centers, delta_hours)
        ):
            ax = fig.add_subplot(grid[row_index, frame_index])
            ax.imshow(frame, **cutout_imshow_kwargs(vmin, vmax))
            draw_gapped_crosshair(
                ax,
                center[0],
                center[1],
                color=example.marker_color,
                **CUTOUT_GAPPED_CROSSHAIR_STYLE,
            )
            add_cutout_label(
                ax,
                8,
                CROP_SIZE - 8,
                f"Exposure {frame_index + 1}\n" + fr"$\Delta t$ = {delta_hour:.2f} h",
                fontsize=11.0,
                ha="left",
                va="top",
            )
            clear_cutout_axis(ax)

    context_axis = fig.add_subplot(grid[:, 3])
    draw_context_flow(context_axis)
    return fig, (vmin, vmax)


def verify_outputs(png_path: Path, pdf_path: Path) -> dict[str, Any]:
    for path in (png_path, pdf_path):
        if not path.is_file() or path.stat().st_size < 20_000:
            raise RuntimeError(f"Fig. 6 output is missing or unexpectedly small: {path}")
    with Image.open(png_path) as image:
        width, height = image.size
        dpi = image.info.get("dpi", (0.0, 0.0))
    if width < 4000 or height < 2200 or width <= height:
        raise RuntimeError(f"unexpected Fig. 6 PNG dimensions: {width}x{height}")
    if len(dpi) != 2 or any(abs(float(value) - 300.0) > 2.0 for value in dpi):
        raise RuntimeError(f"Fig. 6 PNG lacks 300-dpi metadata: {dpi}")
    with pdf_path.open("rb") as handle:
        if handle.read(5) != b"%PDF-":
            raise RuntimeError(f"invalid Fig. 6 PDF signature: {pdf_path}")
    return {
        "width": width,
        "height": height,
        "dpi": tuple(round(float(value), 2) for value in dpi),
        "png_bytes": png_path.stat().st_size,
        "pdf_bytes": pdf_path.stat().st_size,
    }


def main() -> None:
    args = parse_args()
    gif_dir = args.gif_dir.resolve()
    if not gif_dir.is_dir():
        raise FileNotFoundError(f"frozen review-GIF directory not found: {gif_dir}")
    rows = load_selected_rows(args.link_table.resolve())
    examples = (
        make_example(
            gif_dir,
            rows[RETAINED_KEY],
            key=RETAINED_KEY,
            gif_name=RETAINED_GIF,
            status_title="Retained after post-hoc audit",
            marker_color=CUTOUT_CROSSHAIR_COLORS[1],
        ),
        make_example(
            gif_dir,
            rows[REJECTED_KEY],
            key=REJECTED_KEY,
            gif_name=REJECTED_GIF,
            status_title="Post-hoc rejected artifact candidate",
            marker_color=CUTOUT_CROSSHAIR_COLORS[2],
        ),
    )
    figure, limits = make_figure(examples)
    png_path, pdf_path = save_pdf_png(
        figure,
        args.output_stem.resolve(),
        png_dpi=300,
        tight=True,
        pad_inches=0.05,
    )
    validation = verify_outputs(png_path, pdf_path)
    print(f"link table: {args.link_table.resolve()}")
    print(f"GIFs: {RETAINED_GIF}, {REJECTED_GIF}")
    print(f"shared display limits: {limits[0]:.6g}, {limits[1]:.6g}")
    print(f"PNG: {png_path} ({validation['width']}x{validation['height']}, {validation['dpi']} dpi)")
    print(f"PDF: {pdf_path} ({validation['pdf_bytes']} bytes)")


if __name__ == "__main__":
    main()
