#!/usr/bin/env python3
"""Semantic figure styles derived from ``gotta_asteroid_1``.

This module centralizes four visually distinct GOTTA figure families without
reading survey data or producing any paper figure on import:

1. Times New Roman statistics panels and PDF/PNG export;
2. Mollweide/HEALPix maps with the inset horizontal colorbar;
3. background-flattened cutouts with shared stretch and gapped crosshairs;
4. Arial Bold pastel workflow diagrams.

The statistics and map values follow ``gotta_asteroid_1``'s current working
tree.  In particular, the reviewer residual plots use a 0.42 opacity for their
16--84 percentile bands (the committed baseline used 0.22).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


# ---------------------------------------------------------------------------
# 1. Statistical figures
# ---------------------------------------------------------------------------

STATISTICS_FOUR_PANEL_FIGSIZE = (18.0, 14.0)
STATISTICS_TWO_PANEL_FIGSIZE = (18.0, 9.0)
STATISTICS_RATE_PANEL_FIGSIZE = (10.0, 5.3)
STATISTICS_LEGEND_FONTSIZE = 21
STATISTICS_TICK_FONTSIZE = 22
STATISTICS_COLORBAR_LABEL_FONTSIZE = 20
STATISTICS_COLORBAR_TICK_FONTSIZE = 17

STATISTICS_RCPARAMS: dict[str, Any] = {
    "font.family": "Times New Roman",
    "font.size": 30,
    "axes.linewidth": 1.4,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.major.size": 7,
    "ytick.major.size": 7,
    "legend.frameon": False,
    "figure.dpi": 150,
}

STATISTICS_COLORS: dict[str, str] = {
    "blue": "#4c78a8",
    "orange": "#f28e2b",
    "green": "#59a14f",
    "red": "#e15759",
    "cyan": "#76b7b2",
    "gold": "#edc948",
    "purple": "#b07aa1",
    "trend_red": "#d62728",
    "ink": "#2b2b2b",
    "mid_grey": "#7f7f7f",
    "light_grey": "#d9d9d9",
}

STATISTICS_HISTOGRAM_STYLE: dict[str, Any] = {
    "histtype": "stepfilled",
    "alpha": 0.5,
    "edgecolor": STATISTICS_COLORS["ink"],
    "linewidth": 0.55,
}
STATISTICS_DENSITY_STYLE: dict[str, Any] = {
    "cmap": "viridis",
    "marker": "o",
    "s": 3.0,
    "linewidths": 0,
    "alpha": 0.82,
    "rasterized": True,
}
STATISTICS_TREND_STYLE: dict[str, Any] = {
    "color": STATISTICS_COLORS["trend_red"],
    "linewidth": 2.2,
}
STATISTICS_PERCENTILE_BAND_STYLE: dict[str, Any] = {
    "color": STATISTICS_COLORS["trend_red"],
    "alpha": 0.42,
}
STATISTICS_HISTOGRAM_GRID: dict[str, float] = {"alpha": 0.16, "linewidth": 0.6}
STATISTICS_DENSITY_GRID: dict[str, float] = {"alpha": 0.12, "linewidth": 0.55}


def apply_statistics_style() -> None:
    """Apply the core rcParams used by GOTTA's current statistics figures."""

    plt.rcParams.update(STATISTICS_RCPARAMS)


def style_statistics_axis(
    ax,
    *,
    density: bool = False,
    tick_fontsize: float = STATISTICS_TICK_FONTSIZE,
) -> None:
    """Apply GOTTA grid and tick styling to a statistics axis."""

    ax.set_axisbelow(True)
    ax.grid(**(STATISTICS_DENSITY_GRID if density else STATISTICS_HISTOGRAM_GRID))
    ax.tick_params(labelsize=tick_fontsize)


def add_statistics_colorbar(ax, mappable, label: str):
    """Append GOTTA's narrow right-hand colorbar and return it."""

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3.5%", pad=0.06)
    colorbar = ax.figure.colorbar(mappable, cax=cax)
    colorbar.set_label(label, fontsize=STATISTICS_COLORBAR_LABEL_FONTSIZE)
    colorbar.ax.tick_params(
        labelsize=STATISTICS_COLORBAR_TICK_FONTSIZE,
        width=0.9,
        length=4,
    )
    colorbar.outline.set_linewidth(0.8)
    return colorbar


def save_pdf_png(
    fig,
    output_stem: str | Path,
    *,
    png_dpi: int = 300,
    tight: bool = True,
    pad_inches: float | None = None,
    close: bool = True,
) -> tuple[Path, Path]:
    """Save a 300-dpi PNG and vector PDF using the GOTTA export contract.

    ``output_stem`` may have a suffix; it is replaced with ``.png``/``.pdf``.
    No raster ``dpi`` is supplied to the PDF save call, preserving vector
    artists except for artists intentionally marked ``rasterized=True``.
    """

    stem = Path(output_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = stem.with_suffix(".png")
    pdf_path = stem.with_suffix(".pdf")

    save_kwargs: dict[str, Any] = {}
    if tight:
        save_kwargs["bbox_inches"] = "tight"
    if pad_inches is not None:
        save_kwargs["pad_inches"] = pad_inches

    fig.savefig(png_path, dpi=png_dpi, **save_kwargs)
    fig.savefig(pdf_path, **save_kwargs)
    if close:
        plt.close(fig)
    return png_path, pdf_path


# ---------------------------------------------------------------------------
# 2. Mollweide / HEALPix maps
# ---------------------------------------------------------------------------

HEALPIX_NSIDE = 64
HEALPIX_NEST = False
HEALPIX_FIGSIZE = (14.0, 8.0)
HEALPIX_CMAP = "rainbow"
HEALPIX_VMIN = 1.0
HEALPIX_VMAX = 100.0
HEALPIX_POINT_SIZE = 8.0
HEALPIX_GRID_ALPHA = 0.3
HEALPIX_CENTER_LONGITUDE_DEG = 180.0
HEALPIX_TICK_STEP_DEG = 60
HEALPIX_TICK_LABEL_LATITUDE_DEG = -15.0
HEALPIX_TICK_FONTSIZE = 24
HEALPIX_COLORBAR_BOUNDS = (0.30, 0.09, 0.40, 0.045)
HEALPIX_COLORBAR_LABEL_FONTSIZE = 18
HEALPIX_COLORBAR_TICK_FONTSIZE = 16
HEALPIX_COLORBAR_TICK_STYLE: dict[str, Any] = {
    "direction": "in",
    "length": 5,
    "width": 1.1,
}
HEALPIX_CELESTIAL_EQUATOR_STYLE: dict[str, Any] = {
    "color": "#ff8c00",
    "linewidth": 3.0,
    "alpha": 0.95,
    "zorder": 9,
}


def longitude_to_mollweide_rad(
    longitude_deg: Sequence[float] | np.ndarray | float,
    *,
    center_deg: float = HEALPIX_CENTER_LONGITUDE_DEG,
) -> np.ndarray:
    """Wrap 0--360 degree longitudes around ``center_deg`` for Mollweide."""

    longitude = np.asarray(longitude_deg, dtype=float)
    wrapped = (longitude - center_deg + 180.0) % 360.0 - 180.0
    return np.deg2rad(wrapped)


def set_mollweide_longitude_ticks(
    ax,
    *,
    center_deg: float = HEALPIX_CENTER_LONGITUDE_DEG,
    step_deg: int = HEALPIX_TICK_STEP_DEG,
    label_latitude_deg: float = HEALPIX_TICK_LABEL_LATITUDE_DEG,
    fontsize: float = HEALPIX_TICK_FONTSIZE,
):
    """Draw GOTTA's manually positioned 0--360 degree longitude labels."""

    labels = np.arange(step_deg, 360, step_deg, dtype=int)
    ticks = longitude_to_mollweide_rad(labels, center_deg=center_deg)
    ax.set_xticks(ticks)
    ax.set_xticklabels([])
    label_latitude = np.deg2rad(label_latitude_deg)
    return [
        ax.text(
            tick,
            label_latitude,
            f"{value:d}°",
            ha="center",
            va="center",
            fontsize=fontsize,
            color="black",
            zorder=10,
        )
        for value, tick in zip(labels, ticks)
    ]


def style_mollweide_healpix_axis(ax) -> None:
    """Apply GOTTA's map grid, ticks, and Times font to a Mollweide axis."""

    apply_statistics_style()
    ax.grid(True, alpha=HEALPIX_GRID_ALPHA)
    set_mollweide_longitude_ticks(ax)
    ax.tick_params(axis="y", labelsize=HEALPIX_TICK_FONTSIZE)


def healpix_scatter_kwargs(
    *,
    vmin: float = HEALPIX_VMIN,
    vmax: float = HEALPIX_VMAX,
) -> dict[str, Any]:
    """Return the scatter kwargs used by GOTTA's HEALPix Mollweide map."""

    cmap = plt.get_cmap(HEALPIX_CMAP).copy()
    cmap.set_bad(alpha=0.0)
    return {
        "s": HEALPIX_POINT_SIZE,
        "linewidths": 0,
        "cmap": cmap,
        "norm": LogNorm(vmin=vmin, vmax=vmax),
        "rasterized": True,
    }


def add_mollweide_inset_colorbar(
    fig,
    ax,
    mappable,
    *,
    label: str = r"$\mathrm{Counts\ per\ pixel}$",
):
    """Add GOTTA's horizontal in-axis HEALPix colorbar and return it."""

    cax = ax.inset_axes(HEALPIX_COLORBAR_BOUNDS)
    colorbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
    colorbar.set_label(label, fontsize=HEALPIX_COLORBAR_LABEL_FONTSIZE)
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.ax.tick_params(
        labelsize=HEALPIX_COLORBAR_TICK_FONTSIZE,
        **HEALPIX_COLORBAR_TICK_STYLE,
    )
    return colorbar


# ---------------------------------------------------------------------------
# 3. Image cutouts
# ---------------------------------------------------------------------------

CUTOUT_MEDIAN_FILTER_SIZE = 61
CUTOUT_SHARED_LOWER_SIGMA = 1.0
CUTOUT_SHARED_UPPER_SIGMA = 7.0
CUTOUT_CMAP = "gray"
CUTOUT_INTERPOLATION = "nearest"
CUTOUT_CROSSHAIR_COLORS = ("#ffcc33", "#00b4d8", "#fb5607", "#80ed99")
CUTOUT_SUPPLEMENTARY_COLOR = "#2dd4bf"
CUTOUT_GAPPED_CROSSHAIR_STYLE: dict[str, float] = {
    "gap": 10.0,
    "length": 16.0,
    "linewidth": 2.4,
}
FULL_FRAME_GAPPED_CROSSHAIR_STYLE: dict[str, float] = {
    "gap": 72.0,
    "length": 160.0,
    "linewidth": 1.5,
}
CUTOUT_LABEL_STYLE: dict[str, Any] = {
    "color": "white",
    "fontweight": "bold",
    "linespacing": 1.15,
    "bbox": {
        "boxstyle": "round,pad=0.22",
        "facecolor": "black",
        "edgecolor": "none",
        "alpha": 0.64,
    },
}


def robust_display_limits(
    data: np.ndarray,
    *,
    lower_percentile: float = 0.5,
    upper_percentile: float = 99.7,
) -> tuple[float, float]:
    """Return the percentile fallback used by the GOTTA cutout workflow."""

    finite = np.asarray(data, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 0.0, 1.0
    vmin, vmax = np.nanpercentile(finite, [lower_percentile, upper_percentile])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        median = float(np.nanmedian(finite))
        std = float(np.nanstd(finite)) or 1.0
        return median - std, median + 5.0 * std
    return float(vmin), float(vmax)


def flatten_cutout_background(
    cutout: np.ndarray,
    *,
    filter_size: int = CUTOUT_MEDIAN_FILTER_SIZE,
) -> np.ndarray:
    """Subtract a 61-pixel median-filter background while preserving NaNs."""

    try:
        from scipy.ndimage import median_filter
    except ImportError as exc:  # pragma: no cover - clear failure on lean envs
        raise RuntimeError("cutout background flattening requires scipy") from exc

    array = np.asarray(cutout, dtype=float)
    finite = np.isfinite(array)
    if not np.any(finite):
        return array.copy()
    fill = float(np.nanmedian(array[finite]))
    filled = np.where(finite, array, fill)
    background = median_filter(filled, size=filter_size, mode="nearest")
    residual = array - background
    residual[~finite] = np.nan
    return residual


def flatten_cutouts(
    cutouts: Sequence[np.ndarray],
    *,
    filter_size: int = CUTOUT_MEDIAN_FILTER_SIZE,
) -> list[np.ndarray]:
    """Background-flatten every cutout with the same median-filter width."""

    return [
        flatten_cutout_background(cutout, filter_size=filter_size)
        for cutout in cutouts
    ]


def shared_cutout_limits(cutouts: Sequence[np.ndarray]) -> tuple[float, float]:
    """Return one median/MAD stretch for a complete set of display cutouts."""

    if len(cutouts) == 0:
        return 0.0, 1.0
    values = np.concatenate([np.asarray(cutout, dtype=float).ravel() for cutout in cutouts])
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0, 1.0

    center = float(np.nanmedian(values))
    mad = float(np.nanmedian(np.abs(values - center)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.nanstd(values))
    if not np.isfinite(sigma) or sigma <= 0:
        return robust_display_limits(values)
    return (
        center - CUTOUT_SHARED_LOWER_SIGMA * sigma,
        center + CUTOUT_SHARED_UPPER_SIGMA * sigma,
    )


def cutout_imshow_kwargs(vmin: float, vmax: float) -> dict[str, Any]:
    """Return the common GOTTA grayscale cutout display kwargs."""

    return {
        "origin": "lower",
        "cmap": CUTOUT_CMAP,
        "interpolation": CUTOUT_INTERPOLATION,
        "vmin": vmin,
        "vmax": vmax,
    }


def clear_cutout_axis(ax) -> None:
    """Remove ticks and spines from an image-cutout axis."""

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_gapped_crosshair(
    ax,
    x: float,
    y: float,
    color: str,
    *,
    gap: float = 11.0,
    length: float = 17.0,
    linewidth: float = 2.4,
):
    """Draw GOTTA's four-segment crosshair and return the line artists."""

    lines = []
    for x0, x1, y0, y1 in (
        (x - gap - length, x - gap, y, y),
        (x + gap, x + gap + length, y, y),
        (x, x, y - gap - length, y - gap),
        (x, x, y + gap, y + gap + length),
    ):
        (line,) = ax.plot(
            [x0, x1],
            [y0, y1],
            color=color,
            linewidth=linewidth,
            solid_capstyle="butt",
        )
        lines.append(line)
    return lines


def add_cutout_label(
    ax,
    x: float,
    y: float,
    text: str,
    *,
    fontsize: float = 10.0,
    ha: str = "left",
    va: str = "top",
):
    """Draw GOTTA's white-on-translucent-black cutout label."""

    return ax.text(
        x,
        y,
        text,
        ha=ha,
        va=va,
        fontsize=fontsize,
        **CUTOUT_LABEL_STYLE,
    )


# ---------------------------------------------------------------------------
# 4. Workflow / architecture diagrams
# ---------------------------------------------------------------------------

WORKFLOW_CANVAS = (2972, 3246)
WORKFLOW_DPI = (504, 504)
WORKFLOW_FONT_FAMILY = "Arial"
WORKFLOW_FONT_WEIGHT = "bold"
WORKFLOW_MAIN_FONT_SIZE = 41
WORKFLOW_TEXT_SPACING = 8
WORKFLOW_FONT_CANDIDATES = (
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/Library/Fonts/Arial.ttf",
    "/System/Library/Fonts/Helvetica.ttc",
)

WORKFLOW_COLORS: dict[str, tuple[int, int, int, int]] = {
    "ink": (31, 41, 55, 255),
    "arrow": (31, 41, 55, 255),
    "blue_fill": (221, 235, 255, 255),
    "blue_line": (37, 99, 255, 255),
    "green_fill": (233, 255, 244, 255),
    "green_line": (0, 155, 106, 255),
    "orange_fill": (255, 244, 214, 255),
    "orange_line": (242, 140, 0, 255),
    "purple_fill": (241, 228, 255, 255),
    "purple_line": (125, 44, 255, 255),
    "red_fill": (255, 242, 242, 255),
    "red_line": (255, 30, 30, 255),
    "gray_fill": (247, 247, 247, 255),
    "gray_line": (107, 114, 128, 255),
    "white": (255, 255, 255, 255),
}

WORKFLOW_NODE_STYLE: dict[str, Any] = {
    "outline_width": 4,
    "rounded_radius": 28,
    "document_outline_width": 3,
    "document_fold_width": 2,
    "text_spacing": WORKFLOW_TEXT_SPACING,
    "text_align": "center",
    "process_shape": "round",
    "io_shape": "parallelogram",
    "decision_shape": "diamond",
    "document_shape": "rect",
}
WORKFLOW_ARROW_STYLE: dict[str, float | int] = {
    "line_width": 4,
    "head_length": 18.0,
    "head_half_width": 8.0,
    "line_end_inset": 14.0,
}


def workflow_node_colors(kind: str):
    """Return the pastel fill and saturated outline for a workflow node kind."""

    if kind not in {"blue", "green", "orange", "purple", "red", "gray"}:
        raise ValueError(f"unknown workflow color kind: {kind!r}")
    return WORKFLOW_COLORS[f"{kind}_fill"], WORKFLOW_COLORS[f"{kind}_line"]


def find_workflow_font() -> Path | None:
    """Return the first available Arial/Helvetica reference font path."""

    for candidate in WORKFLOW_FONT_CANDIDATES:
        path = Path(candidate)
        if path.exists():
            return path
    return None


def workflow_arrowhead_points(
    tip: tuple[float, float],
    previous: tuple[float, float],
    *,
    head_length: float = float(WORKFLOW_ARROW_STYLE["head_length"]),
    head_half_width: float = float(WORKFLOW_ARROW_STYLE["head_half_width"]),
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    """Return the three polygon vertices for the GOTTA workflow arrowhead."""

    x2, y2 = tip
    x1, y1 = previous
    dx, dy = x2 - x1, y2 - y1
    segment_length = math.hypot(dx, dy)
    if segment_length == 0:
        raise ValueError("arrow tip and previous point must differ")
    ux, uy = dx / segment_length, dy / segment_length
    px, py = -uy, ux
    base_x = x2 - ux * head_length
    base_y = y2 - uy * head_length
    return (
        (base_x + px * head_half_width, base_y + py * head_half_width),
        (base_x - px * head_half_width, base_y - py * head_half_width),
        (x2, y2),
    )


def _self_test() -> None:
    """Run dependency-light checks without writing files or loading live data."""

    apply_statistics_style()
    assert plt.rcParams["font.family"][0] == "Times New Roman"
    assert float(plt.rcParams["font.size"]) == 30.0

    wrapped = longitude_to_mollweide_rad([0.0, 180.0, 359.0])
    assert wrapped.shape == (3,)
    assert np.all(wrapped >= -np.pi) and np.all(wrapped <= np.pi)

    cutout = np.arange(81, dtype=float).reshape(9, 9)
    flattened = flatten_cutout_background(cutout)
    assert flattened.shape == cutout.shape
    vmin, vmax = shared_cutout_limits([flattened, flattened + 1.0])
    assert np.isfinite(vmin) and np.isfinite(vmax) and vmin < vmax

    fill, outline = workflow_node_colors("blue")
    assert fill == WORKFLOW_COLORS["blue_fill"]
    assert outline == WORKFLOW_COLORS["blue_line"]
    points = workflow_arrowhead_points((10.0, 0.0), (0.0, 0.0))
    assert len(points) == 3 and points[-1] == (10.0, 0.0)


if __name__ == "__main__":
    _self_test()
    print("figure_styles self-test: ok")
