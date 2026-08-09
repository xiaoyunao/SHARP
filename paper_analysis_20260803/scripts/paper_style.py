#!/usr/bin/env python3
"""Shared plotting style matched to ``gotta_asteroid_1`` paper figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


FOUR_PANEL_FIGSIZE = (18, 14)
TWO_PANEL_FIGSIZE = (18, 8)
SINGLE_PANEL_FIGSIZE = (10, 7.2)
MEDIAN_LEGEND_FONTSIZE = 21

# The same Tableau-derived colors used by gotta_asteroid_1.
BLUE = "#4c78a8"
ORANGE = "#f28e2b"
GREEN = "#59a14f"
RED = "#e15759"
CYAN = "#76b7b2"
GOLD = "#edc948"
PURPLE = "#b07aa1"
DARK_RED = "#d62728"
INK = "#2b2b2b"
MID_GREY = "#777777"
LIGHT_GREY = "#d9d9d9"
PALETTE = (BLUE, ORANGE, GREEN, RED, CYAN, GOLD, PURPLE)


def setup_style() -> None:
    """Apply the exact core rcParams used by gotta_asteroid_1."""

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.size": 30,
            "axes.linewidth": 1.4,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "legend.frameon": False,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def finish_axis(ax, *, tick_size: float = 22, grid_alpha: float = 0.16) -> None:
    ax.set_axisbelow(True)
    ax.grid(alpha=grid_alpha, linewidth=0.6)
    ax.tick_params(labelsize=tick_size)


def panel_label(ax, label: str, *, x: float = 0.015, y: float = 0.985) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=25,
        fontweight="bold",
        color=INK,
    )


def add_colorbar(ax, mappable, label: str) -> None:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3.5%", pad=0.06)
    colorbar = ax.figure.colorbar(mappable, cax=cax)
    colorbar.set_label(label, fontsize=20)
    colorbar.ax.tick_params(labelsize=17, width=0.9, length=4)
    colorbar.outline.set_linewidth(0.8)


def save_figure(fig: plt.Figure, output_stem: Path, *, tight: bool = True) -> tuple[Path, Path]:
    """Write 300-dpi PNG plus vector PDF, matching gotta's export contract."""

    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    png = output_stem.with_suffix(".png")
    pdf = output_stem.with_suffix(".pdf")
    kwargs = {"bbox_inches": "tight"} if tight else {}
    fig.savefig(png, dpi=300, **kwargs)
    fig.savefig(pdf, **kwargs)
    plt.close(fig)
    return png, pdf
