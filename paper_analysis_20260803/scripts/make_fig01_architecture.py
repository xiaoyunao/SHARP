#!/usr/bin/env python3
"""Generate Fig. 1, the SHARP production architecture diagram.

The script reads only the frozen paper configuration and checks a fixed list of
repository code paths.  It never scans or reads live survey data.  All plotted
text is English, the PDF consists of vector Matplotlib artists, and the PNG is
saved at 300 dpi.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon
from matplotlib.path import Path as MplPath
from PIL import Image

from figure_styles import (
    WORKFLOW_ARROW_STYLE,
    WORKFLOW_COLORS,
    WORKFLOW_FONT_FAMILY,
    WORKFLOW_FONT_WEIGHT,
    WORKFLOW_MAIN_FONT_SIZE,
    WORKFLOW_NODE_STYLE,
    find_workflow_font,
    save_pdf_png,
    workflow_node_colors,
)


SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parent
REPOSITORY_ROOT = ANALYSIS_ROOT.parent
DEFAULT_SNAPSHOT = ANALYSIS_ROOT / "config" / "snapshot.json"
DEFAULT_OUTPUT_STEM = ANALYSIS_ROOT / "figures" / "fig01_system_architecture"

# Fixed paths only: validating these does not inspect any live data directory.
EXPECTED_CODE_PATHS = {
    "scheduler": REPOSITORY_ROOT / "survey" / "scheduler.py",
    "scheduler daily runner": REPOSITORY_ROOT / "survey" / "run_daily.py",
    "follow-up state machine": REPOSITORY_ROOT / "survey" / "followup.py",
    "follow-up plan insertion": REPOSITORY_ROOT / "survey" / "apply_followup.py",
    "known-object matcher": REPOSITORY_ROOT / "known_asteroid" / "match_single_night.py",
    "known-object ADES exporter": REPOSITORY_ROOT / "known_asteroid" / "export_ades.py",
    "daily orchestration": REPOSITORY_ROOT / "heliolincrr" / "run_daily_pipeline.sh",
    "Gaia masking": REPOSITORY_ROOT / "heliolincrr" / "mask_gaia.py",
    "tracklet construction": REPOSITORY_ROOT / "heliolincrr" / "make_tracklet_linreproj.py",
    "shared-endpoint linking": REPOSITORY_ROOT / "heliolincrr" / "run_linear_links_from_tracklets.py",
    "short-arc consistency": REPOSITORY_ROOT / "heliolincrr" / "orbit_confirm_links.py",
    "known-object subtraction": REPOSITORY_ROOT / "heliolincrr" / "remask_unknown_with_known.py",
    "review packaging": REPOSITORY_ROOT / "heliolincrr" / "package_unknown_review.py",
    "reviewed submission": REPOSITORY_ROOT / "heliolincrr" / "submit_reviewed_unknown.py",
}


@dataclass(frozen=True)
class Node:
    """A diagram node in the figure's 14 by 16 coordinate system."""

    key: str
    label: str
    x: float
    y: float
    width: float
    height: float
    category: str
    shape: str
    fontsize: float = 11.5

    def anchor(self, side: str) -> tuple[float, float]:
        if side == "top":
            return self.x, self.y + self.height / 2.0
        if side == "bottom":
            return self.x, self.y - self.height / 2.0
        if side == "left":
            return self.x - self.width / 2.0, self.y
        if side == "right":
            return self.x + self.width / 2.0, self.y
        raise ValueError(f"unknown anchor side: {side!r}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=DEFAULT_SNAPSHOT,
        help="Frozen paper configuration JSON.",
    )
    parser.add_argument(
        "--output-stem",
        type=Path,
        default=DEFAULT_OUTPUT_STEM,
        help="Output path without extension.",
    )
    return parser.parse_args()


def load_snapshot(path: Path) -> dict:
    """Load and minimally validate the frozen paper configuration."""

    if not path.is_file():
        raise FileNotFoundError(f"frozen snapshot configuration not found: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    required = (
        "snapshot_label",
        "observation_start",
        "observation_end",
        "production_algorithm_commit",
    )
    missing = [key for key in required if not payload.get(key)]
    if missing:
        raise ValueError(f"snapshot configuration lacks required fields: {missing}")
    return payload


def validate_code_paths() -> dict[str, Path]:
    """Check only the fixed repository modules represented in the diagram."""

    missing = {label: path for label, path in EXPECTED_CODE_PATHS.items() if not path.is_file()}
    if missing:
        details = "\n".join(f"- {label}: {path}" for label, path in missing.items())
        raise FileNotFoundError(f"architecture code-path validation failed:\n{details}")
    return dict(EXPECTED_CODE_PATHS)


def font(size: float) -> FontProperties:
    """Return the exact Arial Bold workflow font where it is installed."""

    font_path = find_workflow_font()
    if font_path is not None:
        return FontProperties(fname=str(font_path), size=size)
    return FontProperties(
        family=WORKFLOW_FONT_FAMILY,
        weight=WORKFLOW_FONT_WEIGHT,
        size=size,
    )


def node_colors(category: str):
    color_kind = {
        "automatic": "green",
        "human": "orange",
        "external": "gray",
        "persistent": "purple",
    }.get(category)
    if color_kind is None:
        raise ValueError(f"unknown architecture category: {category!r}")
    fill, outline = workflow_node_colors(color_kind)
    return matplotlib_color(fill), matplotlib_color(outline)


def matplotlib_color(rgba: tuple[int, int, int, int]) -> tuple[float, float, float, float]:
    """Convert the shared Pillow-style 0--255 RGBA palette for Matplotlib."""

    return tuple(component / 255.0 for component in rgba)


def workflow_color(name: str) -> tuple[float, float, float, float]:
    return matplotlib_color(WORKFLOW_COLORS[name])


def draw_node(ax, node: Node) -> None:
    """Draw a rounded process, diamond, parallelogram, or document node."""

    fill, outline = node_colors(node.category)
    x0 = node.x - node.width / 2.0
    x1 = node.x + node.width / 2.0
    y0 = node.y - node.height / 2.0
    y1 = node.y + node.height / 2.0
    linewidth = float(WORKFLOW_NODE_STYLE["outline_width"]) / 2.0

    if node.shape == "round":
        patch = FancyBboxPatch(
            (x0, y0),
            node.width,
            node.height,
            boxstyle="round,pad=0.025,rounding_size=0.14",
            facecolor=fill,
            edgecolor=outline,
            linewidth=linewidth,
            joinstyle="round",
            zorder=3,
        )
        ax.add_patch(patch)
    elif node.shape == "parallelogram":
        skew = min(0.28, node.width * 0.08)
        patch = Polygon(
            [(x0 + skew, y1), (x1, y1), (x1 - skew, y0), (x0, y0)],
            closed=True,
            facecolor=fill,
            edgecolor=outline,
            linewidth=linewidth,
            joinstyle="round",
            zorder=3,
        )
        ax.add_patch(patch)
    elif node.shape == "diamond":
        patch = Polygon(
            [(node.x, y1), (x1, node.y), (node.x, y0), (x0, node.y)],
            closed=True,
            facecolor=fill,
            edgecolor=outline,
            linewidth=linewidth,
            joinstyle="round",
            zorder=3,
        )
        ax.add_patch(patch)
    elif node.shape == "document":
        fold = min(0.24, node.height * 0.23)
        points = [(x0, y0), (x1, y0), (x1, y1 - fold), (x1 - fold, y1), (x0, y1)]
        patch = Polygon(
            points,
            closed=True,
            facecolor=fill,
            edgecolor=outline,
            linewidth=linewidth,
            joinstyle="round",
            zorder=3,
        )
        ax.add_patch(patch)
        ax.plot(
            [x1 - fold, x1 - fold, x1],
            [y1, y1 - fold, y1 - fold],
            color=outline,
            linewidth=max(1.0, linewidth * 0.7),
            zorder=4,
        )
    else:
        raise ValueError(f"unknown workflow node shape: {node.shape!r}")

    ax.text(
        node.x,
        node.y,
        node.label,
        ha="center",
        va="center",
        color=workflow_color("ink"),
        fontproperties=font(node.fontsize),
        linespacing=1.12,
        zorder=5,
    )


def draw_arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    dashed: bool = False,
    dotted: bool = False,
    connectionstyle: str = "arc3,rad=0",
    zorder: int = 2,
) -> FancyArrowPatch:
    """Draw a workflow arrow using the shared ink and arrow constants."""

    if dashed and dotted:
        raise ValueError("an arrow cannot be both dashed and dotted")
    linestyle: str | tuple = "solid"
    if dashed:
        linestyle = (0, (6, 4))
    elif dotted:
        linestyle = (0, (1.5, 2.5))

    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=float(WORKFLOW_ARROW_STYLE["head_length"]),
        linewidth=float(WORKFLOW_ARROW_STYLE["line_width"]) / 2.0,
        linestyle=linestyle,
        color=workflow_color("arrow"),
        connectionstyle=connectionstyle,
        shrinkA=0,
        shrinkB=0,
        capstyle="butt",
        joinstyle="round",
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def draw_polyline_arrow(
    ax,
    points: Iterable[tuple[float, float]],
    *,
    dashed: bool = False,
    dotted: bool = False,
) -> FancyArrowPatch:
    """Draw an orthogonal multi-segment arrow, used for the feedback loop."""

    vertices = list(points)
    if len(vertices) < 2:
        raise ValueError("a polyline arrow requires at least two points")
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(vertices) - 1)
    path = MplPath(vertices, codes)
    linestyle: str | tuple = "solid"
    if dashed:
        linestyle = (0, (6, 4))
    elif dotted:
        linestyle = (0, (1.5, 2.5))
    patch = FancyArrowPatch(
        path=path,
        arrowstyle="-|>",
        mutation_scale=float(WORKFLOW_ARROW_STYLE["head_length"]),
        linewidth=float(WORKFLOW_ARROW_STYLE["line_width"]) / 2.0,
        linestyle=linestyle,
        color=workflow_color("arrow"),
        capstyle="butt",
        joinstyle="round",
        zorder=2,
    )
    ax.add_patch(patch)
    return patch


def build_nodes() -> dict[str, Node]:
    """Define the complete production and operational architecture."""

    items = [
        Node(
            "scheduler_state",
            "Plans + exposure history\n+ follow-up state",
            1.75,
            14.25,
            2.8,
            0.95,
            "persistent",
            "document",
            10.5,
        ),
        Node(
            "scheduler",
            "Dynamic ecliptic-plane scheduler",
            7.0,
            14.25,
            4.2,
            0.9,
            "automatic",
            "round",
            12.0,
        ),
        Node(
            "observing",
            "Telescope observing\n30 s survey exposures",
            7.0,
            12.95,
            4.0,
            0.95,
            "external",
            "parallelogram",
            11.5,
        ),
        Node(
            "calibration",
            "External L1/L2 calibration\nWCS + source extraction + photometry",
            7.0,
            11.55,
            5.15,
            1.05,
            "external",
            "parallelogram",
            11.0,
        ),
        Node(
            "nightly_products",
            "Nightly L1 images + L2 catalogs\nwith FITS headers",
            7.0,
            10.15,
            4.7,
            1.0,
            "persistent",
            "document",
            11.5,
        ),
        Node(
            "known_query",
            "Ephemeris query +\nWCS detector footprint",
            3.45,
            8.65,
            4.45,
            0.95,
            "automatic",
            "round",
        ),
        Node(
            "known_match",
            "1 arcsec source association\n+ photometric consistency",
            3.45,
            7.30,
            4.45,
            0.95,
            "automatic",
            "round",
        ),
        Node(
            "known_store",
            "Matched detections\n+ nightly status",
            3.45,
            5.95,
            3.75,
            0.95,
            "persistent",
            "document",
        ),
        Node(
            "known_ades",
            "Generate ADES PSV\n+ duplicate guard",
            3.45,
            4.60,
            3.75,
            0.95,
            "automatic",
            "round",
        ),
        Node(
            "gaia_mask",
            "Quality prefilter\n+ Gaia source mask",
            9.45,
            8.65,
            4.55,
            0.95,
            "automatic",
            "round",
        ),
        Node(
            "tracklets_links",
            "Two-point tracklets\n+ shared-endpoint links",
            9.45,
            7.30,
            4.55,
            0.95,
            "automatic",
            "round",
        ),
        Node(
            "orbit_subtract",
            "Short-arc consistency\n+ 1.5 arcsec known subtraction",
            9.45,
            5.95,
            4.75,
            0.95,
            "automatic",
            "round",
            11.0,
        ),
        Node(
            "review_package",
            "Review package\ncatalog + cutout GIFs",
            9.45,
            4.60,
            3.75,
            0.95,
            "persistent",
            "document",
        ),
        Node(
            "review_state",
            "trkSub history\nreview signatures\nsubmission ledger",
            12.42,
            4.60,
            2.05,
            1.05,
            "persistent",
            "document",
            8.5,
        ),
        Node(
            "human_review",
            "Human candidate\nvetting",
            9.45,
            3.25,
            3.6,
            1.15,
            "human",
            "diamond",
            11.0,
        ),
        Node(
            "unknown_ades",
            "Generate unknown-object ADES PSV\n+ submission guard",
            8.65,
            1.95,
            4.05,
            0.95,
            "automatic",
            "round",
            10.0,
        ),
        Node(
            "followup_state",
            "Conditional follow-up\nqueue + persistent\nstate",
            12.00,
            1.95,
            2.10,
            0.95,
            "persistent",
            "document",
            8.5,
        ),
        Node(
            "mpc",
            "Minor Planet Center\nsubmission + server replies",
            6.25,
            0.85,
            4.65,
            0.95,
            "external",
            "parallelogram",
            11.0,
        ),
    ]
    return {item.key: item for item in items}


def draw_legend(ax) -> None:
    entries = [
        ("automatic", "Automated SHARP step", "round"),
        ("human", "Human decision", "diamond"),
        ("external", "External system/service", "parallelogram"),
        ("persistent", "Persistent product/state", "document"),
    ]
    x_positions = (0.95, 4.45, 7.75, 11.0)
    y = 0.08
    for x, (category, label, shape) in zip(x_positions, entries):
        sample = Node(
            f"legend_{category}",
            "",
            x,
            y + 0.11,
            0.36,
            0.22,
            category,
            shape,
            1,
        )
        draw_node(ax, sample)
        ax.text(
            x + 0.31,
            y + 0.11,
            label,
            ha="left",
            va="center",
            color=workflow_color("ink"),
            fontproperties=font(8.7),
            zorder=5,
        )


def make_figure(snapshot: dict, validated_paths: dict[str, Path]):
    """Build the vector architecture figure in memory."""

    plt.rcParams.update(
        {
            "font.family": WORKFLOW_FONT_FAMILY,
            "font.weight": WORKFLOW_FONT_WEIGHT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(figsize=(14, 16), facecolor="white")
    ax.set_xlim(0, 14)
    ax.set_ylim(-0.05, 16)
    ax.set_aspect("equal")
    ax.axis("off")

    commit = str(snapshot["production_algorithm_commit"])
    subtitle = (
        f"Frozen configuration {snapshot['snapshot_label']}  |  "
        f"observations {snapshot['observation_start']} to {snapshot['observation_end']}  |  "
        f"production algorithm {commit[:7]}  |  modules validated {len(validated_paths)}/{len(EXPECTED_CODE_PATHS)}"
    )
    ax.text(
        7.0,
        15.68,
        "SHARP Production Architecture and Operational Feedback Loop",
        ha="center",
        va="center",
        color=workflow_color("ink"),
        fontproperties=font(20),
    )
    ax.text(
        7.0,
        15.26,
        subtitle,
        ha="center",
        va="center",
        color=workflow_color("gray_line"),
        fontproperties=font(8.8),
    )

    nodes = build_nodes()

    # Main acquisition chain.
    draw_arrow(ax, nodes["scheduler"].anchor("bottom"), nodes["observing"].anchor("top"))
    draw_arrow(ax, nodes["observing"].anchor("bottom"), nodes["calibration"].anchor("top"))
    draw_arrow(ax, nodes["calibration"].anchor("bottom"), nodes["nightly_products"].anchor("top"))
    draw_arrow(
        ax,
        nodes["scheduler_state"].anchor("right"),
        nodes["scheduler"].anchor("left"),
        dotted=True,
    )

    # L1/L2 dispatch to the two automated branches.
    draw_arrow(
        ax,
        (6.75, nodes["nightly_products"].anchor("bottom")[1]),
        nodes["known_query"].anchor("top"),
        connectionstyle="arc3,rad=0.12",
    )
    draw_arrow(
        ax,
        (7.25, nodes["nightly_products"].anchor("bottom")[1]),
        nodes["gaia_mask"].anchor("top"),
        connectionstyle="arc3,rad=-0.12",
    )

    # Known-object branch.
    for source, target in (
        ("known_query", "known_match"),
        ("known_match", "known_store"),
        ("known_store", "known_ades"),
    ):
        draw_arrow(ax, nodes[source].anchor("bottom"), nodes[target].anchor("top"))

    # Unknown-object branch and persistent review state.
    for source, target in (
        ("gaia_mask", "tracklets_links"),
        ("tracklets_links", "orbit_subtract"),
        ("orbit_subtract", "review_package"),
        ("review_package", "human_review"),
        ("human_review", "unknown_ades"),
    ):
        draw_arrow(ax, nodes[source].anchor("bottom"), nodes[target].anchor("top"))
    draw_arrow(
        ax,
        nodes["review_state"].anchor("left"),
        nodes["review_package"].anchor("right"),
        dotted=True,
    )

    # Both report products reach the external MPC service.
    draw_arrow(
        ax,
        nodes["known_ades"].anchor("bottom"),
        nodes["mpc"].anchor("left"),
        connectionstyle="arc3,rad=-0.12",
    )
    draw_arrow(
        ax,
        nodes["unknown_ades"].anchor("bottom"),
        nodes["mpc"].anchor("right"),
        connectionstyle="arc3,rad=0.10",
    )

    # Conditional follow-up is deliberately dashed and loops back to scheduling.
    draw_arrow(
        ax,
        nodes["human_review"].anchor("right"),
        nodes["followup_state"].anchor("top"),
        dashed=True,
        connectionstyle="arc3,rad=-0.18",
    )
    draw_polyline_arrow(
        ax,
        [
            nodes["followup_state"].anchor("right"),
            (13.78, 1.95),
            (13.78, 14.25),
            nodes["scheduler"].anchor("right"),
        ],
        dashed=True,
    )

    # Draw nodes after connectors so boxes mask line endpoints cleanly.
    for node in nodes.values():
        draw_node(ax, node)

    # Branch headers and feedback-loop status annotation.
    for x, text in ((3.45, "KNOWN-OBJECT BRANCH"), (9.45, "UNKNOWN MOVING-OBJECT BRANCH")):
        ax.text(
            x,
            9.40,
            text,
            ha="center",
            va="center",
            color=workflow_color("blue_line"),
            fontproperties=font(11.5),
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
            zorder=5,
        )
        ax.plot(
            [x - 2.05, x + 2.05],
            [9.20, 9.20],
            color=workflow_color("blue_line"),
            linewidth=1.5,
        )

    ax.text(
        13.56,
        9.85,
        "Conditional follow-up request\n(infrastructure; no completed on-sky cycle in this snapshot)",
        rotation=90,
        ha="center",
        va="center",
        color=workflow_color("gray_line"),
        fontproperties=font(8.5),
        linespacing=1.15,
    )
    ax.text(
        11.45,
        2.92,
        "eligible",
        ha="center",
        va="center",
        color=workflow_color("ink"),
        fontproperties=font(8.5),
    )

    draw_legend(ax)
    return fig


def verify_outputs(png_path: Path, pdf_path: Path) -> dict[str, object]:
    """Perform basic size, DPI, signature, and non-empty-file checks."""

    for path in (png_path, pdf_path):
        if not path.is_file() or path.stat().st_size < 10_000:
            raise RuntimeError(f"figure output is missing or unexpectedly small: {path}")

    with Image.open(png_path) as image:
        width, height = image.size
        dpi = image.info.get("dpi", (0.0, 0.0))
    if width < 3000 or height < 3500 or height <= width:
        raise RuntimeError(f"unexpected architecture PNG dimensions: {width}x{height}")
    if len(dpi) != 2 or any(abs(float(value) - 300.0) > 2.0 for value in dpi):
        raise RuntimeError(f"architecture PNG lacks 300-dpi metadata: {dpi}")

    with pdf_path.open("rb") as handle:
        signature = handle.read(5)
    if signature != b"%PDF-":
        raise RuntimeError(f"invalid PDF signature for {pdf_path}")

    return {
        "png_width_px": width,
        "png_height_px": height,
        "png_dpi": tuple(round(float(value), 2) for value in dpi),
        "png_bytes": png_path.stat().st_size,
        "pdf_bytes": pdf_path.stat().st_size,
    }


def main() -> None:
    args = parse_args()
    snapshot = load_snapshot(args.snapshot.resolve())
    validated_paths = validate_code_paths()
    figure = make_figure(snapshot, validated_paths)
    png_path, pdf_path = save_pdf_png(
        figure,
        args.output_stem.resolve(),
        png_dpi=300,
        tight=True,
        pad_inches=0.05,
    )
    validation = verify_outputs(png_path, pdf_path)
    print(f"snapshot: {args.snapshot.resolve()}")
    print(f"validated code paths: {len(validated_paths)}/{len(EXPECTED_CODE_PATHS)}")
    print(f"PNG: {png_path} ({validation['png_width_px']}x{validation['png_height_px']}, {validation['png_dpi']} dpi)")
    print(f"PDF: {pdf_path} ({validation['pdf_bytes']} bytes)")


if __name__ == "__main__":
    main()
