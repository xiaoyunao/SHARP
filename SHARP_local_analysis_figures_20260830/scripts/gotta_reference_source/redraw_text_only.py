#!/usr/bin/env python3
"""Standalone workflow revision with original structure and uniform node text size.

This does not modify the manuscript figure. It generates a clean candidate in
the revision directory, using the approved layout as the coordinate template.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from textwrap import wrap

from PIL import Image, ImageDraw, ImageFont


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = REPO_ROOT / "workflow_figure_revision_20260702"
DPI = (504, 504)
CANVAS = (2972, 3246)
MAIN_FONT_SIZE = 41

COLORS = {
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

FONT_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/Library/Fonts/Arial.ttf",
    "/System/Library/Fonts/Helvetica.ttc",
]


@dataclass
class Node:
    key: str
    text: str
    cx: int
    cy: int
    kind: str
    shape: str
    w: int
    h: int
    skew: int = 0
    font_size: int = MAIN_FONT_SIZE
    wrap_width: int | None = None

    @property
    def bbox(self) -> tuple[int, int, int, int]:
        return (self.cx - self.w // 2, self.cy - self.h // 2, self.cx + self.w // 2, self.cy + self.h // 2)

    @property
    def top(self) -> tuple[int, int]:
        return (self.cx, self.cy - self.h // 2)

    @property
    def bottom(self) -> tuple[int, int]:
        return (self.cx, self.cy + self.h // 2)

    @property
    def left(self) -> tuple[int, int]:
        return (self.cx - self.w // 2, self.cy)

    @property
    def right(self) -> tuple[int, int]:
        return (self.cx + self.w // 2, self.cy)


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for candidate in FONT_CANDIDATES:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size)
    return ImageFont.load_default()


def wrap_text(text: str, width: int | None) -> str:
    if width is None:
        return text
    lines: list[str] = []
    for raw in text.splitlines():
        lines.extend(wrap(raw, width=width, break_long_words=False) or [""])
    return "\n".join(lines)


def text_extent(draw: ImageDraw.ImageDraw, text: str, size: int) -> tuple[int, int]:
    font = load_font(size)
    bbox = draw.multiline_textbbox((0, 0), text, font=font, spacing=8, align="center")
    return bbox[2] - bbox[0], bbox[3] - bbox[1]


def fit_node_to_text(draw: ImageDraw.ImageDraw, node: Node) -> None:
    text = wrap_text(node.text, node.wrap_width)
    tw, th = text_extent(draw, text, node.font_size)
    extra_w = 90 if node.shape == "parallelogram" else 0
    if node.shape == "diamond":
        extra_w = 180
        extra_h = 70
    else:
        extra_h = 32
    node.w = max(node.w, tw + extra_w + 54)
    node.h = max(node.h, th + extra_h)
    # Keep the large bottom parallelograms inside the canvas after automatic growth.
    x1, _, x2, _ = node.bbox
    if x1 < 25:
        node.cx += 25 - x1
    if x2 > CANVAS[0] - 25:
        node.cx -= x2 - (CANVAS[0] - 25)


def style(kind: str) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    return COLORS[f"{kind}_fill"], COLORS[f"{kind}_line"]


def draw_centered_text(draw: ImageDraw.ImageDraw, node: Node, text_box: tuple[int, int, int, int]) -> None:
    text = wrap_text(node.text, node.wrap_width)
    font = load_font(node.font_size)
    bbox = draw.multiline_textbbox((0, 0), text, font=font, spacing=8, align="center")
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    x1, y1, x2, y2 = text_box
    draw.multiline_text(
        (x1 + (x2 - x1 - tw) / 2, y1 + (y2 - y1 - th) / 2 - 1),
        text,
        font=font,
        fill=COLORS["ink"],
        spacing=8,
        align="center",
    )


def draw_node(draw: ImageDraw.ImageDraw, node: Node) -> None:
    x1, y1, x2, y2 = node.bbox
    fill, line = style(node.kind)
    lw = 4
    if node.shape == "round":
        draw.rounded_rectangle((x1, y1, x2, y2), radius=28, fill=fill, outline=line, width=lw)
        text_box = (x1 + 12, y1 + 8, x2 - 12, y2 - 8)
    elif node.shape == "rect":
        draw.rectangle((x1, y1, x2, y2), fill=fill, outline=line, width=3)
        draw.line((x2 - 24, y1, x2, y1 + 24), fill=line, width=2)
        text_box = (x1 + 12, y1 + 8, x2 - 12, y2 - 8)
    elif node.shape == "parallelogram":
        s = node.skew
        pts = [(x1 + s, y1), (x2, y1), (x2 - s, y2), (x1, y2)]
        draw.polygon(pts, fill=fill)
        draw.line(pts + [pts[0]], fill=line, width=lw, joint="curve")
        text_box = (x1 + s // 2 + 10, y1 + 8, x2 - s // 2 - 10, y2 - 8)
    elif node.shape == "diamond":
        pts = [((x1 + x2) // 2, y1), (x2, (y1 + y2) // 2), ((x1 + x2) // 2, y2), (x1, (y1 + y2) // 2)]
        draw.polygon(pts, fill=fill)
        draw.line(pts + [pts[0]], fill=line, width=lw)
        text_box = (x1 + 95, y1 + 32, x2 - 95, y2 - 32)
    else:
        raise ValueError(node.shape)
    draw_centered_text(draw, node, text_box)


def arrowhead(draw: ImageDraw.ImageDraw, tip: tuple[int, int], prev: tuple[int, int], color: tuple[int, int, int, int]) -> None:
    x2, y2 = tip
    x1, y1 = prev
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy)
    if length == 0:
        return
    ux, uy = dx / length, dy / length
    px, py = -uy, ux
    size = 18
    bx, by = x2 - ux * size, y2 - uy * size
    draw.polygon([(bx + px * 8, by + py * 8), (bx - px * 8, by - py * 8), (x2, y2)], fill=color)


def draw_arrow(
    draw: ImageDraw.ImageDraw,
    pts: list[tuple[int, int]],
    color: tuple[int, int, int, int] = COLORS["arrow"],
    width: int = 4,
) -> None:
    if len(pts) < 2:
        return
    end = pts[-1]
    prev = pts[-2]
    dx, dy = end[0] - prev[0], end[1] - prev[1]
    length = math.hypot(dx, dy)
    if length:
        ux, uy = dx / length, dy / length
        line_end = (round(end[0] - ux * 14), round(end[1] - uy * 14))
        draw.line(pts[:-1] + [line_end], fill=color, width=width)
    arrowhead(draw, end, prev, color)


def nodes(draw: ImageDraw.ImageDraw) -> dict[str, Node]:
    items = {
        "input": Node("input", "Nightly calibrated products\n(source catalogs + FITS headers)", 1478, 138, "blue", "parallelogram", 1171, 184, 240),
        "note": Node(
            "note",
            "This paper: the same core logic is re-run in batch mode on the accumulated prototype data set for homogeneous statistics.",
            2190,
            392,
            "gray",
            "rect",
            850,
            150,
            font_size=MAIN_FONT_SIZE,
            wrap_width=42,
        ),
        "dispatch": Node("dispatch", "Parallel dispatcher\n(assign new exposures to workers)", 1133, 382, "green", "round", 588, 98),
        "wcs": Node("wcs", "Build WCS footprint\nand mid-exposure epoch", 1133, 584, "green", "round", 432, 98),
        "orbit": Node("orbit", "Orbit-element screening\nLowell/MPC elements -> predicted RA/Dec", 1133, 778, "green", "round", 692, 98),
        "filter": Node("filter", "Candidate filtering\nfield radius + magnitude + CMOS footprint", 1133, 973, "green", "round", 692, 99),
        "association": Node("association", "Sky-coordinate association\nnearest catalog source + Gaia mask", 1133, 1168, "orange", "round", 614, 98),
        "accepted": Node("accepted", "Accepted asteroid-source\nassociations", 1130, 1461, "purple", "parallelogram", 918, 183, 180),
        "reject": Node("reject", "Reject likely stationary\nGaia-source contaminants", 1936, 1461, "red", "round", 453, 99),
        "local": Node(
            "local",
            "Local ephemeris-table matching\nidentity, class, orbit, geometry, rates",
            1130,
            1699,
            "green",
            "round",
            623,
            98,
            font_size=MAIN_FONT_SIZE,
            wrap_width=44,
        ),
        "catalog": Node("catalog", "Catalog rows + cutout images", 1130, 1893, "purple", "parallelogram", 940, 99, 180),
        "group": Node("group", "Nightly grouping by asteroid", 1130, 2088, "green", "round", 482, 98),
        "store": Node("store", "Store all recovered detections\n(no minimum-count requirement)", 636, 2406, "purple", "parallelogram", 1183, 185, 240),
        "science": Node("science", "Science products\nstatistics + light-curve inputs", 636, 2696, "purple", "parallelogram", 1027, 184, 210),
        "decision": Node("decision", "Single-night detections\nfor one asteroid >= 3?", 1626, 2498, "orange", "diamond", 666, 175, font_size=MAIN_FONT_SIZE, wrap_width=28),
        "nompc": Node("nompc", "No MPC report for this object\nretain in local catalog", 1487, 2825, "gray", "round", 560, 125, font_size=MAIN_FONT_SIZE, wrap_width=27),
        "ades": Node("ades", "Generate ADES PSV astrometry\nfor MPC reporting", 2365, 2826, "blue", "parallelogram", 1123, 184, 230, font_size=MAIN_FONT_SIZE, wrap_width=36),
        "mpc": Node("mpc", "MPC submission\n(single-night report product)", 2365, 3106, "blue", "parallelogram", 1007, 184, 210),
    }
    for item in items.values():
        fit_node_to_text(draw, item)
    return items


def draw_labels(draw: ImageDraw.ImageDraw) -> None:
    label_font = load_font(MAIN_FONT_SIZE)
    for text, xy in [
        ("pass", (1160, 1276)),
        ("Gaia match\nwithin 1 arcsec", (2010, 1256)),
        ("no", (1408, 2650)),
        ("yes", (1930, 2604)),
    ]:
        draw.multiline_text(xy, text, font=label_font, fill=COLORS["ink"], spacing=8, align="left")


def draw_workflow() -> Image.Image:
    img = Image.new("RGBA", CANVAS, COLORS["white"])
    draw = ImageDraw.Draw(img)
    n = nodes(draw)

    # Arrows first, nodes second. This preserves clean node borders.
    draw_arrow(draw, [(1133, n["input"].bottom[1]), n["dispatch"].top])
    chain = ["dispatch", "wcs", "orbit", "filter", "association"]
    for a, b in zip(chain, chain[1:]):
        draw_arrow(draw, [n[a].bottom, n[b].top])
    # Keep the accepted-sample branch visually vertical despite the
    # parallelogram geometry and the slight center offsets in the original.
    draw_arrow(draw, [(1133, n["association"].bottom[1]), (1133, n["accepted"].top[1])])
    for a, b in [("accepted", "local"), ("local", "catalog"), ("catalog", "group")]:
        draw_arrow(draw, [n[a].bottom, n[b].top])

    # Note connector and rejection branch.
    note_x = 1705
    note_y = n["note"].top[1] + 42
    draw.line([(note_x, n["input"].bottom[1]), (note_x, note_y), (n["note"].left[0], note_y)], fill=(156, 163, 175, 255), width=3)
    draw_arrow(draw, [n["association"].right, (1943, n["association"].right[1]), (1943, n["reject"].top[1])], COLORS["red_line"], 4)

    # Bottom split.
    split_y = 2250
    draw_arrow(draw, [n["group"].bottom, (n["group"].bottom[0], split_y), (n["store"].top[0], split_y), n["store"].top])
    draw_arrow(draw, [n["group"].bottom, (n["group"].bottom[0], split_y), (n["decision"].top[0], split_y), n["decision"].top])
    draw_arrow(draw, [n["store"].bottom, n["science"].top])
    draw_arrow(draw, [(1517, 2556), (1517, n["nompc"].top[1])])
    draw_arrow(draw, [(1842, 2554), (2035, n["ades"].top[1])], COLORS["green_line"], 4)
    draw_arrow(draw, [n["ades"].bottom, n["mpc"].top])

    for key in [
        "input",
        "note",
        "dispatch",
        "wcs",
        "orbit",
        "filter",
        "association",
        "accepted",
        "reject",
        "local",
        "catalog",
        "group",
        "store",
        "science",
        "decision",
        "nompc",
        "ades",
        "mpc",
    ]:
        draw_node(draw, n[key])
    draw_labels(draw)
    return img


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate standalone workflow revision.")
    parser.add_argument("--outdir", default=str(OUTDIR))
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    img = draw_workflow()
    png = outdir / "known_object_processing_text_larger_cmos.png"
    pdf = outdir / "known_object_processing_text_larger_cmos.pdf"
    img.save(png, dpi=DPI)
    img.convert("RGB").save(pdf, resolution=DPI[0])
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
