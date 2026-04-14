#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib import gridspec

for candidate in [
    Path(__file__).resolve().parents[1],
    Path("/pipeline/xiaoyunao"),
    Path("/Users/yunaoxiao/Desktop/smt_asteroid"),
]:
    if candidate.exists():
        text = str(candidate)
        if text not in sys.path:
            sys.path.insert(0, text)

from heliolincrr.make_tracklet_linreproj import group_exposures, read_l2_catalog, remove_static_sources


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a flow-style Gaia-mask + static-source subtraction figure.")
    parser.add_argument("night", help="Night in YYYYMMDD")
    parser.add_argument("--processed-root", default="/processed1")
    parser.add_argument("--root-out", default="/pipeline/xiaoyunao/data/heliolincrr")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--group-index", type=int, default=None, help="Optional fixed 3-exposure group index")
    parser.add_argument("--r-static", type=float, default=2.0)
    parser.add_argument("--min-repeat", type=int, default=2)
    return parser.parse_args()


def choose_best_group(mask_dir: Path, group_index: int | None, r_static: float, min_repeat: int):
    groups = group_exposures(sorted(str(p) for p in mask_dir.glob("*.fits.gz")), max_sep_deg=0.5)
    candidates: list[tuple[float, int, list[str], list[int], list[int]]] = []
    for gi, group in enumerate(groups):
        if len(group) != 3:
            continue
        files = [Path(item[0]).name for item in group]
        tabs = [read_l2_catalog(str(mask_dir / name)) for name in files]
        before = [len(t) for t in tabs]
        after_tabs = remove_static_sources(tabs, r_arcsec=r_static, min_repeat=min_repeat)
        after = [len(t) for t in after_tabs]
        removed = sum(before) - sum(after)
        remain = sum(after)
        score = removed + 0.25 * remain
        candidates.append((score, gi, files, before, after))
    if not candidates:
        raise RuntimeError("No 3-exposure groups found.")
    candidates.sort(reverse=True)
    if group_index is None:
        _, gi, files, before, after = candidates[0]
    else:
        found = [item for item in candidates if item[1] == group_index]
        if not found:
            raise RuntimeError(f"group-index={group_index} is not a valid 3-exposure group.")
        _, gi, files, before, after = found[0]
    return gi, files, before, after


def load_stage_tables(before_dir: Path, gaia_dir: Path, files: list[str]):
    raw_tabs = [read_l2_catalog(str(before_dir / name)) for name in files]
    gaia_tabs = [read_l2_catalog(str(gaia_dir / name)) for name in files]
    static_tabs = remove_static_sources(gaia_tabs, r_arcsec=2.0, min_repeat=2)
    return raw_tabs, gaia_tabs, static_tabs


def axis_style(ax):
    ax.minorticks_on()
    ax.tick_params(which="major", direction="in", bottom=True, left=True, labelsize=12, length=6, width=1.2)
    ax.tick_params(which="minor", direction="in", bottom=True, left=True, length=3.5, width=1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.invert_xaxis()


def plot_catalog(ax, tab: Table, color: str, title: str, note: str):
    x = np.asarray(tab["X_Win"], dtype=float)
    y = np.asarray(tab["Y_Win"], dtype=float)
    ax.scatter(x, y, s=1.6, c=color, linewidths=0, alpha=0.45, rasterized=True)
    ax.set_title(title, fontsize=15, pad=6)
    ax.text(0.03, 0.97, note, transform=ax.transAxes, ha="left", va="top", fontsize=12, weight="bold")
    ax.set_xlabel("X_Win (pix)")
    ax.set_ylabel("Y_Win (pix)")
    axis_style(ax)


def plot_overlay(ax, tabs: list[Table], title: str, subtitle: str):
    colors = ["tab:blue", "tab:orange", "tab:green"]
    for idx, (tab, color) in enumerate(zip(tabs, colors), start=1):
        x = np.asarray(tab["X_Win"], dtype=float)
        y = np.asarray(tab["Y_Win"], dtype=float)
        ax.scatter(x, y, s=1.8, c=color, linewidths=0, alpha=0.35, rasterized=True, label=f"exp{idx}")
    ax.text(0.03, 0.97, f"{title}\n{subtitle}", transform=ax.transAxes, ha="left", va="top", fontsize=12, weight="bold")
    ax.set_xlabel("X_Win (pix)")
    ax.set_ylabel("Y_Win (pix)")
    axis_style(ax)
    ax.legend(fontsize=12, loc="lower left", frameon=False)


def main() -> None:
    args = parse_args()
    plt.rc("font", family="Times New Roman", size=18, weight="bold")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    before_dir = Path(args.processed_root) / args.night / "L2"
    gaia_dir = Path(args.root_out) / args.night / "mask_gaia"

    gi, files, before_counts, after_counts = choose_best_group(gaia_dir, args.group_index, args.r_static, args.min_repeat)
    raw_tabs, gaia_tabs, static_tabs = load_stage_tables(before_dir, gaia_dir, files)

    fig = plt.figure(figsize=(18, 12.5))
    spec = gridspec.GridSpec(3, 4, figure=fig, width_ratios=[1, 1, 1, 1.08], height_ratios=[1, 1, 1], wspace=0.22, hspace=0.36)

    raw_axes = [fig.add_subplot(spec[0, i]) for i in range(3)]
    gaia_axes = [fig.add_subplot(spec[1, i]) for i in range(3)]
    overlay_before_ax = fig.add_subplot(spec[2, 0:2])
    overlay_after_ax = fig.add_subplot(spec[2, 2:4])
    note_ax = fig.add_subplot(spec[:, 3])
    note_ax.axis("off")

    for idx in range(3):
        plot_catalog(
            raw_axes[idx],
            raw_tabs[idx],
            color="#7a8fb5",
            title=f"Raw catalog {idx + 1}",
            note=f"N = {len(raw_tabs[idx]):,}",
        )
        plot_catalog(
            gaia_axes[idx],
            gaia_tabs[idx],
            color="#3657a7",
            title=f"After Gaia mask {idx + 1}",
            note=f"N = {len(gaia_tabs[idx]):,}",
        )

    plot_overlay(
        overlay_before_ax,
        gaia_tabs,
        title="Three Gaia-masked catalogs overlaid",
        subtitle=f"Before static-source subtraction  |  total = {sum(len(t) for t in gaia_tabs):,}",
    )
    plot_overlay(
        overlay_after_ax,
        static_tabs,
        title="After static-source subtraction",
        subtitle=f"r_static = {args.r_static:.1f}\"  min_repeat = {args.min_repeat}  |  total = {sum(len(t) for t in static_tabs):,}",
    )

    removed_total = sum(len(t) for t in gaia_tabs) - sum(len(t) for t in static_tabs)
    note_lines = [
        "Masking / Static-source Flow",
        f"Night: {args.night}",
        f"Chosen group: {gi:03d}",
        "",
        "Files:",
        *(f"{idx + 1}. {name}" for idx, name in enumerate(files)),
        "",
        f"Gaia-mask totals: {sum(len(t) for t in gaia_tabs):,}",
        f"After static subtraction: {sum(len(t) for t in static_tabs):,}",
        f"Static sources removed: {removed_total:,}",
    ]
    note_ax.text(0.02, 0.98, "\n".join(note_lines), ha="left", va="top", fontsize=15, weight="bold")

    fig.suptitle(
        "Single-night catalog cleaning flow: raw -> Gaia mask -> static-source subtraction",
        fontsize=22,
        fontweight="bold",
        y=0.98,
    )
    fig.text(0.035, 0.78, "Step 1\nRaw catalogs", ha="center", va="center", fontsize=16, fontweight="bold")
    fig.text(0.035, 0.49, "Step 2\nGaia mask", ha="center", va="center", fontsize=16, fontweight="bold")
    fig.text(0.035, 0.19, "Step 3\nOverlay and\nsubtract static", ha="center", va="center", fontsize=16, fontweight="bold")

    out_png = outdir / f"{args.night}_gaia_static_flow.png"
    fig.savefig(out_png, dpi=250)
    plt.close(fig)

    summary = {
        "night": args.night,
        "group_index": gi,
        "files": files,
        "raw_counts": [len(t) for t in raw_tabs],
        "gaia_counts": [len(t) for t in gaia_tabs],
        "static_counts": [len(t) for t in static_tabs],
        "static_removed_total": removed_total,
        "output_png": str(out_png),
    }
    (outdir / f"{args.night}_gaia_static_flow.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
