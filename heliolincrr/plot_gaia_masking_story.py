#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from PIL import Image


plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.size"] = 13


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a before/after Gaia masking story board for one night.")
    parser.add_argument("night", help="Night in YYYYMMDD")
    parser.add_argument("--processed-root", default="/processed1")
    parser.add_argument("--root-out", default="/pipeline/xiaoyunao/data/heliolincrr")
    parser.add_argument("--filename", default=None, help="Specific catalog filename to visualize")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def load_pair(before_path: Path, after_path: Path) -> tuple[Table, Table]:
    return Table.read(before_path, memmap=True), Table.read(after_path, memmap=True)


def best_filename(before_dir: Path, after_dir: Path) -> str:
    best_name = None
    best_removed = -1
    for before_path in sorted(before_dir.glob("*.fits.gz")):
        after_path = after_dir / before_path.name
        if not after_path.exists():
            continue
        before = Table.read(before_path, memmap=True)
        after = Table.read(after_path, memmap=True)
        removed = len(before) - len(after)
        if removed > best_removed:
            best_removed = removed
            best_name = before_path.name
    if best_name is None:
        raise FileNotFoundError("No matching before/after Gaia masking catalogs found.")
    return best_name


def split_removed(before: Table, after: Table) -> tuple[np.ndarray, np.ndarray]:
    if "objID" in before.colnames and "objID" in after.colnames:
        before_ids = np.asarray(before["objID"])
        after_ids = set(np.asarray(after["objID"]).tolist())
        removed_mask = np.asarray([obj_id not in after_ids for obj_id in before_ids], dtype=bool)
        return removed_mask, ~removed_mask
    before_keys = np.asarray(before["RA_Win"], dtype=float).round(7) * 1e7 + np.asarray(before["DEC_Win"], dtype=float).round(7)
    after_keys = set((np.asarray(after["RA_Win"], dtype=float).round(7) * 1e7 + np.asarray(after["DEC_Win"], dtype=float).round(7)).tolist())
    removed_mask = np.asarray([key not in after_keys for key in before_keys], dtype=bool)
    return removed_mask, ~removed_mask


def render_story(before: Table, after: Table, filename: str, night: str, out_path: Path) -> dict[str, object]:
    removed_mask, kept_mask = split_removed(before, after)
    x_before = np.asarray(before["X_Win"], dtype=float)
    y_before = np.asarray(before["Y_Win"], dtype=float)
    x_after = np.asarray(after["X_Win"], dtype=float)
    y_after = np.asarray(after["Y_Win"], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    panels = [
        ("Before masking", x_before, y_before, "#8d99ae"),
        ("Gaia-matched detections removed", x_before[removed_mask], y_before[removed_mask], "#d1495b"),
        ("After masking", x_after, y_after, "#2b2d7f"),
    ]
    for ax, (title, x, y, color) in zip(axes, panels):
        ax.scatter(x, y, s=2.0, color=color, alpha=0.62, edgecolors="none", rasterized=True)
        ax.set_title(title)
        ax.set_xlabel("X_Win (pix)")
        ax.set_ylabel("Y_Win (pix)")
        ax.set_aspect("equal", adjustable="box")
        ax.invert_xaxis()
    axes[0].text(0.02, 0.98, f"N = {len(before):,}", transform=axes[0].transAxes, va="top", ha="left")
    axes[1].text(0.02, 0.98, f"Removed = {int(removed_mask.sum()):,}", transform=axes[1].transAxes, va="top", ha="left")
    axes[2].text(0.02, 0.98, f"N = {len(after):,}", transform=axes[2].transAxes, va="top", ha="left")
    fig.suptitle(f"Gaia masking story | {night} | {filename}", fontsize=16)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return {
        "night": night,
        "filename": filename,
        "before_rows": int(len(before)),
        "after_rows": int(len(after)),
        "removed_rows": int(removed_mask.sum()),
        "removed_fraction": float(removed_mask.mean()) if len(before) else None,
    }


def render_flip_gif(png_path: Path, out_gif: Path) -> None:
    frame = Image.open(png_path).convert("P", palette=Image.ADAPTIVE)
    frame.save(
        out_gif,
        save_all=True,
        append_images=[frame.copy()],
        duration=[1300, 1300],
        loop=0,
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    before_dir = Path(args.processed_root) / args.night / "L2"
    after_dir = Path(args.root_out) / args.night / "mask_gaia"
    filename = args.filename or best_filename(before_dir, after_dir)
    before, after = load_pair(before_dir / filename, after_dir / filename)
    png_path = outdir / f"{args.night}_gaia_masking_story.png"
    summary = render_story(before, after, filename, args.night, png_path)
    gif_path = outdir / f"{args.night}_gaia_masking_story.gif"
    render_flip_gif(png_path, gif_path)
    (outdir / f"{args.night}_gaia_masking_story.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
