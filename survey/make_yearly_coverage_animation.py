#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timedelta
import sys
from pathlib import Path

import matplotlib
from PIL import Image

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon

for candidate in [
    Path(__file__).resolve().parents[1],
    Path("/pipeline/xiaoyunao"),
    Path("/Users/yunaoxiao/Desktop/smt_asteroid"),
]:
    if candidate.exists():
        candidate_text = str(candidate)
        if candidate_text not in sys.path:
            sys.path.insert(0, candidate_text)

from survey.config import SchedulerConfig, SiteConfig
from survey.history import create_history_template, save_history
from survey.run_daily import prepare_footprints
from survey.scheduler import StripScheduler


plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.size"] = 14


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Simulate 365 survey nights and render a coverage animation.")
    parser.add_argument("--start-date", default="2025-01-01", help="Start night in YYYY-MM-DD")
    parser.add_argument("--n-nights", type=int, default=365)
    parser.add_argument(
        "--footprints",
        default="/pipeline/xiaoyunao/survey/footprints/survey_fov_footprints_with_visibility.fits",
    )
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--history-out", default=None, help="Optional path for the simulated history FITS")
    parser.add_argument("--gif-name", default="survey_yearly_coverage.gif")
    parser.add_argument("--fps", type=float, default=12.0)
    return parser.parse_args()


def wrap_ra_delta_deg(ra_deg: np.ndarray) -> np.ndarray:
    return (180.0 - np.asarray(ra_deg, dtype=float) + 180.0) % 360.0 - 180.0


def split_wrapped_poly(ra_deg: np.ndarray, dec_deg: np.ndarray) -> list[np.ndarray]:
    lon = np.deg2rad(wrap_ra_delta_deg(ra_deg))
    lat = np.deg2rad(dec_deg)
    breaks = np.where(np.abs(np.diff(lon)) > np.pi)[0] + 1
    indices = np.concatenate(([0], breaks, [len(lon)]))
    segments: list[np.ndarray] = []
    for start, end in zip(indices[:-1], indices[1:]):
        if end - start >= 2:
            segments.append(np.column_stack([lon[start:end], lat[start:end]]))
    return segments


def polygon_from_row(row) -> tuple[np.ndarray, np.ndarray]:
    ra = np.array(
        [row["corner_ra_1"], row["corner_ra_2"], row["corner_ra_3"], row["corner_ra_4"], row["corner_ra_1"]],
        dtype=float,
    )
    dec = np.array(
        [row["corner_dec_1"], row["corner_dec_2"], row["corner_dec_3"], row["corner_dec_4"], row["corner_dec_1"]],
        dtype=float,
    )
    return ra, dec


def frame_duration_s(index: int, total: int, fps: float) -> float:
    base = 1.0 / max(fps, 1e-6)
    edge = max(8, total // 12)
    if index < edge or index >= total - edge:
        return base * 2.8
    if index < edge * 2 or index >= total - edge * 2:
        return base * 1.7
    return base


def build_patches(footprints) -> list[list[Polygon]]:
    patches: list[list[Polygon]] = []
    for row in footprints:
        ra, dec = polygon_from_row(row)
        patches.append([Polygon(segment, closed=False) for segment in split_wrapped_poly(ra, dec)])
    return patches


def render_frame(
    footprints,
    footprint_patches,
    exposure_counts: np.ndarray,
    tonight_ids: set[str],
    night_label: str,
    frame_index: int,
    total_frames: int,
    out_path: Path,
) -> None:
    fig = plt.figure(figsize=(15, 8.8))
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(True, alpha=0.25)

    norm = Normalize(vmin=0.0, vmax=max(1.0, float(exposure_counts.max())))
    cmap = plt.get_cmap("Greys")

    base_patches = []
    base_values = []
    tonight_patches = []
    for row, row_patches, count in zip(footprints, footprint_patches, exposure_counts):
        field_id = str(row["field_id"]).strip()
        if field_id in tonight_ids:
            tonight_patches.extend(row_patches)
        elif count > 0:
            base_patches.extend(row_patches)
            base_values.extend([float(count)] * len(row_patches))

    if base_patches:
        collection = PatchCollection(
            base_patches,
            cmap=cmap,
            norm=norm,
            linewidth=0.18,
            edgecolor="#6e6e6e",
            alpha=0.96,
        )
        collection.set_array(np.asarray(base_values))
        ax.add_collection(collection)
        cbar = fig.colorbar(collection, ax=ax, shrink=0.85, pad=0.08)
        cbar.set_label("Historical exposure count")

    if tonight_patches:
        ax.add_collection(
            PatchCollection(
                tonight_patches,
                facecolor="#3847a8",
                edgecolor="#1f2560",
                linewidth=0.55,
                alpha=0.86,
            )
        )

    ax.set_title(
        f"Yearly Survey Simulation | Night {frame_index + 1:03d}/{total_frames} | {night_label}\n"
        f"Grey: cumulative exposure history   Indigo: fields observed tonight",
        y=1.08,
    )
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    frames_dir = outdir / "frames"
    outdir.mkdir(parents=True, exist_ok=True)
    frames_dir.mkdir(parents=True, exist_ok=True)

    footprints = prepare_footprints(Table.read(args.footprints))
    history = create_history_template(footprints)
    scheduler = StripScheduler(footprints, history, SiteConfig(), SchedulerConfig())
    footprint_patches = build_patches(footprints)

    field_index = {str(fid).strip(): idx for idx, fid in enumerate(footprints["field_id"])}
    start_date = datetime.strptime(args.start_date, "%Y-%m-%d").date()
    frame_paths: list[Path] = []
    durations: list[float] = []

    for offset in range(args.n_nights):
        night = start_date + timedelta(days=offset)
        plan = scheduler.build_plan(night)
        tonight_ids = {str(row["field_id"]).strip() for row in plan}
        night_tag = night.strftime("%Y%m%d")
        for field_id in tonight_ids:
            idx = field_index.get(field_id)
            if idx is None:
                continue
            history["exposure_count"][idx] += 1
            history["last_observed_night"][idx] = night_tag
        scheduler.history = history.copy()

        frame_path = frames_dir / f"{offset + 1:03d}_{night_tag}.png"
        render_frame(
            footprints=footprints,
            footprint_patches=footprint_patches,
            exposure_counts=np.asarray(history["exposure_count"], dtype=float),
            tonight_ids=tonight_ids,
            night_label=night.strftime("%Y-%m-%d"),
            frame_index=offset,
            total_frames=args.n_nights,
            out_path=frame_path,
        )
        frame_paths.append(frame_path)
        durations.append(frame_duration_s(offset, args.n_nights, args.fps))
        print(f"[frame] {offset + 1:03d}/{args.n_nights} {night_tag} fields={len(tonight_ids)}", flush=True)

    gif_path = outdir / args.gif_name
    gif_frames: list[Image.Image] = []
    target_size = None
    for path in frame_paths:
        frame = Image.open(path).convert("P", palette=Image.ADAPTIVE)
        if target_size is None:
            target_size = frame.size
        elif frame.size != target_size:
            frame = frame.resize(target_size)
        gif_frames.append(frame)
    gif_frames[0].save(
        gif_path,
        save_all=True,
        append_images=gif_frames[1:],
        duration=[int(round(value * 1000.0)) for value in durations],
        loop=0,
    )
    print(f"[done] gif={gif_path}", flush=True)

    history_out = Path(args.history_out) if args.history_out else outdir / "simulated_exposure_history.fits"
    save_history(history, history_out)
    print(f"[done] history={history_out}", flush=True)


if __name__ == "__main__":
    main()
