#!/usr/bin/env python3
from __future__ import annotations

import argparse
import io
import math
import re
import warnings
from collections import defaultdict
from pathlib import Path

import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astroplan import moon
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.visualization import ZScaleInterval
from astropy.wcs import FITSFixedWarning, WCS
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Polygon


plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.size"] = 13
warnings.filterwarnings("ignore", category=FITSFixedWarning)

FIELD_RE = re.compile(r"OBJ_MP_(\d{4})_\d+\.fits\.gz$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot nightly known asteroid visualizations.")
    parser.add_argument("night", help="Target night in YYYYMMDD format")
    parser.add_argument("--processed-root", default="/processed1")
    parser.add_argument("--plot-root", default="/pipeline/xiaoyunao/known_asteroid/plots")
    parser.add_argument("--survey-history", default="/pipeline/xiaoyunao/survey/runtime/history/exposure_history.fits")
    parser.add_argument("--survey-footprints", default="/pipeline/xiaoyunao/survey/footprints/survey_fov_footprints_with_visibility.fits")
    parser.add_argument("--survey-plan", default=None)
    parser.add_argument("--all-matched-path", default="/pipeline/xiaoyunao/known_asteroid/runtime/history/all_matched_asteroids.fits")
    parser.add_argument("--top-gif-count", type=int, default=5)
    parser.add_argument("--gif-size", type=int, default=280)
    parser.add_argument("--gif-duration", type=float, default=0.5)
    parser.add_argument("--slew-readout-seconds", type=float, default=0.0)
    return parser.parse_args()


def wrap_ra_delta_deg(ra_deg: np.ndarray) -> np.ndarray:
    return (180.0 - np.asarray(ra_deg, dtype=float) + 180.0) % 360.0 - 180.0


def angular_sep_ra_deg(ra_a: np.ndarray, ra_b: np.ndarray) -> np.ndarray:
    return np.abs(((np.asarray(ra_b, dtype=float) - np.asarray(ra_a, dtype=float) + 180.0) % 360.0) - 180.0)


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


def polygon_from_row(row: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    ra = np.array(
        [row["corner_ra_1"], row["corner_ra_2"], row["corner_ra_3"], row["corner_ra_4"], row["corner_ra_1"]],
        dtype=float,
    )
    dec = np.array(
        [row["corner_dec_1"], row["corner_dec_2"], row["corner_dec_3"], row["corner_dec_4"], row["corner_dec_1"]],
        dtype=float,
    )
    return ra, dec


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def clean_text(value: object) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore").strip()
    return str(value).strip()


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("_") or "asteroid"


def matched_paths(processed_root: Path, night: str) -> tuple[Path, Path, Path]:
    night_dir = processed_root / night
    l4_dir = night_dir / "L4"
    return (
        l4_dir / f"{night}_matched_asteroids.fits",
        l4_dir / f"{night}_all_asteroids.fits",
        l4_dir / f"{night}_matched_asteroids_ades.psv",
    )


def build_observed_fields(l1_dir: Path) -> set[str]:
    observed: set[str] = set()
    for path in sorted(l1_dir.glob("OBJ_MP_*.fits.gz")):
        match = FIELD_RE.fullmatch(path.name)
        if match:
            observed.add(match.group(1))
    return observed


def load_planned_fields(plan_path: Path | None) -> set[str]:
    if plan_path is None or not plan_path.exists():
        return set()
    try:
        plan = pd.read_json(plan_path)
    except Exception:
        return set()
    if "field_id" not in plan.columns:
        return set()
    return {str(x).strip().zfill(4) for x in plan["field_id"].tolist()}


def count_submitted_unique(path: Path) -> int:
    if not path.exists():
        return 0
    unique_keys: set[str] = set()
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#") or text.startswith("!") or text.startswith("permID"):
                continue
            parts = [item.strip() for item in text.split("|")]
            if len(parts) < 2:
                continue
            key = parts[0] or parts[1]
            if key:
                unique_keys.add(key)
    return len(unique_keys)


def dynamic_mag_upper_from_moon_phase(dateobs_values: pd.Series) -> float:
    obs_time = pd.to_datetime(dateobs_values, utc=True, errors="coerce").dropna()
    if obs_time.empty:
        return 21.5
    midpoint = obs_time.iloc[len(obs_time) // 2]
    illum = float(moon.moon_illumination(Time(midpoint.to_pydatetime())))
    return 21.5 - illum * (21.5 - 20.5)


def ra_dec_to_aitoff(ra_deg: np.ndarray, dec_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.deg2rad(wrap_ra_delta_deg(ra_deg)), np.deg2rad(np.asarray(dec_deg, dtype=float))


def load_historical_all(all_matched_path: Path, night: str) -> Table | None:
    if not all_matched_path.exists():
        return None
    table = Table.read(all_matched_path, memmap=True)
    if len(table) == 0:
        return None
    if "source_night" in table.colnames:
        keep = np.asarray([clean_text(x) for x in table["source_night"]]) != night
        table = table[keep]
    return table["RA_Win", "DEC_Win"] if len(table) else None


def plot_survey_coverage(
    footprints: pd.DataFrame,
    history: pd.DataFrame,
    planned_fields: set[str],
    observed_fields: set[str],
    out_path: Path,
    night: str,
) -> None:
    fig = plt.figure(figsize=(15, 8.5))
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(True, alpha=0.3)

    exposure_map = history.set_index("field_id")["exposure_count"].astype(float)
    norm = Normalize(vmin=0.0, vmax=max(1.0, float(exposure_map.max())))
    cmap = plt.get_cmap("Greys")

    base_patches = []
    base_values = []
    planned_patches = []
    observed_patches = []
    for _, row in footprints.iterrows():
        field_id = str(row["field_id"])
        ra, dec = polygon_from_row(row)
        for segment in split_wrapped_poly(ra, dec):
            patch = Polygon(segment, closed=False)
            if field_id in observed_fields:
                observed_patches.append(patch)
            elif field_id in planned_fields:
                planned_patches.append(patch)
            else:
                base_patches.append(patch)
                base_values.append(float(exposure_map.get(field_id, 0.0)))

    if base_patches:
        collection = PatchCollection(base_patches, cmap=cmap, norm=norm, linewidth=0.25, edgecolor="#6b6b6b", alpha=0.95)
        collection.set_array(np.asarray(base_values))
        ax.add_collection(collection)
        cbar = fig.colorbar(collection, ax=ax, shrink=0.82, pad=0.08)
        cbar.set_label("Historical exposure count")

    if planned_patches:
        ax.add_collection(
            PatchCollection(planned_patches, facecolor="#8ecae6", edgecolor="#3a86b3", linewidth=0.6, alpha=0.75)
        )
    if observed_patches:
        ax.add_collection(
            PatchCollection(observed_patches, facecolor="#f4a261", edgecolor="#d62828", linewidth=0.9, alpha=0.45)
        )

    ax.set_title(
        f"Survey Coverage Through {night} | Historical + Planned Fields + Last Night Observed Fields",
        y=1.08,
    )
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_known_asteroid_allsky(historical: Table | None, nightly: Table, out_path: Path, night: str) -> None:
    fig = plt.figure(figsize=(15, 8.5))
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(True, alpha=0.3)

    if historical is not None and len(historical) > 0:
        lon_hist, lat_hist = ra_dec_to_aitoff(historical["RA_Win"], historical["DEC_Win"])
        xb = np.linspace(-np.pi, np.pi, 140)
        yb = np.linspace(-np.pi / 2.0, np.pi / 2.0, 90)
        hist2d, _, _ = np.histogram2d(lon_hist, lat_hist, bins=[xb, yb])
        x_idx = np.clip(np.digitize(lon_hist, xb) - 1, 0, hist2d.shape[0] - 1)
        y_idx = np.clip(np.digitize(lat_hist, yb) - 1, 0, hist2d.shape[1] - 1)
        density = hist2d[x_idx, y_idx]
        positive = density[density > 0]
        norm = LogNorm(vmin=max(1.0, float(np.percentile(positive, 5))), vmax=float(np.percentile(positive, 99.5)))
        ax.scatter(
            lon_hist,
            lat_hist,
            c=density,
            cmap="Greys",
            norm=norm,
            s=7.0,
            alpha=0.28,
            edgecolors="none",
            rasterized=True,
        )
        sm = ScalarMappable(norm=norm, cmap="Greys")
        cbar = fig.colorbar(sm, ax=ax, shrink=0.82, pad=0.08)
        cbar.set_label("Historical counts (log scale)")

    lon_night, lat_night = ra_dec_to_aitoff(nightly["RA_Win"], nightly["DEC_Win"])
    ax.scatter(
        lon_night,
        lat_night,
        s=20,
        c="#00d0ff",
        alpha=0.9,
        marker="+",
        linewidths=0.6,
        label=f"{night} detections",
        rasterized=True,
    )
    ax.legend(loc="upper right", framealpha=0.9)
    ax.set_title(f"Known Asteroid All-Sky Distribution Through {night}", y=1.08)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def astrometric_residual_arcsec(tbl: Table) -> np.ndarray:
    dra = (np.asarray(tbl["RA_Win"]) - np.asarray(tbl["ra"])) * np.cos(np.deg2rad(np.asarray(tbl["DEC_Win"])))
    ddec = np.asarray(tbl["DEC_Win"]) - np.asarray(tbl["dec"])
    return 3600.0 * np.hypot(dra, ddec)


def motion_rate_arcsec_per_min(tbl: Table) -> np.ndarray:
    rates = []
    for _, group in tbl.to_pandas().groupby("name"):
        if len(group) < 2:
            continue
        group = group.sort_values("MJD")
        ra = group["RA_Win"].to_numpy(dtype=float)
        dec = group["DEC_Win"].to_numpy(dtype=float)
        mjd = group["MJD"].to_numpy(dtype=float)
        dt_min = np.diff(mjd) * 1440.0
        valid = dt_min > 0
        if not np.any(valid):
            continue
        dra = np.diff(ra) * np.cos(np.deg2rad(0.5 * (dec[:-1] + dec[1:])))
        ddec = np.diff(dec)
        rate = 3600.0 * np.hypot(dra[valid], ddec[valid]) / dt_min[valid]
        rates.extend(rate.tolist())
    return np.asarray(rates, dtype=float)


def plot_statistics_panel(matched_tbl: Table, all_tbl: Table, submitted_unique: int, out_path: Path, night: str) -> None:
    df = matched_tbl.to_pandas().copy()
    for col in ["name", "DATEOBS", "source_file"]:
        if col in df.columns:
            df[col] = df[col].map(clean_text)
    unique_counts = df.groupby("name").size().sort_values(ascending=False)
    residual = astrometric_residual_arcsec(matched_tbl)
    matched_mag = np.asarray(df["Mag_Aper8"], dtype=float)
    matched_magerr = np.asarray(df["MagErr_Aper8"], dtype=float)
    all_pred_mag = np.asarray(all_tbl["mag"], dtype=float)
    matched_pred_mag = np.asarray(matched_tbl["mag"], dtype=float)
    motion = motion_rate_arcsec_per_min(matched_tbl)
    df["obs_time"] = pd.to_datetime(df["DATEOBS"].str.replace("T", " ", regex=False), errors="coerce", utc=True)
    mag_upper = dynamic_mag_upper_from_moon_phase(df["DATEOBS"])
    timeline = df.dropna(subset=["obs_time"]).copy()
    if not timeline.empty:
        timeline["time_bin"] = timeline["obs_time"].dt.floor("20min")
        time_counts = timeline.groupby("time_bin").size()
        first_seen = timeline.sort_values("obs_time").drop_duplicates("name")
        unique_accum = np.arange(1, len(first_seen) + 1)
    else:
        time_counts = pd.Series(dtype=int)
        first_seen = pd.DataFrame(columns=["obs_time"])
        unique_accum = np.array([], dtype=int)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"Known Asteroid Nightly Statistics | {night}", fontsize=18, y=0.98)

    ax = axes[0, 0]
    bins = np.arange(0.5, unique_counts.max() + 1.6, 1.0)
    ax.hist(unique_counts.to_numpy(), bins=bins, color="#457b9d", edgecolor="white", alpha=0.9)
    median_count = float(np.median(unique_counts))
    ax.axvline(median_count, color="#e63946", linestyle="--", linewidth=2, label=f"Median = {median_count:.1f}")
    ax.set_xlabel("Detections per unique asteroid")
    ax.set_ylabel("Number of asteroids")
    ax.set_title("Detection Multiplicity")
    ax.legend(loc="upper right")
    ax.text(
        0.96,
        0.80,
        f"Unique asteroids: {len(unique_counts)}\nSubmitted unique asteroids: {submitted_unique}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "#cccccc"},
    )

    ax = axes[0, 1]
    res_p95 = max(0.5, np.nanpercentile(residual, 97))
    ax.hist(residual[np.isfinite(residual)], bins=np.linspace(0.0, res_p95, 28), color="#8d99ae", edgecolor="white", alpha=0.95)
    res_med = float(np.nanmedian(residual))
    ax.axvline(res_med, color="#d62828", linestyle="--", linewidth=2, label=f"Median = {res_med:.3f}\"")
    ax.set_xlabel("Astrometric residual (arcsec)")
    ax.set_ylabel("Detections")
    ax.set_title("Measured vs Predicted Position")
    ax.legend(loc="upper right")

    ax = axes[0, 2]
    mag_range = (8.0, mag_upper)
    all_in = np.isfinite(all_pred_mag) & (all_pred_mag >= mag_range[0]) & (all_pred_mag <= mag_range[1])
    matched_in = np.isfinite(matched_pred_mag) & (matched_pred_mag >= mag_range[0]) & (matched_pred_mag <= mag_range[1])
    hist_bins = np.arange(mag_range[0], mag_range[1] + 0.25, 0.25)
    ax.hist(all_pred_mag[all_in], bins=hist_bins, color="#b8c0cc", alpha=0.85, edgecolor="white", label=f"All asteroids: {np.count_nonzero(all_in)}")
    ax.hist(matched_pred_mag[matched_in], bins=hist_bins, color="#f4a261", alpha=0.85, edgecolor="white", label=f"Matched asteroids: {np.count_nonzero(matched_in)}")
    ax.set_xlim(mag_range)
    ax.set_xlabel("Predicted V magnitude")
    ax.set_ylabel("Detections")
    ax.set_title("Brightness Distribution: All vs Matched")
    ax.set_yscale("log")
    ax.legend(loc="upper right")

    ax = axes[1, 0]
    valid = np.isfinite(matched_mag) & np.isfinite(matched_magerr)
    h = ax.hist2d(
        matched_mag[valid],
        matched_magerr[valid],
        bins=[70, 60],
        range=[[8.0, mag_upper], [-0.1, 1.0]],
        cmap="turbo",
        cmin=1,
    )
    cbar = fig.colorbar(h[3], ax=ax, pad=0.01)
    cbar.set_label("Counts")
    ax.set_xlim(8.0, mag_upper)
    ax.set_ylim(-0.1, 1.0)
    ax.set_xlabel("Aperture magnitude (Mag_Aper8)")
    ax.set_ylabel("Photometric error (mag)")
    ax.set_title("Photometric Precision")

    ax = axes[1, 1]
    if motion.size > 0:
        motion_hi = max(1.0, np.nanpercentile(motion, 97))
        ax.hist(motion, bins=np.linspace(0.0, motion_hi, 24), color="#2a9d8f", edgecolor="white", alpha=0.95)
        motion_med = float(np.nanmedian(motion))
        ax.axvline(motion_med, color="#d62828", linestyle="--", linewidth=2, label=f"Median = {motion_med:.3f}")
        ax.legend(loc="upper right")
    ax.set_xlabel("Sky motion rate (arcsec / min)")
    ax.set_ylabel("Asteroids")
    ax.set_title("Apparent Motion")

    ax = axes[1, 2]
    if not time_counts.empty:
        width_days = 12.0 / 1440.0
        ax.bar(time_counts.index, time_counts.to_numpy(), width=width_days, color="#6c8ebf", alpha=0.9)
        ax.set_ylabel("Detections per 20 min")
        ax2 = ax.twinx()
        ax2.plot(first_seen["obs_time"], unique_accum, color="#e76f51", linewidth=2.4)
        ax2.set_ylabel("Cumulative unique asteroids")
        ax2.grid(False)
    ax.set_xlabel("UTC time")
    ax.set_title("Nightly Detection Timeline")
    ax.tick_params(axis="x", rotation=30)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_field_yield_map(footprints: pd.DataFrame, nightly: Table, out_path: Path, night: str) -> None:
    df = nightly.to_pandas()
    if "source_file" in df.columns:
        df["source_file"] = df["source_file"].map(clean_text)
    df["field_id"] = df["source_file"].str.extract(r"OBJ_MP_(\d{4})_")
    yield_map = df.dropna(subset=["field_id"]).groupby("field_id").size()

    fig = plt.figure(figsize=(15, 8.5))
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(True, alpha=0.3)

    active = []
    colors = []
    for _, row in footprints.iterrows():
        field_id = str(row["field_id"])
        if field_id not in yield_map.index:
            continue
        ra, dec = polygon_from_row(row)
        for segment in split_wrapped_poly(ra, dec):
            active.append(Polygon(segment, closed=False))
            colors.append(float(yield_map.loc[field_id]))

    if active:
        collection = PatchCollection(
            active,
            cmap="viridis",
            norm=Normalize(vmin=1.0, vmax=float(max(colors))),
            linewidth=0.35,
            edgecolor="#1f1f1f",
            alpha=0.95,
        )
        collection.set_array(np.asarray(colors))
        ax.add_collection(collection)
        cbar = fig.colorbar(collection, ax=ax, shrink=0.82, pad=0.08)
        cbar.set_label("Matched asteroid detections per observed field")

    ax.set_title(f"Known Asteroid Yield by Observed Field | {night}", y=1.08)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def load_slew_series(header_dir: Path, readout_seconds: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows: list[tuple[pd.Timestamp, float, float, float]] = []
    for path in sorted(header_dir.glob("OBJ_MP_*.fits.gz")):
        try:
            hdr = fits.getheader(path, ext=1)
            start = clean_text(hdr.get("DATE-OBS") or hdr.get("OBS_DATE") or hdr.get("DATE"))
            if not start:
                continue
            start_ts = pd.to_datetime(start, utc=True, errors="coerce")
            if pd.isna(start_ts):
                continue
            exptime = float(hdr.get("EXPOSURE", 30.0))
            ra = float(hdr.get("CEN_RA", hdr.get("TEL_RA", 0.0)))
            dec = float(hdr.get("CEN_DEC", hdr.get("TEL_DEC", 0.0)))
            rows.append((start_ts, exptime, ra, dec))
        except Exception:
            continue
    if len(rows) < 2:
        return np.array([]), np.array([]), np.array([])
    rows.sort(key=lambda item: item[0])
    start = np.array([item[0].value for item in rows], dtype=np.int64) / 1e9
    exptime = np.array([item[1] for item in rows], dtype=float)
    ra = np.array([item[2] for item in rows], dtype=float)
    dec = np.array([item[3] for item in rows], dtype=float)
    delta_ra = angular_sep_ra_deg(ra[:-1], ra[1:])
    delta_dec = np.abs(dec[1:] - dec[:-1])
    available = start[1:] - start[:-1] - exptime[:-1] - float(readout_seconds)
    return delta_ra, delta_dec, np.clip(available, 0.0, None)


def plot_axis_hist(ax: plt.Axes, values: np.ndarray, xlabel: str, title: str, color: str, unit_fmt: str) -> None:
    if values.size > 0:
        ax.hist(values, bins=28, color=color, edgecolor="white", alpha=0.9)
        med = float(np.nanmedian(values))
        ax.axvline(med, color="#d62828", linestyle="--", linewidth=2, label=f"Median = {unit_fmt.format(med)}")
        ax.legend(loc="upper right")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Transitions")
    ax.set_title(title)
    ax.set_yscale("log")


def plot_slew_statistics(header_dir: Path, out_path: Path, night: str, readout_seconds: float) -> None:
    delta_ra, delta_dec, available = load_slew_series(header_dir, readout_seconds)
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle(f"Telescope Axis Slew Statistics | {night}", fontsize=18, y=0.98)

    plot_axis_hist(axes[0, 0], delta_ra, "RA-axis slew angle (deg)", "RA-axis Slew Angle", "#457b9d", "{:.3f}")
    plot_axis_hist(axes[0, 1], available[delta_ra > 0], "RA-axis available slew time (s)", "RA-axis Slew Time", "#8ecae6", "{:.2f}")
    plot_axis_hist(axes[1, 0], delta_dec, "Dec-axis slew angle (deg)", "Dec-axis Slew Angle", "#2a9d8f", "{:.3f}")
    plot_axis_hist(axes[1, 1], available[delta_dec > 0], "Dec-axis available slew time (s)", "Dec-axis Slew Time", "#9c89b8", "{:.2f}")
    axes[0, 1].set_xscale("log")
    axes[1, 1].set_xscale("log")

    note = f"Available slew time = next DATE-OBS - previous DATE-OBS - previous exposure - {readout_seconds:.1f}s readout assumption"
    fig.text(0.5, 0.015, note, ha="center", fontsize=11)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


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
    return imageio.imread(buf)


def build_gifs(tbl: Table, l1_dir: Path, out_dir: Path, night: str, top_n: int, gif_size: int, duration: float) -> None:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in tbl:
        name = clean_text(row["name"])
        fits_name = clean_text(row["source_file"]).replace("_cat.fits.gz", ".fits.gz")
        grouped[name].append(
            {
                "mjd": float(row["MJD"]),
                "dateobs": clean_text(row["DATEOBS"]),
                "ra": float(row["RA_Win"]),
                "dec": float(row["DEC_Win"]),
                "fits_name": fits_name,
            }
        )

    top_targets = sorted(grouped.items(), key=lambda item: len(item[1]), reverse=True)[:top_n]
    for rank, (name, rows) in enumerate(top_targets, start=1):
        rows = sorted(rows, key=lambda item: item["mjd"])
        mid_row = rows[len(rows) // 2]
        center_world = np.array([[mid_row["ra"], mid_row["dec"]]], dtype=float)
        frames = []
        for row in rows:
            fits_path = l1_dir / str(row["fits_name"])
            if not fits_path.exists() or not fits_path.name.endswith(".fits.gz"):
                continue
            with fits.open(fits_path, memmap=False) as hdul:
                data = hdul[1].data.astype(float)
                wcs = WCS(hdul[1].header)
                center_pix = wcs.all_world2pix(center_world, 1)[0]
                marker_pix = wcs.all_world2pix(np.array([[row["ra"], row["dec"]]], dtype=float), 1)[0]
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
            local_time = pd.to_datetime(row["dateobs"], utc=True, errors="coerce")
            if pd.notna(local_time):
                stamp = local_time.tz_convert("Asia/Shanghai").strftime("%Y-%m-%d %H:%M:%S UTC+8")
            else:
                stamp = clean_text(row["dateobs"])
            ax.text(
                0.03,
                0.96,
                name,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=11,
                color="white",
                bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none"},
            )
            ax.text(
                0.03,
                0.04,
                stamp,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=10,
                color="#cdeffd",
                bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none"},
            )
            frames.append(buffer_to_array(fig))
            plt.close(fig)
        if frames:
            frames = frames + [frames[-1], frames[-1]]
            imageio.mimsave(out_dir / f"top{rank:02d}_{sanitize_name(name)}_{night}.gif", frames, duration=duration, loop=1)


def main() -> None:
    args = parse_args()
    night = args.night
    processed_root = Path(args.processed_root)
    matched_path, all_path, ades_path = matched_paths(processed_root, night)
    l1_dir = processed_root / night / "L1"
    l2_dir = processed_root / night / "L2"
    out_dir = Path(args.plot_root) / night
    ensure_dir(out_dir)

    if not matched_path.exists():
        raise SystemExit(f"Missing matched FITS: {matched_path}")
    if not all_path.exists():
        raise SystemExit(f"Missing all-asteroids FITS: {all_path}")
    if not l1_dir.exists():
        raise SystemExit(f"Missing L1 directory: {l1_dir}")
    if not l2_dir.exists():
        raise SystemExit(f"Missing L2 directory: {l2_dir}")

    matched_tbl = Table.read(matched_path)
    if len(matched_tbl) == 0:
        raise SystemExit(f"No rows in {matched_path}")
    all_tbl = Table.read(all_path)

    plan_path = Path(args.survey_plan) if args.survey_plan else Path(f"/pipeline/xiaoyunao/survey/runtime/plans/{night}_plan.json")
    footprints = Table.read(Path(args.survey_footprints)).to_pandas()
    history = Table.read(Path(args.survey_history)).to_pandas()
    footprints["field_id"] = footprints["field_id"].astype(str).str.zfill(4)
    history["field_id"] = history["field_id"].astype(str).str.zfill(4)
    planned_fields = load_planned_fields(plan_path)
    observed_fields = build_observed_fields(l1_dir)
    historical_matched = load_historical_all(Path(args.all_matched_path), night)
    submitted_unique = count_submitted_unique(ades_path)

    plot_survey_coverage(footprints, history, planned_fields, observed_fields, out_dir / f"survey_coverage_{night}.png", night)
    plot_known_asteroid_allsky(historical_matched, matched_tbl, out_dir / f"known_asteroid_allsky_{night}.png", night)
    plot_statistics_panel(matched_tbl, all_tbl, submitted_unique, out_dir / f"known_asteroid_statistics_{night}.png", night)
    plot_field_yield_map(footprints, matched_tbl, out_dir / f"known_asteroid_field_yield_{night}.png", night)
    plot_slew_statistics(l2_dir, out_dir / f"telescope_slew_statistics_{night}.png", night, args.slew_readout_seconds)
    build_gifs(matched_tbl, l1_dir, out_dir, night, args.top_gif_count, args.gif_size, args.gif_duration)


if __name__ == "__main__":
    main()
