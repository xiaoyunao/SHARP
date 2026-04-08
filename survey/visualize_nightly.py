#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import astropy.units as u
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz, Angle, SkyCoord, GeocentricTrueEcliptic, get_body
from astropy.table import Table
from astropy.time import Time
from astroplan import moon
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon

from .astro_utils import build_site
from .config import SchedulerConfig, SiteConfig

plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.size"] = 18


def wrap_ra_rad(ra_deg: np.ndarray) -> np.ndarray:
    return Angle(ra_deg * u.deg).wrap_at(180 * u.deg).radian


def plot_line(ax, coords: SkyCoord, **kwargs) -> None:
    ra = wrap_ra_rad(coords.ra.deg)
    dec = coords.dec.radian
    breaks = np.where(np.abs(np.diff(ra)) > np.pi)[0] + 1
    idx = np.concatenate(([0], breaks, [len(ra)]))
    for start, end in zip(idx[:-1], idx[1:]):
        ax.plot(ra[start:end], dec[start:end], **kwargs)


def polygon_segments(ra_deg: np.ndarray, dec_deg: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    ra = wrap_ra_rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))
    ra_closed = np.append(ra, ra[0])
    dec_closed = np.append(dec, dec[0])
    breaks = np.where(np.abs(np.diff(ra_closed)) > np.pi)[0] + 1
    idx = np.concatenate(([0], breaks, [len(ra_closed)]))
    return [(ra_closed[start:end], dec_closed[start:end]) for start, end in zip(idx[:-1], idx[1:]) if end - start >= 3]


def add_polygon(ax, ra_deg: np.ndarray, dec_deg: np.ndarray, **kwargs) -> None:
    for seg_ra, seg_dec in polygon_segments(ra_deg, dec_deg):
        ax.add_patch(Polygon(np.column_stack([seg_ra, seg_dec]), closed=True, **kwargs))


def plot_altitude_30_line(ax, obstime: Time, site, min_altitude_deg: float) -> None:
    grid_ra = np.linspace(0.0, 360.0, 240)
    grid_dec = np.linspace(-30.0, 85.0, 160)
    ra_mesh, dec_mesh = np.meshgrid(grid_ra, grid_dec)
    grid = SkyCoord(ra=ra_mesh.ravel() * u.deg, dec=dec_mesh.ravel() * u.deg)
    alt_mesh = grid.transform_to(AltAz(obstime=obstime, location=site)).alt.deg.reshape(ra_mesh.shape)
    cs = ax.contour(ra_mesh, dec_mesh, alt_mesh, levels=[min_altitude_deg])
    labeled = False
    for seg in cs.allsegs[0]:
        if len(seg) < 3:
            continue
        seg_ra = wrap_ra_rad(seg[:, 0])
        seg_dec = np.deg2rad(seg[:, 1])
        breaks = np.where(np.abs(np.diff(seg_ra)) > np.pi)[0] + 1
        idx = np.concatenate(([0], breaks, [len(seg_ra)]))
        for start, end in zip(idx[:-1], idx[1:]):
            if end - start < 2:
                continue
            ax.plot(
                seg_ra[start:end],
                seg_dec[start:end],
                color="#2a9d8f",
                linewidth=1.2,
                linestyle="dashed",
                label=f"Alt > {min_altitude_deg:.0f}°" if not labeled else None,
            )
            labeled = True
    if hasattr(cs, "collections"):
        for coll in cs.collections:
            coll.remove()
    elif hasattr(cs, "remove"):
        cs.remove()


def polygon_arrays(row) -> tuple[np.ndarray, np.ndarray]:
    ra_deg = np.array([row["corner_ra_1"], row["corner_ra_2"], row["corner_ra_3"], row["corner_ra_4"]], dtype=float)
    dec_deg = np.array([row["corner_dec_1"], row["corner_dec_2"], row["corner_dec_3"], row["corner_dec_4"]], dtype=float)
    return ra_deg, dec_deg


def _dedupe_labels(ax) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    by_label: dict[str, object] = {}
    for handle, label in zip(handles, labels):
        if label not in by_label:
            by_label[label] = handle
    ax.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=10)


def _style_axes(ax) -> None:
    ax.grid(True, linestyle=":", alpha=0.35)
    ax.set_xticks(np.radians([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150]))
    ax.set_xticklabels(["210", "240", "270", "300", "330", "0", "30", "60", "90", "120", "150"])


def _draw_background(ax, footprints: Table, bg_count: np.ndarray, cmap, norm) -> None:
    for idx, row in enumerate(footprints):
        count = bg_count[idx]
        if count <= 0:
            continue
        ra_deg, dec_deg = polygon_arrays(row)
        add_polygon(ax, ra_deg, dec_deg, facecolor=cmap(norm(count)), edgecolor="none", linewidth=0.0, alpha=0.95)


def generate_plan_visualizations(
    plan: list[dict],
    history: Table,
    footprints: Table,
    outdir: Path,
    night_tag: str,
    site_config: SiteConfig | None = None,
    scheduler_config: SchedulerConfig | None = None,
) -> dict[str, object]:
    outdir.mkdir(parents=True, exist_ok=True)
    site_config = site_config or SiteConfig()
    scheduler_config = scheduler_config or SchedulerConfig()
    site = build_site(site_config.lat_deg, site_config.lon_deg, site_config.height_m)
    fp_by_id = {str(row["field_id"]): row for row in footprints}
    counts = {str(fid): int(cnt) for fid, cnt in zip(history["field_id"], history["exposure_count"])}
    bg_count = np.array([counts.get(str(fid), 0) for fid in footprints["field_id"]], dtype=float)
    norm = Normalize(vmin=0.0, vmax=max(float(bg_count.max()), 1.0))
    cmap = cm.get_cmap("Greys")
    lon_vals = np.linspace(0.0, 360.0, 720) * u.deg

    cycle_paths: list[Path] = []
    groups: dict[int, list[dict]] = {}
    for row in plan:
        groups.setdefault(int(row["cycle_id"]), []).append(row)

    for cycle_id in sorted(groups):
        rows = groups[cycle_id]
        times = Time([row["start_utc"] for row in rows], format="isot", scale="utc")
        obstime = times[0] + 0.5 * (times[-1] - times[0])
        moon_coord = get_body("moon", obstime, location=site)
        moon_phase = moon.moon_illumination(obstime)
        fig = plt.figure(figsize=(14, 8))
        ax = fig.add_subplot(111, projection="aitoff")
        plt.subplots_adjust(top=0.88, bottom=0.18)
        _style_axes(ax)
        _draw_background(ax, footprints, bg_count, cmap, norm)

        for idx, row in enumerate(rows):
            sel = fp_by_id.get(str(row["field_id"]))
            if sel is None:
                continue
            ra_deg, dec_deg = polygon_arrays(sel)
            add_polygon(
                ax,
                ra_deg,
                dec_deg,
                facecolor="#d94801",
                edgecolor="#7f2704",
                linewidth=0.6,
                alpha=0.35 if idx else 0.48,
                label="Planned FoVs" if idx == 0 else None,
            )

        ax.scatter(
            wrap_ra_rad(np.array([moon_coord.ra.deg])),
            np.deg2rad(np.array([moon_coord.dec.deg])),
            color="#ffd166",
            s=90,
            marker="o",
            edgecolor="black",
            linewidth=0.5,
            label=f"Moon {moon_phase:.0%}",
        )

        for lat_deg, style in [(-20.0, "--"), (20.0, "--")]:
            band = SkyCoord(lon=lon_vals, lat=np.full(len(lon_vals), lat_deg) * u.deg, frame=GeocentricTrueEcliptic()).icrs
            plot_line(ax, band, color="#1f78b4", linestyle=style, linewidth=1.1, alpha=0.8, label="Ecliptic ±20°" if lat_deg < 0 else None)

        plot_altitude_30_line(ax, obstime, site, scheduler_config.min_altitude_deg)
        _dedupe_labels(ax)
        ax.set_title(
            f"{night_tag} cycle {cycle_id:02d}/{len(groups)}\nmidpoint={obstime.isot}  n_fields={len(rows)}",
            fontsize=13,
            y=1.08,
        )
        dummy = ax.scatter([], [], c=[], cmap=cmap, s=5, marker="s", vmin=0.0, vmax=max(float(bg_count.max()), 1.0))
        cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.02])
        cbar = fig.colorbar(dummy, cax=cbar_ax, orientation="horizontal")
        cbar.set_label("Historical Exposure Count")
        frame_path = outdir / f"{night_tag}_cycle{cycle_id:02d}.png"
        fig.savefig(frame_path, dpi=160, bbox_inches="tight")
        plt.close(fig)
        cycle_paths.append(frame_path)

    summary_path = outdir / f"{night_tag}_all_fields.png"
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111, projection="aitoff")
    plt.subplots_adjust(top=0.88, bottom=0.18)
    _style_axes(ax)
    _draw_background(ax, footprints, bg_count, cmap, norm)

    unique_ids = []
    seen_ids: set[str] = set()
    for row in plan:
        field_id = str(row["field_id"])
        if field_id not in seen_ids:
            seen_ids.add(field_id)
            unique_ids.append(field_id)

    for idx, field_id in enumerate(unique_ids):
        sel = fp_by_id.get(field_id)
        if sel is None:
            continue
        ra_deg, dec_deg = polygon_arrays(sel)
        add_polygon(
            ax,
            ra_deg,
            dec_deg,
            facecolor="#d94801",
            edgecolor="#7f2704",
            linewidth=0.9,
            alpha=0.40,
            label="Tonight Targets" if idx == 0 else None,
        )

    if plan:
        start_times = Time([row["start_utc"] for row in plan], format="isot", scale="utc")
        obstime = start_times[0] + 0.5 * (start_times[-1] - start_times[0])
        moon_coord = get_body("moon", obstime, location=site)
        moon_phase = moon.moon_illumination(obstime)
        ax.scatter(
            wrap_ra_rad(np.array([moon_coord.ra.deg])),
            np.deg2rad(np.array([moon_coord.dec.deg])),
            color="#ffd166",
            s=90,
            marker="o",
            edgecolor="black",
            linewidth=0.5,
            label=f"Moon {moon_phase:.0%}",
        )
        plot_altitude_30_line(ax, obstime, site, scheduler_config.min_altitude_deg)
    for lat_deg, style in [(-20.0, "--"), (20.0, "--")]:
        band = SkyCoord(lon=lon_vals, lat=np.full(len(lon_vals), lat_deg) * u.deg, frame=GeocentricTrueEcliptic()).icrs
        plot_line(ax, band, color="#1f78b4", linestyle=style, linewidth=1.1, alpha=0.8, label="Ecliptic ±20°" if lat_deg < 0 else None)
    _dedupe_labels(ax)
    ax.set_title(f"{night_tag} all planned fields\nunique_fields={len(unique_ids)}  total_exposures={len(plan)}", fontsize=13, y=1.08)
    dummy = ax.scatter([], [], c=[], cmap=cmap, s=5, marker="s", vmin=0.0, vmax=max(float(bg_count.max()), 1.0))
    cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.02])
    cbar = fig.colorbar(dummy, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Historical Exposure Count")
    fig.savefig(summary_path, dpi=160, bbox_inches="tight")
    plt.close(fig)

    gif_path = outdir / f"{night_tag}_nightly_cycles.gif"
    if cycle_paths:
        with imageio.get_writer(gif_path, mode="I", duration=[1.0] * len(cycle_paths)) as writer:
            for frame_path in cycle_paths:
                writer.append_data(imageio.imread(frame_path))

    manifest = {
        "night_tag": night_tag,
        "cycle_pngs": [str(path) for path in cycle_paths],
        "summary_png": str(summary_path),
        "gif": str(gif_path) if cycle_paths else "",
    }
    manifest_path = outdir / f"{night_tag}_plots.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest
