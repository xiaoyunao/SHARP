#!/usr/bin/env python3
"""Generate Fig. 8 from frozen known-object recovery products."""

from __future__ import annotations

import argparse
from pathlib import Path

import healpy as hp
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

from figure_styles import (
    STATISTICS_COLORS,
    add_mollweide_inset_colorbar,
    add_statistics_colorbar,
    apply_statistics_style,
    healpix_scatter_kwargs,
    longitude_to_mollweide_rad,
    save_pdf_png,
    style_mollweide_healpix_axis,
    style_statistics_axis,
)


def panel_label(axis, label: str) -> None:
    axis.text(
        0.03,
        0.96,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        weight="bold",
        bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
        zorder=10,
    )


def recovery_panel(axis, binned: pd.DataFrame, metric: str, xlabel: str, log_x: bool = False) -> None:
    data = binned.loc[binned["metric"].eq(metric)].copy()
    data = data.loc[np.isfinite(pd.to_numeric(data["upper"], errors="coerce"))]
    if data.empty:
        raise ValueError(f"no finite recovery bins for metric={metric}")
    lower = pd.to_numeric(data["lower"], errors="coerce").to_numpy(dtype=float)
    upper = pd.to_numeric(data["upper"], errors="coerce").to_numpy(dtype=float)
    center = np.sqrt(lower * upper) if log_x else (lower + upper) / 2
    fraction = pd.to_numeric(data["recovery_fraction"], errors="coerce").to_numpy(dtype=float)
    low = pd.to_numeric(data["ci95_low"], errors="coerce").to_numpy(dtype=float)
    high = pd.to_numeric(data["ci95_high"], errors="coerce").to_numpy(dtype=float)
    axis.fill_between(center, low, high, color=STATISTICS_COLORS["red"], alpha=0.42, step="mid")
    axis.plot(center, fraction, marker="o", ms=6, color=STATISTICS_COLORS["trend_red"], lw=2.2)
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Nominal recovery fraction")
    axis.set_ylim(0, 1.02)
    if log_x:
        axis.set_xscale("log")
    style_statistics_axis(axis, tick_fontsize=16)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--residuals", type=Path, required=True)
    parser.add_argument("--recovery-binned", type=Path, required=True)
    parser.add_argument("--recovery-by-night", type=Path, required=True)
    parser.add_argument("--night-status", type=Path, required=True)
    parser.add_argument("--output-base", type=Path, required=True)
    parser.add_argument("--figure-data-dir", type=Path, required=True)
    args = parser.parse_args()

    apply_statistics_style()
    residuals = pd.read_parquet(
        args.residuals,
        columns=[
            "night",
            "pred_ra_deg",
            "pred_dec_deg",
            "dra_cosdec_arcsec",
            "ddec_arcsec",
            "radial_residual_arcsec",
        ],
    )
    binned = pd.read_csv(args.recovery_binned)
    nightly = pd.read_csv(args.recovery_by_night, dtype={"night": "string"})
    status = pd.read_csv(args.night_status, dtype={"night": "string"})
    nightly["night"] = nightly["night"].str.zfill(8)
    status["night"] = status["night"].str.zfill(8)
    nightly = nightly.merge(
        status[[column for column in ["night", "known_ades_n", "quality_code", "primary_science_included"] if column in status]],
        on="night",
        how="left",
        validate="one_to_one",
    )
    args.figure_data_dir.mkdir(parents=True, exist_ok=True)

    nside = 64
    finite_sky = np.isfinite(residuals["pred_ra_deg"]) & np.isfinite(residuals["pred_dec_deg"])
    theta = np.deg2rad(90.0 - residuals.loc[finite_sky, "pred_dec_deg"].to_numpy(dtype=float))
    phi = np.deg2rad(residuals.loc[finite_sky, "pred_ra_deg"].to_numpy(dtype=float) % 360.0)
    pixels = hp.ang2pix(nside, theta, phi, nest=False)
    counts = np.bincount(pixels, minlength=hp.nside2npix(nside))
    used = np.flatnonzero(counts)
    pix_theta, pix_phi = hp.pix2ang(nside, used, nest=False)
    sky_data = pd.DataFrame(
        {
            "healpix_nside": nside,
            "healpix_ordering": "RING",
            "pixel": used,
            "ra_deg": np.rad2deg(pix_phi),
            "dec_deg": 90.0 - np.rad2deg(pix_theta),
            "matched_detection_n": counts[used],
        }
    )
    sky_data.to_csv(args.figure_data_dir / "fig08_matched_sky_healpix.csv", index=False)

    edges = np.linspace(-1.0, 1.0, 81)
    signed_counts, _, _ = np.histogram2d(
        pd.to_numeric(residuals["dra_cosdec_arcsec"], errors="coerce"),
        pd.to_numeric(residuals["ddec_arcsec"], errors="coerce"),
        bins=[edges, edges],
    )
    signed_data = pd.DataFrame(
        {
            "dra_bin_center_arcsec": np.repeat((edges[:-1] + edges[1:]) / 2, len(edges) - 1),
            "ddec_bin_center_arcsec": np.tile((edges[:-1] + edges[1:]) / 2, len(edges) - 1),
            "matched_n": signed_counts.astype(int).ravel(),
        }
    )
    signed_data.to_csv(args.figure_data_dir / "fig08_signed_residual_density.csv", index=False)
    binned.to_csv(args.figure_data_dir / "fig08_recovery_binned.csv", index=False)
    nightly.to_csv(args.figure_data_dir / "fig08_known_by_night.csv", index=False)

    figure = plt.figure(figsize=(21, 18.5))
    grid = figure.add_gridspec(3, 3, height_ratios=[1.12, 1.0, 0.72])
    sky = figure.add_subplot(grid[0, 0:2], projection="mollweide")
    signed = figure.add_subplot(grid[0, 2])
    radial = figure.add_subplot(grid[1, 0])
    magnitude = figure.add_subplot(grid[1, 1])
    rate = figure.add_subplot(grid[1, 2])
    timeline = figure.add_subplot(grid[2, :])

    style_mollweide_healpix_axis(sky)
    scatter_kwargs = healpix_scatter_kwargs(vmin=1, vmax=max(10, float(np.percentile(counts[used], 99))))
    scatter = sky.scatter(
        longitude_to_mollweide_rad(sky_data["ra_deg"]),
        np.deg2rad(sky_data["dec_deg"]),
        c=sky_data["matched_detection_n"],
        **scatter_kwargs,
    )
    add_mollweide_inset_colorbar(figure, sky, scatter, label="Matched detections per nside=64 pixel")
    panel_label(sky, f"(a) Matched sky density (N={len(residuals):,}; ICRS)")

    masked = np.ma.masked_where(signed_counts.T <= 0, signed_counts.T)
    mesh = signed.pcolormesh(
        edges,
        edges,
        masked,
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max(2, float(np.nanmax(signed_counts)))),
        shading="auto",
        rasterized=True,
    )
    signed.axhline(0, color="white", linewidth=0.7, alpha=0.65)
    signed.axvline(0, color="white", linewidth=0.7, alpha=0.65)
    signed.set_aspect("equal")
    signed.set_xlabel(r"$\Delta\alpha\cos\delta$ (arcsec)")
    signed.set_ylabel(r"$\Delta\delta$ (arcsec)")
    style_statistics_axis(signed, density=True, tick_fontsize=16)
    add_statistics_colorbar(signed, mesh, "Matched detections per bin")
    panel_label(signed, "(b) Signed residuals")

    radial_values = pd.to_numeric(residuals["radial_residual_arcsec"], errors="coerce").to_numpy(dtype=float)
    radial_values = radial_values[np.isfinite(radial_values)]
    radial.hist(
        radial_values,
        bins=np.linspace(0, 1, 51),
        color=STATISTICS_COLORS["blue"],
        alpha=0.62,
        edgecolor=STATISTICS_COLORS["ink"],
        linewidth=0.45,
    )
    for percentile, color in zip([50, 84, 95], [STATISTICS_COLORS["red"], STATISTICS_COLORS["orange"], STATISTICS_COLORS["green"]]):
        value = float(np.percentile(radial_values, percentile))
        radial.axvline(value, color=color, lw=2.0, linestyle="--", label=f"p{percentile}={value:.3f}″")
    radial.axvline(1.0, color=STATISTICS_COLORS["ink"], lw=2.0, label="1″ truncation")
    radial.set_xlabel("Radial separation (arcsec)")
    radial.set_ylabel("Matched detections")
    radial.legend(fontsize=14, loc="center right")
    style_statistics_axis(radial, tick_fontsize=16)
    panel_label(radial, "(c) Threshold-truncated radial residual")

    recovery_panel(magnitude, binned, "predicted_v_mag", "Predicted V magnitude")
    panel_label(magnitude, "(d) Recovery versus predicted brightness")
    recovery_panel(rate, binned, "sky_rate_arcsec_per_hour", r"Predicted angular rate (arcsec h$^{-1}$)", log_x=True)
    panel_label(rate, "(e) Recovery versus angular rate")

    dates = pd.to_datetime(nightly["night"], format="%Y%m%d")
    matches = pd.to_numeric(nightly["matched_1arcsec_n"], errors="coerce").fillna(0)
    ades = pd.to_numeric(nightly.get("known_ades_n", 0), errors="coerce").fillna(0)
    timeline.plot(dates, matches, color=STATISTICS_COLORS["blue"], lw=1.8, label="1″ matched detections")
    timeline.fill_between(dates, 0, matches, color=STATISTICS_COLORS["blue"], alpha=0.18)
    timeline.plot(dates, ades, color=STATISTICS_COLORS["orange"], lw=1.6, label="ADES observation rows")
    failure = pd.Timestamp("2026-01-11")
    timeline.axvspan(failure - pd.Timedelta(hours=12), failure + pd.Timedelta(hours=12), color=STATISTICS_COLORS["red"], alpha=0.22)
    timeline.annotate(
        "2026-01-11\nupstream WCS failure case",
        xy=(failure, float(matches.loc[dates == failure].iloc[0]) if np.any(dates == failure) else 0),
        xytext=(pd.Timestamp("2026-02-02"), max(matches.max() * 0.78, 1)),
        fontsize=14,
        color=STATISTICS_COLORS["red"],
        arrowprops={"arrowstyle": "->", "color": STATISTICS_COLORS["red"], "lw": 1.4},
    )
    timeline.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    timeline.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    timeline.tick_params(axis="x", rotation=25)
    timeline.set_xlabel("Observing night")
    timeline.set_ylabel("Detection / ADES rows")
    timeline.legend(fontsize=14, loc="upper right")
    style_statistics_axis(timeline, tick_fontsize=15)
    panel_label(timeline, "(f) Nightly production accounting")

    figure.text(
        0.5,
        0.012,
        "Recovery denominator: unique V<22.5 predictions inside the nominal WCS-projected detector rectangle; bad regions, saturation, limiting depth, and ephemeris uncertainty are not yet modeled.",
        ha="center",
        fontsize=15,
        color=STATISTICS_COLORS["mid_grey"],
    )
    figure.tight_layout(rect=[0.015, 0.040, 0.99, 0.99], h_pad=1.7, w_pad=1.4)
    save_pdf_png(figure, args.output_base)


if __name__ == "__main__":
    main()
