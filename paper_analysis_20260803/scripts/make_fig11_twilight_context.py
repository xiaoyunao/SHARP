#!/usr/bin/env python3
"""Generate Fig. 11 with an exposure-level twilight denominator."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

from figure_styles import (
    STATISTICS_COLORS,
    add_statistics_colorbar,
    apply_statistics_style,
    save_pdf_png,
    style_statistics_axis,
)


def finite(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    return out.replace([np.inf, -np.inf], np.nan).dropna(subset=columns)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--links", type=Path, required=True)
    parser.add_argument("--exposures", type=Path, required=True)
    parser.add_argument("--output-base", type=Path, required=True)
    parser.add_argument("--figure-data-dir", type=Path, required=True)
    args = parser.parse_args()

    apply_statistics_style()
    links = pd.read_csv(args.links, dtype={"night": "string"})
    retained = links.loc[links["posthoc_retained"].astype(str).str.lower().isin({"true", "1"})].copy()
    exposures = pd.read_csv(args.exposures, dtype={"night": "string"})
    if len(retained) != 58:
        raise ValueError(f"expected 58 retained linkages, found {len(retained)}")
    if exposures.empty:
        raise ValueError("exposure-level denominator is empty")

    twilight_column = "nearest_twilight_signed_hours"
    retained = finite(
        retained,
        [twilight_column, "solar_elongation_deg", "median_mag_aper4", "speed_arcsec_per_hour"],
    )
    exposures = finite(exposures, [twilight_column, "solar_elongation_deg"])
    args.figure_data_dir.mkdir(parents=True, exist_ok=True)
    retained.to_csv(args.figure_data_dir / "fig11_retained_linkages.csv", index=False)
    exposures.to_csv(args.figure_data_dir / "fig11_exposure_denominator.csv", index=False)

    time_edges = np.linspace(-6, 6, 25)
    elongation_edges = np.linspace(20, 180, 33)
    exposure_time_counts, _ = np.histogram(exposures[twilight_column], bins=time_edges)
    candidate_time_counts, _ = np.histogram(retained[twilight_column], bins=time_edges)
    exposure_elongation_counts, _ = np.histogram(exposures["solar_elongation_deg"], bins=elongation_edges)
    candidate_elongation_counts, _ = np.histogram(retained["solar_elongation_deg"], bins=elongation_edges)
    histogram_data = pd.DataFrame(
        {
            "time_bin_left_hours": time_edges[:-1],
            "time_bin_right_hours": time_edges[1:],
            "exposure_n": exposure_time_counts,
            "retained_linkage_n": candidate_time_counts,
        }
    )
    histogram_data.to_csv(args.figure_data_dir / "fig11_twilight_histogram.csv", index=False)
    pd.DataFrame(
        {
            "elongation_bin_left_deg": elongation_edges[:-1],
            "elongation_bin_right_deg": elongation_edges[1:],
            "exposure_n": exposure_elongation_counts,
            "retained_linkage_n": candidate_elongation_counts,
        }
    ).to_csv(args.figure_data_dir / "fig11_elongation_histogram.csv", index=False)

    figure, axes = plt.subplots(2, 2, figsize=(19.5, 14.5))
    twilight, context, magnitude, rate = axes.ravel()

    exposure_weights = np.ones(len(exposures)) / len(exposures)
    retained_weights = np.ones(len(retained)) / len(retained)
    twilight.hist(
        exposures[twilight_column],
        bins=time_edges,
        weights=exposure_weights,
        histtype="stepfilled",
        color="#b9c0ca",
        alpha=0.60,
        edgecolor=STATISTICS_COLORS["ink"],
        linewidth=0.6,
        label=f"All L2 exposures (N={len(exposures):,})",
    )
    twilight.hist(
        retained[twilight_column],
        bins=time_edges,
        weights=retained_weights,
        histtype="step",
        color=STATISTICS_COLORS["red"],
        linewidth=2.8,
        label=f"Retained linkages (N={len(retained):,})",
    )
    twilight.axvline(0, color=STATISTICS_COLORS["ink"], lw=1.4)
    twilight.text(-0.18, 0.97, "toward dawn", transform=twilight.transAxes, ha="right", va="top", fontsize=16)
    twilight.text(0.18, 0.97, "after dusk", transform=twilight.transAxes, ha="left", va="top", fontsize=16)
    twilight.set_xlabel("Signed time from nearest astronomical twilight (h)")
    twilight.set_ylabel("Fraction per 0.5-h bin")
    twilight.set_title("(a) Twilight context with exposure denominator", loc="left", fontsize=20, weight="bold")
    twilight.legend(fontsize=15, loc="upper right")
    style_statistics_axis(twilight, tick_fontsize=18)

    hb = context.hexbin(
        exposures[twilight_column],
        exposures["solar_elongation_deg"],
        gridsize=(42, 42),
        extent=(-6, 6, 20, 180),
        mincnt=1,
        cmap="Greys",
        norm=LogNorm(),
        linewidths=0,
        rasterized=True,
    )
    context.scatter(
        retained[twilight_column],
        retained["solar_elongation_deg"],
        s=42,
        c=STATISTICS_COLORS["red"],
        edgecolors="white",
        linewidths=0.45,
        alpha=0.88,
        label="Retained linkages",
    )
    context.axvline(0, color=STATISTICS_COLORS["ink"], lw=1.1)
    context.set_xlim(-6, 6)
    context.set_ylim(20, 180)
    context.set_xlabel("Signed time from nearest astronomical twilight (h)")
    context.set_ylabel("Solar elongation (deg)")
    context.set_title("(b) Exposure density and retained linkages", loc="left", fontsize=20, weight="bold")
    context.legend(fontsize=15, loc="lower right")
    style_statistics_axis(context, density=True, tick_fontsize=18)
    context_colorbar = add_statistics_colorbar(context, hb, "L2 exposures per hexbin")
    context_colorbar.set_label("L2 exposures per hexbin", fontsize=15)
    context_colorbar.ax.tick_params(labelsize=13)

    mag_scatter = magnitude.scatter(
        retained["solar_elongation_deg"],
        retained["median_mag_aper4"],
        c=retained["speed_arcsec_per_hour"],
        cmap="viridis",
        s=62,
        edgecolors="white",
        linewidths=0.55,
        alpha=0.90,
    )
    magnitude.invert_yaxis()
    magnitude.set_xlabel("Solar elongation (deg)")
    magnitude.set_ylabel("Median aperture magnitude")
    magnitude.set_title("(c) Retained brightness and rate", loc="left", fontsize=20, weight="bold")
    style_statistics_axis(magnitude, density=True, tick_fontsize=18)
    magnitude_colorbar = add_statistics_colorbar(
        magnitude, mag_scatter, r"Sky-plane rate (arcsec h$^{-1}$)"
    )
    magnitude_colorbar.set_label(r"Sky-plane rate (arcsec h$^{-1}$)", fontsize=15)
    magnitude_colorbar.ax.tick_params(labelsize=13)
    denominator_axis = magnitude.twinx()
    centers = (elongation_edges[:-1] + elongation_edges[1:]) / 2
    denominator_axis.plot(
        centers,
        exposure_elongation_counts / max(exposure_elongation_counts.max(), 1),
        color=STATISTICS_COLORS["mid_grey"],
        alpha=0.48,
        lw=1.6,
        drawstyle="steps-mid",
    )
    denominator_axis.set_ylim(0, 4)
    denominator_axis.set_yticks([])
    magnitude.text(
        0.03,
        0.05,
        "gray step: exposure elongation density (scaled)",
        transform=magnitude.transAxes,
        fontsize=12,
        color=STATISTICS_COLORS["mid_grey"],
    )

    rate_scatter = rate.scatter(
        retained[twilight_column],
        retained["speed_arcsec_per_hour"],
        c=retained["median_mag_aper4"],
        cmap="plasma_r",
        s=62,
        edgecolors="white",
        linewidths=0.55,
        alpha=0.90,
    )
    rate.axvline(0, color=STATISTICS_COLORS["ink"], lw=1.1)
    rate.set_xlabel("Signed time from nearest astronomical twilight (h)")
    rate.set_ylabel(r"Sky-plane rate (arcsec h$^{-1}$)")
    rate.set_title("(d) Retained motion versus twilight", loc="left", fontsize=20, weight="bold")
    style_statistics_axis(rate, density=True, tick_fontsize=18)
    rate_colorbar = add_statistics_colorbar(rate, rate_scatter, "Median aperture magnitude")
    rate_colorbar.set_label("Median aperture magnitude", fontsize=15)
    rate_colorbar.ax.tick_params(labelsize=13)

    for axis in (twilight, context, magnitude, rate):
        axis.xaxis.label.set_fontsize(18)
        axis.yaxis.label.set_fontsize(18)
        axis.title.set_fontsize(18)
        axis.tick_params(labelsize=15)

    figure.text(
        0.5,
        0.012,
        "Context only—not a detection-efficiency measurement. Exposure times use midpoint epochs; astronomical twilight uses provisional scheduler/known coordinates (117.575° E, 40.393° N, 960 m) pending surveyed-site confirmation.",
        ha="center",
        fontsize=12,
        color=STATISTICS_COLORS["mid_grey"],
    )
    figure.tight_layout(rect=[0.035, 0.06, 0.985, 0.99], h_pad=3.4, w_pad=4.6)
    save_pdf_png(figure, args.output_base)


if __name__ == "__main__":
    main()
