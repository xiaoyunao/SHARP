#!/usr/bin/env python3
"""Generate Fig. 10 from frozen automatic and post-audit linkages."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import GeocentricMeanEcliptic, SkyCoord
from astropy.time import Time

from figure_styles import (
    STATISTICS_COLORS,
    STATISTICS_HISTOGRAM_STYLE,
    apply_statistics_style,
    longitude_to_mollweide_rad,
    save_pdf_png,
    style_mollweide_healpix_axis,
    style_statistics_axis,
)


def histogram(axis, all_values, retained_values, bins, xlabel):
    all_values = np.asarray(all_values, dtype=float)
    retained_values = np.asarray(retained_values, dtype=float)
    all_values = all_values[np.isfinite(all_values)]
    retained_values = retained_values[np.isfinite(retained_values)]
    axis.hist(
        all_values,
        bins=bins,
        color="#b9c0ca",
        label=f"Automatic catalog (N={len(all_values):,})",
        **STATISTICS_HISTOGRAM_STYLE,
    )
    axis.hist(
        retained_values,
        bins=bins,
        color=STATISTICS_COLORS["red"],
        alpha=0.66,
        histtype="stepfilled",
        edgecolor=STATISTICS_COLORS["ink"],
        linewidth=0.8,
        label=f"Post-audit retained (N={len(retained_values):,})",
    )
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Single-night linkages")
    style_statistics_axis(axis)
    axis.legend(fontsize=17, loc="upper right")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--links", type=Path, required=True)
    parser.add_argument("--output-base", type=Path, required=True)
    args = parser.parse_args()
    apply_statistics_style()
    links = pd.read_csv(args.links)
    retained = links.loc[links["posthoc_retained"].astype(str).str.lower().isin({"true", "1"})].copy()
    if len(links) != 4762 or len(retained) != 58:
        raise ValueError(f"expected 4762/58 linkages, found {len(links)}/{len(retained)}")

    figure = plt.figure(figsize=(18, 14))
    sky = figure.add_subplot(2, 2, 1, projection="mollweide")
    magnitude = figure.add_subplot(2, 2, 2)
    speed = figure.add_subplot(2, 2, 3)
    rms = figure.add_subplot(2, 2, 4)

    style_mollweide_healpix_axis(sky)
    sky.set_yticks(np.deg2rad([-45, -30, -15, 0, 15, 30, 45]))
    sky.scatter(
        longitude_to_mollweide_rad(links["median_ra_deg"]),
        np.deg2rad(links["median_dec_deg"]),
        s=6,
        color="#aeb5bf",
        alpha=0.24,
        linewidths=0,
        rasterized=True,
        label="Automatic catalog",
    )
    sky.scatter(
        longitude_to_mollweide_rad(retained["median_ra_deg"]),
        np.deg2rad(retained["median_dec_deg"]),
        s=30,
        color=STATISTICS_COLORS["red"],
        alpha=0.92,
        edgecolors="white",
        linewidths=0.45,
        zorder=5,
        label="Post-audit retained",
    )
    longitude = np.linspace(0, 360, 1441, endpoint=False)
    ecliptic = SkyCoord(
        lon=longitude * u.deg,
        lat=np.zeros_like(longitude) * u.deg,
        frame=GeocentricMeanEcliptic(equinox=Time("J2000")),
    ).icrs
    x = longitude_to_mollweide_rad(ecliptic.ra.deg)
    order = np.argsort(x)
    sky.plot(x[order], np.deg2rad(ecliptic.dec.deg[order]), color="#f28e2b", lw=2.3, label="Ecliptic (J2000)")
    sky.text(
        0.04,
        0.93,
        "(a) Sky distribution (ICRS)",
        transform=sky.transAxes,
        ha="left",
        va="top",
        fontsize=21,
        weight="bold",
        bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "none", "alpha": 0.8},
    )
    sky.legend(fontsize=15, loc="lower right", bbox_to_anchor=(0.99, 0.07))

    histogram(
        magnitude,
        links["median_mag_aper4"],
        retained["median_mag_aper4"],
        np.arange(13, 21.5, 0.25),
        "Median aperture magnitude",
    )
    magnitude.text(0.03, 0.95, "(b) Photometry", transform=magnitude.transAxes, ha="left", va="top", fontsize=21, weight="bold")
    histogram(
        speed,
        links["speed_arcsec_per_hour"],
        retained["speed_arcsec_per_hour"],
        np.linspace(3, 63, 31),
        r"Linear sky-plane rate (arcsec h$^{-1}$)",
    )
    speed.text(0.03, 0.95, "(c) Sky-plane motion", transform=speed.transAxes, ha="left", va="top", fontsize=21, weight="bold")
    histogram(
        rms,
        links["rms_arcsec"],
        retained["rms_arcsec"],
        np.linspace(0, 2.5, 31),
        "Short-arc fit RMS (arcsec)",
    )
    rms.text(0.03, 0.95, "(d) Consistency-fit residual", transform=rms.transAxes, ha="left", va="top", fontsize=21, weight="bold")

    figure.text(
        0.5,
        0.012,
        "All counts are single-night linkages, not confirmed objects; short-arc orbital elements are not used for population inference.",
        ha="center",
        fontsize=17,
        color=STATISTICS_COLORS["mid_grey"],
    )
    figure.tight_layout(rect=[0.02, 0.045, 0.99, 0.99], h_pad=2.0, w_pad=1.7)
    save_pdf_png(figure, args.output_base)


if __name__ == "__main__":
    main()
