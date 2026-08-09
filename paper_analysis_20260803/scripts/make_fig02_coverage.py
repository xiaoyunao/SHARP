#!/usr/bin/env python3
"""Generate Fig. 2 from the deployed footprint and frozen raw manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import GeocentricMeanEcliptic, SkyCoord
from astropy.table import Table
from astropy.time import Time

from figure_styles import (
    add_mollweide_inset_colorbar,
    apply_statistics_style,
    longitude_to_mollweide_rad,
    save_pdf_png,
    style_mollweide_healpix_axis,
)


def field_string(value) -> str:
    if isinstance(value, (bytes, bytearray)):
        value = value.decode("utf-8", "replace")
    return str(value).strip().zfill(4)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--footprints", type=Path, required=True)
    parser.add_argument("--field-visits", type=Path, required=True)
    parser.add_argument("--coverage-summary", type=Path, required=True)
    parser.add_argument("--output-base", type=Path, required=True)
    args = parser.parse_args()
    apply_statistics_style()

    footprint = Table.read(args.footprints)
    visits = pd.read_csv(args.field_visits, dtype={"field_id": "string"})
    visits["field_id"] = visits["field_id"].str.zfill(4)
    summary = json.loads(args.coverage_summary.read_text(encoding="utf-8"))
    foot = pd.DataFrame(
        {
            "field_id": [field_string(value) for value in footprint["field_id"]],
            "ra_deg": np.asarray(footprint["center_ra"], dtype=float),
            "dec_deg": np.asarray(footprint["center_dec"], dtype=float),
        }
    )
    observed = foot.merge(visits, on="field_id", how="inner", validate="one_to_one")

    ecliptic_longitude = np.linspace(0, 360, 1441, endpoint=False)
    ecliptic = SkyCoord(
        lon=ecliptic_longitude * u.deg,
        lat=np.zeros_like(ecliptic_longitude) * u.deg,
        frame=GeocentricMeanEcliptic(equinox=Time("J2000")),
    ).icrs
    ex = longitude_to_mollweide_rad(ecliptic.ra.deg)
    order = np.argsort(ex)

    figure = plt.figure(figsize=(18, 8.7))
    axes = [
        figure.add_subplot(1, 2, 1, projection="mollweide"),
        figure.add_subplot(1, 2, 2, projection="mollweide"),
    ]
    for axis in axes:
        style_mollweide_healpix_axis(axis)
        axis.set_yticks(np.deg2rad([-45, -30, -15, 0, 15, 30, 45]))
        axis.plot(
            ex[order],
            np.deg2rad(ecliptic.dec.deg[order]),
            color="#f28e2b",
            lw=2.4,
            zorder=3,
            label="Ecliptic (J2000)",
        )

    axes[0].scatter(
        longitude_to_mollweide_rad(foot["ra_deg"].to_numpy()),
        np.deg2rad(foot["dec_deg"].to_numpy()),
        s=10,
        color="#b9c0ca",
        alpha=0.62,
        linewidths=0,
        rasterized=True,
        label="Deployed footprint centers",
    )
    axes[0].scatter(
        longitude_to_mollweide_rad(observed["ra_deg"].to_numpy()),
        np.deg2rad(observed["dec_deg"].to_numpy()),
        s=8,
        color="#4c78a8",
        alpha=0.70,
        linewidths=0,
        rasterized=True,
        label="Observed fields",
    )
    axes[0].text(
        0.18,
        0.91,
        "(a) Deployed ecliptic footprint",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=20,
        weight="bold",
    )
    axes[0].text(
        0.18,
        0.80,
        f"{len(foot):,} deployed fields\n{len(observed):,} observed fields",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=16,
        color="#2b2b2b",
        bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
    )
    axes[0].legend(loc="lower right", bbox_to_anchor=(0.98, 0.09), fontsize=16, frameon=False)

    norm = mcolors.LogNorm(vmin=max(1, observed["exposure_n"].min()), vmax=observed["exposure_n"].max())
    points = axes[1].scatter(
        longitude_to_mollweide_rad(observed["ra_deg"].to_numpy()),
        np.deg2rad(observed["dec_deg"].to_numpy()),
        c=observed["exposure_n"].to_numpy(),
        s=13,
        cmap="rainbow",
        norm=norm,
        alpha=0.90,
        linewidths=0,
        rasterized=True,
        zorder=2,
    )
    axes[1].text(
        0.18,
        0.91,
        "(b) Acquired-frame visit density",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=20,
        weight="bold",
    )
    axes[1].text(
        0.18,
        0.80,
        f"N = {summary['raw_exposure_n']:,} exposures\nHEALPix union = {summary['unique_area_deg2']:,.2f} deg²",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=16,
        color="#2b2b2b",
        bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
    )
    add_mollweide_inset_colorbar(
        figure,
        axes[1],
        points,
        label="Acquired exposures per field",
    )
    figure.subplots_adjust(left=0.035, right=0.985, top=0.98, bottom=0.02, wspace=0.08)
    save_pdf_png(figure, args.output_base)
    plt.close(figure)


if __name__ == "__main__":
    main()
