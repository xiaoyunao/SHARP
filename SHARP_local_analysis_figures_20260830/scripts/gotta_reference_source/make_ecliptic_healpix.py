#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import healpy as hp
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from matplotlib.colors import LogNorm


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FITS_PATH = ROOT / "gotta_asteroids.fits"
OUT_PNG = HERE / "gotta_ecliptic_healpix_nside64.png"
OUT_PDF = HERE / "gotta_ecliptic_healpix_nside64.pdf"


def load_columns(path: Path) -> dict[str, np.ndarray]:
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        return {name: np.asarray(data[name]) for name in data.names}


def safe_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones_like(arrays[0], dtype=bool)
    for array in arrays:
        mask &= np.isfinite(array)
    return mask


def lonlat_to_pix(lon_deg: np.ndarray, lat_deg: np.ndarray, nside: int = 64, nest: bool = False) -> np.ndarray:
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)
    theta = 0.5 * np.pi - lat
    phi = lon % (2 * np.pi)
    return hp.ang2pix(nside, theta, phi, nest=nest)


def lon_to_mollweide_rad(lon_deg: np.ndarray, center_deg: float = 180.0) -> np.ndarray:
    centered = lon_deg - center_deg
    centered = (centered + 180.0) % 360.0 - 180.0
    return np.deg2rad(centered)


def set_mollweide_lon_ticks_0_360(ax, center_deg: float = 180.0, step: int = 60) -> None:
    lon_labels = np.arange(step, 301, step)
    lon_ticks = lon_to_mollweide_rad(lon_labels, center_deg=center_deg)
    ax.set_xticks(lon_ticks)
    ax.set_xticklabels([])
    label_lat = np.deg2rad(-15.0)
    for value, lon_tick in zip(lon_labels, lon_ticks):
        ax.text(
            lon_tick,
            label_lat,
            f"{value:d}°",
            ha="center",
            va="center",
            fontsize=24,
            color="black",
            zorder=10,
        )


def icrs_to_ecliptic(ra_deg: np.ndarray, dec_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    ecl = coord.barycentrictrueecliptic
    return ecl.lon.to_value(u.deg), ecl.lat.to_value(u.deg)


def plot_celestial_equator(ax) -> None:
    ra = np.linspace(0.0, 360.0, 2000)
    dec = np.zeros_like(ra)
    lon_ecl, lat_ecl = icrs_to_ecliptic(ra, dec)
    x = lon_to_mollweide_rad(lon_ecl, center_deg=180.0)
    y = np.deg2rad(lat_ecl)
    order = np.argsort(x)
    ax.plot(x[order], y[order], color="#ff8c00", linewidth=3.0, alpha=0.95, zorder=9)


def plot_ecliptic_healpix(columns: dict[str, np.ndarray], out_png: Path, out_pdf: Path, nside: int = 64) -> None:
    ra = np.asarray(columns["ra"], dtype=float)
    dec = np.asarray(columns["dec"], dtype=float)
    mask = safe_mask(ra, dec) & (dec > -90) & (dec < 90)
    lon_ecl, lat_ecl = icrs_to_ecliptic(ra[mask], dec[mask])

    nest = False
    npix = hp.nside2npix(nside)
    pix = lonlat_to_pix(lon_ecl, lat_ecl, nside=nside, nest=nest)
    hp_counts = np.bincount(pix, minlength=npix).astype(np.float64)

    ipix = np.arange(npix)
    theta, phi = hp.pix2ang(nside, ipix, nest=nest)
    lat_c = np.rad2deg(0.5 * np.pi - theta)
    lon_c = np.rad2deg(phi)

    lon = lon_to_mollweide_rad(lon_c, center_deg=180.0)
    lat = np.deg2rad(lat_c)

    cmap = plt.get_cmap("rainbow").copy()
    cmap.set_bad(alpha=0.0)
    counts = hp_counts.astype(np.float64)
    counts_plot = np.clip(counts, 1, 100)
    val = np.ma.masked_where(counts == 0, counts_plot)

    font = {"family": "Times New Roman", "size": 30, "weight": "normal"}
    plt.rc("font", **font)

    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111, projection="mollweide")
    im = ax.scatter(
        lon,
        lat,
        c=val,
        s=8,
        linewidths=0,
        cmap=cmap,
        norm=LogNorm(vmin=1, vmax=100),
        rasterized=True,
    )

    ax.grid(True, alpha=0.3)
    set_mollweide_lon_ticks_0_360(ax, center_deg=180.0, step=60)
    ax.tick_params(axis="y", labelsize=24)
    plot_celestial_equator(ax)

    cax = ax.inset_axes([0.30, 0.09, 0.40, 0.045])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label(r"$\mathrm{Counts\ per\ pixel}$", fontsize=18)
    cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=16, direction="in", length=5, width=1.1)

    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    columns = load_columns(FITS_PATH)
    plot_ecliptic_healpix(columns, OUT_PNG, OUT_PDF, nside=64)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
