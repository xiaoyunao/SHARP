#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import healpy as hp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LogNorm

# ===============================
# 1. User config
# ===============================
FITS_PATH = "/Users/yunaoxiao/Desktop/sitian.fits"
OUTDIR = "/Users/yunaoxiao/Desktop"
os.makedirs(OUTDIR, exist_ok=True)

COL_RA  = "ra"
COL_DEC = "dec"

# Healpix
NSIDE = 64
NEST = False
NPIX = hp.nside2npix(NSIDE)

# ===============================
# 0. Font  (完全按你原来的)
# ===============================
font = {'family': 'Times New Roman', 'size': 18, 'weight': 'bold'}
plt.rc('font', **font)

# ===============================
# 2. Helpers  (完全照你给的写法)
# ===============================
def safe_mask(*arrs):
    m = np.ones_like(arrs[0], dtype=bool)
    for a in arrs:
        m &= np.isfinite(a)
    return m

def radec_to_pix(ra_deg, dec_deg, nside=64, nest=False):
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    theta = 0.5*np.pi - dec
    phi = ra % (2*np.pi)
    return hp.ang2pix(nside, theta, phi, nest=nest)

def ra_to_mollweide_lon_rad(ra_deg, center_deg=180.0):
    lon_deg = (ra_deg - center_deg)  # shift
    lon_deg = (lon_deg + 180.0) % 360.0 - 180.0  # wrap to [-180,180)
    return np.deg2rad(lon_deg)

def set_mollweide_ra_ticks_0_360(ax, center_deg=180.0, step=60):
    # 与你贴的保持一致：不包含0、也不包含360（只标 60..300）
    ra_labels = np.arange(step, 301, step)
    lon_ticks = ra_to_mollweide_lon_rad(ra_labels, center_deg=center_deg)
    ax.set_xticks(lon_ticks)
    ax.set_xticklabels([f"{v:d}" for v in ra_labels], fontsize=15)

def apply_specdis_style(ax):
    ax.minorticks_on()
    ax.tick_params(which="major", direction="in", bottom=True, left=True,
                   labelsize=16, length=7, width=1.3)
    ax.tick_params(which="minor", direction="in", bottom=True, left=True,
                   length=4, width=1.0)

def add_inset_colorbar_specdis(fig, ax, mappable, label):
    cax = inset_axes(ax, width="40%", height="4%", loc="upper right", borderpad=1.2)
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal")
    cb.set_label(label, fontsize=15)
    cb.ax.tick_params(labelsize=15, direction="in", length=6, width=1.2)
    return cb

# ===============================
# 3. Read FITS + Healpix counts
# ===============================
t = Table.read(FITS_PATH)
ra  = np.asarray(t[COL_RA], dtype=float)
dec = np.asarray(t[COL_DEC], dtype=float)

m = safe_mask(ra, dec) & (dec > -90) & (dec < 90)
ra = ra[m]
dec = dec[m]

pix = radec_to_pix(ra, dec, nside=NSIDE, nest=NEST)
hp_counts = np.bincount(pix, minlength=NPIX).astype(np.float64)

# ===============================
# 6. Plot: RA-Dec healpix counts (linear, no-color where 0)
# ===============================
ipix = np.arange(NPIX)
theta, phi = hp.pix2ang(NSIDE, ipix, nest=NEST)
dec_c = np.rad2deg(0.5*np.pi - theta)
ra_c  = np.rad2deg(phi)  # [0,360)

lon = ra_to_mollweide_lon_rad(ra_c, center_deg=180.0)
lat = np.deg2rad(dec_c)

cmap = plt.get_cmap("rainbow").copy()
cmap.set_bad(alpha=0.0)

counts = hp_counts.astype(np.float64)

# 对数：0 必须 mask；显示值 clip 到 [1,200]
counts_plot = np.clip(counts, 1, 200)
val = np.ma.masked_where(counts == 0, counts_plot)

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection="mollweide")

im = ax.scatter(
    lon, lat,
    c=val,
    s=8,
    linewidths=0,
    cmap=cmap,
    norm=LogNorm(vmin=1, vmax=200)
)

ax.grid(True, alpha=0.3)
set_mollweide_ra_ticks_0_360(ax, center_deg=180.0, step=60)

ax.set_title(r"$\mathrm{RA\!-\!Dec\ density\ (healpix,\ nside=64)}$", fontsize=15, pad=20)

# 你原程序这里没调用 apply_specdis_style；如果你想完全一致就不调用
# apply_specdis_style(ax)

add_inset_colorbar_specdis(fig, ax, im, r"$\mathrm{Counts\ per\ pixel}$")

plt.tight_layout()
outpng = os.path.join(OUTDIR, "sitian_ra_dec_healpix_nside64_counts_linear_clim0_200.png")
plt.savefig(outpng, dpi=250)
plt.show()

print("Saved:", outpng)
