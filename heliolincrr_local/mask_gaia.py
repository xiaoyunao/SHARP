#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Mask (remove) Gaia sources from nightly L2 catalogs and write new catalogs
with the SAME filenames into /pipeline/xiaoyunao/data/heliolincrr/<night>/mask_gaia/.

Pipeline:
1) Read L2 catalog (OBJ_MP_*_cat.fits.gz)
2) Prefilter: Flag==0 (if exists) and Mag_PSF <= mag_psf_max (if exists)
3) Read CRVAL1/CRVAL2 + DATE-OBS from FITS header (L2 header copied from L1)
4) Load Gaia sources from /pipeline/ref/healpix (healpix-XXXXX.fits, nside=32, ring)
   that intersect a cone around the field center
5) Proper-motion propagate Gaia (pmra/pmdec) from pm_zp (default 2016.0) to DATE-OBS epoch
6) Cross-match within match_radius_arcsec (default 1.5") and remove matched L2 rows
7) Write output catalog (rows removed, columns unchanged) to /pipeline/xiaoyunao/data/heliolincrr/<night>/mask_gaia/<same filename>
8) Write a summary log file

Requirements:
- numpy, astropy
- healpy (for healpix cone selection)

Example:
python mask_gaia.py 20251206 --nproc 6
"""

import os
import glob
import argparse
import numpy as np
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import astropy.units as u
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import FITSFixedWarning

warnings.simplefilter("ignore", FITSFixedWarning)

try:
    import healpy as hp
    _HAS_HEALPY = True
except Exception:
    _HAS_HEALPY = False


# -----------------------------
# Utilities
# -----------------------------
def get_radec_cols(tab: Table):
    """Choose which RA/Dec columns to use in L2 table."""
    if "RA_Win" in tab.colnames and "DEC_Win" in tab.colnames:
        return "RA_Win", "DEC_Win"
    if "RA_PSF" in tab.colnames and "DEC_PSF" in tab.colnames:
        return "RA_PSF", "DEC_PSF"
    if "A" in tab.colnames and "B" in tab.colnames:
        return "A", "B"
    raise KeyError("No recognized RA/Dec columns in L2 catalog (RA_Win/DEC_Win, RA_PSF/DEC_PSF, or A/B).")


def prefilter_l2(tab: Table, mag_psf_max=21.0, require_flag0=True) -> Table:
    """Prefilter L2 sources before Gaia matching."""
    if len(tab) == 0:
        return tab
    m = np.ones(len(tab), dtype=bool)

    if require_flag0 and "Flag" in tab.colnames:
        m &= (tab["Flag"] == 0)

    if (mag_psf_max is not None) and ("Mag_PSF" in tab.colnames):
        m &= np.isfinite(tab["Mag_PSF"])
        m &= (tab["Mag_PSF"] <= float(mag_psf_max))

    return tab[m]


def parse_obsinfo_from_header(l2_fn: str, hdu_header: int = 1):
    """Read CRVAL1/CRVAL2 and DATE-OBS from the FITS header of the catalog file."""
    with fits.open(l2_fn, memmap=False) as hdul:
        hdr = hdul[hdu_header].header if hdu_header < len(hdul) else hdul[0].header

    ra0 = hdr.get("CRVAL1", None)
    dec0 = hdr.get("CRVAL2", None)
    dateobs = hdr.get("DATE-OBS", None)

    if ra0 is None or dec0 is None or dateobs is None:
        raise KeyError(f"Missing CRVAL1/CRVAL2/DATE-OBS in header: {l2_fn}")

    return float(ra0), float(dec0), str(dateobs)


def load_gaia_healpix_cone(ra_deg: float, dec_deg: float, search_radius_deg: float,
                           gaia_dir: str, nside: int = 32) -> Table:
    """
    Load Gaia sources from healpix files that intersect a disc around (ra,dec).
    Files expected: gaia_dir/healpix-XXXXX.fits (ring indexing).
    """
    if not _HAS_HEALPY:
        raise RuntimeError("healpy is required for Gaia healpix cone selection, but it is not installed.")

    c = SkyCoord(ra_deg, dec_deg, unit="deg", frame="icrs")
    phi = np.deg2rad(c.ra.deg)
    theta = np.deg2rad(90.0 - c.dec.deg)
    vec = hp.ang2vec(theta, phi)

    pix = hp.query_disc(nside=nside, vec=vec, radius=np.deg2rad(search_radius_deg))
    pix = np.unique(np.array(pix, dtype=int))
    pix = pix[pix >= 0]

    tabs = []
    for p in pix:
        fname = os.path.join(gaia_dir, f"healpix-{p:05d}.fits")
        if not os.path.isfile(fname):
            continue
        try:
            tabs.append(Table.read(fname))
        except Exception:
            continue

    if len(tabs) == 0:
        return Table()

    return vstack(tabs, metadata_conflicts="silent")


def propagate_gaia_to_epoch(gaia_tab: Table, obstime: Time, pm_zp: float = 2016.0) -> SkyCoord:
    """
    Build Gaia SkyCoord at observation epoch by propagating from pm_zp.
    Uses pmra/pmdec in mas/yr. Gaia pmra is mu_alpha* (includes cos(dec)).
    If pm columns missing, returns coords at catalog positions (no propagation).
    """
    ra = np.array(gaia_tab["ra"], dtype=float)
    dec = np.array(gaia_tab["dec"], dtype=float)

    has_pm = ("pmra" in gaia_tab.colnames) and ("pmdec" in gaia_tab.colnames)
    if not has_pm:
        return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

    pmra = np.array(gaia_tab["pmra"], dtype=float)   # mas/yr, mu_alpha*
    pmdec = np.array(gaia_tab["pmdec"], dtype=float) # mas/yr

    good = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(pmra) & np.isfinite(pmdec)

    ra2 = ra.copy()
    dec2 = dec.copy()

    dt_yr = float(obstime.decimalyear) - float(pm_zp)

    if np.any(good):
        pmra_degyr = (pmra[good] / 1000.0) / 3600.0
        pmdec_degyr = (pmdec[good] / 1000.0) / 3600.0

        dec_rad = np.deg2rad(dec[good])
        cosdec = np.cos(dec_rad)
        cosdec = np.where(cosdec == 0, np.nan, cosdec)

        dra_deg = (pmra_degyr / cosdec) * dt_yr
        ddec_deg = pmdec_degyr * dt_yr

        ra2[good] = (ra[good] + dra_deg) % 360.0
        dec2[good] = dec[good] + ddec_deg

    return SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg, frame="icrs")


# -----------------------------
# Core: mask Gaia from one L2 file
# -----------------------------
def mask_gaia_from_l2_file(
    l2_fn: str,
    out_fn: str,
    gaia_dir: str,
    gaia_cone_deg: float = 2.2,
    match_radius_arcsec: float = 1.5,
    pm_zp: float = 2016.0,
    hdu_header: int = 1,
    mag_psf_max: float = 21.0,
    require_flag0: bool = True,
):
    """
    Returns a dict summary for logging.
    """
    # Read L2 table and prefilter
    tab = Table.read(l2_fn)
    tab = prefilter_l2(tab, mag_psf_max=mag_psf_max, require_flag0=require_flag0)

    n_in = int(len(tab))
    if n_in == 0:
        os.makedirs(os.path.dirname(out_fn), exist_ok=True)
        tab.write(out_fn, overwrite=True)
        return {"file": l2_fn, "status": "empty_after_prefilter", "n_in": 0, "n_out": 0,
                "n_gaia": 0, "n_match": 0}

    ra_col, dec_col = get_radec_cols(tab)

    # Header info for cone + epoch
    ra0, dec0, dateobs = parse_obsinfo_from_header(l2_fn, hdu_header=hdu_header)
    obstime = Time(dateobs, format="isot", scale="utc")

    # Load Gaia in cone
    gaia = load_gaia_healpix_cone(ra0, dec0, gaia_cone_deg, gaia_dir, nside=32)
    n_gaia = int(len(gaia))
    if n_gaia == 0:
        os.makedirs(os.path.dirname(out_fn), exist_ok=True)
        tab.write(out_fn, overwrite=True)
        return {"file": l2_fn, "status": "no_gaia_files", "n_in": n_in, "n_out": n_in,
                "n_gaia": 0, "n_match": 0}

    if ("ra" not in gaia.colnames) or ("dec" not in gaia.colnames):
        os.makedirs(os.path.dirname(out_fn), exist_ok=True)
        tab.write(out_fn, overwrite=True)
        return {"file": l2_fn, "status": "gaia_missing_ra_dec", "n_in": n_in, "n_out": n_in,
                "n_gaia": n_gaia, "n_match": 0}

    # Proper-motion propagation (if pm columns exist)
    c_gaia = propagate_gaia_to_epoch(gaia, obstime, pm_zp=pm_zp)

    # L2 coords
    c_l2 = SkyCoord(np.array(tab[ra_col], dtype=float) * u.deg,
                    np.array(tab[dec_col], dtype=float) * u.deg,
                    frame="icrs")

    # Match L2 -> Gaia and drop hits
    idx, sep, _ = c_l2.match_to_catalog_sky(c_gaia)
    hit = (sep.arcsec <= float(match_radius_arcsec))
    n_match = int(np.sum(hit))

    out_tab = tab[~hit]
    n_out = int(len(out_tab))

    os.makedirs(os.path.dirname(out_fn), exist_ok=True)
    out_tab.write(out_fn, overwrite=True)

    return {"file": l2_fn, "status": "ok", "n_in": n_in, "n_out": n_out,
            "n_gaia": n_gaia, "n_match": n_match}


# -----------------------------
# Top-level worker for multiprocessing (must be picklable on macOS spawn)
# -----------------------------
def mask_gaia_worker(kwargs):
    return mask_gaia_from_l2_file(**kwargs)


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("night", help="e.g. 20251206")
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent / "local_data" / "processed1"), help="input root, contains <night>/L2")
    ap.add_argument("--gaia-dir", default="/pipeline/ref/healpix", help="Gaia healpix directory")
    ap.add_argument("--out-root", default=str(Path(__file__).resolve().parent.parent / "local_data" / "heliolincrr"), help="output root directory")

    ap.add_argument("--gaia-cone-deg", type=float, default=2.2,
                    help="cone radius (deg) for loading Gaia healpix tiles around CRVAL1/2")
    ap.add_argument("--match-arcsec", type=float, default=1.5, help="match radius in arcsec")
    ap.add_argument("--pm-zp", type=float, default=2016.0, help="Gaia reference epoch for proper motion propagation")
    ap.add_argument("--hdu-header", type=int, default=1, help="HDU index to read header keywords from")

    ap.add_argument("--mag-psf-max", type=float, default=21.0, help="prefilter: keep Mag_PSF<=this (if column exists)")
    ap.add_argument("--require-flag0", action="store_true", help="prefilter: require Flag==0 (if column exists)")
    ap.add_argument("--nproc", type=int, default=6, help="parallel processes")
    ap.add_argument("--log", default=None, help="write a summary log file")

    args = ap.parse_args()

    if not _HAS_HEALPY:
        raise SystemExit("ERROR: healpy is not installed, but is required for Gaia healpix cone selection.")

    in_l2 = os.path.join(args.root, args.night, "L2")
    files = sorted(glob.glob(os.path.join(in_l2, "OBJ_MP_*_cat.fits.gz")))

    out_l2 = os.path.join(args.out_root, args.night, "mask_gaia")
    os.makedirs(out_l2, exist_ok=True)

    log_fn = args.log
    if log_fn is None:
        log_fn = os.path.join(args.out_root, args.night, "mask_gaia", "mask_gaia.log")
    os.makedirs(os.path.dirname(log_fn), exist_ok=True)

    packed_tasks = []
    for fn in files:
        out_fn = os.path.join(out_l2, os.path.basename(fn))  # filename unchanged
        packed_tasks.append(dict(
            l2_fn=fn,
            out_fn=out_fn,
            gaia_dir=args.gaia_dir,
            gaia_cone_deg=float(args.gaia_cone_deg),
            match_radius_arcsec=float(args.match_arcsec),
            pm_zp=float(args.pm_zp),
            hdu_header=int(args.hdu_header),
            mag_psf_max=float(args.mag_psf_max),
            require_flag0=bool(args.require_flag0),
        ))

    results = []
    if args.nproc <= 1:
        for d in packed_tasks:
            results.append(mask_gaia_worker(d))
    else:
        with ProcessPoolExecutor(max_workers=args.nproc) as ex:
            futs = [ex.submit(mask_gaia_worker, d) for d in packed_tasks]
            for fut in as_completed(futs):
                results.append(fut.result())

    # Sort by filename for stable logs
    results.sort(key=lambda r: r["file"])

    n_ok = sum(r["status"] == "ok" for r in results)
    n_empty = sum(r["status"] == "empty_after_prefilter" for r in results)

    with open(log_fn, "w", encoding="utf-8") as f:
        f.write(f"night={args.night}\n")
        f.write(f"n_files={len(results)} n_ok={n_ok} n_empty_after_prefilter={n_empty}\n")
        f.write(f"input={in_l2}\noutput={out_l2}\n")
        f.write(f"gaia_dir={args.gaia_dir}\n")
        f.write(f"gaia_cone_deg={args.gaia_cone_deg} match_arcsec={args.match_arcsec} pm_zp={args.pm_zp} hdu_header={args.hdu_header}\n")
        f.write(f"prefilter: Flag==0={args.require_flag0} (if exists), Mag_PSF<={args.mag_psf_max} (if exists)\n")
        f.write(f"nproc={args.nproc}\n")
        f.write("=" * 120 + "\n")
        f.write("status\tin\tout\tgaia\tmatch\tfile\n")
        f.write("-" * 120 + "\n")
        for r in results:
            f.write(f"{r['status']}\t{r['n_in']}\t{r['n_out']}\t{r['n_gaia']}\t{r['n_match']}\t{r['file']}\n")


if __name__ == "__main__":
    main()
