#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import numpy as np
import warnings

import astropy.units as u
from astropy.table import Table, vstack
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs import FITSFixedWarning
from astropy.io.fits.verify import VerifyWarning

warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", VerifyWarning)

try:
    from reproject import reproject_interp
    _HAS_REPROJECT = True
except Exception:
    _HAS_REPROJECT = False

try:
    from scipy.ndimage import binary_erosion
    from scipy.ndimage import generate_binary_structure, iterate_structure
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

# ============================================================
# WCS linearization: remove distortion keywords (SIP/TPV/DIS/PV)
# ============================================================

def linearize_wcs_header(hdr_in):
    """
    Return a COPY of header with distortion-related keywords removed.
    Force CTYPE to TAN (e.g., TPV->TAN, SIP->TAN).
    Intended for robust reproject of coverage masks only.
    """
    hdr = hdr_in.copy()

    # Force TAN in CTYPE (common: -TPV or -SIP)
    ctype1 = hdr.get("CTYPE1", "")
    ctype2 = hdr.get("CTYPE2", "")
    if isinstance(ctype1, str):
        hdr["CTYPE1"] = ctype1.replace("-TPV", "-TAN").replace("-SIP", "-TAN")
    if isinstance(ctype2, str):
        hdr["CTYPE2"] = ctype2.replace("-TPV", "-TAN").replace("-SIP", "-TAN")

    # Remove SIP distortion keywords, PV, and DIS families
    sip_prefixes = ("A_", "B_", "AP_", "BP_")
    sip_exact = ("A_ORDER", "B_ORDER", "AP_ORDER", "BP_ORDER")
    pv_prefixes = ("PV1_", "PV2_")
    dis_prefixes = ("CPDIS", "DPDIS", "CQDIS", "DQDIS", "DIS", "D2IM", "WAT")

    keys = list(hdr.keys())
    for k in keys:
        if k is None:
            continue
        kk = str(k).strip()
        if kk in sip_exact:
            try:
                del hdr[kk]
            except Exception:
                pass
            continue
        if kk.startswith(sip_prefixes) or kk.startswith(pv_prefixes) or kk.startswith(dis_prefixes):
            try:
                del hdr[kk]
            except Exception:
                pass

    return hdr


# ============================================================
# Catalog reading
# ============================================================

def read_l2_catalog(fn):
    tab = Table.read(fn)

    # choose RA/Dec columns
    if "RA_Win" in tab.colnames and "DEC_Win" in tab.colnames:
        ra_col, dec_col = "RA_Win", "DEC_Win"
    elif "RA_PSF" in tab.colnames and "DEC_PSF" in tab.colnames:
        ra_col, dec_col = "RA_PSF", "DEC_PSF"
    else:
        ra_col, dec_col = "A", "B"

    keep = ["objID", "MJD", ra_col, dec_col]

    # keep photometry/quality columns if present (NO filtering here)
    if "FWHM" in tab.colnames:
        keep.append("FWHM")
    if "Mag_PSF" in tab.colnames:
        keep.append("Mag_PSF")
    if "MagErr_PSF" in tab.colnames:
        keep.append("MagErr_PSF")

    tab = tab[keep]
    tab.rename_column(ra_col, "RA")
    tab.rename_column(dec_col, "DEC")

    # rename optional columns to consistent internal names
    if "Mag_PSF" in tab.colnames:
        tab.rename_column("Mag_PSF", "MAG_PSF")
    if "MagErr_PSF" in tab.colnames:
        tab.rename_column("MagErr_PSF", "MAGERR_PSF")

    m = np.isfinite(tab["RA"]) & np.isfinite(tab["DEC"])
    return tab[m]

# ============================================================
# Image center and FoV radius (deg)
# ============================================================

def image_center_and_radius(fn):
    with fits.open(fn, memmap=False) as hdul:
        for h in (1, 0):
            if h >= len(hdul):
                continue
            hdr = hdul[h].header
            if "CRVAL1" not in hdr or "CRVAL2" not in hdr:
                continue
            w = WCS(hdr).celestial
            nx, ny = 9216, 9232
            if nx is None or ny is None:
                continue

            ra0, dec0 = w.wcs_pix2world(nx / 2.0, ny / 2.0, 0)

            # pixel scale (deg/pix)
            if w.wcs.cd is not None:
                cd = w.wcs.cd
                pixscale = np.sqrt(np.abs(np.linalg.det(cd)))
            else:
                pixscale = np.mean(np.abs(w.wcs.cdelt))

            fov_r = 0.5 * np.sqrt(nx**2 + ny**2) * pixscale
            return float(ra0), float(dec0), float(fov_r)

    return None, None, None


# ============================================================
# Group exposures by rough overlap (center distance < r1+r2)
# ============================================================

def group_exposures(cat_files, max_sep_deg=0.5):
    """
    Group exposures by image-center separation < max_sep_deg (deg).
    This prevents chain-merging across neighboring fields.
    """
    exps = []
    for fn in cat_files:
        ra, dec, r = image_center_and_radius(fn)
        if ra is None:
            continue
        exps.append((fn, ra, dec, r))

    # Sort by time if available (optional). If no time, keep filename order.
    # (You can add MJD sorting later; not required for the grouping fix.)
    groups = []

    for fn, ra, dec, r in exps:
        c = SkyCoord(ra, dec, unit="deg")
        placed = False

        for g in groups:
            # use the FIRST exposure center as the group's anchor (prevents chain merge)
            _, ra0, dec0, _ = g[0]
            c0 = SkyCoord(ra0, dec0, unit="deg")
            if c.separation(c0).deg < max_sep_deg:
                g.append((fn, ra, dec, r))
                placed = True
                break

        if not placed:
            groups.append([(fn, ra, dec, r)])

    return groups


# ============================================================
# Common coverage mask via reproject, with linear WCS
# ============================================================

def common_mask_from_group_reproject(files, hdu_index=1, erode_pix=0):
    """
    Returns (ref_wcs_lin, common_mask_bool) using reproject with linearized WCS.
    May raise exceptions if headers are badly broken (caller should catch).
    """
    if not _HAS_REPROJECT:
        raise RuntimeError("reproject is not available")

    ref_fn = files[0]
    with fits.open(ref_fn, memmap=False) as hdul:
        ref_hdr0 = hdul[hdu_index].header
        shape_out = (9232, 9216)
        ref_hdr = linearize_wcs_header(ref_hdr0)
        ref_wcs = WCS(ref_hdr).celestial

    masks = []
    for fn in files:
        with fits.open(fn, memmap=False) as hdul:
            hdr0 = hdul[hdu_index].header
            hdr = linearize_wcs_header(hdr0)
            wcs = WCS(hdr).celestial
            img = np.ones((9232, 9216), dtype=np.float32)

        # order=0 nearest neighbor
        try:
            _, footprint = reproject_interp(
                (img, wcs),
                ref_wcs,
                shape_out=shape_out,
                order=0,
                return_footprint=True
            )
        except TypeError:
            _, footprint = reproject_interp(
                (img, wcs),
                ref_wcs,
                shape_out=shape_out,
                order=0
            )

        masks.append(footprint > 0)

    common = np.logical_and.reduce(masks)

    # Optional conservative erosion to absorb edge errors from linearization
    if erode_pix and erode_pix > 0 and _HAS_SCIPY:
        st = generate_binary_structure(2, 2)      # 3x3 with diagonals
        st = iterate_structure(st, int(erode_pix)) # grow structuring element
        common = binary_erosion(common, structure=st)
    else:
        # No scipy: do nothing (still works)
        pass
    
    return ref_wcs, common


def apply_mask(tab, wcs, mask):
    """Keep sources whose sky position falls inside mask in ref_wcs pixel grid."""
    if len(tab) == 0:
        return tab

    sky = SkyCoord(tab["RA"], tab["DEC"], unit="deg")
    x, y = wcs.world_to_pixel(sky)
    x = x.astype(np.int64)
    y = y.astype(np.int64)

    ok = (
        (x >= 0) & (x < mask.shape[1]) &
        (y >= 0) & (y < mask.shape[0])
    )
    idx = np.where(ok)[0]
    keep = np.zeros(len(tab), dtype=bool)
    keep[idx] = mask[y[idx], x[idx]]
    return tab[keep]


# ============================================================
# Fallback common-sky cut (spherical incircle of FoVs)
# ============================================================

def mask_common_sky_incircle(tab, centers, radii, margin_deg=0.1):
    """Fallback: keep sources inside all FoVs (minus margin) in spherical space."""
    if len(tab) == 0:
        return tab

    c_src = SkyCoord(tab["RA"], tab["DEC"], unit="deg")
    keep = np.ones(len(tab), dtype=bool)
    for (ra0, dec0), R in zip(centers, radii):
        c0 = SkyCoord(ra0, dec0, unit="deg")
        keep &= (c_src.separation(c0).deg < (R - margin_deg))
    return tab[keep]


# ============================================================
# Static source removal (robust to empty catalogs)
# ============================================================

def remove_static_sources(tabs, r_arcsec=0.5, min_repeat=2):
    """
    Remove static sources repeating in >=min_repeat exposures within r_arcsec.
    Robust to empty catalogs. Uses bincount to accumulate match counts onto reference catalog.
    """
    if tabs is None or len(tabs) == 0:
        return tabs

    lens = np.array([len(t) for t in tabs], dtype=int)
    if np.all(lens == 0):
        return tabs

    # choose largest non-empty catalog as reference
    ref_idx = int(np.argmax(lens))
    ref = tabs[ref_idx]
    if len(ref) == 0:
        return tabs

    cref = SkyCoord(ref["RA"], ref["DEC"], unit="deg")
    nref = len(ref)

    # counts per ref source: start at 1 to count the reference exposure itself
    counts = np.ones(nref, dtype=int)

    # accumulate matches from other exposures onto reference indices
    for i, tab in enumerate(tabs):
        if i == ref_idx or len(tab) == 0:
            continue

        c = SkyCoord(tab["RA"], tab["DEC"], unit="deg")
        idx, sep, _ = c.match_to_catalog_sky(cref)

        hit = (sep.arcsec < r_arcsec)
        if np.any(hit):
            # add 1 to the matched reference entries
            counts += np.bincount(idx[hit], minlength=nref)

    static_ref = (counts >= int(min_repeat))

    # remove static sources from each tab by matching to reference and masking
    out = []
    for tab in tabs:
        if len(tab) == 0:
            out.append(tab)
            continue

        c = SkyCoord(tab["RA"], tab["DEC"], unit="deg")
        idx, sep, _ = c.match_to_catalog_sky(cref)
        is_static = (sep.arcsec < r_arcsec) & static_ref[idx]
        out.append(tab[~is_static])

    return out


# ============================================================
# Two-point tracklets (spherical)
# ============================================================

def pair_two_exposures(tab1, tab2, vmin, vmax, dmag_max=1.0, file1=None, file2=None):
    if len(tab1) == 0 or len(tab2) == 0:
        return None

    c1 = SkyCoord(tab1["RA"], tab1["DEC"], unit="deg")
    c2 = SkyCoord(tab2["RA"], tab2["DEC"], unit="deg")

    mjd1 = float(np.median(tab1["MJD"]))
    mjd2 = float(np.median(tab2["MJD"]))
    dt_hr = (mjd2 - mjd1) * 24.0
    if not np.isfinite(dt_hr) or dt_hr <= 0:
        return None

    idx, sep, _ = c1.match_to_catalog_sky(c2)
    v = sep.arcsec / dt_hr  # arcsec / hour

    ok = (v > vmin) & (v < vmax)

    # photometry (optional)
    has_mag = ("MAG_PSF" in tab1.colnames) and ("MAG_PSF" in tab2.colnames)
    if has_mag:
        m1 = np.array(tab1["MAG_PSF"], dtype=float)
        m2 = np.array(tab2["MAG_PSF"], dtype=float)
        dmag = m2[idx] - m1
        if dmag_max is not None:
            ok &= np.isfinite(dmag) & (np.abs(dmag) < float(dmag_max))
    else:
        dmag = None  # do not output column if not available

    sel = np.where(ok)[0]
    if sel.size == 0:
        return None

    # optional per-detection columns
    has_fwhm = ("FWHM" in tab1.colnames) and ("FWHM" in tab2.colnames)
    has_merr = ("MAGERR_PSF" in tab1.colnames) and ("MAGERR_PSF" in tab2.colnames)

    # normalize file paths (absolute)
    file1 = os.path.abspath(file1) if file1 is not None else ""
    file2 = os.path.abspath(file2) if file2 is not None else ""

    rows = []
    for i in sel:
        j = idx[i]
        row = [
            float(tab1["RA"][i]), float(tab1["DEC"][i]),
            float(tab2["RA"][j]), float(tab2["DEC"][j]),
            float(v[i]),
            int(tab1["objID"][i]), int(tab2["objID"][j]),
            mjd1, mjd2,
            file1, file2
        ]

        if has_fwhm:
            row += [float(tab1["FWHM"][i]), float(tab2["FWHM"][j])]
        if has_mag:
            row += [float(tab1["MAG_PSF"][i]), float(tab2["MAG_PSF"][j])]
        if has_merr:
            row += [float(tab1["MAGERR_PSF"][i]), float(tab2["MAGERR_PSF"][j])]
        if dmag is not None:
            row += [float(dmag[i])]

        rows.append(tuple(row))

    names = [
        "ra1", "dec1", "ra2", "dec2",
        "v_arcsec_hr",
        "objID1", "objID2",
        "mjd1", "mjd2",
        "file1", "file2"
    ]

    if has_fwhm:
        names += ["fwhm1", "fwhm2"]
    if has_mag:
        names += ["mag_psf1", "mag_psf2"]
    if has_merr:
        names += ["magerr_psf1", "magerr_psf2"]
    if dmag is not None:
        names += ["dmag"]

    t = Table(rows=rows, names=names)

    # force fixed-width unicode columns for FITS safety
    t["file1"] = np.array(t["file1"], dtype="U512")
    t["file2"] = np.array(t["file2"], dtype="U512")

    return t



from scipy.ndimage import binary_erosion

def filter_tracklets_near_mask_edge(
    trk,
    wcs_ref,
    common_mask,
    edge_pix=50
):
    """
    Remove tracklets whose any endpoint lies outside the
    eroded common mask (i.e. within edge_pix of the mask boundary).
    """
    if len(trk) == 0:
        return trk

    # 1) Erode the common mask to define a safe interior region
    safe_mask = binary_erosion(common_mask, iterations=edge_pix)

    ny, nx = safe_mask.shape

    # 2) Project tracklet endpoints to ref pixel coordinates
    c1 = SkyCoord(trk["ra1"], trk["dec1"], unit="deg")
    c2 = SkyCoord(trk["ra2"], trk["dec2"], unit="deg")

    x1, y1 = wcs_ref.world_to_pixel(c1)
    x2, y2 = wcs_ref.world_to_pixel(c2)

    x1 = np.asarray(x1, dtype=int)
    y1 = np.asarray(y1, dtype=int)
    x2 = np.asarray(x2, dtype=int)
    y2 = np.asarray(y2, dtype=int)

    # 3) Check inside bounds
    def inside(x, y):
        return (x >= 0) & (x < nx) & (y >= 0) & (y < ny)

    ok1 = inside(x1, y1)
    ok2 = inside(x2, y2)

    keep = ok1 & ok2
    keep &= safe_mask[y1, x1]
    keep &= safe_mask[y2, x2]

    return trk[keep]

# ============================================================
# Main
# ============================================================

def process_one_group(payload):
    """
    Top-level worker for ProcessPoolExecutor (must be picklable).
    No prints. Returns dict with log string.
    """
    import traceback

    gi, g, ad, outdir = payload
    log_lines = []

    def log(s):
        log_lines.append(s)

    try:
        if len(g) < 2:
            log(f"[group {gi:03d}] skip: len(g)<2")
            return {"gi": gi, "status": "skip_len<2", "ntrk": 0, "log": "\n".join(log_lines)}

        log(f"[group {gi:03d}] start n_exp={len(g)}")

        fns = [x[0] for x in g]
        centers = [(x[1], x[2]) for x in g]
        radii = [x[3] for x in g]

        # --- read catalogs first ---
        raw_tabs = [read_l2_catalog(fn) for fn in fns]

        # --- common-area selection: prefer reproject mask, fallback to spherical incircle ---
        used_fallback = False
        wcs_ref = None
        common = None

        try:
            wcs_ref, common = common_mask_from_group_reproject(
                fns, hdu_index=ad["hdu"], erode_pix=ad["erode_pix"]
            )
            tabs = [apply_mask(tab, wcs_ref, common) for tab in raw_tabs]
            log(f"[group {gi:03d}] common_area=reproject-lin erode_pix={ad['erode_pix']}")
            log(f"[group {gi:03d}] common_mask True fraction={float(np.mean(common)):.6f}")
        except Exception as e:
            used_fallback = True
            tabs = [mask_common_sky_incircle(tab, centers, radii, margin_deg=ad["fallback_margin"])
                    for tab in raw_tabs]
            log(f"[group {gi:03d}] common_area=fallback margin_deg={ad['fallback_margin']:.3f}  err={repr(e)}")

        raw_ns = [len(t) for t in raw_tabs]
        cut_ns = [len(t) for t in tabs]
        log(f"[group {gi:03d}] raw Ns={raw_ns}")
        log(f"[group {gi:03d}] after common-area Ns={cut_ns} ({'fallback' if used_fallback else 'reproj-lin'})")

        if sum(len(t) > 0 for t in tabs) < 2:
            log(f"[group {gi:03d}] skip: fewer than 2 exposures with sources after common-area cut")
            return {"gi": gi, "status": "skip_after_common", "ntrk": 0,
                    "used_fallback": used_fallback, "log": "\n".join(log_lines)}

        # --- remove static sources ---
        tabs = remove_static_sources(tabs, r_arcsec=ad["r_static"], min_repeat=ad["min_repeat"])
        log(f"[group {gi:03d}] after static removal: " +
            ", ".join([f"exp{i}={len(tabs[i])}" for i in range(len(tabs))]))

        if sum(len(t) > 0 for t in tabs) < 2:
            log(f"[group {gi:03d}] skip: fewer than 2 exposures with sources after static removal")
            return {"gi": gi, "status": "skip_after_static", "ntrk": 0,
                    "used_fallback": used_fallback, "log": "\n".join(log_lines)}

        # --- pair consecutive exposures into 2-point tracklets ---
        all_trk = []
        for i in range(len(tabs) - 1):
            t = pair_two_exposures(
                    tabs[i], tabs[i + 1],
                    ad["vmin"], ad["vmax"],
                    dmag_max=ad["dmag_max"],
                    file1=fns[i],
                    file2=fns[i + 1],
                )
            if t is not None:
                t["group"] = gi
                t["exp_i"] = i
                t["exp_j"] = i + 1
                all_trk.append(t)

        if not all_trk:
            log(f"[group {gi:03d}] skip: no tracklets")
            return {"gi": gi, "status": "no_tracklets", "ntrk": 0,
                    "used_fallback": used_fallback, "log": "\n".join(log_lines)}

        trk = vstack(all_trk)

        # edge filter ONLY when reproject succeeded (we have wcs_ref/common)
        if (not used_fallback) and (wcs_ref is not None) and (common is not None) and ad["edge_pix"] > 0:
            n0 = len(trk)
            trk = filter_tracklets_near_mask_edge(trk, wcs_ref, common, edge_pix=ad["edge_pix"])
            log(f"[group {gi:03d}] edge_veto_pix={ad['edge_pix']} : {n0} -> {len(trk)}")
        else:
            log(f"[group {gi:03d}] edge_veto skipped ({'fallback' if used_fallback else 'no_mask'})")

        if len(trk) == 0:
            log(f"[group {gi:03d}] skip: all tracklets filtered")
            return {"gi": gi, "status": "all_filtered", "ntrk": 0,
                    "used_fallback": used_fallback, "log": "\n".join(log_lines)}

        # tracklet_id: yyyymmdd + group(4) + exp_i(2) + exp_j(2) + seq(4), seq within group
        group_id = gi + 1
        exp_i = np.array(trk["exp_i"], dtype=int)
        exp_j = np.array(trk["exp_j"], dtype=int)
        seq = np.arange(1, len(trk) + 1, dtype=int)
        tids = [
            f"{ad['night']}{group_id:04d}{i:02d}{j:02d}{s:04d}"
            for i, j, s in zip(exp_i, exp_j, seq)
        ]
        trk["tracklet_id"] = np.array(tids, dtype="U32")

        out = os.path.join(outdir, f"tracklets_{ad['night']}_g{gi:03d}.fits")
        trk.meta["used_fallback_common_area"] = bool(used_fallback)
        trk.meta["erode_pix"] = int(ad["erode_pix"])
        trk.meta["r_static_arcsec"] = float(ad["r_static"])
        trk.meta["vmin_arcsec_hr"] = float(ad["vmin"])
        trk.meta["vmax_arcsec_hr"] = float(ad["vmax"])
        trk.meta["edge_pix"] = int(ad["edge_pix"])
        trk.write(out, overwrite=True)

        ra_med = float(np.median(trk["ra1"]))
        dec_med = float(np.median(trk["dec1"]))
        log(f"[group {gi:03d}] DONE tracklets={len(trk)} med(RA,Dec)=({ra_med:.3f},{dec_med:.3f}) -> {out}")

        return {"gi": gi, "status": "ok", "ntrk": int(len(trk)),
                "used_fallback": used_fallback, "out": out,
                "log": "\n".join(log_lines)}

    except Exception as e:
        log(f"[group {gi:03d}] ERROR: {repr(e)}")
        log(traceback.format_exc())
        return {"gi": gi, "status": "error", "ntrk": 0, "log": "\n".join(log_lines)}


def main():
    from concurrent.futures import ProcessPoolExecutor, as_completed

    ap = argparse.ArgumentParser()
    ap.add_argument("night", help="e.g. 20251206")
    ap.add_argument("--root", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--input-subdir", default="mask_gaia",
                    help="input subdirectory under <root>/<night>/ that contains Gaia-masked L2 catalogs")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--mag-max", type=float, default=None)

    ap.add_argument("--vmin", type=float, default=0.5)
    ap.add_argument("--vmax", type=float, default=80.0)
    ap.add_argument("--dmag-max", type=float, default=1.0)

    ap.add_argument("--r-static", type=float, default=0.5, help="arcsec")
    ap.add_argument("--min-repeat", type=int, default=2)

    ap.add_argument("--erode-pix", type=int, default=30,
                    help="Erode common mask by N pixels (0 disables). Requires scipy for effect.")
    ap.add_argument("--fallback-margin", type=float, default=0.1,
                    help="deg, only used in spherical fallback common-area cut")

    ap.add_argument("--hdu", type=int, default=1)

    # NEW: parallel + logging
    ap.add_argument("--nproc", type=int, default=6, help="number of processes for groups")
    ap.add_argument("--edge-pix", type=int, default=500, help="tracklet edge veto in pixels (reproject only)")
    ap.add_argument("--log", default=None, help="log file path (default: outdir/make_tracklet.log)")

    args = ap.parse_args()

    l2 = os.path.join(args.root, args.night, args.input_subdir)
    files = sorted(glob.glob(os.path.join(l2, "OBJ_MP_*_cat.fits.gz")))

    if args.outdir is None:
        outdir = os.path.join(args.root, args.night, "tracklets_linreproj")
    else:
        outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    log_fn = args.log if args.log is not None else os.path.join(outdir, "make_tracklet.log")

    groups = group_exposures(files, max_sep_deg=0.5)

    # pack args into plain dict (safe to pickle)
    ad = dict(
        night=args.night,
        vmin=float(args.vmin),
        vmax=float(args.vmax),
        dmag_max=float(args.dmag_max),
        r_static=float(args.r_static),
        min_repeat=int(args.min_repeat),
        erode_pix=int(args.erode_pix),
        fallback_margin=float(args.fallback_margin),
        hdu=int(args.hdu),
        edge_pix=int(args.edge_pix),
    )

    tasks = [(gi, g, ad, outdir) for gi, g in enumerate(groups) if len(g) >= 2]

    # run in parallel and log in main process
    results = {}
    with ProcessPoolExecutor(max_workers=args.nproc) as ex:
        futs = {ex.submit(process_one_group, t): t[0] for t in tasks}
        for fut in as_completed(futs):
            res = fut.result()
            results[res["gi"]] = res

    # write log in group index order for readability
    with open(log_fn, "w", encoding="utf-8") as f:
        f.write(f"night={args.night} nproc={args.nproc} n_groups={len(tasks)}\n")
        f.write(f"outdir={outdir}\n")
        f.write(f"params: vmin={args.vmin} vmax={args.vmax} r_static={args.r_static} "
                f"min_repeat={args.min_repeat} erode_pix={args.erode_pix} "
                f"edge_pix={args.edge_pix} fallback_margin={args.fallback_margin} hdu={args.hdu}\n")
        f.write("=" * 80 + "\n")

        for gi in sorted(results.keys()):
            f.write(results[gi].get("log", "") + "\n")
            f.write("-" * 80 + "\n")


if __name__ == "__main__":
    main()

