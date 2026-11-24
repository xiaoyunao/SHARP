#!/usr/bin/env python3
"""Read only the specified historical L2 catalogs and L1 pixel sections on smtdata."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import warnings

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.wcs import FITSFixedWarning
import numpy as np

warnings.filterwarnings("ignore", category=FITSFixedWarning)


def name(d, aliases):
    return next((n for n in aliases if n in d.names), None)


def column(d, *aliases):
    key = name(d, aliases)
    return np.asarray(d[key]) if key is not None else np.full(len(d), np.nan)


def l2(job):
    path = Path(job["path"])
    expected = job.get("size")
    if expected and path.stat().st_size != expected:
        raise ValueError("L2 size differs from frozen ledger")
    if job.get("sha256"):
        with path.open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != job["sha256"]:
                raise ValueError("L2 hash differs from previous aperture audit")
    with fits.open(path, memmap=False) as hdus:
        h, d = hdus[1].header, hdus[1].data
        if job.get("mjd") is not None and abs(float(h["MJD"])-job["mjd"])*86400 > .1:
            raise ValueError("L2 epoch differs from the frozen ledger")
        ra, dec = column(d, "RA_Win"), column(d, "DEC_Win")
        idx = np.arange(len(d))
        output = {}
        if job.get("apertures"):
            fwhm = float(h["FWHM"])
            correction = np.array([h[f"APER{k:02d}"] for k in range(2, 20, 2)])
            ff = np.array([column(d, f"Flux_Aper{k}") for k in range(1, 10)]).T
            mag = column(d, "Mag_Aper4")
            valid = ((column(d, "Flag", "FLAGS") == 0) & (column(d, "Flag_ISO", "IMAFLAGS_ISO") == 0)
                     & (column(d, "Class_ANN", "CLASS_STAR") > .9) & (mag >= 14) & (mag < 16)
                     & np.isfinite(ff).all(axis=1) & (ff > 0).all(axis=1))
            back = ff[valid]*10**(.4*(correction-correction[-1]))
            ratios = back/back[:, -1, None]
            q = np.percentile(ratios, [16, 50, 84], axis=0)
            zp = float(h["CCDZP"])
            flux4_adu = ff[valid, 3] * 10**(.4*(zp-22.5+correction[3]))
            m4_calc = zp-2.5*np.log10(flux4_adu)+correction[3]
            result = {k: float(h[k]) for k in ("MLIM", "FWHM", "SKY", "CCDZP", "GAIN", "RON")}
            result.update(APER_REF=int(h["APER_REF"]), APER08=float(correction[3]),
                          APER18=float(correction[-1]), stars=int(valid.sum()),
                          max_abs_mag4_formula_residual=float(np.max(np.abs(m4_calc-mag[valid]))))
            output["apertures"] = result
            output["growth"] = [
                dict(diameter_px=k, p16=float(x[0]), median=float(x[1]), p84=float(x[2]))
                for k, x in zip(range(2, 20, 2), q.T)]
        if job.get("targets"):
            cols = [column(d, *names) for names in [
                ("objID", "NUMBER"), ("Mag_Aper4",), ("MagErr_Aper4",), ("Flag", "FLAGS"),
                ("Flag_ISO", "IMAFLAGS_ISO"), ("FWHM",), ("X_Win",), ("Y_Win",)]]
            collected = []
            coverage = []
            wcs = WCS(h, naxis=2)
            target_coordinates = np.array([(t["ra"], t["dec"]) for t in job["targets"]])
            pixels = wcs.all_world2pix(target_coordinates, 0, quiet=True)
            for t, (x, y) in zip(job["targets"], pixels):
                coverage.append({"member": t["member"], "kind": t["kind"],
                                 "inside": bool(4 <= x < 9212 and 4 <= y < 9228)})
                dra = ((ra-t["ra"]+180)%360-180)*np.cos(np.deg2rad(t["dec"]))
                sep = np.hypot(dra, dec-t["dec"])*3600
                for i in idx[sep <= 7]:
                    row = [float(x[i]) if np.isfinite(x[i]) else None for x in cols]
                    collected.append({**t, "sep_arcsec": float(sep[i]), "src_ra": float(ra[i]),
                                      "src_dec": float(dec[i]), "data": row})
            output["sources"] = collected
            output["coverage"] = coverage
        return output


def stamp(job):
    path = job["path"]
    with fits.open(path, memmap=True) as hdus:
        image = next(hdu for hdu in hdus if hdu.header.get("NAXIS") == 2)
        wcs = WCS(image.header).celestial
        ny, nx = image.shape
        stamps = []
        for t in job["targets"]:
            x, y = wcs.world_to_pixel_values(t["ra"], t["dec"])
            if not np.isfinite([x, y]).all():
                continue
            xc, yc = int(round(float(x))), int(round(float(y)))
            if min(xc, yc, nx-xc-1, ny-yc-1) < 20:
                continue
            # Request one contiguous strip in the gzip stream. Requesting each
            # short row separately repeatedly inflates the file from the start.
            strip = image.section[yc-20:yc+21, :]
            data = np.asarray(strip[:, xc-20:xc+21], dtype=float)
            points = [dict(label=q["label"], ra=q["ra"], dec=q["dec"],
                           pixel=np.array(wcs.world_to_pixel_values(q["ra"],q["dec"])).reshape(2).tolist())
                      for q in t["markers"]]
            stamps.append(dict(label=t["label"], pixels=data.tolist(), x0=xc-20, y0=yc-20,
                               markers=points, date_obs=image.header.get("DATE-OBS", "")))
        return {"stamps": stamps}


def run(job):
    try:
        data = stamp(job) if job.get("image") else l2(job)
        return {"path": job["path"], **data}
    except (OSError, KeyError, TypeError, ValueError, IndexError, EOFError) as err:
        return {"path": job["path"], "error": f"{type(err).__name__}: {err}"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=argparse.FileType("r"))
    parser.add_argument("--workers", type=int, default=5)
    args = parser.parse_args()
    jobs = json.load(args.jobs) if args.jobs else JOBS
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for row in pool.map(run, jobs):
            print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
