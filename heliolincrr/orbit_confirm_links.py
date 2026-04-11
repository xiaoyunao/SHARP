#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Orbit confirmation / validation for RR links (effective multiprocessing).

Key improvement over naive multiprocessing:
  - Large shared objects (idx, tarr, hypos, loc, thresholds) are initialized ONCE per worker process
    via ProcessPoolExecutor(initializer=..., initargs=...).
  - Tasks only send small payloads (lid, tids).
  - Use executor.map(..., chunksize=...) to reduce scheduling overhead.

Parallelism:
  --cores 1  => serial
  --cores N  => multiprocessing with global worker cache (macOS-safe spawn)

Outputs (default: <rr_dir>/orbit_confirm):
  - orbit_links.fits
  - orbit_obs_residuals.fits
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
from numpy.linalg import norm

from astropy.table import Table
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import EarthLocation, get_body_barycentric, CartesianRepresentation

from poliastro.iod import izzo
from poliastro.bodies import Sun
from poliastro.twobody import Orbit

import warnings
from astropy.io.fits.verify import VerifyWarning
from astropy.utils.iers import IERSDegradedAccuracyWarning
warnings.filterwarnings("ignore", category=VerifyWarning)
warnings.filterwarnings("ignore", category=IERSDegradedAccuracyWarning)
warnings.filterwarnings(
    "ignore",
    message="Tried to get polar motions for times after IERS data is valid.*",
)

from astropy.utils import iers
iers.conf.auto_download = False
iers.conf.iers_degraded_accuracy = "warn"

# ---------- default site: MPC 327 ----------
MPC327_LON_DEG = 117.5750
MPC327_LAT_DEG = 40.394239
MPC327_ALT_M = 868.221

PROFILE_DEFAULTS = {
    "single-night": {
        "min_init_earth_au": 0.02,
        "max_v_kms": 30.0,
        "hypos": [
            (1.3, 0.0, 0.0),
            (1.7, 0.0, 0.0),
            (2.1, 0.0, 0.0),
            (2.5, 0.0, 0.0),
            (2.9, 0.0, 0.0),
            (3.3, 0.0, 0.0),
            (3.7, 0.0, 0.0),
            (4.1, 0.0, 0.0),
        ],
    },
    "w15": {
        "min_init_earth_au": 0.02,
        "max_v_kms": 200.0,
        "hypos": [
            (1.3, -0.02, 0.0),
            (1.3, 0.0, 0.0),
            (1.3, 0.02, 0.0),
            (1.7, -0.02, 0.0),
            (1.7, 0.0, 0.0),
            (1.7, 0.02, 0.0),
            (2.1, -0.02, 0.0),
            (2.1, 0.0, 0.0),
            (2.1, 0.02, 0.0),
            (2.5, -0.02, 0.0),
            (2.5, 0.0, 0.0),
            (2.5, 0.02, 0.0),
            (2.9, -0.02, 0.0),
            (2.9, 0.0, 0.0),
            (2.9, 0.02, 0.0),
            (3.3, -0.02, 0.0),
            (3.3, 0.0, 0.0),
            (3.3, 0.02, 0.0),
            (3.7, -0.02, 0.0),
            (3.7, 0.0, 0.0),
            (3.7, 0.02, 0.0),
            (4.1, -0.02, 0.0),
            (4.1, 0.0, 0.0),
            (4.1, 0.02, 0.0),
        ],
    },
}


def earth_location_mpc327() -> EarthLocation:
    return EarthLocation.from_geodetic(
        lon=MPC327_LON_DEG * u.deg,
        lat=MPC327_LAT_DEG * u.deg,
        height=MPC327_ALT_M * u.m,
    )


def heliocentric_observer_xyz_au(times_utc: Time, loc: EarthLocation) -> np.ndarray:
    """Observer position in heliocentric frame (AU), matching run_rr_from_tracklets.py style.
    Always returns shape (N, 3), even if times_utc is scalar.
    """
    earth_bary = get_body_barycentric("earth", times_utc)
    sun_bary = get_body_barycentric("sun", times_utc)
    topo = loc.get_gcrs(times_utc).cartesian
    topo_au = CartesianRepresentation(topo.x.to(u.au), topo.y.to(u.au), topo.z.to(u.au))
    obs_bary = earth_bary + topo_au
    obs_helio = obs_bary - sun_bary

    x = np.atleast_1d(obs_helio.x.to_value(u.au))
    y = np.atleast_1d(obs_helio.y.to_value(u.au))
    z = np.atleast_1d(obs_helio.z.to_value(u.au))

    return np.stack([x, y, z], axis=1).astype(np.float64)


def linear_motion_test(mjd: np.ndarray, ra: np.ndarray, dec: np.ndarray) -> Tuple[float, float, float]:
    """Simple sky-plane linear motion test: return (lin_rms_arcsec, lin_speed_arcsec/day, lin_dir_deg)."""
    t = mjd - np.mean(mjd)

    ra0 = np.mean(ra)
    dra = (ra - ra0) * 3600.0 * np.cos(np.deg2rad(np.mean(dec)))
    ddec = (dec - np.mean(dec)) * 3600.0

    A = np.vstack([t, np.ones_like(t)]).T
    vr, cr = np.linalg.lstsq(A, dra, rcond=None)[0]
    vd, cd = np.linalg.lstsq(A, ddec, rcond=None)[0]

    dra_fit = vr * t + cr
    ddec_fit = vd * t + cd

    resid = np.hypot(dra - dra_fit, ddec - ddec_fit)
    rms = float(np.sqrt(np.mean(resid**2)))

    speed = float(np.hypot(vr, vd))  # arcsec/day
    direction = float(np.degrees(np.arctan2(vd, vr)) % 360.0)
    return rms, speed, direction


def eq_unitvec(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    """ICRS unit vector for given RA/Dec in degrees."""
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    x = np.cos(ra) * np.cos(dec)
    y = np.sin(ra) * np.cos(dec)
    z = np.sin(dec)
    v = np.stack([x, y, z], axis=1)
    v /= np.clip(np.linalg.norm(v, axis=1, keepdims=True), 1e-30, None)
    return v


def angular_sep_arcsec(u1: np.ndarray, u2: np.ndarray) -> np.ndarray:
    """Angular separation between unit vectors, arcsec."""
    dot = np.sum(u1 * u2, axis=1)
    dot = np.clip(dot, -1.0, 1.0)
    ang = np.arccos(dot)
    return ang * (180.0 / math.pi) * 3600.0


def geod2heliod(r_so: np.ndarray, r_ob_unit: np.ndarray, d: float) -> float:
    """Solve for alpha>0 such that |r_so + alpha*r_ob_unit| == d."""
    a = float(np.sum(r_ob_unit ** 2))
    b = float(2.0 * np.sum(r_ob_unit * r_so))
    c = float(-(d ** 2 - np.sum(r_so ** 2)))
    roots = np.roots([a, b, c])
    roots = roots[np.isreal(roots)].real
    roots = roots[roots > 0]
    if roots.size == 0:
        raise ValueError("no positive root")
    return float(np.min(roots))


def resolve_profile(rr_dir: Path, explicit_profile: str | None) -> str:
    if explicit_profile:
        return explicit_profile
    mode = detect_mode_from_rr_dir(rr_dir)
    if mode == "night":
        return "single-night"
    if mode == "w15":
        return "w15"
    return "single-night"


def build_hypos_default(profile: str) -> List[Tuple[float, float, float]]:
    return [tuple(map(float, hypo)) for hypo in PROFILE_DEFAULTS[profile]["hypos"]]


def detect_mode_from_rr_dir(rr_dir: Path) -> str:
    p = str(rr_dir)
    if p.endswith("/rr_links_15") or "/rr_links_15/" in p:
        return "w15"
    if p.endswith("/rr_links") or "/rr_links/" in p:
        return "night"
    return "unknown"


def extract_night_from_rr_dir(rr_dir: Path) -> str:
    if rr_dir.name in {"rr_links", "rr_links_15"}:
        return rr_dir.parent.name
    return rr_dir.name


def resolve_tracklets_path(rr_dir: Path, root: Path, explicit: str | None) -> Path:
    if explicit:
        return Path(explicit).expanduser()

    mode = detect_mode_from_rr_dir(rr_dir)
    night = extract_night_from_rr_dir(rr_dir)

    if mode == "night":
        return root / night / "tracklets_linreproj" / f"tracklets_{night}_ALL.fits"

    if mode == "w15":
        wdir = root / night / "tracklets_windows_15"
        cand = sorted(wdir.glob(f"tracklets_*_{night}_W15.fits"))
        if not cand:
            raise FileNotFoundError(f"Cannot find W15 tracklets in {wdir} matching *_{night}_W15.fits")
        return cand[-1]

    raise FileNotFoundError("Cannot auto-resolve tracklets file from rr_dir. Please pass --tracklets explicitly.")


def load_rr_links(rr_dir: Path) -> Tuple[Table, Table]:
    links = Table.read(rr_dir / "links_tracklets.fits")
    members = Table.read(rr_dir / "linkage_members.fits")
    return links, members


def build_tracklet_index(tracklets: Table) -> Dict[str, int]:
    tids = np.asarray(tracklets["tracklet_id"]).astype("U64")
    return {tid: i for i, tid in enumerate(tids)}


def count_nights_from_mjd(mjd: np.ndarray) -> int:
    if mjd.size == 0:
        return 0
    mjd = np.sort(mjd)
    return int(np.sum(mjd[1:] - mjd[:-1] > 0.5) + 1)


def fit_best_hypo_lambert(
    mjd: np.ndarray,
    ra: np.ndarray,
    dec: np.ndarray,
    loc: EarthLocation,
    hypos: List[Tuple[float, float, float]],
    min_init_earth_au: float,
    max_v_kms: float,
    outlier_clip_arcsec: float | None = None,
) -> Dict:
    t = Time(mjd, format="mjd", scale="utc")
    jd = t.jd.astype(float)
    ut0 = jd[0]
    dt_days = jd - ut0

    r_so = heliocentric_observer_xyz_au(t, loc)
    u_obs = eq_unitvec(ra, dec)

    i0 = 0
    i1 = len(jd) - 1
    t0 = Time(jd[i0], format="jd", scale="utc")
    t1 = Time(jd[i1], format="jd", scale="utc")
    dt_lambert = t1 - t0

    best = {
        "ok": False,
        "rms_arcsec": np.inf,
        "fail_reason": "no_valid_hypo",
        "fail_counts": {},
        "best_v1_kms": np.nan,
        "min_rejected_max_v_kms": np.nan,
        "best_rejected_max_v_kms": np.nan,
    }
    fail_counts: Dict[str, int] = {}

    def note_fail(reason: str) -> None:
        fail_counts[reason] = fail_counts.get(reason, 0) + 1
        best["fail_reason"] = reason

    for (r0, rdot, rdd) in hypos:
        d0 = float(r0 + rdot * dt_days[i0] + 0.5 * rdd * dt_days[i0] ** 2)
        d1 = float(r0 + rdot * dt_days[i1] + 0.5 * rdd * dt_days[i1] ** 2)
        if (d0 <= 0) or (d1 <= 0):
            note_fail("nonpositive_hypo_distance")
            continue

        try:
            alpha0 = geod2heliod(r_so[i0], u_obs[i0], d0)
            alpha1 = geod2heliod(r_so[i1], u_obs[i1], d1)
        except Exception:
            note_fail("geod2heliod")
            continue

        if (alpha0 <= min_init_earth_au) or (alpha1 <= min_init_earth_au):
            note_fail("min_init_earth")
            continue

        r_helio0 = (r_so[i0] + alpha0 * u_obs[i0]) * u.AU
        r_helio1 = (r_so[i1] + alpha1 * u_obs[i1]) * u.AU

        try:
            velocities = izzo.lambert(Sun.k, r_helio0, r_helio1, dt_lambert, rtol=1e-4)
            try:
                v0, v1 = velocities
            except Exception:
                (v0, v1), = velocities
        except Exception:
            note_fail("lambert")
            continue

        v1_kms = float(norm(v1.to_value(u.km / u.s)))
        if v1_kms >= max_v_kms:
            if (not np.isfinite(best.get("min_rejected_max_v_kms", np.nan))) or (v1_kms < float(best["min_rejected_max_v_kms"])):
                best["min_rejected_max_v_kms"] = v1_kms
            note_fail("max_v")
            continue

        try:
            orb = Orbit.from_vectors(Sun, r_helio1, v1, epoch=t1)
        except Exception:
            note_fail("orbit_from_vectors")
            continue

        dt_vec = Time(jd, format="jd", scale="utc") - t1
        try:
            r_obj = []
            for dti in dt_vec:
                oi = orb.propagate(dti)
                r_obj.append(oi.r.to(u.AU).value)
            r_obj = np.asarray(r_obj, dtype=float)
        except Exception:
            note_fail("propagate")
            continue

        rho = r_obj - r_so
        u_pred = rho / np.clip(norm(rho, axis=1, keepdims=True), 1e-30, None)
        resid = angular_sep_arcsec(u_obs, u_pred)

        if outlier_clip_arcsec is not None and np.isfinite(outlier_clip_arcsec):
            keep = resid <= float(outlier_clip_arcsec)
            if np.sum(keep) < 3:
                note_fail("outlier_clip")
                continue
            resid_use = resid[keep]
        else:
            keep = np.ones_like(resid, dtype=bool)
            resid_use = resid

        rms = float(np.sqrt(np.mean(resid_use ** 2)))
        if not np.isfinite(rms):
            note_fail("nonfinite_rms")
            continue

        if rms < best["rms_arcsec"]:
            best = {
                "ok": True,
                "hypo_r_au": float(r0),
                "hypo_rdot_au_day": float(rdot),
                "hypo_rdd_au_day2": float(rdd),
                "rms_arcsec": rms,
                "med_arcsec": float(np.median(resid_use)),
                "max_arcsec": float(np.max(resid_use)),
                "n_obs": int(len(jd)),
                "n_used": int(np.sum(keep)),
                "used_mask": keep,
                "resid_arcsec": resid,
                "t_epoch_jd": float(t1.jd),
                "orbit": orb,
                "fail_reason": "",
                "fail_counts": dict(fail_counts),
                "best_v1_kms": v1_kms,
                "best_rejected_max_v_kms": float(best["min_rejected_max_v_kms"]) if np.isfinite(best.get("min_rejected_max_v_kms", np.nan)) else np.nan,
            }

    best["fail_counts"] = dict(fail_counts)
    if not best.get("ok", False):
        best["best_rejected_max_v_kms"] = float(best["min_rejected_max_v_kms"]) if np.isfinite(best.get("min_rejected_max_v_kms", np.nan)) else np.nan
    return best


def compute_orbit_residuals_arcsec(
    orb: Orbit,
    mjd: np.ndarray,
    ra: np.ndarray,
    dec: np.ndarray,
    loc: EarthLocation,
) -> np.ndarray:
    t = Time(mjd, format="mjd", scale="utc")
    r_so = heliocentric_observer_xyz_au(t, loc)
    u_obs = eq_unitvec(ra, dec)
    dt_vec = t - orb.epoch
    r_obj = []
    for dti in dt_vec:
        oi = orb.propagate(dti)
        r_obj.append(oi.r.to(u.AU).value)
    r_obj = np.asarray(r_obj, dtype=float)
    rho = r_obj - r_so
    u_pred = rho / np.clip(norm(rho, axis=1, keepdims=True), 1e-30, None)
    return angular_sep_arcsec(u_obs, u_pred)


def build_tracklet_obs_masks(src_tid: np.ndarray) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    tids = np.asarray(src_tid).astype("U64")
    for tid in np.unique(tids):
        out[str(tid)] = tids == tid
    return out


def build_subset_tracklet_lists(tracklet_ids: np.ndarray) -> list[list[str]]:
    tids = [str(x) for x in np.asarray(tracklet_ids).astype("U64")]
    uniq = []
    seen = set()
    for tid in tids:
        if tid in seen:
            continue
        seen.add(tid)
        uniq.append(tid)

    subsets: list[list[str]] = []
    n = len(uniq)
    if n >= 3:
        subsets.extend([list(c) for c in itertools.combinations(uniq, 3)])
    subsets.extend([list(c) for c in itertools.combinations(uniq, 2)])
    if not subsets and uniq:
        subsets.append(uniq[:])
    return subsets[:24]


def fit_single_night_robust(
    tids: np.ndarray,
    mjd: np.ndarray,
    ra: np.ndarray,
    dec: np.ndarray,
    src_tid: np.ndarray,
    loc: EarthLocation,
    hypos: List[Tuple[float, float, float]],
    min_init_earth_au: float,
    max_v_kms: float,
    outlier_clip_arcsec: float | None,
    inlier_arcsec: float,
) -> Dict:
    best = fit_best_hypo_lambert(
        mjd=mjd,
        ra=ra,
        dec=dec,
        loc=loc,
        hypos=hypos,
        min_init_earth_au=min_init_earth_au,
        max_v_kms=max_v_kms,
        outlier_clip_arcsec=outlier_clip_arcsec,
    )

    tracklet_masks = build_tracklet_obs_masks(src_tid)
    subset_lists = build_subset_tracklet_lists(tids)
    candidate_fits: list[Dict] = []
    if best.get("ok", False):
        candidate_fits.append(best)

    for subset in subset_lists:
        subset_mask = np.zeros(len(mjd), dtype=bool)
        for tid in subset:
            subset_mask |= tracklet_masks.get(str(tid), np.zeros(len(mjd), dtype=bool))
        if int(np.sum(subset_mask)) < 3:
            continue
        fit = fit_best_hypo_lambert(
            mjd=mjd[subset_mask],
            ra=ra[subset_mask],
            dec=dec[subset_mask],
            loc=loc,
            hypos=hypos,
            min_init_earth_au=min_init_earth_au,
            max_v_kms=max_v_kms,
            outlier_clip_arcsec=outlier_clip_arcsec,
        )
        if fit.get("ok", False):
            fit["seed_tracklets"] = ";".join(subset)
            candidate_fits.append(fit)

    for seed_fit in candidate_fits:
        try:
            resid_all = compute_orbit_residuals_arcsec(seed_fit["orbit"], mjd, ra, dec, loc)
        except Exception:
            continue
        inlier_mask = resid_all <= float(inlier_arcsec)
        if int(np.sum(inlier_mask)) < 3:
            continue
        inlier_tracklets = {
            str(tid)
            for tid, keep in zip(np.asarray(src_tid).astype("U64"), inlier_mask)
            if keep
        }
        if len(inlier_tracklets) < 2:
            continue
        refit = fit_best_hypo_lambert(
            mjd=mjd[inlier_mask],
            ra=ra[inlier_mask],
            dec=dec[inlier_mask],
            loc=loc,
            hypos=hypos,
            min_init_earth_au=min_init_earth_au,
            max_v_kms=max_v_kms,
            outlier_clip_arcsec=outlier_clip_arcsec,
        )
        if not refit.get("ok", False):
            continue

        try:
            final_resid = compute_orbit_residuals_arcsec(refit["orbit"], mjd, ra, dec, loc)
        except Exception:
            continue
        final_mask = final_resid <= float(inlier_arcsec)
        if int(np.sum(final_mask)) < 3:
            continue
        refit["used_mask"] = final_mask
        refit["resid_arcsec"] = final_resid
        refit["n_obs"] = int(len(mjd))
        refit["n_used"] = int(np.sum(final_mask))
        refit["med_arcsec"] = float(np.median(final_resid[final_mask]))
        refit["max_arcsec"] = float(np.max(final_resid[final_mask]))
        refit["rms_arcsec"] = float(np.sqrt(np.mean(final_resid[final_mask] ** 2)))
        refit["seed_tracklets"] = seed_fit.get("seed_tracklets", "")

        if (not best.get("ok", False)) or (refit["rms_arcsec"] < best["rms_arcsec"]):
            best = refit

    return best


def format_fail_counts(fail_counts: Dict[str, int]) -> str:
    if not fail_counts:
        return ""
    parts = [f"{key}:{fail_counts[key]}" for key in sorted(fail_counts)]
    return ";".join(parts)


def orbit_elements_from_poliastro(orb: Orbit) -> Dict[str, float]:
    return {
        "a_au": float(orb.a.to_value(u.AU)),
        "ecc": float(orb.ecc.value),
        "inc_deg": float(orb.inc.to_value(u.deg)),
        "raan_deg": float(orb.raan.to_value(u.deg)),
        "argp_deg": float(orb.argp.to_value(u.deg)),
        "nu_deg": float(orb.nu.to_value(u.deg)),
    }


def predict_radec(orb: Orbit, t_pred: Time, loc: EarthLocation) -> Tuple[float, float]:
    """
    Predict topocentric RA/Dec by propagating 2-body orbit and subtracting heliocentric observer vector.

    Robust across poliastro/astropy versions:
      - observer computed in UTC (needed by EarthLocation.get_gcrs)
      - propagation dt computed in seconds, in TDB to match typical poliastro epoch scale
    """
    # Observer heliocentric position: use UTC for EarthLocation.get_gcrs stability
    t_pred_utc = t_pred.utc
    r_so = heliocentric_observer_xyz_au(t_pred_utc, loc)[0]

    # Propagation: use a Quantity time-of-flight (seconds), in TDB to match poliastro internals
    # poliastro Orbit.epoch is often in TDB
    dt = (t_pred.tdb - orb.epoch).to(u.s)  # Quantity

    oi = orb.propagate(dt)  # works with Quantity in modern poliastro
    r_obj = oi.r.to(u.AU).value

    rho = r_obj - r_so
    uvec = rho / np.clip(norm(rho), 1e-30, None)

    x, y, z = uvec
    ra = math.degrees(math.atan2(y, x)) % 360.0
    dec = math.degrees(math.asin(np.clip(z, -1.0, 1.0)))
    return ra, dec


# ---------------------- packed arrays ----------------------

def _pack_tracklets_arrays(tracklets: Table) -> dict:
    tid = np.asarray(tracklets["tracklet_id"]).astype("U64")

    out = {
        "tid": tid,
        "mjd1": np.asarray(tracklets["mjd1"], dtype=np.float64),
        "ra1": np.asarray(tracklets["ra1"], dtype=np.float64),
        "dec1": np.asarray(tracklets["dec1"], dtype=np.float64),
        "mjd2": np.asarray(tracklets["mjd2"], dtype=np.float64),
        "ra2": np.asarray(tracklets["ra2"], dtype=np.float64),
        "dec2": np.asarray(tracklets["dec2"], dtype=np.float64),
    }

    # detection identity (preferred)
    if "file1" in tracklets.colnames and "file2" in tracklets.colnames:
        out["file1"] = np.asarray(tracklets["file1"]).astype("U512")
        out["file2"] = np.asarray(tracklets["file2"]).astype("U512")
    else:
        out["file1"] = None
        out["file2"] = None

    if "objID1" in tracklets.colnames and "objID2" in tracklets.colnames:
        out["objID1"] = np.asarray(tracklets["objID1"], dtype=np.int64)
        out["objID2"] = np.asarray(tracklets["objID2"], dtype=np.int64)
    else:
        out["objID1"] = None
        out["objID2"] = None

    return out


def _expand_link_to_observations_fast(link_tids: np.ndarray, tarr: dict, idx: Dict[str, int]):
    """
    Expand link tracklets into per-detection observations and de-duplicate by detection key.

    Preferred de-dup key: (file, objID)
      - For consecutive tracklets AB and BC, the shared middle detection B will be deduped.
    Fallback (if no file/objID): round(mjd, ra, dec) and unique.
    """
    mjd_list, ra_list, dec_list = [], [], []
    file_list, obj_list = [], []
    obs_key_list, src_tid_list = [], []

    has_detkey = (tarr.get("file1") is not None) and (tarr.get("objID1") is not None)

    for tid in link_tids:
        j = idx.get(tid, None)
        if j is None:
            continue

        # endpoint 1
        mjd_list.append(float(tarr["mjd1"][j]))
        ra_list.append(float(tarr["ra1"][j]))
        dec_list.append(float(tarr["dec1"][j]))
        src_tid_list.append(str(tid))

        if has_detkey:
            f1 = str(tarr["file1"][j])
            o1 = int(tarr["objID1"][j])
            file_list.append(f1)
            obj_list.append(o1)
            obs_key_list.append(f"{Path(f1).name}:{o1}")
        else:
            file_list.append("")
            obj_list.append(-1)
            obs_key_list.append(f"{tid}:1")

        # endpoint 2
        mjd_list.append(float(tarr["mjd2"][j]))
        ra_list.append(float(tarr["ra2"][j]))
        dec_list.append(float(tarr["dec2"][j]))
        src_tid_list.append(str(tid))

        if has_detkey:
            f2 = str(tarr["file2"][j])
            o2 = int(tarr["objID2"][j])
            file_list.append(f2)
            obj_list.append(o2)
            obs_key_list.append(f"{Path(f2).name}:{o2}")
        else:
            file_list.append("")
            obj_list.append(-1)
            obs_key_list.append(f"{tid}:2")

    if not mjd_list:
        return None

    mjd = np.asarray(mjd_list, dtype=np.float64)
    ra = np.asarray(ra_list, dtype=np.float64)
    dec = np.asarray(dec_list, dtype=np.float64)
    obs_key = np.asarray(obs_key_list, dtype="U128")
    src_tid = np.asarray(src_tid_list, dtype="U64")

    # sort by time
    order = np.argsort(mjd, kind="mergesort")
    mjd, ra, dec, obs_key, src_tid = mjd[order], ra[order], dec[order], obs_key[order], src_tid[order]

    if has_detkey:
        file_arr = np.asarray(file_list, dtype="U512")[order]
        obj_arr = np.asarray(obj_list, dtype=np.int64)[order]

        # build structured key for fast unique
        key = np.empty(len(mjd), dtype=[("f", "U512"), ("o", "i8")])
        key["f"] = file_arr
        key["o"] = obj_arr

        _, uidx = np.unique(key, return_index=True)
        uidx = np.sort(uidx)
        return mjd[uidx], ra[uidx], dec[uidx], obs_key[uidx], src_tid[uidx]

    # fallback: unique by rounded mjd/ra/dec (coarse)
    key = np.empty(len(mjd), dtype=[("t", "i8"), ("ra", "i8"), ("dec", "i8")])
    key["t"] = np.round(mjd * 86400.0).astype(np.int64)  # 1-second bins
    key["ra"] = np.round(ra * 3600.0 * 100.0).astype(np.int64)   # 0.01 arcsec bins
    key["dec"] = np.round(dec * 3600.0 * 100.0).astype(np.int64) # 0.01 arcsec bins
    _, uidx = np.unique(key, return_index=True)
    uidx = np.sort(uidx)
    return mjd[uidx], ra[uidx], dec[uidx], obs_key[uidx], src_tid[uidx]


# ---------------------- worker globals + initializer ----------------------

G_IDX: Optional[Dict[str, int]] = None
G_TARR: Optional[dict] = None
G_LOC: Optional[EarthLocation] = None
G_HYPOS: Optional[List[Tuple[float, float, float]]] = None
G_MIN_INIT: float = 0.0
G_MAX_V: float = 200.0
G_OUTLIER_CLIP: Optional[float] = None
G_SINGLE_THR: float = 3.0
G_MULTI_THR: float = 1.2
G_SINGLE_MAX: float = 6.0
G_MULTI_MAX: float = 3.0
G_MIN_USED_FRAC: float = 0.9
G_ECC_MAX_MULTI: float = 1.05
G_PRED_DAYS: float = 1.0


def _init_worker(
    idx: Dict[str, int],
    tarr: dict,
    loc: EarthLocation,
    hypos: List[Tuple[float, float, float]],
    min_init: float,
    max_v: float,
    outlier_clip: Optional[float],
    single_thr: float,
    multi_thr: float,
    single_max: float,
    multi_max: float,
    min_used_frac: float,
    ecc_max_multi: float,
    pred_days: float,
):
    global G_IDX, G_TARR, G_LOC, G_HYPOS
    global G_MIN_INIT, G_MAX_V, G_OUTLIER_CLIP
    global G_SINGLE_THR, G_MULTI_THR, G_SINGLE_MAX, G_MULTI_MAX
    global G_MIN_USED_FRAC, G_ECC_MAX_MULTI, G_PRED_DAYS

    G_IDX = idx
    G_TARR = tarr
    G_LOC = loc
    G_HYPOS = hypos

    G_MIN_INIT = float(min_init)
    G_MAX_V = float(max_v)
    G_OUTLIER_CLIP = outlier_clip

    G_SINGLE_THR = float(single_thr)
    G_MULTI_THR = float(multi_thr)
    G_SINGLE_MAX = float(single_max)
    G_MULTI_MAX = float(multi_max)
    G_MIN_USED_FRAC = float(min_used_frac)
    G_ECC_MAX_MULTI = float(ecc_max_multi)
    G_PRED_DAYS = float(pred_days)


def _process_one_link_small(task: Tuple[int, np.ndarray]) -> Tuple[Optional[tuple], List[tuple]]:
    """
    Task only sends (lid, tids). Everything else comes from worker globals.
    Returns: (summary_row or None, resid_rows list)
    """
    lid, tids = task
    assert G_IDX is not None and G_TARR is not None and G_LOC is not None and G_HYPOS is not None

    obs = _expand_link_to_observations_fast(tids, G_TARR, G_IDX)
    if obs is None:
        return None, []

    mjd, ra, dec, obs_key, src_tid = obs
    lin_rms, lin_speed, lin_dir = linear_motion_test(mjd, ra, dec)
    n_nights = count_nights_from_mjd(mjd)
    thr_rms = G_SINGLE_THR if n_nights <= 1 else G_MULTI_THR

    if n_nights <= 1:
        fit = fit_single_night_robust(
            tids=tids,
            mjd=mjd,
            ra=ra,
            dec=dec,
            src_tid=src_tid,
            loc=G_LOC,
            hypos=G_HYPOS,
            min_init_earth_au=G_MIN_INIT,
            max_v_kms=G_MAX_V,
            outlier_clip_arcsec=G_OUTLIER_CLIP,
            inlier_arcsec=min(G_SINGLE_MAX, max(G_SINGLE_THR * 1.5, 3.0)),
        )
    else:
        fit = fit_best_hypo_lambert(
            mjd=mjd,
            ra=ra,
            dec=dec,
            loc=G_LOC,
            hypos=G_HYPOS,
            min_init_earth_au=G_MIN_INIT,
            max_v_kms=G_MAX_V,
            outlier_clip_arcsec=G_OUTLIER_CLIP,
        )

    if not fit.get("ok", False):
        summary_row = (
            int(lid), int(len(tids)), int(n_nights), int(len(mjd)),
            False,
            np.nan, np.nan, np.nan,
            np.inf, np.nan, np.nan,
            False,
            np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
            np.nan, np.nan,
            float(lin_rms), float(lin_speed), float(lin_dir),
            float(fit.get("best_v1_kms", np.nan)),
            float(fit.get("min_rejected_max_v_kms", np.nan)),
            float(fit.get("best_rejected_max_v_kms", np.nan)),
            str(fit.get("fail_reason", "no_valid_hypo")),
            format_fail_counts(fit.get("fail_counts", {})),
        )
        return summary_row, []

    orb: Orbit = fit["orbit"]
    elems = orbit_elements_from_poliastro(orb)

    t_last = Time(float(mjd[-1] + G_PRED_DAYS), format="mjd", scale="utc")
    try:
        ra_pred, dec_pred = predict_radec(orb, t_last, G_LOC)
    except Exception as e:
        # print("[pred error]", repr(e))
        ra_pred, dec_pred = (np.nan, np.nan)

    # quality gating
    if n_nights <= 1:
        lin_ok = lin_rms <= 10
        is_good = (
            np.isfinite(fit["rms_arcsec"]) and (fit["rms_arcsec"] <= thr_rms) and
            np.isfinite(fit["max_arcsec"]) and (fit["max_arcsec"] <= G_SINGLE_MAX) and
            (fit["n_used"] >= 3) and
            (len(tids) >= 2) and
            lin_ok
        )
    else:
        used_frac = fit["n_used"] / max(1, fit["n_obs"])
        ecc_ok = np.isfinite(elems["ecc"]) and (elems["ecc"] <= G_ECC_MAX_MULTI)
        a_ok = np.isfinite(elems["a_au"]) and (0.5 <= elems["a_au"] <= 50.0)
        lin_ok = lin_rms <= 10
        is_good = (
            np.isfinite(fit["rms_arcsec"]) and (fit["rms_arcsec"] <= thr_rms) and
            np.isfinite(fit["max_arcsec"]) and (fit["max_arcsec"] <= G_MULTI_MAX) and
            (used_frac >= G_MIN_USED_FRAC) and
            ecc_ok and a_ok and lin_ok
        )

    summary_row = (
        int(lid), int(len(tids)), int(n_nights), int(fit["n_obs"]),
        True,
        float(fit["hypo_r_au"]),
        float(fit["hypo_rdot_au_day"]),
        float(fit["hypo_rdd_au_day2"]),
        float(fit["rms_arcsec"]),
        float(fit["med_arcsec"]),
        float(fit["max_arcsec"]),
        bool(is_good),
        float(elems["a_au"]),
        float(elems["ecc"]),
        float(elems["inc_deg"]),
        float(elems["raan_deg"]),
        float(elems["argp_deg"]),
        float(elems["nu_deg"]),
        float(ra_pred),
        float(dec_pred),
        float(lin_rms),
        float(lin_speed),
        float(lin_dir),
        float(fit.get("best_v1_kms", np.nan)),
        float(fit.get("min_rejected_max_v_kms", np.nan)),
        float(fit.get("best_rejected_max_v_kms", np.nan)),
        "",
        format_fail_counts(fit.get("fail_counts", {})),
    )

    resid = np.asarray(fit["resid_arcsec"], dtype=np.float64)
    used = np.asarray(fit["used_mask"], dtype=bool)
    resid_rows = []
    for i in range(len(mjd)):
        resid_rows.append(
            (
                int(lid),
                str(obs_key[i]),
                str(src_tid[i]),
                float(mjd[i]),
                float(ra[i]),
                float(dec[i]),
                float(resid[i]),
                bool(used[i]),
            )
        )

    return summary_row, resid_rows


def _choose_chunksize(n_links: int, cores: int) -> int:
    """
    Heuristic chunksize: avoid too many tiny tasks; keep good load balance.
    For 50k links and 8 cores, ~100-300 is usually fine.
    """
    if cores <= 1:
        return 1
    # target ~ 1000-5000 chunks total
    target_chunks = min(max(cores * 400, 1000), 5000)
    cs = max(10, int(math.ceil(n_links / target_chunks)))
    # cap to avoid too large chunks hurting balance
    return min(cs, 500)


# ---------------------- main ----------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rr-dir", required=True)
    ap.add_argument("--profile", choices=sorted(PROFILE_DEFAULTS.keys()), default=None)
    ap.add_argument("--root", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--tracklets", default=None)
    ap.add_argument("--outdir", default=None)

    ap.add_argument("--min-init-earth-au", type=float, default=None)
    ap.add_argument("--max-v-kms", type=float, default=None)
    ap.add_argument("--hypos", type=str, default=None)

    ap.add_argument("--single-rms-arcsec", type=float, default=5.0)
    ap.add_argument("--multi-rms-arcsec", type=float, default=3.0)
    ap.add_argument("--single-max-arcsec", type=float, default=8.0)
    ap.add_argument("--multi-max-arcsec", type=float, default=5.0)
    ap.add_argument("--min-used-frac", type=float, default=0.8)
    ap.add_argument("--ecc-max-multi", type=float, default=1.1)
    ap.add_argument("--outlier-clip-arcsec", type=float, default=5.0)
    ap.add_argument("--limit-links", type=int, default=0)

    ap.add_argument("--pred-days", type=float, default=1.0)

    ap.add_argument("--cores", type=int, default=1)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--log-every", type=int, default=0)

    ap.add_argument("--set-threads-1", action="store_true", default=True,
                    help="set OMP/MKL/OPENBLAS/NUMEXPR threads to 1 (recommended for multiprocessing)")
    ap.add_argument("--no-set-threads-1", dest="set_threads_1", action="store_false")

    args = ap.parse_args()

    # avoid BLAS oversubscription
    if args.set_threads_1:
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        os.environ.setdefault("MKL_NUM_THREADS", "1")
        os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
        os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    rr_dir = Path(args.rr_dir).expanduser()
    profile = resolve_profile(rr_dir, args.profile)
    root = Path(args.root).expanduser()
    outdir = Path(args.outdir).expanduser() if args.outdir else (rr_dir / "orbit_confirm")
    outdir.mkdir(parents=True, exist_ok=True)

    def log(msg: str):
        if not args.quiet:
            print(msg)

    mode = detect_mode_from_rr_dir(rr_dir)
    night = extract_night_from_rr_dir(rr_dir)
    tracklets_path = resolve_tracklets_path(rr_dir, root, args.tracklets)

    log(f"[info] rr_dir={rr_dir}")
    log(f"[info] mode={mode}  night={night}")
    log(f"[info] tracklets={tracklets_path}")

    links, members = load_rr_links(rr_dir)

    tracklets = Table.read(tracklets_path)
    idx = build_tracklet_index(tracklets)
    tarr = _pack_tracklets_arrays(tracklets)

    loc = earth_location_mpc327()

    # hypos
    if args.hypos:
        hypos: List[Tuple[float, float, float]] = []
        with open(args.hypos, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if (not s) or s.startswith("#"):
                    continue
                r, rd, rdd = [float(x) for x in s.split()[:3]]
                hypos.append((r, rd, rdd))
    else:
        hypos = build_hypos_default(profile)

    if args.min_init_earth_au is None:
        min_init = float(PROFILE_DEFAULTS[profile]["min_init_earth_au"])
    else:
        min_init = float(args.min_init_earth_au)

    max_v_kms = float(PROFILE_DEFAULTS[profile]["max_v_kms"]) if args.max_v_kms is None else float(args.max_v_kms)

    outlier_clip = None if args.outlier_clip_arcsec <= 0 else float(args.outlier_clip_arcsec)

    # linkage_id -> list(tracklet_id)
    mem: Dict[int, List[str]] = {}
    for lid, tid in zip(np.asarray(members["linkage_id"], dtype=int),
                        np.asarray(members["tracklet_id"]).astype("U64")):
        mem.setdefault(int(lid), []).append(str(tid))

    lids = np.asarray(links["linkage_id"], dtype=int)
    if args.limit_links and args.limit_links > 0:
        lids = lids[: int(args.limit_links)]

    tasks: List[Tuple[int, np.ndarray]] = []
    for lid in lids:
        tids = np.array(mem.get(int(lid), []), dtype="U64")
        tasks.append((int(lid), tids))

    summary_rows: List[tuple] = []
    resid_rows: List[tuple] = []

    n_links = len(tasks)
    if args.cores <= 1:
        # serial: use same initializer to reuse logic
        _init_worker(
            idx=idx,
            tarr=tarr,
            loc=loc,
            hypos=hypos,
            min_init=min_init,
            max_v=max_v_kms,
            outlier_clip=outlier_clip,
            single_thr=float(args.single_rms_arcsec),
            multi_thr=float(args.multi_rms_arcsec),
            single_max=float(args.single_max_arcsec),
            multi_max=float(args.multi_max_arcsec),
            min_used_frac=float(args.min_used_frac),
            ecc_max_multi=float(args.ecc_max_multi),
            pred_days=float(args.pred_days),
        )

        for k, task in enumerate(tasks, start=1):
            srow, rrows = _process_one_link_small(task)
            if srow is not None:
                summary_rows.append(srow)
            if rrows:
                resid_rows.extend(rrows)

            if (not args.quiet) and args.log_every and (k % args.log_every == 0):
                print(f"[info] processed {k}/{n_links} links")

    else:
        nproc = int(args.cores)
        chunksize = _choose_chunksize(n_links, nproc)
        log(f"[info] multiprocessing: cores={nproc}  chunksize={chunksize}")

        with ProcessPoolExecutor(
            max_workers=nproc,
            initializer=_init_worker,
            initargs=(
                idx,
                tarr,
                loc,
                hypos,
                float(min_init),
                max_v_kms,
                outlier_clip,
                float(args.single_rms_arcsec),
                float(args.multi_rms_arcsec),
                float(args.single_max_arcsec),
                float(args.multi_max_arcsec),
                float(args.min_used_frac),
                float(args.ecc_max_multi),
                float(args.pred_days),
            ),
        ) as ex:
            # map returns results in input order (good for reproducibility)
            for k, (srow, rrows) in enumerate(ex.map(_process_one_link_small, tasks, chunksize=chunksize), start=1):
                if srow is not None:
                    summary_rows.append(srow)
                if rrows:
                    resid_rows.extend(rrows)

                if (not args.quiet) and args.log_every and (k % args.log_every == 0):
                    print(f"[info] processed {k}/{n_links} links")

    # sort
    if summary_rows:
        summary_rows.sort(key=lambda r: int(r[0]))
    if resid_rows:
        resid_rows.sort(key=lambda r: (int(r[0]), float(r[3])))

    summary = Table(
        rows=summary_rows,
        names=[
            "linkage_id",
            "n_tracklets",
            "n_nights",
            "n_obs",
            "fit_ok",
            "hypo_r_au",
            "hypo_rdot_au_day",
            "hypo_rdd_au_day2",
            "rms_arcsec",
            "med_arcsec",
            "max_arcsec",
            "is_good",
            "a_au",
            "ecc",
            "inc_deg",
            "raan_deg",
            "argp_deg",
            "nu_deg",
            "pred_ra_deg",
            "pred_dec_deg",
            "lin_rms_arcsec",
            "lin_speed_arcsec_per_day",
            "lin_dir_deg",
            "best_v1_kms",
            "min_rejected_max_v_kms",
            "best_rejected_max_v_kms",
            "fail_reason",
            "fail_counts",
        ],
    )

    summary.meta["rr_dir"] = rr_dir.name
    summary.meta["tracklets"] = tracklets_path.name
    summary.meta["mode"] = mode
    summary.meta["profile"] = profile
    summary.meta["night"] = night
    summary.meta["n_hypos"] = int(len(hypos))
    summary.meta["min_init_earth_au"] = float(min_init)
    summary.meta["max_v_kms"] = float(max_v_kms)
    summary.meta["single_rms_arcsec"] = float(args.single_rms_arcsec)
    summary.meta["multi_rms_arcsec"] = float(args.multi_rms_arcsec)
    summary.meta["pred_days"] = float(args.pred_days)
    summary.meta["cores"] = int(args.cores)

    out_links = outdir / "orbit_links.fits"
    summary.write(out_links, overwrite=True)
    log(f"[write] {out_links}  n_links={len(summary)}")

    if resid_rows:
        resid_tab = Table(
            rows=resid_rows,
            names=[
                "linkage_id",
                "obs_key",
                "tracklet_id",
                "mjd",
                "ra_deg",
                "dec_deg",
                "resid_arcsec",
                "used",
            ],
        )
        out_resid = outdir / "orbit_obs_residuals.fits"
        resid_tab.write(out_resid, overwrite=True)
        log(f"[write] {out_resid}  n_rows={len(resid_tab)}")

    try:
        (outdir / "provenance.txt").write_text(
            f"rr_dir: {rr_dir}\ntracklets: {tracklets_path}\nmode: {mode}\nnight: {night}\n",
            encoding="utf-8",
        )
    except Exception:
        pass

    if len(summary) > 0:
        n_ok = int(np.sum(np.asarray(summary["fit_ok"], dtype=bool)))
        n_good = int(np.sum(np.asarray(summary["is_good"], dtype=bool)))
        log(f"[done] fit_ok={n_ok}/{len(summary)}  is_good={n_good}/{len(summary)}")
        fail_mask = ~np.asarray(summary["fit_ok"], dtype=bool)
        if np.any(fail_mask):
            reasons = np.asarray(summary["fail_reason"]).astype("U64")[fail_mask]
            uniq, cnt = np.unique(reasons, return_counts=True)
            pairs = sorted(zip(cnt, uniq), reverse=True)
            joined = ", ".join(f"{reason}={int(n)}" for n, reason in pairs)
            log(f"[done] fit_fail_reasons: {joined}")


if __name__ == "__main__":
    main()
