#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path
from itertools import chain

import numpy as np
from numpy import linalg as LA
from astropy.table import Table
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import EarthLocation, get_body_barycentric, CartesianRepresentation
from sklearn.neighbors import KDTree

from poliastro.iod import izzo
from poliastro.bodies import Sun
from poliastro.twobody import Orbit

import warnings
from astropy.io.fits.verify import VerifyWarning

warnings.filterwarnings(
    "ignore",
    category=VerifyWarning
)

try:
    from poliastro.twobody.propagation.farnocchia import FarnocchiaPropagator
    POLINEW = True
except Exception:
    from poliastro.twobody.propagation import propagate
    POLINEW = False


# MPC 327 default (you can override in CLI later if you want)
MPC327_LON_DEG = 117.5750
MPC327_LAT_DEG = 40.394239
MPC327_ALT_M = 868.221


def earth_location_mpc327() -> EarthLocation:
    return EarthLocation.from_geodetic(
        lon=MPC327_LON_DEG * u.deg,
        lat=MPC327_LAT_DEG * u.deg,
        height=MPC327_ALT_M * u.m,
    )


def heliocentric_observer_xyz_au(times_utc: Time, loc: EarthLocation) -> np.ndarray:
    earth_bary = get_body_barycentric("earth", times_utc)
    sun_bary = get_body_barycentric("sun", times_utc)
    topo = loc.get_gcrs(times_utc).cartesian
    topo_au = CartesianRepresentation(topo.x.to(u.au), topo.y.to(u.au), topo.z.to(u.au))
    obs_bary = earth_bary + topo_au
    obs_helio = obs_bary - sun_bary
    return np.stack(
        [obs_helio.x.to_value(u.au), obs_helio.y.to_value(u.au), obs_helio.z.to_value(u.au)],
        axis=1,
    ).astype(np.float64)


def propagate_positions_compat(orb: Orbit, timedelta_vector) -> np.ndarray:
    if POLINEW:
        xyz, _ = FarnocchiaPropagator.propagate_many(orb, orb._state, timedelta_vector)
        return np.asarray(np.transpose(xyz.to(u.AU).value), dtype=np.float64)

    rows = []
    for dt in timedelta_vector:
        propagated = propagate(orb, dt, rtol=1e-4)
        if hasattr(propagated, "xyz"):
            xyz = propagated.xyz.to(u.AU).value
        else:
            xyz = propagated.r.to(u.AU).value
        rows.append(np.asarray(xyz, dtype=np.float64))
    return np.asarray(rows, dtype=np.float64)


def _fmt_elapsed(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, sec = divmod(int(seconds), 60)
    if minutes < 60:
        return f"{minutes}m{sec:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m{sec:02d}s"


class ProgressLogger:
    def __init__(self, outdir: Path):
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.log_path = self.outdir / "progress.log"

    def log(self, message: str) -> None:
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}"
        print(line, flush=True)
        with open(self.log_path, "a", encoding="utf-8") as f:
            f.write(line + "\n")


class HeliolincRR_Tracklets:
    def __init__(self, ra1, dec1, ra2, dec2, jd1, jd2, r_so_1, r_so_2, tracklet_ids):
        self.tracklet_ids = np.asarray(tracklet_ids).astype("U64")
        n = len(ra1)

        # Build 2N observations in original order
        ra = np.empty(2 * n, dtype=np.float64)
        dec = np.empty(2 * n, dtype=np.float64)
        jdates = np.empty(2 * n, dtype=np.float64)
        r_so = np.empty((2 * n, 3), dtype=np.float64)

        ra[0::2] = ra1
        dec[0::2] = dec1
        jdates[0::2] = jd1
        r_so[0::2, :] = r_so_1

        ra[1::2] = ra2
        dec[1::2] = dec2
        jdates[1::2] = jd2
        r_so[1::2, :] = r_so_2

        # Sort chronologically
        order = np.argsort(jdates, kind="mergesort")
        ra = ra[order]
        dec = dec[order]
        jdates = jdates[order]
        r_so = r_so[order]

        inv = np.empty_like(order)
        inv[order] = np.arange(len(order), dtype=np.int64)

        # combs in sorted obs indices
        self.combs = []
        for i in range(n):
            a0 = int(inv[2 * i])
            a1 = int(inv[2 * i + 1])
            if jdates[a0] <= jdates[a1]:
                self.combs.append([a0, a1])
            else:
                self.combs.append([a1, a0])

        # unique times and fnums
        utimes, inv_t = np.unique(jdates, return_inverse=True)
        self.utimes = utimes.astype(np.float64)

        first_idx = np.full(len(utimes), -1, dtype=np.int64)
        for i_obs, fnum in enumerate(inv_t):
            if first_idx[fnum] < 0:
                first_idx[fnum] = i_obs
        self.r_so = r_so[first_idx, :].astype(np.float64)

        # bbf: [ra, dec, orig_obs_id, fnum]
        orig_obs_id_sorted = order.astype(np.int64)
        self.bbf = np.empty((len(jdates), 4), dtype=object)
        self.bbf[:, 0] = ra
        self.bbf[:, 1] = dec
        self.bbf[:, 2] = orig_obs_id_sorted
        self.bbf[:, 3] = inv_t.astype(np.int64)

        self.r_ob_unit = self.eq2cart(ra, dec, 1.0)[0].astype(np.float64)
        self.times = Time(self.utimes, format="jd")
        self.dts = np.asarray([(t - self.times[0]).value for t in self.times], dtype=np.float64)
        self.fnums = np.asarray(self.bbf[:, 3], dtype=np.int64)
        self.combs_arr = np.asarray(self.combs, dtype=np.int64)
        self.comb_f0 = self.fnums[self.combs_arr[:, 0]]
        self.comb_f1 = self.fnums[self.combs_arr[:, 1]]
        self.comb_r_ob0 = self.r_ob_unit[self.combs_arr[:, 0], :]
        self.comb_r_ob1 = self.r_ob_unit[self.combs_arr[:, 1], :]
        self.comb_r_so0 = self.r_so[self.comb_f0, :]
        self.comb_r_so1 = self.r_so[self.comb_f1, :]

    @staticmethod
    def eq2cart(ra, dec, d=1.0):
        ra_rad = np.deg2rad(ra)
        dec_rad = np.deg2rad(dec)
        ooe = 0.40909262968940374
        rot_ooe = np.matrix([[1, 0, 0], [0, np.cos(ooe), np.sin(ooe)], [0, -np.sin(ooe), np.cos(ooe)]])
        x = d * np.cos(ra_rad) * np.cos(dec_rad)
        y = d * np.sin(ra_rad) * np.cos(dec_rad)
        z = d * np.sin(dec_rad)
        eq_cart = np.column_stack((x, y, z))
        ec_cart = np.asarray((rot_ooe * np.matrix(eq_cart).T).T)
        return eq_cart, ec_cart

    @staticmethod
    def geodist_from_heliodist(r_so, r_ob, d):
        a = np.sum(r_ob**2)
        b = 2 * np.sum(r_ob * r_so)
        c = -1 * (d**2 - np.sum(r_so**2))
        r = np.roots([a, b, c])
        return r[r > 0][0]

    def propagator(self, hypo, epochs_jd, max_v_kms, min_init_earth_au):
        r, rdot, rdotdot = hypo

        epochs = Time(epochs_jd, format="jd")

        prvs = []
        combs_used = []

        for idx, c in enumerate(self.combs):
            f0 = self.comb_f0[idx]
            f1 = self.comb_f1[idx]
            t0 = self.times[f0]
            t1 = self.times[f1]

            r_ob0 = self.comb_r_ob0[idx]
            r_ob1 = self.comb_r_ob1[idx]
            r_so0 = self.comb_r_so0[idx]
            r_so1 = self.comb_r_so1[idx]

            range_dt0 = r + rdot * self.dts[f0] + 0.5 * rdotdot * self.dts[f0] ** 2
            range_dt1 = r + rdot * self.dts[f1] + 0.5 * rdotdot * self.dts[f1] ** 2

            try:
                alpha0 = self.geodist_from_heliodist(r_so0, r_ob0, range_dt0)
                alpha1 = self.geodist_from_heliodist(r_so1, r_ob1, range_dt1)
            except Exception:
                continue

            if (alpha0 <= min_init_earth_au) or (alpha1 <= min_init_earth_au):
                continue

            hr0 = (r_so0 + alpha0 * r_ob0) * u.AU
            hr1 = (r_so1 + alpha1 * r_ob1) * u.AU

            try:
                velocities = izzo.lambert(Sun.k, hr0, hr1, t1 - t0, rtol=1e-4)
                if POLINEW:
                    _, v1 = velocities
                else:
                    (_, v1), = velocities
            except Exception:
                continue

            if LA.norm(v1.value) >= max_v_kms:
                continue

            try:
                orb = Orbit.from_vectors(Sun, hr1, v1)
                timedelta_vector = epochs - t1
                prv = np.reshape(
                    propagate_positions_compat(orb, timedelta_vector),
                    (1, 6),
                )[0].tolist()
            except Exception:
                continue

            if not np.any(np.isnan(prv)):
                prvs.append(prv)
                combs_used.append(c)

        return np.asarray(prvs, dtype=np.float64), combs_used, hypo

    def propagate_all(self, hypos, epochs_jd, cores, max_v_kms, min_init_earth_au, logger: ProgressLogger | None = None):
        hypos = list(hypos)

        results = []
        t0 = time.perf_counter()
        n_hypos = len(hypos)
        if cores <= 1:
            for i, hypo in enumerate(hypos, start=1):
                prvs, combs_used, _ = self.propagator(hypo, epochs_jd, max_v_kms, min_init_earth_au)
                results.append((prvs, combs_used, hypo))
                if logger is not None:
                    logger.log(
                        f"[propagate] {i}/{n_hypos} hypo={hypo} n_prvs={len(prvs)} elapsed={_fmt_elapsed(time.perf_counter() - t0)}"
                    )
        else:
            with mp.Pool(cores) as p:
                jobs = [p.apply_async(self.propagator, args=(h, epochs_jd, max_v_kms, min_init_earth_au)) for h in hypos]
                for i, (hypo, job) in enumerate(zip(hypos, jobs), start=1):
                    prvs, combs_used, _ = job.get()
                    results.append((prvs, combs_used, hypo))
                    if logger is not None:
                        logger.log(
                            f"[propagate] {i}/{n_hypos} hypo={hypo} n_prvs={len(prvs)} elapsed={_fmt_elapsed(time.perf_counter() - t0)}"
                        )
        if logger is not None:
            logger.log(f"[propagate] done n_hypos={n_hypos} elapsed={_fmt_elapsed(time.perf_counter() - t0)}")
        return results

    @staticmethod
    def _min_nights_filter(links_obs, bbf, utimes, min_nights):
        accepted = []
        for link in links_obs:
            fnums = np.asarray([int(bbf[i, 3]) for i in link], dtype=np.int64)
            ltimes = np.unique(utimes[fnums])
            nights = int(np.sum(ltimes[1:] - ltimes[:-1] > 0.5) + 1)
            if nights >= min_nights:
                accepted.append(link)
        return accepted

    @staticmethod
    def _unique_links(links_obs):
        return [list(l) for l in sorted(set(tuple(l) for l in links_obs))]

    def cluster_one(self, prvs, combs_used, hypo, tol, min_len_obs, min_nights, k_neighbors_cap):
        if prvs.size == 0:
            return []

        drvs = prvs.copy()
        drvs[:, 0:3] /= hypo[0]
        drvs[:, 3:6] /= hypo[0]

        tree = KDTree(drvs)
        num_tracklets = drvs.shape[0]
        k_query = min(int(k_neighbors_cap), num_tracklets)
        fnums = np.asarray(self.bbf[:, 3], dtype=np.int64)

        comb_obs = [tuple(int(x) for x in comb) for comb in combs_used]
        comb_fnums = [tuple(int(fnums[idx]) for idx in comb) for comb in combs_used]

        edges = []
        for i, rv in enumerate(drvs):
            dist, kidx = tree.query(np.asarray([rv]), k=k_query)
            for d, j in zip(dist[0], kidx[0]):
                j = int(j)
                if j <= i or d > tol:
                    continue
                edges.append((float(d), i, j))
        edges.sort(key=lambda x: (x[0], x[1], x[2]))

        parent = list(range(num_tracklets))
        rank = [0] * num_tracklets
        comp_obs = [set(obs) for obs in comb_obs]
        comp_fnums = [set(fs) for fs in comb_fnums]

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra = find(a)
            rb = find(b)
            if ra == rb:
                return
            if comp_fnums[ra] & comp_fnums[rb]:
                return
            if rank[ra] < rank[rb]:
                ra, rb = rb, ra
            parent[rb] = ra
            comp_obs[ra].update(comp_obs[rb])
            comp_fnums[ra].update(comp_fnums[rb])
            if rank[ra] == rank[rb]:
                rank[ra] += 1

        for _, i, j in edges:
            union(i, j)

        links = []
        seen = set()
        for i in range(num_tracklets):
            root = find(i)
            if root in seen:
                continue
            seen.add(root)
            link = sorted(comp_obs[root])
            if len(link) >= min_len_obs:
                links.append(link)

        links = self._min_nights_filter(links, self.bbf, self.utimes, min_nights)
        return self._unique_links(links)

    def cluster_all(self, results, tol, min_len_obs, min_nights, k_neighbors_cap, logger: ProgressLogger | None = None):
        all_links = []
        t0 = time.perf_counter()
        n_hypos = len(results)
        for i, (prvs, combs_used, hypo) in enumerate(results, start=1):
            links = self.cluster_one(prvs, combs_used, hypo, tol, min_len_obs, min_nights, k_neighbors_cap)
            all_links.extend(links)
            if logger is not None:
                logger.log(
                    f"[cluster] {i}/{n_hypos} hypo={hypo} n_prvs={len(prvs)} n_links_added={len(links)} elapsed={_fmt_elapsed(time.perf_counter() - t0)}"
                )
        if logger is not None:
            logger.log(f"[cluster] done n_hypos={n_hypos} elapsed={_fmt_elapsed(time.perf_counter() - t0)}")
        return self._unique_links(all_links)


def build_hypos_default():
    rs = np.array([1.3, 1.7, 2.1, 2.5, 2.9, 3.3, 3.7, 4.1], dtype=np.float64)
    rdots = np.array([-0.02, 0.0, +0.02], dtype=np.float64)
    rdd = np.array([0.0], dtype=np.float64)
    return [(float(r), float(rd), float(rdd0)) for r in rs for rd in rdots for rdd0 in rdd]


def choose_reference_epochs(jd1, jd2, mode: str, ref_dt_days: float):
    all_jd = np.concatenate([jd1, jd2])
    if mode == "mid":
        jd_ref0 = float(np.median(all_jd))
    elif mode == "start":
        jd_ref0 = float(np.min(all_jd))
    else:
        raise ValueError(mode)
    return np.array([jd_ref0, jd_ref0 + float(ref_dt_days)], dtype=np.float64)


def summarize_int_column(arr) -> dict[str, object]:
    arr = np.asarray(arr, dtype=np.int64)
    if arr.size == 0:
        return {
            "min": 0,
            "median": 0.0,
            "max": 0,
            "counts": {},
        }
    vals, cnts = np.unique(arr, return_counts=True)
    return {
        "min": int(arr.min()),
        "median": float(np.median(arr)),
        "max": int(arr.max()),
        "counts": {str(int(v)): int(c) for v, c in zip(vals, cnts)},
    }


def write_rr_summary_json(
    outdir: Path,
    infile: Path,
    args,
    hypos,
    epochs,
    results,
    links_tab: Table,
    members_tab: Table,
    elapsed_propagation_s: float,
    elapsed_clustering_s: float,
):
    summary = {
        "infile": str(infile),
        "outdir": str(outdir),
        "params": {
            "cores": int(args.cores),
            "max_v_kms": float(args.max_v_kms),
            "min_init_earth_au": float(args.min_init_earth_au),
            "ref_epoch_mode": str(args.ref_epoch_mode),
            "ref_dt_days": float(args.ref_dt_days),
            "tol": float(args.tol),
            "min_len_obs": int(args.min_len_obs),
            "min_nights": int(args.min_nights),
            "k_neighbors_cap": int(args.k_neighbors_cap),
            "hypos_path": str(args.hypos) if args.hypos else None,
        },
        "reference_epochs_jd": [float(x) for x in epochs],
        "n_hypos": int(len(hypos)),
        "propagation_total_elapsed_s": float(elapsed_propagation_s),
        "clustering_total_elapsed_s": float(elapsed_clustering_s),
        "per_hypo": [
            {
                "hypo": {
                    "r": float(hypo[0]),
                    "rdot": float(hypo[1]),
                    "rddot": float(hypo[2]),
                },
                "n_prvs": int(len(prvs)),
            }
            for prvs, _, hypo in results
        ],
        "outputs": {
            "n_links": int(len(links_tab)),
            "n_member_rows": int(len(members_tab)),
            "n_tracklets_per_link": summarize_int_column(links_tab["n_tracklets"] if len(links_tab) else []),
            "n_nights_per_link": summarize_int_column(links_tab["n_nights"] if len(links_tab) else []),
        },
    }
    out_path = outdir / "rr_summary.json"
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out_path


def write_links_fits_with_endpoints(
    outdir: Path,
    links_obs,
    bbf,
    utimes,
    tracklets_tab: Table,
):
    """
    Writes:
      - links_tracklets.fits (with endpoints ra/dec/mjd columns)
      - linkage_members.fits
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # obs -> tracklet index
    orig_obs_id = bbf[:, 2].astype(np.int64)
    obs_to_tracklet = orig_obs_id // 2

    tracklet_ids = np.asarray(tracklets_tab["tracklet_id"]).astype("U64")
    trk_index = {tid: i for i, tid in enumerate(tracklet_ids)}

    link_rows = []
    mem_linkage = []
    mem_tracklet = []

    # endpoint columns buffers
    n_links = len(links_obs)
    def nancol():
        return np.full(n_links, np.nan, dtype=np.float64)

    ra1 = nancol(); dec1 = nancol(); ra2 = nancol(); dec2 = nancol()
    ra3 = nancol(); dec3 = nancol(); ra4 = nancol(); dec4 = nancol()
    mjd1 = nancol(); mjd2 = nancol(); mjd3 = nancol(); mjd4 = nancol()

    for lid, obs_link in enumerate(links_obs):
        obs_link = np.asarray(obs_link, dtype=np.int64)

        tids = np.unique(tracklet_ids[obs_to_tracklet[obs_link]])

        fnums = bbf[obs_link, 3].astype(np.int64)
        ltimes = np.unique(utimes[fnums])
        nights = int(np.sum(ltimes[1:] - ltimes[:-1] > 0.5) + 1)

        link_rows.append((lid, int(len(tids)), int(nights)))

        mem_linkage.extend([lid] * len(tids))
        mem_tracklet.extend(tids.tolist())

        # ---- endpoints (first/last tracklet by start time) ----
        idxs = [trk_index[tid] for tid in tids if tid in trk_index]
        if not idxs:
            continue
        subt = tracklets_tab[idxs]
        t0 = np.minimum(subt["mjd1"], subt["mjd2"])
        order = np.argsort(t0)
        subt = subt[order]
        first = subt[0]
        last = subt[-1]

        ra1[lid] = float(first["ra1"]); dec1[lid] = float(first["dec1"]); mjd1[lid] = float(first["mjd1"])
        ra2[lid] = float(first["ra2"]); dec2[lid] = float(first["dec2"]); mjd2[lid] = float(first["mjd2"])

        ra3[lid] = float(last["ra1"]);  dec3[lid] = float(last["dec1"]);  mjd3[lid] = float(last["mjd1"])
        ra4[lid] = float(last["ra2"]);  dec4[lid] = float(last["dec2"]);  mjd4[lid] = float(last["mjd2"])

    links_tab = Table(rows=link_rows, names=["linkage_id", "n_tracklets", "n_nights"])
    links_tab["ra1"] = ra1; links_tab["dec1"] = dec1; links_tab["ra2"] = ra2; links_tab["dec2"] = dec2
    links_tab["ra3"] = ra3; links_tab["dec3"] = dec3; links_tab["ra4"] = ra4; links_tab["dec4"] = dec4
    links_tab["mjd1"] = mjd1; links_tab["mjd2"] = mjd2; links_tab["mjd3"] = mjd3; links_tab["mjd4"] = mjd4

    members_tab = Table(
        [np.array(mem_linkage, dtype=np.int64), np.array(mem_tracklet, dtype="U64")],
        names=["linkage_id", "tracklet_id"],
    )

    out_links = outdir / "links_tracklets.fits"
    out_members = outdir / "linkage_members.fits"
    links_tab.write(out_links, overwrite=True)
    members_tab.write(out_members, overwrite=True)

    print(f"[write] {out_links}  n_links={len(links_tab)}")
    print(f"[write] {out_members}  n_rows={len(members_tab)}")
    return links_tab, members_tab


def validate_args(args):
    if args.cores < 1:
        raise ValueError("--cores must be >= 1")
    if args.max_v_kms <= 0:
        raise ValueError("--max-v-kms must be > 0")
    if args.min_init_earth_au <= 0:
        raise ValueError("--min-init-earth-au must be > 0")
    if args.ref_dt_days <= 0:
        raise ValueError("--ref-dt-days must be > 0")
    if args.tol <= 0:
        raise ValueError("--tol must be > 0")
    if args.min_len_obs < 2:
        raise ValueError("--min-len-obs must be >= 2")
    if args.min_nights < 1:
        raise ValueError("--min-nights must be >= 1")
    if args.k_neighbors_cap < 1:
        raise ValueError("--k-neighbors-cap must be >= 1")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--outdir", required=True)

    ap.add_argument("--cores", type=int, default=1)
    ap.add_argument("--max-v-kms", type=float, default=30.0)
    ap.add_argument("--min-init-earth-au", type=float, default=0.01)

    ap.add_argument("--ref-epoch-mode", choices=["mid", "start"], default="mid")
    ap.add_argument("--ref-dt-days", type=float, default=0.05)

    ap.add_argument("--tol", type=float, default=0.03)
    ap.add_argument("--min-len-obs", type=int, default=4)
    ap.add_argument("--min-nights", type=int, default=2)
    ap.add_argument("--k-neighbors-cap", type=int, default=300)

    ap.add_argument("--hypos", type=str, default=None)

    args = ap.parse_args()
    validate_args(args)

    infile = Path(args.infile).expanduser()
    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)
    logger = ProgressLogger(outdir)

    tab = Table.read(infile)
    need = ["ra1", "dec1", "ra2", "dec2", "mjd1", "mjd2", "tracklet_id"]
    miss = [c for c in need if c not in tab.colnames]
    if miss:
        raise RuntimeError(f"Missing columns in {infile}: {miss}")

    t1 = Time(tab["mjd1"], format="mjd", scale="utc")
    t2 = Time(tab["mjd2"], format="mjd", scale="utc")
    jd1 = t1.jd.astype(np.float64)
    jd2 = t2.jd.astype(np.float64)

    loc = earth_location_mpc327()
    r_so_1 = heliocentric_observer_xyz_au(t1, loc)
    r_so_2 = heliocentric_observer_xyz_au(t2, loc)

    epochs = choose_reference_epochs(jd1, jd2, args.ref_epoch_mode, args.ref_dt_days)
    logger.log(f"[info] infile={infile} n_tracklets={len(tab)} cores={args.cores} POLINEW={POLINEW}")
    logger.log(f"[info] reference epochs (JD): {epochs[0]:.6f}, {epochs[1]:.6f}")

    if args.hypos:
        hypos = []
        with open(args.hypos, "r") as f:
            for line in f:
                if not line.strip() or line.strip().startswith("#"):
                    continue
                r, rd, rdd = [float(x) for x in line.split()[:3]]
                hypos.append((r, rd, rdd))
    else:
        hypos = build_hypos_default()
    logger.log(f"[info] n_hypos = {len(hypos)}")

    hl = HeliolincRR_Tracklets(
        ra1=np.asarray(tab["ra1"], dtype=np.float64),
        dec1=np.asarray(tab["dec1"], dtype=np.float64),
        ra2=np.asarray(tab["ra2"], dtype=np.float64),
        dec2=np.asarray(tab["dec2"], dtype=np.float64),
        jd1=jd1,
        jd2=jd2,
        r_so_1=r_so_1,
        r_so_2=r_so_2,
        tracklet_ids=np.asarray(tab["tracklet_id"]).astype("U64"),
    )

    t_prop0 = time.perf_counter()
    results = hl.propagate_all(
        hypos=hypos,
        epochs_jd=epochs,
        cores=args.cores,
        max_v_kms=args.max_v_kms,
        min_init_earth_au=args.min_init_earth_au,
        logger=logger,
    )
    elapsed_propagation_s = time.perf_counter() - t_prop0
    logger.log(f"[summary] propagation_total_elapsed={_fmt_elapsed(elapsed_propagation_s)}")

    t_cluster0 = time.perf_counter()
    links_obs = hl.cluster_all(
        results=results,
        tol=args.tol,
        min_len_obs=args.min_len_obs,
        min_nights=args.min_nights,
        k_neighbors_cap=args.k_neighbors_cap,
        logger=logger,
    )
    elapsed_clustering_s = time.perf_counter() - t_cluster0
    logger.log(f"[summary] clustering_total_elapsed={_fmt_elapsed(elapsed_clustering_s)}")
    logger.log(f"[ok] n_links (obs-level) = {len(links_obs)}")

    # keep debug arrays (useful)
    np.save(outdir / "bbf.npy", hl.bbf, allow_pickle=True)
    np.save(outdir / "utimes.npy", hl.utimes)

    # FITS outputs + endpoints in links table
    links_tab, members_tab = write_links_fits_with_endpoints(outdir, links_obs, hl.bbf, hl.utimes, tab)
    summary_path = write_rr_summary_json(
        outdir=outdir,
        infile=infile,
        args=args,
        hypos=hypos,
        epochs=epochs,
        results=results,
        links_tab=links_tab,
        members_tab=members_tab,
        elapsed_propagation_s=elapsed_propagation_s,
        elapsed_clustering_s=elapsed_clustering_s,
    )
    logger.log(f"[write] {summary_path}")
    logger.log(f"[done] outputs written to {outdir}")


if __name__ == "__main__":
    main()
