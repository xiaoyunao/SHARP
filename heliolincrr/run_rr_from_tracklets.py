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

PROFILE_DEFAULTS = {
    "single-night": {
        "cores": 16,
        "max_v_kms": 30.0,
        "min_init_earth_au": 0.02,
        "ref_epoch_mode": "mid",
        "ref_dt_days": 0.05,
        "tol": 0.03,
        "min_len_obs": 3,
        "min_nights": 1,
        "k_neighbors_cap": 300,
    },
    "w15": {
        "cores": 20,
        "max_v_kms": 200.0,
        "min_init_earth_au": 0.02,
        "ref_epoch_mode": "mid",
        "ref_dt_days": 0.50,
        "tol": 0.02,
        "min_len_obs": 4,
        "min_nights": 2,
        "k_neighbors_cap": 200,
    },
}


def earth_location_mpc327() -> EarthLocation:
    return EarthLocation.from_geodetic(
        lon=MPC327_LON_DEG * u.deg,
        lat=MPC327_LAT_DEG * u.deg,
        height=MPC327_ALT_M * u.m,
    )


def build_detection_key_array(file_col, objid_col) -> np.ndarray:
    return np.asarray(
        [f"{Path(str(f)).name}:{int(o)}" for f, o in zip(file_col, objid_col)],
        dtype="U256",
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
    def __init__(
        self,
        ra1,
        dec1,
        ra2,
        dec2,
        jd1,
        jd2,
        r_so_1,
        r_so_2,
        tracklet_ids,
        det_key1=None,
        det_key2=None,
    ):
        self.tracklet_ids = np.asarray(tracklet_ids).astype("U64")
        n = len(ra1)
        self.start_det_key = None
        self.end_det_key = None

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
        self.comb_u_mean = self.comb_r_ob0 + self.comb_r_ob1
        self.comb_u_mean /= np.clip(LA.norm(self.comb_u_mean, axis=1, keepdims=True), 1e-12, None)

        if det_key1 is not None and det_key2 is not None:
            det_key1 = np.asarray(det_key1).astype("U256")
            det_key2 = np.asarray(det_key2).astype("U256")
            start_keys = np.empty(n, dtype="U256")
            end_keys = np.empty(n, dtype="U256")
            first_is_start = np.asarray(jd1, dtype=np.float64) <= np.asarray(jd2, dtype=np.float64)
            start_keys[first_is_start] = det_key1[first_is_start]
            end_keys[first_is_start] = det_key2[first_is_start]
            start_keys[~first_is_start] = det_key2[~first_is_start]
            end_keys[~first_is_start] = det_key1[~first_is_start]
            self.start_det_key = start_keys
            self.end_det_key = end_keys

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

    @staticmethod
    def _connected_components_from_adjacency(adj: list[list[int]]) -> list[list[int]]:
        comps = []
        seen = set()
        for start in range(len(adj)):
            if start in seen:
                continue
            stack = [start]
            seen.add(start)
            comp = []
            while stack:
                cur = stack.pop()
                comp.append(cur)
                for nxt in adj[cur]:
                    if nxt in seen:
                        continue
                    seen.add(nxt)
                    stack.append(nxt)
            comps.append(sorted(comp))
        return comps

    @staticmethod
    def _pairwise_distances(points: np.ndarray) -> np.ndarray:
        if len(points) == 0:
            return np.empty((0, 0), dtype=np.float64)
        diff = points[:, None, :] - points[None, :, :]
        return np.sqrt(np.sum(diff * diff, axis=2))

    @staticmethod
    def _pairwise_angles_deg(unit_vectors: np.ndarray) -> np.ndarray:
        if len(unit_vectors) == 0:
            return np.empty((0, 0), dtype=np.float64)
        dot = np.clip(unit_vectors @ unit_vectors.T, -1.0, 1.0)
        return np.rad2deg(np.arccos(dot))

    @staticmethod
    def _compact_core_indices(
        points: np.ndarray,
        sky_units: np.ndarray,
        diameter_limit: float,
        centroid_limit: float,
        motion_limit: float,
        sky_pair_soft_deg: float,
        sky_centroid_soft_deg: float,
    ) -> list[int]:
        active = list(range(len(points)))
        if not active:
            return []

        while len(active) >= 2:
            sub = points[active]
            sub_sky = sky_units[active]
            pos0 = sub[:, 0:3]
            pos1 = sub[:, 3:6]
            disp = pos1 - pos0

            centroid = np.mean(sub, axis=0)
            centroid_dist = np.sqrt(np.sum((sub - centroid) ** 2, axis=1))
            disp_centroid = np.mean(disp, axis=0)
            motion_dist = np.sqrt(np.sum((disp - disp_centroid) ** 2, axis=1))
            sky_centroid = np.mean(sub_sky, axis=0)
            sky_centroid /= max(float(LA.norm(sky_centroid)), 1e-12)
            sky_dot = np.clip(sub_sky @ sky_centroid, -1.0, 1.0)
            sky_centroid_dist = np.rad2deg(np.arccos(sky_dot))
            pairwise = HeliolincRR_Tracklets._pairwise_distances(sub)
            sky_pairwise = HeliolincRR_Tracklets._pairwise_angles_deg(sub_sky)
            max_pair = float(np.max(pairwise)) if pairwise.size else 0.0
            max_centroid = float(np.max(centroid_dist)) if centroid_dist.size else 0.0
            max_motion = float(np.max(motion_dist)) if motion_dist.size else 0.0
            if (
                max_pair <= diameter_limit
                and max_centroid <= centroid_limit
                and max_motion <= motion_limit
            ):
                return active

            mean_pair = np.mean(pairwise, axis=1) if pairwise.size else np.zeros(len(active), dtype=np.float64)
            score = np.maximum(
                centroid_dist / max(centroid_limit, 1e-12),
                mean_pair / max(diameter_limit, 1e-12),
            )
            score = np.maximum(score, motion_dist / max(motion_limit, 1e-12))

            # Treat sky pointing as a soft outlier signal: large-angle members are clipped first,
            # but a mostly good >=3-tracklet core should survive instead of dropping the whole link.
            if len(active) >= 3:
                score = np.maximum(score, sky_centroid_dist / max(sky_centroid_soft_deg, 1e-12))
                if sky_pairwise.size:
                    score = np.maximum(score, np.max(sky_pairwise, axis=1) / max(sky_pair_soft_deg, 1e-12))
            else:
                # With only two tracklets left, a large sky separation means the pair itself is not viable.
                max_sky_pair = float(np.max(sky_pairwise)) if sky_pairwise.size else 0.0
                max_sky_centroid = float(np.max(sky_centroid_dist)) if sky_centroid_dist.size else 0.0
                if max_sky_pair <= sky_pair_soft_deg and max_sky_centroid <= sky_centroid_soft_deg:
                    return active
                score = np.maximum(score, sky_centroid_dist / max(sky_centroid_soft_deg, 1e-12))

            worst_local = int(np.argmax(score))
            del active[worst_local]

        return active

    def _component_links_from_members(
        self,
        member_indices: list[int],
        drvs: np.ndarray,
        comb_obs: list[tuple[int, ...]],
        tol: float,
        min_len_obs: int,
        profile: str,
    ) -> list[list[int]]:
        if not member_indices:
            return []

        if profile == "single-night":
            split_edge = float(tol * 0.88)
            diameter_limit = float(tol * 1.60)
            centroid_limit = float(tol * 0.78)
            motion_limit = np.inf
            sky_pair_soft_deg = 10.0
            sky_centroid_soft_deg = 6.0
        else:
            split_edge = float(tol * 0.90)
            diameter_limit = float(tol * 1.60)
            centroid_limit = float(tol * 0.80)
            motion_limit = np.inf
            sky_pair_soft_deg = 180.0
            sky_centroid_soft_deg = 180.0

        points = drvs[np.asarray(member_indices, dtype=np.int64)]
        sky_units = self.comb_u_mean[np.asarray(member_indices, dtype=np.int64)]
        pairwise = self._pairwise_distances(points)
        adjacency = [[] for _ in range(len(member_indices))]
        for i in range(len(member_indices)):
            for j in range(i + 1, len(member_indices)):
                if pairwise[i, j] <= split_edge:
                    adjacency[i].append(j)
                    adjacency[j].append(i)

        subcomponents = self._connected_components_from_adjacency(adjacency)
        links = []
        for local_members in subcomponents:
            if not local_members:
                continue
            sub_points = points[np.asarray(local_members, dtype=np.int64)]
            sub_sky_units = sky_units[np.asarray(local_members, dtype=np.int64)]
            kept_local = self._compact_core_indices(
                sub_points,
                sub_sky_units,
                diameter_limit,
                centroid_limit,
                motion_limit,
                sky_pair_soft_deg,
                sky_centroid_soft_deg,
            )
            if len(kept_local) == 0:
                continue
            kept_members = [member_indices[local_members[i]] for i in kept_local]
            link_obs = sorted({obs for idx in kept_members for obs in comb_obs[idx]})
            if len(link_obs) >= min_len_obs:
                links.append(link_obs)

        return links

    def cluster_one(self, prvs, combs_used, hypo, tol, min_len_obs, min_nights, k_neighbors_cap, profile):
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
                if profile == "single-night" and self.start_det_key is not None and self.end_det_key is not None:
                    can_chain = (
                        self.end_det_key[i] == self.start_det_key[j]
                        or self.end_det_key[j] == self.start_det_key[i]
                    )
                    if not can_chain:
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

        components = []
        seen = set()
        for i in range(num_tracklets):
            root = find(i)
            if root in seen:
                continue
            seen.add(root)
            members = [j for j in range(num_tracklets) if find(j) == root]
            components.append(members)

        links = []
        for members in components:
            links.extend(
                self._component_links_from_members(
                    member_indices=members,
                    drvs=drvs,
                    comb_obs=comb_obs,
                    tol=tol,
                    min_len_obs=min_len_obs,
                    profile=profile,
                )
            )

        links = self._min_nights_filter(links, self.bbf, self.utimes, min_nights)
        return self._unique_links(links)

    def cluster_all(self, results, tol, min_len_obs, min_nights, k_neighbors_cap, profile, logger: ProgressLogger | None = None):
        all_links = []
        t0 = time.perf_counter()
        n_hypos = len(results)
        for i, (prvs, combs_used, hypo) in enumerate(results, start=1):
            links = self.cluster_one(prvs, combs_used, hypo, tol, min_len_obs, min_nights, k_neighbors_cap, profile)
            all_links.extend(links)
            if logger is not None:
                logger.log(
                    f"[cluster] {i}/{n_hypos} hypo={hypo} n_prvs={len(prvs)} n_links_added={len(links)} elapsed={_fmt_elapsed(time.perf_counter() - t0)}"
                )
        if logger is not None:
            logger.log(f"[cluster] done n_hypos={n_hypos} elapsed={_fmt_elapsed(time.perf_counter() - t0)}")
        return self._unique_links(all_links)


def build_hypos_default():
    # (r, rdot, rddot) are heliocentric range hypotheses in AU, AU/day, AU/day^2.
    rs = np.array([1.3, 1.7, 2.1, 2.5, 2.9, 3.3, 3.7, 4.1], dtype=np.float64)
    rdots = np.array([0.0], dtype=np.float64)
    rdd = np.array([0.0], dtype=np.float64)
    return [(float(r), float(rd), float(rdd0)) for r in rs for rd in rdots for rdd0 in rdd]


def apply_profile_defaults(args):
    defaults = PROFILE_DEFAULTS[args.profile]
    for key, value in defaults.items():
        if getattr(args, key) is None:
            setattr(args, key, value)
    return args


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
    ap.add_argument("--profile", choices=sorted(PROFILE_DEFAULTS.keys()), default="single-night")

    ap.add_argument("--cores", type=int, default=None)
    ap.add_argument("--max-v-kms", type=float, default=None)
    ap.add_argument("--min-init-earth-au", type=float, default=None)

    ap.add_argument("--ref-epoch-mode", choices=["mid", "start"], default=None)
    ap.add_argument("--ref-dt-days", type=float, default=None)

    ap.add_argument("--tol", type=float, default=None)
    ap.add_argument("--min-len-obs", type=int, default=None)
    ap.add_argument("--min-nights", type=int, default=None)
    ap.add_argument("--k-neighbors-cap", type=int, default=None)

    ap.add_argument("--hypos", type=str, default=None)

    args = ap.parse_args()
    args = apply_profile_defaults(args)
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
        det_key1=build_detection_key_array(tab["file1"], tab["objID1"]) if {"file1", "objID1"}.issubset(tab.colnames) else None,
        det_key2=build_detection_key_array(tab["file2"], tab["objID2"]) if {"file2", "objID2"}.issubset(tab.colnames) else None,
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
        profile=args.profile,
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
