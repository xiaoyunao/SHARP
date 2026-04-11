#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
import math
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.table import Table


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


def build_detection_keys(files, objids) -> np.ndarray:
    return np.asarray(
        [f"{Path(str(f)).name}:{int(o)}" for f, o in zip(files, objids)],
        dtype="U256",
    )


def wrap_angle_diff_deg(a: float, b: float) -> float:
    diff = abs(float(a) - float(b)) % 360.0
    return min(diff, 360.0 - diff)


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


def tracklet_nights(mjds: np.ndarray) -> int:
    if mjds.size == 0:
        return 0
    mjds = np.sort(np.asarray(mjds, dtype=np.float64))
    return int(np.sum(mjds[1:] - mjds[:-1] > 0.5) + 1)


def load_tracklets(infile: Path) -> Table:
    tab = Table.read(infile)
    need = [
        "tracklet_id",
        "ra1",
        "dec1",
        "ra2",
        "dec2",
        "mjd1",
        "mjd2",
        "file1",
        "file2",
        "objID1",
        "objID2",
    ]
    miss = [c for c in need if c not in tab.colnames]
    if miss:
        raise RuntimeError(f"Missing columns in {infile}: {miss}")
    return tab


def prepare_tracklet_rows(tab: Table):
    det1 = build_detection_keys(tab["file1"], tab["objID1"])
    det2 = build_detection_keys(tab["file2"], tab["objID2"])
    rows = []
    for i in range(len(tab)):
        t1 = float(tab["mjd1"][i])
        t2 = float(tab["mjd2"][i])
        ra1 = float(tab["ra1"][i])
        ra2 = float(tab["ra2"][i])
        dec1 = float(tab["dec1"][i])
        dec2 = float(tab["dec2"][i])
        if t1 <= t2:
            start_mjd, end_mjd = t1, t2
            start_ra, end_ra = ra1, ra2
            start_dec, end_dec = dec1, dec2
            start_key, end_key = str(det1[i]), str(det2[i])
        else:
            start_mjd, end_mjd = t2, t1
            start_ra, end_ra = ra2, ra1
            start_dec, end_dec = dec2, dec1
            start_key, end_key = str(det2[i]), str(det1[i])

        dt_hours = max((end_mjd - start_mjd) * 24.0, 1e-12)
        dra_arcsec = (end_ra - start_ra) * math.cos(math.radians(0.5 * (start_dec + end_dec))) * 3600.0
        ddec_arcsec = (end_dec - start_dec) * 3600.0
        speed = float(math.hypot(dra_arcsec, ddec_arcsec) / dt_hours)
        direction = float(math.degrees(math.atan2(ddec_arcsec, dra_arcsec)) % 360.0)
        rows.append(
            {
                "tracklet_id": str(tab["tracklet_id"][i]),
                "start_key": start_key,
                "end_key": end_key,
                "start_mjd": start_mjd,
                "end_mjd": end_mjd,
                "speed_arcsec_per_hour": speed,
                "direction_deg": direction,
            }
        )
    return rows


def build_edges(
    rows,
    speed_thresh_arcsec_per_hour: float,
    direction_thresh_deg: float,
    require_shared_endpoint: bool,
):
    by_start = defaultdict(list)
    for idx, row in enumerate(rows):
        by_start[row["start_key"]].append(idx)

    outgoing = defaultdict(list)
    incoming = defaultdict(list)
    edge_rows = []
    for i, row in enumerate(rows):
        if require_shared_endpoint:
            candidates = by_start.get(row["end_key"], [])
        else:
            candidates = range(len(rows))
        for j in candidates:
            if i == j:
                continue
            other = rows[j]
            shared_key = row["end_key"] if row["end_key"] == other["start_key"] else ""
            if require_shared_endpoint and not shared_key:
                continue
            if other["end_mjd"] <= row["end_mjd"]:
                continue
            speed_diff = abs(row["speed_arcsec_per_hour"] - other["speed_arcsec_per_hour"])
            dir_diff = wrap_angle_diff_deg(row["direction_deg"], other["direction_deg"])
            if speed_diff > float(speed_thresh_arcsec_per_hour):
                continue
            if dir_diff > float(direction_thresh_deg):
                continue
            outgoing[i].append(j)
            incoming[j].append(i)
            edge_rows.append(
                (
                    rows[i]["tracklet_id"],
                    rows[j]["tracklet_id"],
                    shared_key,
                    float(speed_diff),
                    float(dir_diff),
                )
            )
    return outgoing, incoming, edge_rows


def enumerate_maximal_paths(rows, outgoing, incoming):
    starts = [i for i in range(len(rows)) if len(incoming.get(i, [])) == 0 and len(outgoing.get(i, [])) > 0]
    if not starts:
        starts = [i for i in range(len(rows)) if len(outgoing.get(i, [])) > 0]

    paths = set()

    def dfs(path):
        cur = path[-1]
        nxts = outgoing.get(cur, [])
        if not nxts:
            if len(path) >= 2:
                paths.add(tuple(path))
            return
        extended = False
        for nxt in nxts:
            if nxt in path:
                continue
            extended = True
            dfs(path + [nxt])
        if not extended and len(path) >= 2:
            paths.add(tuple(path))

    for start in starts:
        dfs([start])

    if not paths:
        for src, dsts in outgoing.items():
            for dst in dsts:
                paths.add((src, dst))

    return [list(p) for p in sorted(paths)]


def build_link_obs(rows, path):
    obs = []
    for i, idx in enumerate(path):
        row = rows[idx]
        if i == 0:
            obs.append((row["start_mjd"], row["start_key"]))
        obs.append((row["end_mjd"], row["end_key"]))
    return obs


def write_outputs(outdir: Path, tab: Table, rows, paths, edge_rows, speed_thresh, dir_thresh):
    outdir.mkdir(parents=True, exist_ok=True)
    tracklet_index = {str(tab["tracklet_id"][i]): i for i in range(len(tab))}

    link_rows = []
    member_linkage = []
    member_tracklet = []
    n_links = len(paths)

    def nancol():
        return np.full(n_links, np.nan, dtype=np.float64)

    ra1 = nancol(); dec1 = nancol(); ra2 = nancol(); dec2 = nancol()
    ra3 = nancol(); dec3 = nancol(); ra4 = nancol(); dec4 = nancol()
    mjd1 = nancol(); mjd2 = nancol(); mjd3 = nancol(); mjd4 = nancol()

    for lid, path in enumerate(paths):
        tids = [rows[idx]["tracklet_id"] for idx in path]
        idxs = [tracklet_index[tid] for tid in tids]
        sub = tab[idxs]
        member_linkage.extend([lid] * len(tids))
        member_tracklet.extend(tids)

        obs = build_link_obs(rows, path)
        link_rows.append((lid, len(tids), tracklet_nights(np.asarray([x[0] for x in obs], dtype=np.float64))))

        first = sub[0]
        last = sub[-1]
        ra1[lid] = float(first["ra1"]); dec1[lid] = float(first["dec1"]); mjd1[lid] = float(first["mjd1"])
        ra2[lid] = float(first["ra2"]); dec2[lid] = float(first["dec2"]); mjd2[lid] = float(first["mjd2"])
        ra3[lid] = float(last["ra1"]); dec3[lid] = float(last["dec1"]); mjd3[lid] = float(last["mjd1"])
        ra4[lid] = float(last["ra2"]); dec4[lid] = float(last["dec2"]); mjd4[lid] = float(last["mjd2"])

    links_tab = Table(rows=link_rows, names=["linkage_id", "n_tracklets", "n_nights"])
    links_tab["ra1"] = ra1; links_tab["dec1"] = dec1; links_tab["ra2"] = ra2; links_tab["dec2"] = dec2
    links_tab["ra3"] = ra3; links_tab["dec3"] = dec3; links_tab["ra4"] = ra4; links_tab["dec4"] = dec4
    links_tab["mjd1"] = mjd1; links_tab["mjd2"] = mjd2; links_tab["mjd3"] = mjd3; links_tab["mjd4"] = mjd4

    members_tab = Table(
        [np.asarray(member_linkage, dtype=np.int64), np.asarray(member_tracklet, dtype="U64")],
        names=["linkage_id", "tracklet_id"],
    )
    edges_tab = Table(
        rows=edge_rows,
        names=["src_tracklet_id", "dst_tracklet_id", "shared_det_key", "speed_diff_arcsec_per_hour", "direction_diff_deg"],
    )

    out_links = outdir / "links_tracklets.fits"
    out_members = outdir / "linkage_members.fits"
    out_edges = outdir / "link_edges.fits"
    links_tab.write(out_links, overwrite=True)
    members_tab.write(out_members, overwrite=True)
    edges_tab.write(out_edges, overwrite=True)

    summary = {
        "params": {
            "speed_thresh_arcsec_per_hour": float(speed_thresh),
            "direction_thresh_deg": float(dir_thresh),
        },
        "outputs": {
            "n_links": int(len(links_tab)),
            "n_member_rows": int(len(members_tab)),
            "n_edges": int(len(edges_tab)),
            "n_tracklets_per_link": summarize_int_column(links_tab["n_tracklets"] if len(links_tab) else []),
            "n_nights_per_link": summarize_int_column(links_tab["n_nights"] if len(links_tab) else []),
        },
    }
    out_summary = outdir / "rr_summary.json"
    out_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return out_links, out_members, out_edges, out_summary, links_tab, members_tab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--speed-thresh-arcsec-per-hour", type=float, default=5.0)
    ap.add_argument("--direction-thresh-deg", type=float, default=10.0)
    ap.add_argument("--require-shared-endpoint", action="store_true", default=True)
    ap.add_argument("--no-require-shared-endpoint", dest="require_shared_endpoint", action="store_false")
    args = ap.parse_args()

    infile = Path(args.infile).expanduser()
    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)
    logger = ProgressLogger(outdir)

    t0 = time.perf_counter()
    tab = load_tracklets(infile)
    logger.log(f"[info] infile={infile} n_tracklets={len(tab)}")
    rows = prepare_tracklet_rows(tab)
    logger.log(f"[build] prepared tracklets elapsed={_fmt_elapsed(time.perf_counter() - t0)}")

    outgoing, incoming, edge_rows = build_edges(
        rows=rows,
        speed_thresh_arcsec_per_hour=float(args.speed_thresh_arcsec_per_hour),
        direction_thresh_deg=float(args.direction_thresh_deg),
        require_shared_endpoint=bool(args.require_shared_endpoint),
    )
    logger.log(
        f"[edges] n_src_with_outgoing={sum(1 for k in outgoing if outgoing[k])} n_edges={len(edge_rows)} elapsed={_fmt_elapsed(time.perf_counter() - t0)}"
    )

    paths = enumerate_maximal_paths(rows, outgoing, incoming)
    logger.log(f"[paths] n_links={len(paths)} elapsed={_fmt_elapsed(time.perf_counter() - t0)}")

    out_links, out_members, out_edges, out_summary, links_tab, members_tab = write_outputs(
        outdir=outdir,
        tab=tab,
        rows=rows,
        paths=paths,
        edge_rows=edge_rows,
        speed_thresh=float(args.speed_thresh_arcsec_per_hour),
        dir_thresh=float(args.direction_thresh_deg),
    )

    logger.log(f"[write] {out_links} n_links={len(links_tab)}")
    logger.log(f"[write] {out_members} n_rows={len(members_tab)}")
    logger.log(f"[write] {out_edges} n_edges={len(edge_rows)}")
    logger.log(f"[write] {out_summary}")
    logger.log(f"[done] elapsed={_fmt_elapsed(time.perf_counter() - t0)}")


if __name__ == "__main__":
    main()
