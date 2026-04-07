#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Dump all detections referenced by links/tracklets into a single FITS table.

For every (linkage_id, tracklet_id, endpoint in {1,2}) that appears in RR links,
this script fetches the corresponding row from the Gaia-removed L2 source catalog
using (file1/objID1) and (file2/objID2) stored in tracklets FITS.

It writes a merged FITS containing:
  - leading metadata columns: linkage_id, fit_ok, is_good, n_nights, n_tracklets,
    tracklet_id, endpoint, l2_file, objID, mjd, ra, dec
  - plus ALL columns from that L2 catalog row (verbatim)

Duplicates across tracklets/links are kept (no dedup).

Typical usage:
  python dump_link_detections_from_L2.py \
    --rr-dir /pipeline/xiaoyunao/data/heliolincrr/20251206/rr_links \
    --root   /pipeline/xiaoyunao/data/heliolincrr \
    --orbit-confirm /pipeline/xiaoyunao/data/heliolincrr/20251206/rr_links/orbit_confirm \
    --out /pipeline/xiaoyunao/data/heliolincrr/20251206/rr_links/orbit_confirm/link_detections_fullL2.fits \
    --which all

Options:
  --which all|good|fitok   choose which links to include using orbit_links.fits
  --include-bad-tracklets  keep endpoints even if tracklet row missing (filled with NaNs) [default: drop]
  --max-links N            debug limit links
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
from astropy.table import Table, vstack, hstack


# ----------------- helpers: resolve paths -----------------

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


def resolve_tracklets_path(rr_dir: Path, root: Path, explicit: Optional[str]) -> Path:
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


def load_rr_members(rr_dir: Path) -> Table:
    return Table.read(rr_dir / "linkage_members.fits")


def load_orbit_links(orbit_confirm_dir: Path) -> Table:
    return Table.read(orbit_confirm_dir / "orbit_links.fits")


# ----------------- helpers: tracklet index & endpoint extraction -----------------

def build_tracklet_index(tracklets: Table) -> Dict[str, int]:
    tids = np.asarray(tracklets["tracklet_id"]).astype("U64")
    return {tid: i for i, tid in enumerate(tids)}


def require_cols(tab: Table, cols: List[str], what: str):
    missing = [c for c in cols if c not in tab.colnames]
    if missing:
        raise RuntimeError(f"{what} missing required columns: {missing}")


def get_tracklet_endpoints(tracklets: Table, idx: Dict[str, int], tid: str):
    """
    Return two endpoints for a tracklet as:
      [(endpoint, file, objID, mjd, ra, dec), ...]
    """
    j = idx.get(tid, None)
    if j is None:
        return None

    tr = tracklets[j]
    return [
        (1, str(tr["file1"]), int(tr["objID1"]), float(tr["mjd1"]), float(tr["ra1"]), float(tr["dec1"])),
        (2, str(tr["file2"]), int(tr["objID2"]), float(tr["mjd2"]), float(tr["ra2"]), float(tr["dec2"])),
    ]


# ----------------- helpers: L2 row fetch with caching -----------------

class L2Cache:
    """
    Cache L2 catalog tables by filename to avoid re-reading FITS many times.
    Keeps 1 table in memory per unique file encountered.
    """
    def __init__(self, root: Path, assume_paths_are_absolute: bool = False):
        self.root = root
        self.assume_abs = assume_paths_are_absolute
        self._cache: Dict[str, Table] = {}

    def resolve_path(self, file_str: str) -> Path:
        p = Path(file_str)
        if p.is_absolute() or self.assume_abs:
            return p
        # Most pipelines store file paths as full paths; but if not, try root-relative.
        # You can customize this if needed.
        return (self.root / p).expanduser()

    def get_table(self, file_str: str) -> Table:
        key = file_str
        if key not in self._cache:
            path = self.resolve_path(file_str)
            if not path.exists():
                raise FileNotFoundError(f"L2 catalog not found: {path}")
            self._cache[key] = Table.read(path)
        return self._cache[key]

    def fetch_row(self, file_str: str, obj_id: int) -> Table:
        tab = self.get_table(file_str)
        # objID usually refers to row index in that catalog (0-based). If yours is 1-based, change here.
        if obj_id < 0 or obj_id >= len(tab):
            raise IndexError(f"objID out of range: objID={obj_id}  nrow={len(tab)}  file={file_str}")
        return tab[obj_id:obj_id+1]  # keep as 1-row Table


# ----------------- main dump logic -----------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rr-dir", required=True, help="RR output dir containing linkage_members.fits")
    ap.add_argument("--root", required=True, help="pipeline root (used to resolve tracklets + possibly L2 relative paths)")
    ap.add_argument("--orbit-confirm", default=None, help="orbit_confirm dir containing orbit_links.fits (default: <rr_dir>/orbit_confirm)")
    ap.add_argument("--tracklets", default=None, help="explicit tracklets FITS (optional)")
    ap.add_argument("--out", required=True, help="output FITS file path")
    ap.add_argument("--which", choices=["all", "good", "fitok"], default="all",
                    help="select links using orbit_links.fits: all / is_good / fit_ok")
    ap.add_argument("--max-links", type=int, default=0, help="debug: only first N links after filtering")
    ap.add_argument("--assume-l2-paths-absolute", action="store_true",
                    help="assume file1/file2 stored in tracklets are absolute paths; otherwise try root-relative")
    ap.add_argument("--keep-missing", action="store_true",
                    help="keep endpoints even if L2 row cannot be loaded (fill L2 columns with masked)")

    args = ap.parse_args()

    rr_dir = Path(args.rr_dir).expanduser()
    root = Path(args.root).expanduser()
    orbit_confirm_dir = Path(args.orbit_confirm).expanduser() if args.orbit_confirm else (rr_dir / "orbit_confirm")
    out_path = Path(args.out).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Load orbit link summary (for fit_ok/is_good flags)
    orbit_links = load_orbit_links(orbit_confirm_dir)
    require_cols(orbit_links, ["linkage_id", "fit_ok", "is_good", "n_nights", "n_tracklets"], "orbit_links")

    # Build selection set
    if args.which == "all":
        sel = np.ones(len(orbit_links), dtype=bool)
    elif args.which == "good":
        sel = np.asarray(orbit_links["is_good"], dtype=bool)
    else:  # fitok
        sel = np.asarray(orbit_links["fit_ok"], dtype=bool)

    orbit_links_sel = orbit_links[sel]
    if args.max_links and args.max_links > 0:
        orbit_links_sel = orbit_links_sel[: int(args.max_links)]

    # Map linkage_id -> flags
    lid_arr = np.asarray(orbit_links_sel["linkage_id"], dtype=int)
    fit_ok_map = {int(l): bool(v) for l, v in zip(lid_arr, np.asarray(orbit_links_sel["fit_ok"], dtype=bool))}
    is_good_map = {int(l): bool(v) for l, v in zip(lid_arr, np.asarray(orbit_links_sel["is_good"], dtype=bool))}
    n_nights_map = {int(l): int(v) for l, v in zip(lid_arr, np.asarray(orbit_links_sel["n_nights"], dtype=int))}
    n_trk_map = {int(l): int(v) for l, v in zip(lid_arr, np.asarray(orbit_links_sel["n_tracklets"], dtype=int))}

    # Load RR members and filter by selected linkage_ids
    members = load_rr_members(rr_dir)
    require_cols(members, ["linkage_id", "tracklet_id"], "linkage_members")

    members_lid = np.asarray(members["linkage_id"], dtype=int)
    keep_mem = np.isin(members_lid, lid_arr)
    members = members[keep_mem]

    # Load tracklets file
    tracklets_path = resolve_tracklets_path(rr_dir, root, args.tracklets)
    tracklets = Table.read(tracklets_path)
    require_cols(tracklets, ["tracklet_id", "file1", "file2", "objID1", "objID2", "mjd1", "mjd2", "ra1", "ra2", "dec1", "dec2"], "tracklets")
    tindex = build_tracklet_index(tracklets)

    # Cache for L2
    l2 = L2Cache(root=root, assume_paths_are_absolute=args.assume_l2_paths_absolute)

    # Collect rows grouped by L2 file (so vstack doesn’t fight schema too much)
    # But different files might have different columns; we will vstack with join_type='outer'.
    out_tables: List[Table] = []

    n_total = len(members)
    for k, (lid, tid) in enumerate(zip(np.asarray(members["linkage_id"], dtype=int),
                                       np.asarray(members["tracklet_id"]).astype("U64")), start=1):
        endpoints = get_tracklet_endpoints(tracklets, tindex, str(tid))
        if endpoints is None:
            # missing tracklet
            if not args.keep_missing:
                continue
            # emit meta-only row(s) with no L2 cols
            for ep in (1, 2):
                meta = Table(
                    rows=[(
                        int(lid),
                        bool(fit_ok_map.get(int(lid), False)),
                        bool(is_good_map.get(int(lid), False)),
                        int(n_nights_map.get(int(lid), -1)),
                        int(n_trk_map.get(int(lid), -1)),
                        str(tid),
                        int(ep),
                        "",
                        -1,
                        np.nan, np.nan, np.nan,
                    )],
                    names=["linkage_id", "fit_ok", "is_good", "n_nights", "n_tracklets",
                           "tracklet_id", "endpoint", "l2_file", "objID", "mjd", "ra_deg", "dec_deg"],
                )
                out_tables.append(meta)
            continue

        for ep, f, oid, mjd, ra, dec in endpoints:
            try:
                row = l2.fetch_row(f, oid)  # 1-row Table with ALL L2 columns
            except Exception:
                if not args.keep_missing:
                    continue
                row = Table()  # empty, will be outer-joined later

            # prepend meta columns
            meta = Table(
                rows=[(
                    int(lid),
                    bool(fit_ok_map.get(int(lid), False)),
                    bool(is_good_map.get(int(lid), False)),
                    int(n_nights_map.get(int(lid), -1)),
                    int(n_trk_map.get(int(lid), -1)),
                    str(tid),
                    int(ep),
                    str(f),
                    int(oid),
                    float(mjd),
                    float(ra),
                    float(dec),
                )],
                names=["linkage_id", "fit_ok", "is_good", "n_nights", "n_tracklets",
                       "tracklet_id", "endpoint", "l2_file", "objID", "mjd", "ra_deg", "dec_deg"],
            )

            if len(row.colnames) == 0:
                out_tables.append(meta)
            else:
                out_tables.append(hstack_safe(meta, row))

        if (k % 5000 == 0):
            print(f"[info] processed members {k}/{n_total}")

    if not out_tables:
        raise RuntimeError("No rows collected. Check inputs/filters.")

    # vstack with outer join (catalog schemas might differ)
    out = vstack(out_tables, join_type="outer", metadata_conflicts="silent")

    # provenance
    out.meta["rr_dir"] = str(rr_dir)
    out.meta["orbit_confirm_dir"] = str(orbit_confirm_dir)
    out.meta["tracklets"] = str(tracklets_path)
    out.meta["which"] = str(args.which)

    out.write(out_path, overwrite=True)
    print(f"[write] {out_path}  n_rows={len(out)}  n_cols={len(out.colnames)}")


def hstack_safe(a: Table, b: Table) -> Table:
    """
    Horizontal stack but avoid name collisions.
    If L2 catalog already has columns with same names as meta (unlikely), we suffix them.
    """
    b2 = b.copy()
    for c in a.colnames:
        if c in b2.colnames:
            b2.rename_column(c, f"{c}_l2")
    return hstack([a, b2], join_type="outer", metadata_conflicts="silent")


if __name__ == "__main__":
    main()
