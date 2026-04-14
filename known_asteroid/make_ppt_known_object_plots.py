#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import healpy as hp
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from aleph.Query import Query
from astropy.table import Table
from astroquery.jplsbdb import SBDB


plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["font.size"] = 13


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build PPT-ready known-object statistics plots.")
    parser.add_argument(
        "--all-matched-path",
        default="/pipeline/xiaoyunao/known_asteroid/runtime/history/all_matched_asteroids.fits",
    )
    parser.add_argument("--astorb", default="/pipeline/xiaoyunao/known_asteroid/astorb.dat")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--nside", type=int, default=64)
    parser.add_argument("--sbdb-cache", default=None)
    parser.add_argument("--max-sbdb", type=int, default=250)
    return parser.parse_args()


def normalize_name(value: object) -> str:
    text = str(value).strip()
    if text.startswith("b'") and text.endswith("'"):
        return text[2:-1].strip()
    return text


def object_key(number: float, name: str) -> str:
    if np.isfinite(number):
        return f"num:{int(number)}"
    return f"name:{name}"


def derive_orbit_class(a: float, e: float) -> str:
    if not np.isfinite(a) or not np.isfinite(e):
        return "Unclassified"
    q = a * (1.0 - e)
    q_big = a * (1.0 + e)
    if a < 1.0 and q_big < 0.983:
        return "Atira"
    if a < 1.0 and q_big >= 0.983:
        return "Aten"
    if a >= 1.0 and q <= 1.017:
        return "Apollo"
    if 1.017 < q < 1.3:
        return "Amor"
    if q < 1.666 and a > 1.3:
        return "Mars-crosser"
    if 1.78 <= a <= 2.0 and e < 0.18:
        return "Hungaria"
    if 3.7 <= a <= 4.2 and e < 0.35:
        return "Hilda"
    if 5.05 <= a <= 5.35 and e < 0.25:
        return "Jupiter Trojan"
    if 2.06 <= a < 2.5:
        return "Inner Main Belt"
    if 2.5 <= a < 2.82:
        return "Middle Main Belt"
    if 2.82 <= a < 3.28:
        return "Outer Main Belt"
    if a >= 30.0:
        return "Trans-Neptunian"
    if a >= 5.5:
        return "Centaur"
    return "Other"


def build_astorb_maps(astorb_path: str) -> tuple[dict[int, dict[str, float]], dict[str, dict[str, float]]]:
    q = Query(service="Lowell", filename=astorb_path)
    by_number: dict[int, dict[str, float]] = {}
    by_name: dict[str, dict[str, float]] = {}
    for row in q.asts:
        entry = {
            "orbit_elements_a": float(row["a"]),
            "orbit_elements_e": float(row["e"]),
            "orbit_elements_i": float(row["incl"]),
        }
        num_text = str(row["num"]).strip()
        if num_text.isdigit():
            by_number[int(num_text)] = entry
        name = normalize_name(row["name"])
        if name:
            by_name[name.lower()] = entry
    return by_number, by_name


def load_sbdb_cache(cache_path: Path) -> dict[str, dict[str, object]]:
    if not cache_path.exists():
        return {}
    return json.loads(cache_path.read_text(encoding="utf-8"))


def save_sbdb_cache(cache_path: Path, cache: dict[str, dict[str, object]]) -> None:
    cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def query_sbdb_entry(target: str) -> dict[str, object] | None:
    try:
        result = SBDB.query(target, full_precision=True)
    except Exception:
        return None
    orbit = result.get("orbit", {})
    elements = orbit.get("elements", [])
    values = {item.get("name"): item.get("value") for item in elements if isinstance(item, dict)}
    try:
        a = float(values["a"])
        e = float(values["e"])
        inc = float(values["i"])
    except Exception:
        return None
    orbit_class = result.get("object", {}).get("orbit_class", {}).get("name")
    return {
        "orbit_elements_a": a,
        "orbit_elements_e": e,
        "orbit_elements_i": inc,
        "object_orbit_class_name": orbit_class or derive_orbit_class(a, e),
    }


def enrich_objects(
    unique_rows: list[dict[str, object]],
    astorb_by_number: dict[int, dict[str, float]],
    astorb_by_name: dict[str, dict[str, float]],
    cache: dict[str, dict[str, object]],
    max_sbdb: int,
) -> tuple[list[dict[str, object]], Counter]:
    stats = Counter()
    enriched: list[dict[str, object]] = []
    sbdb_budget = max_sbdb
    for row in unique_rows:
        number = row["number"]
        name = row["name"]
        orbit = None
        if np.isfinite(number):
            orbit = astorb_by_number.get(int(number))
            if orbit is not None:
                stats["astorb_number"] += 1
        if orbit is None and name:
            orbit = astorb_by_name.get(name.lower())
            if orbit is not None:
                stats["astorb_name"] += 1
        cache_key = str(int(number)) if np.isfinite(number) else name
        if orbit is None and cache_key in cache:
            orbit = cache[cache_key]
            stats["sbdb_cache"] += 1
        if orbit is None and sbdb_budget > 0:
            queried = query_sbdb_entry(cache_key)
            if queried is not None:
                cache[cache_key] = queried
                orbit = queried
                sbdb_budget -= 1
                stats["sbdb_live"] += 1
        if orbit is None:
            stats["missing"] += 1
            continue
        a = float(orbit["orbit_elements_a"])
        e = float(orbit["orbit_elements_e"])
        inc = float(orbit["orbit_elements_i"])
        orbit_class = str(orbit.get("object_orbit_class_name") or derive_orbit_class(a, e))
        enriched.append(
            {
                **row,
                "orbit_elements_a": a,
                "orbit_elements_e": e,
                "orbit_elements_i": inc,
                "object_orbit_class_name": orbit_class,
            }
        )
    return enriched, stats


def make_unique_objects(table: Table) -> list[dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    names = np.asarray(table["name"])
    numbers = np.asarray(table["number"], dtype=float)
    mjd = np.asarray(table["MJD"], dtype=float)
    for idx in range(len(table)):
        name = normalize_name(names[idx])
        number = float(numbers[idx])
        key = object_key(number, name)
        prev = rows.get(key)
        current = {"key": key, "name": name, "number": number, "last_mjd": float(mjd[idx])}
        if prev is None or current["last_mjd"] > prev["last_mjd"]:
            rows[key] = current
    return list(rows.values())


def plot_orbits(enriched: list[dict[str, object]], outdir: Path) -> None:
    classes = [
        "Atira",
        "Aten",
        "Apollo",
        "Amor",
        "Mars-crosser",
        "Hungaria",
        "Inner Main Belt",
        "Middle Main Belt",
        "Outer Main Belt",
        "Hilda",
        "Jupiter Trojan",
        "Other",
    ]
    colors = {
        "Atira": "#9c6644",
        "Aten": "#b56576",
        "Apollo": "#e76f51",
        "Amor": "#f4a261",
        "Mars-crosser": "#f6bd60",
        "Hungaria": "#84a59d",
        "Inner Main Belt": "#277da1",
        "Middle Main Belt": "#4d908e",
        "Outer Main Belt": "#577590",
        "Hilda": "#8d99ae",
        "Jupiter Trojan": "#6d597a",
        "Other": "#6c757d",
    }
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    for orbit_class in classes:
        subset = [row for row in enriched if row["object_orbit_class_name"] == orbit_class]
        if not subset:
            continue
        a = np.asarray([row["orbit_elements_a"] for row in subset], dtype=float)
        e = np.asarray([row["orbit_elements_e"] for row in subset], dtype=float)
        inc = np.asarray([row["orbit_elements_i"] for row in subset], dtype=float)
        axes[0].scatter(a, e, s=8, alpha=0.55, color=colors[orbit_class], edgecolors="none", label=orbit_class)
        axes[1].scatter(a, inc, s=8, alpha=0.55, color=colors[orbit_class], edgecolors="none", label=orbit_class)
    axes[0].set_xlabel("Semi-major axis a (au)")
    axes[0].set_ylabel("Eccentricity e")
    axes[0].set_title("Known-object orbital distribution: a vs e")
    axes[1].set_xlabel("Semi-major axis a (au)")
    axes[1].set_ylabel("Inclination i (deg)")
    axes[1].set_title("Known-object orbital distribution: a vs i")
    axes[0].set_xlim(0, 6)
    axes[0].set_ylim(0, 1)
    axes[1].set_xlim(0, 6)
    axes[1].set_ylim(0, 40)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.subplots_adjust(top=0.90, bottom=0.12, left=0.07, right=0.82, wspace=0.10)
    fig.legend(handles, labels, loc="center left", ncol=1, frameon=False, bbox_to_anchor=(0.835, 0.5))
    fig.savefig(outdir / "known_object_orbits.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_healpix_counts(table: Table, nside: int, outdir: Path) -> None:
    ra = np.asarray(table["RA_Win"], dtype=float)
    dec = np.asarray(table["DEC_Win"], dtype=float)
    valid = np.isfinite(ra) & np.isfinite(dec)
    theta = np.deg2rad(90.0 - dec[valid])
    phi = np.deg2rad(ra[valid] % 360.0)
    pix = hp.ang2pix(nside, theta, phi, nest=False)
    counts = np.bincount(pix, minlength=hp.nside2npix(nside)).astype(float)

    fig = plt.figure(figsize=(13, 8))
    hp.mollview(
        counts,
        fig=fig.number,
        title=f"Known-object detections on sky (HEALPix nside={nside})",
        unit="Detections per pixel",
        cmap="viridis",
        cbar=True,
        notext=True,
    )
    hp.graticule(color="white", alpha=0.25)
    fig.savefig(outdir / f"known_object_nside{nside}_counts.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    table_out = Table()
    table_out["hpix"] = np.arange(len(counts), dtype=np.int64)
    table_out["count"] = counts.astype(np.int32)
    table_out.write(outdir / f"known_object_nside{nside}_counts.fits", overwrite=True)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cache_path = Path(args.sbdb_cache) if args.sbdb_cache else outdir / "sbdb_cache.json"

    table = Table.read(args.all_matched_path, memmap=True)
    unique_rows = make_unique_objects(table)
    astorb_by_number, astorb_by_name = build_astorb_maps(args.astorb)
    cache = load_sbdb_cache(cache_path)
    enriched, stats = enrich_objects(unique_rows, astorb_by_number, astorb_by_name, cache, args.max_sbdb)
    save_sbdb_cache(cache_path, cache)

    enriched_table = Table(rows=enriched)
    enriched_table.write(outdir / "known_object_orbits_enriched.fits", overwrite=True)
    plot_orbits(enriched, outdir)
    plot_healpix_counts(table, args.nside, outdir)

    summary = {
        "all_matched_rows": int(len(table)),
        "unique_objects": int(len(unique_rows)),
        "enriched_objects": int(len(enriched)),
        "class_counts": Counter([row["object_orbit_class_name"] for row in enriched]),
        "enrichment_stats": dict(stats),
    }
    (outdir / "known_object_plot_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=int) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=int))


if __name__ == "__main__":
    main()
