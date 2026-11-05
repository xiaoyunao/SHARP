#!/usr/bin/env python3
"""Plot actual SHARP plan metadata exported by export_plot_data.py; no network.

Creates pointings_aitoff.png, monthly_pointings_YYYY-MM.png, and (when --night
and --iers are provided) night_YYYY-MM-DD.png. Input coordinates/timestamps are
never generated or replaced. Field borders are spherical FOV outlines. These
figures describe source plans, not successful observations or scientific yield.

Dependencies: numpy, matplotlib, astropy. Run as a module from the SHARP repository.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sys
import textwrap
from zoneinfo import ZoneInfo
from concurrent.futures import ProcessPoolExecutor
import multiprocessing


BJ = ZoneInfo("Asia/Shanghai")
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from targeted_observing.plot_style import MODES as COLORS, apply_style, save_pair

LINESTYLES = dict(EV="solid", AM="dashed", TJ="solid")
LABELS = {"EV": "Evening twilight", "AM": "Morning twilight", "TJ": "Trojan"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def utc(value):
    t = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if t.tzinfo is None or t.utcoffset().total_seconds() != 0:
        raise ValueError(f"Expected explicit UTC timestamp, got {value!r}")
    return t


def category(row):
    if row["mode"] == "trojan" or row["block_id"] == "TJ":
        return "TJ"
    for key in ("EV", "AM"):
        if row["block_id"].startswith(key):
            return key
    raise ValueError(f"Unknown source block/mode: {row['block_id']}/{row['mode']}")


def unique_fields(rows):
    result = {}
    for row in rows:
        key = (row["field_id"], row["field_version"])
        if key in result:
            old = result[key]
            if any(old[k] != row[k] for k in ("ra_deg", "dec_deg", "pa_deg", "fov_deg")):
                raise ValueError(f"Same field/version has inconsistent geometry: {key}")
        else:
            result[key] = row
    return list(result.values())


def counts(rows):
    ids = {str(k) for row in rows for name in ("core_targets", "secondary_targets") for k in row[name]}
    return {"night_block_field_rows": len(rows), "unique_field_versions": len(unique_fields(rows)),
            "unique_pointing_geometries": len({(r['ra_deg'], r['dec_deg'], r['pa_deg'], r['fov_deg']) for r in rows}),
            "unique_planned_target_ids": len(ids), "exposure_rows": sum(r["visits"] for r in rows)}


def load_bundle(path):
    required = ["manifest.json", "fields.csv", "exposures.csv.gz", "night_index.json"]
    missing = [k for k in required if not (path / k).is_file()]
    if missing:
        raise ValueError("Actual exported bundle is missing " + ", ".join(missing) +
                         "; run export_plot_data.py on real selected_nominal plans first. No substitute coordinates will be generated.")
    manifest = read_json(path / "manifest.json")
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported plot-bundle schema_version")
    for name in required[1:]:
        info = manifest.get("outputs", {}).get(name)
        if not info or sha(path / name) != info.get("sha256"):
            raise ValueError(f"Bundle hash missing or mismatched: {name}")
    index = read_json(path / "night_index.json")
    nights = {r["date"]: r for r in index["nights"]}
    if len(nights) != manifest["counts"]["nights"]:
        raise ValueError("Night index count disagrees with manifest")
    with (path / "fields.csv").open(encoding="utf-8", newline="") as handle:
        fields = list(csv.DictReader(handle))
    for row in fields:
        night = nights[row["date"]]
        conf = manifest["configurations"][night["config_id"]]
        if "fov_deg" not in conf:
            raise ValueError(f"Missing explicit field width for {row['date']}; not assuming 3 degrees")
        for key in ("ra_deg", "dec_deg", "pa_deg"):
            row[key] = float(row[key])
        row["fov_deg"] = float(conf["fov_deg"])
        if not 0 < row["fov_deg"] < 180:
            raise ValueError("Invalid source FOV width")
        row["visits"] = int(row["visits"])
        for key in ("core_targets", "secondary_targets"):
            row[key] = json.loads(row[key + "_json"])
        row["category"] = category(row)
    with gzip.open(path / "exposures.csv.gz", "rt", encoding="utf-8", newline="") as handle:
        exposures = list(csv.DictReader(handle))
    if len(fields) != manifest["counts"]["field_visits"] or len(exposures) != manifest["counts"]["exposures"]:
        raise ValueError("Exported table counts disagree with manifest")
    if sum(r["visits"] for r in fields) != len(exposures):
        raise ValueError("Field visit totals disagree with exposure rows")
    return manifest, nights, fields, exposures


def footprint(row):
    """Densely sampled square boundary in a spherical sky-offset frame."""
    h = row["fov_deg"] / 2
    line = np.linspace(-h, h, 17, endpoint=False)
    x = np.r_[line, np.full(17, h), -line, np.full(17, -h), -h]
    y = np.r_[np.full(17, -h), line, np.full(17, h), -line, -h]
    origin = SkyCoord(row["ra_deg"] * u.deg, row["dec_deg"] * u.deg, frame="icrs")
    # SHARP's current products use PA=0. The frame rotation is retained explicitly.
    frame = origin.skyoffset_frame(rotation=row["pa_deg"] * u.deg)
    border = SkyCoord(lon=x * u.deg, lat=y * u.deg, frame=frame).icrs
    return border.ra.deg, border.dec.deg


def map_segments(ra, dec):
    """Astronomical RA orientation, split exactly at the Aitoff map seam."""
    lon = (np.deg2rad(180 - ra) + np.pi) % (2 * np.pi) - np.pi
    lat = np.deg2rad(dec)
    segments = []; current = [(lon[0], lat[0])]
    for a, b, ya, yb in zip(lon[:-1], lon[1:], lat[:-1], lat[1:]):
        if abs(b - a) <= np.pi:
            current.append((b, yb)); continue
        adjusted = b - np.sign(b - a) * 2 * np.pi
        edge = np.pi if adjusted > a else -np.pi
        yedge = ya + (yb - ya) * (edge - a) / (adjusted - a)
        current.append((edge, yedge)); segments.append(current)
        current = [(-edge, yedge), (b, yb)]
    segments.append(current)
    return segments


def status_note(nights):
    reconstructed = any(n.get("reconstruction", {}).get("is_original_annual_plan") is False or
                        "TEST" in n.get("status", "") for n in nights.values())
    if reconstructed:
        return "RECONSTRUCTED TEST INPUT - not the original annual schedule"
    if all(n.get("observations_completed") is False for n in nights.values()):
        return "EXPORTED PLANNED POINTINGS - not executed observations"
    return "SOURCE PLAN METADATA - observation status retained in night_index.json"


def save_figure(fig, path, note, extra=""):
    footer = note + (" | " + extra if extra else "")
    save_pair(fig, path, footer)
    plt.close(fig)


def overview(fields, nights, output, note):
    fig = plt.figure(figsize=(13, 7.5)); ax = fig.add_subplot(111, projection="aitoff")
    legend = []
    for key in ("EV", "AM", "TJ"):
        part = [r for r in fields if r["category"] == key]
        segments = [s for r in unique_fields(part) for s in map_segments(*footprint(r))]
        if segments:
            ax.add_collection(LineCollection(segments, colors=COLORS[key], linestyles=LINESTYLES[key], linewidths=.65, alpha=.45, transform=ax.transData))
        legend.append(Line2D([], [], color=COLORS[key], ls=LINESTYLES[key], label=LABELS[key]))
    ticks = np.arange(-150, 180, 30)
    ax.set_xticks(np.deg2rad(ticks)); ax.set_xticklabels([f"{int((180-x)/15)}h" for x in ticks])
    ax.grid(alpha=.3); ax.set_xlabel("Right ascension (increases to the left)", labelpad=22)
    dates = sorted(nights); total = counts(fields)
    ax.set_title(f"Nominal pointings · {dates[0]} to {dates[-1]}", pad=20)
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(.5, -.08), frameon=False, ncol=3)
    fig.subplots_adjust(bottom=.22, top=.87)
    save_figure(fig, output / "pointings_aitoff.png", note,
                f"All modes: {total['unique_field_versions']} distinct field/version keys, {total['unique_planned_target_ids']} target IDs, {total['exposure_rows']} exposure rows. Overlapping fields are not independent sky area; target IDs are not recovered lightcurves.")


def monthly_maps(fields, nights, output, note):
    for month in sorted({day[:7] for day in nights}):
        fig, axes = plt.subplots(1, 3, figsize=(16, 5.4), layout="constrained")
        for ax, key in zip(axes, ("EV", "AM", "TJ")):
            part = [r for r in fields if r["date"].startswith(month) and r["category"] == key]
            unique = unique_fields(part)
            ax.set_title(LABELS[key])
            if not unique:
                ax.text(.5, .5, "No scheduled fields", transform=ax.transAxes, ha="center", color="#777777")
                ax.set_xticks([]); ax.set_yticks([]); continue
            borders = [footprint(r) for r in unique]
            ras = np.sort(np.concatenate([r for r, _ in borders]) % 360)
            gaps = np.diff(np.r_[ras, ras[0] + 360]); start = ras[(int(np.argmax(gaps)) + 1) % len(ras)]
            unwrap = lambda ra: np.mod(ra - start, 360) + start
            xs = []; ys = []
            for row, (ra, dec) in zip(unique, borders):
                x = unwrap(ra); xs.extend(x); ys.extend(dec)
                ax.plot(x, dec, color=COLORS[key], ls=LINESTYLES[key], alpha=.55, lw=.7)
            xp = max(1., np.ptp(xs) * .06); yp = max(1., np.ptp(ys) * .06)
            ax.set_xlim(max(xs) + xp, min(xs) - xp); ax.set_ylim(max(-90, min(ys) - yp), min(90, max(ys) + yp))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value % 360:g}"))
            ax.set_xlabel("RA (deg; increases left)"); ax.set_ylabel("Dec (deg)"); ax.grid(alpha=.2)
        fig.suptitle(month)
        save_figure(fig, output / f"monthly_pointings_{month}.png", note,
                    "Each panel has its own coordinate limits. Target counts use field membership IDs, not inferred target positions or scientific yield.")


def get_site(manifest, night, args):
    explicit = (args.site_lat, args.site_lon, args.site_height)
    if any(v is not None for v in explicit):
        if not all(v is not None for v in explicit):
            raise ValueError("Provide --site-lat, --site-lon and --site-height together")
        site = dict(zip(("lat_deg", "lon_deg", "height_m"), explicit)); basis = "explicit CLI override"
    else:
        config = manifest["configurations"].get(night["config_id"], {})
        site = config.get("site")
        if not site or not all(k in site for k in ("lat_deg", "lon_deg", "height_m")):
            raise ValueError("No explicit site in bundle; provide all three --site-* arguments. No default observatory is assumed.")
        site = {k: float(site[k]) for k in ("lat_deg", "lon_deg", "height_m")}; basis = "bundle manifest/config_id"
    if not -90 <= site["lat_deg"] <= 90 or not -180 <= site["lon_deg"] <= 180:
        raise ValueError("Invalid site latitude/longitude")
    return site, basis


def night_card(day, manifest, nights, fields, exposures, args, output, note, table):
    if day not in nights:
        raise ValueError(f"Night {day} is absent from the actual bundle")
    night = nights[day]; site, site_basis = get_site(manifest, night, args)
    config = manifest["configurations"][night["config_id"]]
    altitude_floor = float(config["min_alt_deg"])
    center_floor = float(config["center_alt_deg"])
    iers.earth_orientation_table.set(table)
    loc = EarthLocation.from_geodetic(site["lon_deg"] * u.deg, site["lat_deg"] * u.deg, site["height_m"] * u.m)
    lo, hi = float(night["night"]["evening"][0]), float(night["night"]["morning"][1])
    seconds = np.linspace(lo, hi, max(2, int(np.ceil((hi-lo)/60)) + 1)); t = Time(seconds, format="unix")
    if t.mjd.min() < table["MJD"][0].value or t.mjd.max() > table["MJD"][-1].value:
        raise ValueError(f"Supplied IERS does not cover the full plotted night {day}")
    frame = AltAz(obstime=t, location=loc, pressure=0 * u.hPa)
    sun = get_sun(t).transform_to(frame); moon = get_body("moon", t, loc).transform_to(frame)
    x = mdates.date2num([datetime.fromtimestamp(s, timezone.utc) for s in seconds])
    current = [r for r in fields if r["date"] == day]
    blocks = {b["block_id"]: b for b in night["blocks"]}
    rows = [r for r in exposures if r["date"] == day]
    lanes = sorted(current, key=lambda r: (utc(blocks[r["block_id"]]["not_before_utc"]), r["first_start_utc"]))
    h = max(8.5, 6 + len(lanes) * .20)
    fig, (ax, gantt) = plt.subplots(2, 1, figsize=(15, h), sharex=True, gridspec_kw={"height_ratios": [3, max(1.4, len(lanes)*.13)]})
    ax.plot(x, sun.alt.deg, color="#9d6815", lw=1.5, label="Sun (computed, airless)")
    ax.plot(x, moon.alt.deg, color="#555555", lw=1.2, ls="--", label="Moon (computed, airless)")
    for level, label, color in ((altitude_floor, "", "#b94c45"), (center_floor, "", "#d28c86"), (-12, "Sun -12 deg", "#9d6815"), (-18, "Sun -18 deg", "#c49b51")):
        ax.axhline(level, color=color, lw=.75, ls=":")
        if label:
            ax.text(1.002, level, label, transform=ax.get_yaxis_transform(), va="center", fontsize=8, color=color)
    ax.text(1.002, (altitude_floor + center_floor) / 2, f"{altitude_floor:g} / {center_floor:g} deg", transform=ax.get_yaxis_transform(), va="center", fontsize=8, color="#b94c45")
    markers = ("o", "s", "^", "D", "v", "P")
    ticks = []; labels = []
    for lane, f in enumerate(lanes):
        color = COLORS[f["category"]]; block = blocks[f["block_id"]]
        rr = [r for r in rows if r["block_id"] == f["block_id"] and r["field_id"] == f["field_id"]]
        point = SkyCoord(f["ra_deg"] * u.deg, f["dec_deg"] * u.deg, frame="icrs")
        alt = point.transform_to(frame).alt.deg
        begin, end = utc(block["not_before_utc"]).timestamp(), utc(block["hard_deadline_utc"]).timestamp()
        use = (seconds >= begin) & (seconds <= end)
        ax.plot(x[use], alt[use], color=color, alpha=.40, lw=.75)
        ts = Time([utc(r["mid_utc"]).timestamp() for r in rr], format="unix")
        mids = point.transform_to(AltAz(obstime=ts, location=loc, pressure=0*u.hPa)).alt.deg
        xm = mdates.date2num([utc(r["mid_utc"]) for r in rr])
        ax.scatter(xm, mids, s=16, color=color, marker=markers[lane % len(markers)], zorder=5)
        gantt.broken_barh([(mdates.date2num(utc(block["not_before_utc"])), (end-begin)/86400)], (lane-.36, .72), facecolors=color, alpha=.065)
        for r in rr:
            a, b = mdates.date2num(utc(r["start_utc"])), mdates.date2num(utc(r["end_utc"]))
            release = mdates.date2num(utc(r["release_utc"])); overhead = mdates.date2num(utc(r["overhead_start_utc"]))
            gantt.broken_barh([(a, b-a)], (lane-.30, .60), facecolors=color)
            if release > b:
                gantt.broken_barh([(b, release-b)], (lane-.30, .60), facecolors="#aaaaaa")
            if a > overhead:
                gantt.broken_barh([(overhead, a-overhead)], (lane-.30, .60), facecolors="#eeeeee", edgecolors="#999999", hatch="///", linewidth=.2)
        ticks.append(lane); labels.append(f"{lane+1:02d} {f['block_id']} / {f['field_id'][-6:]}")
    for item in night["free_intervals"]:
        a, b = mdates.date2num(utc(item["start_utc"])), mdates.date2num(utc(item["end_utc"]))
        ax.axvspan(a, b, color="#e5e5e5", alpha=.22, zorder=0)
    ax.set_ylim(-90, 90); ax.set_ylabel("Altitude (deg)"); ax.grid(alpha=.2)
    ax.legend(loc="lower left", fontsize=8, framealpha=.85)
    tj_count = sum(r["mode"] == "trojan" for r in rows)
    suffix = f" · {tj_count} Trojan exposures" if tj_count else " · no Trojan block"
    ax.set_title(day + suffix)
    gantt.set_yticks(ticks); gantt.set_yticklabels(labels, fontsize=8)
    if lanes:
        gantt.set_ylim(len(lanes)-.5, -.5)
    else:
        reasons = "; ".join(str(x) for x in night.get("skips", [])) or "Source plan contains no blocks."
        gantt.text(.5, .5, "Explicit zero-plan night: no exposure rows\n" + textwrap.fill(reasons, 105),
                   transform=gantt.transAxes, ha="center", va="center", fontsize=9)
    gantt.set_ylabel("Exported block / field"); gantt.grid(axis="x", alpha=.2)
    gantt.xaxis.set_major_locator(mdates.HourLocator(interval=1, tz=BJ))
    gantt.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=BJ)); gantt.set_xlim(x[0], x[-1])
    gantt.set_xlabel("Beijing time [UTC+08]")
    fig.legend(handles=[Patch(color="#3883b9", label="Exposure"), Patch(color="#aaaaaa", label="Readout"),
                          Patch(facecolor="#eeeeee", edgecolor="#999999", hatch="///", label="Slew / filter overhead")], loc="lower center", bbox_to_anchor=(.5, .015), ncol=3, fontsize=11, frameon=False)
    fig.subplots_adjust(top=.91, bottom=.15, left=.16, right=.91, hspace=.16)
    save_figure(fig, output / f"night_{day}.png", note,
                f"Site lat={site['lat_deg']:.4f}, lon={site['lon_deg']:.4f}, {site['height_m']:g} m ({site_basis}); local IERS {args.iers.name}. Altitude lines are config floors before sampling guards; no photometric success inferred.")
    return {"date": day, "exposure_rows": len(rows), "field_lanes": len(lanes), "site": site, "site_basis": site_basis,
            "iers_sha256": sha(args.iers), "airless": True, "environment_grid_step_s_max": float(np.diff(seconds).max())}


def _plot_night(day):
    return night_card(day, *_NIGHT_CONTEXT)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bundle", type=Path, help="output directory of export_plot_data.py")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--night", action="append", default=[], help="night YYYY-MM-DD; repeat to plot more than one")
    parser.add_argument("--all-nights", action="store_true", help="plot every night in the actual night index, including zero-plan nights")
    parser.add_argument("--workers", type=int, default=1, help="local numerical plotting processes; shared read-only input")
    parser.add_argument("--iers", type=Path, help="local covering IERS-A file, required for --night or --all-nights")
    parser.add_argument("--deps", type=Path, help="optional existing Python package directory")
    parser.add_argument("--no-monthly", action="store_true", help="skip month-by-month local views")
    parser.add_argument("--site-lat", type=float); parser.add_argument("--site-lon", type=float)
    parser.add_argument("--site-height", type=float, help="meters; all three --site-* values are required for an explicit override")
    args = parser.parse_args()
    if args.all_nights and args.night:
        parser.error("Use --all-nights or repeated --night, not both")
    if (args.night or args.all_nights) and (args.iers is None or not args.iers.is_file()):
        parser.error("Night plots require an existing local --iers file; this tool will not download one")
    if args.deps:
        sys.path.insert(0, str(args.deps.resolve()))
    global np, plt, u, SkyCoord, AltAz, EarthLocation, get_body, get_sun, Time, iers, mdates, LineCollection, Line2D, Patch, FuncFormatter
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter
    import astropy.units as u
    from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_body, get_sun, solar_system_ephemeris
    from astropy.time import Time
    from astropy.utils import iers
    iers.conf.auto_download = False; iers.conf.iers_degraded_accuracy = "error"
    solar_system_ephemeris.set("builtin")
    apply_style()
    manifest, nights, fields, exposures = load_bundle(args.bundle)
    if args.all_nights:
        args.night = sorted(nights)
    table = iers.IERS_A.open(str(args.iers)) if args.night else None
    for day in args.night:
        if day not in nights:
            raise ValueError(f"Requested night {day} is absent from the actual bundle")
        get_site(manifest, nights[day], args)
        bounds = Time([nights[day]["night"]["evening"][0], nights[day]["night"]["morning"][1]], format="unix").mjd
        if bounds.min() < table["MJD"][0].value or bounds.max() > table["MJD"][-1].value:
            raise ValueError(f"Supplied IERS does not cover requested night {day}; no automatic download or extrapolation")
    args.output.mkdir(parents=True, exist_ok=True)
    note = status_note(nights)
    overview(fields, nights, args.output, note)
    if not args.no_monthly:
        monthly_maps(fields, nights, args.output, note)
    global _NIGHT_CONTEXT
    _NIGHT_CONTEXT=(manifest, nights, fields, exposures, args, args.output, note, table)
    night_results=[]
    if args.workers>1 and args.night:
        # The initialized scientific imports are inherited read-only; each process
        # owns its figures and writes a distinct night. Never share pyplot threads.
        with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('fork')) as pool:
            for i,result in enumerate(pool.map(_plot_night,args.night),1):
                night_results.append(result)
                if i%30==0:print(f"Plotted {i}/{len(args.night)} nights",flush=True)
    else:
        night_results=[_plot_night(day) for day in args.night]
    result = {"source_manifest_sha256": sha(args.bundle / "manifest.json"), "script_sha256": sha(Path(__file__)),
              "source_commit": manifest["source"].get("commit"), "source_statuses": manifest["source"].get("statuses"),
              "source_note": note, "counts": counts(fields), "counts_by_mode": {k: counts([r for r in fields if r['category']==k]) for k in COLORS},
              "nights_plotted": night_results, "night_selection": "all_exported_nights" if args.all_nights else "explicit_dates",
              "scope": "Coordinates and times copied from actual export. Spherical outlines computed from explicit FOV/PA. Altitudes computed offline from explicit site/IERS. Unique planned IDs are not successful lightcurves; overlapping fields are not independent area."}
    (args.output / "plot_provenance.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "source_note": note, "counts": result["counts"], "nights_plotted": args.night}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, ImportError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
