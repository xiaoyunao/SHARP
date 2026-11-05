#!/usr/bin/env python3
"""Offline representative-night geometry for a supplied SHARP Horizons index.

Uses existing snapshot geometry/catalog code and cached inputs only. This is not
an annual scheduler or a prediction of photometric/scientific yield. By default
it evaluates the 15th of each month, November 2026 through October 2027, and saves
five-minute curves for December 15 and April 15. See output provenance/rules.

Example:
  python3 compute_special_visibility.py --snapshot /path/to/snapshot \
    --deps /path/to/python-packages --index /path/to/special_ephemeris_index.json \
    --iers /path/to/snapshot/targeted_observing/cache/iers.all \
    --output /path/to/plot_data
"""
from __future__ import annotations

import argparse
import csv
from datetime import date, datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import sys
from zoneinfo import ZoneInfo


BEIJING = ZoneInfo("Asia/Shanghai")


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def iso_utc(seconds):
    return datetime.fromtimestamp(float(seconds), timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def iso_bj(seconds):
    return datetime.fromtimestamp(float(seconds), BEIJING).isoformat(timespec="seconds")


def intervals(times, series, condition):
    """Union of passing intervals of the piecewise-linear sampled quantities.

    Every threshold crossing is a subdivision, including moonrise and the 75%
    illumination switch. Returned boundaries are interpolated, not certified
    continuous astronomical roots.
    """
    import numpy as np
    result = []
    for i in range(len(times) - 1):
        left, right = series[:, i], series[:, i + 1]
        delta = right - left
        cross = (left * right < 0) & (delta != 0)
        cuts = np.unique(np.r_[0., -left[cross] / delta[cross], 1.])
        for a, b in zip(cuts[:-1], cuts[1:]):
            if condition(left + delta * ((a + b) / 2)):
                lo = float(times[i] + a * (times[i + 1] - times[i]))
                hi = float(times[i] + b * (times[i + 1] - times[i]))
                if result and abs(result[-1][1] - lo) < .001:
                    result[-1][1] = hi
                else:
                    result.append([lo, hi])
    return result


def interval_records(values):
    return [{"start_utc": iso_utc(a), "end_utc": iso_utc(b),
             "start_beijing": iso_bj(a), "end_beijing": iso_bj(b),
             "duration_h": round((b - a) / 3600, 6)} for a, b in values]


def month_dates(first, months, day):
    year, month = map(int, first.split("-"))
    output = []
    for i in range(months):
        n = year * 12 + month - 1 + i
        output.append(date(n // 12, n % 12 + 1, day).isoformat())
    return output


def write_json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--deps", type=Path, help="optional existing Python package directory")
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--iers", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--first-month", default="2026-11")
    parser.add_argument("--months", type=int, default=12)
    parser.add_argument("--day", type=int, default=15)
    parser.add_argument("--step-s", type=int, default=60, choices=[15, 30, 60, 300])
    parser.add_argument("--night-dates", nargs="+", default=["2026-12-15", "2027-04-15"])
    parser.add_argument("--source-commit", default="a5cb9454e56771ebe4211b7a544412056cce4c40")
    args = parser.parse_args()
    if args.deps:
        sys.path.insert(0, str(args.deps.resolve()))
    sys.path.insert(0, str(args.snapshot.resolve()))
    import numpy as np
    import astropy
    from astropy.time import Time
    from astropy.coordinates import solar_system_ephemeris
    from targeted_observing.geometry import Sky, setup_iers, coord
    from targeted_observing.catalog import interpolate

    config_path = args.snapshot / "targeted_observing/defaults.json"
    config = json.loads(config_path.read_text())
    iers_info = setup_iers(args.iers)
    solar_system_ephemeris.set("builtin")
    sky = Sky(config)
    index = json.loads(args.index.read_text())
    if index.get("failures"):
        raise ValueError(f"Input ephemeris index has failures: {index['failures']}")
    tracks = {}; source_files = {}
    for key, raw_path in index["paths"].items():
        path = Path(raw_path)
        if not path.exists():
            path = args.snapshot / "targeted_observing/cache/annual" / path.name
        track = json.loads(path.read_text())
        if str(track["pdes"]) != str(key) or track["site"] != config["site"]:
            raise ValueError(f"Ephemeris ID or site differs from input configuration: {path}")
        tracks[str(key)] = track
        source_files[str(key)] = {"file": path.name, "sha256": hash_file(path),
                                  "targetname": track["targetname"], "retrieved_utc": track.get("retrieved_utc"),
                                  "source": track.get("source"), "query": track.get("query")}
    if not tracks:
        raise ValueError("No ephemeris tracks in input index")
    ids = list(tracks)
    monthly_dates = month_dates(args.first_month, args.months, args.day)
    all_dates = sorted(set(monthly_dates) | set(args.night_dates))
    rules = {
        "representative_day": args.day, "sampling_step_s": args.step_s,
        "sun_altitude_limit_deg": -18., "target_altitude_limit_deg": 30.,
        "field_center_altitude_limit_deg": 32.8, "field_circumscribed_cap_deg": 2.122,
        "moon_separation_deg": config["moon_sep_deg"], "bright_moon_separation_deg": config["bright_moon_sep_deg"],
        "bright_moon_fraction": config["bright_moon_fraction"],
        "moon_rule": "Only when lunar center altitude >0: center-Moon separation exceeds 30 or 40 deg plus 2.122 deg cap; illumination >=0.75 selects 40 deg.",
        "field_center": "Fixed at each object's interpolated RA/Dec at 00:00 Beijing following the observing-night date. This is a representative target-centered field, not a selected or optimized pointing.",
        "magnitude": "Horizons predicted V at that same Beijing midnight; V is not a calibrated r magnitude.",
        "continuous_window_method": "Threshold crossings interpolated between samples; includes lunar horizon and illumination switches. No extra between-grid angular guard, no independent fine-grid execution validation.",
        "shared_time": "Astronomical darkness only; twilight hour reservations and other programmes are NOT subtracted. No weather, stellar masks, saturation, photometric threshold, or readout/slew costs applied.",
        "night_curve_grid": "Every 300 s from 16:00 Beijing to next-day 08:00, inclusive.",
    }
    provenance = {
        "scope": "17-object input subset (actual count below), 12 monthly representative nights; not a full-year extrema search, annual schedule, or scientific-yield prediction.",
        "source_commit_declared": args.source_commit,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": hash_file(Path(__file__)), "ephemeris_index_sha256": hash_file(args.index),
        "configuration_sha256": hash_file(config_path), "site": config["site"], "iers": iers_info,
        "geometry_py_sha256": hash_file(args.snapshot / "targeted_observing/geometry.py"),
        "catalog_py_sha256": hash_file(args.snapshot / "targeted_observing/catalog.py"),
        "ephemerides": source_files, "object_count": len(ids),
        "interpolation": "Existing catalog.interpolate: RA unwrapped, RA/Dec linearly interpolated in JD; V linearly interpolated from the same six-hour Horizons table; no extrapolation.",
        "sun_moon_model": "Astropy builtin solar-system ephemerides; airless AltAz, supplied frozen IERS, all network downloads disabled.",
        "packages": {"numpy": np.__version__, "astropy": astropy.__version__},
    }
    monthly = []; night_curves = []
    for day in all_dates:
        noon = datetime.fromisoformat(day).replace(hour=12, tzinfo=BEIJING)
        seconds = noon.timestamp() + np.arange(0, 86400 + 1, args.step_s, dtype=float)
        ref = noon.timestamp() + 12 * 3600
        frame, sun, moon, fraction = sky.environment(seconds)
        sun_alt = sun.alt.deg; moon_alt = moon.alt.deg
        ras = []; decs = []; mags = []; centers_ra = []; centers_dec = []; ref_mag = []
        jd = Time(seconds, format="unix").jd
        ref_jd = float(Time(ref, format="unix").jd)
        for key in ids:
            t = tracks[key]
            ra, dec = interpolate(t, seconds)
            r0, d0 = interpolate(t, [ref])
            v = np.interp(jd, t["datetime_jd"], t["V"])
            if not np.isfinite(np.r_[ra, dec, v]).all():
                raise ValueError(f"Non-finite ephemeris samples for {key} on {day}")
            ras.append(ra); decs.append(dec); mags.append(v)
            centers_ra.append(r0[0]); centers_dec.append(d0[0])
            ref_mag.append(float(np.interp(ref_jd, t["datetime_jd"], t["V"])))
        ras, decs, mags = map(np.array, (ras, decs, mags))
        targets = coord(ras, decs).transform_to(frame[None, :])
        centers = coord(centers_ra, centers_dec)[:, None].transform_to(frame[None, :])
        alt = targets.alt.deg; sep = targets.separation(moon[None, :]).deg
        center_alt = centers.alt.deg; center_sep = centers.separation(moon[None, :]).deg
        dark = -18. - sun_alt
        dark_intervals = intervals(seconds, dark[None, :], lambda q: q[0] > 0)
        dark_hours = sum(b - a for a, b in dark_intervals) / 3600
        moon_low = config["moon_sep_deg"] + 2.122
        moon_high = config["bright_moon_sep_deg"] + 2.122
        for j, key in enumerate(ids):
            visible = intervals(seconds, np.array([dark, alt[j] - 30.]), lambda q: bool(np.all(q > 0)))
            series = np.array([dark, center_alt[j] - 32.8, moon_alt,
                               fraction - config["bright_moon_fraction"], center_sep[j] - moon_low, center_sep[j] - moon_high])
            eligible = intervals(seconds, series, lambda q: bool(q[0] > 0 and q[1] > 0 and
                                (q[2] <= 0 or (q[5] > 0 if q[3] >= 0 else q[4] > 0))))
            two_hour = [[a, b] for a, b in eligible if b - a >= 7200]
            start_ranges = [[a, b - 7200] for a, b in two_hour]
            dark_indices = np.flatnonzero(dark > 0)
            peak = int(dark_indices[np.argmax(alt[j, dark_indices])]) if len(dark_indices) else None
            record = {
                "date": day, "object_id": key, "targetname": tracks[key]["targetname"],
                "reference_beijing": iso_bj(ref), "V_midnight": round(ref_mag[j], 4),
                "ra_midnight_deg": round(float(centers_ra[j]), 6), "dec_midnight_deg": round(float(centers_dec[j]), 6),
                "astronomical_dark_h": round(dark_hours, 6),
                "dark_alt30_h": round(sum(b - a for a, b in visible) / 3600, 6),
                "night_max_alt_deg": round(float(alt[j, peak]), 4) if peak is not None else None,
                "night_max_alt_at_beijing": iso_bj(seconds[peak]) if peak is not None else None,
                "center_moon_eligible_h": round(sum(b - a for a, b in eligible) / 3600, 6),
                "longest_center_moon_window_h": round(max((b - a for a, b in eligible), default=0.) / 3600, 6),
                "has_continuous_2h_window": bool(two_hour), "continuous_2h_window_count": len(two_hour),
                "dark_alt30_intervals": interval_records(visible),
                "center_moon_eligible_intervals": interval_records(eligible),
                "continuous_2h_windows": interval_records(two_hour),
                "two_hour_start_ranges": interval_records(start_ranges),
            }
            if day in monthly_dates:
                monthly.append(record)
        if day in args.night_dates:
            sel = np.flatnonzero((seconds >= noon.timestamp() + 4 * 3600) &
                                 (seconds <= noon.timestamp() + 20 * 3600) &
                                 (np.remainder(seconds - noon.timestamp(), 300) == 0))
            rounded = lambda x: np.round(np.asarray(x), 5).tolist()
            night_curves.append({
                "date": day, "utc": [iso_utc(s) for s in seconds[sel]],
                "beijing": [iso_bj(s) for s in seconds[sel]], "unix_utc": seconds[sel].tolist(),
                "sun_alt_deg": rounded(sun_alt[sel]), "moon_alt_deg": rounded(moon_alt[sel]),
                "moon_illuminated_fraction": rounded(fraction[sel]),
                "targets": {key: {"targetname": tracks[key]["targetname"],
                                  "ra_deg": rounded(ras[j, sel]), "dec_deg": rounded(decs[j, sel]),
                                  "V": rounded(mags[j, sel]), "alt_deg": rounded(alt[j, sel]),
                                  "moon_sep_deg": rounded(sep[j, sel]),
                                  "field_center_ra_deg": float(centers_ra[j]), "field_center_dec_deg": float(centers_dec[j]),
                                  "field_center_alt_deg": rounded(center_alt[j, sel]),
                                  "field_center_moon_sep_deg": rounded(center_sep[j, sel])}
                            for j, key in enumerate(ids)},
            })
        count = sum(r["has_continuous_2h_window"] for r in monthly if r["date"] == day)
        print(f"{day}: {len(ids)} targets; {dark_hours:.2f} h astronomical darkness; {count} representative 2 h geometrically eligible fields", flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / "special_visibility_monthly.json", {"schema_version": 1, "provenance": provenance,
                "rules": rules, "representative_nights": monthly_dates, "records": monthly})
    interval_keys = ("dark_alt30_intervals", "center_moon_eligible_intervals", "continuous_2h_windows", "two_hour_start_ranges")
    table = [{k: json.dumps(v, separators=(",", ":"), ensure_ascii=False) if k in interval_keys else v
              for k, v in r.items()} for r in monthly]
    with (args.output / "special_visibility_monthly.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0])); writer.writeheader(); writer.writerows(table)
    write_json(args.output / "special_night_tracks.json", {"schema_version": 1, "provenance": provenance,
                "rules": rules, "nights": night_curves})
    print(json.dumps({"output": str(args.output), "monthly_rows": len(monthly), "night_curve_count": len(night_curves),
                      "objects": len(ids), "representative_dates_only": True}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
