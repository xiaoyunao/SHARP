#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

from astropy.io import fits
from astropy.table import Table
from astropy.time import Time


PROCESSED = Path("/processed1")
TZ = ZoneInfo("Asia/Shanghai")


def clean_permid(value) -> str:
    s = str(value).strip()
    if s == "" or s.lower() in {"nan", "none", "null", "--", "-", "—"}:
        return ""
    try:
        return str(int(float(s)))
    except Exception:
        return ""


def obs_time_from_mjd(value) -> str:
    isot = Time(float(value), format="mjd", scale="utc").isot.replace(" ", "T")
    if "." in isot:
        return isot + "Z"
    return isot + ".000Z"


def assign_observing_night(date_obs: str | None) -> str | None:
    if not date_obs:
        return None
    dt = Time(date_obs, format="isot", scale="utc").to_datetime(timezone=TZ)
    if dt.hour < 12:
        dt = dt - timedelta(days=1)
    return dt.strftime("%Y%m%d")


def read_ades_duplicates() -> dict[tuple[str, str], list[str]]:
    seen: dict[tuple[str, str], list[str]] = defaultdict(list)
    for path in sorted(PROCESSED.glob("20*/L4/*_matched_asteroids_ades.psv")):
        night = path.parts[2]
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or line.startswith("#") or line.startswith("!"):
                    continue
                if "permID" in line and "obsTime" in line:
                    continue
                parts = [p.strip() for p in line.rstrip("\n").split("|")]
                if len(parts) < 11:
                    continue
                obj = parts[0] or parts[1]
                obs_time = parts[4]
                if obj and obs_time:
                    seen[(obj, obs_time)].append(night)
    return {k: v for k, v in seen.items() if len(v) > 1}


def load_matched_lookup(night: str, cache: dict[str, tuple[Table, dict[tuple[str, str], list[int]]]]):
    if night in cache:
        return cache[night]
    path = PROCESSED / night / "L4" / f"{night}_matched_asteroids.fits"
    table = Table.read(path, memmap=True)
    lookup: dict[tuple[str, str], list[int]] = defaultdict(list)
    for idx, row in enumerate(table):
        obj = clean_permid(row["number"])
        if obj:
            lookup[(obj, obs_time_from_mjd(row["epoch"]))].append(idx)
    cache[night] = (table, lookup)
    return cache[night]


def read_header_info(
    night: str,
    source_file: str,
    cache: dict[tuple[str, str], tuple[str | None, str | None]],
) -> tuple[str | None, str | None]:
    key = (night, source_file)
    if key in cache:
        return cache[key]
    path = PROCESSED / night / "L2" / source_file
    date_obs = None
    try:
        with fits.open(path, memmap=False) as hdul:
            hdr = hdul[1].header if len(hdul) > 1 else hdul[0].header
            date_obs = hdr.get("DATE-OBS") or hdr.get("DATE_OBS") or hdr.get("OBS_DATE")
    except Exception:
        date_obs = None
    observing_night = assign_observing_night(date_obs)
    cache[key] = (date_obs, observing_night)
    return cache[key]


def main() -> None:
    duplicates = read_ades_duplicates()
    print(f"duplicate_keys {len(duplicates)}")

    summary: Counter[str] = Counter()
    examples: dict[str, list[tuple]] = defaultdict(list)
    matched_cache: dict[str, tuple[Table, dict[tuple[str, str], list[int]]]] = {}
    header_cache: dict[tuple[str, str], tuple[str | None, str | None]] = {}

    for (obj, obs_time), nights in duplicates.items():
        unique_nights = sorted(set(nights))
        if len(unique_nights) > 2:
            bucket = "duplicate_in_more_than_two_files"
            summary[bucket] += 1
            if len(examples[bucket]) < 5:
                examples[bucket].append((obj, obs_time, unique_nights))
            continue
        if len(unique_nights) < 2:
            bucket = "duplicate_within_single_file_set"
            summary[bucket] += 1
            continue

        n1, n2 = unique_nights
        table1, lookup1 = load_matched_lookup(n1, matched_cache)
        table2, lookup2 = load_matched_lookup(n2, matched_cache)
        idxs1 = lookup1.get((obj, obs_time), [])
        idxs2 = lookup2.get((obj, obs_time), [])
        if len(idxs1) != 1 or len(idxs2) != 1:
            bucket = "missing_or_multiple_matched_rows"
            summary[bucket] += 1
            if len(examples[bucket]) < 5:
                examples[bucket].append((obj, obs_time, n1, len(idxs1), n2, len(idxs2)))
            continue

        row1 = table1[idxs1[0]]
        row2 = table2[idxs2[0]]
        src1 = str(row1["source_file"]).strip() if "source_file" in table1.colnames else ""
        src2 = str(row2["source_file"]).strip() if "source_file" in table2.colnames else ""
        date_obs1, obs_night1 = read_header_info(n1, src1, header_cache) if src1 else (None, None)
        date_obs2, obs_night2 = read_header_info(n2, src2, header_cache) if src2 else (None, None)

        same_source = bool(src1) and src1 == src2
        same_dateobs = date_obs1 is not None and date_obs1 == date_obs2
        same_observing_night = obs_night1 is not None and obs_night1 == obs_night2
        adjacent_nights = abs(int(n2) - int(n1)) == 1

        if same_source and same_dateobs and same_observing_night:
            bucket = "same_source_same_dateobs_same_observing_night"
            summary[bucket] += 1
            if adjacent_nights:
                summary["adjacent_night_overlap_same_source"] += 1
            else:
                summary["non_adjacent_same_source"] += 1
            if len(examples[bucket]) < 5:
                examples[bucket].append((obj, obs_time, n1, n2, src1, date_obs1, obs_night1))
        elif same_observing_night:
            bucket = "same_observing_night_but_different_source_or_dateobs"
            summary[bucket] += 1
            if len(examples[bucket]) < 5:
                examples[bucket].append((obj, obs_time, n1, n2, src1, src2, date_obs1, date_obs2, obs_night1, obs_night2))
        else:
            bucket = "other_pattern"
            summary[bucket] += 1
            if len(examples[bucket]) < 10:
                examples[bucket].append((obj, obs_time, n1, n2, src1, src2, date_obs1, date_obs2, obs_night1, obs_night2))

    print("summary")
    for key, value in summary.most_common():
        print(key, value)
    for bucket in sorted(examples):
        print(f"EXAMPLES {bucket}")
        for item in examples[bucket]:
            print(item)


if __name__ == "__main__":
    main()
