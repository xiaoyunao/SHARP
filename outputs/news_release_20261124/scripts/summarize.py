#!/usr/bin/env python3
"""Build the physical-exposure and known-object ledgers and annual HEALPix summaries."""
from collections import defaultdict
import csv
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import re

from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.utils import iers
import healpy as hp
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence"
PAPER = EVIDENCE / "paper_revision_20261015"
BASE_END = "20260715"
END = "20261123"
TZ = "Asia/Shanghai"
CORNERS = [f"l2_corner_{i}_{axis}" for i in range(1, 5) for axis in ("ra", "dec")]
MODES = ("MP", "twilight", "TJN", "followup")
NUMERIC_FIELD = re.compile(r"^MP_\d{4}$")
FILENAME = re.compile(r"(_cat)?\.fits(?:\.gz)?$")
FILTERED_SUFFIX = re.compile(r"_(?:[grizw]_)?\d{4,6}$", re.I)
iers.conf.auto_download = False


def read_stream(kind):
    with gzip.open(EVIDENCE / f"server_{kind}.jsonl.gz", "rt") as handle:
        return [json.loads(x) for x in handle]


def stem(filename):
    return FILENAME.sub("", Path(filename).name)


def timestamp(value):
    return pd.Timestamp(value, tz="UTC") if value else pd.NaT


def start_and_stop(header):
    start = timestamp(header.get("EXPSTA") or header.get("OBS_DATE") or header.get("DATE-OBS"))
    stop = timestamp(header.get("EXPEND"))
    elapsed = (stop - start).total_seconds() if pd.notna(stop) and pd.notna(start) else np.nan
    seconds = float(header.get("EXPTIME") or header.get("EXPOSURE") or np.nan)
    if pd.notna(start) and pd.isna(stop) and seconds > 0:
        stop = start + pd.to_timedelta(seconds, unit="s")
    consistent = not (elapsed <= 0 or (np.isfinite(elapsed) and abs(seconds-elapsed) > 1))
    return start, stop, seconds if seconds > 0 else elapsed, consistent


def noon_night(mid):
    if pd.isna(mid):
        return None
    return (mid.tz_convert(TZ) - pd.Timedelta(hours=12)).strftime("%Y%m%d")


def task_name(key):
    return FILTERED_SUFFIX.sub("", key.removeprefix("OBJ_"))


def main_mode(name):
    upper = name.upper()
    if upper.startswith(("MP_FU", "FOLLOW", "MP_FOLLOW")):
        return "followup"
    if upper.startswith(("TJ", "MP_TJ")):
        return "TJN"
    if upper.startswith(("TW_", "MP_TW_")):
        return "twilight"
    if NUMERIC_FIELD.fullmatch(upper) or upper.startswith("MP_"):
        return "MP"
    return "other"


def plan_mode(plan):
    mode = str(plan.get("mode", "")).lower()
    block = str(plan.get("block", "")).upper()
    if mode in {"ev", "am", "near_sun", "twilight"} or block in {"EV", "AM"}:
        return "twilight"
    if block == "FILL" or mode == "fill":
        return "MP"
    return None


def assignment(plans, frame):
    """Use observed field/filter/duration to label a matching planned activity."""
    field = frame["field"]
    filt = frame["filter"]
    options = [r for r in plans.get(field, []) if r["_mode"] and
               (not r.get("filter") or r["filter"].upper() == filt) and
               abs(float(r.get("exptime_s", frame["seconds"]))-frame["seconds"]) <= 0.11]
    if not options:
        return
    modes = {r["_mode"] for r in options}
    if len(modes) == 1:
        frame["task"] = modes.pop()
        frame["task_source"] = "same-night field, filter and duration in observing list"
    elif pd.notna(frame["start_utc"]):
        best = min(options, key=lambda r: abs((timestamp(r["start_utc"]) - frame["start_utc"]).total_seconds()))
        if abs((timestamp(best["start_utc"]) - frame["start_utc"]).total_seconds()) <= 1800:
            frame["task"] = best["_mode"]
            frame["task_source"] = "nearest actual start among matching same-night plan entries"


def dimensional_valid(corners):
    xyz = hp.ang2vec(corners[:, 0], corners[:, 1], lonlat=True)
    sides = np.rad2deg(np.arccos(np.einsum("ij,ij->i", xyz, np.roll(xyz, 1, axis=0)).clip(-1, 1)))
    return bool(np.isfinite(sides).all() and 2.80 <= sides.max() <= 3.10)


def ledger(nights):
    old = pd.read_csv(PAPER / "tables/exposure_ledger.csv.gz", low_memory=False)
    timed = ["start_utc", "mid_utc", "raw_start_utc"]
    for col in timed:
        old[col] = pd.to_datetime(old[col], format="mixed", utc=True)
    rows = []
    inputs = {entry["night"]: entry for entry in nights}
    old_manifest = {n: set(row["manifest"]) for n, row in inputs.items() if n <= BASE_END}
    for a in old.itertuples(index=False):
        n = str(a.night)
        key = stem(a.file_name)
        seconds = a.actual_exposure_s if pd.notna(a.raw_start_utc) else np.nan
        stop = a.start_utc + pd.to_timedelta(seconds, unit="s") if pd.notna(seconds) else pd.NaT
        rows.append(dict(directory_night=n, night=noon_night(a.mid_utc), stem=key,
                         field=task_name(key), task="MP", task_source="frozen MP cohort",
                         raw_path=f"/raw1/{n}/OBJ/{a.file_name}",
                         l2_path=a.path if isinstance(a.path, str) else None,
                         time_source=a.time_source, start_utc=a.start_utc, end_utc=stop,
                         mid_utc=a.mid_utc, seconds=seconds, filter=a.filter_reconciled,
                         filter_source="paper frozen filter reconciliation",
                         quality=bool(a.quality), quality_reason=a.quality_exclusion,
                         input_known=bool(pd.notna(a.prediction_rows) or (inputs.get(n, {}).get("known_exists") and
                                          key in old_manifest.get(n, set()))),
                         input_unknown=bool(inputs.get(n, {}).get("unknown_exists") and
                                            key in old_manifest.get(n, set())),
                         source="paper_e71f9f6a", invalid_directory=False,
                         ra=a.center_ra_deg, dec=a.center_dec_deg,
                         **{k: getattr(a, k) for k in CORNERS}))
    for entry in nights:
        n = entry["night"]
        accepted = entry["science_accepted"]
        used = set(accepted) if accepted else set(entry["manifest"])
        plans = defaultdict(list)
        for plan in entry["plan"]:
            field = str(plan.get("field_id", ""))
            field = f"MP_{field.zfill(4)}" if field.isdigit() else field
            plans[field].append({**plan, "_mode": plan_mode(plan)})
        for r in entry["rows"]:
            h = r.get("raw_header", {})
            l2 = r.get("l2_header", {})
            field = task_name(r["stem"])
            if re.match(r"^(?:test|start)(?:_|\d|$)", field, re.I):
                continue
            if str(h.get("IMAGETYP", "")).strip().upper() not in {"", "LIGHT", "0"}:
                continue
            start, end, seconds, consistent = start_and_stop(h)
            mid = start + (end - start) / 2
            physical = noon_night(mid)
            corners = np.array(l2.get("corners", [[np.nan, np.nan]] * 4))
            ra, dec = l2.get("center", (np.nan, np.nan))
            filter_name = re.search(r"_([grizw])_\d+$", r["stem"], re.I)
            actual_filter = (filter_name[1] if filter_name else h.get("FILTER") or l2.get("FILTER") or "")
            filter_source = "controller filename" if filter_name else ("raw FILTER" if h.get("FILTER") else
                             "L2 FILTER" if l2.get("FILTER") else "missing")
            l2_start = start_and_stop(l2)[0]
            mismatch = pd.notna(start) and pd.notna(l2_start) and abs((start-l2_start).total_seconds()) > 1
            l2_mjd = l2.get("MJD", np.nan)
            if np.isfinite(l2_mjd) and pd.notna(l2_start):
                l2_isot_delta = abs((Time(l2_mjd, format="mjd", scale="utc") -
                                     Time(l2_start.to_pydatetime(), scale="utc")).sec)
            else:
                l2_isot_delta = np.nan
            reasons = {
                "marked_invalid_directory": bool(entry["invalid"]),
                "raw_header_unreadable": "raw_error" in r,
                "missing_or_conflicting_shutter_time": pd.isna(mid) or not (seconds > 0) or not consistent,
                "directory_night_mismatch": physical != n,
                "missing_L2": r["l2"] is None,
                "missing_or_abnormal_WCS": not np.isfinite(corners).all() or not dimensional_valid(corners),
                "raw_L2_start_disagreement_gt_1s": mismatch,
                "L2_ISO_MJD_disagreement_gt_1s": l2_isot_delta > 1,
            }
            reason = next((k for k, v in reasons.items() if v), "qualified")
            row = dict(directory_night=n, night=physical, stem=r["stem"], field=field,
                       task=main_mode(field), task_source="actual target name", raw_path=r["raw"],
                       l2_path=r["l2"], time_source="raw EXPSTA/EXPEND" if h.get("EXPEND") else
                       "raw start and EXPTIME/EXPOSURE", start_utc=start, end_utc=end, mid_utc=mid,
                       seconds=seconds, filter=actual_filter.upper(), filter_source=filter_source,
                       quality=reason == "qualified", quality_reason=reason,
                       input_known=entry["known_exists"] and r["stem"] in used,
                       input_unknown=entry["unknown_exists"] and r["stem"] in used,
                       source="server_header", invalid_directory=bool(entry["invalid"]),
                       ra=ra, dec=dec, **dict(zip(CORNERS, corners.ravel())))
            if row["task"] == "MP" and physical == n:
                # Earlier trial nights have per-block executed counts in the on-site night audit.
                if n in {"20261111", "20261112", "20261113"} and NUMERIC_FIELD.fullmatch(field):
                    row["task"] = "twilight"
                    row["task_source"] = "November 11–13 site execution audit"
                elif plans:
                    assignment(plans, row)
            rows.append(row)
    data = pd.DataFrame(rows)
    possibly_visible = data.quality & data.source.eq("server_header")
    f = data.loc[possibly_visible]
    site = EarthLocation.from_geodetic(117.575 * u.deg, 40.3967 * u.deg, 960 * u.m)
    frame = AltAz(obstime=Time(f.mid_utc.to_numpy(dtype="datetime64[ns]"), scale="utc"),
                  location=site, pressure=0 * u.hPa)
    alt = SkyCoord(f.ra.to_numpy() * u.deg, f.dec.to_numpy() * u.deg, frame="icrs").transform_to(frame).alt.deg
    data.loc[f.index, "altitude_deg"] = alt
    bad = f.index[alt < 0]
    data.loc[bad, ["quality", "quality_reason"]] = [False, "WCS_center_below_horizon"]
    identity = data.start_utc.dt.floor("ms").astype("string")
    # One camera, one physical shutter opening. For missing raw times keep the inventory row.
    data["exposure_id"] = identity.fillna(data.directory_night + "/" + data.stem)
    data["canonical_directory"] = data.night.eq(data.directory_night)
    sort = data.sort_values(["invalid_directory", "quality", "canonical_directory", "input_known"],
                            ascending=[True, False, False, False], kind="stable")
    data["counted_raw"] = ~data.index.isin(sort.index[sort.exposure_id.duplicated()])
    data["main_sample"] = data.task.isin(MODES)
    data["used_in_analysis"] = data.input_known | data.input_unknown
    data.loc[data.invalid_directory, ["input_known", "input_unknown", "used_in_analysis"]] = False
    data["month"] = data.night.fillna(data.directory_night).str[:6]
    data["night_label"] = data.night.fillna(data.directory_night)
    return data


def norm(name):
    return re.sub("[^a-z0-9]", "", str(name).casefold())


def raw_key(number, name):
    return "N:" + number if number and number != "0" else "P:" + norm(name)


def known_matches(data):
    old = pd.read_parquet(PAPER / "inputs/known_matched.parquet",
                          columns=["night", "asteroid_name", "asteroid_number", "source_file", "epoch_mjd", "obj_id"])
    old.columns = ["directory_night", "name", "number", "source_file", "epoch", "objid"]
    old.number = old.number.fillna("")
    old["cohort"] = "paper"
    rows = [(rec["night"], name, str(number or ""), src, t, obj, "server")
            for rec in read_stream("known") for name, number, src, t, obj in rec["rows"]]
    new = pd.DataFrame(rows, columns=list(old.columns))
    det = pd.concat([old, new], ignore_index=True)
    det["stem"] = det.source_file.map(stem)
    exposure_cols = ["directory_night", "stem", "night", "task", "main_sample", "quality",
                     "counted_raw", "exposure_id", "start_utc", "mid_utc"]
    matched = det.merge(data[exposure_cols], on=["directory_night", "stem"], how="left",
                        validate="many_to_one", indicator=True)
    assert len(matched[matched.cohort.eq("paper") & matched._merge.ne("both")]) == 0
    quality = matched.quality.eq(True) & matched.counted_raw.eq(True)
    epoch_from_start = (matched.epoch - (matched.start_utc.astype("int64") / 1e9 / 86400 + 40587)) * 86400
    epoch_from_mid = (matched.epoch - (matched.mid_utc.astype("int64") / 1e9 / 86400 + 40587)) * 86400
    epoch_ok = np.minimum(abs(epoch_from_start), abs(epoch_from_mid)) <= 1
    # Frozen rows were already reconciled through the independent paper validation.
    q = matched.loc[quality].copy()
    q["epoch_verified"] = (epoch_ok | matched.cohort.eq("paper")).loc[q.index]
    q["original_identity"] = [raw_key(n, s) for n, s in zip(q.number, q.name)]
    paper = q[q.cohort.eq("paper")]
    assert (len(paper), paper.original_identity.nunique()) == (531612, 58281)
    timing = matched.loc[quality & ~epoch_ok & matched.cohort.eq("server"),
                         ["directory_night", "stem", "task"]].copy()
    timing["offset_from_start_s"] = epoch_from_start.loc[timing.index]
    timing["offset_from_mid_s"] = epoch_from_mid.loc[timing.index]
    timing.groupby(["directory_night", "task"], as_index=False).agg(
        measurements=("stem", "size"), exposures=("stem", "nunique"),
        median_offset_from_mid_s=("offset_from_mid_s", "median"),
        median_offset_from_start_s=("offset_from_start_s", "median"),
    ).to_csv(EVIDENCE / "legacy_known_epoch_audit.csv", index=False)
    info = dict(server_rows=len(new), server_matched_input=int(new.shape[0] -
                matched.loc[matched.cohort.eq("server"), "_merge"].ne("both").sum()),
                panel_rows_rejected_epoch=int((quality & ~epoch_ok & matched.cohort.eq("server")).sum()))
    return q, info


def classify_known(d):
    c = pd.read_csv(EVIDENCE / "astorb_elements.csv.gz", dtype={"number": str, "name": str},
                    keep_default_na=False)
    alias_evidence = json.loads((EVIDENCE / "jpl_2021vt99.json").read_text())["object"]
    assert "2021 VT99" in alias_evidence["fullname"]
    primary = alias_evidence["des"]
    names = {norm(s) for s in d.name.unique()} | {norm(primary)}
    numbers = set(d.number.unique()) - {""}
    norm_c = c.name.map(norm)
    c = c.loc[c.number.isin(numbers) | norm_c.isin(names)].copy()
    c["alias"] = norm_c.loc[c.index]
    c = c[(c.alias != "") | (c.number != "")]
    conflicts = c[c.alias.duplicated(keep=False) & c.alias.ne("")]
    safe_names = c[~c.alias.isin(conflicts.alias)].set_index("alias")
    num_data = c[c.number.ne("")].set_index("number", verify_integrity=True)
    identities = d[["number", "name", "original_identity"]].drop_duplicates().copy()
    identities["alias"] = identities.name.map(norm)
    alias_numbers = identities.alias.map(safe_names.number).fillna("")
    canonical_number = identities.number.where(identities.number.ne(""), alias_numbers)
    identities["identity"] = [raw_key(n, s) for n, s in zip(canonical_number, identities.name)]
    assert norm(primary) in safe_names.index
    identities.loc[identities.identity.eq(raw_key("", "2021 VT99")), "identity"] = raw_key("", primary)
    by_raw = identities[["original_identity", "identity"]].drop_duplicates()
    assert by_raw.original_identity.is_unique
    d = d.merge(by_raw, how="left", on="original_identity", validate="many_to_one")
    d = d.sort_values("epoch_verified", ascending=False, kind="stable").drop_duplicates(
        ["exposure_id", "identity"]).copy()
    jpl = {}
    for group in ("neo", "trojan"):
        with (EVIDENCE / f"jpl_{group}.json").open() as handle:
            ref = json.load(handle)
        assert ref["fields"] == ["pdes", "full_name", "class"]
        jpl[group] = {raw_key(x[0], "") if x[0].isdigit() else raw_key("", x[0]): x[2] for x in ref["data"]}
    objrows = []
    for identity, hits in d[d.epoch_verified & d.main_sample.eq(True)].groupby("identity", sort=False):
        key = identity[2:]
        orbit = num_data.loc[key] if identity.startswith("N:") and key in num_data.index else (
            safe_names.loc[key] if key in safe_names.index else None)
        row = dict(identity=identity, name=hits.name.iloc[0], detections=len(hits),
                   nights=hits.night.nunique(), jpl_neo=identity in jpl["neo"],
                   jpl_trojan=identity in jpl["trojan"],
                   jpl_special_class=jpl["neo"].get(identity) or jpl["trojan"].get(identity))
        if orbit is not None:
            a, e = float(orbit.a_au), float(orbit.ecc)
            q = a * (1 - e)
            categories = ["NEO (q < 1.3 au)", "Mars-crosser", "Inner belt / Hungaria",
                          "Main belt", "Outer belt / Hilda", "Jupiter Trojan region", "Distant (a ≥ 5.4 au)"]
            idx = next(i for i, yes in enumerate((q < 1.3, q < 1.666, a < 2.06, a < 3.28,
                                                   a < 5.05, a < 5.4, True)) if yes)
            row.update(a_au=a, q_au=q, ecc=e, inc_deg=float(orbit.inc_deg),
                       astorb_group=categories[idx], orbit_name=orbit["name"])
        else:
            row["astorb_group"] = "Comet" if re.match(r"^(?:[CPDXAI]/|\d+P(?:/|$))", hits.name.iloc[0]) else "Unclassified"
        objrows.append(row)
    return d, pd.DataFrame(objrows)


def coverage(data, nsides=(512, 1024)):
    cohorts = dict(paper=data[data.source.eq("paper_e71f9f6a")], full=data)
    results = []
    for nside in nsides:
        n = hp.nside2npix(nside)
        total = np.zeros(n, dtype=np.uint32)
        triple = np.zeros(n, dtype=bool)
        area = hp.nside2pixarea(nside, degrees=True)
        for scope, cohort in cohorts.items():
            section = cohort if scope == "paper" else cohort[cohort.source.ne("paper_e71f9f6a")]
            for _, night in section.groupby("night"):
                pix = [hp.query_polygon(nside, hp.ang2vec(row[:, 0], row[:, 1], lonlat=True),
                                        inclusive=False)
                       for row in night[CORNERS].to_numpy().reshape(-1, 4, 2)]
                pixels, count = np.unique(np.concatenate(pix), return_counts=True)
                total[pixels] += count.astype(np.uint32)
                triple[pixels[count >= 3]] = True
            summary = dict(cohort=scope, nside=nside, exposures=len(cohort),
                           area_once_deg2=np.count_nonzero(total)*area,
                           area_3_same_night_deg2=np.count_nonzero(triple)*area,
                           max_count=int(total.max()))
            results.append(summary)
            print(summary, flush=True)
        if nside == 512:
            active = np.flatnonzero(total)
            np.savez_compressed(EVIDENCE / "wcs_coverage.npz", ring_pixel=active,
                                count=total[active], same_night_ge3=triple[active], nside=nside)
    sums = pd.DataFrame(results)
    sums.to_csv(ROOT / "coverage_summary.csv", index=False)
    return sums


def counts(raw, q):
    return {"raw_exposures": len(raw), "raw_nights": raw.night_label.nunique(),
            "raw_timed_exposures": raw.start_utc.notna().sum(),
            "raw_shutter_hours_timed": raw.seconds.sum()/3600,
            "quality_exposures": len(q), "quality_nights": q.night.nunique(),
            "quality_shutter_hours": q.seconds.sum()/3600, "known_input_exposures": int(q.input_known.sum()),
            "unknown_manifest_exposures": int(q.input_unknown.sum()),
            "analysis_exposures": int(q.used_in_analysis.sum())}


def iso(x):
    return x.tz_convert(TZ).isoformat(timespec="milliseconds")


def main():
    nights = read_stream("exposures")
    assert nights[-1]["night"] == END
    data = ledger(nights)
    persist = data.drop(columns="altitude_deg")
    persist.to_csv(ROOT / "exposure_ledger.csv.gz", index=False)
    inventory = data.loc[data.counted_raw & data.main_sample]
    raw = inventory.loc[~inventory.invalid_directory]
    q = raw.loc[raw.quality]
    start = q.start_utc.min()
    assert len(q[q.source.eq("paper_e71f9f6a")]) == 39590
    assert q[CORNERS].notna().all().all()
    assert q.exposure_id.is_unique
    qother = data.loc[data.counted_raw & ~data.main_sample & data.quality]
    cov = coverage(q)
    detection, audit = known_matches(data)
    all_dedup, obj = classify_known(detection)
    archive_main = all_dedup.loc[all_dedup.main_sample.eq(True)]
    dedup = archive_main.loc[archive_main.epoch_verified]
    archive_asteroids = archive_main.loc[~archive_main.name.str.match(r"^(?:[CPDXAI]/|\d+P(?:/|$))")]
    asteroid = obj[obj.astorb_group.ne("Comet")]
    obj.to_csv(ROOT / "known_object_summary.csv.gz", index=False)
    dedup.to_csv(EVIDENCE / "known_measurement_crosswalk.csv.gz", index=False, columns=[
        "identity", "directory_night", "night", "source_file", "exposure_id", "task", "epoch", "objid", "cohort"])
    scopes = [("all_main", raw, q)] + [(t, raw[raw.task.eq(t)], q[q.task.eq(t)]) for t in MODES]
    month_rows = []
    for m in pd.period_range("2025-11", "2026-11", freq="M"):
        for name, rr, qq in scopes:
            month = m.strftime("%Y%m")
            month_rows.append(dict(month=str(m), task=name, **counts(rr[rr.month.eq(month)], qq[qq.month.eq(month)])))
    pd.DataFrame(month_rows).to_csv(ROOT / "monthly_summary.csv", index=False)
    area = cov.query("cohort == 'full' and nside == 1024").iloc[0]
    assert len(asteroid) == obj.astorb_group.ne("Comet").sum()
    ext = q.loc[q.source.eq("server_header")]
    supplement = all_dedup.loc[all_dedup.main_sample.eq(False) & all_dedup.epoch_verified]
    duplicates_n17 = data[data.directory_night.eq("20261117") & ~data.counted_raw]
    excluded_n17 = data[data.directory_night.eq("20261117") & data.counted_raw]
    summary = dict(
        data_through_observing_night=END, qualification_start=iso(start),
        qualification_last_end=iso(q.end_utc.max()), main_raw_last_end=iso(raw.end_utc.max()),
        main_raw_last_mid=iso(raw.mid_utc.max()), last_qualified_mid=iso(q.mid_utc.max()),
        any_science_last_end=iso(data.loc[data.counted_raw, "end_utc"].max()),
        collected_at_beijing=datetime.now(TZINFO).isoformat(timespec="seconds"),
        **{k: int(v) if isinstance(v, np.integer) else v for k, v in counts(raw, q).items()},
        area_once_deg2=float(area.area_once_deg2), area_3_same_night_deg2=float(area.area_3_same_night_deg2),
        asteroids=int(len(asteroid)), asteroid_measurements=int(asteroid.detections.sum()),
        asteroid_multinight=int(asteroid.nights.ge(2).sum()),
        archive_asteroids=int(archive_asteroids.identity.nunique()),
        archive_asteroid_measurements=len(archive_asteroids),
        known_comets=int(obj.astorb_group.eq("Comet").sum()),
        known_comet_measurements=int(obj.loc[obj.astorb_group.eq("Comet"), "detections"].sum()),
        asteroid_catalog_neo=int(asteroid.jpl_neo.sum()), asteroid_catalog_trojan=int(asteroid.jpl_trojan.sum()),
        asteroid_classes=asteroid.groupby("astorb_group").agg(objects=("identity", "size"),
                               measurements=("detections", "sum")).to_dict("index"),
        other_quality_exposures=len(qother), supplementary_known_measurements=len(supplement),
        supplementary_known_identities=int(supplement.identity.nunique()),
        supplementary_known_unshared_identities=int((~supplement.identity.isin(asteroid.identity)
                                                      ).groupby(supplement.identity).any().sum()),
        duplicate_nov17_rows=len(duplicates_n17), unique_invalid_nov17_rows=len(excluded_n17),
        raw_archive_first_directory=raw.directory_night.min(),
        raw_inventory_including_invalid=len(inventory),
        raw_invalid_nov17_main=int(excluded_n17.main_sample.sum()),
        raw_first_timed_noon_night=raw.night.dropna().min(),
        raw_archive_nights=raw.directory_night.nunique(),
        raw_noon_labels=raw.night.dropna().nunique(),
        raw_verified_same_night_count=raw.loc[raw.canonical_directory, "night"].nunique(),
        ext_exclusions=data.loc[data.source.eq("server_header") & data.counted_raw & data.main_sample, "quality_reason"]
                          .value_counts().to_dict(),
        qualified_filter_counts=q["filter"].value_counts().to_dict(),
        ext_filter_counts=ext["filter"].value_counts().to_dict(),
        analyzed_last_end=iso(q.loc[q.used_in_analysis, "end_utc"].max()),
        recent_night_counts={n: {k: v.item() if isinstance(v, np.generic) else v
                               for k, v in counts(raw[raw.night.eq(n)], q[q.night.eq(n)]).items()}
                             for n in ["20261118", "20261120", "20261121", "20261122", END]},
        nov18_tjn_fields=int(q.loc[q.night.eq("20261118") & q.task.eq("TJN"), "field"].nunique()),
        nov18_tjn_exposures=int((q.night.eq("20261118") & q.task.eq("TJN")).sum()),
        nov18_twilight_exposures=int((q.night.eq("20261118") & q.task.eq("twilight")).sum()),
        nov23_pipeline_completion=next(x["status"]["completed_at"] for x in nights if x["night"] == END),
        known_audit=audit,
    )
    (EVIDENCE / "annual_metrics.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    columns = ["metric", "value", "unit", "scope", "start_night", "end_night", "source"]
    units = {"exposures": "frames", "nights": "nights", "hours": "hours", "area": "deg²", "night_count": "nights", "labels": "nights",
             "asteroids": "objects", "measurements": "measurements", "neo": "objects", "trojan": "objects",
             "multinight": "objects", "comets": "objects", "rows": "frames"}
    annual = []
    def add(key, val, scope="all_main", source="physical exposure ledger", unit=None):
        unit = unit or next((v for k, v in units.items() if k in key), "value")
        annual.append([key, val, unit, scope, "20251115" if key.startswith("raw") else "20251116", END, source])
    for n, rr, qq in scopes:
        for key, val in counts(rr, qq).items():
            add(key, val, n)
    for name in ("area_once_deg2", "area_3_same_night_deg2"):
        add(name, summary[name], source="ICRS TPV polygon, ring HEALPix NSIDE=1024")
    for name in ("asteroids", "asteroid_measurements", "asteroid_multinight", "known_comets", "known_comet_measurements",
                 "asteroid_catalog_neo", "asteroid_catalog_trojan", "archive_asteroids", "archive_asteroid_measurements"):
        add(name, summary[name], source="deduplicated 1-arcsec L4 associations, ASTORB, JPL SBDB")
    for name, values in summary["asteroid_classes"].items():
        add("asteroids", values["objects"], scope=name, source="ASTORB orbital-element bins")
        add("asteroid_measurements", values["measurements"], scope=name, source="ASTORB orbital-element bins")
    for name in ("other_quality_exposures", "supplementary_known_measurements", "supplementary_known_identities",
                 "supplementary_known_unshared_identities", "raw_archive_nights", "raw_noon_labels",
                 "raw_verified_same_night_count", "duplicate_nov17_rows", "unique_invalid_nov17_rows",
                 "raw_inventory_including_invalid", "raw_invalid_nov17_main"):
        add(name, summary[name], scope="supplementary" if name.startswith(("other", "supplementary")) else "all_main",
            unit="objects" if name.endswith("identities") else "frames" if name in {"raw_inventory_including_invalid", "raw_invalid_nov17_main"} else None)
    with (ROOT / "annual_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(annual)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    from zoneinfo import ZoneInfo
    TZINFO = ZoneInfo(TZ)
    main()
