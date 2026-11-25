#!/usr/bin/env python3
"""Reconstruct the exposure ledger and same-night geometric revisit coverage."""
from pathlib import Path

import astropy.units as u
import healpy as hp
import json
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

ROOT = Path(__file__).resolve().parents[1]
INPUT, TABLES, DATA = (ROOT / d for d in ("inputs", "tables", "figure_data"))
SITE = EarthLocation.from_geodetic(117.575*u.deg, 40.393*u.deg, 960*u.m)
iers.conf.auto_download = False
CUTOFF = pd.Timestamp("2026-10-15T16:00:00Z")
OLD_END = "20260715"
CORNERS = [f"l2_corner_{i}_{a}" for i in range(1, 5) for a in ("ra", "dec")]


def read(name, **kwargs):
    return pd.read_csv(INPUT / name, low_memory=False, **kwargs)


def mjd_utc(mjd):
    return pd.to_datetime(np.asarray(mjd), unit="D", origin="1858-11-17", utc=True)


def iso_time(series):
    # Camera software occasionally rounds x:xx:59.9995 to x:xx:60.000.
    # Preserve the carry to the next minute, tied to the independent MJD check.
    text = series.astype("string")
    carry = text.str.contains(r"T\d\d:\d\d:60(?:\.|$)", na=False)
    text = text.str.replace(r"(T\d\d:\d\d:)60", r"\g<1>59", regex=True)
    return pd.to_datetime(text, format="mixed", utc=True, errors="coerce") + pd.to_timedelta(carry.astype(int), unit="s")


def noon_night(utc):
    return (utc + pd.Timedelta(hours=-4)).dt.strftime("%Y%m%d")


def fits_name(stems):
    return stems.str.replace(r"_cat\.fits(?:\.gz)?$", ".fits", regex=True)


def aperture_altitude(frame):
    use = frame.mid_utc.notna() & frame.center_ra_deg.notna() & frame.center_dec_deg.notna()
    time = Time(frame.loc[use, "mid_utc"].to_numpy(dtype="datetime64[ns]"), scale="utc")
    altaz = AltAz(obstime=time, location=SITE, pressure=0*u.hPa)
    target = SkyCoord(frame.loc[use, "center_ra_deg"].to_numpy()*u.deg,
                      frame.loc[use, "center_dec_deg"].to_numpy()*u.deg)
    frame.loc[use, "field_altitude_deg"] = target.transform_to(altaz).alt.deg
    frame.loc[use, "sun_altitude_deg"] = get_sun(time).transform_to(altaz).alt.deg


def positions(ra, dec):
    return SkyCoord(np.asarray(ra, dtype=float)*u.deg, np.asarray(dec, dtype=float)*u.deg)


def pointing_metrics(frame):
    where = frame.wcs_verified.eq(True)
    maps = frame.loc[where, CORNERS].to_numpy(dtype=float).reshape(-1, 4, 2)
    corners = [positions(maps[:, i, 0], maps[:, i, 1]) for i in range(4)]
    sizes = [corners[i].separation(corners[(i+1)%4]).deg for i in range(4)]
    frame.loc[where, "wcs_max_side_deg"] = np.max(sizes, axis=0)
    source = positions(frame.center_ra_deg, frame.center_dec_deg)
    nominal = positions(frame.center_ra, frame.center_dec)
    telemetry = positions(pd.to_numeric(frame.tel_ra, errors="coerce")*15,
                          pd.to_numeric(frame.tel_dec, errors="coerce"))
    command = positions(pd.to_numeric(frame.obj_ra, errors="coerce")*15,
                        pd.to_numeric(frame.obj_dec, errors="coerce"))
    frame["wcs_nominal_offset_deg"] = source.separation(nominal).deg
    frame["wcs_telemetry_offset_deg"] = source.separation(telemetry).deg
    frame["wcs_command_offset_deg"] = source.separation(command).deg


def ledger():
    raw = read("raw_manifest.csv.gz", dtype={"night": str, "field_id": str, "sequence_id": str})
    l2 = read("l2_manifest.csv.gz", dtype={"night": str, "field_id": str, "sequence_id": str})
    server = read("server_headers.csv.gz", dtype={"directory_night": str})
    nm = ["night", "file_name"]
    remote = server.loc[server.kind.eq("raw")].rename(columns={"directory_night": "night", "source_file": "file_name"})
    l2["l2_file"] = l2.file_name
    l2["file_name"] = fits_name(l2.file_name)
    footprint = Table.read(INPUT / "survey_fov_footprints_with_visibility.fits").to_pandas()
    footprint.field_id = footprint.field_id.str.decode("ascii") if isinstance(footprint.field_id.iloc[0], bytes) else footprint.field_id
    raw_cols = ["night", "file_name", "field_id", "sequence_id", "size_bytes", "exposure_s"]
    old = raw[raw_cols].rename(columns={"exposure_s": "nominal_exposure_s", "size_bytes": "frozen_raw_bytes"})
    remote_cols = ["bytes", "error", "l1_present", "l2_present", "obs_date", "expsta", "exposure", "expreq",
                   "exptime", "filter", "naxis1", "naxis2", "tel_ra", "tel_dec", "obj_ra", "obj_dec", "expstop"]
    old = old.merge(remote[nm+remote_cols], on=nm, how="left", validate="one_to_one")
    old = old.rename(columns={"bytes": "current_raw_bytes", "filter": "raw_filter",
                              "l2_present": "current_l2_present"})
    wcs = server.loc[server.kind.eq("l2_wcs")].rename(
        columns={"directory_night": "night", "source_file": "l2_file", "mtime_ns": "current_l2_mtime_ns",
                 "bytes": "current_l2_bytes", "error": "l2_wcs_error"})
    l2 = l2.merge(wcs[["night", "l2_file", "current_l2_bytes", "current_l2_mtime_ns",
                       "l2_wcs_error", *CORNERS]], on=["night", "l2_file"], how="left", validate="one_to_one")
    frozen_mtime = pd.to_datetime(l2.mtime_utc, format="mixed", utc=True)
    # FITS readers in the 20260803 collector formatted timestamps to microseconds.
    mtime_ns = l2.current_l2_mtime_ns
    l2["wcs_verified"] = (l2.current_l2_bytes.eq(l2.size_bytes)
                          & (np.abs(mtime_ns - frozen_mtime.astype("int64")) < 3000)
                          & l2[CORNERS].notna().all(axis=1))
    old = old.merge(l2.drop(columns=["field_id", "sequence_id"]), on=nm, how="left", validate="one_to_one")
    old["l2_available"] = old.l2_file.notna()
    old["raw_start_utc"] = iso_time(old.obs_date)
    old["l2_start_utc"] = mjd_utc(old.mjd)
    old["l2_header_date_error_s"] = (iso_time(old.obs_time_utc) - old.l2_start_utc).dt.total_seconds()
    old["raw_l2_time_error_s"] = (old.raw_start_utc - old.l2_start_utc).dt.total_seconds()
    old["start_utc"] = old.raw_start_utc.combine_first(old.l2_start_utc)
    old["time_source"] = np.select([old.raw_start_utc.notna(), old.l2_start_utc.notna()],
                                  ["raw OBS_DATE", "frozen L2 MJD"], default="missing")
    old["actual_exposure_s"] = pd.to_numeric(old.exposure, errors="coerce").combine_first(old.exposure_s)
    old["mid_utc"] = old.start_utc + pd.to_timedelta(old.actual_exposure_s, unit="s")/2
    old["noon_night"] = noon_night(old.mid_utc)
    old["before_cutoff"] = old.mid_utc.lt(CUTOFF)
    old["exposure_key"] = old.file_name + "@" + old.start_utc.dt.floor("ms").astype("string")
    old = old.merge(footprint[["field_id", "center_ra", "center_dec"]], how="left", on="field_id", validate="many_to_one")
    old.center_ra_deg = old.center_ra_deg.fillna(old.center_ra)
    old.center_dec_deg = old.center_dec_deg.fillna(old.center_dec)
    with iers.earth_orientation_table.set(iers.IERS_A.open(str(INPUT / "paper_iers.all"))):
        aperture_altitude(old)
    pointing_metrics(old)
    old["duplicate_time_and_name"] = old.exposure_key.notna() & old.exposure_key.duplicated(keep=False)
    old["canonical_directory_copy"] = old.noon_night.eq(old.night)
    # Prefer the correctly named raw directory. A proven duplicate can supply the
    # sole frozen catalog when the canonical raw directory never received L2.
    better_l2 = old.loc[old.l2_available & old.canonical_directory_copy, "exposure_key"]
    canonical_keys = old.loc[old.canonical_directory_copy & old.raw_start_utc.notna(), "exposure_key"]
    old["l2_reassigned_from_copy"] = (old.l2_available & ~old.canonical_directory_copy
                                      & old.exposure_key.isin(canonical_keys)
                                      & ~old.exposure_key.isin(better_l2))
    old["noncanonical_duplicate"] = (old.duplicate_time_and_name & ~old.canonical_directory_copy
                                     & ~old.l2_reassigned_from_copy)
    old["historical_wcs_exclusion"] = old.night.eq("20260111")
    old["physical_time_conflict"] = old.field_altitude_deg.lt(0)
    old["abnormal_wcs_perimeter"] = old.wcs_verified.eq(True) & ~old.wcs_max_side_deg.between(2.80, 3.10)
    old["quality_exclusion"] = np.select(
        [~old.l2_available, old.historical_wcs_exclusion, old.mid_utc.isna(),
         old.physical_time_conflict, old.abnormal_wcs_perimeter, old.noncanonical_duplicate,
         old.noon_night.ne(old.night) & ~old.l2_reassigned_from_copy,
         old.raw_l2_time_error_s.abs().gt(1), old.l2_header_date_error_s.abs().gt(1)],
        ["missing_frozen_l2", "20260111_astrometric_WCS_offset", "missing_exposure_time",
         "field_below_horizon_at_header_time", "abnormal_WCS_detector_perimeter", "noncanonical_duplicate",
         "header_night_disagrees_with_directory",
         "raw_L2_start_disagreement_gt_1s", "L2_ISO_MJD_disagreement_gt_1s"],
        default="qualified")
    old["quality"] = old.quality_exclusion.eq("qualified") & old.before_cutoff
    old["observation_mode"] = "standard_MP"
    old["filter_reconciled"] = old.raw_filter.combine_first(old["filter"]).str.upper()
    old["historical_config"] = "20260803 snapshot; production baseline d2f0057; exact per-exposure version unavailable"
    night_status = read("nightly_products.csv.gz", dtype={"night": str})
    status = night_status[["night", "known_match1_exists", "known_mask15_exists", "unknown_links_exists", "unknown_links_rows"]]
    old = old.merge(status, on="night", how="left", validate="many_to_one")
    f = read("known_l4_frames.csv.gz", dtype={"night": str})
    f["file_name"] = fits_name(f.source_file)
    old = old.merge(f.drop(columns="source_file"), on=nm, how="left", validate="one_to_one")
    old["known_prediction_state"] = np.select(
        [old.prediction_rows.notna(), old.known_match1_exists.eq(True)],
        ["prediction_rows_present", "night_product_present_without_frame_rows"], default="night_product_absent")
    old.to_csv(TABLES / "exposure_ledger.csv.gz", index=False)
    return old, raw, footprint, server


def extension(server):
    frame = server.loc[server.kind.eq("raw") & server.directory_night.gt(OLD_END)].copy()
    frame["start_utc"] = iso_time(frame.obs_date).combine_first(iso_time(frame.expsta))
    frame["mid_utc"] = frame.start_utc + pd.to_timedelta(pd.to_numeric(frame.exposure, errors="coerce"), unit="s")/2
    frame["noon_night"] = noon_night(frame.mid_utc)
    frame["time_status"] = np.select([frame.mid_utc.isna(), frame.mid_utc.lt(CUTOFF)],
                                     ["unknown", "before_cutoff"], default="after_cutoff")
    frame["entered_extension"] = frame.time_status.eq("before_cutoff") & frame.strict_name.eq(True) & frame.noon_night.gt(OLD_END)
    keep = frame.loc[frame.mid_utc.isna() | frame.noon_night.le("20261015")]
    keep.to_csv(TABLES / "extension_exposure_inventory.csv.gz", index=False)
    group = keep.groupby(["directory_night", "time_status", "strict_name"], as_index=False).agg(
        raw_n=("source_file", "size"), in_extension_n=("entered_extension", "sum"),
        l1_available_n=("l1_present", "sum"), l2_available_n=("l2_present", "sum"),
        first_mid_utc=("mid_utc", "min"), last_mid_utc=("mid_utc", "max"),
        any_known_night_product=("known_night_product_present", "max"),
        any_unknown_night_product=("unknown_night_product_present", "max"))
    group.to_csv(TABLES / "extension_by_night.csv", index=False)
    return keep


def nominal_polygons(footprint, fields, nside, inclusive):
    pixels = {}
    idx = footprint.set_index("field_id")
    for field in sorted(fields):
        row = idx.loc[field]
        ra = [row[f"corner_ra_{i}"] for i in range(1, 5)]
        dec = [row[f"corner_dec_{i}"] for i in range(1, 5)]
        pixels[field] = hp.query_polygon(nside, hp.ang2vec(ra, dec, lonlat=True), inclusive=inclusive, fact=4)
    return pixels


def cadence_pointings(intervals):
    keys = ["noon_night", "field_id"]
    intervals["pointing_group"] = 1
    intervals["wcs_max_same_label_separation_deg"] = 0.0
    thresholds = (.25, .3, .5)
    alternative_counts = {t: {"groups": 0, "ge3_groups": 0, "split_label_groups": 0} for t in thresholds}
    for _, indices in intervals.groupby(keys).groups.items():
        group = intervals.loc[indices]
        if len(group) > 1:
            xyz = positions(group.center_ra_deg, group.center_dec_deg).cartesian.xyz.value.T
            distances = pdist(xyz)
            max_sep = np.rad2deg(2*np.arcsin(np.clip(distances.max()/2, 0, 1)))
        else:
            max_sep = 0.0
        intervals.loc[indices, "wcs_max_same_label_separation_deg"] = max_sep
        tree = linkage(distances, method="complete") if max_sep > min(thresholds) else None
        for limit in thresholds:
            if tree is None or max_sep <= limit:
                counts = np.array([len(group)])
            else:
                clusters = fcluster(tree, 2*np.sin(np.deg2rad(limit)/2), criterion="distance")
                counts = np.unique(clusters, return_counts=True)[1]
                if limit == .3:
                    intervals.loc[indices, "pointing_group"] = clusters
            row = alternative_counts[limit]
            row["groups"] += len(counts)
            row["ge3_groups"] += int((counts >= 3).sum())
            row["split_label_groups"] += int(len(counts) > 1)
    by_label = intervals.groupby(keys, as_index=False).agg(exposures=("mid_utc", "size"),
        first=("mid_utc", "min"), last=("mid_utc", "max"),
        max_wcs_center_separation_deg=("wcs_max_same_label_separation_deg", "first"),
        pointing_groups=("pointing_group", "nunique"))
    by_label.to_csv(TABLES / "field_label_night_cadence.csv", index=False)
    pd.DataFrame([{"max_center_separation_deg": t, **v} for t,v in alternative_counts.items()]).to_csv(
        TABLES / "cadence_pointing_tolerance.csv", index=False)
    return by_label


def run_coverage(old, raw, footprint):
    nominal = nominal_polygons(footprint, raw.field_id.unique(), 512, True)
    baseline_area = len(np.unique(np.concatenate(list(nominal.values()))))*hp.nside2pixarea(512, degrees=True)
    good = old.loc[old.quality].copy()
    intervals = good.sort_values(["noon_night", "field_id", "mid_utc"])
    cadence_pointings(intervals)
    by_group = intervals.groupby(["noon_night", "field_id", "pointing_group"])
    intervals["adjacent_min"] = by_group.mid_utc.diff().dt.total_seconds()/60
    cadence = by_group.agg(exposures=("file_name", "size"), first=("mid_utc", "min"), last=("mid_utc", "max")).reset_index()
    cadence["span_min"] = (cadence["last"]-cadence["first"]).dt.total_seconds()/60
    cadence.to_csv(TABLES / "field_night_cadence.csv", index=False)
    intervals[["noon_night", "field_id", "pointing_group", "file_name", "mid_utc", "adjacent_min"]].to_csv(
        TABLES / "adjacent_visit_intervals.csv.gz", index=False)
    activity = old.groupby("night", as_index=False).agg(raw=("file_name", "size"),
        l2=("l2_available", "sum"), qualified=("quality", "sum"))
    activity["month"] = activity.night.str[:6]
    activity.to_csv(DATA / "fig04_nightly_activity.csv", index=False)
    activity.groupby("month", as_index=False).agg(raw=("raw", "sum"), l2=("l2", "sum"),
        qualified=("qualified", "sum"), raw_nights=("night", "size")).to_csv(DATA / "fig04_monthly_activity.csv", index=False)
    cadence.exposures.clip(upper=3).value_counts().sort_index().rename_axis("visits_3_means_ge3").reset_index(name="field_nights").to_csv(
        DATA / "fig04_field_night_counts.csv", index=False)

    good = good.loc[good.wcs_verified.eq(True)]
    results = []
    for ns in (512, 1024):
        all_counts = np.zeros(hp.nside2npix(ns), dtype="uint16")
        thr_union = np.zeros(len(all_counts), dtype=bool)
        for _, night in good.groupby("noon_night"):
            coverage = []
            for row in night[CORNERS].to_numpy(dtype=float).reshape(-1, 4, 2):
                vectors = hp.ang2vec(row[:, 0], row[:, 1], lonlat=True)
                coverage.append(hp.query_polygon(ns, vectors, inclusive=False))
            pixels, count = np.unique(np.concatenate(coverage), return_counts=True)
            all_counts[pixels] += count.astype("uint16")
            thr_union[pixels[count >= 3]] = True
        pix_area = hp.nside2pixarea(ns, degrees=True)
        results.append({"nside": ns, "cohort": "qualified_with_verified_frozen_WCS",
                        "exposures": len(good), "geometric_once_deg2": np.count_nonzero(all_counts)*pix_area,
                        "geometric_ge3_same_night_deg2": thr_union.sum()*pix_area})
        if ns == 512:
            active = np.flatnonzero(all_counts)
            pd.DataFrame({"ring_pixel": active, "qualified_exposures": all_counts[active],
                          "ge3_same_night": thr_union[active]}).to_csv(DATA / "fig03_wcs_visits_icrs_nside512.csv.gz", index=False)
    pd.DataFrame(results).to_csv(TABLES / "geometric_coverage.csv", index=False)
    widths = {"adjacent": np.arange(0, 125, 5), "span_ge3": np.arange(0, 245, 10)}
    for key, series in [("adjacent", intervals.adjacent_min.dropna()),
                        ("span_ge3", cadence.loc[cadence.exposures.ge(3), "span_min"])]:
        bins = np.r_[widths[key], np.inf]
        counts = np.histogram(series, bins=bins)[0]
        assert counts.sum() == len(series) and series.gt(0).all()
        pd.DataFrame({"left_min": bins[:-1], "right_min": bins[1:], "n": counts}).to_csv(DATA / f"fig04_{key}.csv", index=False)
    print("baseline_nominal_area_deg2", baseline_area, "verified_quality_wcs", len(good))
    return baseline_area, cadence, intervals


if __name__ == "__main__":
    old, raw, footprint, server = ledger()
    print(old.groupby(["night", "quality_exclusion"]).size().loc[lambda x: x.index.get_level_values(1) != "qualified"].to_string())
    ext = extension(server)
    print("extension", ext.entered_extension.sum(), "boundary night",
          ext.loc[ext.noon_night.eq("20261015")].groupby(["strict_name", "time_status"]).size().to_dict())
    baseline_area, cadence, intervals = run_coverage(old, raw, footprint)
    interval = intervals.adjacent_min.dropna()
    span = cadence.loc[cadence.exposures.ge(3), "span_min"]
    q = old.loc[old.quality]
    assert q.exposure_key.is_unique
    l2 = old.loc[old.l2_available]
    ext_selected = ext.loc[ext.entered_extension]
    duplicate = old.loc[old.duplicate_time_and_name]
    summary = {
        "baseline_raw_file_inventory": int(len(raw)), "baseline_raw_directory_nights": int(raw.night.nunique()),
        "baseline_raw_nominal_shutter_hours": raw.exposure_s.sum()/3600, "baseline_raw_fields": int(raw.field_id.nunique()),
        "baseline_nominal_coverage_deg2": baseline_area,
        "frozen_l2_catalogs": len(l2), "frozen_l2_directory_nights": l2.night.nunique(),
        "quality_exposures": len(q), "quality_noon_nights": q.noon_night.nunique(),
        "quality_first_night": q.noon_night.min(), "quality_last_night": q.noon_night.max(),
        "quality_fields": q.field_id.nunique(), "quality_shutter_hours": q.actual_exposure_s.sum()/3600,
        "qualified_with_verified_wcs": int(q.wcs_verified.eq(True).sum()),
        "frozen_legacy_exposures_without_parseable_current_raw_time": int(old.raw_start_utc.isna().sum()),
        "relocated_sole_l2_copy": int(old.l2_reassigned_from_copy.sum()),
        "duplicate_key_file_rows": len(duplicate), "duplicate_physical_keys": duplicate.exposure_key.nunique(),
        "exclusion_counts": old.quality_exclusion.value_counts().to_dict(),
        "raw_wcs_scale_arcsec_median": float(np.sqrt(np.abs(q.cd1_1*q.cd2_2-q.cd1_2*q.cd2_1)).median()*3600),
        "field_nights": len(cadence),
        "field_nights_visits": {str(k):int(v) for k,v in cadence.exposures.clip(upper=3).value_counts().items()},
        "adjacent_pairs": len(interval), "adjacent_min_p16_p50_p84": interval.quantile([.16,.5,.84]).tolist(),
        "ge3_span_n": len(span), "ge3_span_min_p16_p50_p84": span.quantile([.16,.5,.84]).tolist(),
        "extension_strict_name": int(len(ext_selected)), "extension_strict_nights": ext_selected.noon_night.nunique(),
        "extension_strict_with_l2_now": int(ext_selected.l2_present.sum()),
        "extension_strict_filters": ext_selected["filter"].value_counts().to_dict(),
        "extension_boundary_counts": ext.loc[ext.noon_night.eq("20261015") & ext.strict_name.eq(True)].time_status.value_counts().to_dict()
    }
    (TABLES / "observing_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    point = old.loc[old.wcs_nominal_offset_deg.gt(1) | old.abnormal_wcs_perimeter]
    point[["night", "noon_night", "file_name", "quality", "center_ra_deg", "center_dec_deg",
           "wcs_nominal_offset_deg", "wcs_command_offset_deg", "wcs_telemetry_offset_deg",
           "wcs_max_side_deg", "n_wcs", "ra_rms_header", "dec_rms_header",
           "prediction_rows"]].to_csv(TABLES / "pointing_and_perimeter_audit.csv", index=False)
    old.loc[~old.quality, ["night", "noon_night", "file_name", "exposure_key",
                           "field_altitude_deg", "sun_altitude_deg", "wcs_max_side_deg", "quality_exclusion"]].to_csv(
        TABLES / "quality_exclusions.csv", index=False)
    old.groupby(["night", "quality_exclusion"], as_index=False).size().to_csv(
        TABLES / "quality_by_directory_night.csv", index=False)
