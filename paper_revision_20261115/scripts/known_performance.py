#!/usr/bin/env python3
"""Measure the frozen program's known associations without re-identifying objects."""
import json
from pathlib import Path

import astropy.units as u
import healpy as hp
import numpy as np
import pandas as pd
from astropy import constants
from astropy.coordinates import (CartesianRepresentation, BarycentricMeanEcliptic, SkyCoord,
                                 angular_separation, get_body_barycentric, get_sun,
                                 solar_system_ephemeris)
from astropy.time import Time
from astropy.utils import iers
from exposures import DATA, INPUT, SITE, TABLES, fits_name, read
from replay_frozen_recovery import binned_rows, object_key, wilson

RADIUS_COLS = ["dra_arcsec", "ddec_arcsec", "radial_arcsec"]
MATCH_KEYS = ["night", "source_file", "obj_id"]


def describe(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if not x.size:
        return {"n": 0}
    q = np.percentile(x, [0, 5, 10, 16, 25, 50, 75, 84, 90, 95, 100])
    labels = ["min", "p05", "p10", "p16", "p25", "median", "p75", "p84", "p90", "p95", "max"]
    return {"n": len(x), **dict(zip(labels, q)), "sigma_iqr": 0.7413*(q[6]-q[4])}


def angular(ra1, dec1, ra2, dec2):
    return angular_separation(*np.deg2rad([ra1, dec1, ra2, dec2])) * (180/np.pi)*3600


def canon_orbit_id(s):
    return s.str.replace("NUMBER:", "N:", regex=False).str.replace("NAME:", "S:", regex=False)


def with_residuals(frame):
    wrap = ((frame.obs_ra_win_deg-frame.pred_ra_deg+180) % 360) - 180
    frame["dra_arcsec"] = wrap*np.cos(np.deg2rad(frame.pred_dec_deg))*3600
    frame["ddec_arcsec"] = (frame.obs_dec_win_deg-frame.pred_dec_deg)*3600
    frame["radial_arcsec"] = angular(frame.pred_ra_deg, frame.pred_dec_deg,
                                     frame.obs_ra_win_deg, frame.obs_dec_win_deg)
    return frame


def read_known(ledger, name):
    d = pd.read_parquet(INPUT / f"{name}.parquet")
    d["object_key"] = object_key(d)
    d["file_name"] = fits_name(d.source_file.astype("string"))
    d = d.merge(ledger[["night", "file_name", "noon_night", "quality", "quality_exclusion",
                        "exposure_key", "mid_utc", "raw_start_utc", "l2_available"]],
                on=["night", "file_name"], how="left", validate="many_to_one")
    if d.quality.isna().any():
        raise ValueError("Unmapped frozen known catalog row")
    return with_residuals(d)


def exposure_pairs(frame):
    d = frame.sort_values(["object_key", "noon_night", "epoch_mjd", "exposure_key", "radial_arcsec"])
    duplicate_epoch = d.duplicated(["object_key", "noon_night", "epoch_mjd"], keep=False)
    duplicate_exposure = d.duplicated(["object_key", "exposure_key"], keep="first")
    d.loc[duplicate_epoch | duplicate_exposure].to_csv(TABLES / "motion_repeated_epochs.csv.gz", index=False)
    d = d.loc[~duplicate_exposure].copy()
    time_col = "epoch_mjd"
    keys = ["object_key", "noon_night"]
    cols = [time_col, "exposure_key", "pred_ra_deg", "pred_dec_deg", "obs_ra_win_deg", "obs_dec_win_deg",
            "source_file", "radial_arcsec"]
    prev = d.groupby(keys)[cols].shift()
    seconds = (d[time_col]-prev[time_col])*86400
    positive = seconds.gt(0)
    pairs = d.loc[positive, [*keys, "asteroid_number", "asteroid_name", *cols, "mag_aper4", "pred_v_mag"]].copy()
    pairs = pairs.join(prev.loc[positive].add_prefix("previous_"))
    pairs["separation_s"] = seconds.loc[positive]
    for label, ra, dec in [("predicted", "pred_ra_deg", "pred_dec_deg"),
                            ("measured", "obs_ra_win_deg", "obs_dec_win_deg")]:
        dist = angular(pairs[ra], pairs[dec], pairs["previous_"+ra], pairs["previous_"+dec])
        pairs[f"{label}_rate_arcsec_h"] = dist/(pairs.separation_s/3600)
    dra = (pairs.pred_ra_deg-pairs.previous_pred_ra_deg)*3600*np.cos(np.deg2rad(pairs.pred_dec_deg))
    ddec = (pairs.pred_dec_deg-pairs.previous_pred_dec_deg)*3600
    pairs["legacy_unwrapped_rate_arcsec_h"] = np.hypot(dra, ddec)/(pairs.separation_s/3600)
    pairs["ra_wrap"] = (pairs.pred_ra_deg-pairs.previous_pred_ra_deg).abs() > 180
    tail = pairs.loc[(pairs.filter(regex="_rate_arcsec_h$") > 1000).any(axis=1) | pairs.separation_s.lt(30)].copy()
    tail["audit_flag"] = np.select([tail.ra_wrap, tail.separation_s.lt(30),
                                   tail.predicted_rate_arcsec_h.gt(1000) & tail.measured_rate_arcsec_h.gt(1000)],
                                  ["RA_wrap", "very_short_pair", "prediction_and_measurement_fast"],
                                  default="inspect_identity_or_pair")
    tail.to_csv(TABLES / "motion_tail_audit.csv", index=False)
    return pairs, {"duplicate_epoch_rows_flagged": int(duplicate_epoch.sum()),
                   "duplicate_exposure_rows_removed": int(duplicate_exposure.sum()),
                   "nonpositive_interval_rows_omitted": int(seconds.le(0).sum()),
                   "legacy_gt1000_n": int(pairs.legacy_unwrapped_rate_arcsec_h.gt(1000).sum()),
                   "corrected_gt1000_n": int(pairs.predicted_rate_arcsec_h.gt(1000).sum()),
                   "ra_wrap_rows": int(pairs.ra_wrap.sum()),
                   "short_lt30s_pairs": int(pairs.separation_s.lt(30).sum())}


def vector_angles(a, b):
    return np.rad2deg(np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1),
                                np.einsum("ij,ij->i", a, b)))


def geometry(frame):
    mjds, index = np.unique(frame.epoch_mjd, return_inverse=True)
    t = Time(mjds, format="mjd", scale="utc")
    with iers.earth_orientation_table.set(iers.IERS_A.open(str(INPUT / "paper_iers.all"))), \
         solar_system_ephemeris.set(str(INPUT / "de432s.bsp")):
        position, _ = SITE.get_gcrs_posvel(t)
        higher = SITE.from_geodetic(SITE.lon, SITE.lat, SITE.height+1000*u.m)
        position_high, _ = higher.get_gcrs_posvel(t)
        zenith = (position_high-position).xyz.to_value(u.m).T
        zenith /= np.linalg.norm(zenith, axis=1)[:, None]
        # Aleph projects the heliocentric retarded asteroid minus the site.
        # Its observer HCRS vector uses this very same Astropy transform.
        obs_hcrs = SkyCoord(position, frame="gcrs", obstime=t).hcrs.cartesian.xyz.to_value(u.au).T
        same_kernel_sun = -obs_hcrs
        # Independently construct a retarded Sun-site vector from SPK barycentric states.
        observer = get_body_barycentric("earth", t) + position
        lighttime = np.linalg.norm(same_kernel_sun, axis=1)*u.au / constants.c
        sun_spk = (get_body_barycentric("sun", t-lighttime)-observer).xyz.to_value(u.au).T
        # Original error translated the finite-distance Sun to the solar-system barycentre.
        old_sun = get_sun(t).icrs.cartesian.xyz.to_value(u.au).T
    a = SkyCoord(frame.pred_ra_deg.to_numpy()*u.deg, frame.pred_dec_deg.to_numpy()*u.deg).cartesian.xyz.value.T
    sun = same_kernel_sun[index]
    elong = vector_angles(a, sun)
    independent = vector_angles(a, sun_spk[index])
    old = vector_angles(a, old_sun[index])
    z = zenith[index]
    target_alt = np.rad2deg(np.arcsin(np.clip(np.einsum("ij,ij->i", a, z), -1, 1)))
    sun_alt = np.rad2deg(np.arcsin(np.clip(np.einsum("ij,ij->i", sun, z)/np.linalg.norm(sun, axis=1), -1, 1)))
    bound = elong-np.abs(target_alt-sun_alt)
    assert bound.min() > -1e-9
    columns = ["night", "noon_night", "file_name", "source_file", "object_key", "asteroid_name",
               "epoch_mjd", "mid_utc", "raw_start_utc", "pred_ra_deg", "pred_dec_deg",
               "quality", "quality_exclusion"]
    geom = frame[columns].assign(
        elongation_deg=elong, legacy_elongation_deg=old, spk_retarded_elongation_deg=independent,
        sun_alt_deg=sun_alt, target_alt_deg=target_alt, altitude_bound_margin_deg=bound)
    minimum = pd.concat([geom.nsmallest(30, "elongation_deg").assign(audit="minimum_all_corrected"),
                         geom.loc[geom.quality].nsmallest(30, "elongation_deg").assign(audit="minimum_quality"),
                         geom.nsmallest(30, "legacy_elongation_deg").assign(audit="minimum_legacy")])
    minimum.to_csv(TABLES / "minimum_solar_elongations.csv", index=False)
    return geom


def quality_tables(d, mask15):
    rows = []
    for metric in RADIUS_COLS:
        for label, sub in [("baseline_1arcsec", d), ("qualified_1arcsec", d.loc[d.quality]),
                           ("qualified_1p5arcsec", mask15.loc[mask15.quality])]:
            rows.append({"sample": label, "metric": metric, **describe(sub[metric])})
    pd.DataFrame(rows).to_csv(TABLES / "astrometry_summary.csv", index=False)
    binned = []
    nightly = []
    qualified = d.loc[d.quality].copy()
    for metric in RADIUS_COLS:
        for night, group in qualified.groupby("noon_night"):
            nightly.append({"night": night, "metric": metric, "qualified_n": len(group), **describe(group[metric])})
        for radius, frame in [(1, qualified), (1.5, mask15.loc[mask15.quality])]:
            for low in np.arange(11, 22, .5):
                sample = frame.loc[frame.mag_aper4.ge(low) & frame.mag_aper4.lt(low+.5), metric]
                binned.append({"radius_arcsec": radius, "g_aper_low": low, "g_aper_high": low+.5,
                               "metric": metric, **describe(sample)})
    pd.DataFrame(nightly).to_csv(DATA / "fig05_astrometry_by_night.csv", index=False)
    pd.DataFrame([{"directory_night": night, "metric": col, **describe(sub[col])}
                  for night, sub in d.groupby("night") for col in RADIUS_COLS]).to_csv(
        TABLES / "baseline_astrometry_by_directory_night.csv", index=False)
    pd.DataFrame(binned).to_csv(DATA / "fig05_astrometry_vs_mag.csv", index=False)
    wide = mask15.loc[mask15.quality].copy()
    wide["matched_within_1"] = wide.radial_arcsec.le(1)
    by_bin = wide.assign(g_aper_bin=np.floor(wide.mag_aper4*2)/2).groupby("g_aper_bin").matched_within_1.agg(["count", "sum"]).reset_index()
    by_bin.to_csv(TABLES / "astrometry_radius_tail_vs_mag.csv", index=False)
    flux_ok = qualified.flag.eq(0) & qualified.magerr_aper4.gt(0) & qualified.flux_aper4.gt(0) & qualified.fluxerr_aper4.gt(0)
    phot = qualified.loc[flux_ok].assign(snr_aper=lambda x: x.flux_aper4/x.fluxerr_aper4)
    mag_rows = []
    for low in np.arange(11, 22, .5):
        sample = phot.loc[phot.mag_aper4.ge(low) & phot.mag_aper4.lt(low+.5)]
        for metric in ["magerr_aper4", "snr_aper"]:
            mag_rows.append({"g_aper_low": low, "g_aper_high": low+.5, "metric": metric, **describe(sample[metric])})
    pd.DataFrame(mag_rows).to_csv(DATA / "fig05_photometry_vs_mag.csv", index=False)
    return {"quality_mask15_gt1": int((~wide.matched_within_1).sum()), "quality_mask15_n": len(wide)}


def sample_tables(d, pairs):
    classes = read("fig04_unique_object_orbits.csv.gz", dtype={"asteroid_number": str, "catalog_number": str})
    class_keys = set(classes.object_key)
    mapping = d.object_key.drop_duplicates().to_frame()
    standard = canon_orbit_id(mapping.object_key)
    # August joins normalize provisional designations to lowercase without spaces.
    alias = "S:" + mapping.object_key.str.removeprefix("NAME:").str.replace(r"\s+", "", regex=True).str.lower()
    values = standard.where(standard.isin(class_keys), alias)
    if not values.isin(class_keys).all():
        unmatched = values.loc[~values.isin(class_keys)]
        raise ValueError(f"Identity differences: {unmatched.head().tolist()}, n={len(unmatched)}")
    d["orbit_key"] = d.object_key.map(dict(zip(mapping.object_key, values)))
    if d.orbit_key.nunique() != d.object_key.nunique():
        raise ValueError("Orbit identity mapping is not one-to-one")
    q = d.loc[d.quality].copy()
    orbits = classes.drop(columns=["first_mjd", "last_mjd", "detection_count", "night_count", "time_baseline_days"])
    objects = q.groupby("orbit_key", as_index=False).agg(detections=("epoch_mjd", "size"),
        nights=("noon_night", "nunique"), first=("epoch_mjd", "min"), last=("epoch_mjd", "max"))
    objects["span_days"] = objects["last"]-objects["first"]
    objects = objects.merge(orbits, left_on="orbit_key", right_on="object_key", validate="one_to_one")
    objects.orbit_class = objects.orbit_class.replace({"Distant / Centaur-TNO": "Distant (a ≥ 5.4 au)"})
    objects.to_csv(DATA / "fig07_objects_orbits_revisits.csv.gz", index=False)
    class_summary = []
    for scope, src in [("baseline", classes.rename(columns={"detection_count": "detections"})), ("quality", objects)]:
        for name, sub in src.groupby("orbit_class"):
            class_summary.append({"cohort": scope, "class": name.replace("Distant / Centaur-TNO", "Distant (a ≥ 5.4 au)"),
                                  "unique_objects": len(sub), "detections": int(sub.detections.sum())})
    pd.DataFrame(class_summary).to_csv(TABLES / "orbital_classes.csv", index=False)

    rare = classes.loc[classes.orbit_class.isin(["Distant / Centaur-TNO", "NEO (q < 1.3 AU)"])].copy()
    rare_keys = mapping.loc[values.isin(rare.object_key), "object_key"]
    rr = d.loc[d.object_key.isin(rare_keys)].copy()
    cols = ["orbit_key", "night", "source_file", "epoch_mjd", "quality", "pred_ra_deg", "pred_dec_deg",
            "obs_ra_win_deg", "obs_dec_win_deg", "pred_v_mag", "mag_aper4", "magerr_aper4", "flag",
            "dra_arcsec", "ddec_arcsec", "radial_arcsec"]
    rr[cols].to_csv(TABLES / "rare_object_detections.csv", index=False)
    stats = rr.groupby("orbit_key").agg(n_quality=("quality", "sum"), night_n=("night", "nunique"),
        max_radial_arcsec=("radial_arcsec", "max"), median_pred_v=("pred_v_mag", "median"),
        median_g_aper=("mag_aper4", "median"), min_epoch_mjd=("epoch_mjd", "min"), max_epoch_mjd=("epoch_mjd", "max"))
    motion = pairs.loc[pairs.object_key.isin(rare_keys)].copy()
    motion["rate_abs_delta"] = (motion.measured_rate_arcsec_h-motion.predicted_rate_arcsec_h).abs()
    motion["orbit_key"] = motion.object_key.map(dict(zip(mapping.object_key, values)))
    stats = stats.join(motion.groupby("orbit_key").agg(positive_same_night_pairs=("separation_s", "size"),
                       med_predicted_rate=("predicted_rate_arcsec_h", "median"),
                       max_meas_pred_rate_delta=("rate_abs_delta", "max")))
    stats["positive_same_night_pairs"] = stats.positive_same_night_pairs.fillna(0).astype(int)
    rare["perihelion_au"] = rare.a_au*(1-rare.ecc)
    rare.merge(stats, left_on="object_key", right_index=True, validate="one_to_one").to_csv(TABLES / "rare_object_audit.csv", index=False)
    columns = {"baseline_detected": len(d), "baseline_identified": d.object_key.nunique(),
               "qualified_detected": len(q), "qualified_identified": q.object_key.nunique(),
               "qualified_detected_nights": q.noon_night.nunique(), "qualified_object_nights": len(q[["object_key","noon_night"]].drop_duplicates()),
               "revisits": {col: describe(objects[col]) for col in ["detections", "nights", "span_days"]}}
    return columns


def retention(ledger):
    obj = pd.read_parquet(INPUT / "blinded_known_recovery_by_object.parquet")
    e = obj.loc[obj.eligible_ge3_and_rate].copy()
    e["span_min"] = (e.last_mjd-e.first_mjd)*1440
    exact = {"n": len(e), "k": int(e.strict_linked.sum())}
    assert exact == {"n": 102347, "k": 75899}
    exact["fraction"] = exact["k"]/exact["n"]
    exact["wilson_interval"] = wilson(exact["k"], exact["n"])
    rows = []
    for col, bins in [("detected_1p5arcsec_n", np.r_[3:10, np.inf]),
                      ("span_min", np.array([0, 20, 30, 40, 50, 60, 75, 90, 120, np.inf]))]:
        rows.extend(binned_rows(e, col, bins, col))
    pd.DataFrame(rows).to_csv(DATA / "fig06_retention_sampling.csv", index=False)
    frozen = read("blinded_known_recovery_binned.csv")
    frozen.to_csv(DATA / "fig06_retention_frozen_bins.csv", index=False)
    night_group = e.groupby("night").strict_linked.agg(["sum", "count"]).reset_index(names="night")
    states = read("night_status.csv", dtype={"night": str})
    night_group = night_group.merge(states[["night", "quality_code", "unknown_mpc_state", "summary_exists"]],
                                    on="night", validate="one_to_one")
    night_group["fraction"] = night_group["sum"]/night_group["count"]
    night_group["development_night"] = night_group.night.eq("20260220")
    night_group.to_csv(TABLES / "conditional_retention_by_night.csv", index=False)
    exact["development_night"] = {k: int(v) for k,v in night_group.loc[
        night_group.development_night, ["sum", "count"]].iloc[0].items()}
    exact["without_saved_stage_summary_n"] = int(night_group.loc[
        night_group.summary_exists.ne(True), "count"].sum())
    exact["historical_quarantine_n"] = int(night_group.loc[
        night_group.quality_code.eq("quarantine_unknown"), "count"].sum())
    e["bright_bin"] = e.median_measured_mag.between(10, 13, inclusive="left")
    e["slow_bin"] = e.first_last_rate_arcsec_per_hour.between(3, 5, inclusive="left")
    e["edge_lt500"] = e.median_detector_edge_px.lt(500)
    e.groupby(["bright_bin", "slow_bin", "edge_lt500"], as_index=False).strict_linked.agg(
        eligible_n="size", retained_n="sum").to_csv(TABLES / "conditional_retention_bright_slow_edge.csv", index=False)
    return exact


def candidate_tables():
    u = pd.read_parquet(INPUT / "unknown_catalog.parquet")
    raw = u.assign(image=u.image_names.str.split(";"), obj=u.objids.str.split(";")).explode(["image", "obj"])
    status = read("review_and_mpc_status.csv", dtype={"origin_night": str, "trk_sub": str})
    retained = status.loc[status.final_paper_status.eq("retained_after_posthoc_audit")].copy()
    retained = retained.merge(u, left_on=["origin_night", "linkage_id", "trk_sub"],
                              right_on=["night", "linkage_id", "trk_sub"], how="left", validate="one_to_one")
    jpl = read("jpl_second_pass_confirmations.csv", dtype={"night": str, "trk_sub": str})
    retained = retained.merge(jpl[["night", "trk_sub", "linkage_id", "second_pass_object_name",
                                  "numerically_confirmed_candidate", "second_pass_separation_arcsec",
                                  "speed_difference_arcsec_per_hour", "direction_difference_deg"]],
                              on=["night", "trk_sub", "linkage_id"], how="left", validate="one_to_one")
    retained["is_c2025_y1"] = retained.numerically_confirmed_candidate.eq(True)
    retained["frozen_identity_status"] = np.where(
        retained.is_c2025_y1, "ephemeris_consistent_C/2025_Y1", "unresolved_in_frozen_audit")
    midpoint = read("unknown_ades_midpoint_correction_179.csv", dtype={"night": str, "trk_sub": str})
    spans = midpoint.groupby(["night", "linkage_id"]).corrected_midpoint_mjd_utc.agg(lambda a: (a.max()-a.min())*1440)
    retained = retained.merge(spans.rename("arc_min_midpoint"), on=["night", "linkage_id"], validate="one_to_one")
    retained[["night", "trk_sub", "linkage_id", "n_obs_y", "rms_arcsec", "max_arcsec", "lin_rms_arcsec",
              "is_c2025_y1", "frozen_identity_status", "arc_min_midpoint", "mpc_identification", "second_pass_separation_arcsec",
              "speed_difference_arcsec_per_hour", "direction_difference_deg"]].to_csv(TABLES / "retained_linkages.csv", index=False)
    midpoint.to_csv(TABLES / "retained_measurements_midpoint.csv", index=False)
    night = read("night_status.csv", dtype={"night": str})
    stages = night[["night", "unknown_mpc_state", "unknown_n", "review_real_n", "submit_real_n",
                    "audit_initial_n", "audit_real_n", "quality_code", "quality_reason"]].copy()
    stages.to_csv(TABLES / "candidate_status_by_night.csv", index=False)
    return {"auto_links": len(u), "memberships": len(raw), "unique_detections": len(raw[["night", "image", "obj"]].drop_duplicates()),
            "fit_ok_and_is_good": int((u.fit_ok & u.is_good).sum()), "marked_real": int(stages.review_real_n.sum()),
            "selected": len(status), "retained": len(retained), "retained_measurements": len(midpoint),
            "retained_nights": retained.night.nunique(), "c2025y1_linkages": int(retained.is_c2025_y1.sum())}


def histogram(data, edges, file):
    pd.DataFrame({"lower": edges[:-1], "upper": edges[1:],
                  "count": np.histogram(data, bins=edges)[0]}).to_csv(DATA / file, index=False)


def main():
    ledger = pd.read_csv(TABLES / "exposure_ledger.csv.gz", dtype={"night": str, "noon_night": str}, low_memory=False)
    d = read_known(ledger, "known_matched")
    m = read_known(ledger, "known_mask15")
    pairs, audit = exposure_pairs(d.loc[d.quality])
    geom = geometry(d)
    stats = {"known": sample_tables(d, pairs), "astrometry": quality_tables(d, m),
             "retention": retention(ledger), "candidates": candidate_tables(),
             "motion_pair_audit": audit,
             "predicted_rate_arcsec_h": describe(pairs.predicted_rate_arcsec_h),
             "measured_rate_arcsec_h": describe(pairs.measured_rate_arcsec_h),
             "legacy_elong_deg": describe(geom.legacy_elongation_deg),
             "corrected_elong_all_deg": describe(geom.elongation_deg),
             "corrected_elong_quality_deg": describe(geom.loc[geom.quality, "elongation_deg"]),
             "spk_check_abs_difference_arcsec": describe((geom.elongation_deg-geom.spk_retarded_elongation_deg).abs()*3600),
             "minimum_altitude_bound_deg": geom.altitude_bound_margin_deg.min()}
    q = d.loc[d.quality]
    c = SkyCoord(q.obs_ra_win_deg.to_numpy()*u.deg, q.obs_dec_win_deg.to_numpy()*u.deg)
    ecl = c.transform_to(BarycentricMeanEcliptic(equinox=Time("J2000")))
    # A fixed obliquity rotation on unit vectors is recorded explicitly in the sky table.
    sky_pixel = hp.ang2pix(128, ecl.lon.deg, ecl.lat.deg, lonlat=True)
    pd.Series(sky_pixel).value_counts().sort_index().rename_axis("ring_pixel").reset_index(name="known_detections").to_csv(
        DATA / "fig03_known_ecliptic_nside128.csv", index=False)
    l2 = ledger.loc[ledger.quality, ["limiting_mag", "seeing_arcsec", "ccd_zeropoint", "ccd_zeropoint_rms"]].copy()
    stats["l2_quality_header_metrics"] = {col: describe(l2[col]) for col in l2}
    histogram(l2.limiting_mag, np.r_[0, np.arange(14, 24, .2), 99], "fig05_limiting_magnitude.csv")
    histogram(geom.loc[geom.quality, "elongation_deg"], np.linspace(0, 180, 61), "figA1_elongation.csv")
    histogram(pairs.predicted_rate_arcsec_h, np.r_[0, np.geomspace(.1, 10000, 76), np.inf], "figA1_predicted_rate.csv")
    histogram(pairs.measured_rate_arcsec_h, np.r_[0, np.geomspace(.1, 10000, 76), np.inf], "figA1_measured_rate.csv")
    assert stats["known"]["baseline_detected"] == 534780
    assert stats["known"]["baseline_identified"] == 58482
    assert list(stats["candidates"].values()) == [4762, 14299, 14159, 4762, 68, 67, 58, 179, 34, 6]
    (TABLES / "performance_summary.json").write_text(json.dumps(stats, indent=2) + "\n")


if __name__ == "__main__":
    main()
