"""Prepare an exploratory period sample from immutable annual compact inputs."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import angular_separation
from astropy.time import Time


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def best_window(frame, max_days=30.):
    """Densest fixed-duration window, selected by times only, never fit quality."""
    time = frame.t_mjd.to_numpy()
    ends = np.searchsorted(time, time + max_days, side="right")
    start = int(np.argmax(ends - np.arange(len(time))))
    return frame.iloc[start:ends[start]]


def prepare(annual, output):
    output.mkdir(parents=True, exist_ok=True)
    source = annual / "inputs/known_matched.parquet"
    ledger_path = annual / "tables/exposure_ledger.csv.gz"
    d = pd.read_parquet(source)
    ledger = pd.read_csv(ledger_path, usecols=["night", "file_name", "quality", "noon_night",
                                             "exposure_key", "mid_utc", "ccd_zeropoint_rms"])
    d["input_row"] = np.arange(len(d))
    d["night"] = d.night.astype(str)
    ledger["night"] = ledger.night.astype(str)
    d["file_name"] = d.source_file.str.replace(r"_cat\.fits(?:\.gz)?$", ".fits", regex=True)
    d = d.merge(ledger, on=["night", "file_name"], how="left", validate="many_to_one")
    if d.quality.isna().any():
        raise ValueError("Annual input contains rows absent from the exposure ledger")
    number = pd.to_numeric(d.asteroid_number, errors="coerce").astype("Int64")
    d["object_key"] = "N:" + number.astype(str)
    d["match_arcsec"] = np.rad2deg(angular_separation(
        *np.deg2rad([d.pred_ra_deg, d.pred_dec_deg, d.obs_ra_win_deg, d.obs_dec_win_deg]))) * 3600
    d["rejection"] = ""

    def reject(mask, reason):
        d.loc[mask & d.rejection.eq(""), "rejection"] = reason

    reject(~d.quality, "not_annual_quality")
    reject(number.isna(), "unnumbered_identity_deferred")
    reject(~d.flag.eq(0), "flag_nonzero")
    reject(~d.match_arcsec.le(1), "match_over_1arcsec_or_invalid")
    reject(~(d.flux_aper4.gt(0) & d.fluxerr_aper4.gt(0)
             & (d.flux_aper4 / d.fluxerr_aper4).ge(10)), "invalid_or_low_snr")
    reject(~(d.mag_aper4.between(9.5, 22) & d.magerr_aper4.gt(0)
             & d.magerr_aper4.le(.15)), "invalid_or_imprecise_magnitude")
    reject(~d.ccd_zeropoint_rms.between(.001, .1), "invalid_or_large_frame_zp_scatter")
    # A detection mapped to multiple catalog identities is not an independent sample.
    ambiguous = d.groupby(["night", "source_file", "obj_id"]).object_key.transform("nunique").gt(1)
    reject(ambiguous, "ambiguous_detection_identity")
    duplicate = d.loc[d.rejection.eq("")].duplicated(["object_key", "exposure_key"], keep=False)
    reject(d.index.isin(duplicate.index[duplicate]), "multiple_matches_same_object_exposure")
    counts = d.rejection.value_counts().rename_axis("first_rejection").reset_index(name="rows")
    counts.loc[counts.first_rejection.eq(""), "first_rejection"] = "accepted_exploratory"
    counts.to_csv(output / "selection_counts.csv", index=False)
    # Every rejected row is recoverable by its index into the hashed immutable input.
    d.loc[d.rejection.ne(""), ["input_row", "rejection"]].to_csv(
        output / "rejected_rows.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    good = d.loc[d.rejection.eq("")].copy()
    times = pd.to_datetime(good.mid_utc, utc=True, format="ISO8601").dt.tz_localize(None).to_numpy()
    good["t_mjd"] = Time(times, format="datetime64").mjd
    good["physical_night"] = good.noon_night.astype(int).astype(str)
    good["mag"] = good.mag_aper4
    good["mag_err"] = good.magerr_aper4
    good["fit_err"] = np.hypot(good.mag_err, .03)
    good["iso_flags_available"] = False
    good = good.sort_values(["object_key", "t_mjd", "exposure_key"])
    rows, windows = [], []
    for key, group in good.groupby("object_key", sort=True):
        window = best_window(group)
        night_counts = window.groupby("physical_night").size()
        eligible_nights = night_counts.index[night_counts.ge(3)]
        window = window.loc[window.physical_night.isin(eligible_nights)]
        n, nights = len(window), window.physical_night.nunique()
        rows.append(dict(object_key=key, annual_points=len(group), annual_nights=group.physical_night.nunique(),
                         annual_span_days=float(np.ptp(group.t_mjd)), window_points=n, window_nights=nights,
                         window_start_mjd=float(window.t_mjd.min()) if n else None,
                         window_end_mjd=float(window.t_mjd.max()) if n else None,
                         median_err_mag=float(window.mag_err.median()) if n else None,
                         pilot_eligible=n >= 40 and nights >= 4 and n >= nights + 12))
        if rows[-1]["pilot_eligible"]:
            windows.append(window)
    catalog = pd.DataFrame(rows).sort_values(["window_points", "object_key"], ascending=[False, True])
    catalog.to_csv(output / "sampling_catalog.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    columns = ["object_key", "asteroid_name", "night", "source_file", "obj_id", "input_row", "exposure_key",
               "physical_night", "t_mjd", "mag", "mag_err", "fit_err", "match_arcsec", "flag",
               "iso_flags_available", "ccd_zeropoint_rms"]
    selected = pd.concat(windows, ignore_index=True)[columns] if windows else pd.DataFrame(columns=columns)
    selected.to_csv(output / "eligible_window_points.csv.gz", index=False,
                    compression={"method": "gzip", "mtime": 0})
    manifest = dict(input_rows=len(d), accepted_rows=len(good), numbered_objects=len(catalog),
                    pilot_eligible_objects=int(catalog.pilot_eligible.sum()),
                    eligible_window_points=len(selected),
                    inputs={str(p.relative_to(annual)): sha256(p) for p in [source, ledger_path]},
                    selection="40+ points, 4+ nights (each 3+ points), n >= nights+12; densest <=30d window",
                    search_window_selection="timestamps only, before looking at magnitudes or periods",
                    error_floor_mag=.03, min_snr=10, max_magerr=.15, max_match_arcsec=1.,
                    max_frame_zp_scatter=.1, time_basis="UTC exposure midpoint MJD from annual ledger",
                    missing=["Flag_ISO", "Flag_ISO_Num", "distance and phase corrections", "light-time correction",
                             "differential comparison-star photometry", "image inspection", "external period novelty check"],
                    status="exploratory_only_not_new_period_measurements")
    (output / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ["input_rows", "accepted_rows", "numbered_objects",
                                              "pilot_eligible_objects", "eligible_window_points"]}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annual", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.annual, args.output)
