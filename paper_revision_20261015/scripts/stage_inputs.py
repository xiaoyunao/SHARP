#!/usr/bin/env python3
"""Extract compact analysis inputs from the checksum-preserved August snapshot."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd
from pyarrow import parquet as pq

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "inputs"
KNOWN_COLUMNS = ["night", "source_file", "asteroid_name", "asteroid_number", "epoch_mjd",
                 "pred_ra_deg", "pred_dec_deg", "pred_v_mag", "obs_ra_win_deg", "obs_dec_win_deg",
                 "x_win_px", "y_win_px", "obj_id", "mag_aper4", "magerr_aper4", "flux_aper4",
                 "fluxerr_aper4", "flag", "pred_geocentric_distance_au"]


def sha(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("august_figures", type=Path)
    parser.add_argument("replayed_recovery", type=Path)
    parser.add_argument("de432s", type=Path)
    args = parser.parse_args()
    snapshot, aug, replay = args.snapshot, args.august_figures, args.replayed_recovery
    sources = {}
    copies = [
        (snapshot / "config/snapshot.json", "snapshot_config.json"),
        *[(snapshot / f"snapshot/{rel}", Path(rel).name) for rel in [
            "inventory/raw_manifest.csv", "inventory/l2_manifest.csv", "tables/quality_mask.csv",
            "tables/night_status.csv", "coverage/coverage_summary.json",
            "inventory/nightly_products.csv", "provenance/survey_fov_footprints_with_visibility.fits",
            "provenance/production_file_hashes.csv", "scheduler/scheduler_mode_summary.json",
            "scheduler/scheduler_plan_realization_by_night.csv",
            "blinded_known/blinded_known_recovery_binned.csv",
            "blinded_known/blinded_known_summary.json",
            "jpl_identification/second_pass/jpl_second_pass_confirmations.csv",
            "review_sample/review_and_mpc_status.csv",
            "unknown_ades_time_correction/unknown_ades_midpoint_correction_179.csv",
            "unknown_ades_time_correction/unknown_ades_midpoint_correction_summary.json",
            "orbit_site_comparison/orbit_site_sensitivity_summary.json",
            "random_shift/known_random_shift_summary.json",
        ]],
        *[(aug / f"figure_data/{name}", name) for name in [
            "fig04_unique_object_orbits.csv", "fig01_sources.csv", "fig06_astrometry.csv",
        ]],
        (aug / "audit/image_source_manifest.csv", "image_source_manifest.csv"),
        (args.de432s, "de432s.bsp"),
    ]
    for src, name in copies:
        dest = INPUT / name
        if name in {"raw_manifest.csv", "l2_manifest.csv", "nightly_products.csv", "fig04_unique_object_orbits.csv"}:
            import gzip
            dest = dest.with_suffix(".csv.gz")
            with src.open("rb") as file, gzip.open(dest, "wb", compresslevel=5) as output:
                shutil.copyfileobj(file, output)
        else:
            shutil.copyfile(src, dest)
        sources[dest.name] = {"original": src.name, "original_sha256": sha(src)}

    inputs = {
        "known_matched": KNOWN_COLUMNS,
        "known_mask15": ["night", "source_file", "obj_id", "asteroid_number",
                         "asteroid_name", "epoch_mjd", "pred_ra_deg", "pred_dec_deg",
                         "obs_ra_win_deg", "obs_dec_win_deg", "mag_aper4",
                         "pred_v_mag", "x_win_px", "y_win_px", "flag"],
        "orbit_links": ["night", "linkage_id", "fit_ok", "is_good"],
        "orbit_obs_residuals": ["night", "linkage_id", "obs_key", "used"],
        "unknown_catalog": ["night", "trk_sub", "linkage_id", "n_obs", "fit_ok", "is_good",
                            "image_names", "objids", "mjds", "rms_arcsec", "max_arcsec",
                            "lin_rms_arcsec", "lin_speed_arcsec_per_day"],
    }
    for name, columns in inputs.items():
        src = snapshot / "snapshot/frozen_products" / f"{name}.parquet"
        pq.write_table(pq.read_table(src, columns=columns), INPUT / src.name, compression="zstd")
        sources[src.name] = {"original": src.name, "original_sha256": sha(src), "selected_columns": columns}

    src = replay / "blinded_known_recovery_by_object.csv"
    frame = pd.read_csv(src, dtype={"night": str}, low_memory=False)
    frame.to_parquet(INPUT / "blinded_known_recovery_by_object.parquet", compression="zstd", index=False)
    sources["blinded_known_recovery_by_object.parquet"] = {"original": src.name, "original_sha256": sha(src)}
    (INPUT / "original_hashes.json").write_text(json.dumps(sources, indent=2) + "\n")


if __name__ == "__main__":
    main()
