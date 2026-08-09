#!/usr/bin/env python3
"""Freeze the locally audited 67 -> 58 unknown-link review sample.

This script normalizes the existing review exports into paper-facing tables.  It
does not infer MPC ingestion, identity, designation, or discovery status; those
fields remain explicitly unresolved until supported by external evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


KEY = ["night", "trk_sub", "linkage_id"]
ORBIT_LINK_COLUMNS = [
    "n_obs",
    "n_tracklets",
    "rms_arcsec",
    "med_arcsec",
    "max_arcsec",
    "a_au",
    "ecc",
    "inc_deg",
    "raan_deg",
    "argp_deg",
    "nu_deg",
    "best_v1_kms",
    "lin_rms_arcsec",
    "lin_speed_arcsec_per_day",
    "lin_dir_deg",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"night": "string", "trk_sub": "string"})


def finite_median(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--copy-example-gifs", action="store_true")
    args = parser.parse_args()

    source = args.source_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    manifest_path = source / "submitted_unknown_gif_manifest.csv"
    detections_path = source / "submitted_unknown_true_detections_photometry.csv"
    twilight_path = (
        source
        / "analysis_after_false_removal"
        / "submitted_unknown_link_twilight_summary.csv"
    )
    required = [manifest_path, detections_path, twilight_path]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing required inputs: {missing}")

    manifest = read_csv(manifest_path)
    detections = read_csv(detections_path)
    twilight = read_csv(twilight_path)
    for frame in (manifest, detections, twilight):
        frame["night"] = frame["night"].str.zfill(8)
        frame["linkage_id"] = pd.to_numeric(frame["linkage_id"], errors="raise").astype(int)

    initial_keys = manifest[KEY].drop_duplicates()
    retained_keys = detections[KEY].drop_duplicates()
    if len(manifest) != 67 or len(initial_keys) != 67:
        raise ValueError(f"expected 67 unique initial links, found {len(initial_keys)}")
    if len(retained_keys) != 58:
        raise ValueError(f"expected 58 retained links, found {len(retained_keys)}")
    if len(detections) != 179:
        raise ValueError(f"expected 179 retained member rows, found {len(detections)}")

    key_check = retained_keys.merge(initial_keys, on=KEY, how="left", indicator=True)
    if not (key_check["_merge"] == "both").all():
        raise ValueError("retained sample contains links absent from the 67-link manifest")

    rejected = manifest.merge(retained_keys, on=KEY, how="left", indicator=True)
    rejected = rejected.loc[rejected["_merge"] == "left_only"].drop(columns="_merge")
    if len(rejected) != 9:
        raise ValueError(f"expected 9 post-hoc rejects, found {len(rejected)}")
    rejected["review_status"] = "rejected_posthoc"
    rejected["reason_class"] = "manual_post_audit_artifact"
    rejected["reason_detail"] = "specific artifact subtype not recorded in source export"

    # A linkage-member row is not automatically a globally unique detection.
    detection_key_cols = ["catalog_path", "source_objID"]
    if not set(detection_key_cols).issubset(detections.columns):
        detection_key_cols = ["catalog_name", "source_objID"]
    detections["detection_key"] = (
        detections[detection_key_cols[0]].astype(str)
        + "::"
        + detections[detection_key_cols[1]].astype(str)
    )
    detections["link_member_row_id"] = np.arange(len(detections), dtype=int)
    detections["final_paper_status"] = "retained_after_posthoc_audit"

    first_values = (
        detections.sort_values(KEY + ["detection_index"])
        .groupby(KEY, as_index=False)
        .first()
    )
    aggregate = (
        detections.groupby(KEY, as_index=False)
        .agg(
            n_link_member_rows=("link_member_row_id", "size"),
            n_unique_detections=("detection_key", "nunique"),
            first_mjd=("MJD", "min"),
            median_mjd=("MJD", finite_median),
            last_mjd=("MJD", "max"),
            median_ra_deg=("RA_Win", finite_median),
            median_dec_deg=("DEC_Win", finite_median),
            median_mag_aper4=("Mag_Aper4", finite_median),
            median_mag_psf=("Mag_PSF", finite_median),
        )
    )
    for column in ORBIT_LINK_COLUMNS:
        if column in first_values.columns:
            aggregate = aggregate.merge(
                first_values[KEY + [column]], on=KEY, how="left", validate="one_to_one"
            )
    aggregate["speed_arcsec_per_hour"] = (
        pd.to_numeric(aggregate.get("lin_speed_arcsec_per_day"), errors="coerce") / 24.0
    )
    aggregate["final_paper_status"] = "retained_after_posthoc_audit"
    aggregate = aggregate.merge(
        twilight,
        on=KEY,
        how="left",
        validate="one_to_one",
        suffixes=("", "_legacy_twilight"),
    )
    aggregate["twilight_metrics_status"] = "legacy_site_coordinates_do_not_use_for_final_figure"

    review = manifest.copy()
    review = review.merge(
        aggregate[KEY + ["n_link_member_rows", "n_unique_detections"]],
        on=KEY,
        how="left",
        validate="one_to_one",
    )
    review["origin_night"] = review["night"]
    review["n_obs"] = review["n_link_member_rows"]
    review["submission_id"] = ""
    review["mpc_ingest_state"] = "unverified"
    review["mpc_identification"] = "unverified"
    review["provisional_designation"] = ""
    review["known_object_id"] = ""
    review["cross_night_group_id"] = "unresolved"
    retained_tuple = set(map(tuple, retained_keys[KEY].itertuples(index=False, name=None)))
    review["final_paper_status"] = [
        "retained_after_posthoc_audit"
        if tuple(row) in retained_tuple
        else "rejected_posthoc"
        for row in review[KEY].itertuples(index=False, name=None)
    ]
    review["notes"] = np.where(
        review["final_paper_status"].eq("rejected_posthoc"),
        "manual post-audit artifact; subtype absent from source export",
        "identity/designation/ingestion not yet reconciled with MPC/JPL",
    )
    review_columns = [
        "trk_sub",
        "origin_night",
        "linkage_id",
        "n_obs",
        "submission_id",
        "mpc_ingest_state",
        "mpc_identification",
        "provisional_designation",
        "known_object_id",
        "cross_night_group_id",
        "final_paper_status",
        "notes",
        "source_gif",
        "local_gif_name",
    ]

    outputs = {
        "unknown_high_confidence_links.csv": aggregate,
        "unknown_high_confidence_detections.csv": detections,
        "unknown_posthoc_rejects.csv": rejected,
        "review_and_mpc_status.csv": review[review_columns],
    }
    for name, frame in outputs.items():
        frame.to_csv(output / name, index=False)

    copied_gifs: list[str] = []
    if args.copy_example_gifs:
        gif_dir = output / "review_gifs"
        gif_dir.mkdir(exist_ok=True)
        example_names = {
            "20251121_000004r_link0598.gif",  # retained
            "20260207_00000Mi_link0178.gif",  # retained, dense field
            "20260529_00001gB_link0128.gif",  # post-hoc reject
            "20260530_00001gK_link0064.gif",  # post-hoc reject
        }
        for name in sorted(example_names):
            src = source / "gifs" / name
            if src.exists():
                shutil.copy2(src, gif_dir / name)
                copied_gifs.append(name)

    cross_link_duplicates = (
        detections.groupby("detection_key")[KEY]
        .apply(lambda frame: len(frame.drop_duplicates()))
        .astype(int)
    )
    summary = {
        "schema_version": "1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_directory": str(source),
        "initial_linkages": int(len(initial_keys)),
        "retained_linkages": int(len(retained_keys)),
        "posthoc_rejected_linkages": int(len(rejected)),
        "retained_link_member_rows": int(len(detections)),
        "retained_unique_detection_keys": int(detections["detection_key"].nunique()),
        "detection_keys_shared_across_links": int((cross_link_duplicates > 1).sum()),
        "retained_nights": int(retained_keys["night"].nunique()),
        "copied_example_gifs": copied_gifs,
        "identity_reconciliation_complete": False,
        "twilight_metrics_final": False,
        "input_hashes_sha256": {path.name: sha256(path) for path in required},
    }
    (output / "review_status_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    hash_lines = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "hashes.sha256":
            hash_lines.append(f"{sha256(path)}  {path.relative_to(output)}")
    (output / "hashes.sha256").write_text("\n".join(hash_lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
