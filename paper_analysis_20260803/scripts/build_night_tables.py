#!/usr/bin/env python3
"""Build the frozen per-night accounting and unknown-stage tables."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path, required: bool = True) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    frame = pd.read_csv(path, dtype={"night": "string", "field_id": "string"})
    if "night" in frame:
        frame["night"] = frame["night"].astype("string").str.zfill(8)
    if "field_id" in frame:
        frame["field_id"] = frame["field_id"].astype("string").str.zfill(4)
    return frame


def numeric(frame: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
    if column not in frame:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").fillna(default)


def bool_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(False, index=frame.index)
    values = frame[column]
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y", "t"})


def per_night_counts(frame: pd.DataFrame, value_name: str) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=["night", value_name])
    return frame.groupby("night", as_index=False).size().rename(columns={"size": value_name})


def multiset_plan_match(raw: pd.DataFrame, plans: pd.DataFrame) -> pd.DataFrame:
    rows = []
    nights = sorted(set(raw.get("night", [])) | set(plans.get("night", [])))
    for night in nights:
        actual = Counter(
            raw.loc[raw["night"] == night, "field_id"].dropna().astype(str)
        )
        planned = Counter(
            plans.loc[plans["night"] == night, "field_id"].dropna().astype(str)
        )
        matched = sum(min(count, planned.get(field, 0)) for field, count in actual.items())
        actual_n = sum(actual.values())
        plan_n = sum(planned.values())
        rows.append(
            {
                "night": night,
                "plan_n": plan_n,
                "plan_matched_raw_n": matched,
                "raw_not_in_plan_n": actual_n - matched,
                "plan_not_acquired_n": plan_n - matched,
                "acquired_frame_plan_compliance": matched / actual_n if actual_n else np.nan,
                "planned_frame_realization": matched / plan_n if plan_n else np.nan,
            }
        )
    return pd.DataFrame(rows)


def distribution_summary(values: pd.Series) -> dict[str, float | int | None]:
    array = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    array = array[np.isfinite(array)]
    if not array.size:
        return {"n": 0, "median": None, "p16": None, "p84": None, "p90": None, "p95": None}
    q = np.percentile(array, [16, 50, 84, 90, 95])
    return {
        "n": int(array.size),
        "median": float(q[1]),
        "p16": float(q[0]),
        "p84": float(q[2]),
        "p90": float(q[3]),
        "p95": float(q[4]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory-dir", type=Path, required=True)
    parser.add_argument("--review-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--tracklet-stage-counts", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    inventory = args.inventory_dir.resolve()
    review_dir = args.review_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = json.loads(args.config.read_text(encoding="utf-8"))

    raw = read_csv(inventory / "raw_manifest.csv")
    engineering = read_csv(inventory / "raw_engineering_manifest.csv")
    l2 = read_csv(inventory / "l2_manifest.csv")
    plans = read_csv(inventory / "plan_rows.csv", required=False)
    products = read_csv(inventory / "nightly_products.csv")
    mask = read_csv(inventory / "mask_gaia_stage_counts.csv", required=False)
    tracklet_stage = read_csv(args.tracklet_stage_counts, required=False) if args.tracklet_stage_counts else pd.DataFrame()
    review_status = read_csv(review_dir / "review_and_mpc_status.csv")
    review_status["origin_night"] = (
        review_status["origin_night"].astype("string").str.zfill(8)
    )

    start = pd.Timestamp(config["observation_start"].replace("-", ""))
    end = pd.Timestamp(config["observation_end"].replace("-", ""))
    nights = pd.date_range(start, end, freq="D").strftime("%Y%m%d")
    status = pd.DataFrame({"night": nights})

    strict_raw = raw.loc[bool_series(raw, "strict_standard_science")].copy()
    strict_l2 = l2.loc[bool_series(l2, "strict_standard_catalog")].copy()
    for frame, name in (
        (strict_raw, "raw_science_n"),
        (engineering, "raw_all_mp_n"),
        (strict_l2, "l2_mp_n"),
    ):
        status = status.merge(per_night_counts(frame, name), on="night", how="left")
    status = status.merge(multiset_plan_match(strict_raw, plans), on="night", how="left")

    product_rename = {
        "known_all_rows": "known_predicted_n",
        "known_match1_rows": "known_match1_n",
        "known_mask15_rows": "known_mask15_n",
        "known_ades_rows": "known_ades_n",
        "unknown_links_rows": "unknown_n",
        "summary_matched_detections_total": "summary_known_mask15_n",
        "summary_mask_gaia_total_detections": "masked_catalog_detection_n",
        "summary_tracklets_total": "tracklet_n",
        "summary_links_total": "link_n",
        "summary_member_rows_total": "link_member_n",
        "summary_orbit_fit_ok": "orbit_fit_n",
        "summary_orbit_is_good": "orbit_is_good_n",
        "link_class_counts_all_non_asteroid": "link_all_non_known_n",
        "fit_ok_link_class_counts_all_non_asteroid": "orbit_fit_non_known_n",
        "is_good_link_class_counts_all_non_asteroid": "orbit_is_good_non_known_n",
        "requested_mask_gaia_known_detections_after_mask": "known_after_gaia_n",
        "requested_known_objects_with_ge_3_detections_total_objects": "known_ge3_object_n",
        "requested_known_objects_with_ge_3_detections_linked_objects": "known_ge3_linked_object_n",
        "review_real": "review_real_n",
        "submit_real": "submit_real_n",
    }
    keep = ["night"] + [column for column in product_rename if column in products]
    if "summary_exists" in products:
        keep.append("summary_exists")
    status = status.merge(
        products[keep].rename(columns=product_rename), on="night", how="left", validate="one_to_one"
    )
    if not mask.empty:
        mask_rename = {
            "prefilter_rows": "mag_flag_prefilter_detection_n",
            "gaia_survivor_rows": "gaia_survivor_n",
            "gaia_matched_rows": "gaia_removed_n",
            "status_ok_files": "gaia_ok_file_n",
            "parsed_rows": "gaia_log_file_n",
            "status_counts": "gaia_status_counts",
        }
        keep = ["night"] + [column for column in mask_rename if column in mask]
        status = status.merge(
            mask[keep].rename(columns=mask_rename), on="night", how="left", validate="one_to_one"
        )
    if not tracklet_stage.empty:
        tracklet_rename = {
            "parse_status": "tracklet_stage_state",
            "grouped_gaia_input_detection_n": "grouped_gaia_input_detection_n",
            "common_area_survivor_detection_n": "common_area_survivor_detection_n",
            "edge_shell_survivor_detection_n": "edge_shell_survivor_detection_n",
            "static_survivor_detection_n": "static_survivor_detection_n",
            "tracklet_n": "tracklet_log_n",
            "group_started_n": "tracklet_group_n",
            "skip_reason_counts": "tracklet_skip_reason_counts",
        }
        keep = ["night"] + [column for column in tracklet_rename if column in tracklet_stage]
        status = status.merge(
            tracklet_stage[keep].rename(columns=tracklet_rename),
            on="night",
            how="left",
            validate="one_to_one",
        )

    frozen_initial = review_status.loc[
        review_status["final_paper_status"].isin(
            ["retained_after_posthoc_audit", "rejected_posthoc"]
        )
    ]
    initial_audit = frozen_initial.groupby("origin_night", as_index=False).size().rename(
        columns={"origin_night": "night", "size": "audit_initial_n"}
    )
    status = status.merge(initial_audit, on="night", how="left")
    retained = review_status.loc[
        review_status["final_paper_status"].eq("retained_after_posthoc_audit")
    ]
    audit = retained.groupby("origin_night", as_index=False).size().rename(
        columns={"origin_night": "night", "size": "audit_real_n"}
    )
    status = status.merge(audit, on="night", how="left")

    count_columns = [
        column
        for column in status.columns
        if column.endswith("_n") and column not in {"plan_compliance_n"}
    ]
    for column in count_columns:
        status[column] = pd.to_numeric(status[column], errors="coerce").fillna(0).astype(int)

    excluded = config["quality_rules"]["exclude_from_primary_science"]
    quarantine = config["quality_rules"]["quarantine_unknown"]
    raw_without_l2 = config["quality_rules"]["raw_without_l2"]
    codes = []
    reasons = []
    for row in status.itertuples(index=False):
        row_codes = []
        row_reasons = []
        if row.night in excluded:
            row_codes.append("excluded_primary")
            row_reasons.append(excluded[row.night])
        if row.night in quarantine:
            row_codes.append("quarantine_unknown")
            row_reasons.append(quarantine[row.night])
        if row.night in raw_without_l2:
            row_codes.append("raw_without_l2")
            row_reasons.append(raw_without_l2[row.night])
        if row.raw_science_n == 0:
            row_codes.append("no_strict_raw")
        if row.raw_science_n > 0 and row.l2_mp_n == 0:
            row_codes.append("l2_missing")
        codes.append(";".join(row_codes) if row_codes else "quality_ok")
        reasons.append(";".join(row_reasons))
    status["quality_code"] = codes
    status["quality_reason"] = reasons
    status["primary_science_included"] = ~status["night"].isin(excluded)
    status["unknown_science_included"] = (
        status["primary_science_included"] & ~status["night"].isin(quarantine)
    )

    def known_state(row: pd.Series) -> str:
        source = products.loc[products["night"].eq(row["night"])]
        if source.empty:
            return "not_run"
        source = source.iloc[0]
        if bool_series(source.to_frame().T, "known_reply_exists").iloc[0]:
            return "reply_received"
        if row["known_ades_n"] > 0:
            return "ades_without_reply"
        if row["known_match1_n"] > 0:
            return "matched_without_ades"
        if "known_all_exists" in source and bool_series(source.to_frame().T, "known_all_exists").iloc[0]:
            return "zero_match_or_no_ades"
        return "not_run"

    def unknown_state(row: pd.Series) -> str:
        source = products.loc[products["night"].eq(row["night"])]
        if source.empty:
            return "not_run"
        source = source.iloc[0]
        if "unknown_reply_exists" in source and bool_series(source.to_frame().T, "unknown_reply_exists").iloc[0]:
            return "reply_received"
        if row["submit_real_n"] > 0:
            return "selected_without_reply"
        if row["unknown_n"] > 0:
            return "catalog_nonempty"
        if "unknown_links_exists" in source and bool_series(source.to_frame().T, "unknown_links_exists").iloc[0]:
            return "true_zero_catalog"
        return "not_run"

    status["known_mpc_state"] = status.apply(known_state, axis=1)
    status["unknown_mpc_state"] = status.apply(unknown_state, axis=1)
    status["l2_detection_n"] = status["night"].map(
        strict_l2.groupby("night")["n_sources"].apply(
            lambda values: pd.to_numeric(values, errors="coerce").sum()
        )
    ).fillna(0).astype(int)
    if "tracklet_log_n" in status:
        summary_available = status["summary_exists"].astype(str).str.lower().isin(
            {"1", "true", "yes", "y", "t"}
        )
        comparable = status["tracklet_stage_state"].isin(["ok", "zero_groups"]) & summary_available
        status["tracklet_count_log_minus_summary"] = np.where(
            comparable,
            status["tracklet_log_n"] - status["tracklet_n"],
            np.nan,
        )
    status["code_commit"] = config["production_algorithm_commit"]
    source_manifest_paths = [
        inventory / "raw_manifest.csv",
        inventory / "l2_manifest.csv",
        inventory / "nightly_products.csv",
    ]
    if args.tracklet_stage_counts and args.tracklet_stage_counts.exists():
        source_manifest_paths.append(args.tracklet_stage_counts)
    combined_hash = hashlib.sha256(
        "".join(sha256(path) for path in source_manifest_paths).encode("ascii")
    ).hexdigest()
    status["input_manifest_hash"] = combined_hash

    priority = [
        "night",
        "raw_science_n",
        "l2_mp_n",
        "plan_n",
        "plan_matched_raw_n",
        "known_predicted_n",
        "known_match1_n",
        "known_mask15_n",
        "known_ades_n",
        "known_mpc_state",
        "gaia_survivor_n",
        "grouped_gaia_input_detection_n",
        "edge_shell_survivor_detection_n",
        "static_survivor_detection_n",
        "tracklet_n",
        "link_n",
        "orbit_fit_n",
        "orbit_is_good_n",
        "unknown_n",
        "review_real_n",
        "submit_real_n",
        "audit_initial_n",
        "audit_real_n",
        "unknown_mpc_state",
        "quality_code",
        "quality_reason",
        "code_commit",
        "input_manifest_hash",
    ]
    ordered = priority + [column for column in status if column not in priority]
    status = status[ordered]
    status.to_csv(output / "night_status.csv", index=False)

    field_visits = (
        strict_raw.groupby("field_id", as_index=False)
        .agg(
            exposure_n=("file_name", "size"),
            night_n=("night", "nunique"),
            first_night=("night", "min"),
            last_night=("night", "max"),
            open_shutter_s=("exposure_s", lambda values: pd.to_numeric(values, errors="coerce").sum()),
        )
        .sort_values(["exposure_n", "field_id"], ascending=[False, True])
    )
    field_visits.to_csv(output / "field_visit_counts.csv", index=False)
    status[["night", "quality_code", "quality_reason", "primary_science_included", "unknown_science_included"]].to_csv(
        output / "quality_mask.csv", index=False
    )

    stage_map = {
        "l2_detection_n": "l2_source_detection",
        "mag_flag_prefilter_detection_n": "magnitude_flag_survivor_detection",
        "gaia_survivor_n": "gaia_survivor_detection",
        "grouped_gaia_input_detection_n": "grouped_gaia_input_detection",
        "common_area_survivor_detection_n": "common_area_survivor_detection",
        "edge_shell_survivor_detection_n": "edge_shell_survivor_detection",
        "static_survivor_detection_n": "reference_anchored_static_survivor_detection",
        "tracklet_n": "two_point_tracklet",
        "link_n": "shared_endpoint_linkage",
        "orbit_fit_n": "orbit_fit_ok_linkage",
        "orbit_is_good_n": "orbit_is_good_linkage",
        "orbit_fit_non_known_n": "fit_ok_all_non_known_linkage",
        "orbit_is_good_non_known_n": "is_good_all_non_known_linkage",
        "unknown_n": "post_known_catalog_linkage",
        "review_real_n": "human_review_marked_real_linkage",
        "submit_real_n": "submission_selected_linkage",
        "audit_initial_n": "frozen_initial_audit_set_linkage",
        "audit_real_n": "posthoc_retained_linkage",
    }
    stage_columns = [column for column in stage_map if column in status]
    stage_context = ["night", "quality_code", "unknown_science_included", "unknown_mpc_state"]
    if "tracklet_stage_state" in status:
        stage_context.append("tracklet_stage_state")
    stage = status[
        stage_context
        + stage_columns
    ].copy()
    stage["night_group"] = np.select(
        [
            stage["quality_code"].str.contains("excluded_primary"),
            stage["quality_code"].str.contains("quarantine_unknown"),
            stage["unknown_mpc_state"].eq("not_run"),
            stage["unknown_mpc_state"].eq("true_zero_catalog"),
        ],
        ["excluded", "quarantine", "not_run", "true_zero_catalog"],
        default="normal",
    )
    stage.to_csv(output / "unknown_stage_counts_by_night.csv", index=False)
    definitions = pd.DataFrame(
        [
            {
                "column": column,
                "stage": stage_name,
                "unit": "detection" if "detection" in stage_name else ("tracklet" if "tracklet" in stage_name else "linkage"),
                "definition_status": "as_executed",
            }
            for column, stage_name in stage_map.items()
            if column in stage
        ]
    )
    definitions.to_csv(output / "unknown_stage_definitions.csv", index=False)

    totals: dict[str, object] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "manifest_hash": combined_hash,
        "all_nights": {},
        "included_unknown_nights": {},
        "nightly_distributions": {},
    }
    included_stage = stage.loc[stage["unknown_science_included"]]
    for column in stage_columns:
        totals["all_nights"][column] = int(stage[column].sum())
        totals["included_unknown_nights"][column] = int(included_stage[column].sum())
        totals["nightly_distributions"][column] = {
            group: distribution_summary(frame[column])
            for group, frame in stage.groupby("night_group")
        }
    (output / "unknown_stage_counts_total.json").write_text(
        json.dumps(totals, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    initial_links = int((review_status["final_paper_status"].isin(["retained_after_posthoc_audit", "rejected_posthoc"])).sum())
    review_marked_real_links = int(status["review_real_n"].sum())
    submission_selected_links = int(status["submit_real_n"].sum())
    retained_links = int((review_status["final_paper_status"] == "retained_after_posthoc_audit").sum())
    taxonomy = pd.DataFrame(
        [
            {
                "stage": "initial_human_review",
                "taxonomy": "not_selected_by_initial_human_review",
                "count": max(int(status["unknown_n"].sum()) - review_marked_real_links, 0),
                "unit": "linkage",
                "evidence_status": "aggregate_only_no_subtype_labels",
            },
            {
                "stage": "submission_selection",
                "taxonomy": "review_real_reverted_before_submission",
                "count": review_marked_real_links - submission_selected_links,
                "unit": "linkage",
                "evidence_status": "row_identified_20260507_00001et",
            },
            {
                "stage": "posthoc_audit",
                "taxonomy": "manual_artifact_subtype_unrecorded",
                "count": initial_links - retained_links,
                "unit": "linkage",
                "evidence_status": "row_identified_subtype_missing",
            },
            {
                "stage": "posthoc_audit",
                "taxonomy": "retained_high_confidence",
                "count": retained_links,
                "unit": "linkage",
                "evidence_status": "row_identified",
            },
        ]
    )
    taxonomy.to_csv(output / "unknown_false_positive_taxonomy.csv", index=False)

    headline = {
        "strict_raw_frames": int(len(strict_raw)),
        "strict_raw_nights": int(strict_raw["night"].nunique()),
        "raw_fields": int(strict_raw["field_id"].nunique()),
        "open_shutter_hours": float(pd.to_numeric(strict_raw["exposure_s"], errors="coerce").sum() / 3600),
        "all_mp_fits": int(len(engineering)),
        "l2_catalogs": int(len(strict_l2)),
        "known_predictions": int(status["known_predicted_n"].sum()),
        "known_matches_1arcsec": int(status["known_match1_n"].sum()),
        "known_matches_1p5arcsec": int(status["known_mask15_n"].sum()),
        "unknown_catalog_linkages": int(status["unknown_n"].sum()),
        "human_review_marked_real_linkages": review_marked_real_links,
        "submission_selected_linkages": submission_selected_links,
        "initial_review_selected_linkages": initial_links,
        "posthoc_retained_linkages": retained_links,
        "posthoc_retained_detection_memberships": int(retained["n_obs"].sum()),
        "quality_mask_status": config["quality_rules"].get("status", "unknown"),
    }
    (output / "snapshot_summary.json").write_text(
        json.dumps(headline, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(headline, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
