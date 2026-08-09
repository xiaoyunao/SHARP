#!/usr/bin/env python3
"""Generate Fig. 9: the frozen unknown-candidate accounting funnel.

The figure deliberately separates source detections, two-point tracklets,
single-night linkages, and review decisions.  Counts at those grains are not
joined by a continuous line.  The complete as-executed source-detection path
comes from ``unknown_stage_counts_by_night.csv``; independent frozen orbit,
post-known, membership, unique-detection, review, and false-positive tables
must reconcile before any figure is written.

The script reads only frozen analysis products.  It never scans production
directories and never substitutes headline constants for missing inputs.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

from figure_styles import (
    STATISTICS_COLORS,
    apply_statistics_style,
    save_pdf_png,
    style_statistics_axis,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = PROJECT_ROOT / "snapshot"
DEFAULT_STAGE_COUNTS = SNAPSHOT / "tables" / "unknown_stage_counts_by_night.csv"
DEFAULT_FALSE_POSITIVE_TAXONOMY = (
    SNAPSHOT / "tables" / "unknown_false_positive_taxonomy.csv"
)
DEFAULT_ORBIT_LINKS = SNAPSHOT / "frozen_products" / "orbit_links.parquet"
DEFAULT_UNKNOWN_LINKS = SNAPSHOT / "derived_unknown" / "unknown_all_links.csv"
DEFAULT_MEMBERSHIPS = SNAPSHOT / "derived_unknown" / "unknown_all_link_memberships.csv"
DEFAULT_UNIQUE_DETECTIONS = SNAPSHOT / "derived_unknown" / "unknown_unique_detections.csv"
DEFAULT_REVIEW_STATUS = SNAPSHOT / "review_sample" / "review_and_mpc_status.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "figures" / "fig09_unknown_funnel"
DEFAULT_FIGURE_DATA = PROJECT_ROOT / "figure_data" / "fig09_unknown_funnel.csv"


DETECTION_COLUMNS = [
    "l2_detection_n",
    "mag_flag_prefilter_detection_n",
    "gaia_survivor_n",
    "grouped_gaia_input_detection_n",
    "common_area_survivor_detection_n",
    "edge_shell_survivor_detection_n",
    "static_survivor_detection_n",
]
LINKAGE_COLUMNS = [
    "link_n",
    "orbit_fit_n",
    "orbit_is_good_n",
    "orbit_fit_non_known_n",
    "orbit_is_good_non_known_n",
    "unknown_n",
    "review_real_n",
    "submit_real_n",
    "audit_initial_n",
    "audit_real_n",
]
REQUIRED_STAGE_COLUMNS = {
    "night",
    "quality_code",
    "unknown_science_included",
    "unknown_mpc_state",
    "tracklet_stage_state",
    "night_group",
    *DETECTION_COLUMNS,
    "tracklet_n",
    *LINKAGE_COLUMNS,
}

UNKNOWN_MPC_STATES = {
    "reply_received",
    "selected_without_reply",
    "catalog_nonempty",
    "true_zero_catalog",
    "not_run",
}
TRACKLET_STAGE_STATES = {"ok", "zero_groups", "missing"}
NIGHT_GROUPS = {"normal", "true_zero_catalog", "not_run", "quarantine", "excluded"}
GROUP_ORDER = ["included", "primary_excluded", "unknown_quarantine"]
GROUP_STYLES = {
    "primary_excluded": {
        "color": "#c8c8c8",
        "hatch": "xx",
        "label": "Primary-science excluded",
    },
    "unknown_quarantine": {
        "color": "#9f9f9f",
        "hatch": "////",
        "label": "Unknown-branch quarantine",
    },
}


def parse_bool(values: pd.Series, name: str) -> pd.Series:
    normalized = values.astype("string").str.strip().str.lower()
    mapping = {
        "true": True,
        "1": True,
        "yes": True,
        "y": True,
        "t": True,
        "false": False,
        "0": False,
        "no": False,
        "n": False,
        "f": False,
    }
    invalid = normalized.isna() | ~normalized.isin(mapping)
    if invalid.any():
        raise ValueError(f"{name} contains invalid booleans: {values[invalid].head().tolist()}")
    return normalized.map(mapping).astype(bool)


def normalize_night(values: pd.Series, name: str) -> pd.Series:
    nights = values.astype("string").str.strip().str.replace(r"\.0$", "", regex=True).str.zfill(8)
    invalid = nights.isna() | ~nights.str.fullmatch(r"\d{8}", na=False)
    if invalid.any():
        raise ValueError(f"{name} contains invalid night keys: {values[invalid].head().tolist()}")
    return nights


def integer_counts(frame: pd.DataFrame, columns: Iterable[str], source: Path) -> None:
    for column in columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        invalid = values.isna() | (values < 0) | ~np.isclose(values, np.round(values))
        if invalid.any():
            examples = frame.loc[invalid, ["night", column]].head().to_dict("records")
            raise ValueError(f"{source}: {column} is not a non-negative integer count: {examples}")
        frame[column] = np.round(values).astype("int64")


def require_file(path: Path, purpose: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing frozen {purpose}: {path}")


def assert_equal(label: str, left: int, right: int) -> None:
    if int(left) != int(right):
        raise ValueError(f"Frozen-source reconciliation failed for {label}: {left:,} != {right:,}")


def assert_relation(
    frame: pd.DataFrame,
    left: str,
    right: str,
    relation: str,
    source: Path,
) -> None:
    if relation == ">=":
        valid = frame[left].ge(frame[right])
    elif relation == "==":
        valid = frame[left].eq(frame[right])
    else:  # pragma: no cover - internal programming guard
        raise ValueError(f"Unsupported relation: {relation}")
    if not valid.all():
        examples = frame.loc[~valid, ["night", left, right]].head().to_dict("records")
        raise ValueError(f"{source}: expected {left} {relation} {right}; examples: {examples}")


def read_stage_counts(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Frozen unknown stage table not found: {path}\n"
            "Run build_night_tables.py on the frozen manifests/products first. "
            "Fig. 9 will not infer missing stages from live data."
        )
    frame = pd.read_csv(path, dtype={"night": "string"}, low_memory=False)
    missing = sorted(REQUIRED_STAGE_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"{path} lacks required Fig. 9 columns: {', '.join(missing)}")

    frame["night"] = normalize_night(frame["night"], "stage-count night")
    if frame["night"].duplicated().any():
        raise ValueError(f"{path}: night must be unique")
    integer_counts(frame, [*DETECTION_COLUMNS, "tracklet_n", *LINKAGE_COLUMNS], path)
    frame["unknown_science_included"] = parse_bool(
        frame["unknown_science_included"], "unknown_science_included"
    )

    for column, allowed in [
        ("unknown_mpc_state", UNKNOWN_MPC_STATES),
        ("tracklet_stage_state", TRACKLET_STAGE_STATES),
        ("night_group", NIGHT_GROUPS),
    ]:
        values = frame[column].fillna("").astype(str).str.strip()
        unexpected = sorted(set(values) - allowed)
        if unexpected:
            raise ValueError(f"{path}: unexpected {column} values: {unexpected}")
        frame[column] = values

    quality = frame["quality_code"].fillna("").astype(str)
    frame["quality_group"] = np.select(
        [
            quality.str.contains("excluded_primary", regex=False),
            quality.str.contains("quarantine_unknown", regex=False),
        ],
        ["primary_excluded", "unknown_quarantine"],
        default="included",
    )
    inferred_included = frame["quality_group"].eq("included")
    mismatch = inferred_included.ne(frame["unknown_science_included"])
    if mismatch.any():
        examples = frame.loc[
            mismatch, ["night", "quality_code", "unknown_science_included"]
        ].head().to_dict("records")
        raise ValueError(f"quality_code and unknown_science_included disagree: {examples}")

    expected_night_group = np.select(
        [
            frame["quality_group"].eq("primary_excluded"),
            frame["quality_group"].eq("unknown_quarantine"),
            frame["unknown_mpc_state"].eq("not_run"),
            frame["unknown_mpc_state"].eq("true_zero_catalog"),
        ],
        ["excluded", "quarantine", "not_run", "true_zero_catalog"],
        default="normal",
    )
    group_mismatch = frame["night_group"].ne(expected_night_group)
    if group_mismatch.any():
        examples = frame.loc[
            group_mismatch,
            ["night", "quality_code", "unknown_mpc_state", "night_group"],
        ].head().to_dict("records")
        raise ValueError(f"night_group is inconsistent with state/quality precedence: {examples}")

    for left, right, relation in [
        ("l2_detection_n", "mag_flag_prefilter_detection_n", ">="),
        ("mag_flag_prefilter_detection_n", "gaia_survivor_n", ">="),
        ("gaia_survivor_n", "grouped_gaia_input_detection_n", ">="),
        # The production command used --skip-common-area.  A difference here
        # means the plotted no-op label would be false and must stop the build.
        ("grouped_gaia_input_detection_n", "common_area_survivor_detection_n", "=="),
        ("common_area_survivor_detection_n", "edge_shell_survivor_detection_n", ">="),
        ("edge_shell_survivor_detection_n", "static_survivor_detection_n", ">="),
        ("link_n", "orbit_fit_n", ">="),
        ("orbit_fit_n", "orbit_is_good_n", ">="),
        ("orbit_fit_n", "orbit_fit_non_known_n", ">="),
        ("orbit_fit_non_known_n", "orbit_is_good_non_known_n", ">="),
        ("orbit_fit_non_known_n", "unknown_n", ">="),
        ("orbit_is_good_non_known_n", "unknown_n", ">="),
        ("unknown_n", "review_real_n", ">="),
        ("review_real_n", "submit_real_n", ">="),
        ("submit_real_n", "audit_initial_n", "=="),
        ("audit_initial_n", "audit_real_n", ">="),
    ]:
        assert_relation(frame, left, right, relation, path)

    zero_catalog = frame["unknown_mpc_state"].eq("true_zero_catalog")
    not_run = frame["unknown_mpc_state"].eq("not_run")
    nonempty = frame["unknown_mpc_state"].isin(
        {"reply_received", "selected_without_reply", "catalog_nonempty"}
    )
    if not frame.loc[zero_catalog | not_run, "unknown_n"].eq(0).all():
        raise ValueError(f"{path}: true-zero/not-run nights must have unknown_n == 0")
    if not frame.loc[nonempty, "unknown_n"].gt(0).all():
        raise ValueError(f"{path}: catalog-nonempty/reply nights must have unknown_n > 0")

    return frame.sort_values("night").reset_index(drop=True)


def attach_per_night(
    stage: pd.DataFrame,
    values: pd.Series,
    column: str,
    source_name: str,
) -> None:
    values = values.copy()
    values.index = normalize_night(
        pd.Series(values.index, dtype="string"), source_name
    ).to_numpy()
    outside = sorted(set(values.index) - set(stage["night"]))
    if outside:
        raise ValueError(f"{source_name} contains nights absent from stage table: {outside[:8]}")
    mapping = pd.to_numeric(values, errors="raise").astype("int64").to_dict()
    stage[column] = stage["night"].map(mapping).fillna(0).astype("int64")


def unique_token_counts(
    frame: pd.DataFrame,
    token_column: str,
    mask: pd.Series | None = None,
) -> pd.Series:
    if mask is not None:
        frame = frame.loc[mask]
    counts: dict[str, int] = {}
    for night, group in frame.groupby("night"):
        tokens: set[str] = set()
        for value in group[token_column].dropna().astype(str):
            tokens.update(token.strip() for token in value.split(";") if token.strip())
        counts[str(night)] = len(tokens)
    return pd.Series(counts, dtype="int64")


def read_false_positive_taxonomy(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    require_file(path, "unknown false-positive taxonomy")
    frame = pd.read_csv(path, low_memory=False)
    required = {"stage", "taxonomy", "count", "unit", "evidence_status"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{path} lacks required columns: {', '.join(missing)}")
    counts = pd.to_numeric(frame["count"], errors="coerce")
    invalid = counts.isna() | (counts < 0) | ~np.isclose(counts, np.round(counts))
    if invalid.any():
        raise ValueError(f"{path}: count must contain non-negative integers")
    frame["count"] = np.round(counts).astype("int64")

    withdrawal = frame.loc[
        frame["taxonomy"].eq("review_real_reverted_before_submission")
    ]
    if len(withdrawal) != 1:
        raise ValueError(
            f"{path}: expected exactly one review_real_reverted_before_submission taxonomy row"
        )
    evidence_status = str(withdrawal.iloc[0]["evidence_status"])
    match = re.fullmatch(r"row_identified_(\d{8})_(.+)", evidence_status)
    identifier = (
        f"{match.group(1)} / {match.group(2)}" if match else evidence_status.replace("_", " ")
    )
    metadata = {
        "review_to_submission_withdrawal_n": int(withdrawal.iloc[0]["count"]),
        "review_to_submission_identifier": identifier,
    }
    return frame, metadata


def read_and_reconcile(
    stage_path: Path,
    taxonomy_path: Path,
    orbit_path: Path,
    unknown_path: Path,
    memberships_path: Path,
    unique_path: Path,
    review_path: Path,
) -> pd.DataFrame:
    stage = read_stage_counts(stage_path)
    taxonomy, metadata = read_false_positive_taxonomy(taxonomy_path)
    for path, purpose in [
        (orbit_path, "orbit-link table"),
        (unknown_path, "post-known link table"),
        (memberships_path, "link-membership table"),
        (unique_path, "unique-detection table"),
        (review_path, "review-status table"),
    ]:
        require_file(path, purpose)

    orbit = pd.read_parquet(orbit_path, columns=["night", "fit_ok", "is_good"])
    orbit["night"] = normalize_night(orbit["night"], "orbit-link night")
    orbit["fit_ok"] = parse_bool(orbit["fit_ok"], "orbit fit_ok")
    orbit["is_good"] = parse_bool(orbit["is_good"], "orbit is_good")
    if (orbit["is_good"] & ~orbit["fit_ok"]).any():
        raise ValueError(f"{orbit_path}: is_good contains a row that is not fit_ok")

    unknown_required = {
        "night",
        "n_obs",
        "tracklet_ids",
        "initial_human_selected",
        "posthoc_retained",
        "fit_ok",
        "is_good",
    }
    unknown = pd.read_csv(unknown_path, dtype={"night": "string"}, low_memory=False)
    missing = sorted(unknown_required - set(unknown.columns))
    if missing:
        raise ValueError(f"{unknown_path} lacks required columns: {', '.join(missing)}")
    unknown["night"] = normalize_night(unknown["night"], "post-known link night")
    for column in ["initial_human_selected", "posthoc_retained", "fit_ok", "is_good"]:
        unknown[column] = parse_bool(unknown[column], f"unknown links {column}")
    if not (unknown["fit_ok"] & unknown["is_good"]).all():
        raise ValueError("post-known frozen catalog contains a row that is not fit_ok and is_good")

    memberships = pd.read_csv(
        memberships_path,
        usecols=lambda column: column in {"night", "detection_key", "final_paper_status"},
        dtype={"night": "string"},
        low_memory=False,
    )
    if not {"night", "detection_key"}.issubset(memberships.columns):
        raise ValueError(f"{memberships_path} lacks night/detection_key")
    memberships["night"] = normalize_night(memberships["night"], "membership night")

    unique = pd.read_csv(
        unique_path,
        usecols=lambda column: column in {"night", "detection_key", "final_paper_status"},
        dtype={"night": "string"},
        low_memory=False,
    )
    if not {"night", "detection_key", "final_paper_status"}.issubset(unique.columns):
        raise ValueError(f"{unique_path} lacks night/detection_key/final_paper_status")
    unique["night"] = normalize_night(unique["night"], "unique-detection night")
    if unique["detection_key"].duplicated().any():
        raise ValueError(f"{unique_path}: detection_key is not globally unique")

    review = pd.read_csv(review_path, dtype={"origin_night": "string"}, low_memory=False)
    review_required = {"origin_night", "final_paper_status"}
    if not review_required.issubset(review.columns):
        raise ValueError(f"{review_path} lacks origin_night/final_paper_status")
    review["origin_night"] = normalize_night(review["origin_night"], "review origin_night")
    allowed_status = {"retained_after_posthoc_audit", "rejected_posthoc"}
    unexpected = sorted(set(review["final_paper_status"].dropna().astype(str)) - allowed_status)
    if unexpected:
        raise ValueError(f"Unexpected final_paper_status values in review table: {unexpected}")

    assert_equal("all shared-endpoint links", stage["link_n"].sum(), len(orbit))
    assert_equal("numerical fit_ok links", stage["orbit_fit_n"].sum(), orbit["fit_ok"].sum())
    assert_equal(
        "thresholded is_good links", stage["orbit_is_good_n"].sum(), orbit["is_good"].sum()
    )
    assert_equal("post-known catalog links", stage["unknown_n"].sum(), len(unknown))

    n_obs = pd.to_numeric(unknown["n_obs"], errors="coerce")
    if n_obs.isna().any() or (n_obs < 0).any():
        raise ValueError(f"{unknown_path}: invalid n_obs")
    assert_equal("link-detection memberships", int(n_obs.sum()), len(memberships))
    assert_equal("unique post-known detections", memberships["detection_key"].nunique(), len(unique))

    initial_mask = unknown["initial_human_selected"]
    retained_mask = unknown["posthoc_retained"]
    review_retained = review["final_paper_status"].eq("retained_after_posthoc_audit")
    assert_equal("submission set (link/review table)", initial_mask.sum(), len(review))
    assert_equal("post-audit retained (link/review table)", retained_mask.sum(), review_retained.sum())
    assert_equal("submission selected (stage/review table)", stage["submit_real_n"].sum(), len(review))
    assert_equal("frozen audit initial set", stage["audit_initial_n"].sum(), len(review))
    assert_equal("post-audit retained (stage/review table)", stage["audit_real_n"].sum(), review_retained.sum())

    withdrawal_n = int(metadata["review_to_submission_withdrawal_n"])
    assert_equal(
        "review-real minus submission-selected withdrawal",
        stage["review_real_n"].sum() - stage["submit_real_n"].sum(),
        withdrawal_n,
    )
    taxonomy_counts = taxonomy.set_index("taxonomy")["count"]
    if "not_selected_by_initial_human_review" in taxonomy_counts:
        assert_equal(
            "not selected by initial review taxonomy",
            stage["unknown_n"].sum() - stage["review_real_n"].sum(),
            taxonomy_counts["not_selected_by_initial_human_review"],
        )
    if "manual_artifact_subtype_unrecorded" in taxonomy_counts:
        assert_equal(
            "post-hoc rejected taxonomy",
            stage["audit_initial_n"].sum() - stage["audit_real_n"].sum(),
            taxonomy_counts["manual_artifact_subtype_unrecorded"],
        )
    if "retained_high_confidence" in taxonomy_counts:
        assert_equal(
            "retained taxonomy", stage["audit_real_n"].sum(), taxonomy_counts["retained_high_confidence"]
        )

    attach_per_night(stage, unique.groupby("night").size(), "postknown_unique_detection_n", "unique detections")
    retained_unique = unique["final_paper_status"].eq("retained_after_posthoc_audit")
    attach_per_night(
        stage,
        unique.loc[retained_unique].groupby("night").size(),
        "retained_unique_detection_n",
        "retained unique detections",
    )
    attach_per_night(
        stage,
        memberships.groupby("night").size(),
        "postknown_membership_n",
        "link memberships",
    )
    attach_per_night(
        stage,
        unique_token_counts(unknown, "tracklet_ids"),
        "postknown_unique_tracklet_n",
        "post-known tracklets",
    )
    attach_per_night(
        stage,
        unique_token_counts(unknown, "tracklet_ids", retained_mask),
        "retained_unique_tracklet_n",
        "retained tracklets",
    )
    attach_per_night(
        stage,
        review.groupby("origin_night").size(),
        "submission_derived_n",
        "review submission set",
    )
    attach_per_night(
        stage,
        review.loc[review_retained].groupby("origin_night").size(),
        "postaudit_derived_n",
        "post-audit review",
    )
    if not stage["submit_real_n"].eq(stage["submission_derived_n"]).all():
        raise ValueError("Per-night submission-selected counts disagree with review_and_mpc_status.csv")
    if not stage["audit_initial_n"].eq(stage["submission_derived_n"]).all():
        raise ValueError("Per-night audit-initial counts disagree with review_and_mpc_status.csv")
    if not stage["audit_real_n"].eq(stage["postaudit_derived_n"]).all():
        raise ValueError("Per-night retained counts disagree with review_and_mpc_status.csv")

    stage.attrs.update(metadata)
    return stage


def distribution(values: pd.Series) -> dict[str, float]:
    array = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    array = array[np.isfinite(array)]
    if not len(array):
        return {"median": np.nan, "p16": np.nan, "p84": np.nan, "p95": np.nan}
    p16, median, p84, p95 = np.percentile(array, [16, 50, 84, 95])
    return {"median": median, "p16": p16, "p84": p84, "p95": p95}


def stage_availability(frame: pd.DataFrame, column: str) -> pd.Series:
    if column in {"l2_detection_n", "mag_flag_prefilter_detection_n", "gaia_survivor_n"}:
        return frame["l2_detection_n"].gt(0)
    if column in {
        "grouped_gaia_input_detection_n",
        "common_area_survivor_detection_n",
        "edge_shell_survivor_detection_n",
        "static_survivor_detection_n",
        "tracklet_n",
        "link_n",
        "orbit_fit_n",
        "orbit_is_good_n",
        "orbit_fit_non_known_n",
        "orbit_is_good_non_known_n",
    }:
        return frame["tracklet_stage_state"].isin({"ok", "zero_groups"})
    if column in {"unknown_n", "review_real_n", "submit_real_n", "audit_initial_n", "audit_real_n"}:
        return ~frame["unknown_mpc_state"].eq("not_run")
    return pd.Series(True, index=frame.index, dtype=bool)


def aggregate_row(
    frame: pd.DataFrame,
    *,
    panel: str,
    grain: str,
    order: int,
    stage_key: str,
    label: str,
    column: str,
    semantics: str,
    plot_role: str,
    source: str,
    note: str = "",
) -> dict[str, object]:
    counts = {
        group: int(frame.loc[frame["quality_group"].eq(group), column].sum())
        for group in GROUP_ORDER
    }
    all_count = int(frame[column].sum())
    available = stage_availability(frame, column)
    analysis_mask = frame["quality_group"].eq("included") & available
    stats = distribution(frame.loc[analysis_mask, column])
    return {
        "panel": panel,
        "grain": grain,
        "order": order,
        "stage_key": stage_key,
        "stage_label": label,
        "unit": grain,
        "stage_semantics": semantics,
        "plot_role": plot_role,
        "source": source,
        "source_column": column,
        "note": note,
        "all_count": all_count,
        "included_count": counts["included"],
        "primary_excluded_count": counts["primary_excluded"],
        "unknown_quarantine_count": counts["unknown_quarantine"],
        "closure_delta": all_count - sum(counts.values()),
        "source_evidence_nights": int(available.sum()),
        "source_evidence_zero_nights": int((available & frame[column].eq(0)).sum()),
        "all_nonzero_nights": int(frame[column].gt(0).sum()),
        "included_nonzero_nights": int(
            (frame["quality_group"].eq("included") & frame[column].gt(0)).sum()
        ),
        "distribution_nights": int(analysis_mask.sum()),
        "night_median": stats["median"],
        "night_p16": stats["p16"],
        "night_p84": stats["p84"],
        "night_p95": stats["p95"],
    }


def status_row(
    frame: pd.DataFrame,
    *,
    order: int,
    key: str,
    label: str,
    mask: pd.Series,
    note: str,
) -> dict[str, object]:
    counts = {
        group: int((mask & frame["quality_group"].eq(group)).sum())
        for group in GROUP_ORDER
    }
    count = int(mask.sum())
    return {
        "panel": "night_status",
        "grain": "night",
        "order": order,
        "stage_key": key,
        "stage_label": label,
        "unit": "night",
        "stage_semantics": "status_category",
        "plot_role": "night_status_bar",
        "source": "unknown_stage_counts_by_night.csv",
        "source_column": "night_group",
        "note": note,
        "all_count": count,
        "included_count": counts["included"],
        "primary_excluded_count": counts["primary_excluded"],
        "unknown_quarantine_count": counts["unknown_quarantine"],
        "closure_delta": count - sum(counts.values()),
        "source_evidence_nights": len(frame),
        "source_evidence_zero_nights": np.nan,
        "all_nonzero_nights": count,
        "included_nonzero_nights": counts["included"],
        "distribution_nights": np.nan,
        "night_median": np.nan,
        "night_p16": np.nan,
        "night_p84": np.nan,
        "night_p95": np.nan,
    }


def build_figure_data(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    detection_specs = [
        ("l2_detection", "L2 source detections", "l2_detection_n", "production_path"),
        (
            "mag_flag_prefilter",
            "Magnitude/flag survivors",
            "mag_flag_prefilter_detection_n",
            "production_path",
        ),
        ("gaia_survivor", "Gaia-mask survivors", "gaia_survivor_n", "production_path"),
        (
            "grouped_gaia_input",
            "Grouped Gaia-filtered input",
            "grouped_gaia_input_detection_n",
            "production_path",
        ),
        (
            "common_area_skipped",
            "Common-area survivors (skipped/no-op)",
            "common_area_survivor_detection_n",
            "no_op_as_executed",
        ),
        (
            "edge_shell_survivor",
            "Edge-shell survivors",
            "edge_shell_survivor_detection_n",
            "production_path",
        ),
        (
            "static_survivor",
            "Reference-anchored static survivors",
            "static_survivor_detection_n",
            "production_path",
        ),
    ]
    for order, (key, label, column, semantics) in enumerate(detection_specs):
        rows.append(
            aggregate_row(
                frame,
                panel="source_detection",
                grain="detection",
                order=order,
                stage_key=key,
                label=label,
                column=column,
                semantics=semantics,
                plot_role="detection_bar",
                source="unknown_stage_counts_by_night.csv",
                note=(
                    "Production used --skip-common-area; grouped input and common-area "
                    "survivor counts are required to be identical."
                    if semantics == "no_op_as_executed"
                    else ""
                ),
            )
        )

    rows.append(
        aggregate_row(
            frame,
            panel="tracklet",
            grain="tracklet",
            order=0,
            stage_key="two_point_tracklet",
            label="Two-point tracklets",
            column="tracklet_n",
            semantics="production_path",
            plot_role="tracklet_bar",
            source="unknown_stage_counts_by_night.csv",
            note="Unit changes from source detections to tracklets at two-point pairing.",
        )
    )

    link_specs = [
        ("shared_link", "Shared-endpoint links", "link_n", "production_path"),
        ("fit_ok", "Numerical fit_ok", "orbit_fit_n", "production_path"),
        (
            "is_good",
            "Thresholded is_good (diagnostic)",
            "orbit_is_good_n",
            "diagnostic_parallel",
        ),
        (
            "fit_ok_non_known",
            "fit_ok + all members non-known",
            "orbit_fit_non_known_n",
            "production_path",
        ),
        (
            "is_good_non_known",
            "is_good + all members non-known (diagnostic)",
            "orbit_is_good_non_known_n",
            "diagnostic_parallel",
        ),
        (
            "post_known_unknown",
            "All-non-known / post-known catalog",
            "unknown_n",
            "production_path",
        ),
    ]
    for order, (key, label, column, semantics) in enumerate(link_specs):
        rows.append(
            aggregate_row(
                frame,
                panel="linkage",
                grain="linkage",
                order=order,
                stage_key=key,
                label=label,
                column=column,
                semantics=semantics,
                plot_role="linkage_bar",
                source="unknown_stage_counts_by_night.csv",
                note=(
                    "Diagnostic threshold; not the selector used to define the production fit_ok path."
                    if semantics == "diagnostic_parallel"
                    else ""
                ),
            )
        )

    withdrawal_identifier = str(frame.attrs["review_to_submission_identifier"])
    withdrawal_n = int(frame.attrs["review_to_submission_withdrawal_n"])
    review_specs = [
        (
            "review_marked_real",
            "Human review marked real",
            "review_real_n",
            "human_review",
            "Unknown-review decisions before submission reconciliation.",
        ),
        (
            "submission_selected",
            "Submission / frozen initial set",
            "submit_real_n",
            "submission_selection",
            f"{withdrawal_n:,} review-real row withdrawn before submission: {withdrawal_identifier}.",
        ),
        (
            "posthoc_retained",
            "Post-hoc retained",
            "audit_real_n",
            "posthoc_audit",
            "Retained high-confidence linkages; not confirmed independent objects.",
        ),
    ]
    for order, (key, label, column, semantics, note) in enumerate(review_specs):
        rows.append(
            aggregate_row(
                frame,
                panel="review_zoom",
                grain="linkage",
                order=order,
                stage_key=key,
                label=label,
                column=column,
                semantics=semantics,
                plot_role="review_bar",
                source="unknown_stage_counts_by_night.csv",
                note=note,
            )
        )
    rows.append(
        aggregate_row(
            frame,
            panel="review_zoom",
            grain="linkage",
            order=20,
            stage_key="audit_initial_crosscheck",
            label="Frozen audit initial set (cross-check)",
            column="audit_initial_n",
            semantics="crosscheck",
            plot_role="context_only",
            source="unknown_stage_counts_by_night.csv",
            note="Required to equal submission_real_n per night and in total.",
        )
    )

    context_specs = [
        (
            "postknown_membership",
            "Post-known link-detection memberships",
            "postknown_membership_n",
            "link-detection membership",
            "unknown_all_link_memberships.csv",
        ),
        (
            "postknown_unique_detection",
            "Unique detections in post-known links",
            "postknown_unique_detection_n",
            "detection",
            "unknown_unique_detections.csv",
        ),
        (
            "retained_unique_detection",
            "Unique detections in retained links",
            "retained_unique_detection_n",
            "detection",
            "unknown_unique_detections.csv",
        ),
        (
            "postknown_unique_tracklet",
            "Unique tracklet IDs in post-known links",
            "postknown_unique_tracklet_n",
            "tracklet",
            "unknown_all_links.csv: tracklet_ids",
        ),
        (
            "retained_unique_tracklet",
            "Unique tracklet IDs in retained links",
            "retained_unique_tracklet_n",
            "tracklet",
            "unknown_all_links.csv: tracklet_ids",
        ),
    ]
    for order, (key, label, column, grain, source) in enumerate(context_specs):
        rows.append(
            aggregate_row(
                frame,
                panel="unit_context",
                grain=grain,
                order=order,
                stage_key=key,
                label=label,
                column=column,
                semantics="context",
                plot_role="context_only",
                source=source,
            )
        )

    quarantine_nights = ", ".join(
        frame.loc[frame["night_group"].eq("quarantine"), "night"].astype(str)
    )
    quarantine_note = (
        "Processed contributions retained for accounting but excluded from unknown science. "
        f"Quarantine nights: {quarantine_nights}."
        if quarantine_nights
        else "No nights are classified as unknown-branch quarantine."
    )
    status_specs = [
        (
            "normal_catalog_nonempty",
            "Catalog non-empty (normal)",
            frame["night_group"].eq("normal"),
            "Unknown branch executed and final catalog was non-empty.",
        ),
        (
            "true_zero_catalog",
            "True zero (executed)",
            frame["night_group"].eq("true_zero_catalog"),
            "Unknown branch executed; final unknown catalog contained zero rows.",
        ),
        (
            "not_run",
            "Not run (quality-included)",
            frame["night_group"].eq("not_run"),
            "Unknown branch was not run; these are not zero-yield executions.",
        ),
        (
            "quarantine",
            "Unknown-branch quarantine",
            frame["night_group"].eq("quarantine"),
            quarantine_note,
        ),
        (
            "primary_excluded",
            "Primary-science excluded",
            frame["night_group"].eq("excluded"),
            "Night excluded by the primary-science quality mask.",
        ),
    ]
    for order, (key, label, mask, note) in enumerate(status_specs):
        rows.append(
            status_row(frame, order=order, key=key, label=label, mask=mask, note=note)
        )

    branch_specs = [
        (
            "branch_catalog_nonempty",
            frame["unknown_mpc_state"].isin(
                {"reply_received", "selected_without_reply", "catalog_nonempty"}
            ),
            "Raw unknown-branch state: catalog non-empty.",
        ),
        (
            "branch_true_zero",
            frame["unknown_mpc_state"].eq("true_zero_catalog"),
            "Raw unknown-branch state: executed true zero.",
        ),
        (
            "branch_not_run",
            frame["unknown_mpc_state"].eq("not_run"),
            "Raw unknown-branch state: not run; includes quarantined/excluded nights.",
        ),
    ]
    for order, (key, mask, note) in enumerate(branch_specs, start=20):
        rows.append(
            status_row(
                frame,
                order=order,
                key=key,
                label=note.split(": ", 1)[-1].rstrip("."),
                mask=mask,
                note=note,
            )
        )
        rows[-1]["plot_role"] = "context_only"

    result = pd.DataFrame(rows).sort_values(["panel", "order", "stage_key"]).reset_index(drop=True)
    finite_closure = result["closure_delta"].dropna()
    if not finite_closure.eq(0).all():
        raise AssertionError("Fig. 9 quality-group accounting failed to close")
    return result


def compact_number(value: float) -> str:
    return f"{int(round(value)):,}" if np.isfinite(value) else "NA"


def count_axis(value: float, _: int) -> str:
    absolute = abs(value)
    if absolute >= 1_000_000_000:
        return f"{value / 1_000_000_000:.1f}B"
    if absolute >= 1_000_000:
        scaled = value / 1_000_000
        return f"{scaled:.1f}M" if absolute < 10_000_000 else f"{scaled:.0f}M"
    if absolute >= 1_000:
        return f"{value / 1_000:.0f}k"
    if np.isclose(value, round(value)):
        return f"{value:,.0f}"
    return f"{value:.2g}"


def included_bar_style(semantics: str) -> dict[str, object]:
    if semantics == "diagnostic_parallel":
        return {
            "color": "white",
            "edgecolor": STATISTICS_COLORS["orange"],
            "hatch": "..",
            "alpha": 1.0,
            "linewidth": 1.25,
        }
    if semantics == "no_op_as_executed":
        return {
            "color": "#dce9f5",
            "edgecolor": STATISTICS_COLORS["blue"],
            "hatch": "///",
            "alpha": 1.0,
            "linewidth": 0.9,
        }
    return {
        "color": STATISTICS_COLORS["blue"],
        "edgecolor": STATISTICS_COLORS["ink"],
        "hatch": "",
        "alpha": 0.78,
        "linewidth": 0.55,
    }


def draw_stage_bars(
    axis,
    rows: pd.DataFrame,
    *,
    title: str = "",
    xlabel: str,
    reserve_rows: float = 0.45,
    tick_fontsize: float = 11.5,
) -> None:
    rows = rows.sort_values("order").reset_index(drop=True)
    maximum = max(float(rows["all_count"].max()), 1.0)
    for index, row in rows.iterrows():
        included = float(row["included_count"])
        primary_excluded = float(row["primary_excluded_count"])
        quarantine = float(row["unknown_quarantine_count"])
        style = included_bar_style(str(row["stage_semantics"]))
        axis.barh(index, included, height=0.64, zorder=3, **style)
        left = included
        for group, value in [
            ("primary_excluded", primary_excluded),
            ("unknown_quarantine", quarantine),
        ]:
            if value <= 0:
                continue
            group_style = GROUP_STYLES[group]
            axis.barh(
                index,
                value,
                left=left,
                height=0.64,
                color=group_style["color"],
                edgecolor=STATISTICS_COLORS["ink"],
                linewidth=0.55,
                hatch=group_style["hatch"],
                zorder=3,
            )
            left += value
        axis.text(
            float(row["all_count"]) + maximum * 0.017,
            index,
            compact_number(float(row["all_count"])),
            ha="left",
            va="center",
            fontsize=11.5,
            fontweight="bold",
        )
    axis.set_yticks(np.arange(len(rows)))
    axis.set_yticklabels(rows["stage_label"], fontsize=tick_fontsize)
    axis.set_ylim(len(rows) - 0.35 + reserve_rows, -0.65)
    axis.set_xlim(0, maximum * 1.27)
    axis.set_xlabel(xlabel, fontsize=13.5)
    axis.xaxis.set_major_formatter(FuncFormatter(count_axis))
    if title:
        axis.set_title(title, loc="left", fontsize=18.5, fontweight="bold", pad=10)
    style_statistics_axis(axis, tick_fontsize=11.5)


def text_box(axis, text: str, *, x: float, y: float, ha: str = "center", fontsize: float = 11.5) -> None:
    axis.text(
        x,
        y,
        text,
        transform=axis.transAxes,
        ha=ha,
        va="center",
        fontsize=fontsize,
        bbox={
            "boxstyle": "round,pad=0.42",
            "facecolor": "white",
            "edgecolor": STATISTICS_COLORS["mid_grey"],
            "linewidth": 0.8,
            "alpha": 0.96,
        },
        zorder=10,
    )


def draw_detection_panel(figure, spec, figure_data: pd.DataFrame) -> None:
    subgrid = spec.subgridspec(2, 1, height_ratios=[0.76, 1.55], hspace=0.40)
    high_axis = figure.add_subplot(subgrid[0, 0])
    lower_axis = figure.add_subplot(subgrid[1, 0])
    rows = figure_data.loc[figure_data["plot_role"].eq("detection_bar")].sort_values("order")
    draw_stage_bars(
        high_axis,
        rows.iloc[:2],
        title="(a) Source-detection funnel",
        xlabel="Source detections — high-volume scale",
        reserve_rows=0.10,
        tick_fontsize=11.5,
    )
    draw_stage_bars(
        lower_axis,
        rows.iloc[2:],
        xlabel="Source detections — filtered-stage scale",
        reserve_rows=0.10,
        tick_fontsize=10.6,
    )
    lower_axis.text(
        0.99,
        0.02,
        "Common-area filtering was skipped as executed: grouped input = survivors.",
        transform=lower_axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.8,
        color=STATISTICS_COLORS["mid_grey"],
    )


def draw_tracklet_panel(axis, figure_data: pd.DataFrame) -> None:
    rows = figure_data.loc[figure_data["plot_role"].eq("tracklet_bar")]
    draw_stage_bars(
        axis,
        rows,
        title="(b) Tracklet grain",
        xlabel="Two-point tracklets",
        reserve_rows=1.90,
        tick_fontsize=12.2,
    )
    lookup = figure_data.set_index("stage_key")
    static_count = float(lookup.loc["static_survivor", "all_count"])
    tracklet_count = float(lookup.loc["two_point_tracklet", "all_count"])
    link_count = float(lookup.loc["shared_link", "all_count"])
    text_box(
        axis,
        f"{compact_number(static_count)}\nreference-anchored static survivors\n(unit: source detections)",
        x=0.27,
        y=0.25,
        fontsize=10.8,
    )
    axis.annotate(
        "two-point pairing\nchanges the unit",
        xy=(0.50, 0.34),
        xytext=(0.50, 0.34),
        xycoords="axes fraction",
        textcoords="axes fraction",
        ha="center",
        va="center",
        fontsize=10.2,
        color=STATISTICS_COLORS["mid_grey"],
    )
    axis.annotate(
        "",
        xy=(0.59, 0.25),
        xytext=(0.41, 0.25),
        xycoords="axes fraction",
        arrowprops={"arrowstyle": "->", "color": STATISTICS_COLORS["mid_grey"], "lw": 1.2},
    )
    text_box(
        axis,
        f"{compact_number(tracklet_count)}\ntwo-point tracklets\n→ {compact_number(link_count)} shared links",
        x=0.74,
        y=0.25,
        fontsize=10.8,
    )
    axis.text(
        0.5,
        0.055,
        "No continuous funnel line crosses detection, tracklet, and linkage grains.",
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=10.5,
        color=STATISTICS_COLORS["mid_grey"],
        style="italic",
    )


def draw_linkage_panel(figure, spec, figure_data: pd.DataFrame) -> None:
    subgrid = spec.subgridspec(2, 1, height_ratios=[2.15, 1.08], hspace=0.58)
    link_axis = figure.add_subplot(subgrid[0, 0])
    review_axis = figure.add_subplot(subgrid[1, 0])
    link_rows = figure_data.loc[figure_data["plot_role"].eq("linkage_bar")]
    review_rows = figure_data.loc[figure_data["plot_role"].eq("review_bar")]
    draw_stage_bars(
        link_axis,
        link_rows,
        title="(c) Linkage and orbit gates",
        xlabel="Single-night linkages",
        reserve_rows=0.15,
        tick_fontsize=9.8,
    )
    draw_stage_bars(
        review_axis,
        review_rows,
        title="Human-review zoom (separate x-axis)",
        xlabel="Reviewed linkages",
        reserve_rows=0.18,
        tick_fontsize=10.1,
    )
    lookup = review_rows.set_index("stage_key")
    review_n = int(lookup.loc["review_marked_real", "all_count"])
    submit_n = int(lookup.loc["submission_selected", "all_count"])
    retained_n = int(lookup.loc["posthoc_retained", "all_count"])
    submission_note = str(lookup.loc["submission_selected", "note"])
    review_axis.text(
        0.99,
        -0.42,
        f"{review_n:,} → {submit_n:,}: {submission_note}  "
        f"{submit_n:,} → {retained_n:,}: {submit_n - retained_n:,} rejected in post-hoc audit.",
        transform=review_axis.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
        color=STATISTICS_COLORS["mid_grey"],
    )


def draw_night_status_panel(axis, figure_data: pd.DataFrame) -> None:
    rows = figure_data.loc[figure_data["plot_role"].eq("night_status_bar")].sort_values("order")
    style_map = {
        "normal_catalog_nonempty": (STATISTICS_COLORS["blue"], ""),
        "true_zero_catalog": (STATISTICS_COLORS["gold"], ""),
        "not_run": (STATISTICS_COLORS["light_grey"], ".."),
        "quarantine": (STATISTICS_COLORS["orange"], "////"),
        "primary_excluded": ("#8d8d8d", "xx"),
    }
    maximum = max(float(rows["all_count"].max()), 1.0)
    for index, row in rows.reset_index(drop=True).iterrows():
        color, hatch = style_map[str(row["stage_key"])]
        axis.barh(
            index,
            float(row["all_count"]),
            height=0.64,
            color=color,
            alpha=0.80,
            hatch=hatch,
            edgecolor=STATISTICS_COLORS["ink"],
            linewidth=0.65,
            zorder=3,
        )
        axis.text(
            float(row["all_count"]) + maximum * 0.025,
            index,
            compact_number(float(row["all_count"])),
            ha="left",
            va="center",
            fontsize=12,
            fontweight="bold",
        )
    axis.set_yticks(np.arange(len(rows)))
    axis.set_yticklabels(rows["stage_label"], fontsize=11.4)
    axis.set_ylim(len(rows) + 1.55, -0.65)
    axis.set_xlim(0, maximum * 1.30)
    axis.set_xlabel("Frozen calendar nights", fontsize=13.5)
    axis.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    axis.set_title("(d) True zero, not run, and quality disposition", loc="left", fontsize=18.5, fontweight="bold")
    style_statistics_axis(axis, tick_fontsize=11.5)

    context = figure_data.set_index("stage_key")
    raw_nonempty = int(context.loc["branch_catalog_nonempty", "all_count"])
    raw_zero = int(context.loc["branch_true_zero", "all_count"])
    raw_not_run = int(context.loc["branch_not_run", "all_count"])
    quarantine_rows = rows.loc[rows["stage_key"].eq("quarantine")]
    quarantine_note = str(quarantine_rows.iloc[0]["note"]).replace(
        ". Quarantine nights:", ".\nQuarantine nights:"
    )
    axis.text(
        0.01,
        0.18,
        f"Raw branch state: {raw_nonempty:,} catalog non-empty; {raw_zero:,} executed true zero; "
        f"{raw_not_run:,} not run.\n"
        "Night-group bars are exhaustive: quarantine/excluded nights are separated from quality-included not-run nights.\n"
        f"{quarantine_note}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=10.0,
        color=STATISTICS_COLORS["mid_grey"],
        linespacing=1.25,
    )


def make_figure(figure_data: pd.DataFrame):
    apply_statistics_style()
    figure = plt.figure(figsize=(18.0, 16.2))
    outer = figure.add_gridspec(
        2,
        2,
        height_ratios=[1.0, 1.12],
        left=0.14,
        right=0.98,
        top=0.895,
        bottom=0.095,
        hspace=0.34,
        wspace=0.38,
    )
    draw_detection_panel(figure, outer[0, 0], figure_data)
    tracklet_axis = figure.add_subplot(outer[0, 1])
    draw_tracklet_panel(tracklet_axis, figure_data)
    draw_linkage_panel(figure, outer[1, 0], figure_data)
    night_axis = figure.add_subplot(outer[1, 1])
    draw_night_status_panel(night_axis, figure_data)

    figure.suptitle(
        "Unknown-candidate accounting by statistical grain",
        fontsize=25,
        fontweight="bold",
        y=0.975,
    )
    figure.text(
        0.5,
        0.934,
        "Full as-executed source funnel → two-point tracklets → linkage/orbit gates → review; exact counts include all frozen nights.",
        ha="center",
        fontsize=13.0,
        color=STATISTICS_COLORS["mid_grey"],
    )
    legend_handles = [
        Patch(
            facecolor=STATISTICS_COLORS["blue"],
            alpha=0.78,
            edgecolor=STATISTICS_COLORS["ink"],
            label="Unknown-science included contribution",
        ),
        Patch(
            facecolor=GROUP_STYLES["primary_excluded"]["color"],
            edgecolor=STATISTICS_COLORS["ink"],
            hatch=GROUP_STYLES["primary_excluded"]["hatch"],
            label=GROUP_STYLES["primary_excluded"]["label"],
        ),
        Patch(
            facecolor=GROUP_STYLES["unknown_quarantine"]["color"],
            edgecolor=STATISTICS_COLORS["ink"],
            hatch=GROUP_STYLES["unknown_quarantine"]["hatch"],
            label=GROUP_STYLES["unknown_quarantine"]["label"],
        ),
        Patch(
            facecolor="#dce9f5",
            edgecolor=STATISTICS_COLORS["blue"],
            hatch="///",
            label="Common-area stage skipped/no-op",
        ),
        Patch(
            facecolor="white",
            edgecolor=STATISTICS_COLORS["orange"],
            hatch="..",
            label="is_good diagnostic (parallel)",
        ),
    ]
    figure.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.018),
        ncol=5,
        fontsize=10.9,
        frameon=False,
    )
    return figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-counts", type=Path, default=DEFAULT_STAGE_COUNTS)
    parser.add_argument(
        "--false-positive-taxonomy", type=Path, default=DEFAULT_FALSE_POSITIVE_TAXONOMY
    )
    parser.add_argument("--orbit-links", type=Path, default=DEFAULT_ORBIT_LINKS)
    parser.add_argument("--unknown-links", type=Path, default=DEFAULT_UNKNOWN_LINKS)
    parser.add_argument("--memberships", type=Path, default=DEFAULT_MEMBERSHIPS)
    parser.add_argument("--unique-detections", type=Path, default=DEFAULT_UNIQUE_DETECTIONS)
    parser.add_argument("--review-status", type=Path, default=DEFAULT_REVIEW_STATUS)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure-data", type=Path, default=DEFAULT_FIGURE_DATA)
    args = parser.parse_args()

    by_night = read_and_reconcile(
        args.stage_counts,
        args.false_positive_taxonomy,
        args.orbit_links,
        args.unknown_links,
        args.memberships,
        args.unique_detections,
        args.review_status,
    )
    figure_data = build_figure_data(by_night)
    args.figure_data.parent.mkdir(parents=True, exist_ok=True)
    figure_data.to_csv(args.figure_data, index=False)
    by_night_path = args.figure_data.with_name(f"{args.figure_data.stem}_by_night.csv")
    by_night.to_csv(by_night_path, index=False)
    figure = make_figure(figure_data)
    png_path, pdf_path = save_pdf_png(figure, args.output_base)
    print(f"wrote {pdf_path}")
    print(f"wrote {png_path}")
    print(f"wrote {args.figure_data}")
    print(f"wrote {by_night_path}")


if __name__ == "__main__":
    main()
