#!/usr/bin/env python3
"""Screen the 58 retained unknown links for repeated linear-motion candidates.

This is deliberately a *candidate-screening* analysis.  A connected component
from this script is not a confirmed independent Solar-system body, an orbital
identity, or an MPC/JPL identification.  It only records links whose local
linear sky-plane motion is mutually compatible under fixed, documented cuts.

For every pair with ``0 < delta_t <= 7 days``, the earlier link is propagated
forward with its own fitted linear velocity and the later link is propagated
backward with its own fitted linear velocity.  Both residuals are evaluated in
the tangent plane of the propagation origin.  The production convention is
used: ``lin_dir_deg = atan2(north, east)``, so 0 deg is increasing
``RA*cos(Dec)`` and 90 deg is increasing Dec.

The primary graph definition is intentionally fixed in code:

* max(forward residual, backward residual) <= 0.03 deg;
* absolute speed difference <= 2 arcsec/hour; and
* wrapped direction difference <= 5 deg.

Outputs are written atomically and input products are never modified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


SCHEMA_VERSION = "1.0"
EXPECTED_LINK_COUNT = 58
MAX_DELTA_DAYS = 7.0
KEY_COLUMNS = ("night", "trk_sub", "linkage_id")
REQUIRED_COLUMNS = (
    "night",
    "trk_sub",
    "linkage_id",
    "median_mjd",
    "median_ra_deg",
    "median_dec_deg",
    "speed_arcsec_per_hour",
    "lin_dir_deg",
)
CANDIDATE_DISCLAIMER = (
    "Linear-motion repeat candidate only; not a confirmed independent object, "
    "orbital identity, or MPC/JPL identification."
)


@dataclass(frozen=True)
class Threshold:
    name: str
    residual_deg: float
    speed_arcsec_per_hour: float
    direction_deg: float

    @property
    def edge_column(self) -> str:
        return f"is_edge_{self.name}"

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "maximum_two_way_residual_deg": self.residual_deg,
            "maximum_speed_difference_arcsec_per_hour": self.speed_arcsec_per_hour,
            "maximum_direction_difference_deg": self.direction_deg,
            "comparison_operator": "all metrics <= their thresholds",
        }


THRESHOLDS: tuple[Threshold, ...] = (
    Threshold("strict", 0.01, 1.0, 3.0),
    Threshold("conservative_primary", 0.03, 2.0, 5.0),
    Threshold("relaxed", 0.05, 3.0, 7.0),
)
PRIMARY_THRESHOLD = THRESHOLDS[1]


class UnionFind:
    """Small deterministic disjoint-set implementation."""

    def __init__(self, items: Iterable[str]) -> None:
        self.parent = {item: item for item in items}
        self.rank = {item: 0 for item in items}

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        rank_left = self.rank[root_left]
        rank_right = self.rank[root_right]
        if rank_left < rank_right:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if rank_left == rank_right:
            self.rank[root_left] += 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--links",
        type=Path,
        required=True,
        help="The frozen 58-row unknown_high_confidence_links_recomputed.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for pair, assignment, and JSON summary products.",
    )
    parser.add_argument(
        "--review-status",
        type=Path,
        help=(
            "Optional review_and_mpc_status.csv to merge into a new copy. "
            "The source file is never modified."
        ),
    )
    parser.add_argument(
        "--review-status-output",
        type=Path,
        help=(
            "Optional path for the merged review-status copy.  Requires "
            "--review-status."
        ),
    )
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_night(value: Any) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    digits = "".join(character for character in text if character.isdigit())
    if len(digits) < 8:
        raise ValueError(f"Cannot normalize night value {value!r} to YYYYMMDD")
    return digits[:8]


def parse_truth(value: Any) -> bool | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n"}:
        return False
    return None


def link_key(night: Any, trk_sub: Any, linkage_id: Any) -> str:
    return f"{normalize_night(night)}|{str(trk_sub).strip()}|{int(linkage_id)}"


def wrapped_ra_delta_deg(target_ra_deg: float, origin_ra_deg: float) -> float:
    """Shortest signed target-origin RA difference in [-180, 180) degrees."""

    return float((target_ra_deg - origin_ra_deg + 180.0) % 360.0 - 180.0)


def wrapped_direction_difference_deg(left_deg: float, right_deg: float) -> float:
    """Smallest absolute direction separation in [0, 180] degrees."""

    return float(abs((left_deg - right_deg + 180.0) % 360.0 - 180.0))


def velocity_components(
    speed_arcsec_per_hour: float, direction_deg: float
) -> tuple[float, float]:
    angle = math.radians(direction_deg)
    return (
        speed_arcsec_per_hour * math.cos(angle),
        speed_arcsec_per_hour * math.sin(angle),
    )


def propagate_and_residual(
    *,
    origin_ra_deg: float,
    origin_dec_deg: float,
    target_ra_deg: float,
    target_dec_deg: float,
    speed_arcsec_per_hour: float,
    direction_deg: float,
    elapsed_hours: float,
) -> dict[str, float]:
    """Propagate in the origin tangent plane and compare with a target point."""

    cos_dec = math.cos(math.radians(origin_dec_deg))
    if abs(cos_dec) < 1.0e-8:
        raise ValueError(
            f"Tangent-plane RA propagation is ill-conditioned at Dec={origin_dec_deg}"
        )

    velocity_east, velocity_north = velocity_components(
        speed_arcsec_per_hour, direction_deg
    )
    predicted_east_deg = velocity_east * elapsed_hours / 3600.0
    predicted_north_deg = velocity_north * elapsed_hours / 3600.0
    observed_east_deg = (
        wrapped_ra_delta_deg(target_ra_deg, origin_ra_deg) * cos_dec
    )
    observed_north_deg = target_dec_deg - origin_dec_deg
    residual_east_deg = observed_east_deg - predicted_east_deg
    residual_north_deg = observed_north_deg - predicted_north_deg
    residual_deg = math.hypot(residual_east_deg, residual_north_deg)

    predicted_ra_deg = (origin_ra_deg + predicted_east_deg / cos_dec) % 360.0
    predicted_dec_deg = origin_dec_deg + predicted_north_deg
    return {
        "predicted_ra_deg": predicted_ra_deg,
        "predicted_dec_deg": predicted_dec_deg,
        "observed_east_offset_deg": observed_east_deg,
        "observed_north_offset_deg": observed_north_deg,
        "predicted_east_offset_deg": predicted_east_deg,
        "predicted_north_offset_deg": predicted_north_deg,
        "residual_east_deg": residual_east_deg,
        "residual_north_deg": residual_north_deg,
        "residual_deg": residual_deg,
    }


def validate_links(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Links input does not exist: {path}")
    frame = pd.read_csv(path, dtype={"night": "string", "trk_sub": "string"})
    missing = sorted(set(REQUIRED_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"Links input is missing required columns: {missing}")
    if len(frame) != EXPECTED_LINK_COUNT:
        raise AssertionError(
            f"Expected exactly {EXPECTED_LINK_COUNT} retained links, found {len(frame)}"
        )

    result = frame.copy()
    result["source_row_order"] = np.arange(len(result), dtype=int)
    result["night"] = result["night"].map(normalize_night)
    result["trk_sub"] = result["trk_sub"].astype("string").str.strip()
    if result["trk_sub"].isna().any() or (result["trk_sub"] == "").any():
        raise ValueError("trk_sub contains missing or empty values")

    result["linkage_id"] = pd.to_numeric(result["linkage_id"], errors="raise")
    if not np.allclose(result["linkage_id"], np.round(result["linkage_id"])):
        raise ValueError("linkage_id contains non-integer values")
    result["linkage_id"] = result["linkage_id"].astype(int)

    numeric_columns = (
        "median_mjd",
        "median_ra_deg",
        "median_dec_deg",
        "speed_arcsec_per_hour",
        "lin_dir_deg",
    )
    for column in numeric_columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")
        if not np.isfinite(result[column].to_numpy(dtype=float)).all():
            raise ValueError(f"{column} contains missing or non-finite values")

    if not result["median_ra_deg"].between(0.0, 360.0, inclusive="left").all():
        raise ValueError("median_ra_deg must be in [0, 360)")
    if not result["median_dec_deg"].between(-90.0, 90.0).all():
        raise ValueError("median_dec_deg must be in [-90, 90]")
    if (result["speed_arcsec_per_hour"] < 0.0).any():
        raise ValueError("speed_arcsec_per_hour cannot be negative")
    result["lin_dir_deg"] = result["lin_dir_deg"] % 360.0

    if result.duplicated(list(KEY_COLUMNS)).any():
        duplicates = result.loc[
            result.duplicated(list(KEY_COLUMNS), keep=False), list(KEY_COLUMNS)
        ]
        raise ValueError(
            "Duplicate link keys found:\n" + duplicates.to_string(index=False)
        )

    if "posthoc_retained" in result.columns:
        retained = result["posthoc_retained"].map(parse_truth)
        if retained.isna().any() or not retained.all():
            raise ValueError("All 58 input links must have posthoc_retained=True")
    if "final_paper_status" in result.columns:
        status = result["final_paper_status"].astype("string").str.strip()
        if not status.eq("retained_after_posthoc_audit").all():
            raise ValueError(
                "All 58 input links must have final_paper_status="
                "retained_after_posthoc_audit"
            )
    for column in ("fit_ok", "is_good"):
        if column in result.columns:
            parsed = result[column].map(parse_truth)
            if parsed.isna().any() or not parsed.all():
                raise ValueError(f"All retained links must have {column}=True")

    result["link_key"] = [
        link_key(night, trk_sub, linkage_id)
        for night, trk_sub, linkage_id in result[list(KEY_COLUMNS)].itertuples(
            index=False, name=None
        )
    ]
    if result["link_key"].nunique() != EXPECTED_LINK_COUNT:
        raise AssertionError("Normalized link_key is not unique across the 58 rows")
    return result


def edge_passes(
    maximum_residual_deg: float,
    speed_difference_arcsec_per_hour: float,
    direction_difference_deg: float,
    threshold: Threshold,
) -> bool:
    return bool(
        maximum_residual_deg <= threshold.residual_deg
        and speed_difference_arcsec_per_hour <= threshold.speed_arcsec_per_hour
        and direction_difference_deg <= threshold.direction_deg
    )


def primary_failure_reasons(
    maximum_residual_deg: float,
    speed_difference_arcsec_per_hour: float,
    direction_difference_deg: float,
) -> str:
    reasons: list[str] = []
    if maximum_residual_deg > PRIMARY_THRESHOLD.residual_deg:
        reasons.append("two_way_residual")
    if speed_difference_arcsec_per_hour > PRIMARY_THRESHOLD.speed_arcsec_per_hour:
        reasons.append("speed_difference")
    if direction_difference_deg > PRIMARY_THRESHOLD.direction_deg:
        reasons.append("direction_difference")
    return ";".join(reasons)


def build_pair_metrics(links: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    ordered = links.sort_values(
        ["median_mjd", "night", "trk_sub", "linkage_id"], kind="stable"
    ).reset_index(drop=True)
    rows: list[dict[str, Any]] = []
    total_possible = len(ordered) * (len(ordered) - 1) // 2
    zero_delta = 0
    beyond_window = 0

    for left_index in range(len(ordered) - 1):
        earlier = ordered.iloc[left_index]
        for right_index in range(left_index + 1, len(ordered)):
            later = ordered.iloc[right_index]
            delta_days = float(later["median_mjd"] - earlier["median_mjd"])
            if delta_days <= 0.0:
                zero_delta += 1
                continue
            if delta_days > MAX_DELTA_DAYS:
                beyond_window += len(ordered) - right_index
                break

            delta_hours = delta_days * 24.0
            forward = propagate_and_residual(
                origin_ra_deg=float(earlier["median_ra_deg"]),
                origin_dec_deg=float(earlier["median_dec_deg"]),
                target_ra_deg=float(later["median_ra_deg"]),
                target_dec_deg=float(later["median_dec_deg"]),
                speed_arcsec_per_hour=float(earlier["speed_arcsec_per_hour"]),
                direction_deg=float(earlier["lin_dir_deg"]),
                elapsed_hours=delta_hours,
            )
            backward = propagate_and_residual(
                origin_ra_deg=float(later["median_ra_deg"]),
                origin_dec_deg=float(later["median_dec_deg"]),
                target_ra_deg=float(earlier["median_ra_deg"]),
                target_dec_deg=float(earlier["median_dec_deg"]),
                speed_arcsec_per_hour=float(later["speed_arcsec_per_hour"]),
                direction_deg=float(later["lin_dir_deg"]),
                elapsed_hours=-delta_hours,
            )
            maximum_residual = max(
                forward["residual_deg"], backward["residual_deg"]
            )
            speed_difference = abs(
                float(earlier["speed_arcsec_per_hour"])
                - float(later["speed_arcsec_per_hour"])
            )
            direction_difference = wrapped_direction_difference_deg(
                float(earlier["lin_dir_deg"]), float(later["lin_dir_deg"])
            )
            same_night = bool(earlier["night"] == later["night"])

            row: dict[str, Any] = {
                "pair_id": f"PAIR{len(rows) + 1:06d}",
                "pair_class": "same_night" if same_night else "cross_night",
                "same_night": same_night,
                "earlier_link_key": earlier["link_key"],
                "earlier_night": earlier["night"],
                "earlier_trk_sub": earlier["trk_sub"],
                "earlier_linkage_id": int(earlier["linkage_id"]),
                "earlier_median_mjd": float(earlier["median_mjd"]),
                "earlier_median_ra_deg": float(earlier["median_ra_deg"]),
                "earlier_median_dec_deg": float(earlier["median_dec_deg"]),
                "earlier_speed_arcsec_per_hour": float(
                    earlier["speed_arcsec_per_hour"]
                ),
                "earlier_direction_deg": float(earlier["lin_dir_deg"]),
                "later_link_key": later["link_key"],
                "later_night": later["night"],
                "later_trk_sub": later["trk_sub"],
                "later_linkage_id": int(later["linkage_id"]),
                "later_median_mjd": float(later["median_mjd"]),
                "later_median_ra_deg": float(later["median_ra_deg"]),
                "later_median_dec_deg": float(later["median_dec_deg"]),
                "later_speed_arcsec_per_hour": float(
                    later["speed_arcsec_per_hour"]
                ),
                "later_direction_deg": float(later["lin_dir_deg"]),
                "delta_t_days": delta_days,
                "delta_t_hours": delta_hours,
            }
            row.update({f"forward_{key}": value for key, value in forward.items()})
            row.update({f"backward_{key}": value for key, value in backward.items()})
            row.update(
                {
                    "max_two_way_residual_deg": maximum_residual,
                    "speed_difference_arcsec_per_hour": speed_difference,
                    "direction_difference_deg": direction_difference,
                }
            )
            for threshold in THRESHOLDS:
                row[threshold.edge_column] = edge_passes(
                    maximum_residual,
                    speed_difference,
                    direction_difference,
                    threshold,
                )
            row["is_primary_conservative_edge"] = row[
                PRIMARY_THRESHOLD.edge_column
            ]
            row["primary_edge_failure_reasons"] = primary_failure_reasons(
                maximum_residual, speed_difference, direction_difference
            )
            rows.append(row)

    pair_frame = pd.DataFrame(rows)
    accounted = len(pair_frame) + zero_delta + beyond_window
    if accounted != total_possible:
        raise AssertionError(
            f"Pair accounting failed: {accounted} != {total_possible}"
        )
    return pair_frame, {
        "all_unordered_pairs": total_possible,
        "eligible_positive_pairs_within_7_days": len(pair_frame),
        "zero_or_negative_delta_pairs_excluded": zero_delta,
        "pairs_beyond_7_days_excluded": beyond_window,
    }


def graph_components(
    links: pd.DataFrame, pairs: pd.DataFrame, edge_column: str
) -> tuple[list[list[str]], dict[str, int]]:
    keys = links["link_key"].tolist()
    union_find = UnionFind(keys)
    edge_rows = pairs[pairs[edge_column].astype(bool)]
    for row in edge_rows.itertuples(index=False):
        union_find.union(row.earlier_link_key, row.later_link_key)

    groups: dict[str, list[str]] = defaultdict(list)
    for key in keys:
        groups[union_find.find(key)].append(key)

    lookup = links.set_index("link_key")

    def member_sort_key(key: str) -> tuple[float, str, str, int]:
        row = lookup.loc[key]
        return (
            float(row["median_mjd"]),
            str(row["night"]),
            str(row["trk_sub"]),
            int(row["linkage_id"]),
        )

    components = [sorted(members, key=member_sort_key) for members in groups.values()]
    components.sort(key=lambda members: member_sort_key(members[0]))
    component_index = {
        member: index for index, members in enumerate(components) for member in members
    }
    return components, component_index


def threshold_summary(
    links: pd.DataFrame, pairs: pd.DataFrame, threshold: Threshold
) -> dict[str, Any]:
    components, _ = graph_components(links, pairs, threshold.edge_column)
    edge_rows = pairs[pairs[threshold.edge_column].astype(bool)]
    component_sizes = [len(component) for component in components]
    size_counts = Counter(component_sizes)
    return {
        "threshold": threshold.as_dict(),
        "eligible_pair_count": int(len(pairs)),
        "edge_count": int(len(edge_rows)),
        "same_night_edge_count": int(edge_rows["same_night"].sum()),
        "cross_night_edge_count": int((~edge_rows["same_night"]).sum()),
        "component_count": len(components),
        "non_singleton_component_count": sum(size > 1 for size in component_sizes),
        "singleton_component_count": sum(size == 1 for size in component_sizes),
        "links_in_non_singleton_components": sum(
            size for size in component_sizes if size > 1
        ),
        "largest_component_size": max(component_sizes) if component_sizes else 0,
        "component_size_distribution": {
            str(size): count for size, count in sorted(size_counts.items())
        },
    }


def metric_distribution(series: pd.Series) -> dict[str, float | int | None]:
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "n": 0,
            "minimum": None,
            "p50": None,
            "p90": None,
            "p95": None,
            "maximum": None,
        }
    return {
        "n": int(values.size),
        "minimum": float(np.min(values)),
        "p50": float(np.percentile(values, 50.0)),
        "p90": float(np.percentile(values, 90.0)),
        "p95": float(np.percentile(values, 95.0)),
        "maximum": float(np.max(values)),
    }


def pair_class_summary(pairs: pd.DataFrame) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for pair_class in ("same_night", "cross_night", "all"):
        subset = pairs if pair_class == "all" else pairs[pairs["pair_class"] == pair_class]
        primary_edges = subset[subset["is_primary_conservative_edge"].astype(bool)]
        summary[pair_class] = {
            "eligible_pair_count": int(len(subset)),
            "primary_edge_count": int(len(primary_edges)),
            "primary_edge_fraction": (
                float(len(primary_edges) / len(subset)) if len(subset) else None
            ),
            "max_two_way_residual_deg": metric_distribution(
                subset["max_two_way_residual_deg"]
            ),
            "speed_difference_arcsec_per_hour": metric_distribution(
                subset["speed_difference_arcsec_per_hour"]
            ),
            "direction_difference_deg": metric_distribution(
                subset["direction_difference_deg"]
            ),
        }
    return summary


def build_primary_assignments(
    links: pd.DataFrame, pairs: pd.DataFrame
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    components, component_index = graph_components(
        links, pairs, "is_primary_conservative_edge"
    )
    group_ids = {
        index: f"LMG{index + 1:03d}" for index in range(len(components))
    }
    edge_rows = pairs[pairs["is_primary_conservative_edge"].astype(bool)].copy()
    links_by_key = links.set_index("link_key", drop=False)
    group_summaries: list[dict[str, Any]] = []

    for index, members in enumerate(components):
        member_set = set(members)
        group_edges = edge_rows[
            edge_rows["earlier_link_key"].isin(member_set)
            & edge_rows["later_link_key"].isin(member_set)
        ]
        member_rows = links_by_key.loc[members]
        group_summaries.append(
            {
                "linear_motion_candidate_group_id": group_ids[index],
                "link_count": len(members),
                "distinct_night_count": int(member_rows["night"].nunique()),
                "direct_primary_edge_count": int(len(group_edges)),
                "same_night_direct_primary_edge_count": int(
                    group_edges["same_night"].sum()
                ),
                "cross_night_direct_primary_edge_count": int(
                    (~group_edges["same_night"]).sum()
                ),
                "has_cross_night_primary_edge": bool(
                    (~group_edges["same_night"]).any()
                ),
                "is_singleton": len(members) == 1,
                "minimum_mjd": float(member_rows["median_mjd"].min()),
                "maximum_mjd": float(member_rows["median_mjd"].max()),
                "time_span_days": float(
                    member_rows["median_mjd"].max()
                    - member_rows["median_mjd"].min()
                ),
                "member_nights": member_rows["night"].astype(str).tolist(),
                "member_link_keys": members,
                "interpretation": CANDIDATE_DISCLAIMER,
            }
        )

    assignment = links.copy()
    assignment["linear_motion_candidate_group_id"] = assignment["link_key"].map(
        lambda key: group_ids[component_index[key]]
    )
    summary_by_id = {
        group["linear_motion_candidate_group_id"]: group for group in group_summaries
    }
    assignment["candidate_group_link_count"] = assignment[
        "linear_motion_candidate_group_id"
    ].map(lambda value: summary_by_id[value]["link_count"])
    assignment["candidate_group_distinct_night_count"] = assignment[
        "linear_motion_candidate_group_id"
    ].map(lambda value: summary_by_id[value]["distinct_night_count"])
    assignment["candidate_group_direct_primary_edge_count"] = assignment[
        "linear_motion_candidate_group_id"
    ].map(lambda value: summary_by_id[value]["direct_primary_edge_count"])
    assignment["candidate_group_same_night_primary_edge_count"] = assignment[
        "linear_motion_candidate_group_id"
    ].map(lambda value: summary_by_id[value]["same_night_direct_primary_edge_count"])
    assignment["candidate_group_cross_night_primary_edge_count"] = assignment[
        "linear_motion_candidate_group_id"
    ].map(lambda value: summary_by_id[value]["cross_night_direct_primary_edge_count"])
    assignment["candidate_group_has_cross_night_primary_edge"] = assignment[
        "linear_motion_candidate_group_id"
    ].map(lambda value: summary_by_id[value]["has_cross_night_primary_edge"])
    assignment["candidate_group_is_singleton"] = assignment[
        "linear_motion_candidate_group_id"
    ].map(lambda value: summary_by_id[value]["is_singleton"])
    assignment["candidate_interpretation"] = CANDIDATE_DISCLAIMER

    leading = [
        "night",
        "trk_sub",
        "linkage_id",
        "link_key",
        "median_mjd",
        "median_ra_deg",
        "median_dec_deg",
        "speed_arcsec_per_hour",
        "lin_dir_deg",
        "linear_motion_candidate_group_id",
        "candidate_group_link_count",
        "candidate_group_distinct_night_count",
        "candidate_group_direct_primary_edge_count",
        "candidate_group_same_night_primary_edge_count",
        "candidate_group_cross_night_primary_edge_count",
        "candidate_group_has_cross_night_primary_edge",
        "candidate_group_is_singleton",
        "candidate_interpretation",
        "source_row_order",
    ]
    remainder = [column for column in assignment.columns if column not in leading]
    assignment = assignment[leading + remainder].sort_values(
        "source_row_order", kind="stable"
    )
    if len(assignment) != EXPECTED_LINK_COUNT:
        raise AssertionError("Primary assignment output must contain exactly 58 rows")
    return assignment, group_summaries


def atomic_csv(frame: pd.DataFrame, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".inprogress")
    frame.to_csv(temporary, index=False)
    temporary.replace(destination)


def atomic_json(payload: dict[str, Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".inprogress")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(destination)


def merge_review_status(
    review_path: Path, output_path: Path, assignment: pd.DataFrame
) -> dict[str, Any]:
    if not review_path.is_file():
        raise FileNotFoundError(f"Review-status input does not exist: {review_path}")
    if review_path.resolve() == output_path.resolve():
        raise ValueError("Merged review output must not overwrite its source input")

    source_hash_before = sha256_file(review_path)
    review = pd.read_csv(
        review_path, dtype={"origin_night": "string", "trk_sub": "string"}
    )
    required = {"origin_night", "trk_sub", "linkage_id"}
    missing = sorted(required - set(review.columns))
    if missing:
        raise ValueError(f"Review-status input is missing key columns: {missing}")

    review_copy = review.copy()
    review_copy["_review_row_order"] = np.arange(len(review_copy), dtype=int)
    review_copy["origin_night"] = review_copy["origin_night"].map(normalize_night)
    review_copy["trk_sub"] = review_copy["trk_sub"].astype("string").str.strip()
    review_copy["linkage_id"] = pd.to_numeric(
        review_copy["linkage_id"], errors="raise"
    ).astype(int)
    if review_copy.duplicated(["origin_night", "trk_sub", "linkage_id"]).any():
        raise ValueError("Review-status input contains duplicate link keys")

    merge_columns = [
        "linear_motion_candidate_group_id",
        "candidate_group_link_count",
        "candidate_group_distinct_night_count",
        "candidate_group_direct_primary_edge_count",
        "candidate_group_same_night_primary_edge_count",
        "candidate_group_cross_night_primary_edge_count",
        "candidate_group_has_cross_night_primary_edge",
        "candidate_group_is_singleton",
        "candidate_interpretation",
    ]
    collisions = sorted(set(merge_columns) & set(review_copy.columns))
    if collisions:
        raise ValueError(
            "Review-status input already contains output columns; refusing to "
            f"overwrite them: {collisions}"
        )
    right = assignment.rename(columns={"night": "origin_night"})[
        ["origin_night", "trk_sub", "linkage_id"] + merge_columns
    ]
    merged = review_copy.merge(
        right,
        how="left",
        on=["origin_night", "trk_sub", "linkage_id"],
        validate="one_to_one",
        indicator=True,
        sort=False,
    )
    matched = merged["_merge"].eq("both")
    if int(matched.sum()) != EXPECTED_LINK_COUNT:
        raise AssertionError(
            "Review-status merge must match every one of the 58 retained links "
            f"exactly once; matched {int(matched.sum())}"
        )
    merged["linear_motion_candidate_group_membership_status"] = np.where(
        matched, "retained_58_assigned", "not_in_retained_58"
    )
    merged = merged.sort_values("_review_row_order", kind="stable").drop(
        columns=["_review_row_order", "_merge"]
    )
    atomic_csv(merged, output_path)
    source_hash_after = sha256_file(review_path)
    if source_hash_before != source_hash_after:
        raise AssertionError("Review-status source changed during read-only merge")
    return {
        "source_path": str(review_path.resolve()),
        "source_sha256": source_hash_before,
        "source_row_count": int(len(review)),
        "output_path": str(output_path.resolve()),
        "output_sha256": sha256_file(output_path),
        "output_row_count": int(len(merged)),
        "retained_58_rows_assigned": int(matched.sum()),
        "rows_not_in_retained_58": int((~matched).sum()),
        "source_preserved": True,
    }


def build_summary(
    *,
    links_path: Path,
    links: pd.DataFrame,
    pair_path: Path,
    pairs: pd.DataFrame,
    assignment_path: Path,
    assignments: pd.DataFrame,
    pair_accounting: dict[str, int],
    group_summaries: list[dict[str, Any]],
    review_merge: dict[str, Any] | None,
) -> dict[str, Any]:
    component_sizes = [group["link_count"] for group in group_summaries]
    non_singletons = [group for group in group_summaries if group["link_count"] > 1]
    primary_edges = pairs[pairs["is_primary_conservative_edge"].astype(bool)]
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": utc_now(),
        "analysis_name": "retained_unknown_linear_motion_repeat_candidate_screen",
        "candidate_semantics": {
            "disclaimer": CANDIDATE_DISCLAIMER,
            "permitted_interpretation": (
                "A component is a screening candidate for repeated observations "
                "under a local linear-motion approximation."
            ),
            "prohibited_interpretations": [
                "confirmed independent object count",
                "confirmed cross-night orbital linkage",
                "MPC/JPL identity or designation",
                "replacement for orbit fitting or external identification",
            ],
            "transitive_closure_note": (
                "Connected components use transitive closure; two members in the "
                "same component need not themselves share a direct qualifying edge."
            ),
        },
        "inputs": {
            "links_path": str(links_path.resolve()),
            "links_sha256": sha256_file(links_path),
            "links_row_count": int(len(links)),
            "expected_links_asserted": EXPECTED_LINK_COUNT,
            "review_status_merge": review_merge,
        },
        "method": {
            "eligible_pair_rule": "all unordered pairs with 0 < delta_t <= 7 days",
            "maximum_delta_days": MAX_DELTA_DAYS,
            "time_coordinate": "median_mjd",
            "position_coordinate": "median_ra_deg, median_dec_deg (ICRS degrees)",
            "speed_unit": "arcsec per hour",
            "direction_convention": (
                "lin_dir_deg = atan2(north, east) modulo 360 degrees; 0 degrees "
                "is increasing RA*cos(Dec), 90 degrees is increasing Dec"
            ),
            "production_convention_source": (
                "heliolincrr/orbit_confirm_links.py linear_motion_test"
            ),
            "ra_wrap_handling": (
                "RA differences use the shortest signed wrap in [-180, 180) degrees"
            ),
            "forward_test": (
                "Propagate the earlier link by +delta_t using its own speed and "
                "direction in the tangent plane centered at the earlier position."
            ),
            "backward_test": (
                "Propagate the later link by -delta_t using its own speed and "
                "direction in the tangent plane centered at the later position."
            ),
            "residual_metric": (
                "Euclidean east/north tangent-plane separation in degrees; the "
                "edge test uses max(forward residual, backward residual)."
            ),
        },
        "primary_graph_definition": {
            "fixed_not_cli_overridable": True,
            "threshold": PRIMARY_THRESHOLD.as_dict(),
            "edge_rule": (
                "max(two-way residual) <= 0.03 deg AND speed difference <= "
                "2 arcsec/hour AND direction difference <= 5 deg"
            ),
            "node_definition": "one retained link row (58 nodes total)",
            "component_definition": "undirected connected component of primary edges",
            "includes_same_night_edges": True,
        },
        "pair_accounting": pair_accounting,
        "same_night_vs_cross_night": pair_class_summary(pairs),
        "threshold_sensitivity": [
            threshold_summary(links, pairs, threshold) for threshold in THRESHOLDS
        ],
        "primary_results": {
            "node_count": EXPECTED_LINK_COUNT,
            "primary_edge_count": int(len(primary_edges)),
            "same_night_primary_edge_count": int(primary_edges["same_night"].sum()),
            "cross_night_primary_edge_count": int(
                (~primary_edges["same_night"]).sum()
            ),
            "component_count": len(group_summaries),
            "non_singleton_component_count": len(non_singletons),
            "singleton_component_count": sum(size == 1 for size in component_sizes),
            "links_in_non_singleton_components": sum(
                size for size in component_sizes if size > 1
            ),
            "largest_component_size": max(component_sizes) if component_sizes else 0,
            "component_size_distribution": {
                str(size): count
                for size, count in sorted(Counter(component_sizes).items())
            },
            "groups": group_summaries,
            "non_singleton_groups": non_singletons,
        },
        "outputs": {
            "pair_metrics_csv": {
                "path": str(pair_path.resolve()),
                "sha256": sha256_file(pair_path),
                "row_count": int(len(pairs)),
            },
            "group_assignments_csv": {
                "path": str(assignment_path.resolve()),
                "sha256": sha256_file(assignment_path),
                "row_count": int(len(assignments)),
            },
        },
    }


def main() -> None:
    args = parse_args()
    if args.review_status_output is not None and args.review_status is None:
        raise ValueError("--review-status-output requires --review-status")

    links_path = args.links.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    pair_path = output_dir / "unknown_cross_night_pair_metrics.csv"
    assignment_path = (
        output_dir / "unknown_linear_motion_candidate_group_assignments.csv"
    )
    summary_path = output_dir / "unknown_linear_motion_candidate_group_summary.json"

    links = validate_links(links_path)
    pairs, pair_accounting = build_pair_metrics(links)
    assignments, group_summaries = build_primary_assignments(links, pairs)

    atomic_csv(pairs, pair_path)
    atomic_csv(assignments, assignment_path)

    review_merge: dict[str, Any] | None = None
    if args.review_status is not None:
        review_output = (
            args.review_status_output.resolve()
            if args.review_status_output is not None
            else output_dir
            / "review_and_mpc_status_with_linear_motion_candidate_groups.csv"
        )
        review_merge = merge_review_status(
            args.review_status.resolve(), review_output, assignments
        )

    summary = build_summary(
        links_path=links_path,
        links=links,
        pair_path=pair_path,
        pairs=pairs,
        assignment_path=assignment_path,
        assignments=assignments,
        pair_accounting=pair_accounting,
        group_summaries=group_summaries,
        review_merge=review_merge,
    )
    atomic_json(summary, summary_path)

    primary = summary["primary_results"]
    print(
        "Completed linear-motion candidate screen: "
        f"{pair_accounting['eligible_positive_pairs_within_7_days']} eligible pairs, "
        f"{primary['primary_edge_count']} primary edges, "
        f"{primary['component_count']} components "
        f"({primary['non_singleton_component_count']} non-singletons)."
    )
    print(CANDIDATE_DISCLAIMER)
    print(f"Pair metrics: {pair_path}")
    print(f"Assignments: {assignment_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
