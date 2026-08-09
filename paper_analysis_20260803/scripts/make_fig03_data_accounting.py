#!/usr/bin/env python3
"""Generate Fig. 3 from the frozen per-night accounting table.

The figure deliberately does not draw a width-scaled Sankey.  Adjacent stages
mix files, detections, tracklets, linkages, ADES rows, and nights, so implying
mass conservation between every pair would be incorrect.  Instead, the top
row gives explicit, closed night/file accounting and the lower panels give
unit-labelled processing-order cards.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch

from figure_styles import (
    STATISTICS_COLORS,
    apply_statistics_style,
    save_pdf_png,
    style_statistics_axis,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NIGHT_STATUS = PROJECT_ROOT / "snapshot" / "tables" / "night_status.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "figures" / "fig03_data_accounting"
DEFAULT_FIGURE_DATA = PROJECT_ROOT / "figure_data" / "fig03_data_accounting.csv"


COUNT_COLUMNS = {
    "raw_science_n",
    "l2_mp_n",
    "known_predicted_n",
    "known_match1_n",
    "known_mask15_n",
    "known_ades_n",
    "l2_detection_n",
    "gaia_survivor_n",
    "tracklet_n",
    "link_n",
    "orbit_fit_n",
    "orbit_is_good_n",
    "unknown_n",
    "review_real_n",
    "submit_real_n",
    "audit_real_n",
}
REQUIRED_COLUMNS = COUNT_COLUMNS | {
    "night",
    "primary_science_included",
    "unknown_science_included",
    "known_mpc_state",
    "unknown_mpc_state",
}


UNIT_COLORS = {
    "file": "#d8e5f2",
    "detection": "#f8ddbe",
    "tracklet": "#d9ead3",
    "linkage": "#e7d6e8",
    "ADES row": "#f8efc7",
    "night": "#d9e7e5",
}


def parse_bool(values: pd.Series, name: str) -> pd.Series:
    """Parse a required boolean column without silently treating unknowns false."""

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
        examples = values.loc[invalid].head(5).tolist()
        raise ValueError(f"{name} contains non-boolean values, for example {examples!r}")
    return normalized.map(mapping).astype(bool)


def validate_count(frame: pd.DataFrame, column: str) -> None:
    values = pd.to_numeric(frame[column], errors="coerce")
    invalid = values.isna() | (values < 0) | ~np.isclose(values, np.round(values))
    if invalid.any():
        examples = frame.loc[invalid, ["night", column]].head(5).to_dict("records")
        raise ValueError(f"{column} must contain non-negative integer counts: {examples}")
    frame[column] = np.round(values).astype("int64")


def read_night_status(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Frozen night table not found: {path}\n"
            "Run build_night_tables.py on the frozen manifests first; Fig. 3 "
            "will not substitute live data or hard-coded headline counts."
        )
    frame = pd.read_csv(path, dtype={"night": "string"}, low_memory=False)
    missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(
            f"{path} is missing required columns for Fig. 3: {', '.join(missing)}"
        )
    frame["night"] = frame["night"].str.strip().str.zfill(8)
    if frame["night"].isna().any() or frame["night"].duplicated().any():
        duplicates = frame.loc[frame["night"].duplicated(False), "night"].head().tolist()
        raise ValueError(f"night must be non-null and unique; duplicates={duplicates}")
    for column in sorted(COUNT_COLUMNS):
        validate_count(frame, column)
    frame["primary_science_included"] = parse_bool(
        frame["primary_science_included"], "primary_science_included"
    )
    frame["unknown_science_included"] = parse_bool(
        frame["unknown_science_included"], "unknown_science_included"
    )
    impossible = frame["unknown_science_included"] & ~frame["primary_science_included"]
    if impossible.any():
        raise ValueError(
            "unknown_science_included cannot be true on a primary-science-excluded night"
        )
    return frame.sort_values("night").reset_index(drop=True)


def append_row(rows: list[dict[str, object]], **kwargs: object) -> None:
    base = {
        "panel": "",
        "order": 0,
        "stage_key": "",
        "stage_label": "",
        "unit": "",
        "scope": "",
        "count": 0,
        "parent_total": np.nan,
        "closure_delta": np.nan,
        "source_column": "",
    }
    base.update(kwargs)
    rows.append(base)


def build_figure_data(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[dict], list[dict]]:
    rows: list[dict[str, object]] = []
    primary = frame["primary_science_included"]
    unknown = frame["unknown_science_included"]

    night_masks = [
        ("primary_excluded", "Primary-science excluded", ~primary),
        ("unknown_quarantine", "Unknown-branch quarantine", primary & ~unknown),
        ("no_strict_raw", "No strict raw frame", unknown & frame["raw_science_n"].eq(0)),
        (
            "raw_without_l2",
            "Raw present, L2 absent",
            unknown & frame["raw_science_n"].gt(0) & frame["l2_mp_n"].eq(0),
        ),
        (
            "analysis_ready",
            "Unknown scope with raw + L2",
            unknown & frame["raw_science_n"].gt(0) & frame["l2_mp_n"].gt(0),
        ),
    ]
    night_total = len(frame)
    night_sum = sum(int(mask.sum()) for _, _, mask in night_masks)
    if night_sum != night_total:
        raise ValueError(
            f"Night-accounting categories do not close: {night_sum} != {night_total}"
        )
    for order, (key, label, mask) in enumerate(night_masks):
        append_row(
            rows,
            panel="closure_nights",
            order=order,
            stage_key=key,
            stage_label=label,
            unit="night",
            scope="all calendar nights in frozen night_status",
            count=int(mask.sum()),
            parent_total=night_total,
            closure_delta=night_total - night_sum,
            source_column=(
                "night;raw_science_n;l2_mp_n;primary_science_included;"
                "unknown_science_included"
            ),
        )

    for order, (column, label) in enumerate(
        [("raw_science_n", "Strict raw science frames"), ("l2_mp_n", "Strict L2 catalogs")]
    ):
        total = int(frame[column].sum())
        included = int(frame.loc[primary, column].sum())
        excluded = int(frame.loc[~primary, column].sum())
        if included + excluded != total:
            raise AssertionError(f"{column} inclusion accounting failed")
        for segment_order, (key, segment_label, value) in enumerate(
            [
                ("primary_included", "Primary-science included", included),
                ("primary_excluded", "Primary-science excluded", excluded),
            ]
        ):
            append_row(
                rows,
                panel=f"closure_{column}",
                order=segment_order,
                stage_key=f"{column}_{key}",
                stage_label=f"{label}: {segment_label}",
                unit="file",
                scope="all frozen nights",
                count=value,
                parent_total=total,
                closure_delta=total - included - excluded,
                source_column=column,
            )

    known_specs = [
        ("raw_science", "Strict raw science", "raw_science_n", "file"),
        ("l2_catalog", "Strict L2 catalogs", "l2_mp_n", "file"),
        ("known_prediction", "Known-object predictions", "known_predicted_n", "detection"),
        ("known_match1", "1 arcsec known matches", "known_match1_n", "detection"),
        ("known_mask15", "1.5 arcsec subtraction mask", "known_mask15_n", "detection"),
        ("known_ades", "Known ADES observations", "known_ades_n", "ADES row"),
    ]
    known_cards: list[dict] = []
    for order, (key, label, column, unit) in enumerate(known_specs):
        included = int(frame.loc[primary, column].sum())
        excluded = int(frame.loc[~primary, column].sum())
        card = {
            "key": key,
            "label": label,
            "unit": unit,
            "count": included,
            "excluded": excluded,
        }
        known_cards.append(card)
        append_row(
            rows,
            panel="known_branch",
            order=order,
            stage_key=key,
            stage_label=label,
            unit=unit,
            scope="primary-science-included nights",
            count=included,
            parent_total=included + excluded,
            closure_delta=0,
            source_column=column,
        )
        append_row(
            rows,
            panel="known_branch_excluded",
            order=order,
            stage_key=f"{key}_excluded",
            stage_label=f"{label}: primary-science excluded",
            unit=unit,
            scope="primary-science-excluded nights",
            count=excluded,
            parent_total=included + excluded,
            closure_delta=0,
            source_column=column,
        )
    known_reply = int(
        (primary & frame["known_mpc_state"].astype("string").eq("reply_received")).sum()
    )
    known_reply_excluded = int(
        (~primary & frame["known_mpc_state"].astype("string").eq("reply_received")).sum()
    )
    known_cards.append(
        {
            "key": "known_reply",
            "label": "Known MPC reply received",
            "unit": "night",
            "count": known_reply,
            "excluded": known_reply_excluded,
        }
    )
    append_row(
        rows,
        panel="known_branch",
        order=len(known_specs),
        stage_key="known_reply",
        stage_label="Known MPC reply received",
        unit="night",
        scope="primary-science-included nights",
        count=known_reply,
        parent_total=known_reply + known_reply_excluded,
        closure_delta=0,
        source_column="known_mpc_state",
    )

    unknown_specs = [
        ("l2_detection", "L2 catalog sources", "l2_detection_n", "detection"),
        ("gaia_survivor", "Gaia-mask survivors", "gaia_survivor_n", "detection"),
        ("tracklet", "Two-point tracklets", "tracklet_n", "tracklet"),
        ("link", "Shared-endpoint links", "link_n", "linkage"),
        ("orbit_fit", "Numerical fit_ok", "orbit_fit_n", "linkage"),
        (
            "orbit_is_good",
            "Thresholded is_good (diagnostic)",
            "orbit_is_good_n",
            "linkage",
        ),
        ("post_known", "Post-known catalog", "unknown_n", "linkage"),
        ("initial_review", "Human review marked real", "review_real_n", "linkage"),
        ("submission_selected", "Submission-selected", "submit_real_n", "linkage"),
        ("posthoc_audit", "Post-audit retained", "audit_real_n", "linkage"),
    ]
    unknown_cards: list[dict] = []
    for order, (key, label, column, unit) in enumerate(unknown_specs):
        included = int(frame.loc[unknown, column].sum())
        primary_excluded = int(frame.loc[~primary, column].sum())
        quarantine = int(frame.loc[primary & ~unknown, column].sum())
        excluded = primary_excluded + quarantine
        unknown_cards.append(
            {
                "key": key,
                "label": label,
                "unit": unit,
                "count": included,
                "excluded": excluded,
            }
        )
        for group_order, (group, mask) in enumerate(
            [
                ("unknown_included", unknown),
                ("primary_excluded", ~primary),
                ("unknown_quarantine", primary & ~unknown),
            ]
        ):
            append_row(
                rows,
                panel="unknown_branch" if group_order == 0 else "unknown_branch_excluded",
                order=order,
                stage_key=key if group_order == 0 else f"{key}_{group}",
                stage_label=label if group_order == 0 else f"{label}: {group.replace('_', ' ')}",
                unit=unit,
                scope=group.replace("_", " "),
                count=int(frame.loc[mask, column].sum()),
                parent_total=included + excluded,
                closure_delta=0,
                source_column=column,
            )
    unknown_reply = int(
        (unknown & frame["unknown_mpc_state"].astype("string").eq("reply_received")).sum()
    )
    unknown_reply_excluded = int(
        (~unknown & frame["unknown_mpc_state"].astype("string").eq("reply_received")).sum()
    )
    unknown_cards.append(
        {
            "key": "unknown_reply",
            "label": "Unknown MPC reply received",
            "unit": "night",
            "count": unknown_reply,
            "excluded": unknown_reply_excluded,
        }
    )
    append_row(
        rows,
        panel="unknown_branch",
        order=len(unknown_specs),
        stage_key="unknown_reply",
        stage_label="Unknown MPC reply received",
        unit="night",
        scope="unknown-science-included nights",
        count=unknown_reply,
        parent_total=unknown_reply + unknown_reply_excluded,
        closure_delta=0,
        source_column="unknown_mpc_state",
    )
    return pd.DataFrame(rows), known_cards, unknown_cards


def draw_closure_bar(
    ax,
    values: Iterable[int],
    labels: Iterable[str],
    colors: Iterable[str],
    hatches: Iterable[str],
    title: str,
    unit: str,
) -> None:
    values = list(values)
    labels = list(labels)
    left = 0
    for value, label, color, hatch in zip(values, labels, colors, hatches):
        ax.barh(
            [0],
            [value],
            left=[left],
            height=0.48,
            color=color,
            edgecolor=STATISTICS_COLORS["ink"],
            linewidth=0.65,
            hatch=hatch,
            label=f"{label} ({value:,})",
        )
        if value and value >= max(values) * 0.08:
            ax.text(left + value / 2, 0, f"{value:,}", ha="center", va="center", fontsize=14)
        left += value
    ax.set_xlim(0, max(left * 1.02, 1))
    ax.set_yticks([])
    ax.set_xlabel(unit, fontsize=15)
    ax.set_title(f"{title}\nTotal = {left:,} {unit}", fontsize=18, fontweight="bold")
    ax.legend(
        fontsize=11,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.30),
        ncol=2 if len(labels) > 2 else 1,
    )
    style_statistics_axis(ax, tick_fontsize=13)


def unit_display(unit: str, count: int) -> str:
    if unit == "ADES row":
        return "ADES row" if count == 1 else "ADES rows"
    return unit if count == 1 else f"{unit}s"


def draw_stage_cards(ax, cards: list[dict], panel: str, scope: str) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.0, 1.02, panel, fontsize=22, fontweight="bold", ha="left", va="bottom")
    ax.text(
        0.0,
        0.982,
        scope,
        fontsize=12.5,
        color=STATISTICS_COLORS["mid_grey"],
        ha="left",
        va="top",
    )
    n = len(cards)
    top, bottom = 0.905, 0.035
    box_h = min(0.082, (top - bottom) / (n * 1.18))
    centers = np.linspace(top - box_h / 2, bottom + box_h / 2, n)
    center_by_key = {card["key"]: center for card, center in zip(cards, centers)}
    for index, (card, center_y) in enumerate(zip(cards, centers)):
        x, width = 0.08, 0.84
        diagnostic = card["key"] == "orbit_is_good"
        patch = FancyBboxPatch(
            (x, center_y - box_h / 2),
            width,
            box_h,
            boxstyle="round,pad=0.010,rounding_size=0.018",
            facecolor=UNIT_COLORS[card["unit"]],
            edgecolor=STATISTICS_COLORS["ink"],
            linewidth=1.25 if diagnostic else 1.0,
            linestyle="--" if diagnostic else "-",
        )
        ax.add_patch(patch)
        count = int(card["count"])
        ax.text(
            x + 0.02,
            center_y,
            card["label"],
            fontsize=14.5,
            fontweight="bold",
            ha="left",
            va="center",
        )
        ax.text(
            x + width - 0.02,
            center_y + (0.010 if card["excluded"] else 0),
            f"{count:,} {unit_display(card['unit'], count)}",
            fontsize=14.5,
            ha="right",
            va="center",
        )
        if card["excluded"]:
            ax.text(
                x + width - 0.02,
                center_y - 0.020,
                f"outside analysis scope: {int(card['excluded']):,}",
                fontsize=9.8,
                color=STATISTICS_COLORS["mid_grey"],
                ha="right",
                va="center",
            )
        if index < n - 1:
            next_y = centers[index + 1]
            if card["key"] == "orbit_is_good":
                continue
            next_is_diagnostic = cards[index + 1]["key"] == "orbit_is_good"
            ax.add_patch(
                FancyArrowPatch(
                    (0.5, center_y - box_h / 2 - 0.004),
                    (0.5, next_y + box_h / 2 + 0.004),
                    arrowstyle="-|>",
                    mutation_scale=12,
                    linewidth=1.05 if next_is_diagnostic else 1.15,
                    linestyle="--" if next_is_diagnostic else "-",
                    color=STATISTICS_COLORS["mid_grey"],
                )
            )
    if {"orbit_fit", "orbit_is_good", "post_known"}.issubset(center_by_key):
        # Production selected fit_ok + all-non-known; is_good is a parallel
        # threshold diagnostic and must not be drawn as the applied gate.
        ax.add_patch(
            FancyArrowPatch(
                (0.079, center_by_key["orbit_fit"]),
                (0.079, center_by_key["post_known"]),
                connectionstyle="arc3,rad=0.34",
                arrowstyle="-|>",
                mutation_scale=12,
                linewidth=1.35,
                color=STATISTICS_COLORS["blue"],
                zorder=5,
            )
        )
        ax.text(
            0.018,
            (center_by_key["orbit_fit"] + center_by_key["post_known"]) / 2,
            "production selector path",
            fontsize=8.3,
            color=STATISTICS_COLORS["blue"],
            rotation=90,
            ha="center",
            va="center",
        )


def make_figure(frame: pd.DataFrame, known_cards: list[dict], unknown_cards: list[dict]):
    apply_statistics_style()
    figure = plt.figure(figsize=(18, 14.5))
    grid = figure.add_gridspec(
        2,
        6,
        height_ratios=[1.05, 4.0],
        left=0.045,
        right=0.985,
        top=0.875,
        bottom=0.075,
        hspace=0.50,
        wspace=0.70,
    )
    ax_night = figure.add_subplot(grid[0, 0:2])
    ax_raw = figure.add_subplot(grid[0, 2:4])
    ax_l2 = figure.add_subplot(grid[0, 4:6])
    ax_known = figure.add_subplot(grid[1, 0:3])
    ax_unknown = figure.add_subplot(grid[1, 3:6])

    primary = frame["primary_science_included"]
    night_values = [
        int((~primary).sum()),
        int((primary & ~frame["unknown_science_included"]).sum()),
        int((frame["unknown_science_included"] & frame["raw_science_n"].eq(0)).sum()),
        int(
            (
                frame["unknown_science_included"]
                & frame["raw_science_n"].gt(0)
                & frame["l2_mp_n"].eq(0)
            ).sum()
        ),
        int(
            (
                frame["unknown_science_included"]
                & frame["raw_science_n"].gt(0)
                & frame["l2_mp_n"].gt(0)
            ).sum()
        ),
    ]
    draw_closure_bar(
        ax_night,
        night_values,
        ["Primary excluded", "Unknown quarantine", "No strict raw", "Raw; no L2", "Raw + L2"],
        ["#9c9c9c", "#b8b8b8", "#eeeeee", "#cfcfcf", STATISTICS_COLORS["blue"]],
        ["xx", "\\\\", "", "////", ""],
        "(a) Night closure",
        "nights",
    )
    for axis, column, title in [
        (ax_raw, "raw_science_n", "(b) Raw-file closure"),
        (ax_l2, "l2_mp_n", "(c) L2-file closure"),
    ]:
        draw_closure_bar(
            axis,
            [int(frame.loc[primary, column].sum()), int(frame.loc[~primary, column].sum())],
            ["Primary included", "Primary excluded"],
            [STATISTICS_COLORS["blue"], "#9c9c9c"],
            ["", "xx"],
            title,
            "files",
        )

    draw_stage_cards(
        ax_known,
        known_cards,
        "(d) Known-object branch",
        "Counts aggregated over primary-science-included nights",
    )
    draw_stage_cards(
        ax_unknown,
        unknown_cards,
        "(e) Unknown-object branch",
        "Counts aggregated over unknown-science-included nights",
    )
    figure.suptitle(
        "Frozen survey data accounting",
        fontsize=27,
        fontweight="bold",
        y=0.978,
    )
    figure.legend(
        handles=[
            Patch(facecolor=color, edgecolor=STATISTICS_COLORS["ink"], label=unit)
            for unit, color in UNIT_COLORS.items()
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=len(UNIT_COLORS),
        fontsize=12,
        frameon=False,
    )
    figure.text(
        0.5,
        0.006,
        "Arrows show processing order only; the blue bypass is the executed fit_ok selector path, while is_good is a parallel diagnostic. Box size does not encode count.",
        ha="center",
        va="bottom",
        fontsize=13,
        color=STATISTICS_COLORS["mid_grey"],
    )
    return figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--night-status", type=Path, default=DEFAULT_NIGHT_STATUS)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure-data", type=Path, default=DEFAULT_FIGURE_DATA)
    args = parser.parse_args()

    frame = read_night_status(args.night_status)
    figure_data, known_cards, unknown_cards = build_figure_data(frame)
    args.figure_data.parent.mkdir(parents=True, exist_ok=True)
    figure_data.to_csv(args.figure_data, index=False)
    figure = make_figure(frame, known_cards, unknown_cards)
    png_path, pdf_path = save_pdf_png(figure, args.output_base)
    print(f"wrote {pdf_path}")
    print(f"wrote {png_path}")
    print(f"wrote {args.figure_data}")


if __name__ == "__main__":
    main()
