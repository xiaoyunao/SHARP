#!/usr/bin/env python3
"""Generate Fig. 7 from frozen nightly status and coverage products.

The required L2 and quality information comes only from ``night_status.csv``;
the script refuses to infer it from the live filesystem.  Raw counts, unique
fields, cumulative open-shutter time, and cumulative HEALPix area are
cross-checked against the independently frozen coverage table.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from figure_styles import (
    STATISTICS_COLORS,
    apply_statistics_style,
    save_pdf_png,
    style_statistics_axis,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NIGHT_STATUS = PROJECT_ROOT / "snapshot" / "tables" / "night_status.csv"
DEFAULT_COVERAGE = PROJECT_ROOT / "snapshot" / "coverage" / "coverage_by_night.csv"
DEFAULT_COVERAGE_SUMMARY = PROJECT_ROOT / "snapshot" / "coverage" / "coverage_summary.json"
DEFAULT_OUTPUT = PROJECT_ROOT / "figures" / "fig07_nightly_exposure_history"
DEFAULT_FIGURE_DATA = PROJECT_ROOT / "figure_data" / "fig07_nightly_exposure_history.csv"


STATUS_COLUMNS = {
    "night",
    "raw_science_n",
    "l2_mp_n",
    "quality_code",
    "quality_reason",
    "primary_science_included",
    "unknown_science_included",
}
COVERAGE_COLUMNS = {
    "night",
    "raw_exposure_n",
    "unique_field_n",
    "cumulative_exposure_n",
    "cumulative_open_shutter_hours",
    "cumulative_unique_field_n",
    "cumulative_healpix_pixel_n",
    "cumulative_area_deg2",
}


QUALITY_STYLES = {
    "no_strict_raw": {
        "label": "No strict raw frame",
        "facecolor": "#eeeeee",
        "edgecolor": "#d0d0d0",
        "hatch": "",
        "alpha": 0.62,
        "zorder": -8,
    },
    "l2_incomplete": {
        "label": "Raw present; L2 missing/incomplete",
        "facecolor": "#c9c9c9",
        "edgecolor": "#8f8f8f",
        "hatch": "////",
        "alpha": 0.48,
        "zorder": -6,
    },
    "quality_flagged": {
        "label": "Other quality/quarantine flag",
        "facecolor": "#dedede",
        "edgecolor": "#999999",
        "hatch": "\\\\",
        "alpha": 0.34,
        "zorder": -5,
    },
    "primary_excluded": {
        "label": "Primary-science excluded",
        "facecolor": "#9c9c9c",
        "edgecolor": "#5f5f5f",
        "hatch": "xx",
        "alpha": 0.48,
        "zorder": -4,
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
        raise ValueError(f"{name} contains invalid boolean values: {values[invalid].head().tolist()}")
    return normalized.map(mapping).astype(bool)


def validate_integer_count(frame: pd.DataFrame, column: str) -> None:
    values = pd.to_numeric(frame[column], errors="coerce")
    invalid = values.isna() | (values < 0) | ~np.isclose(values, np.round(values))
    if invalid.any():
        examples = frame.loc[invalid, ["night", column]].head().to_dict("records")
        raise ValueError(f"{column} must be a non-negative integer count: {examples}")
    frame[column] = np.round(values).astype("int64")


def read_inputs(
    night_status_path: Path,
    coverage_path: Path,
    coverage_summary_path: Path,
) -> pd.DataFrame:
    if not night_status_path.exists():
        raise FileNotFoundError(
            f"Frozen night table not found: {night_status_path}\n"
            "Fig. 7 requires per-night L2 counts and quality codes. Build the "
            "frozen night table first; a coverage-only substitute is not publication-safe."
        )
    if not coverage_path.exists():
        raise FileNotFoundError(f"Frozen coverage table not found: {coverage_path}")
    if not coverage_summary_path.exists():
        raise FileNotFoundError(f"Frozen coverage summary not found: {coverage_summary_path}")

    status = pd.read_csv(night_status_path, dtype={"night": "string"}, low_memory=False)
    coverage = pd.read_csv(coverage_path, dtype={"night": "string"}, low_memory=False)
    missing_status = sorted(STATUS_COLUMNS - set(status.columns))
    missing_coverage = sorted(COVERAGE_COLUMNS - set(coverage.columns))
    if missing_status:
        raise ValueError(
            f"{night_status_path} lacks Fig. 7 columns: {', '.join(missing_status)}"
        )
    if missing_coverage:
        raise ValueError(
            f"{coverage_path} lacks Fig. 7 columns: {', '.join(missing_coverage)}"
        )

    for frame, name in [(status, "night_status"), (coverage, "coverage")]:
        frame["night"] = frame["night"].str.strip().str.zfill(8)
        if frame["night"].isna().any() or frame["night"].duplicated().any():
            duplicates = frame.loc[frame["night"].duplicated(False), "night"].head().tolist()
            raise ValueError(f"{name} night keys must be non-null and unique: {duplicates}")

    for column in ["raw_science_n", "l2_mp_n"]:
        validate_integer_count(status, column)
    for column in [
        "raw_exposure_n",
        "unique_field_n",
        "cumulative_exposure_n",
        "cumulative_unique_field_n",
        "cumulative_healpix_pixel_n",
    ]:
        validate_integer_count(coverage, column)
    for column in ["cumulative_open_shutter_hours", "cumulative_area_deg2"]:
        coverage[column] = pd.to_numeric(coverage[column], errors="coerce")
        if coverage[column].isna().any() or (coverage[column] < 0).any():
            raise ValueError(f"{coverage_path}: invalid non-negative numeric column {column}")
        if (coverage.sort_values("night")[column].diff().dropna() < -1e-9).any():
            raise ValueError(f"{coverage_path}: {column} is not monotonic non-decreasing")

    status["primary_science_included"] = parse_bool(
        status["primary_science_included"], "primary_science_included"
    )
    status["unknown_science_included"] = parse_bool(
        status["unknown_science_included"], "unknown_science_included"
    )
    status["date"] = pd.to_datetime(status["night"], format="%Y%m%d", errors="raise")
    status = status.sort_values("date").reset_index(drop=True)
    expected = pd.date_range(status["date"].min(), status["date"].max(), freq="D")
    if len(status) != len(expected) or not status["date"].reset_index(drop=True).equals(
        pd.Series(expected, name="date")
    ):
        present = set(status["date"])
        missing = [date.strftime("%Y%m%d") for date in expected if date not in present]
        raise ValueError(
            "night_status must contain one row for every calendar night; "
            f"missing examples={missing[:8]}"
        )

    coverage["date"] = pd.to_datetime(coverage["night"], format="%Y%m%d", errors="raise")
    merged = status.merge(
        coverage.drop(columns="night"),
        on="date",
        how="left",
        validate="one_to_one",
    )
    nightly_zero = ["raw_exposure_n", "unique_field_n"]
    cumulative = [
        "cumulative_exposure_n",
        "cumulative_open_shutter_hours",
        "cumulative_unique_field_n",
        "cumulative_healpix_pixel_n",
        "cumulative_area_deg2",
    ]
    merged[nightly_zero] = merged[nightly_zero].fillna(0)
    merged[cumulative] = merged[cumulative].ffill().fillna(0)
    for column in nightly_zero + [
        "cumulative_exposure_n",
        "cumulative_unique_field_n",
        "cumulative_healpix_pixel_n",
    ]:
        merged[column] = np.round(merged[column]).astype("int64")

    mismatch = merged["raw_science_n"].ne(merged["raw_exposure_n"])
    if mismatch.any():
        examples = merged.loc[
            mismatch, ["night", "raw_science_n", "raw_exposure_n"]
        ].head().to_dict("records")
        raise ValueError(
            "night_status raw counts do not match the frozen coverage raw counts: "
            f"{examples}"
        )
    if int(merged["raw_science_n"].sum()) != int(merged["cumulative_exposure_n"].iloc[-1]):
        raise ValueError("final cumulative exposure count does not close to nightly raw counts")

    summary = json.loads(coverage_summary_path.read_text(encoding="utf-8"))
    summary_required = {
        "healpix_nside",
        "raw_exposure_n",
        "open_shutter_hours",
        "unique_area_deg2",
    }
    missing_summary = sorted(summary_required - set(summary))
    if missing_summary:
        raise ValueError(
            f"{coverage_summary_path} lacks Fig. 7 keys: {', '.join(missing_summary)}"
        )
    assert_pairs = [
        ("raw_exposure_n", float(merged["cumulative_exposure_n"].iloc[-1]), 0.0),
        (
            "open_shutter_hours",
            float(merged["cumulative_open_shutter_hours"].iloc[-1]),
            1e-8,
        ),
        ("unique_area_deg2", float(merged["cumulative_area_deg2"].iloc[-1]), 1e-8),
    ]
    for key, table_value, tolerance in assert_pairs:
        summary_value = float(summary[key])
        if not np.isclose(summary_value, table_value, rtol=0.0, atol=tolerance):
            raise ValueError(
                f"coverage summary/table mismatch for {key}: {summary_value} != {table_value}"
            )
    healpix_nside = int(summary["healpix_nside"])
    if healpix_nside <= 0:
        raise ValueError(f"{coverage_summary_path}: healpix_nside must be positive")
    merged.attrs["healpix_nside"] = healpix_nside

    merged["l2_gap_n"] = merged["raw_science_n"] - merged["l2_mp_n"]
    # A negative gap would mean more strict L2 catalogs than strict raw files,
    # which needs upstream reconciliation rather than graphical suppression.
    if (merged["l2_gap_n"] < 0).any():
        examples = merged.loc[
            merged["l2_gap_n"] < 0, ["night", "raw_science_n", "l2_mp_n"]
        ].head().to_dict("records")
        raise ValueError(f"strict L2 count exceeds strict raw count: {examples}")

    primary = merged["primary_science_included"]
    quality_code = merged["quality_code"].fillna("").astype(str)
    merged["quality_plot_class"] = "quality_ok"
    merged.loc[primary & merged["raw_science_n"].eq(0), "quality_plot_class"] = "no_strict_raw"
    merged.loc[primary & merged["l2_gap_n"].gt(0), "quality_plot_class"] = "l2_incomplete"
    other_flag = (
        primary
        & merged["quality_plot_class"].eq("quality_ok")
        & ~quality_code.eq("quality_ok")
    )
    merged.loc[other_flag, "quality_plot_class"] = "quality_flagged"
    merged.loc[~primary, "quality_plot_class"] = "primary_excluded"
    return merged


def shade_quality_intervals(axis, frame: pd.DataFrame) -> None:
    """Shade contiguous status runs, retaining hatches for grayscale legibility."""

    dates = frame["date"]
    categories = frame["quality_plot_class"]
    for category, style in QUALITY_STYLES.items():
        mask = categories.eq(category).to_numpy()
        if not mask.any():
            continue
        starts = np.flatnonzero(mask & np.r_[True, ~mask[:-1]])
        stops = np.flatnonzero(mask & np.r_[~mask[1:], True])
        for start, stop in zip(starts, stops):
            axis.axvspan(
                dates.iloc[start] - pd.Timedelta(hours=12),
                dates.iloc[stop] + pd.Timedelta(hours=12),
                facecolor=style["facecolor"],
                edgecolor=style["edgecolor"],
                hatch=style["hatch"],
                alpha=style["alpha"],
                linewidth=0.0 if not style["hatch"] else 0.5,
                zorder=style["zorder"],
            )


def panel_label(axis, text: str) -> None:
    axis.text(
        0.012,
        0.965,
        text,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=21,
        fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.76, "pad": 2.0},
        zorder=20,
    )


def make_figure(frame: pd.DataFrame):
    apply_statistics_style()
    figure, axes = plt.subplots(
        3,
        1,
        figsize=(18, 13.5),
        sharex=True,
        gridspec_kw={"height_ratios": [1.25, 1.0, 1.25], "hspace": 0.08},
    )
    raw_axis, field_axis, cumulative_axis = axes
    dates = frame["date"]
    for axis in axes:
        shade_quality_intervals(axis, frame)
        style_statistics_axis(axis, tick_fontsize=16)

    raw_axis.bar(
        dates,
        frame["raw_science_n"],
        width=0.82,
        color=STATISTICS_COLORS["blue"],
        alpha=0.74,
        edgecolor="none",
        label="Strict raw science files",
        zorder=3,
    )
    raw_axis.plot(
        dates,
        frame["l2_mp_n"],
        color=STATISTICS_COLORS["orange"],
        linewidth=1.35,
        marker="o",
        markersize=2.4,
        markeredgewidth=0,
        label="Strict L2 catalogs",
        zorder=5,
    )
    raw_axis.set_ylabel("Files per night", fontsize=18)
    raw_axis.set_ylim(bottom=0)
    raw_axis.legend(loc="upper right", fontsize=14, ncol=2)
    panel_label(raw_axis, "(a) Acquired raw and L2 products")

    field_axis.bar(
        dates,
        frame["unique_field_n"],
        width=0.82,
        color=STATISTICS_COLORS["green"],
        alpha=0.75,
        edgecolor="none",
        zorder=3,
    )
    field_axis.set_ylabel("Distinct fields\nper night", fontsize=18)
    field_axis.set_ylim(bottom=0)
    panel_label(field_axis, "(b) Nightly field coverage")

    time_line = cumulative_axis.plot(
        dates,
        frame["cumulative_open_shutter_hours"],
        color=STATISTICS_COLORS["blue"],
        linewidth=2.6,
        label="Cumulative open-shutter time",
        zorder=6,
    )[0]
    cumulative_axis.set_ylabel("Open-shutter time (h)", fontsize=18, color=STATISTICS_COLORS["blue"])
    cumulative_axis.tick_params(axis="y", colors=STATISTICS_COLORS["blue"])
    cumulative_axis.set_ylim(bottom=0)
    area_axis = cumulative_axis.twinx()
    area_line = area_axis.plot(
        dates,
        frame["cumulative_area_deg2"],
        color=STATISTICS_COLORS["red"],
        linewidth=2.4,
        linestyle="--",
        label="Cumulative unique HEALPix area",
        zorder=7,
    )[0]
    area_axis.set_ylabel(
        r"Unique area (deg$^2$)", fontsize=18, color=STATISTICS_COLORS["red"]
    )
    area_axis.tick_params(axis="y", colors=STATISTICS_COLORS["red"], labelsize=16)
    area_axis.set_ylim(bottom=0)
    area_axis.grid(False)
    cumulative_axis.legend(
        [time_line, area_line],
        [time_line.get_label(), area_line.get_label()],
        loc="upper left",
        bbox_to_anchor=(0.40, 0.98),
        fontsize=14,
        ncol=1,
    )
    panel_label(cumulative_axis, "(c) Cumulative exposure and sky coverage")

    locator = mdates.MonthLocator(interval=1)
    cumulative_axis.xaxis.set_major_locator(locator)
    cumulative_axis.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    cumulative_axis.tick_params(axis="x", rotation=32, labelsize=15)
    cumulative_axis.set_xlabel("Observing night (local-date label)", fontsize=19)
    cumulative_axis.set_xlim(
        dates.min() - pd.Timedelta(days=2), dates.max() + pd.Timedelta(days=2)
    )

    quality_handles = []
    for category, style in QUALITY_STYLES.items():
        if frame["quality_plot_class"].eq(category).any():
            quality_handles.append(
                Patch(
                    facecolor=style["facecolor"],
                    edgecolor=style["edgecolor"],
                    hatch=style["hatch"],
                    label=style["label"],
                    alpha=max(style["alpha"], 0.65),
                )
            )
    if quality_handles:
        figure.legend(
            handles=quality_handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.015),
            ncol=min(4, len(quality_handles)),
            fontsize=13,
            frameon=False,
        )
    figure.suptitle(
        "Nightly survey exposure history",
        fontsize=27,
        fontweight="bold",
        y=0.992,
    )
    figure.text(
        0.5,
        0.947,
        "Night grain; raw/L2 are files, fields are distinct field IDs per night, "
        f"and area is the cumulative NSIDE={frame.attrs['healpix_nside']} HEALPix union.",
        ha="center",
        fontsize=14,
        color=STATISTICS_COLORS["mid_grey"],
    )
    figure.subplots_adjust(left=0.075, right=0.915, top=0.900, bottom=0.115)
    return figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--night-status", type=Path, default=DEFAULT_NIGHT_STATUS)
    parser.add_argument("--coverage", type=Path, default=DEFAULT_COVERAGE)
    parser.add_argument("--coverage-summary", type=Path, default=DEFAULT_COVERAGE_SUMMARY)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure-data", type=Path, default=DEFAULT_FIGURE_DATA)
    args = parser.parse_args()

    figure_data = read_inputs(args.night_status, args.coverage, args.coverage_summary)
    args.figure_data.parent.mkdir(parents=True, exist_ok=True)
    export_columns = [
        "night",
        "date",
        "raw_science_n",
        "l2_mp_n",
        "l2_gap_n",
        "raw_exposure_n",
        "unique_field_n",
        "cumulative_exposure_n",
        "cumulative_open_shutter_hours",
        "cumulative_unique_field_n",
        "cumulative_healpix_pixel_n",
        "cumulative_area_deg2",
        "quality_code",
        "quality_reason",
        "quality_plot_class",
        "primary_science_included",
        "unknown_science_included",
    ]
    figure_data[export_columns].to_csv(args.figure_data, index=False)
    figure = make_figure(figure_data)
    png_path, pdf_path = save_pdf_png(figure, args.output_base)
    print(f"wrote {pdf_path}")
    print(f"wrote {png_path}")
    print(f"wrote {args.figure_data}")


if __name__ == "__main__":
    main()
