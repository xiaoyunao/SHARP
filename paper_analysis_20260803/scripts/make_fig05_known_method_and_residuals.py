#!/usr/bin/env python3
"""Generate Fig. 5: known-object query geometry and match diagnostics."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle, Polygon, Rectangle

from figure_styles import (
    STATISTICS_COLORS,
    add_statistics_colorbar,
    apply_statistics_style,
    save_pdf_png,
    style_statistics_axis,
)


def wilson(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return np.nan, np.nan
    fraction = success / total
    denominator = 1 + z * z / total
    center = (fraction + z * z / (2 * total)) / denominator
    half = z * np.sqrt(fraction * (1 - fraction) / total + z * z / (4 * total * total)) / denominator
    return center - half, center + half


def geometry_panel(axis) -> None:
    axis.set_aspect("equal")
    axis.set_xlim(-3.05, 3.05)
    axis.set_ylim(-3.05, 3.05)
    detector = np.array([[-2.1, -1.65], [-1.85, 1.72], [2.0, 1.52], [2.18, -1.48]])
    axis.add_patch(
        Polygon(detector, closed=True, facecolor="#dce8f4", edgecolor=STATISTICS_COLORS["blue"], linewidth=2.4)
    )
    axis.add_patch(Circle((0, 0), 2.72, fill=False, linestyle="--", linewidth=2.0, color=STATISTICS_COLORS["orange"]))
    axis.scatter(detector[:, 0], detector[:, 1], marker="s", s=55, color=STATISTICS_COLORS["blue"], zorder=5)
    axis.scatter([0], [0], marker="*", s=180, color=STATISTICS_COLORS["red"], edgecolor="white", linewidth=0.6, zorder=6)
    axis.axvline(0, color=STATISTICS_COLORS["mid_grey"], linewidth=1.0, alpha=0.45)
    axis.text(0.10, -1.05, r"RA = 0$^\circ$ seam", ha="left", va="center", rotation=90, fontsize=13, color=STATISTICS_COLORS["mid_grey"])
    axis.annotate(
        "nominal detector\nrectangle",
        xy=(1.65, 1.22),
        xytext=(1.28, 1.92),
        fontsize=15,
        arrowprops={"arrowstyle": "->", "lw": 1.4, "color": STATISTICS_COLORS["ink"]},
    )
    axis.annotate(
        "query radius = farthest\ncorner + 0.05°",
        xy=(-1.82, -1.75),
        xytext=(-2.86, -0.82),
        fontsize=14,
        arrowprops={"arrowstyle": "->", "lw": 1.4, "color": STATISTICS_COLORS["ink"]},
    )
    axis.scatter([-0.28, 0.28], [0.35, 0.35], s=85, facecolor="white", edgecolor=STATISTICS_COLORS["green"], linewidth=2)
    axis.text(0, 0.62, "extra query centers\n359° and 0°", ha="center", fontsize=15, color=STATISTICS_COLORS["green"])
    axis.text(-2.90, 2.90, "(a) Production query geometry (schematic)", fontsize=19, weight="bold", va="top")
    axis.text(
        -2.90,
        -2.88,
        "Ephemerides are deduplicated, projected through WCS,\nand retained only inside 1≤x≤9216, 1≤y≤9232.",
        fontsize=12.5,
        va="bottom",
        color=STATISTICS_COLORS["mid_grey"],
        bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
    )
    axis.axis("off")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--residuals", type=Path, required=True)
    parser.add_argument("--random-shift", type=Path, required=True)
    parser.add_argument("--output-base", type=Path, required=True)
    parser.add_argument("--figure-data-dir", type=Path, required=True)
    args = parser.parse_args()

    apply_statistics_style()
    residuals = pd.read_parquet(
        args.residuals,
        columns=["dra_cosdec_arcsec", "ddec_arcsec", "radial_residual_arcsec"],
    )
    for column in residuals:
        residuals[column] = pd.to_numeric(residuals[column], errors="coerce")
    residuals = residuals.replace([np.inf, -np.inf], np.nan).dropna()
    shifts = pd.read_csv(args.random_shift)
    args.figure_data_dir.mkdir(parents=True, exist_ok=True)

    x_edges = np.linspace(-1.0, 1.0, 101)
    y_edges = np.linspace(-1.0, 1.0, 101)
    counts, _, _ = np.histogram2d(
        residuals["dra_cosdec_arcsec"], residuals["ddec_arcsec"], bins=[x_edges, y_edges]
    )
    grid = pd.DataFrame(
        {
            "dra_bin_center_arcsec": np.repeat((x_edges[:-1] + x_edges[1:]) / 2, len(y_edges) - 1),
            "ddec_bin_center_arcsec": np.tile((y_edges[:-1] + y_edges[1:]) / 2, len(x_edges) - 1),
            "matched_n": counts.astype(int).ravel(),
        }
    )
    grid.to_csv(args.figure_data_dir / "fig05_signed_residual_density.csv", index=False)

    radial = np.sort(residuals["radial_residual_arcsec"].to_numpy(dtype=float))
    radial_cdf = pd.DataFrame(
        {
            "radial_residual_arcsec": radial,
            "empirical_cdf": np.arange(1, len(radial) + 1) / len(radial),
        }
    )
    radial_cdf.iloc[:: max(1, len(radial_cdf) // 5000)].to_csv(
        args.figure_data_dir / "fig05_radial_residual_cdf.csv", index=False
    )

    categories = [
        ("Unshifted", shifts["test"].eq("unshifted_baseline")),
        ("Random 30″ shift", shifts["test"].eq("random_position_shift")),
    ]
    chance_rows = []
    for label, mask in categories:
        sample = shifts.loc[mask]
        total = int(pd.to_numeric(sample["predicted_n"], errors="coerce").fillna(0).sum())
        success = int(pd.to_numeric(sample["match_lt_1arcsec_n"], errors="coerce").fillna(0).sum())
        low, high = wilson(success, total)
        chance_rows.append(
            {
                "test": label,
                "trial_rows": int(len(sample)),
                "prediction_trials_n": total,
                "matches_lt_1arcsec_n": success,
                "match_fraction": success / total if total else np.nan,
                "ci95_low": low,
                "ci95_high": high,
            }
        )
    chance = pd.DataFrame(chance_rows)
    chance.to_csv(args.figure_data_dir / "fig05_random_shift_summary.csv", index=False)

    figure = plt.figure(figsize=(18, 14))
    geometry = figure.add_subplot(2, 2, 1)
    density = figure.add_subplot(2, 2, 2)
    cdf = figure.add_subplot(2, 2, 3)
    random = figure.add_subplot(2, 2, 4)
    geometry_panel(geometry)

    masked = np.ma.masked_where(counts.T <= 0, counts.T)
    mesh = density.pcolormesh(
        x_edges,
        y_edges,
        masked,
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max(2, float(np.nanmax(counts)))),
        shading="auto",
        rasterized=True,
    )
    density.axhline(0, color="white", linewidth=0.7, alpha=0.65)
    density.axvline(0, color="white", linewidth=0.7, alpha=0.65)
    density.add_patch(Rectangle((-1, -1), 2, 2, fill=False, edgecolor=STATISTICS_COLORS["ink"], linewidth=1.0))
    density.set_xlim(-1, 1)
    density.set_ylim(-1, 1)
    density.set_aspect("equal")
    density.set_xlabel(r"$\Delta\alpha\cos\delta$ (arcsec)")
    density.set_ylabel(r"$\Delta\delta$ (arcsec)")
    density.set_title(f"(b) Signed residuals (N={len(residuals):,})", loc="left", fontsize=20, weight="bold")
    style_statistics_axis(density, density=True, tick_fontsize=18)
    add_statistics_colorbar(density, mesh, "Matched detections per bin")

    cdf.plot(radial, np.arange(1, len(radial) + 1) / len(radial), color=STATISTICS_COLORS["blue"], lw=2.7)
    quantiles = np.percentile(radial, [50, 84, 95])
    for value, percentile, color in zip(quantiles, [50, 84, 95], [STATISTICS_COLORS["red"], STATISTICS_COLORS["orange"], STATISTICS_COLORS["green"]]):
        cdf.axvline(value, color=color, linestyle="--", linewidth=1.9, label=f"p{percentile} = {value:.3f}″")
    cdf.axvline(1.0, color=STATISTICS_COLORS["ink"], linewidth=2.0, label="1″ match truncation")
    cdf.set_xlim(0, 1.02)
    cdf.set_ylim(0, 1.01)
    cdf.set_xlabel("Radial predicted–measured separation (arcsec)")
    cdf.set_ylabel("Empirical cumulative fraction")
    cdf.set_title("(c) Radial separation is threshold-truncated", loc="left", fontsize=20, weight="bold")
    cdf.legend(fontsize=16, loc="lower right")
    style_statistics_axis(cdf, tick_fontsize=18)

    y = np.arange(len(chance))
    fractions = chance["match_fraction"].to_numpy(dtype=float) * 100
    lower = (chance["match_fraction"] - chance["ci95_low"]).to_numpy(dtype=float) * 100
    upper = (chance["ci95_high"] - chance["match_fraction"]).to_numpy(dtype=float) * 100
    random.barh(y, fractions, color=[STATISTICS_COLORS["blue"], STATISTICS_COLORS["orange"]], alpha=0.78)
    random.errorbar(fractions, y, xerr=np.vstack([lower, upper]), fmt="none", ecolor=STATISTICS_COLORS["ink"], capsize=5, lw=1.5)
    random.set_yticks(y, chance["test"], fontsize=18)
    random.invert_yaxis()
    random.set_xlabel("Matches within 1″ (%)")
    random.set_title("(d) Position-shift chance-match control", loc="left", fontsize=20, weight="bold")
    for index, row in chance.iterrows():
        random.text(
            row["match_fraction"] * 100 + max(fractions) * 0.025,
            index,
            f"{int(row['matches_lt_1arcsec_n']):,}/{int(row['prediction_trials_n']):,}\n({row['match_fraction']*100:.3f}%)",
            va="center",
            fontsize=16,
        )
    random.set_xlim(0, max(fractions) * 1.34)
    style_statistics_axis(random, tick_fontsize=18)
    random.text(
        0.02,
        0.04,
        "Control: one sampled exposure per production night;\n8 independent 30″ position shifts per exposure.",
        transform=random.transAxes,
        fontsize=14,
        color=STATISTICS_COLORS["mid_grey"],
        va="bottom",
    )

    figure.text(
        0.5,
        0.012,
        "Production association uses nearest-neighbor separation <1″ and a finite Mag_Kron schema gate; it is not one-to-one matching.",
        ha="center",
        fontsize=16,
        color=STATISTICS_COLORS["mid_grey"],
    )
    figure.tight_layout(rect=[0.02, 0.045, 0.985, 0.99], h_pad=2.0, w_pad=1.7)
    save_pdf_png(figure, args.output_base)


if __name__ == "__main__":
    main()
