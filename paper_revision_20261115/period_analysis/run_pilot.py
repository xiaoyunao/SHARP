"""Run a bounded, offline pilot. Candidate periods are not physical spin claims."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from period_core import search, fit, design, phase_gap, held_night_test, half_cycle_test


def diagnostic(points, candidates, spectrum, output, key):
    t = points.t_mjd.to_numpy() - points.t_mjd.min()
    y, err = points.mag.to_numpy(), points.fit_err.to_numpy()
    night = points.physical_night.to_numpy()
    best = candidates[0][0]
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.7), constrained_layout=True)
    axes[0].errorbar(t, y, yerr=err, fmt=".", ms=3, color="#235b85", alpha=.7)
    axes[0].invert_yaxis()
    axes[0].set(xlabel="UTC midpoint − first midpoint (d)", ylabel="Calibrated aperture magnitude")
    freq, chi2 = spectrum["frequency_cpd"], spectrum["fourier_chi2"]
    finite = np.isfinite(chi2)
    axes[1].plot(freq[finite], (chi2[finite] - min(chi2[finite])) / len(t), lw=.7, color="#235b85")
    for candidate, _, _, _ in candidates[:3]:
        axes[1].axvline(candidate.frequency, color="#b55432", lw=.7, alpha=.6)
    axes[1].set(xlabel="Trial frequency (cycles/day)", ylabel="Grid Δχ² / N")
    offsets = best.coefficients[:len(best.levels)]
    level_index = np.searchsorted(best.levels, night)
    relative = y - offsets[level_index]
    phase = np.mod(t * best.frequency, 1)
    _, colors = np.unique(night, return_inverse=True)
    for shift in (0, 1):
        axes[2].errorbar(phase + shift, relative, yerr=err, fmt="none", ecolor=".7", lw=.5)
        axes[2].scatter(phase + shift, relative, c=colors, cmap="tab20", s=10)
    grid = np.linspace(0, 2, 500)
    x, _ = design(grid / best.frequency, np.zeros(len(grid)), best.frequency)
    axes[2].plot(grid, x[:, -4:] @ best.coefficients[-4:], color="black", lw=1)
    axes[2].invert_yaxis()
    axes[2].set(xlabel=f"Phase at candidate {24 / best.frequency:.5f} h", ylabel="Magnitude − fitted night mean")
    fig.suptitle(f"{key} | {len(t)} points / {len(best.levels)} nights | exploratory; ISO flags unavailable", fontsize=10)
    fig.savefig(output / f"{key.replace(':', '_')}.png", dpi=160)
    fig.savefig(output / f"{key.replace(':', '_')}.pdf")
    plt.close(fig)


def run(inputs, output, limit, plots_limit=10, cv_targets=10, objects=None):
    output.mkdir(parents=True, exist_ok=True)
    catalog = pd.read_csv(inputs / "sampling_catalog.csv.gz")
    selected = catalog.loc[catalog.pilot_eligible]
    if objects:
        selected = selected.loc[selected.object_key.isin(objects)]
        if set(objects) != set(selected.object_key):
            raise ValueError("Requested object missing from eligible window catalog")
    selected = selected.head(limit)
    points = pd.read_csv(inputs / "eligible_window_points.csv.gz", dtype={"physical_night": str})
    all_rows, summaries = [], []
    for target_index, key in enumerate(selected.object_key):
        d = points.loc[points.object_key.eq(key)].sort_values("t_mjd")
        t = d.t_mjd.to_numpy() - d.t_mjd.min()
        y, err, night = d.mag.to_numpy(), d.fit_err.to_numpy(), d.physical_night.to_numpy()
        candidates, spectrum = search(t, y, err, night)
        if not candidates:
            summaries.append(dict(object_key=key, status="rank_failure"))
            continue
        best = candidates[0][0]
        target_rows = []
        for rank, (result, edge, origin, width) in enumerate(candidates, 1):
            phase = np.linspace(0, 1, 1000)
            x, _ = design(phase / result.frequency, np.zeros(len(phase)), result.frequency)
            model = x[:, -4:] @ result.coefficients[-4:]
            shared = fit(t, y, err, night, result.frequency, mode="shared")
            row = dict(object_key=key, rank=rank, candidate_period_h=24 / result.frequency,
                       frequency_cpd=result.frequency, seed_origin=origin, delta_bic=result.bic - best.bic,
                       chi2=result.chi2, dof=len(t) - len(result.coefficients) - 1,
                       fourier_order=2, amplitude_mag=float(np.ptp(model)),
                       rms_mag=float(np.sqrt(np.mean(result.residual ** 2))),
                       max_phase_gap=phase_gap(t, result.frequency),
                       refinement_boundary=edge, night_model_bic=result.bic,
                       shared_model_bic=shared.bic if shared else None,
                       condition=result.condition, n_points=len(t), n_nights=len(np.unique(night)),
                       baseline_days=float(np.ptp(t)), status="exploratory_missing_iso_flags_and_geometry")
            if rank <= 6 and target_index < cv_targets:
                row.update(held_night_test(t, y, err, night, result.frequency, width))
            else:
                row.update(cv_nights=0, cv_conditional_rmse_mag=None, cv_null_rmse_mag=None)
            all_rows.append(row)
            target_rows.append(row)
        near = [r for r in target_rows[1:] if r["delta_bic"] < 10]
        symmetry = half_cycle_test(t, y, err, night, best.frequency / 2)
        summaries.append(dict(object_key=key, candidate_period_h=24 / best.frequency,
                              doubled_candidate_period_h=48 / best.frequency,
                              doubled_period_odd_even_delta_bic=symmetry["half_cycle_delta_bic"],
                              doubled_period_status=symmetry["half_cycle_status"],
                              n_points=len(t), n_nights=len(best.levels),
                              baseline_days=float(np.ptp(t)), competing_delta_bic_lt10=len(near),
                              **{k: target_rows[0][k] for k in ["rms_mag", "amplitude_mag", "cv_nights",
                                  "cv_conditional_rmse_mag", "cv_null_rmse_mag", "max_phase_gap"]},
                              status="exploratory_not_confirmed_not_novelty_checked"))
        if target_index < plots_limit:
            pd.DataFrame(spectrum).to_csv(output / f"{key.replace(':', '_')}_spectrum.csv.gz", index=False,
                                         compression={"method": "gzip", "mtime": 0})
            diagnostic(d, candidates, spectrum, output, key)
        print(key, "points", len(d), "candidate_h", round(24 / best.frequency, 5),
              "close_competitors", len(near), flush=True)
    pd.DataFrame(all_rows).to_csv(output / "period_hypotheses.csv", index=False)
    pd.DataFrame(summaries).to_csv(output / "pilot_summary.csv", index=False)
    points.loc[points.object_key.isin(selected.object_key)].to_csv(output / "pilot_points.csv.gz", index=False,
                                                                 compression={"method": "gzip", "mtime": 0})
    (output / "run_config.json").write_text(json.dumps(dict(
        limit=limit, plots_limit=plots_limit, cv_targets=cv_targets, objects=objects,
        ranking="most points in time-selected <=30d window; no reference periods",
        min_period_h=1, max_period_h=120, grid_samples_per_peak=8, fourier_order=2,
        error_floor_mag=.03, clipping=False, night_intercepts=True,
        cv="top6 basins; train-only local frequency fit; held-night mean fitted conditionally",
        doubled_period_test="order4 vs even harmonics2,4 at 2P on identical points; scaled errors",
        interpretation="all periods are search hypotheses, not confirmed spin or new periods"), indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=10)
    parser.add_argument("--plots-limit", type=int, default=10)
    parser.add_argument("--cv-targets", type=int, default=10)
    parser.add_argument("--objects", nargs="+")
    args = parser.parse_args()
    run(args.inputs, args.output, args.limit, args.plots_limit, args.cv_targets, args.objects)
