"""Numerical checks; optional read-only equivalence check against delivery source."""
import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from period_core import fit, search, half_cycle_test


def validate(delivery=None):
    rng = np.random.default_rng(27)
    t = np.concatenate([n + np.sort(rng.uniform(0, .25, 45)) for n in [0, 1, 3, 4, 6]])
    night = np.repeat(np.arange(5), 45)
    true_f = 3.712
    y = (17 + .16 * np.cos(2 * np.pi * true_f * t) + .09 * np.sin(4 * np.pi * true_f * t)
         + np.repeat([0, .1, -.05, .03, .1], 45) + rng.normal(0, .012, len(t)))
    err = np.full(len(t), .02)
    candidates, _ = search(t, y, err, night)
    measured_f = candidates[0][0].frequency
    assert abs(measured_f - true_f) < .002
    assert any(abs(c[0].frequency - true_f / 2) < .02 for c in candidates)
    assert any(abs(c[0].frequency - (true_f - 1)) < .04 for c in candidates)
    # One observation per night cannot constrain a model with free nightly means.
    assert fit(t[:12], y[:12], err[:12], np.arange(12), true_f) is None
    doubled = half_cycle_test(t, y, err, night, measured_f / 2)
    assert doubled["half_cycle_status"] == "half_double_ambiguous"
    report = dict(synthetic_seed=27, injected_frequency_cpd=true_f, recovered_frequency_cpd=measured_f,
                  relative_period_error=abs(true_f / measured_f - 1),
                  half_frequency_and_daily_alias_retained=True, unconstrained_night_model_rejected=True,
                  doubled_period_capacity_test=doubled)
    if delivery:
        code = delivery / "sitian_period_analysis_code"
        checked = []
        for line in (code / "SHA256SUMS.txt").read_text().splitlines():
            if not line.strip():
                continue
            expected, name = line.split(None, 1)
            actual = hashlib.sha256((code / name.strip()).read_bytes()).hexdigest()
            assert actual == expected, name
            checked.append(dict(path=name.strip(), sha256=actual))
        source = (code / "upstream/peroid_ana_update.py").read_text()
        # Execute only the small, reviewed pure numerical definition, not the
        # upstream module (which contains server paths and publication code).
        outer = next(n for n in ast.parse(source).body
                     if isinstance(n, ast.FunctionDef) and n.name == "phase_fold_and_plot_P_only")
        node = next(n for n in outer.body if isinstance(n, ast.FunctionDef) and n.name == "fourier_fit_fixed")
        namespace = {"np": np}
        exec(compile(ast.Module(body=[node], type_ignores=[]), "<reviewed_fourier_fit_fixed>", "exec"), namespace)
        results = pd.read_csv(delivery / "sitian_verified_full_results.csv")
        maximum_error, count, points = 0., 0, 0
        for row in results.itertuples():
            d = pd.read_csv(delivery / "sitian_verified_folded_points/points" / f"{row.number}.csv")
            d = d.loc[d.in_fit.eq(True)]
            phase = d.phase.to_numpy()
            mag, error = d.mag.to_numpy(), d.mag_err.to_numpy()
            order = int(row.verified_fit_order)
            old, _ = namespace["fourier_fit_fixed"](phase, mag, error, order)
            new = fit(phase, mag, error, np.zeros(len(d)), 1., order, "shared")
            assert new is not None
            diff = float(np.max(np.abs(old - new.coefficients)))
            assert diff < 1e-9
            maximum_error = max(maximum_error, diff)
            count += 1
            points += len(d)
        report.update(delivery_files_sha256=checked, delivery_fixed_sample_targets=count,
                      delivery_fixed_sample_points=points, max_coefficient_difference=maximum_error,
                      equivalence_scope="weighted fixed-period Fourier solver only; search/clipping intentionally changed")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--delivery", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = validate(args.delivery)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
