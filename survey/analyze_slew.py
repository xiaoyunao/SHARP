#!/usr/bin/env python3
"""Fit effective two-axis survey transition times from frozen raw headers.

This is an offline calibration tool, not a scheduler configuration writer.
Angles are coordinate axis angles (RA hours * 15, without cos(Dec)).
All errors and fits are in seconds after subtracting exposure and readout.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


def wrap(angle):
    return (angle + 180.0) % 360.0 - 180.0


def ramp_time(angle, velocity, acceleration):
    """Symmetric acceleration/deceleration, with an optional cruise segment."""
    return np.where(angle <= velocity**2 / acceleration,
                    2 * np.sqrt(angle / acceleration),
                    angle / velocity + velocity / acceleration)


MODELS = {
    "constant": ([35.], [0.], [300.], ["constant_s"]),
    "parallel_linear": ([20., 2., 2.], [0., 0., 0.], [60., 100., 100.],
                        ["overhead_s", "ra_s_per_deg", "dec_s_per_deg"]),
    "sequential_linear": ([20., 2., 2.], [0., 0., 0.], [60., 100., 100.],
                          ["overhead_s", "ra_s_per_deg", "dec_s_per_deg"]),
    "parallel_sqrt": ([9., 12., 12.], [0., 0., 0.], [60., 100., 100.],
                      ["overhead_s", "ra_s_per_sqrt_deg", "dec_s_per_sqrt_deg"]),
    "sequential_sqrt": ([9., 3., 12.], [0., 0., 0.], [60., 100., 100.],
                        ["overhead_s", "ra_s_per_sqrt_deg", "dec_s_per_sqrt_deg"]),
    "parallel_power": ([10., 10., .55, 10., .55], [0., .001, .2, .001, .2],
                       [60., 100., 1.5, 100., 1.5],
                       ["overhead_s", "ra_coefficient", "ra_exponent",
                        "dec_coefficient", "dec_exponent"]),
    "parallel_trapezoid": ([9., 1., .027, .7, .027],
                           [0., .03, .0001, .03, .0001], [60., 20., 20., 20., 20.],
                           ["overhead_s", "ra_vmax_deg_s", "ra_accel_deg_s2",
                            "dec_vmax_deg_s", "dec_accel_deg_s2"]),
}


def predict(model, p, x, y):
    if model == "constant":
        return np.full_like(x, p[0], dtype=float)
    if model == "parallel_linear":
        return p[0] + np.maximum(p[1] * x, p[2] * y)
    if model == "sequential_linear":
        return p[0] + p[1] * x + p[2] * y
    if model == "parallel_sqrt":
        return p[0] + np.maximum(p[1] * np.sqrt(x), p[2] * np.sqrt(y))
    if model == "sequential_sqrt":
        return p[0] + p[1] * np.sqrt(x) + p[2] * np.sqrt(y)
    if model == "parallel_power":
        return p[0] + np.maximum(p[1] * x**p[2], p[3] * y**p[4])
    if model == "parallel_trapezoid":
        return p[0] + np.maximum(ramp_time(x, p[1], p[2]),
                                 ramp_time(y, p[3], p[4]))
    raise ValueError(model)


def fit(model, x, y, time, start=None):
    initial, lower, upper, _ = MODELS[model]
    return least_squares(lambda p: predict(model, p, x, y) - time,
                         initial if start is None else start,
                         bounds=(lower, upper), loss="soft_l1", f_scale=1.,
                         x_scale="jac", max_nfev=2000)


def error_metrics(observed, predicted):
    residual = observed - predicted
    return {"n": int(len(residual)), "mae_s": float(np.mean(abs(residual))),
            "median_abs_error_s": float(np.median(abs(residual))),
            "p90_abs_error_s": float(np.quantile(abs(residual), .9)),
            "rmse_s": float(np.sqrt(np.mean(residual**2))),
            "median_residual_s": float(np.median(residual)),
            "abs_error_gt5_s": int(np.sum(abs(residual) > 5))}


def make_pairs(snapshots, readout):
    frames, metadata = [], []
    for path in snapshots:
        source = json.loads(path.read_text())
        metadata.append({key: source[key] for key in
                         ("night", "directory", "retrieved_utc", "method")})
        metadata[-1]["snapshot_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        for row in source["rows"]:
            frames.append({**{k: v for k, v in row.items() if k not in ("header", "comments")},
                           **row.get("header", {}), "night": source["night"]})
    frames = pd.DataFrame(frames)
    frames["timestamp"] = pd.to_datetime(frames["OBS_DATE"], errors="coerce")
    frames["date_timestamp"] = pd.to_datetime(frames["DATE"], errors="coerce")
    # Sort all OBJ exposures before selecting science pairs: never bridge a test frame.
    frames = frames.sort_values(["night", "timestamp", "name"]).reset_index(drop=True)
    for key in ["EXPOSURE", "OBJ_RA", "OBJ_DEC", "TEL_RA", "TEL_DEC", "HOURANG", "OBLST"]:
        frames[key] = pd.to_numeric(frames[key], errors="coerce")
    previous = frames.groupby("night", sort=False).shift(1)
    pairs = frames[["night", "name", "timestamp", "sequence", "TEL_EL"]].copy()
    pairs["previous_name"] = previous["name"]
    pairs["sequence_step"] = frames.sequence - previous.sequence
    pairs["dt_obs_s"] = (frames.timestamp - previous.timestamp).dt.total_seconds()
    pairs["dt_date_s"] = (frames.date_timestamp - previous.date_timestamp).dt.total_seconds()
    pairs["previous_exposure_s"] = previous.EXPOSURE
    pairs["time_s"] = pairs.dt_obs_s - previous.EXPOSURE - readout
    pairs["time_date_s"] = pairs.dt_date_s - previous.EXPOSURE - readout
    pairs["signed_ra_deg"] = wrap(15 * (frames.OBJ_RA - previous.OBJ_RA))
    pairs["signed_dec_deg"] = frames.OBJ_DEC - previous.OBJ_DEC
    pairs["ra_deg"] = abs(pairs.signed_ra_deg)
    pairs["dec_deg"] = abs(pairs.signed_dec_deg)
    pairs["tel_ra_deg"] = abs(wrap(15 * (frames.TEL_RA - previous.TEL_RA)))
    pairs["tel_dec_deg"] = abs(frames.TEL_DEC - previous.TEL_DEC)
    pairs["hourangle_step_deg"] = wrap(15 * (frames.HOURANG - previous.HOURANG))
    pairs["lst_step_deg"] = wrap(15 * (frames.OBLST - previous.OBLST))
    pairs["offset_obs_minus_date_s"] = (frames.timestamp - frames.date_timestamp).dt.total_seconds()
    # Use a disjoint exclusion label, preserving all rows for audit.
    reason = np.full(len(pairs), "eligible", dtype=object)
    reason[~np.isfinite(pairs.time_s)] = "no_valid_previous_time"
    valid = reason == "eligible"
    science = frames.standard_mp.eq(True) & previous.standard_mp.eq(True)
    reason[valid & ~science] = "nonstandard_exposure_boundary"
    valid = reason == "eligible"
    reason[valid & (pairs.sequence_step != 1)] = "sequence_not_consecutive"
    valid = reason == "eligible"
    finite = np.isfinite(pairs[["ra_deg", "dec_deg", "previous_exposure_s"]]).all(axis=1)
    reason[valid & ~finite] = "invalid_coordinate_or_exposure"
    valid = reason == "eligible"
    reason[valid & ((pairs.time_s <= 0) | (pairs.time_s > 180))] = "time_outside_0_180_s"
    pairs["selection"] = reason
    return frames, pairs, metadata


def cross_validate(model, x, y, time, folds):
    predicted = np.empty(len(time))
    for fold in np.unique(folds):
        train, test = folds != fold, folds == fold
        result = fit(model, x[train], y[train], time[train])
        predicted[test] = predict(model, result.x, x[test], y[test])
    return predicted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshots", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--readout", type=float, default=14.47)
    parser.add_argument("--bootstrap", type=int, default=300)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    frames, pairs, metadata = make_pairs(args.snapshots, args.readout)
    eligible = pairs[pairs.selection == "eligible"].copy().reset_index(drop=True)
    x, y, time = [eligible[key].to_numpy(float) for key in ("ra_deg", "dec_deg", "time_s")]
    # Keep repeated field transitions together in the geometry validation.
    transition = eligible.previous_name.str.extract(r"OBJ_MP_(\d+)_", expand=False) + "_" + eligible.name.str.extract(r"OBJ_MP_(\d+)_", expand=False)
    rng = np.random.default_rng(20261022)
    unique = np.unique(transition)
    shuffled = rng.permutation(unique)
    mapping = {key: i % 5 for i, key in enumerate(shuffled)}
    geometry_folds = transition.map(mapping).to_numpy()
    time_folds = np.minimum(np.arange(len(time)) * 5 // len(time), 4)
    eligible["geometry_fold"] = geometry_folds
    eligible["time_fold"] = time_folds
    results, fitted = {}, {}
    for model in MODELS:
        result = fit(model, x, y, time)
        # Multiple deterministic starts for the more flexible motion model.
        if model == "parallel_trapezoid":
            for vr, vd in ((.5, .5), (2., 2.), (5., 5.)):
                candidate = fit(model, x, y, time, [9., vr, .027, vd, .027])
                if candidate.cost < result.cost:
                    result = candidate
        fitted[model] = result.x
        pred = predict(model, result.x, x, y)
        geo_pred = cross_validate(model, x, y, time, geometry_folds)
        temporal_pred = cross_validate(model, x, y, time, time_folds)
        eligible[model + "_pred_s"] = pred
        eligible[model + "_geometry_oof_s"] = geo_pred
        eligible[model + "_temporal_oof_s"] = temporal_pred
        results[model] = {"parameters": dict(zip(MODELS[model][3], result.x.tolist())),
                          "robust_cost": float(result.cost),
                          "in_sample": error_metrics(time, pred),
                          "geometry_grouped_cv": error_metrics(time, geo_pred),
                          "contiguous_time_cv": error_metrics(time, temporal_pred)}
        print(model, json.dumps(results[model]), flush=True)
    chosen = "parallel_sqrt"
    p = fitted[chosen]
    ra_component, dec_component = p[1] * np.sqrt(x), p[2] * np.sqrt(y)
    eligible["dominant_axis"] = np.where(ra_component > dec_component, "RA", "Dec")
    eligible["axis_margin_s"] = abs(ra_component - dec_component)
    eligible["residual_s"] = time - predict(chosen, p, x, y)
    eligible["anomaly_abs_residual_gt5_s"] = abs(eligible.residual_s) > 5
    # Bootstrapping blocks of 15 adjacent transitions retains short-time dependence.
    blocks = np.arange(len(time)) // 15
    unique_blocks = np.unique(blocks)
    bootstrap = []
    for _ in range(args.bootstrap):
        sample = rng.choice(unique_blocks, len(unique_blocks), replace=True)
        indices = np.concatenate([np.flatnonzero(blocks == b) for b in sample])
        b = fit(chosen, x[indices], y[indices], time[indices], p).x
        bootstrap.append(b)
    bootstrap = np.asarray(bootstrap)
    np.savetxt(out / "bootstrap_parameters.csv", bootstrap, delimiter=",",
               header=",".join(MODELS[chosen][3]), comments="")
    intervals = {name: {"value": float(p[i]), "block_bootstrap_95pct":
                        np.quantile(bootstrap[:, i], [.025, .975]).tolist()}
                 for i, name in enumerate(MODELS[chosen][3])}
    sensitivity = {}
    for label, sx, sy, st in (
        ("DATE_integer_seconds", x, y, eligible.time_date_s.to_numpy()),
        ("TEL_coordinates", eligible.tel_ra_deg.to_numpy(), eligible.tel_dec_deg.to_numpy(), time),
        ("readout_plus_1_second", x, y, time - 1),
        ("readout_minus_1_second", x, y, time + 1),
    ):
        estimate = fit(chosen, sx, sy, st).x
        sensitivity[label] = {"parameters": estimate.tolist(),
                              "in_sample": error_metrics(st, predict(chosen, estimate, sx, sy))}
    # Profile cruise speed while re-fitting all remaining parameters. A flat
    # high-speed profile explicitly demonstrates absence of an upper constraint.
    profile = []
    for axis, parameter_index in (("RA", 1), ("Dec", 3)):
        grid = np.geomspace(.2, 10., 60)
        for velocity in grid:
            free = [i for i in range(5) if i != parameter_index]
            initial, lower, upper, _ = MODELS["parallel_trapezoid"]
            def full(q):
                v = fitted["parallel_trapezoid"].copy()
                v[parameter_index], v[free] = velocity, q
                return v
            result = least_squares(
                lambda q: predict("parallel_trapezoid", full(q), x, y) - time,
                fitted["parallel_trapezoid"][free],
                bounds=(np.array(lower)[free], np.array(upper)[free]),
                loss="soft_l1", f_scale=1., x_scale="jac", max_nfev=1000)
            profile.append({"axis": axis, "vmax_deg_s": velocity,
                            "robust_cost": result.cost})
    pd.DataFrame(profile).to_csv(out / "vmax_profile.csv", index=False)
    audit = {
        "source": metadata, "readout_s": args.readout,
        "frame_count": len(frames), "standard_mp_frames": int(frames.standard_mp.sum()),
        "header_error_count": int(frames.get("error", pd.Series(dtype=str)).notna().sum()),
        "pair_selection_counts": pairs.selection.value_counts().to_dict(),
        "date_equals_obs_floor_second": bool((frames.timestamp.dt.floor("s") == frames.date_timestamp).all()),
        "timestamp_difference_disagreement_std_s": float((eligible.dt_obs_s - eligible.dt_date_s).std()),
        "unique_directed_field_transitions": len(unique),
        "ra_range_deg": [float(min(x)), float(max(x))],
        "dec_range_deg": [float(min(y)), float(max(y))],
        "dominant_axis_counts": eligible.dominant_axis.value_counts().to_dict(),
        "dominant_angle_ranges_deg": {
            axis: [float(eligible.loc[eligible.dominant_axis == axis, column].min()),
                   float(eligible.loc[eligible.dominant_axis == axis, column].max())]
            for axis, column in (("RA", "ra_deg"), ("Dec", "dec_deg"))},
        "confident_axis_counts_margin_gt2s": eligible.loc[eligible.axis_margin_s > 2, "dominant_axis"].value_counts().to_dict(),
        "anomalies": eligible.loc[eligible.anomaly_abs_residual_gt5_s,
                    ["night", "previous_name", "name", "time_s", "ra_deg", "dec_deg", "residual_s"]].to_dict("records"),
        "models": results, "chosen_model": chosen, "parameter_intervals": intervals,
        "sensitivity": sensitivity,
        "interpretation": [
            "OBS_DATE differences cancel any constant timestamp offset; variable delay does not cancel.",
            "No fit residual is used to exclude an eligible pair. Robust loss includes all pairs.",
            "The chosen square-root model is a parsimonious effective interpolation over this sample.",
            "The slightly more accurate free-exponent model does not materially improve contiguous-time validation.",
            "Its coefficients have a kinematic acceleration interpretation only under the symmetric-ramp assumption.",
            "Cruise velocity is not measured if the high-speed profile is flat; do not report a fitted plateau as vmax.",
            "Block bootstrap describes within-night statistical variation, not cross-night or hardware systematics.",
            "Coordinate changes describe motion relative to tracking, not calibrated encoder trajectories.",
        ],
    }
    (out / "results.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False, allow_nan=False))
    pairs.to_csv(out / "all_pair_audit.csv", index=False)
    eligible.to_csv(out / "fitted_pairs.csv", index=False)
    frames.drop(columns=["timestamp", "date_timestamp"]).to_csv(out / "frames.csv", index=False)
    draw_figures(eligible, fitted, intervals, pd.DataFrame(profile), out)
    print("FINAL", json.dumps({"intervals": intervals, "axes": audit["dominant_axis_counts"],
                               "anomalies": len(audit["anomalies"])}), flush=True)


def draw_figures(data, fitted, intervals, profile, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": .18})
    blue, orange = "#246A9E", "#CA722C"
    p = fitted["parallel_sqrt"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    for axis, col, coef, ax in (("RA", "ra_deg", p[1], axes[0, 0]),
                                ("Dec", "dec_deg", p[2], axes[0, 1])):
        ax.scatter(data[col], data.time_s, color="#AAAEB3", s=17, alpha=.5,
                   label="All eligible pairs")
        use = (data.dominant_axis == axis) & (data.axis_margin_s > 2)
        ax.scatter(data.loc[use, col], data.loc[use, "time_s"], s=27,
                   color=blue if axis == "RA" else orange,
                   marker="o" if axis == "RA" else "^", label=f"{axis} slower by >2 s (model)")
        limits = (float(data.loc[use, col].min()), float(data.loc[use, col].max()))
        angle = np.linspace(limits[0], limits[1], 300)
        ax.plot(angle, p[0] + coef * np.sqrt(angle), color="#252A30", lw=1.8,
                label=f"{p[0]:.2f} + {coef:.2f} sqrt(angle)")
        ax.set(xlabel=f"|Delta {axis}| [deg]", ylabel="Gap minus exposure and readout [s]",
               title=f"{axis} branch: inferred from the joint fit")
        ax.legend(fontsize=9, loc="upper left")
    ax = axes[1, 0]
    for model, color, label in (("parallel_linear", orange, "Constant axis speed"),
                                 ("parallel_sqrt", blue, "Square-root motion")):
        pred = data[model + "_geometry_oof_s"]
        ax.scatter(pred, data.time_s, s=16, color=color, alpha=.55, label=label)
    ax.plot([0, 175], [0, 175], "--", color="#252A30", lw=1)
    ax.set(xlabel="Held-out prediction [s]", ylabel="Observed transition time [s]",
           title="5-fold validation: repeated field pairs stay together", xlim=(10, 175), ylim=(10, 175))
    ax.legend(fontsize=9)
    ax = axes[1, 1]
    elapsed = (pd.to_datetime(data.timestamp) - pd.to_datetime(data.timestamp.iloc[0])).dt.total_seconds() / 3600
    ax.scatter(elapsed, data.time_s - data.parallel_sqrt_geometry_oof_s, s=19, color=blue)
    ax.axhline(0, color="#252A30", lw=1)
    ax.set(xlabel="Hours since first eligible pair", ylabel="Observed - held-out prediction [s]",
           title="Residuals across the night (all eligible pairs)")
    fig.suptitle(f"SMT survey slew calibration | {data.night.min()} | {len(data)} consecutive science pairs", fontsize=15)
    fig.savefig(out / "slew_fit.png", dpi=180)
    fig.savefig(out / "slew_fit.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), layout="constrained")
    angle = np.geomspace(.1, 70, 300)
    for axis, coef, color, style in (("RA", p[1], blue, "-"), ("Dec", p[2], orange, "--")):
        use = (data.dominant_axis == axis) & (data.axis_margin_s > 2)
        col = "ra_deg" if axis == "RA" else "dec_deg"
        lo, hi = data.loc[use, col].min(), data.loc[use, col].max()
        # Draw only the observed confidently dominant angle range.
        a = angle[(angle >= lo) & (angle <= hi)]
        axes[0].plot(a, np.sqrt(a) / coef, color=color, ls=style, lw=2,
                     label=f"{axis}: motion only")
        axes[0].plot(a, a / (p[0] + coef * np.sqrt(a)), color=color, ls=":", lw=1.7,
                     label=f"{axis}: including fixed overhead")
        sub = profile[profile.axis == axis]
        axes[1].plot(sub.vmax_deg_s, sub.robust_cost - sub.robust_cost.min(),
                     color=color, ls=style, label=axis)
    axes[0].set(xscale="log", xlabel="Axis move angle [deg]", ylabel="Effective average speed [deg/s]",
                title="Speed varies with move angle")
    axes[0].legend(fontsize=8)
    axes[1].set(xscale="log", xlabel="Assumed maximum axis speed [deg/s]",
                ylabel="Increase in robust fitting objective", title="Maximum-speed identifiability")
    axes[1].set_ylim(0, 30)
    axes[1].legend()
    residual = data.time_s - data.parallel_sqrt_geometry_oof_s
    axes[2].hist(residual, bins=np.arange(-4, 4.001, .25), color=blue, alpha=.8)
    outside = int(np.sum(abs(residual) > 4))
    axes[2].set(xlabel="Held-out residual [s]", ylabel="Number of transitions",
                title=f"Residual core; {outside} pairs outside +/-4 s")
    fig.suptitle("One-night limits: effective timing is measurable; cruise speed may not be", fontsize=14)
    fig.savefig(out / "slew_diagnostics.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
