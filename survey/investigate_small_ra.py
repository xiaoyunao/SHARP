#!/usr/bin/env python3
"""Audit archived small-RA transitions from frozen inventory and header probes.

Keep direct small-Dec pairs, conditional survey delays, and ancillary programs
separate. A long exposure-to-exposure gap is not an axis fault diagnosis.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from . import analyze_slew as slew
except ImportError:
    import analyze_slew as slew


def pair_values(previous, current):
    a, b = previous, current
    ra = slew.wrap(15 * (float(b["OBJ_RA"]) - float(a["OBJ_RA"])))
    dec = float(b["OBJ_DEC"]) - float(a["OBJ_DEC"])
    gap = (pd.Timestamp(b["OBS_DATE"]) - pd.Timestamp(a["OBS_DATE"])).total_seconds()
    return {"signed_ra_deg": float(ra), "ra_deg": abs(float(ra)),
            "dec_deg": abs(dec), "time_s": gap - float(a["EXPOSURE"]) - 14.47,
            "previous_exposure_s": float(a["EXPOSURE"]), "current_exposure_s": float(b["EXPOSURE"])}


def ancillary_pairs(path):
    records = [json.loads(line) for line in path.read_text().splitlines() if line]
    lookup = {(r["night"], r["name"]): r for r in records}
    pairs = []
    for row in records:
        previous = lookup.get((row["night"], row["inventory_previous_name"]))
        if previous is None or "header" not in previous or "header" not in row:
            continue
        if not (row["inventory_sequence_unique"] and previous["inventory_sequence_unique"]):
            continue
        if row["sequence"] - previous["sequence"] != 1:
            continue
        names = previous["name"] + row["name"]
        if not (row["name"].startswith("OBJ_") and previous["name"].startswith("OBJ_")):
            continue
        if any(key in names.lower() for key in ("test", "stellar", "start")):
            continue
        try:
            values = pair_values(previous["header"], row["header"])
        except (ValueError, TypeError, KeyError):
            continue
        if not all(np.isfinite(v) for v in values.values()):
            continue
        pairs.append({"night": row["night"], "previous_name": previous["name"],
                      "name": row["name"], **values})
    return pd.DataFrame(pairs), len(records), sum("error" in r for r in records)


def archive_pairs(directory):
    rows = [json.loads(line) for line in (directory / "archive_mp_headers.jsonl").read_text().splitlines() if line]
    expected = json.loads((directory / "archive_collection.json").read_text())["expected_mp_headers"]
    if len(rows) != expected:
        raise ValueError(f"Archive collection incomplete: {len(rows)}/{expected}")
    lookup = {(r["night"], r["name"]): r for r in rows}
    if len(lookup) != len(rows):
        raise ValueError("Duplicate archive frame identities")
    footprint = json.loads((directory / "server_footprints.json").read_text())
    grid = {str(int(r["field_id"])): r for r in footprint["rows"]}
    mismatch = []
    valid, errors = 0, 0
    for row in rows:
        if "header" not in row:
            errors += 1
            continue
        valid += 1
        field = str(int(row["name"].split("_")[2]))
        target = grid[field]
        header = row["header"]
        try:
            dx = abs(float(slew.wrap(15 * float(header["OBJ_RA"]) - target["ra_deg"])))
            dy = abs(float(header["OBJ_DEC"]) - target["dec_deg"])
            if dx > .001 or dy > .001:
                mismatch.append({"night": row["night"], "name": row["name"], "ra_offset_deg": dx, "dec_offset_deg": dy})
        except (ValueError, KeyError, TypeError):
            mismatch.append({"night": row["night"], "name": row["name"], "invalid_coordinate": True})
    mismatch_keys = {(r["night"], r["name"]) for r in mismatch}
    exposure_groups = {}
    for row in rows:
        if "header" not in row:
            continue
        header = row["header"]
        key = tuple(header.get(k) for k in ("OBS_DATE", "OBJ_RA", "OBJ_DEC", "EXPOSURE", "CHANNEL"))
        exposure_groups.setdefault(key, []).append(row)
    aliases, duplicate_rows = set(), []
    for group in exposure_groups.values():
        if len(group) < 2:
            continue
        def rank(row):
            date = str(row["header"].get("OBS_DATE", ""))[:10].replace("-", "")
            return row["night"] != date, row["night"], row["name"]
        canonical = min(group, key=rank)
        for row in group:
            if row is canonical:
                continue
            aliases.add((row["night"], row["name"]))
            duplicate_rows.append({"night": row["night"], "name": row["name"],
                                   "canonical_night": canonical["night"], "canonical_name": canonical["name"],
                                   "header_hash_identical": row.get("header_sha256") == canonical.get("header_sha256")})
    pd.DataFrame(duplicate_rows).to_csv(directory / "duplicate_exposure_aliases.csv", index=False)
    pairs = []
    for row in rows:
        previous = lookup.get((row["night"], row["inventory_previous_name"]))
        if previous is None or "header" not in row or "header" not in previous:
            continue
        if not (row["inventory_sequence_unique"] and previous["inventory_sequence_unique"]):
            continue
        if row["sequence"] - previous["sequence"] != 1:
            continue
        try:
            values = pair_values(previous["header"], row["header"])
        except (ValueError, TypeError, KeyError):
            continue
        if not all(np.isfinite(v) for v in values.values()):
            continue
        pairs.append({"night": row["night"], "previous_name": previous["name"],
                      "name": row["name"], **values,
                      "both_grid_consistent": ((row["night"], row["name"]) not in mismatch_keys and
                                                (previous["night"], previous["name"]) not in mismatch_keys),
                      "both_unique_exposure": ((row["night"], row["name"]) not in aliases and
                                                (previous["night"], previous["name"]) not in aliases)})
    pairs = pd.DataFrame(pairs)
    pairs.to_csv(directory / "archive_adjacent_mp_pairs.csv", index=False)
    pd.DataFrame(mismatch).to_csv(directory / "archive_grid_mismatches.csv", index=False)
    zero = pairs[(pairs.ra_deg <= 1e-5) & (pairs.dec_deg <= 1e-5) & (pairs.time_s > 0) & pairs.both_unique_exposure]
    return pairs, {"mp_frames": len(rows), "readable_headers": valid, "header_errors": errors,
                   "mp_nights": len({row["night"] for row in rows}),
                   "duplicate_exposure_aliases": len(aliases), "unique_readable_exposure_keys": len(exposure_groups),
                   "pairs_touching_duplicate_alias": int((~pairs.both_unique_exposure).sum()),
                   "adjacent_mp_pairs": len(pairs), "grid_mismatch_frames": len(mismatch),
                   "pairs_touching_grid_mismatch": int((~pairs.both_grid_consistent).sum()),
                   "zero_recorded_pointing_change": {"pairs": len(zero),
                       "median_time_s": float(zero.time_s.median()) if len(zero) else None,
                       "time_gt30s": int((zero.time_s > 30).sum()),
                       "interpretation": "Zero recorded target displacement does not establish whether a go-to correction was reissued."}}


def matched_comparison(events, small_mask, reference_mask, output):
    rows = []
    for (night, dec_bin), group in events.groupby(["night", "dec_bin"]):
        small = group[small_mask.loc[group.index]]
        reference = group[reference_mask.loc[group.index]]
        if len(small) < 5 or len(reference) < 5:
            continue
        rows.append({"night": str(night), "dec_bin": int(dec_bin),
                     "small_n": len(small), "reference_n": len(reference),
                     "small_rate": float(small.extra_gt5s.mean()),
                     "reference_rate": float(reference.extra_gt5s.mean()),
                     # The common nightly offset cancels from this difference.
                     "median_adjusted_time_difference_s": float(small.extra_time_s.median() - reference.extra_time_s.median())})
    matched = pd.DataFrame(rows)
    matched.to_csv(output, index=False)
    if not len(matched):
        return None
    aggregates = matched.assign(
        risk_numerator=(matched.small_rate - matched.reference_rate) * matched.small_n,
        time_numerator=matched.median_adjusted_time_difference_s * matched.small_n,
    ).groupby("night").agg(risk_numerator=("risk_numerator", "sum"),
                           time_numerator=("time_numerator", "sum"), denominator=("small_n", "sum"))
    rng = np.random.default_rng(20261022)
    sample = rng.integers(0, len(aggregates), size=(500, len(aggregates)))
    denominator = aggregates.denominator.to_numpy()[sample].sum(axis=1)
    risks = aggregates.risk_numerator.to_numpy()[sample].sum(axis=1) / denominator
    times = aggregates.time_numerator.to_numpy()[sample].sum(axis=1) / denominator
    any_events = bool((matched.small_rate > 0).any() or (matched.reference_rate > 0).any())
    return {"strata": len(matched), "nights": len(aggregates),
            "small_pairs": int(matched.small_n.sum()), "reference_pairs": int(matched.reference_n.sum()),
            "standardized_risk_difference": float(aggregates.risk_numerator.sum() / aggregates.denominator.sum()),
            "risk_night_bootstrap_95pct": np.quantile(risks, [.025, .975]).tolist() if any_events else None,
            "mean_matched_median_time_difference_s": float(aggregates.time_numerator.sum() / aggregates.denominator.sum()),
            "time_night_bootstrap_95pct_s": np.quantile(times, [.025, .975]).tolist(),
            "zero_event_caveat": "If both groups have zero observed events, an empirical bootstrap supplies no risk confidence bound.",
            "interpretation": "Association of total waiting after night/Dec matching, not an RA fault rate. Residual hour, target and control-mode confounding remains."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--recent-start", default="20260829")
    parser.add_argument("--use-archive", action="store_true",
                        help="Audit conditional delays over all collected MP headers; retain recent calibration separately")
    args = parser.parse_args()
    source = args.source
    out = source / "expanded"
    inventory = json.loads((out / "archive_inventory.json").read_text())
    geometry = json.loads((out / "geometry_inventory.json").read_text())
    direct = json.loads((out / "isolated_headers.json").read_text())
    direct_rows = []
    for row in direct["rows"]:
        values = pair_values(row["previous"]["header"], row["current"]["header"])
        assert abs(values["ra_deg"] - row["ra_deg_from_grid"]) < .0001
        assert abs(values["dec_deg"] - row["dec_deg_from_grid"]) < .0001
        direct_rows.append({"night": row["night"], "previous_name": row["previous_name"],
                            "name": row["name"], **values})
    direct_pairs = pd.DataFrame(direct_rows)
    full_archive, archive_summary = None, None
    if args.use_archive:
        full_archive, archive_summary = archive_pairs(out)
        all_direct = full_archive[(full_archive.ra_deg > 1e-5) & (full_archive.ra_deg < 1) &
                                  (full_archive.dec_deg < 2) & (full_archive.time_s > 0)].copy()
        all_direct.to_csv(out / "archive_small_ra_all_exposures.csv", index=False)
        direct_pairs = all_direct[all_direct.previous_exposure_s.between(29.9, 30.1) &
                                  all_direct.current_exposure_s.between(29.9, 30.1) &
                                  all_direct.both_grid_consistent & all_direct.both_unique_exposure].copy()
    direct_pairs.to_csv(out / "isolated_small_ra_pairs.csv", index=False)
    paths = [p for p in sorted(source.glob("headers_*.json"))
             if args.recent_start <= p.stem.removeprefix("headers_") <= "20261021"]
    frames, pairs, metadata = slew.make_pairs(paths, 14.47)
    moving = (pairs.ra_deg > 1e-5) | (pairs.dec_deg > 1e-5)
    science_exposure = pairs.previous_exposure_s.between(29.9, 30.1) & pairs.current_exposure_s.between(29.9, 30.1)
    usable = pairs.selection.eq("eligible") & moving & science_exposure
    fit_data = pairs[usable].copy().reset_index(drop=True)
    x, y, t = [fit_data[c].to_numpy(float) for c in ("ra_deg", "dec_deg", "time_s")]
    p = slew.fit("parallel_sqrt", x, y, t).x
    fit_data["dominant_axis"] = np.where(p[1] * np.sqrt(x) > p[2] * np.sqrt(y), "RA", "Dec")
    fit_data["axis_margin_s"] = abs(p[1] * np.sqrt(x) - p[2] * np.sqrt(y))
    prediction = np.empty_like(t)
    nightly = []
    for night, indices in fit_data.groupby("night").groups.items():
        test = np.zeros(len(t), dtype=bool)
        test[np.asarray(indices)] = True
        estimate = slew.fit("parallel_sqrt", x[~test], y[~test], t[~test]).x
        prediction[test] = slew.predict("parallel_sqrt", estimate, x[test], y[test])
        nightly.append({"night": str(night), **slew.error_metrics(t[test], prediction[test])})
    fit_data["night_oof_s"] = prediction
    fit_data["residual_s"] = t - slew.predict("parallel_sqrt", p, x, y)
    fit_data.to_csv(out / "recent_fitted_pairs.csv", index=False)
    slew.draw_axis_fits_only(fit_data, {"parallel_sqrt": p}, out)

    # Include long positive intervals in the delay audit. A pre-fit 180 s cap
    # must not remove exactly the events sought by this investigation.
    selected = pairs.selection.isin(["eligible", "time_outside_0_180_s"]) & moving & (pairs.time_s > 0)
    events = (full_archive[(full_archive.time_s > 0) & full_archive.both_grid_consistent & full_archive.both_unique_exposure].copy()
              if args.use_archive else pairs[selected].copy())
    events = events[(events.ra_deg > 1e-5) & (events.ra_deg < 3) &
                    (events.dec_deg >= 2) & (events.dec_deg <= 20) &
                    (events.ra_deg < .7 * events.dec_deg) &
                    events.previous_exposure_s.between(29.9, 30.1) &
                    events.current_exposure_s.between(29.9, 30.1)]
    control_offsets = {}
    for night, group in events.groupby("night"):
        controls = group[group.ra_deg.between(.3, 1., inclusive="left")]
        if len(controls) >= 20:
            control_offsets[night] = float(np.median(controls.time_s - p[2] * np.sqrt(controls.dec_deg)))
    events["night_offset_s"] = events.night.map(control_offsets)
    events = events[events.night_offset_s.notna()].copy()
    events["extra_time_s"] = events.time_s - events.night_offset_s - p[2] * np.sqrt(events.dec_deg)
    events["extra_gt5s"] = events.extra_time_s > 5
    events["extra_gt10s"] = events.extra_time_s > 10
    events["extra_gt30s"] = events.extra_time_s > 30
    events["ra_bin"] = pd.cut(events.ra_deg, [0, .01, .03, .1, .3, 1, 3])
    bins = events.groupby("ra_bin", observed=True).agg(
        pairs=("time_s", "size"), nights=("night", "nunique"),
        median_extra_s=("extra_time_s", "median"),
        p90_extra_s=("extra_time_s", lambda a: a.quantile(.9)),
        extra_gt5s=("extra_gt5s", "sum"), extra_gt10s=("extra_gt10s", "sum"),
        extra_gt30s=("extra_gt30s", "sum"))
    bins["extra_gt5s_fraction"] = bins.extra_gt5s / bins.pairs
    bins.to_csv(out / "conditional_ra_bins.csv")
    events.to_csv(out / "conditional_ra_events.csv", index=False)

    # Compare low and reference RA angles only within the same night and a
    # narrow Dec-angle bin, requiring observations from both groups.
    events["dec_bin"] = np.floor(events.dec_deg / .5).astype(int)
    comparisons = {}
    for label, low, reference in (
        ("lt0.1_vs_0.3to1", events.ra_deg < .1, events.ra_deg.between(.3, 1., inclusive="left")),
        ("lt1_vs_1to3", events.ra_deg < 1., events.ra_deg.between(1., 3., inclusive="left")),
        ("lt0.03_vs_0.3to1", events.ra_deg < .03, events.ra_deg.between(.3, 1., inclusive="left")),
    ):
        comparisons[label] = matched_comparison(events, low, reference, out / f"matched_{label}.csv")

    ancillary, probe_count, probe_errors = ancillary_pairs(out / "other_science_headers.jsonl")
    ancillary.to_csv(out / "ancillary_adjacent_pairs.csv", index=False)
    small_ancillary = ancillary[(ancillary.ra_deg > 1e-5) & (ancillary.ra_deg < 1) &
                                (ancillary.dec_deg < 2) & (ancillary.time_s > 0)] if len(ancillary) else ancillary
    small_ancillary.to_csv(out / "ancillary_small_ra_pairs.csv", index=False)
    result = {
        "inventory_nights": len(inventory["rows"]), "inventory_frames": sum(r["frames"] for r in inventory["rows"]),
        "inventory_start": inventory["rows"][0]["night"], "inventory_end": inventory["rows"][-1]["night"],
        "archive_header_audit": archive_summary,
        "isolated_mp_pairs": len(direct_pairs), "isolated_mp_nights": direct_pairs.night.nunique(),
        "isolated_mp_summary": direct_pairs.groupby("ra_deg").time_s.agg(["size", "median", "min", "max"]).reset_index().to_dict("records"),
        "same_field_inventory_pairs": len(geometry["same_field_transitions"]),
        "recent_source": metadata, "recent_frames": len(frames), "recent_fit_pairs": len(fit_data),
        "recent_nights": fit_data.night.nunique(),
        "normal_parameters": dict(zip(slew.MODELS["parallel_sqrt"][3], p.tolist())),
        "leave_night_out": slew.error_metrics(t, prediction),
        "matched_small_ra_comparisons": comparisons,
        "ancillary_probe_headers": probe_count, "ancillary_probe_errors": probe_errors,
        "ancillary_adjacent_science_pairs": len(ancillary),
        "ancillary_nonzero_small_ra_pairs": len(small_ancillary),
        "limits": ["A full filename/grid search locates potential field transitions; candidate movements are verified against raw headers.",
                   "Two isolated RA angles do not identify a general nonmonotone small-angle curve.",
                   "Repeated same-field exposures do not establish commanded small RA motion.",
                   "Ancillary observations have different scheduling modes and are not pooled with MP survey calibration.",
                   "Bootstrap intervals condition on the sampled nights, geometry and matched cohorts; they are not mechanical-axis response confidence bands.",
                   "Total excess waiting cannot identify which motor or software stage caused it."]}
    (out / "small_ra_results.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))
    draw_small_ra(direct_pairs, comparisons, out)
    print(json.dumps(result, ensure_ascii=False, indent=2))


def draw_small_ra(direct, comparisons, out):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
    for sign, marker, color, label in ((-1, "o", "#246A9E", "Decreasing RA"),
                                       (1, "^", "#CA722C", "Increasing RA")):
        use = np.sign(direct.signed_ra_deg) == sign
        axes[0].scatter(direct.loc[use, "ra_deg"], direct.loc[use, "time_s"],
                         marker=marker, color=color, s=40, label=label)
    axes[0].set(xlabel="|Delta RA| [deg]", ylabel="Gap minus exposure and readout [s]",
                title="Direct small-RA cases with small Dec movement", ylim=(0, None))
    axes[0].legend(fontsize=9)
    labels, positions = [], []
    for i, (key, label) in enumerate((("lt1_vs_1to3", "<1 vs 1-3 deg"),
                                     ("lt0.1_vs_0.3to1", "<0.1 vs 0.3-1 deg"),
                                     ("lt0.03_vs_0.3to1", "<0.03 vs 0.3-1 deg"))):
        row = comparisons.get(key)
        if row is None:
            continue
        mid = row["mean_matched_median_time_difference_s"]
        low, high = row["time_night_bootstrap_95pct_s"]
        axes[1].errorbar(mid, i, xerr=np.array([[mid-low], [high-mid]]),
                         fmt="o", color="#246A9E", capsize=4, lw=1.7)
        labels.append(f"{label}\n{row['small_pairs']:,} / {row['reference_pairs']:,} pairs")
        positions.append(i)
    axes[1].axvline(0, color="#555C63", ls="--", lw=1)
    axes[1].set_yticks(positions, labels, fontsize=9)
    axes[1].invert_yaxis()
    axes[1].set(xlabel="Matched median timing difference [s]",
                title="Same night and similar Dec angle\n95% whole-night bootstrap intervals")
    axes[1].margins(y=.3)
    fig.suptitle("Archived small-RA investigation | net transition times, not motor telemetry", fontsize=12)
    fig.savefig(out / "small_ra_evidence.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
