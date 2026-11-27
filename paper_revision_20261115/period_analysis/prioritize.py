"""Post-search triage only. Ranking is a review priority, not a reliability class."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def prioritize(results):
    summary = pd.read_csv(results / "pilot_summary.csv")
    hypotheses = pd.read_csv(results / "period_hypotheses.csv")
    points = pd.read_csv(results / "pilot_points.csv.gz", dtype={"physical_night": str})
    refs = pd.read_csv(results / "reference_comparison.csv")
    rows = []
    for item in summary.itertuples():
        d = points.loc[points.object_key.eq(item.object_key)]
        null_chi2 = 0.
        for _, group in d.groupby("physical_night"):
            w = 1 / group.fit_err.to_numpy() ** 2
            mean = np.average(group.mag, weights=w)
            null_chi2 += np.sum(w * (group.mag.to_numpy() - mean) ** 2)
        choices = hypotheses.loc[hypotheses.object_key.eq(item.object_key)].sort_values("rank")
        best = choices.iloc[0]
        improvement = null_chi2 - best.chi2 - 5 * np.log(len(d))
        margin = choices.iloc[1].delta_bic if len(choices) > 1 else np.nan
        rows.append(dict(object_key=item.object_key, delta_bic_vs_night_constant=improvement,
                         next_search_basin_delta_bic=margin,
                         amplitude_over_rms=item.amplitude_mag / item.rms_mag,
                         at_refinement_edge=bool(best.refinement_boundary)))
    triage = summary.merge(pd.DataFrame(rows), on="object_key", validate="one_to_one").merge(
        refs[["object_key", "fullname", "sbdb_period_h", "interpretation"]],
        on="object_key", validate="one_to_one")
    triage["review_gate"] = ((triage.delta_bic_vs_night_constant > 25)
                             & (triage.next_search_basin_delta_bic >= 10)
                             & (triage.amplitude_over_rms >= 3)
                             & (triage.max_phase_gap < .15)
                             & triage.candidate_period_h.between(1.1, 100)
                             & ~triage.at_refinement_edge)
    triage["sbdb_period_absent"] = triage.interpretation.eq("no_sbdb_period_is_not_evidence_of_novelty")
    triage = triage.sort_values(["review_gate", "delta_bic_vs_night_constant"], ascending=False)
    triage.to_csv(results / "review_priority.csv", index=False)
    shortlist = triage.loc[triage.review_gate & triage.sbdb_period_absent]
    shortlist.to_csv(results / "missing_sbdb_review_queue.csv", index=False)
    info = dict(screened=len(triage), sbdb_period_present=int(triage.sbdb_period_h.notna().sum()),
                sbdb_period_absent=int(triage.sbdb_period_absent.sum()),
                lookup_failed=int(triage.interpretation.eq("lookup_failed").sum()),
                heuristic_review_gate=int(triage.review_gate.sum()),
                missing_sbdb_review_queue=len(shortlist),
                gate="deltaBIC vs night constant>25; next basin deltaBIC>=10; amplitude/RMS>=3; phase gap<0.15; 1.1<=P<=100h; no refinement edge",
                warning="Heuristic triage without a global false-positive calibration; ISO flags/geometry/images/LCDB still required",
                confirmed_new_periods=0)
    (results / "screening_summary.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2))
    print(shortlist.head(10)[["object_key", "candidate_period_h", "n_points", "delta_bic_vs_night_constant"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    prioritize(parser.parse_args().results)
