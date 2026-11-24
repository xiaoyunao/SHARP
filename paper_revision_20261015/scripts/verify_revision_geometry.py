#!/usr/bin/env python3
"""Recompute geometry from frozen inputs and prove the panels retain their sample."""
import argparse
import json
import re
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from pypdf import PdfReader

import known_performance as known
from exposures import DATA, TABLES


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("draft_root", type=Path)
    args = parser.parse_args()
    # The original period analysis writes small diagnostic tables; keep those
    # working outputs away from the frozen evidence used in the manuscript.
    ledger = pd.read_csv(TABLES/"exposure_ledger.csv.gz", dtype={"night":str,"noon_night":str}, low_memory=False)
    frame = known.read_known(ledger, "known_matched")
    with tempfile.TemporaryDirectory(prefix="sharp-geometry-") as scratch:
        known.TABLES = Path(scratch)
        g = known.geometry(frame)
        pairs, audit = known.exposure_pairs(frame.loc[frame.quality])
    p = json.loads((TABLES/"performance_summary.json").read_text())
    q = g.loc[g.quality]
    assert len(g)==534780 and len(q)==531612
    assert np.isclose(g.elongation_deg.median(),p["corrected_elong_all_deg"]["median"],atol=1e-9)
    assert np.isclose(q.elongation_deg.median(),p["corrected_elong_quality_deg"]["median"],atol=1e-9)
    assert np.isclose(q.elongation_deg.min(),p["corrected_elong_quality_deg"]["min"],atol=1e-9)
    assert len(pairs)==304909 and pairs.separation_s.gt(0).all()
    assert pairs.exposure_key.ne(pairs.previous_exposure_key).all()
    for name, source in [("elongation",q.elongation_deg),
                          ("predicted_rate",pairs.predicted_rate_arcsec_h),
                          ("measured_rate",pairs.measured_rate_arcsec_h)]:
        b = pd.read_csv(DATA/f"figA1_{name}.csv")
        edges = np.r_[b.lower,b.upper.iloc[-1]]
        assert np.array_equal(np.histogram(source, bins=edges)[0],b["count"]),name
    fast=pairs.loc[pairs.predicted_rate_arcsec_h.gt(1000)]
    assert len(fast)==5 and fast.measured_rate_arcsec_h.gt(1000).all()
    assert set(fast.asteroid_name)=={"2026 AA","2026 CA"}
    path=args.draft_root/"manuscript/manuscript.pdf"
    reader=PdfReader(path)
    assert len(reader.pages)==16
    rendered=" ".join(page.extract_text() for page in reader.pages)
    rendered = re.sub(r"(?<=\d)\s+(?=\.\d)", "", rendered)
    for s in ["39,590","531,612","74.16","149.34","304,909","20.25","0.263"]:
        assert s in rendered,s
    return {
        "manuscript_pages":len(reader.pages),
        "baseline_elongation_median_deg":float(g.elongation_deg.median()),
        "qualified_elongation_median_deg":float(q.elongation_deg.median()),
        "qualified_elongation_min_deg":float(q.elongation_deg.min()),
        "rate_pair_audit":audit,"angular_rate_pair_count":len(pairs),"fast_pairs_agree":len(fast),
        "appendix_figure_histograms_match":True,
    }


if __name__=="__main__":
    result=main()
    (TABLES/"revision_geometry_verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
