#!/usr/bin/env python3
"""Record the cited source and statistic definition for each draft claim and caption."""
from pathlib import Path
import argparse
import hashlib
import json

import pandas as pd

from exposures import ROOT, TABLES


def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manuscript_root", type=Path, help="Extracted GPT Pro project with manuscript/ and evidence/")
    args = parser.parse_args()
    src = args.manuscript_root
    astro = pd.read_csv(TABLES/"astrometry_summary.csv")
    qualified = astro.loc[astro["sample"].eq("qualified_1arcsec")].set_index("metric")
    baseline = astro.loc[astro["sample"].eq("baseline_1arcsec")].set_index("metric")
    perf = json.loads((TABLES/"performance_summary.json").read_text())
    obs = json.loads((TABLES/"observing_summary.json").read_text())
    areas = pd.read_csv(TABLES/"geometric_coverage.csv").set_index("nside")
    mod = perf["l2_quality_header_metrics"]["limiting_mag"]
    coeff = "e71f9f6a2ea95205e0a30230039ea65b50949914"
    raw = perf["known"]
    ph = perf["retention"]
    can = perf["candidates"]
    s = "qualified 20260803 commissioning exposures"
    h = "historical 20260803 processing cohort"
    rows = []

    def row(position, quantity, unit, value, observation, criteria, scope, file, generator):
        rows.append(dict(manuscript_location=position, quantity=quantity, units=unit, value=value,
                         statistical_unit=observation, denominator_or_filter=criteria, cohort=scope,
                         source_file_and_key=file, generation_script=generator))

    exposure_script = "scripts/exposures.py"
    known_script = "scripts/known_performance.py"
    row("Abstract; Sections 2, 4, 7; Table 2","Exposure count / physical nights","count",
        f"{obs['quality_exposures']} / {obs['quality_noon_nights']}", "exposure / physical night",
        "valid L2 WCS and midpoint; unique exposure_key",s,"tables/observing_summary.json:quality_exposures,quality_noon_nights",exposure_script)
    row("Abstract; Sections 4, 7; Table 2","Open-shutter time","h",obs["quality_shutter_hours"],"exposure duration",
        "sum actual exposure_s / 3600",s,"tables/observing_summary.json:quality_shutter_hours",exposure_script)
    for col, label, subset in [("geometric_once_deg2","Geometric area with ≥1 exposure","sky union, any exposure"),
                               ("geometric_ge3_same_night_deg2","Geometric area with ≥3 same-night visits",
                                "union over nights after thresholding each physical night independently")]:
        row("Abstract; Section 4; Table 2",label,"deg2",areas.loc[1024,col],"HEALPix area",
            subset,s,f"tables/geometric_coverage.csv:nside=1024,{col}",exposure_script)
    row("Section 4; Table 2; Fig. 4c","At least three visits","groups; %","9710 / 14400 = 67.43056",
        "pointing group within field label and physical night","qualified ≥3 / all qualified groups; 0.30-deg split",
        s,"tables/observing_summary.json:field_nights_visits,field_nights; tables/cadence_pointing_tolerance.csv:0.3",exposure_script)
    for label, key, flt in [("Median adjacent-visit interval","adjacent_min_p16_p50_p84","adjacent pair in a qualified pointing group"),
                            ("Median first-last span","ge3_span_min_p16_p50_p84","pointing group with ≥3 visits")]:
        row("Abstract; Sections 4, 7; Table 2; Fig. 4",label,"min",obs[key][1],flt,"median physical midpoint difference",
            s,f"tables/observing_summary.json:{key}[1]",exposure_script)
    row("Abstract; Sections 4, 7; Table 2; Fig. 5a","Median model depth and IQR","mag",
        f"{mod['median']}; {mod['p25']}–{mod['p75']}","exposure",
        "historical MLIM with reference-aperture CCDZP and a FWHM/2 background-noise radius",
        s,"tables/performance_summary.json:l2_quality_header_metrics.limiting_mag; inputs/depth_definition_note.md",known_script)
    sees = perf["l2_quality_header_metrics"]["seeing_arcsec"]
    row("Section 4 image quality; Table 2","Exposure seeing median and IQR","arcsec",
        f"{sees['median']}; {sees['p25']}–{sees['p75']}","exposure","L2 seeing_arcsec",s,
        "tables/performance_summary.json:l2_quality_header_metrics.seeing_arcsec",known_script)
    row("Abstract; Sections 4, 5, 7; Table 2","Primary associations / identities","count",
        f"{raw['qualified_detected']} / {raw['qualified_identified']}","detection / catalog identity",
        "qualified exposures with archived 1-arcsec nearest matches",s,
        "tables/performance_summary.json:known.qualified_detected,qualified_identified",known_script)
    cls=pd.read_csv(TABLES/"orbital_classes.csv")
    mb=int(cls.loc[cls.cohort.eq("quality") & cls["class"].eq("Main belt"),"unique_objects"].iloc[0])
    row("Abstract; Sections 4, 5; Table 2; Fig. 7","Main-belt-region catalog identities","identities; %",f"{mb}; {100*mb/raw['qualified_identified']:.5f}",
        "catalog identity","main-belt region and unique identities",s,"tables/orbital_classes.csv:quality,Main belt",known_script)
    for label, val, crit, source in [
        ("Median radial O−C",qualified.loc["radial_arcsec","median"],"qualified 1 arcsec, n=531612","qualified_1arcsec,radial_arcsec,median"),
        ("Historical radial O−C comparison",baseline.loc["radial_arcsec","median"],"historical 1 arcsec, n=534780","baseline_1arcsec,radial_arcsec,median"),
        ("Median signed RA*cos(Dec) / Dec",
         f"{qualified.loc['dra_arcsec','median']} / {qualified.loc['ddec_arcsec','median']}","qualified 1 arcsec","qualified_1arcsec,dra_arcsec/ddec_arcsec,median"),
        ("Robust signed RA*cos(Dec) / Dec scatter",
         f"{qualified.loc['dra_arcsec','sigma_iqr']} / {qualified.loc['ddec_arcsec','sigma_iqr']}",
         "qualified 1 arcsec; 0.7413 IQR","qualified_1arcsec,dra_arcsec/ddec_arcsec,sigma_iqr"),
    ]:
        row("Abstract; Section 4; Table 2; Fig. 5c–d",label,"arcsec",val,"association",crit,
            h if label.startswith("Historical") else s,f"tables/astrometry_summary.csv:{source}",known_script)
    row("Section 4 astrometry","Auxiliary count / outside 1 arcsec","associations",
        f"{perf['astrometry']['quality_mask15_n']} / {perf['astrometry']['quality_mask15_gt1']}",
        "archived auxiliary association","quality selected with 1.5 arcsec radius; outside means r>1 arcsec",s,
        "tables/performance_summary.json:astrometry; tables/astrometry_summary.csv:qualified_1p5arcsec",known_script)
    row("Abstract; Sections 3, 4, 7; Table 2; Fig. 6","Conditional connection retention","object-nights; %",
        f"{ph['k']} / {ph['n']} = {ph['fraction']*100:.5f}","historical catalog identity-night",
        "≥3 detections; 3–63 arcsec/h; ≥1 archived pure pre-subtraction link",h,
        "inputs/blinded_known_summary.json:survival_for_ge3_and_rate; figure_data/fig06_retention_frozen_bins.csv", "scripts/replay_frozen_recovery.py")
    row("Sections 4, 5; Tables 2–3; Fig. A2","Candidate review sequence","linkages",
        f"{can['auto_links']} / {can['marked_real']} / {can['selected']} / {can['retained']}",
        "single-night linkage","automatic / initially marked real / selected / retained",h,
        "tables/performance_summary.json:candidates; tables/candidate_status_by_night.csv",known_script)
    row("Sections 4, 5; Tables 2–3; Appendix B","Retained measurements / nights / comet-consistent linkages","count",
        f"{can['retained_measurements']} / {can['retained_nights']} / {can['c2025y1_linkages']}",
        "measurement / physical night / linkage","58 reviewed historical linkages",h,
        "tables/retained_measurements_midpoint.csv; tables/retained_linkages.csv:is_c2025_y1",known_script)
    for key, clabel, loc in [("corrected_elong_quality_deg",s,"Section 5; Appendix A; Fig. A1b"),
                             ("corrected_elong_all_deg",h,"Historical sample cross-check only")]:
        val=perf[key]
        row(loc,"Corrected median solar elongation","deg",val["median"],"association",
            f"topocentric, {val['n']} archived identities per exposure",clabel,
            f"tables/performance_summary.json:{key}.median",known_script)
    pred, meas = perf["predicted_rate_arcsec_h"],perf["measured_rate_arcsec_h"]
    row("Appendix A; Fig. A1a","Angular rate pair count / median predicted and measured","pairs; arcsec/h",
        f"{pred['n']}; {pred['median']} / {meas['median']}","adjacent detection pair",
        "same physical night and identity; different exposures; positive elapsed time",s,
        "tables/performance_summary.json:predicted_rate_arcsec_h,measured_rate_arcsec_h",known_script)
    row("Appendix A; Fig. A1a","Predicted speed above 1000","pairs",5,"adjacent detection pair",
        "all five supported by measured displacement; 2026 AA/2026 CA",s,
        "tables/motion_tail_audit.csv:predicted_rate_arcsec_h>1000",known_script)

    fig_checks = json.loads((src/"evidence/figure_source_sha256.json").read_text())["figures"]
    caption = [
        ("Fig. 1","raw field and labeled cutouts", "inputs/image_source_manifest.csv; figure_data/fig01_sources.csv", "historical L1; time annotation is exposure start"),
        ("Fig. 2","historical processing and review", "figure_data/fig02_workflow_nodes.csv; fig02_workflow_edges.csv", "known and search branches; pre-subtraction orbit fits"),
        ("Fig. 3","WCS visits and detection-weighted density", "figure_data/fig03_wcs_visits_icrs_nside512.csv.gz; fig03_known_ecliptic_nside128.csv", "WCS visit pixels vs primary association pixels"),
        ("Fig. 4","counts and cadence", "figure_data/fig04_nightly_activity.csv; fig04_monthly_activity.csv; fig04_field_night_counts.csv; fig04_adjacent.csv; fig04_span_ge3.csv", "panels a,b raw and qualified by directory date; c–e midpoint physical nights"),
        ("Fig. 5","measurement diagnostics", "figure_data/fig05_limiting_magnitude.csv; fig05_photometry_vs_mag.csv; fig05_astrometry_vs_mag.csv; fig05_astrometry_by_night.csv; tables/stellar_repeatability_vs_magnitude.csv", "b known Flag=0, positive aperture flux and error; stellar sample 18 exposures in six sequences; c,d qualified 1 arcsec"),
        ("Fig. 6","historical conditional retention", "figure_data/fig06_retention_frozen_bins.csv", "historical object-night, unchanged sample and denominator"),
        ("Fig. 7","catalog orbit regions and repeat sampling", "figure_data/fig07_objects_orbits_revisits.csv.gz", "one record per catalog identity; distant class is geometric orbital region"),
        ("Fig. 8","20260210 link 35","figure_data/fig06_astrometry.csv; tables/retained_measurements_midpoint.csv:trk_sub=00000P2", "three original frames and midpoint cadence"),
        ("Fig. A1","known angular rates and solar elongation","figure_data/figA1_predicted_rate.csv; figA1_measured_rate.csv; figA1_elongation.csv", "rates weighted per adjacent pair; elongation per quality-selected association"),
        ("Fig. A2","retained linkage residuals and arc spans","figure_data/figA2_candidate_rms_arc.csv", "historical retained linkage, arc span from corrected midpoints"),
    ]
    evidence = []
    assert len(fig_checks)==len(caption)==10
    for c, verified in zip(caption,fig_checks):
        label, name, file, flt = c
        fig = verified["file"]
        assert sha(ROOT/"figures"/fig)==sha(src/"manuscript/figures"/fig)==verified["sha256"], fig
        row(label,name,"caption","verified","panel-specific",flt,"caption-defined",file,"scripts/make_figures.py or archived original image")
        evidence.append({"figure":label,"pdf":fig,"sha256":verified["sha256"],"pinned_commit":coeff})
    pd.DataFrame(evidence).to_csv(TABLES/"manuscript_figure_identity.csv",index=False)
    pd.DataFrame(rows).to_csv(TABLES/"manuscript_number_source_map.csv",index=False)
    print(len(rows),"manuscript entries; verified all 10 manuscript PDF hashes")


if __name__=="__main__":
    main()
