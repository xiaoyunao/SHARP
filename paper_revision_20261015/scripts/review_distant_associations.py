#!/usr/bin/env python3
"""Compare the moving and stationary hypotheses around frozen distant-region matches."""
import gzip
import json

import numpy as np
import pandas as pd

from exposures import INPUT, TABLES
from known_performance import angular

FIELDS = ("obj_id", "g_aper", "g_aper_error", "flag", "iso_flag", "fwhm", "x_win", "y_win")


def check_depth(rows):
    checks, growth = [], []
    for r in rows:
        if "apertures" not in r:
            continue
        a = r["apertures"]
        k = 5*np.sqrt((a["SKY"]*a["GAIN"]+a["RON"]**2)*np.pi*(a["FWHM"]/2)**2)/a["GAIN"]
        per = pd.DataFrame(r["growth"])
        # FWHM diameters rarely coincide with a stored fixed aperture. Report bracketing
        # growth-curve points rather than imposing a PSF shape.
        lower = per.loc[per.diameter_px.le(a["FWHM"])].tail(1)
        upper = per.loc[per.diameter_px.ge(a["FWHM"])].head(1)
        clip = lambda z, col: float(z[col].iloc[0]) if len(z) else np.nan
        checks.append({"path":r["path"], **a, "mlim_replayed":a["CCDZP"]-2.5*np.log10(k),
                       "background_5sigma_flux_adu":k,
                       "bracket_low_diameter_px":clip(lower,"diameter_px"),
                       "bracket_high_diameter_px":clip(upper,"diameter_px"),
                       "bracket_low_median_F_over_F18":clip(lower,"median"),
                       "bracket_high_median_F_over_F18":clip(upper,"median")})
        growth.extend([{"path":r["path"], **p} for p in r["growth"]])
    checks = pd.DataFrame(checks)
    assert len(checks)==18 and checks.stars.min()>0
    assert (checks.mlim_replayed-checks.MLIM).abs().max()<.006
    assert checks.max_abs_mag4_formula_residual.max()<.00005
    checks.to_csv(TABLES / "depth_aperture_reconciliation_18.csv",index=False)
    pd.DataFrame(growth).to_csv(TABLES / "depth_aperture_growth_curves_18.csv",index=False)


def classify_members(members, results, jobs):
    targets = []
    files = {r["path"]:r for r in results}
    for job in jobs:
        for t in job["targets"]:
            if t["kind"] == "prediction":
                targets.append({"path":job["path"], **t})
    comparisons = pd.DataFrame(targets).drop_duplicates(["path","member"])
    hit_rows = [{"path":r["path"], **{k:v for k,v in t.items() if k!="data"},
                 **dict(zip(FIELDS,t["data"]))} for r in results for t in r.get("sources",[])]
    hits = pd.DataFrame(hit_rows)
    hits.to_csv(TABLES / "distant_nearby_l2_sources.csv",index=False)
    anchor = hits.loc[hits.kind.eq("anchor")]
    prediction = hits.loc[hits.kind.eq("prediction")]
    nearest_a = anchor.sort_values("sep_arcsec").drop_duplicates(["path","member"])
    nearest_p = prediction.sort_values("sep_arcsec").drop_duplicates(["path","member"])
    compare = comparisons.merge(nearest_a[["path","member","sep_arcsec","g_aper","src_ra","src_dec"]].rename(
                columns={k:f"fixed_{k}" for k in ("sep_arcsec","g_aper","src_ra","src_dec")}),
                on=["path","member"], how="left", validate="one_to_one")
    compare = compare.merge(nearest_p[["path","member","sep_arcsec","g_aper","src_ra","src_dec"]].rename(
                columns={k:f"prediction_{k}" for k in ("sep_arcsec","g_aper","src_ra","src_dec")}),
                on=["path","member"], how="left", validate="one_to_one")
    coverage = pd.DataFrame([{"path": r["path"], **c} for r in results for c in r.get("coverage",[])])
    both = coverage.pivot(index=["path","member"],columns="kind",values="inside").all(axis=1)
    compare = compare.merge(both.rename("both_positions_on_detector"), on=["path","member"], validate="one_to_one")
    compare["catalog_read"] = compare.path.map(lambda p:p in files and "error" not in files[p])
    compare = compare.merge(members[["member","path","pred_ra_deg","pred_dec_deg","obs_ra_win_deg",
                                     "obs_dec_win_deg","mag_aper4","epoch_mjd","noon_night"]].rename(columns={"path":"matched_path"}),
                            on="member",validate="many_to_one")
    compare["other_exposure"] = compare.path.ne(compare.matched_path)
    compare["elapsed_min"] = (compare.frozen_mjd-compare.epoch_mjd)*1440
    compare["prediction_motion_arcsec"] = angular(compare.pred_ra_deg,compare.pred_dec_deg,compare.ra,compare.dec)
    compare["prediction_to_anchor_arcsec"] = angular(compare.obs_ra_win_deg,compare.obs_dec_win_deg,compare.ra,compare.dec)
    compare["predicted_source_displacement_arcsec"] = angular(compare.obs_ra_win_deg,compare.obs_dec_win_deg,
                                                            compare.prediction_src_ra,compare.prediction_src_dec)
    dx = lambda a,b: (a-b+180)%360-180
    error_x = (dx(compare.prediction_src_ra,compare.obs_ra_win_deg)-dx(compare.ra,compare.pred_ra_deg))
    error_y = (compare.prediction_src_dec-compare.obs_dec_win_deg)-(compare.dec-compare.pred_dec_deg)
    compare["displacement_resid_arcsec"] = 3600*np.hypot(error_x*np.cos(np.deg2rad(compare.dec)),error_y)
    compare["same_physical_night"] = compare.noon_night.eq(compare.comparison_night)
    sibling = (compare.other_exposure & compare.catalog_read & compare.both_positions_on_detector &
               compare.elapsed_min.abs().gt(0.5))
    compare["stationary_support"] = (sibling & compare.prediction_to_anchor_arcsec.ge(1.5)
                                      & compare.fixed_sep_arcsec.le(.7)
                                      & (compare.fixed_g_aper-compare.mag_aper4).abs().le(.65))
    sufficient_motion = np.where(compare.same_physical_night, 1.5, 5.0)
    compare["moving_support"] = (sibling & compare.prediction_motion_arcsec.ge(sufficient_motion)
                                   & compare.prediction_sep_arcsec.le(1.05)
                                   & compare.predicted_source_displacement_arcsec.ge(1.1)
                                   & compare.displacement_resid_arcsec.le(1.2)
                                   & (compare.prediction_g_aper-compare.mag_aper4).abs().le(1.0)
                                   & ~compare.stationary_support)
    compare.to_csv(TABLES / "distant_l2_trial_exposures.csv",index=False)
    own = compare.loc[~compare.other_exposure].set_index("member")
    match_ok = own.catalog_read & own.fixed_sep_arcsec.le(.02)
    assert own.index.is_unique and len(own)==119
    count = compare.loc[compare.other_exposure].groupby("member").agg(
        sibling_catalogs=("catalog_read","sum"),
        stationary_support_frames=("stationary_support","sum"),
        moving_support_frames=("moving_support","sum"),
        max_predicted_motion_arcsec=("prediction_motion_arcsec","max"))
    out = members.merge(count,on="member",how="left",validate="one_to_one")
    out["verified_l2_source"] = out.member.map(match_ok).eq(True)
    for c in count.columns[:3]:
        out[c] = out[c].fillna(0).astype(int)
    out["pred_V_minus_g_aper"] = out.pred_v_mag-out.mag_aper4
    out["magnitude_discrepancy_gt_2"] = out.pred_V_minus_g_aper.abs().gt(2)
    s=out.stationary_support_frames.gt(0)
    m=out.moving_support_frames.gt(0)
    out["association_state"] = np.select(
        [~out.verified_l2_source, s & m, s, m & ~out.magnitude_discrepancy_gt_2],
        ["key_unmapped", "mixed_local_evidence", "stationary_coincidence",
         "moving_source_supported"], default="insufficient")
    return out, compare


def summarize(members, audit):
    results = []
    for key, d in members.groupby("orbit_key"):
        states = d.association_state.value_counts()
        nmove = int(states.get("moving_source_supported",0))
        nstatic = int(states.get("stationary_coincidence",0))
        state = ("moving_source_supported" if nmove else
                 "stationary_coincidence" if nstatic else
                 "insufficient" if d.verified_l2_source.any() else "key_unmapped")
        row = audit.loc[key]
        results.append(dict(orbit_key=key,identity=row.asteroid_name, a_au=row.a_au,
                            historical_associations=len(d), matched_images=d.path.nunique(),
                            moving_supported_associations=nmove, stationary_coincidences=nstatic,
                            mixed_evidence=int(states.get("mixed_local_evidence",0)),
                            insufficient=int(states.get("insufficient",0)),
                            key_unmapped=int(states.get("key_unmapped",0)),
                            brightness_flags=int(d.magnitude_discrepancy_gt_2.sum()),
                            maximum_abs_pred_V_minus_g_aper=d.pred_V_minus_g_aper.abs().max(),
                            identity_evidence=state))
    return pd.DataFrame(results)


def main():
    jobs=json.loads((INPUT/"followup_l2_jobs.json").read_text())
    with gzip.open(INPUT/"followup_l2_evidence.jsonl.gz","rt") as handle:
        rows=[json.loads(line) for line in handle]
    assert len(rows)==len(jobs) and {r["path"] for r in rows}=={j["path"] for j in jobs}
    check_depth(rows)
    base=pd.read_csv(TABLES/"distant_matched_members.csv",dtype={"night":str,"noon_night":str,"field_id":str})
    audit=pd.read_csv(TABLES/"rare_object_audit.csv").set_index("object_key")
    members, comparisons=classify_members(base,rows,jobs)
    members.to_csv(TABLES/"distant_association_local_checks.csv",index=False)
    identities=summarize(members,audit)
    assert len(identities)==17 and identities.historical_associations.sum()==119
    identities.to_csv(TABLES/"rare_object_local_checks.csv",index=False)
    print(identities[["identity","historical_associations","moving_supported_associations",
                      "stationary_coincidences","insufficient","identity_evidence"]].to_string(index=False))
    print(members.association_state.value_counts().to_string())


if __name__=="__main__":
    main()
