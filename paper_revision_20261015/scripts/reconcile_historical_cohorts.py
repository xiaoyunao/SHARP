#!/usr/bin/env python3
"""Join frozen detection and linkage members to physical, qualified exposures."""
import json

import numpy as np
import pandas as pd

from exposures import INPUT, TABLES
from replay_frozen_recovery import normalize_obj_id, object_key, wilson

PRODUCTS = ["orbit_links", "orbit_obs_residuals"]
HIST = ["night", "object_key"]
PHYSICAL = ["noon_night", "object_key"]
LKEY = ["night", "linkage_id"]
DKEY = ["night", "obs_key"]
MID_KEYS = ["night", "trk_sub", "linkage_id"]


def night_scope():
    base = pd.read_csv(INPUT / "night_status.csv", dtype={"night": str})
    status = pd.read_csv(INPUT / "frozen_pre_subtraction_product_status.csv", dtype={"night": str})
    status["readable"] = status.status.str.startswith("ok", na=False) & status["exists"].eq(True)
    state = status.pivot(index="night", columns="product", values="readable").reindex(columns=PRODUCTS).fillna(False)
    both = state.all(axis=1)
    base["both_orbit_products"] = base.night.isin(both.index[both])
    base["stage_state"] = np.select(
        [base.quality_code.eq("quarantine_unknown"),
         base.unknown_mpc_state.eq("not_run") & base.summary_exists.ne(True),
         ~base.both_orbit_products],
        ["historical_quarantine", "not_run_no_stage_summary", "missing_pre_subtraction_product"],
        default="stage_products_present")
    return base[["night", "both_orbit_products", "stage_state", "summary_exists", "unknown_mpc_state"]]


def exposure_catalog():
    cols = ["night", "l2_file", "exposure_key", "noon_night", "mid_utc",
            "quality", "quality_exclusion", "duplicate_time_and_name"]
    raw = pd.read_csv(TABLES / "exposure_ledger.csv.gz", usecols=cols,
                      dtype={"night": str, "noon_night": str})
    good_exposure = set(raw.loc[raw.quality, "exposure_key"])
    led = raw.loc[raw.l2_file.notna(), cols].copy().rename(columns={"l2_file": "source_file"})
    led["qualified_physical_exposure_exists"] = led.exposure_key.isin(good_exposure)
    led["midpoint_mjd"] = (pd.to_datetime(led.mid_utc, format="mixed", utc=True).astype("int64")
                           /86400e9+40587)
    assert not led.duplicated(["night", "source_file"]).any()
    return led


def attach(frame, ledger, epoch):
    d = frame.merge(ledger, on=["night", "source_file"], how="left",
                    validate="many_to_one", indicator=True)
    d["key_mapped"] = d.pop("_merge").eq("both")
    d["midpoint_delta_s"] = (d[epoch]-d.midpoint_mjd)*86400
    d["retained"] = d.key_mapped & d.quality.eq(True) & d.midpoint_delta_s.abs().le(1.0)
    d["member_state"] = np.select(
        [~d.key_mapped, d.quality.ne(True), d.midpoint_delta_s.abs().gt(1.0)],
        ["catalog_key_unmapped", d.quality_exclusion.fillna("ledger_quality_unavailable"),
         "midpoint_disagrees_gt_1s"], default="qualified")
    return d


def known_members(ledger):
    d = pd.read_parquet(INPUT / "known_mask15.parquet")
    d["night"] = d.night.astype(str).str.zfill(8)
    d["object_key"] = object_key(d)
    d["obs_key"] = d.source_file + ":" + d.obj_id.map(normalize_obj_id)
    d = d.drop_duplicates([*HIST, "obs_key"])
    return attach(d, ledger, "epoch_mjd")


def historical_links(d):
    obs = pd.read_parquet(INPUT / "orbit_obs_residuals.parquet")
    obs.night = obs.night.astype(str).str.zfill(8)
    obs["source_file"] = obs.obs_key.str.rsplit(":", n=1).str[0]
    mapping = d[[*DKEY, "object_key"]].drop_duplicates()
    counts = mapping.groupby(DKEY).object_key.transform("nunique")
    ambiguous_n = mapping.loc[counts.gt(1), DKEY].drop_duplicates().shape[0]
    mapping = mapping.loc[counts.eq(1)]
    obs = obs.merge(mapping, how="left", on=DKEY, validate="many_to_one")
    # Use the frozen algorithm's exact membership and identity rule, independent of fit flags.
    classified = (obs.groupby(LKEY, as_index=False)
                  .agg(link_obs_n=("obs_key","nunique"), mapped_n=("object_key","count"),
                       identity_n=("object_key","nunique")))
    keep = classified.identity_n.eq(1) & classified.mapped_n.eq(classified.link_obs_n) & classified.link_obs_n.ge(3)
    pure = obs.merge(classified.loc[keep, LKEY], on=LKEY, validate="many_to_one")
    unique_good = d.loc[d.retained, [*DKEY, "object_key", "retained", "noon_night", "exposure_key"]]
    pure = pure.merge(unique_good, on=[*DKEY, "object_key"], how="left", validate="many_to_one")
    return pure, ambiguous_n, int(keep.sum())


def membership_state(frame):
    return np.select([frame["unmapped"].gt(0), frame["retained"].eq(frame["members"]),
                      frame["retained"].gt(0)],
                     ["key_unmapped", "all_retained", "partially_retained"], default="all_excluded")


def group_state(df, keys):
    z = df.groupby(keys, as_index=False).agg(members=("retained", "size"),
              retained=("retained", "sum"), unmapped=("key_mapped", lambda x: int((~x).sum())))
    z["membership_state"] = membership_state(z)
    return z


def historic_cohort(d, pure, scope):
    old = pd.read_parquet(INPUT / "blinded_known_recovery_by_object.parquet")
    old = old.loc[old.eligible_ge3_and_rate, [*HIST, "strict_linked"]].copy()
    old.night = old.night.astype(str).str.zfill(8)
    hist_members = d.merge(old[HIST], on=HIST, validate="many_to_one")
    membership = group_state(hist_members, HIST)
    linked_members = pure.merge(old.loc[old.strict_linked, HIST], on=HIST, validate="many_to_one")
    linked_members = linked_members.drop_duplicates([*HIST, "obs_key"]).copy()
    linked_members["key_mapped"] = linked_members.set_index(DKEY).index.isin(
        ledger_mapped_keys(hist_members))
    linked_members["retained"] = linked_members.retained.eq(True)
    linked_state = group_state(linked_members, HIST).rename(columns={
        name: "pure_link_member_"+name for name in ("members", "retained", "unmapped", "membership_state")})
    migration = (hist_members.dropna(subset=["noon_night"]).groupby([*HIST,"noon_night"], as_index=False)
                 .agg(historical_members=("obs_key","nunique"), qualified_members=("retained","sum")))
    moved = migration.night.ne(migration.noon_night)
    a = migration.loc[moved,HIST].drop_duplicates()
    selected = migration.merge(a, on=HIST, validate="many_to_one")
    selected["archive_group_splits"] = selected.groupby(HIST).noon_night.transform("nunique").gt(1)
    inbound = migration.groupby(PHYSICAL).night.nunique().rename("archive_nights_in_physical_group")
    selected = selected.merge(inbound, on=PHYSICAL, validate="many_to_one")
    selected.to_csv(TABLES / "retention_historical_night_migration.csv",index=False)
    physical = (hist_members.dropna(subset=["noon_night"]).groupby(HIST).noon_night
                .agg(lambda x: "|".join(sorted(set(x)))).rename("mapped_physical_nights"))
    old = (old.merge(membership, on=HIST, validate="one_to_one").merge(
           linked_state, on=HIST, how="left", validate="one_to_one")
           .merge(physical, on=HIST, how="left", validate="one_to_one")
           .merge(scope, on="night", how="left", validate="many_to_one"))
    old["connection_state"] = np.select(
        [old.strict_linked, old.stage_state.eq("stage_products_present")],
        ["historically_linked", "recorded_stage_without_pure_link"], default=old.stage_state)
    stats = []
    for label, table, col in [
        ("historical_denominator", old, "membership_state"),
        ("historical_numerator_all_auxiliary_detections", old.loc[old.strict_linked], "membership_state"),
        ("historical_numerator_union_of_pure_link_members", old.loc[old.strict_linked], "pure_link_member_membership_state"),
    ]:
        counts = table.groupby([col, "stage_state"], dropna=False, as_index=False).size()
        for _, row in counts.iterrows():
            stats.append({"cohort": label, "member_outcome": row[col], "stage_state": row.stage_state,
                          "object_nights": int(row["size"])})
    pd.DataFrame(stats).to_csv(TABLES / "retention_cohort_reconciliation.csv", index=False)
    old.to_csv(TABLES / "retention_cohort_object_nights.csv.gz", index=False)
    reason = d.merge(old[HIST], on=HIST).groupby("member_state").size().to_dict()
    return old, reason


def ledger_mapped_keys(d):
    return d.loc[d.key_mapped, DKEY].drop_duplicates().set_index(DKEY).index


def physical_cohort(d, pure, scope):
    eligible_members = d.loc[d.retained].sort_values(
        ["midpoint_mjd", "epoch_mjd", "obs_key"]).drop_duplicates([*PHYSICAL, "exposure_key"])
    grouped = eligible_members.groupby(PHYSICAL)
    new = grouped.agg(exposures=("exposure_key","size"), start=("midpoint_mjd","min"), end=("midpoint_mjd","max"),
                      first_ra=("pred_ra_deg","first"), last_ra=("pred_ra_deg","last"),
                      first_dec=("pred_dec_deg","first"), last_dec=("pred_dec_deg","last")).reset_index()
    hours = (new.end-new.start)*24
    dra = (new.last_ra-new.first_ra+180)%360-180
    dx = dra*np.cos(np.deg2rad((new.first_dec+new.last_dec)/2))
    new["rate"] = np.hypot(dx, new.last_dec-new.first_dec)*3600/hours.where(hours.gt(0))
    new = new.loc[new.exposures.ge(3) & new.rate.between(3, 63)].copy()
    pure = pure.loc[pure.retained.eq(True)].copy()
    support = pure.groupby([*LKEY, *PHYSICAL]).exposure_key.nunique().rename("qualified_distinct_measurements")
    supported = support.loc[support.ge(3)].reset_index()[PHYSICAL].drop_duplicates().assign(strict_linked=True)
    new = new.merge(supported, on=PHYSICAL, how="left", validate="one_to_one")
    new["strict_linked"] = new.strict_linked.eq(True)
    stages = (eligible_members.merge(scope, on="night", validate="many_to_one").groupby(PHYSICAL)
              .agg(directory_nights=("night", lambda x: "|".join(sorted(set(x)))),
                   stage_states=("stage_state", lambda x: "|".join(sorted(set(x)))),
                   all_orbit_products=("both_orbit_products","all")).reset_index())
    new = new.merge(stages, on=PHYSICAL, validate="one_to_one")
    new["connection_state"] = np.select(
        [new.strict_linked, new.stage_states.eq("stage_products_present")],
        ["supported_by_preserved_pure_link", "recorded_stage_without_qualified_pure_link"],
        default=new.stage_states)
    details = []
    for label, cohort in [
        ("historical_stage_policy_all_archival_states", new),
        ("recorded_nonquarantined_stage_only", new.loc[new.stage_states.eq("stage_products_present")]),
    ]:
        k, n = int(cohort.strict_linked.sum()), len(cohort)
        details.append(dict(cohort=label, numerator=k, denominator=n, fraction=k/n,
                            wilson95=wilson(k,n), physical_nights=cohort.noon_night.nunique(),
                            stage_counts=cohort.stage_states.value_counts().to_dict(),
                            status_counts=cohort.connection_state.value_counts().to_dict()))
    new.to_csv(TABLES / "retention_qualified_object_nights.csv.gz", index=False)
    night_map = (eligible_members[["night", *PHYSICAL]].drop_duplicates()
                 .merge(new[PHYSICAL], on=PHYSICAL, validate="many_to_one"))
    moved = night_map.night.ne(night_map.noon_night)
    if moved.any():
        moved_nights = night_map.loc[moved, "noon_night"]
        night_map[night_map.noon_night.isin(moved_nights)].to_csv(
            TABLES / "retention_reassigned_object_nights.csv", index=False)
    return details, int(moved.sum())


def candidate_members(ledger):
    d = pd.read_csv(TABLES / "retained_measurements_midpoint.csv", dtype={"night": str, "trk_sub": str})
    catalog = pd.read_parquet(INPUT / "unknown_catalog.parquet")
    catalog["night"] = catalog.night.astype(str).str.zfill(8)
    catalog = catalog.merge(d[MID_KEYS].drop_duplicates(), on=MID_KEYS, validate="one_to_one")
    members = catalog.assign(image_member=catalog.image_names.str.split(";"),
                             obj=catalog.objids.str.split(";")).explode(["image_member","obj"])
    d["original_id"] = d.detection_key.str.rsplit("::", n=1).str[-1].map(normalize_obj_id)
    members["original_id"] = members.obj.map(normalize_obj_id)
    keys = [*MID_KEYS, "image_name", "original_id"]
    d = d.merge(members.rename(columns={"image_member":"image_name"})[keys],
                on=keys, how="left", validate="one_to_one", indicator=True)
    assert d.pop("_merge").eq("both").all() and len(members)==len(d)==179
    d["source_file"] = d.catalog_name
    d = attach(d, ledger, "corrected_midpoint_mjd_utc")
    good=d.loc[d.retained].groupby(MID_KEYS).agg(
        qualified_distinct_exposures=("exposure_key","nunique"),
        retained_physical_nights=("noon_night","nunique"),
        physical_night_labels=("noon_night",lambda x:"|".join(sorted(set(x)))),
        qualified_exposure_keys=("exposure_key",lambda x:"|".join(sorted(set(x)))))
    z=group_state(d, MID_KEYS).merge(good,on=MID_KEYS,how="left",validate="one_to_one")
    z.qualified_distinct_exposures=z.qualified_distinct_exposures.fillna(0).astype(int)
    z.retained_physical_nights=z.retained_physical_nights.fillna(0).astype(int)
    z["keeps_single_night_minimum"]=z.qualified_distinct_exposures.ge(3) & z.retained_physical_nights.eq(1)
    reason=(d.loc[~d.retained].groupby(MID_KEYS).member_state
            .agg(lambda x:"|".join(sorted(set(x)))).rename("exclusion_reasons"))
    z=z.merge(reason,on=MID_KEYS,how="left",validate="one_to_one")
    z.to_csv(TABLES / "linkage_cohort_reconciliation.csv",index=False)
    d[["night", "trk_sub", "linkage_id", "detection_index","image_name","original_id","source_file",
       "exposure_key","noon_night","mid_utc","midpoint_delta_s","retained","member_state"]].to_csv(
       TABLES / "linkage_cohort_members.csv",index=False)
    assert len(z)==58 and z.members.sum()==179
    return {"retained_measurements":int(d.retained.sum()), "retains_single_night_minimum":int(z.keeps_single_night_minimum.sum()),
            "outcomes":z.membership_state.value_counts().to_dict(),"member_states":d.member_state.value_counts().to_dict(),
            "max_midpoint_delta_s":float(d.midpoint_delta_s.abs().max())}


def main():
    led, scope = exposure_catalog(), night_scope()
    scope.to_csv(TABLES / "retention_archive_night_stage_states.csv",index=False)
    d = known_members(led)
    pure, ambiguous, strict_links = historical_links(d)
    old, reasons = historic_cohort(d, pure, scope)
    assert len(old)==102347 and old.strict_linked.sum()==75899 and ambiguous==13 and strict_links==79621
    new, moved = physical_cohort(d, pure, scope)
    result = {"historical": {"denominator":len(old),"numerator":int(old.strict_linked.sum()),
               "ambiguous_detection_keys":ambiguous,"pre_subtraction_pure_links":strict_links,
               "auxiliary_member_reasons":reasons, "object_nights":old.membership_state.value_counts().to_dict(),
               "linked_object_nights_aux_members":old.loc[old.strict_linked].membership_state.value_counts().to_dict(),
               "linked_object_nights_link_members": old.loc[old.strict_linked].pure_link_member_membership_state.value_counts().to_dict()},
              "physical_quality_cohorts":new,"reassigned_eligible_archive_object_nights":moved,
              "reviewed_linkages":candidate_members(led)}
    (TABLES / "retention_reconciliation_summary.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__=="__main__":
    main()
