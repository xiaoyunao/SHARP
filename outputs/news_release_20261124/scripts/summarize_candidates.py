#!/usr/bin/env python3
"""Join single-night linkages to qualified raw frames and time/position-matched MPC records."""
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone
import json
import re

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
import numpy as np
import pandas as pd

from summarize import EVIDENCE, ROOT, read_stream, stem


def parse_observation(line):
    year, month, day = line[15:32].split()
    utc = datetime(int(year), int(month), 1, tzinfo=timezone.utc) + timedelta(days=float(day)-1)
    coords = SkyCoord(line[32:44].strip(), line[44:56].strip(), unit=(u.hourangle, u.deg))
    return Time(utc).mjd, coords


def mpc_identity(row, mpc):
    ras = np.fromstring(row["ras_deg"], sep=";")
    decs = np.fromstring(row["decs_deg"], sep=";")
    mjds = np.fromstring(row["mjds"], sep=";")
    source_positions = SkyCoord(ras * u.deg, decs * u.deg)
    hits = defaultdict(dict)
    for obs in mpc:
        line = obs.get("obs80") or ""
        if obs.get("status") not in {"p", "P"} or len(line) != 80 or line[77:] != "327":
            continue
        if not obs.get("iau_desig"):
            continue
        t, coord = parse_observation(line)
        match = np.flatnonzero((abs(mjds - t) * 86400 <= 20) &
                               (source_positions.separation(coord).arcsec <= 1))
        if len(match) == 1:
            hits[obs["iau_desig"]][int(match[0])] = obs
    valid = {name: matched for name, matched in hits.items() if len(matched) == len(ras)}
    if not valid:
        return "", "", "", False, 0
    assert len(valid) == 1, (row["trk_sub"], list(valid))
    name, match = next(iter(valid.items()))
    records = list(match.values())
    statuses = "".join(sorted({o["status"] for o in records}))
    refs = "|".join(sorted({o["ref"] for o in records if o.get("ref")}))
    discovery = any(o["status"] == "P" and o["obs80"][12:13] == "*" for o in records)
    return name, statuses, refs, discovery, len(match)


def main():
    keys = ["directory_night", "stem", "night", "main_sample", "quality", "counted_raw", "exposure_id"]
    e = pd.read_csv(ROOT / "exposure_ledger.csv.gz", usecols=keys, dtype={"directory_night": str, "night": str})
    frames = {(r.directory_night, r.stem): r for r in e.itertuples(index=False)}
    artifacts = set(pd.read_csv(EVIDENCE / "nov09/candidate_posthoc_artifacts.csv").trk_sub)
    false_ids = set(pd.read_csv(EVIDENCE / "nov09/candidate_known_false_sources.csv")
                    .query("is_real == 0").tracklet_id)
    comet_check = pd.read_csv(EVIDENCE / "nov09/candidate_comet_identifications.csv", dtype={"night": str})
    comet_keys = set(zip(comet_check.night, comet_check.trk_sub, comet_check.linkage_id))
    bundle = json.loads((EVIDENCE / "wamo.json").read_text())
    mpc = {k: v for entry in bundle["found"] for k, v in entry.items()}
    linkages, nightly = [], []
    for block in read_stream("candidates"):
        n = block["night"]
        if n == "20261117":
            continue
        human = {r["tracklet_id"] for r in block["human"] or [] if r.get("is_real", "").strip() == "1"}
        screened = {r["trk_sub"]: r.get("status") for r in block["morphology"] if "trk_sub" in r}
        rset = {r["trk_sub"]: r for r in block["raw"] or []}
        assert human <= set(rset), (n, human - set(rset))
        nightly.append(dict(night=n, saved_catalog=block["raw"] is not None,
                            pre_gate_count=block["pre_gate_count"], review_package=block["package_n"],
                            review_package_measurements=block["package_obs"], positive_labels=len(human),
                            review_csv_mtime=block["human_mtime"], morphology=block["morphology_counts"]))
        for r in rset.values():
            tid = r["trk_sub"]
            images = r["image_names"].split(";")
            refs = [frames.get((n, stem(img))) for img in images]
            valid = [rr for rr in refs if rr is not None and rr.quality and rr.counted_raw and rr.night == n]
            all_main = len(valid) == len(refs) and all(rr.main_sample for rr in valid)
            all_good = len(valid) == len(refs) and len({rr.exposure_id for rr in valid}) == len(refs)
            scope = "all_main_quality" if all_main and all_good else (
                "other_or_mixed_quality" if all_good else "incomplete_quality_join")
            decision = screened.get(tid, "not_screened")
            exclusion = ("archived_posthoc_artifact" if tid in artifacts else
                         "archived_human_false_source" if tid in false_ids else
                         "morphology_reject" if decision == "reject" else "")
            accepted = tid in human
            retained = accepted and not exclusion
            desig = status = refs_pub = ""
            discovery = False
            count = 0
            if retained:
                desig, status, refs_pub, discovery, count = mpc_identity(r, mpc.get(f"{tid} 327", []))
            if (n, tid, r["linkage_id"]) in comet_keys and retained and not desig:
                desig, status = "C/2025 Y1", "frozen_ephemeris_match"
            if not retained:
                identity = ""
            elif discovery:
                identity = "published_discovery_flag"
            elif not desig:
                identity = "identity_pending"
            elif re.match(r"^(?:[CPDXAI]/|\d+P(?:/|$))", desig):
                identity = "identified_comet"
            elif desig.startswith(("Jupiter ", "Saturn ")):
                identity = "identified_natural_satellite"
            else:
                identity = "identified_minor_planet"
            linkages.append(dict(night=n, trk_sub=tid, linkage_id=r["linkage_id"], scope=scope,
                                 observations=r["n_obs"], quality_member_exposures=len(valid),
                                 auto_screen=decision, human_positive=accepted, retained=retained,
                                 exclusion=exclusion if accepted else "", identity_status=identity,
                                 iau_designation=desig, mpc_status=status, mpc_measurements=count,
                                 mpc_publication=refs_pub, discovery_flag=discovery,
                                 image_names=";".join(images),
                                 identity_checked_utc=bundle["completed_utc"] if retained else ""))
    all_links = pd.DataFrame(linkages)
    all_links[all_links.human_positive].to_csv(ROOT / "candidate_review_crosswalk.csv", index=False)
    pd.DataFrame(nightly).to_csv(EVIDENCE / "candidate_stage_by_night.csv", index=False)
    records = []
    for scope, x in [("all_projects", all_links), *list(all_links.groupby("scope"))]:
        for stage, selected in [
            ("saved_automatic", x), ("human_positive", x[x.human_positive]),
            ("post_review_retained", x[x.retained]),
            *list(x[x.retained].groupby("identity_status"))]:
            records.append(dict(scope=scope, stage=stage, single_night_linkages=len(selected),
                                measurements=int(selected.observations.sum()), nights=selected.night.nunique(),
                                first_night=selected.night.min(), last_night=selected.night.max(),
                                identified_objects=selected.iau_designation.replace("", np.nan).nunique()
                                if stage.startswith(("identified_", "published_")) else None,
                                evaluated_utc=bundle["completed_utc"]))
        if not x[x.retained].discovery_flag.any():
            records.append(dict(scope=scope, stage="published_discovery_flag", single_night_linkages=0,
                                measurements=0, nights=0, identified_objects=0, evaluated_utc=bundle["completed_utc"]))
    result = pd.DataFrame(records)
    result.to_csv(ROOT / "candidate_status_summary.csv", index=False)
    print(result.to_string(index=False))
    x = all_links[all_links.retained & all_links.scope.eq("all_main_quality")]
    print("retained identity breakdown", x.groupby("identity_status").size().to_dict())
    print("comets", x[x.identity_status.eq("identified_comet")].iau_designation.value_counts().to_dict())
    print("unassigned reason", x[x.iau_designation.eq("")].mpc_status.value_counts().to_dict())
    print("v2 review packages", {r["night"]: r["review_package"] for r in nightly if r["night"] >= "20261118"})
    print("retained newly reviewed", x[x.night.ge("20261110")][["night", "trk_sub", "iau_designation"]].to_dict("records"))
    print("saved nightly gate", sum(r["pre_gate_count"] or 0 for r in nightly if not r["saved_catalog"]))
    print("fig3 tracking", x[x.trk_sub.isin(["00000P2", "00001mX"])][
        ["night", "trk_sub", "identity_status", "iau_designation", "mpc_publication"]].to_dict("records"))


if __name__ == "__main__":
    main()
