#!/usr/bin/env python3
"""Choose only matched exposure and same-night, same-field comparison catalogs."""
import argparse
import json

import pandas as pd

from exposures import INPUT, TABLES


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("jobs")
    args = parser.parse_args()
    calib = pd.read_csv(TABLES / "sample_photometric_header_and_depth_check.csv")
    jobs = {row.path: dict(path=row.path, apertures=True, sha256=row.sha256, targets=[])
            for row in calib.itertuples()}
    dtype = {"night": str, "noon_night": str, "field_id": str}
    ledger = pd.read_csv(TABLES / "exposure_ledger.csv.gz", dtype=dtype, low_memory=False)
    target = pd.read_csv(INPUT / "distant_frozen_predictions.csv.gz", dtype=dtype)
    ids = set(target.orbit_key)
    members = pd.read_csv(TABLES / "rare_object_detections.csv", dtype=dtype)
    members = members.loc[members.orbit_key.isin(ids)].copy()
    guide = ledger.loc[ledger.l2_file.notna(), ["night", "l2_file", "path", "size_bytes", "field_id", "mjd", "quality", "noon_night"]]
    target = target.merge(guide, left_on=["night", "source_file"], right_on=["night", "l2_file"],
                          how="left", validate="many_to_one")
    members = members.merge(guide, left_on=["night", "source_file"], right_on=["night", "l2_file"],
                            how="left", validate="many_to_one")
    assert len(members) == 119 and members.path.notna().all()
    members["member"] = members.index
    for row in members.itertuples():
        # The exact frozen ephemeris for a sibling must correspond to the same field
        # and verified physical night. Includes catalogs with no historical 1-arcsec match.
        siblings = target.loc[target.orbit_key.eq(row.orbit_key) & target.noon_night.eq(row.noon_night) &
                              target.field_id.eq(row.field_id) & target.quality.eq(True)].copy()
        # Eris has little displacement within several of its same-night sequences.
        # Follow its nearest same-field sequences on other physical nights using the
        # already-frozen ephemeris. Also check the three unresolved singletons.
        if row.orbit_key in {"N:136199", "N:523716", "N:562274", "N:84922"}:
            delta = target.mjd - row.mjd
            nearby = target.loc[target.orbit_key.eq(row.orbit_key) & target.field_id.eq(row.field_id) &
                                target.noon_night.ne(row.noon_night) & target.quality.eq(True) &
                                delta.abs().between(.25, 14)].copy()
            if len(nearby):
                nearby["dt"] = (nearby.mjd-row.mjd).abs()
                nearby["side"] = nearby.mjd.gt(row.mjd)
                siblings = pd.concat([siblings, nearby.sort_values("dt").drop_duplicates("side")])
        siblings = siblings.sort_values(["mjd", "path"]).drop_duplicates("path")
        if row.path not in set(siblings.path):
            siblings = pd.concat([siblings, target.loc[target.orbit_key.eq(row.orbit_key) & target.path.eq(row.path)]])
        for s in siblings.itertuples():
            job = jobs.setdefault(s.path, dict(path=s.path, size=int(s.size_bytes), mjd=s.mjd, targets=[]))
            for kind, ra, dec in [("prediction", s.pred_ra_deg, s.pred_dec_deg),
                                  ("anchor", row.obs_ra_win_deg, row.obs_dec_win_deg)]:
                job["targets"].append(dict(member=row.member, kind=kind, ra=ra, dec=dec,
                                           frozen_mjd=s.epoch_mjd, comparison_night=s.noon_night))
    with open(args.jobs, "w") as handle:
        json.dump(list(jobs.values()), handle, indent=2)
    members.to_csv(TABLES / "distant_matched_members.csv", index=False)
    print(len(jobs), "catalogs,", sum(bool(j["targets"]) for j in jobs.values()), "distant comparison catalogs")


if __name__ == "__main__":
    main()
