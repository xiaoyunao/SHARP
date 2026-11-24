#!/usr/bin/env python3
"""Extract frozen predictions only for the paper's 17 distant-region identities."""
import argparse
import hashlib

import pandas as pd
import pyarrow.parquet as pq

from exposures import INPUT, TABLES
from known_performance import canon_orbit_id
from replay_frozen_recovery import object_key


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("known_all")
    args = parser.parse_args()
    a = pd.read_csv(TABLES / "rare_object_audit.csv")
    wanted = set(a.loc[a.orbit_class.str.startswith("Distant"), "object_key"])
    assert len(wanted) == 17
    columns = ["night", "source_file", "asteroid_number", "asteroid_name",
               "epoch_mjd", "pred_ra_deg", "pred_dec_deg", "pred_v_mag"]
    results = []
    with open(args.known_all, "rb") as handle:
        print("Archive SHA-256", hashlib.file_digest(handle, "sha256").hexdigest())
    file = pq.ParquetFile(args.known_all)
    for batch in file.iter_batches(columns=columns, batch_size=250_000):
        rows = batch.to_pandas()
        key = canon_orbit_id(object_key(rows))
        # Unnumbered catalog identifiers are stored normalized as S:<designation>.
        missing = key.str.startswith("S:")
        key.loc[missing] = key.loc[missing].str.replace(r"\s+", "", regex=True).str.lower().str.replace("s:", "S:", n=1)
        use = key.isin(wanted)
        if use.any():
            results.append(rows.loc[use].assign(orbit_key=key.loc[use]))
    result = pd.concat(results).drop_duplicates()
    assert set(result.orbit_key) == wanted
    result.to_csv(INPUT / "distant_frozen_predictions.csv.gz", index=False)
    print(result.groupby("orbit_key").size().to_string(), "total", len(result),
          "distinct exposures", len(result[["night","source_file"]].drop_duplicates()))


if __name__ == "__main__":
    main()
