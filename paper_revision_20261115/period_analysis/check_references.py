"""Optional POST-fit SBDB comparison; never feeds the blind search or rankings."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import urllib.parse
import urllib.request

import pandas as pd


def check(results):
    summary = pd.read_csv(results / "pilot_summary.csv")
    output = results / "reference_cache"
    output.mkdir(exist_ok=True)
    rows = []
    # Batch queries retain stable numbered identities and avoid one call per body.
    for start in range(0, len(summary), 40):
        chunk = summary.iloc[start:start + 40]
        numbers = [key.removeprefix("N:") for key in chunk.object_key]
        filters = {"OR": ["pdes|EQ|" + number for number in numbers]}
        url = "https://ssd-api.jpl.nasa.gov/sbdb_query.api?" + urllib.parse.urlencode({
            "fields": "pdes,full_name,rot_per", "sb-cdata": json.dumps(filters, separators=(",", ":"))})
        path = output / f"sbdb_query_{start:04d}.json"
        try:
            if path.exists():
                cache = json.loads(path.read_text())
                if cache["query_url"] != url:
                    raise ValueError("Cached query differs; use a new results directory")
            else:
                with urllib.request.urlopen(url, timeout=20) as response:
                    data = json.load(response)
                if data.get("fields") != ["pdes", "full_name", "rot_per"]:
                    raise ValueError("Unexpected SBDB fields or error response")
                cache = dict(query_url=url, retrieved_utc=datetime.now(timezone.utc).isoformat(), response=data)
                path.write_text(json.dumps(cache, indent=2) + "\n")
            data = cache["response"]
            mapping = {str(values[0]): values for values in data.get("data", [])}
            if set(mapping) != set(numbers):
                raise ValueError("SBDB identities missing or unexpected")
            for item in chunk.itertuples():
                number = item.object_key.removeprefix("N:")
                _, fullname, raw_period = mapping[number]
                period = float(raw_period) if raw_period is not None else None
                rows.append(dict(object_key=item.object_key, fullname=fullname.strip(), sbdb_period_h=period,
                                 candidate_period_h=item.candidate_period_h,
                                 candidate_over_reference=item.candidate_period_h / period if period else None,
                                 doubled_candidate_over_reference=2 * item.candidate_period_h / period if period else None,
                                 interpretation="reference_exists_not_new_by_this_check" if period else
                                 "no_sbdb_period_is_not_evidence_of_novelty", query_url=url,
                                 retrieved_utc=cache["retrieved_utc"],
                                 api_version=data.get("signature", {}).get("version")))
                print(number, "SBDB_period_h", period, flush=True)
        except Exception as exc:
            for item in chunk.itertuples():
                rows.append(dict(object_key=item.object_key, query_url=url, interpretation="lookup_failed", error=str(exc)))
    pd.DataFrame(rows).to_csv(results / "reference_comparison.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    check(parser.parse_args().results)
