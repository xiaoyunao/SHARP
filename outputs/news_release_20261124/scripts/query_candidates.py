#!/usr/bin/env python3
"""Read the MPC WAMO status for positively reviewed SHARP tracklet identifiers."""
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path

import requests

EVIDENCE = Path(__file__).resolve().parents[1] / "evidence"
URL = "https://data.minorplanetcenter.net/api/wamo"
QUERY_LIMIT = 70


def main():
    query = set()
    with gzip.open(EVIDENCE / "server_candidates.jsonl.gz", "rt") as handle:
        for raw in handle:
            night = json.loads(raw)
            for row in night["human"] or []:
                if row.get("is_real", "").strip() == "1":
                    query.add(f"{row['tracklet_id']} 327")
    query = sorted(query)
    bundle = {"url": URL, "started_utc": datetime.now(timezone.utc).isoformat(),
              "query": query, "found": [], "not_found": [], "malformed": []}
    with requests.Session() as session:
        session.trust_env = False
        for index in range(0, len(query), QUERY_LIMIT):
            response = session.get(URL, json=query[index:index+QUERY_LIMIT], timeout=45)
            response.raise_for_status()
            result = response.json()
            for key in ("found", "not_found", "malformed"):
                bundle[key].extend(result.get(key, []))
    bundle["completed_utc"] = datetime.now(timezone.utc).isoformat()
    (EVIDENCE / "wamo.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
