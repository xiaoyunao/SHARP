#!/usr/bin/env python3
"""Archive JPL SBDB classifications and the displayed asteroid's identity."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

EVIDENCE = Path(__file__).resolve().parents[1] / "evidence"
ENDPOINT = "https://ssd-api.jpl.nasa.gov"


def get(path):
    return json.loads(subprocess.check_output(
        ["curl", "-fsS", "--retry", "2", "--max-time", "40", ENDPOINT + path]))


for name, group in [("neo", "sb-group=neo"), ("trojan", "sb-class=TJN")]:
    path = f"/sbdb_query.api?{group}&sb-kind=a&fields=pdes,full_name,class&limit=100000"
    data = get(path)
    assert int(data["count"]) == len(data["data"]), (name, data.keys())
    data["source"] = ENDPOINT + path
    data["checked_utc"] = datetime.now(timezone.utc).isoformat()
    (EVIDENCE / f"jpl_{name}.json").write_text(json.dumps(data) + "\n")
    print(name, len(data["data"]))
(EVIDENCE / "jpl_282917.json").write_text(json.dumps(get("/sbdb.api?sstr=282917")) + "\n")

# An historical provisional designation absent from the current ASTORB snapshot.
(EVIDENCE / "jpl_2021vt99.json").write_text(
    json.dumps(get("/sbdb.api?sstr=2021%20VT99")) + "\n")
