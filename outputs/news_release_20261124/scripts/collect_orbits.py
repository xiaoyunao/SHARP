#!/usr/bin/env python3
"""Read ASTORB on smtdata and stream compact elements and alias evidence."""
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from aleph.DataBase import get_astdata

path = Path("/pipeline/xiaoyunao/known_asteroid/astorb.dat")
before = path.stat()
digest = hashlib.sha256()
writer = csv.writer(sys.stdout)
writer.writerow(["number", "name", "a_au", "ecc", "inc_deg", "epoch"])
rows = 0
with path.open("rb") as handle:
    for raw in handle:
        digest.update(raw)
        line = raw.decode("ascii")
        rec = get_astdata(line, source="astorb")
        writer.writerow([rec["num"].strip(), rec["name"].strip(), rec["a"], rec["e"],
                         rec["incl"], rec["epoch"]])
        rows += 1
assert path.stat().st_mtime_ns == before.st_mtime_ns
print(json.dumps({"path": str(path), "sha256": digest.hexdigest(), "rows": rows,
                  "bytes": before.st_size,
                  "mtime_utc": datetime.fromtimestamp(before.st_mtime, timezone.utc).isoformat()}),
      file=sys.stderr)
