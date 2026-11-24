#!/usr/bin/env python3
"""Join independently streamed date shards, rejecting conflicting duplicate rows."""
import argparse
import csv
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("shards", nargs="+", help="path:first_directory_night:last_directory_night")
    args = parser.parse_args()
    rows = {}
    fields = None
    for shard in args.shards:
        name, first, last = shard.rsplit(":", 2)
        with gzip.open(name, "rt", newline="") as file:
            reader = csv.DictReader(file)
            if fields is not None and reader.fieldnames != fields:
                raise ValueError(f"Header schema differs in {name}")
            fields = reader.fieldnames
            try:
                for row in reader:
                    night = row["directory_night"]
                    if night > last:
                        break
                    if night < first:
                        continue
                    if None in row or any(val is None for val in row.values()):
                        raise ValueError(f"Incomplete source row in {name}: {night}")
                    key = (night, row["kind"], row["source_file"])
                    if key in rows and rows[key] != row:
                        raise ValueError(f"Inconsistent source row: {key}")
                    rows[key] = row
            except EOFError:
                raise ValueError(f"Stream ended before requested end {last}: {name}") from None
    path = ROOT / "inputs" / "server_headers.csv.gz"
    with gzip.open(path, "wt", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows[key] for key in sorted(rows))
    print(f"Combined {len(rows)} source rows")


if __name__ == "__main__":
    main()
