#!/usr/bin/env python3
"""Collect raw MP headers over a frozen filename inventory, reusing snapshots.

Only primary FITS header blocks are read on the server. Image data are not read.
Run locally; output and resume data stay in the supplied analysis directory.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess


REMOTE = r'''
from pathlib import Path
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip, hashlib, json, sys
tasks = json.loads(PAYLOAD)
def read(row):
    row = dict(row)
    path = Path(row['path'])
    try:
        opener = gzip.open if path.name.endswith('.gz') else open
        blocks = []
        with opener(path, 'rb') as stream:
            for _ in range(128):
                block = stream.read(2880)
                if len(block) != 2880:
                    raise ValueError('Truncated primary FITS header')
                blocks.append(block)
                if any(block[i:i+8] == b'END     ' for i in range(0,2880,80)):
                    break
            else:
                raise ValueError('No END card within 128 header blocks')
        raw = b''.join(blocks)
        header = fits.Header.fromstring(raw.decode('ascii'), sep='')
        row['header'] = {k:header[k] for k in ['DATE','OBS_DATE','EXPOSURE','OBJ_RA','OBJ_DEC','TEL_RA','TEL_DEC','OBJECT','READMODE','CHANNEL','TELESCOP'] if k in header}
        row['header_sha256'] = hashlib.sha256(raw).hexdigest()
    except Exception as exc:
        row['error'] = f'{type(exc).__name__}: {exc}'
    return row
with ThreadPoolExecutor(max_workers=WORKERS) as pool:
    futures = [pool.submit(read, row) for row in tasks]
    for i, future in enumerate(as_completed(futures)):
        print(json.dumps(future.result()), flush=True)
        if (i+1)%1000 == 0:
            print('HEADERS',i+1,'/',len(tasks),file=sys.stderr,flush=True)
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    base = args.source
    out = base / "expanded"
    inventory = json.loads((out / "archive_inventory.json").read_text())
    cache, reused_sources = {}, []
    for path in sorted(base.glob("headers_*.json")):
        snapshot = json.loads(path.read_text())
        reused_sources.append({"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        for row in snapshot["rows"]:
            if "header" in row:
                cache[(snapshot["night"], row["name"])] = row
    target = out / "archive_mp_headers.jsonl"
    if target.exists():
        for line in target.read_text().splitlines():
            try:
                row = json.loads(line)
                if "header" in row:
                    cache[(row["night"], row["name"])] = row
            except json.JSONDecodeError:
                pass
    tasks, reused = [], []
    for night in inventory["rows"]:
        parsed = []
        for name in night["files"]:
            match = re.search(r"_(\d+)\.fits(?:\.gz)?$", name)
            if match:
                parsed.append((int(match[1]), name))
        parsed.sort()
        counts = Counter(sequence for sequence, _ in parsed)
        for i, (sequence, name) in enumerate(parsed):
            if not re.fullmatch(r"OBJ_MP_\d+_\d+\.fits(?:\.gz)?", name):
                continue
            row = {"night": night["night"], "name": name, "sequence": sequence,
                   "path": str(Path(night["directory"]) / name),
                   "inventory_previous_name": parsed[i-1][1] if i else None,
                   "inventory_sequence_unique": counts[sequence] == 1}
            key = (night["night"], name)
            if key in cache:
                cached = cache[key]
                row.update({field: cached[field] for field in ("header", "header_sha256") if field in cached})
                reused.append(row)
            else:
                tasks.append(row)
    metadata = {"expected_mp_headers": len(tasks) + len(reused), "reused_headers": len(reused),
                "requested_headers": len(tasks), "workers": args.workers,
                "inventory_sha256": hashlib.sha256((out / "archive_inventory.json").read_bytes()).hexdigest(),
                "reused_sources": reused_sources,
                "method": "Read 2880-byte primary-header blocks through END; parse with astropy.Header; no image data"}
    (out / "archive_collection.json").write_text(json.dumps(metadata, indent=2))
    print(json.dumps(metadata | {"reused_sources": len(reused_sources)}), flush=True)
    with target.open("wb") as stream:
        for row in reused:
            stream.write((json.dumps(row) + "\n").encode())
        stream.flush()
        code = REMOTE.replace("PAYLOAD", repr(json.dumps(tasks))).replace("WORKERS", str(args.workers))
        subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "-p", "20093",
                        "smtpipeline@www.xinglong-naoc.cn", "/home/smtpipeline/Softwares/miniconda3/bin/python -"],
                       input=code.encode(), stdout=stream, check=True)
    print("COMPLETE", len(tasks) + len(reused), flush=True)


if __name__ == "__main__":
    main()
