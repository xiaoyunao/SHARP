#!/usr/bin/env python3
"""Offline digest2 rehearsal on an immutable copy of unknown ADES PSV.

No production entry point imports this module. This program cannot submit to MPC.
The catalog must describe exactly the objects represented by the input PSV (use
the matching submit_masked catalog when replaying a historical submission).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

VERSION = "2.8.2"
THRESHOLD = 65.0
REQUIRED = {"trkSub", "mode", "stn", "obsTime", "ra", "dec", "rmsRA", "rmsDec", "mag", "band"}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def json_hash(value: object) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode())


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def read_psv(text: str) -> tuple[list[str], list[str], dict[str, list[dict[str, str]]]]:
    """Strict single-block parser; never silently discard a measurement."""
    metadata, columns = [], []
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        if line.startswith(("#", "!")):
            if columns:
                raise ValueError("Multiple ADES blocks are not supported in this rehearsal")
            metadata.append(line)
            continue
        values = [v.strip() for v in line.split("|")]
        if not columns:
            columns = values
            if len(set(columns)) != len(columns) or not REQUIRED.issubset(columns):
                raise ValueError("Missing or duplicate PSV columns")
            continue
        if len(values) != len(columns):
            raise ValueError(f"PSV line {number}: wrong column count")
        row = dict(zip(columns, values))
        key = row["trkSub"]
        if not re.fullmatch(r"[-A-Za-z0-9_]{1,7}", key):
            raise ValueError(f"PSV line {number}: invalid trkSub {key!r}")
        groups[key].append(row)
    if not columns or not metadata or metadata[0].replace(" ", "") != "#version=2017":
        raise ValueError("Expected ADES 2017 header and observation columns")
    return metadata, columns, dict(groups)


def utc_time(value: str) -> datetime:
    if not value.endswith("Z"):
        raise ValueError("obsTime must explicitly end in Z (UTC)")
    return datetime.fromisoformat(value[:-1] + "+00:00")


def check_rows(rows: list[dict[str, str]], expected: int, sites: set[str]) -> None:
    if len(rows) != expected or expected < 2:
        raise ValueError(f"Incomplete tracklet: {len(rows)} PSV rows, {expected} catalog observations")
    epochs = []
    for row in rows:
        if row.get("provID") or row.get("permID"):
            raise ValueError("Designated object requires identity routing, not NEW NEO routing")
        if any(row.get(k) for k in ("sys", "pos1", "pos2", "pos3")) or row["mode"] != "CCD":
            raise ValueError("Only ground-based CCD observations are supported")
        if row["stn"] not in sites:
            raise ValueError(f"Unknown/invalid observatory {row['stn']}")
        if row["band"] != "G":
            raise ValueError("This pipeline rehearsal expects its documented G photometric band")
        epoch = utc_time(row["obsTime"])
        epochs.append(epoch)
        numbers = {k: float(row[k]) for k in ("ra", "dec", "rmsRA", "rmsDec", "mag")}
        if not all(math.isfinite(v) for v in numbers.values()):
            raise ValueError("Non-finite astrometry/photometry")
        if not 0 <= numbers["ra"] < 360 or not -90 <= numbers["dec"] <= 90:
            raise ValueError("RA/Dec outside angular bounds")
        if numbers["rmsRA"] <= 0 or numbers["rmsDec"] <= 0:
            raise ValueError("Non-positive astrometric uncertainty")
    if len({r["stn"] for r in rows}) != 1:
        raise ValueError("Mixed observatories require separate tracklets")
    if len(set(epochs)) != len(epochs):
        raise ValueError("Duplicate epoch: resolve duplicate/conflicting detections before scoring")
    if (max(epochs) - min(epochs)).total_seconds() > 86400:
        raise ValueError("Only single-night tracklets are supported")


def read_catalog(text: str) -> dict[str, int]:
    catalog = json.loads(text)
    if not isinstance(catalog, list):
        raise ValueError("Expected a list of unknown catalog records")
    expected = {}
    for row in catalog:
        key = str(row.get("trk_sub") or row.get("trkSub") or "")
        if not key or key in expected:
            raise ValueError("Missing or duplicate catalog trkSub")
        count = int(row["n_obs"])
        # Catalog count alone must not hide an inconsistent member list.
        members = []
        for name in ("image_names", "objids"):
            value = row.get(name)
            if value is not None:
                values = value if isinstance(value, list) else str(value).split(";")
                if len(values) != count:
                    raise ValueError(f"{key}: inconsistent catalog {name}")
                members.append(values)
        if len(members) == 2 and len(set(zip(*members))) != count:
            raise ValueError(f"{key}: duplicate catalog detection")
        expected[key] = count
    return expected


def read_decisions(text: str) -> dict[str, int | None]:
    result = {}
    reader = csv.DictReader(text.splitlines())
    if not {"tracklet_id", "is_real"}.issubset(reader.fieldnames or []):
        raise ValueError("Review CSV requires tracklet_id,is_real")
    for row in reader:
        key, value = row["tracklet_id"].strip(), row["is_real"].strip()
        if not key or key in result or value not in ("", "0", "1"):
            raise ValueError(f"Invalid or duplicate review decision: {key!r}")
        result[key] = int(value) if value else None
    return result


def route(decision: int | None, score: dict) -> str:
    if decision == 0:
        return "rejected"
    if decision is None:
        return "pending_review"
    if score["status"] != "ok":
        return "blocked"
    return "neo candidate" if score["neo_noid"] >= THRESHOLD else "unclassified"


def ground_sites(path: Path) -> set[str]:
    sites = set()
    for line in path.read_text().splitlines():
        try:
            constants = [float(line[a:b]) for a, b in ((4, 13), (13, 21), (21, 30))]
            if all(math.isfinite(v) for v in constants) and re.fullmatch(r"[A-Za-z0-9]{3}", line[:3]):
                sites.add(line[:3])
        except ValueError:
            continue
    return sites


def prepare(args: argparse.Namespace) -> dict:
    from digest2 import Digest2
    from digest2.model import find_config_path, find_model_path
    from digest2.observation import parse_ades_psv

    if importlib.metadata.version("digest2") != VERSION:
        raise ValueError(f"Rehearsal requires digest2=={VERSION}")
    psv_bytes = args.psv.read_bytes()
    catalog_bytes = args.catalog.read_bytes()
    review_bytes = args.review_csv.read_bytes() if args.review_csv else None
    metadata, columns, groups = read_psv(psv_bytes.decode())
    expected = read_catalog(catalog_bytes.decode())
    if set(groups) - set(expected):
        raise ValueError("PSV contains objects absent from the supplied catalog")
    decisions = read_decisions(review_bytes.decode()) if review_bytes is not None else {}
    model_path = Path(args.model or find_model_path())
    config_path = Path(args.config or find_config_path())
    if "random" in [s.strip() for s in config_path.read_text().splitlines()]:
        raise ValueError("Config overrides repeatable mode with random")
    sites = ground_sites(args.obscodes)
    if "327" not in sites:
        raise ValueError("ObsCodes does not contain a valid ground station 327")
    # New directory only: never replace existing review, score or submission files.
    args.out_dir.mkdir(parents=True, exist_ok=False)
    inputs = args.out_dir / "inputs"
    inputs.mkdir()
    (inputs / "observations.psv").write_bytes(psv_bytes)
    (inputs / "catalog.json").write_bytes(catalog_bytes)
    if review_bytes is not None:
        (inputs / "review.csv").write_bytes(review_bytes)
    provenance = {}
    for name, path in (("model", model_path), ("config", config_path), ("obscodes", args.obscodes)):
        data = path.read_bytes()
        copy = inputs / path.name
        if copy.exists():
            raise ValueError("Colliding dependency snapshot filenames")
        copy.write_bytes(data)
        provenance[name] = {"source": str(path.resolve()), "sha256": sha256(data), "snapshot": str(copy.relative_to(args.out_dir))}
    provenance.update({"digest2": VERSION, "repeatable": True, "is_ades": True,
                       "no_threshold": False, "input_psv_sha256": sha256(psv_bytes),
                       "catalog_sha256": sha256(catalog_bytes),
                       "review_sha256": sha256(review_bytes) if review_bytes is not None else None})
    # The official parser applies the G -> V magnitude correction. Constructing
    # Observation(mag=G) directly would bypass it and change the resulting scores.
    parsed = parse_ades_psv(str(inputs / "observations.psv"))
    scored = {}
    with Digest2(model_path=str(args.out_dir / provenance["model"]["snapshot"]),
                 config_path=str(args.out_dir / provenance["config"]["snapshot"]),
                 obscodes_path=str(args.out_dir / provenance["obscodes"]["snapshot"]),
                 repeatable=True, no_threshold=False) as engine:
        for key in sorted(expected):
            rows = sorted(groups.get(key, []), key=lambda r: r["obsTime"])
            item = {"status": "error", "n_obs": len(rows), "expected_n_obs": expected[key],
                    "observations_sha256": json_hash(rows), "neo_noid": None, "neo_raw": None}
            try:
                check_rows(rows, expected[key], sites)
                obs = parsed.get(key, [])
                if len(obs) != len(rows):
                    raise ValueError("Official parser dropped observations")
                result = engine.classify_tracklet(obs, classes=["NEO"], is_ades=True)
                noid, raw = float(result.noid.NEO), float(result.raw.NEO)
                if not all(math.isfinite(v) and 0 <= v <= 100 for v in (noid, raw)):
                    raise ValueError("digest2 returned an invalid score")
                item.update(status="ok", neo_noid=noid, neo_raw=raw)
            except (ValueError, RuntimeError, KeyError, TypeError, OverflowError) as exc:
                item["error"] = str(exc)
            item["decision"] = decisions.get(key)
            item["route"] = route(item["decision"], item)
            scored[key] = item
    batches = []
    for obj_type, filename in (("neo candidate", "neo_candidate.psv"), ("unclassified", "unclassified.psv")):
        keys = sorted(k for k, v in scored.items() if v["route"] == obj_type)
        if not keys:
            continue
        rows = [r for k in keys for r in sorted(groups[k], key=lambda r: r["obsTime"])]
        payload = "\n".join(metadata + [" | ".join(columns)] + ["|".join(r[c] for c in columns) for r in rows]) + "\n"
        digest = sha256(payload.encode())
        (args.out_dir / filename).write_text(payload)
        batches.append({"id": json_hash({"obj_type": obj_type, "psv_sha256": digest}),
                        "obj_type": obj_type, "file": filename, "sha256": digest,
                        "trkSub": keys, "n_obs": len(rows)})
    report = {"schema_version": 1, "mode": "offline_rehearsal", "production_enabled": False,
              "created_utc": datetime.now(timezone.utc).isoformat(), "threshold": THRESHOLD,
              "provenance": provenance, "scores": scored, "batches": batches,
              "blocked": sorted(k for k, v in scored.items() if v["route"] == "blocked")}
    write_json(args.out_dir / "routing.json", report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psv", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--review-csv", type=Path)
    parser.add_argument("--obscodes", type=Path, required=True)
    parser.add_argument("--model", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True, help="Must not already exist")
    args = parser.parse_args()
    report = prepare(args)
    print(json.dumps({"output": str(args.out_dir), "objects": len(report["scores"]),
                      "batches": [(b["obj_type"], len(b["trkSub"])) for b in report["batches"]],
                      "blocked": report["blocked"]}, indent=2))
    if report["blocked"] or any(v["status"] != "ok" for v in report["scores"].values()):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
