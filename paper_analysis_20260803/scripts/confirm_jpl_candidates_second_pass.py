#!/usr/bin/env python3
"""Numerically confirm first-pass JPL candidate associations."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from query_jpl_identifications import API_DOC, API_URL, safe_float, sexagesimal_dec, sexagesimal_ra


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def angle_difference(first: float, second: float) -> float:
    return abs((first - second + 180.0) % 360.0 - 180.0)


def query(row: dict, args) -> tuple[dict, dict]:
    params = {
        "mpc-code": args.mpc_code,
        "obs-time": f"{safe_float(row['median_mjd']) + 2400000.5:.8f}",
        "fov-ra-center": sexagesimal_ra(safe_float(row["median_ra_deg"])),
        "fov-dec-center": sexagesimal_dec(safe_float(row["median_dec_deg"])),
        "fov-ra-hwidth": f"{args.fov_halfwidth_deg:.5f}",
        "fov-dec-hwidth": f"{args.fov_halfwidth_deg:.5f}",
        "vmag-lim": f"{args.vmag_limit:.2f}",
        "mag-required": "true",
        "two-pass": "true",
        "suppress-first-pass": "true",
    }
    url = API_URL + "?" + urllib.parse.urlencode(params)
    started = time.monotonic()
    request = urllib.request.Request(url, headers={"User-Agent": "SHARP-PASP-frozen-analysis/1.0"})
    with urllib.request.urlopen(request, timeout=args.timeout_s) as response:
        payload = json.loads(response.read().decode("utf-8"))
    return {
        "night": row["night"],
        "trk_sub": row["trk_sub"],
        "linkage_id": int(row["linkage_id"]),
        "query_url": url,
        "elapsed_s": time.monotonic() - started,
    }, payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reconciliation", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mpc-code", default="327")
    parser.add_argument("--fov-halfwidth-deg", type=float, default=0.03)
    parser.add_argument("--vmag-limit", type=float, default=22.5)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout-s", type=float, default=180.0)
    args = parser.parse_args()

    input_path = args.reconciliation.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    raw_dir = output / "raw"
    raw_dir.mkdir()
    frame = pd.read_csv(input_path, dtype={"night": "string", "trk_sub": "string"})
    strict = frame.loc[frame["strict_candidate_association"].astype(str).str.lower().isin({"true", "1"})].copy()
    if strict.empty:
        raise ValueError("no strict first-pass candidates")

    metadata = []
    payloads = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        future_map = {executor.submit(query, row, args): row for row in strict.to_dict("records")}
        for index, future in enumerate(concurrent.futures.as_completed(future_map), 1):
            meta, payload = future.result()
            metadata.append(meta)
            payloads.append((future_map[future], payload))
            raw_path = raw_dir / f"{meta['night']}_{meta['trk_sub']}_{meta['linkage_id']}.json"
            raw_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            print(f"[jpl-second] {index}/{len(strict)} {meta['night']} {meta['trk_sub']}", flush=True)

    rows = []
    for source, payload in payloads:
        fields = payload.get("fields_second", [])
        records = [dict(zip(fields, values)) for values in payload.get("data_second_pass", [])]
        for record in records:
            separation = safe_float(record.get('Dist. from center Norm (")'))
            ra_rate = safe_float(record.get('RA rate ("/h)'))
            dec_rate = safe_float(record.get('Dec rate ("/h)'))
            speed = float(np.hypot(ra_rate, dec_rate))
            direction = float(np.mod(np.rad2deg(np.arctan2(dec_rate, ra_rate)), 360.0))
            speed_difference = abs(speed - safe_float(source["candidate_speed_arcsec_per_hour"]))
            direction_difference = angle_difference(direction, safe_float(source["candidate_direction_deg"]))
            object_name = record.get("Object name", "")
            rows.append(
                {
                    "night": source["night"],
                    "trk_sub": source["trk_sub"],
                    "linkage_id": int(source["linkage_id"]),
                    "first_pass_object_name": source["jpl_object_name"],
                    "second_pass_object_name": object_name,
                    "second_pass_separation_arcsec": separation,
                    "second_pass_ra_rate_arcsec_per_hour": ra_rate,
                    "second_pass_dec_rate_arcsec_per_hour": dec_rate,
                    "second_pass_speed_arcsec_per_hour": speed,
                    "candidate_speed_arcsec_per_hour": safe_float(source["candidate_speed_arcsec_per_hour"]),
                    "speed_difference_arcsec_per_hour": speed_difference,
                    "direction_difference_deg": direction_difference,
                    "object_name_agrees": object_name == source["jpl_object_name"],
                    "numerically_confirmed_candidate": bool(
                        object_name == source["jpl_object_name"]
                        and separation <= 1.0
                        and speed_difference <= 1.0
                        and direction_difference <= 5.0
                    ),
                    "api_version": payload.get("signature", {}).get("version", ""),
                }
            )
    result = pd.DataFrame(rows).sort_values(["night", "trk_sub", "second_pass_separation_arcsec"])
    result.to_csv(output / "jpl_second_pass_confirmations.csv", index=False)
    pd.DataFrame(metadata).sort_values(["night", "trk_sub"]).to_csv(output / "jpl_second_pass_query_status.csv", index=False)
    confirmed_links = result.loc[result["numerically_confirmed_candidate"], ["night", "trk_sub", "linkage_id"]].drop_duplicates()
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "api": API_URL,
        "api_documentation": API_DOC,
        "first_pass_strict_linkages": int(len(strict)),
        "second_pass_result_rows": int(len(result)),
        "numerically_confirmed_linkages": int(len(confirmed_links)),
        "confirmed_object_names": sorted(result.loc[result["numerically_confirmed_candidate"], "second_pass_object_name"].unique().tolist()),
        "guardrail": "numerically confirmed JPL candidate association; MPC ingest/designation status still requires submission records",
        "input_sha256": sha256(input_path),
    }
    (output / "jpl_second_pass_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (output / "hashes.sha256").open("w", encoding="utf-8") as handle:
        handle.write(f"{sha256(input_path)}  {input_path}\n")
        for path in sorted(raw_dir.glob("*.json")):
            handle.write(f"{sha256(path)}  {path.relative_to(output)}\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
