#!/usr/bin/env python3
"""Freeze JPL SB Identification queries for the retained 58 linkages.

This is a current-epoch cross-identification aid, not an MPC ingest-status
lookup.  Any positional/rate agreement is labelled a candidate association
until the authors reconcile it with MPC submission records.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


API_URL = "https://ssd-api.jpl.nasa.gov/sb_ident.api"
API_DOC = "https://ssd-api.jpl.nasa.gov/doc/sb_ident.html"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sexagesimal_ra(ra_deg: float) -> str:
    total_seconds = (ra_deg % 360.0) / 15.0 * 3600.0
    hours = int(total_seconds // 3600)
    minutes = int((total_seconds - hours * 3600) // 60)
    seconds = total_seconds - hours * 3600 - minutes * 60
    if seconds >= 59.995:
        seconds = 0.0
        minutes += 1
    if minutes >= 60:
        minutes = 0
        hours = (hours + 1) % 24
    return f"{hours:02d}-{minutes:02d}-{seconds:05.2f}"


def sexagesimal_dec(dec_deg: float) -> str:
    prefix = "M" if dec_deg < 0 else ""
    value = abs(dec_deg)
    degrees = int(value)
    minutes_float = (value - degrees) * 60.0
    minutes = int(minutes_float)
    seconds = (minutes_float - minutes) * 60.0
    if seconds >= 59.995:
        seconds = 0.0
        minutes += 1
    if minutes >= 60:
        minutes = 0
        degrees += 1
    return f"{prefix}{degrees:02d}-{minutes:02d}-{seconds:05.2f}"


def safe_float(value) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else np.nan
    except (TypeError, ValueError):
        return np.nan


def angular_difference(first: float, second: float) -> float:
    return abs((first - second + 180.0) % 360.0 - 180.0)


def fetch_one(row: dict, args) -> tuple[dict, dict | None]:
    params = {
        "mpc-code": args.mpc_code,
        "obs-time": f"{safe_float(row['median_mjd']):.8f}",
        "fov-ra-center": sexagesimal_ra(safe_float(row["median_ra_deg"])),
        "fov-dec-center": sexagesimal_dec(safe_float(row["median_dec_deg"])),
        "fov-ra-hwidth": f"{args.fov_halfwidth_deg:.5f}",
        "fov-dec-hwidth": f"{args.fov_halfwidth_deg:.5f}",
        "vmag-lim": f"{args.vmag_limit:.2f}",
        "mag-required": "true",
        "two-pass": "false",
    }
    # SB Identification accepts JD, whereas the frozen links store MJD.
    params["obs-time"] = f"{safe_float(row['median_mjd']) + 2400000.5:.8f}"
    url = API_URL + "?" + urllib.parse.urlencode(params)
    error = ""
    started = time.monotonic()
    for attempt in range(1, args.retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "SHARP-PASP-frozen-analysis/1.0"})
            with urllib.request.urlopen(request, timeout=args.timeout_s) as response:
                payload = json.loads(response.read().decode("utf-8"))
            elapsed = time.monotonic() - started
            return (
                {
                    "night": row["night"],
                    "trk_sub": row["trk_sub"],
                    "linkage_id": int(row["linkage_id"]),
                    "status": "ok",
                    "attempts": attempt,
                    "elapsed_s": elapsed,
                    "url": url,
                    "api_version": payload.get("signature", {}).get("version", ""),
                    "returned_first_pass_n": int(payload.get("n_first_pass", 0)),
                    "error": "",
                },
                payload,
            )
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
            error = f"{type(exc).__name__}: {exc}"
            if attempt < args.retries:
                time.sleep(min(2 ** attempt, 8))
    return (
        {
            "night": row["night"],
            "trk_sub": row["trk_sub"],
            "linkage_id": int(row["linkage_id"]),
            "status": "error",
            "attempts": args.retries,
            "elapsed_s": time.monotonic() - started,
            "url": url,
            "api_version": "",
            "returned_first_pass_n": 0,
            "error": error,
        },
        None,
    )


def parse_candidates(row: dict, payload: dict, maximum_separation_arcsec: float) -> list[dict]:
    fields = payload.get("fields_first", [])
    data = payload.get("data_first_pass", [])
    candidates = []
    candidate_speed = safe_float(row["speed_arcsec_per_hour"])
    candidate_direction = safe_float(row["lin_dir_deg"])
    for values in data:
        record = dict(zip(fields, values))
        separation = safe_float(record.get('Dist. from center Norm (")'))
        if not np.isfinite(separation) or separation > maximum_separation_arcsec:
            continue
        ra_rate = safe_float(record.get('RA rate ("/h)'))
        dec_rate = safe_float(record.get('Dec rate ("/h)'))
        jpl_speed = float(np.hypot(ra_rate, dec_rate)) if np.isfinite(ra_rate) and np.isfinite(dec_rate) else np.nan
        jpl_direction = float(np.mod(np.rad2deg(np.arctan2(dec_rate, ra_rate)), 360.0)) if np.isfinite(jpl_speed) else np.nan
        speed_difference = abs(jpl_speed - candidate_speed) if np.isfinite(jpl_speed) and np.isfinite(candidate_speed) else np.nan
        direction_difference = angular_difference(jpl_direction, candidate_direction) if np.isfinite(jpl_direction) and np.isfinite(candidate_direction) else np.nan
        strict = bool(
            separation <= 5.0
            and np.isfinite(speed_difference)
            and speed_difference <= 2.0
            and np.isfinite(direction_difference)
            and direction_difference <= 10.0
        )
        plausible = bool(
            separation <= 30.0
            and np.isfinite(speed_difference)
            and speed_difference <= 5.0
            and np.isfinite(direction_difference)
            and direction_difference <= 20.0
        )
        score_terms = [separation / 5.0]
        if np.isfinite(speed_difference):
            score_terms.append(speed_difference / 2.0)
        if np.isfinite(direction_difference):
            score_terms.append(direction_difference / 10.0)
        candidates.append(
            {
                "night": row["night"],
                "trk_sub": row["trk_sub"],
                "linkage_id": int(row["linkage_id"]),
                "jpl_object_name": record.get("Object name", ""),
                "jpl_ra": record.get("Astrometric RA (hh:mm:ss)", ""),
                "jpl_dec": record.get('Astrometric Dec (dd mm\'ss")', ""),
                "separation_arcsec": separation,
                "jpl_visual_magnitude": safe_float(record.get("Visual magnitude (V)")),
                "jpl_ra_rate_arcsec_per_hour": ra_rate,
                "jpl_dec_rate_arcsec_per_hour": dec_rate,
                "jpl_speed_arcsec_per_hour": jpl_speed,
                "jpl_direction_deg": jpl_direction,
                "candidate_speed_arcsec_per_hour": candidate_speed,
                "candidate_direction_deg": candidate_direction,
                "speed_difference_arcsec_per_hour": speed_difference,
                "direction_difference_deg": direction_difference,
                "estimated_error_ra_arcsec": safe_float(record.get('Est. error RA (")')),
                "estimated_error_dec_arcsec": safe_float(record.get('Est. error Dec (")')),
                "joint_score": float(np.hypot.reduce(score_terms)),
                "strict_candidate_association": strict,
                "plausible_candidate_association": plausible,
            }
        )
    return sorted(candidates, key=lambda item: item["joint_score"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--links", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mpc-code", default="327")
    parser.add_argument("--fov-halfwidth-deg", type=float, default=0.03)
    parser.add_argument("--vmag-limit", type=float, default=22.5)
    parser.add_argument("--candidate-separation-arcsec", type=float, default=180.0)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout-s", type=float, default=90.0)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()

    links_path = args.links.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    raw_dir = output / "raw"
    raw_dir.mkdir()
    links = pd.read_csv(links_path, dtype={"night": "string", "trk_sub": "string"})
    links["night"] = links["night"].str.zfill(8)
    if len(links) != 58:
        raise ValueError(f"expected 58 retained linkages, found {len(links)}")
    required = ["night", "trk_sub", "linkage_id", "median_mjd", "median_ra_deg", "median_dec_deg", "speed_arcsec_per_hour", "lin_dir_deg"]
    missing = [column for column in required if column not in links]
    if missing:
        raise ValueError(f"missing columns: {missing}")

    rows = links[required].to_dict("records")
    statuses = []
    candidate_rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {executor.submit(fetch_one, row, args): row for row in rows}
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            row = futures[future]
            status, payload = future.result()
            statuses.append(status)
            if payload is not None:
                raw_path = raw_dir / f"{row['night']}_{row['trk_sub']}_{int(row['linkage_id'])}.json"
                raw_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                candidate_rows.extend(parse_candidates(row, payload, args.candidate_separation_arcsec))
            print(f"[jpl] {index}/{len(rows)} {row['night']} {row['trk_sub']} status={status['status']}", flush=True)

    status_frame = pd.DataFrame(statuses).sort_values(["night", "trk_sub"])
    candidate_frame = pd.DataFrame(candidate_rows)
    if candidate_frame.empty:
        candidate_frame = pd.DataFrame(columns=["night", "trk_sub", "linkage_id", "jpl_object_name", "separation_arcsec", "joint_score", "strict_candidate_association", "plausible_candidate_association"])
    else:
        candidate_frame = candidate_frame.sort_values(["night", "trk_sub", "joint_score"])
    best = (
        candidate_frame.groupby(["night", "trk_sub", "linkage_id"], as_index=False).first()
        if not candidate_frame.empty
        else candidate_frame.copy()
    )
    reconciliation = links.merge(best, on=["night", "trk_sub", "linkage_id"], how="left", suffixes=("", "_jpl"))
    reconciliation["jpl_query_state"] = np.where(
        reconciliation["jpl_object_name"].notna(), "candidate_rows_within_search_radius", "no_candidate_within_search_radius"
    )
    reconciliation["identification_authority_state"] = "unverified_requires_mpc_submission_reconciliation"
    status_frame.to_csv(output / "jpl_query_status.csv", index=False)
    candidate_frame.to_csv(output / "jpl_candidate_matches.csv", index=False)
    best.to_csv(output / "jpl_best_match_by_link.csv", index=False)
    reconciliation.to_csv(output / "mpc_reconciliation.csv", index=False)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "api": API_URL,
        "api_documentation": API_DOC,
        "mpc_observatory_code": args.mpc_code,
        "query_linkages": int(len(links)),
        "successful_queries": int(status_frame["status"].eq("ok").sum()),
        "failed_queries": int(status_frame["status"].eq("error").sum()),
        "candidate_rows_within_180_arcsec": int(len(candidate_frame)),
        "links_with_any_candidate_within_180_arcsec": int(candidate_frame[["night", "trk_sub", "linkage_id"]].drop_duplicates().shape[0]),
        "strict_candidate_association_links": int(candidate_frame.loc[candidate_frame.get("strict_candidate_association", False).fillna(False), ["night", "trk_sub", "linkage_id"]].drop_duplicates().shape[0]) if "strict_candidate_association" in candidate_frame else 0,
        "plausible_candidate_association_links": int(candidate_frame.loc[candidate_frame.get("plausible_candidate_association", False).fillna(False), ["night", "trk_sub", "linkage_id"]].drop_duplicates().shape[0]) if "plausible_candidate_association" in candidate_frame else 0,
        "guardrail": "JPL positional/rate candidates are not MPC ingest or designation states; row-level MPC records remain author-supplied",
        "input_sha256": sha256(links_path),
    }
    (output / "jpl_identification_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (output / "hashes.sha256").open("w", encoding="utf-8") as handle:
        handle.write(f"{sha256(links_path)}  {links_path}\n")
        for path in sorted(raw_dir.glob("*.json")):
            handle.write(f"{sha256(path)}  {path.relative_to(output)}\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
