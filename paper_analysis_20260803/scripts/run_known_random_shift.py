#!/usr/bin/env python3
"""Run a deterministic random-position shift test for known-object matches.

The script samples source catalogs across every available night, repeats the
production nearest-neighbour geometry after shifting all predicted coordinates
by 30 or 60 arcsec in a deterministic random direction, and writes only to an
isolated output directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree


def iter_nights(start: str, end: str):
    value = datetime.strptime(start, "%Y-%m-%d").date()
    finish = datetime.strptime(end, "%Y-%m-%d").date()
    while value <= finish:
        yield value.strftime("%Y%m%d")
        value += timedelta(days=1)


def text(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            item.decode("utf-8", "replace").strip()
            if isinstance(item, (bytes, bytearray))
            else str(item).strip()
            for item in values
        ]
    )


def unit_vectors(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    cos_dec = np.cos(dec)
    return np.column_stack((cos_dec * np.cos(ra), cos_dec * np.sin(ra), np.sin(dec)))


def nearest_separation_arcsec(
    predicted_ra: np.ndarray,
    predicted_dec: np.ndarray,
    detected_ra: np.ndarray,
    detected_dec: np.ndarray,
) -> np.ndarray:
    tree = cKDTree(unit_vectors(detected_ra, detected_dec))
    chord, _ = tree.query(unit_vectors(predicted_ra, predicted_dec), k=1, workers=1)
    radians = 2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))
    return np.rad2deg(radians) * 3600.0


def shifted_coordinates(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    distance_arcsec: float,
    position_angle_rad: float,
) -> tuple[np.ndarray, np.ndarray]:
    distance_deg = distance_arcsec / 3600.0
    dec_shift = distance_deg * math.cos(position_angle_rad)
    cos_dec = np.cos(np.deg2rad(dec_deg))
    ra_shift = distance_deg * math.sin(position_angle_rad) / np.clip(cos_dec, 0.1, None)
    return (ra_deg + ra_shift) % 360.0, np.clip(dec_deg + dec_shift, -89.999, 89.999)


def deterministic_angle(night: str, source_file: str, replicate: int) -> float:
    digest = hashlib.sha256(f"{night}|{source_file}|{replicate}".encode()).digest()
    integer = int.from_bytes(digest[:8], "big")
    return 2.0 * math.pi * integer / 2**64


def select_sources(source_files: np.ndarray, per_night: int) -> list[str]:
    unique = sorted(set(text(source_files)))
    ranked = sorted(
        unique,
        key=lambda name: hashlib.sha256(name.encode()).hexdigest(),
    )
    return ranked[:per_night]


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--processed-root", type=Path, default=Path("/processed1"))
    parser.add_argument("--start", default="2025-11-15")
    parser.add_argument("--end", default="2026-07-15")
    parser.add_argument("--per-night", type=int, default=1)
    parser.add_argument("--replicates", type=int, default=4)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    processed = args.processed_root.resolve()
    if output == processed or processed in output.parents or output in processed.parents:
        raise ValueError("output must not overlap processed root")
    output.mkdir(parents=True, exist_ok=False)

    rows: list[dict] = []
    errors: list[dict] = []
    sampled_nights = 0
    for night in iter_nights(args.start, args.end):
        all_path = processed / night / "L4" / f"{night}_all_asteroids.fits"
        if not all_path.exists():
            continue
        with fits.open(all_path, memmap=True) as hdul:
            data = hdul[1].data
            sources = select_sources(data["source_file"], args.per_night)
            source_values = text(data["source_file"])
            for source_file in sources:
                selection = source_values == source_file
                predicted_ra = np.asarray(data["ra"][selection], dtype=float)
                predicted_dec = np.asarray(data["dec"][selection], dtype=float)
                predicted_mag = np.asarray(data["mag"][selection], dtype=float)
                predicted_ok = np.isfinite(predicted_ra) & np.isfinite(predicted_dec) & np.isfinite(predicted_mag)
                predicted_ra = predicted_ra[predicted_ok]
                predicted_dec = predicted_dec[predicted_ok]
                l2_path = processed / night / "L2" / source_file
                try:
                    with fits.open(l2_path, memmap=True) as l2_hdul:
                        table = l2_hdul[1].data
                        names = set(l2_hdul[1].columns.names)
                        if not {"RA_Win", "DEC_Win", "Mag_Kron"}.issubset(names):
                            raise KeyError("RA_Win/DEC_Win/Mag_Kron schema gate failed")
                        detected_ra = np.asarray(table["RA_Win"], dtype=float)
                        detected_dec = np.asarray(table["DEC_Win"], dtype=float)
                        detected_mag = np.asarray(table["Mag_Kron"], dtype=float)
                    detected_ok = np.isfinite(detected_ra) & np.isfinite(detected_dec) & np.isfinite(detected_mag)
                    detected_ra = detected_ra[detected_ok]
                    detected_dec = detected_dec[detected_ok]
                    if not len(predicted_ra) or not len(detected_ra):
                        raise ValueError("empty finite prediction or detection table")
                    baseline = nearest_separation_arcsec(
                        predicted_ra, predicted_dec, detected_ra, detected_dec
                    )
                    rows.append(
                        {
                            "night": night,
                            "source_file": source_file,
                            "test": "unshifted_baseline",
                            "shift_arcsec": 0.0,
                            "replicate": 0,
                            "position_angle_deg": np.nan,
                            "predicted_n": len(predicted_ra),
                            "detected_finite_mag_n": len(detected_ra),
                            "match_lt_1arcsec_n": int((baseline < 1.0).sum()),
                            "match_fraction": float((baseline < 1.0).mean()),
                            "median_nearest_arcsec": float(np.median(baseline)),
                        }
                    )
                    for distance in (30.0, 60.0):
                        for replicate in range(args.replicates):
                            angle = deterministic_angle(night, source_file, replicate + int(distance))
                            shifted_ra, shifted_dec = shifted_coordinates(
                                predicted_ra, predicted_dec, distance, angle
                            )
                            separation = nearest_separation_arcsec(
                                shifted_ra, shifted_dec, detected_ra, detected_dec
                            )
                            rows.append(
                                {
                                    "night": night,
                                    "source_file": source_file,
                                    "test": "random_position_shift",
                                    "shift_arcsec": distance,
                                    "replicate": replicate,
                                    "position_angle_deg": math.degrees(angle),
                                    "predicted_n": len(predicted_ra),
                                    "detected_finite_mag_n": len(detected_ra),
                                    "match_lt_1arcsec_n": int((separation < 1.0).sum()),
                                    "match_fraction": float((separation < 1.0).mean()),
                                    "median_nearest_arcsec": float(np.median(separation)),
                                }
                            )
                except Exception as exc:
                    errors.append(
                        {
                            "night": night,
                            "source_file": source_file,
                            "error": f"{type(exc).__name__}: {exc}",
                        }
                    )
        sampled_nights += 1
        if sampled_nights % 10 == 0:
            print(f"[shift] nights={sampled_nights} rows={len(rows)} errors={len(errors)}", flush=True)

    write_csv(output / "known_random_shift_test.csv", rows)
    write_csv(output / "known_random_shift_errors.csv", errors)
    shift_rows = [row for row in rows if row["test"] == "random_position_shift"]
    total_pred = sum(row["predicted_n"] for row in shift_rows)
    total_match = sum(row["match_lt_1arcsec_n"] for row in shift_rows)
    baseline_rows = [row for row in rows if row["test"] == "unshifted_baseline"]
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sampled_nights": sampled_nights,
        "sampled_exposures": len(baseline_rows),
        "random_shift_trials": len(shift_rows),
        "random_shift_prediction_trials": total_pred,
        "random_shift_matches_lt_1arcsec": total_match,
        "random_shift_chance_fraction": total_match / total_pred if total_pred else None,
        "baseline_prediction_n": sum(row["predicted_n"] for row in baseline_rows),
        "baseline_match_n": sum(row["match_lt_1arcsec_n"] for row in baseline_rows),
        "errors": len(errors),
        "method": "deterministic exposure-stratified 30/60 arcsec random-direction position shifts; finite Mag_Kron schema gate; nearest neighbor; strict <1 arcsec",
    }
    (output / "known_random_shift_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
