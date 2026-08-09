#!/usr/bin/env python3
"""Freeze field visits and the nside=512 inclusive HEALPix sky coverage."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import healpy as hp
import numpy as np
import pandas as pd
from astropy.table import Table


NSIDE = 512
INCLUSIVE = True
FACT = 4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def field_string(value) -> str:
    if isinstance(value, (bytes, bytearray)):
        value = value.decode("utf-8", "replace")
    return str(value).strip().zfill(4)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-manifest", type=Path, required=True)
    parser.add_argument("--footprints", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(
        args.raw_manifest,
        dtype={"night": "string", "field_id": "string", "sequence_id": "string"},
    )
    raw["night"] = raw["night"].str.zfill(8)
    raw["field_id"] = raw["field_id"].str.zfill(4)
    strict = raw["strict_standard_science"].astype(str).str.lower().isin({"true", "1"})
    raw = raw.loc[strict].copy()
    raw["exposure_s"] = pd.to_numeric(raw["exposure_s"], errors="coerce").fillna(30.0)
    footprints = Table.read(args.footprints)
    index = {field_string(value): i for i, value in enumerate(footprints["field_id"])}
    missing_fields = sorted(set(raw["field_id"]) - set(index))
    if missing_fields:
        raise ValueError(f"raw fields absent from deployed footprint: {missing_fields[:20]}")

    field_pixels: dict[str, np.ndarray] = {}
    for field_id in sorted(raw["field_id"].unique()):
        row = index[field_id]
        ra = np.array([footprints[f"corner_ra_{corner}"][row] for corner in range(1, 5)], dtype=float)
        dec = np.array([footprints[f"corner_dec_{corner}"][row] for corner in range(1, 5)], dtype=float)
        field_pixels[field_id] = hp.query_polygon(
            NSIDE,
            hp.ang2vec(ra, dec, lonlat=True),
            inclusive=INCLUSIVE,
            fact=FACT,
            nest=False,
        )

    coverage_map = np.zeros(hp.nside2npix(NSIDE), dtype=np.uint8)
    for pixels in field_pixels.values():
        coverage_map[pixels] = 1
    coverage_path = output / "healpix_coverage.fits"
    hp.write_map(
        coverage_path,
        coverage_map,
        nest=False,
        dtype=np.uint8,
        overwrite=True,
        column_names=["OBSERVED"],
        extra_header=[
            ("METHOD", "query_polygon"),
            ("INCLUSIV", True),
            ("FACT", FACT),
            ("FRAME", "ICRS"),
        ],
    )

    visits = (
        raw.groupby("field_id", as_index=False)
        .agg(
            exposure_n=("file_name", "size"),
            night_n=("night", "nunique"),
            first_night=("night", "min"),
            last_night=("night", "max"),
            open_shutter_s=("exposure_s", "sum"),
        )
        .sort_values(["exposure_n", "field_id"], ascending=[False, True])
    )
    center_rows = []
    for field_id in visits["field_id"]:
        row = index[field_id]
        center_rows.append(
            {
                "field_id": field_id,
                "center_ra_deg": float(footprints["center_ra"][row]),
                "center_dec_deg": float(footprints["center_dec"][row]),
                "nominal_pixel_n_nside512": int(len(field_pixels[field_id])),
            }
        )
    visits = visits.merge(pd.DataFrame(center_rows), on="field_id", validate="one_to_one")
    visits.to_csv(output / "field_visit_counts.csv", index=False)

    cumulative = np.zeros_like(coverage_map)
    cumulative_exposure_s = 0.0
    history_rows = []
    for night, frame in raw.groupby("night", sort=True):
        for field_id in frame["field_id"].unique():
            cumulative[field_pixels[field_id]] = 1
        cumulative_exposure_s += float(frame["exposure_s"].sum())
        history_rows.append(
            {
                "night": night,
                "raw_exposure_n": int(len(frame)),
                "unique_field_n": int(frame["field_id"].nunique()),
                "cumulative_exposure_n": int(raw.loc[raw["night"] <= night].shape[0]),
                "cumulative_open_shutter_hours": cumulative_exposure_s / 3600.0,
                "cumulative_unique_field_n": int(raw.loc[raw["night"] <= night, "field_id"].nunique()),
                "cumulative_healpix_pixel_n": int(cumulative.sum()),
                "cumulative_area_deg2": float(cumulative.sum() * hp.nside2pixarea(NSIDE, degrees=True)),
            }
        )
    pd.DataFrame(history_rows).to_csv(output / "coverage_by_night.csv", index=False)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "raw_exposure_n": int(len(raw)),
        "raw_night_n": int(raw["night"].nunique()),
        "observed_field_n": int(raw["field_id"].nunique()),
        "open_shutter_hours": float(raw["exposure_s"].sum() / 3600.0),
        "healpix_nside": NSIDE,
        "healpix_ordering": "RING",
        "healpix_frame": "ICRS",
        "query_polygon_inclusive": INCLUSIVE,
        "query_polygon_fact": FACT,
        "observed_pixel_n": int(coverage_map.sum()),
        "unique_area_deg2": float(coverage_map.sum() * hp.nside2pixarea(NSIDE, degrees=True)),
        "input_hashes": {
            "raw_manifest": sha256(args.raw_manifest),
            "footprints": sha256(args.footprints),
        },
        "output_hashes": {"healpix_coverage.fits": sha256(coverage_path)},
    }
    (output / "coverage_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
