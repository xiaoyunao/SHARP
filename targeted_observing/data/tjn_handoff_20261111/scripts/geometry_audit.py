#!/usr/bin/env python3
"""Audit the proposed SHARP TJN fixed grid using nominal, undistorted TAN WCS.

This is an independent analysis script, not production scheduler code.  It uses
only Python's standard library and NumPy.  It DOES NOT use a measured image WCS,
fit TPV distortion, or establish an operationally safe crop/pointing tolerance.

Reproduce from this directory:
    python3 geometry_audit.py

The default input is ../sharp_source/targeted_observing/data/mp_footprints.json.
Output is geometry_audit.json beside this script; both paths can be overridden.

Proposed grid (angles in degrees):
    dec_j = j * 2.67496
    N_j = ceil(360 * cos(dec_j) / 2.67264)
    ra_jk = k * 360 / N_j,  k = 0, ..., N_j - 1

Each row owns a stripe dec_j +/- 2.67496/2.  Each center owns a cell spanning
half the RA row pitch to either side.  To test the nominal geometrical coverage,
we project that cell's densely sampled boundary to its center's TAN plane.
Positive x is east, positive y is north; a north-up detector may reverse x.

For d=Dec, d0=center Dec, a=wrapped RA difference, all in radians:
    D = sin(d0)*sin(d) + cos(d0)*cos(d)*cos(a)
    x = cos(d)*sin(a)/D
    y = (cos(d0)*sin(d) - sin(d0)*cos(d)*cos(a))/D

The TAN plane coordinates x,y are dimensionless radians.  Multiplication by
180/pi gives intermediate projection-plane degrees, to compare with N_pix *
nominal arcsec_per_pix / 3600.  Do NOT replace x with atan(x) for that comparison.
Margins are distances from the worst cell point to nominal detector edges in
this projection plane; they must jointly accommodate clipping, pointing, PA,
distortion and other operational effects.  The comparison to 0.15 deg means a
literal 0.15-deg inward crop on EACH detector edge, not a total overlap width.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


S_X_DEG = 2.67264
S_DEC_DEG = 2.67496
NX, NY = 9216, 9232
DEC_MIN, DEC_MAX = -26.0, 57.0
EDGE_SAMPLES = 2049
INSET_COMPARISON_DEG = 0.15


def wrap_deg(angle):
    return (np.asarray(angle) + 180.0) % 360.0 - 180.0


def audit_mp(path: Path) -> dict:
    raw = path.read_bytes()
    source = json.loads(raw)
    rows = source["rows"]
    ra = np.array([r["center_ra"] for r in rows], dtype=float)
    dec = np.array([r["center_dec"] for r in rows], dtype=float)
    width = np.abs(wrap_deg(np.array([r["corner_ra_2"] for r in rows]) -
                            np.array([r["corner_ra_1"] for r in rows])))
    width *= np.cos(np.deg2rad(dec))
    height = np.abs(np.array([r["corner_dec_3"] for r in rows]) -
                    np.array([r["corner_dec_1"] for r in rows]))
    r, d = np.deg2rad(ra), np.deg2rad(dec)
    xyz = np.column_stack((np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)))
    dots = xyz @ xyz.T
    np.fill_diagonal(dots, -1.0)
    nearest = np.rad2deg(np.arccos(np.clip(dots.max(axis=1), -1.0, 1.0)))
    wx, hy = float(np.median(width)), float(np.median(height))
    return {
        "local_path": str(path.resolve()),
        "local_sha256": hashlib.sha256(raw).hexdigest(),
        "upstream_source_path": source.get("source_path"),
        "upstream_source_sha256": source.get("source_sha256"),
        "upstream_retrieved_at": source.get("retrieved_at"),
        "field_count": len(rows),
        "dec_min_deg": float(dec.min()),
        "dec_max_deg": float(dec.max()),
        "distinct_dec_count": int(len(np.unique(dec))),
        "rectangle_width_cosdec_median_deg": wx,
        "rectangle_height_median_deg": hy,
        "inferred_nominal_scale_x_arcsec_per_pix": wx*3600/NX,
        "inferred_nominal_scale_y_arcsec_per_pix": hy*3600/NY,
        "equivalent_width_pixels_at_1p16_arcsec": wx*3600/1.16,
        "equivalent_height_pixels_at_1p16_arcsec": hy*3600/1.16,
        "height_difference_from_9232_pixels_at_1p16_arcsec": NY-hy*3600/1.16,
        "ninety_percent_width_deg": wx*0.9,
        "ninety_percent_height_deg": hy*0.9,
        "nearest_center_separation_min_deg": float(nearest.min()),
        "nearest_center_separation_median_deg": float(np.median(nearest)),
        "nearest_center_separation_max_deg": float(nearest.max()),
        "caution": "MP corners are RA-Dec rectangle approximations; not measured WCS. RA differences were wrapped before measuring widths. At 1.16 arcsec/pixel the stored height implies 9224 pixels, eight fewer than the assumed 9232-pixel nominal model; the reason for this discrepancy is not inferred.",
    }


def boundary_projection(dec_deg: float, ra_pitch_deg: float) -> dict:
    hdec, hra = S_DEC_DEG/2, ra_pitch_deg/2
    ras = np.linspace(-hra, hra, EDGE_SAMPLES)
    decs = np.linspace(dec_deg-hdec, dec_deg+hdec, EDGE_SAMPLES)
    da = np.deg2rad(np.concatenate((ras, ras, np.full_like(decs, -hra),
                                    np.full_like(decs, hra))))
    d = np.deg2rad(np.concatenate((np.full_like(ras, dec_deg-hdec),
                                  np.full_like(ras, dec_deg+hdec), decs, decs)))
    d0 = math.radians(dec_deg)
    denominator = np.sin(d0)*np.sin(d) + np.cos(d0)*np.cos(d)*np.cos(da)
    if not np.all(denominator > 0):
        raise ValueError("Cell crosses the TAN projection horizon")
    x = np.rad2deg(np.cos(d)*np.sin(da)/denominator)
    y = np.rad2deg((np.cos(d0)*np.sin(d) - np.sin(d0)*np.cos(d)*np.cos(da))/denominator)
    ix, iy = int(np.argmax(np.abs(x))), int(np.argmax(np.abs(y)))
    # Corners occur at the first/last entry of each edge. Check our dense-boundary
    # extrema against those corners, rather than assuming a rectangular sky map.
    ci = np.array([0, EDGE_SAMPLES-1, EDGE_SAMPLES, 2*EDGE_SAMPLES-1])
    x_extreme, y_extreme = float(np.max(np.abs(x))), float(np.max(np.abs(y)))
    return {
        "maximum_abs_tan_x_deg": x_extreme,
        "maximum_abs_tan_y_deg": y_extreme,
        "minimum_denominator": float(denominator.min()),
        "x_extreme_sky_offset": {"ra_offset_deg": float(np.rad2deg(da[ix])),
                                 "dec_deg": float(np.rad2deg(d[ix]))},
        "y_extreme_sky_offset": {"ra_offset_deg": float(np.rad2deg(da[iy])),
                                 "dec_deg": float(np.rad2deg(d[iy]))},
        "boundary_extrema_equal_corner_extrema": bool(
            np.isclose(x_extreme, np.max(np.abs(x[ci])), rtol=0, atol=1e-12) and
            np.isclose(y_extreme, np.max(np.abs(y[ci])), rtol=0, atol=1e-12)),
    }


def audit_row(dec_deg: float, row_index: int | None = None) -> dict:
    n_ra = math.ceil(360*math.cos(math.radians(dec_deg))/S_X_DEG)
    pitch = 360.0/n_ra
    projected = boundary_projection(dec_deg, pitch)
    separation = math.degrees(2*math.asin(math.cos(math.radians(dec_deg))*
                                         math.sin(math.radians(pitch)/2)))
    models = {}
    for scale in (1.15, 1.16):
        width, height = NX*scale/3600, NY*scale/3600
        mx = width/2 - projected["maximum_abs_tan_x_deg"]
        my = height/2 - projected["maximum_abs_tan_y_deg"]
        models[f"nominal_{scale:.2f}_arcsec_per_pix"] = {
            "scale_arcsec_per_pix": scale,
            "nominal_width_tan_deg": width,
            "nominal_height_tan_deg": height,
            "minimum_x_edge_margin_deg": mx,
            "minimum_y_edge_margin_deg": my,
            "minimum_any_edge_margin_deg": min(mx, my),
            "minimum_x_edge_margin_arcmin": mx*60,
            "minimum_y_edge_margin_arcmin": my*60,
            "minimum_x_edge_margin_pixels": mx*3600/scale,
            "minimum_y_edge_margin_pixels": my*3600/scale,
            "nominal_uncropped_cell_covered": bool(mx >= 0 and my >= 0),
            "cell_covered_after_0p15_deg_per_edge_crop": bool(
                mx >= INSET_COMPARISON_DEG and my >= INSET_COMPARISON_DEG),
            "minimum_remaining_margin_after_0p15_deg_crop_deg":
                min(mx, my) - INSET_COMPARISON_DEG,
        }
    return {
        "row_index": row_index,
        "center_dec_deg": dec_deg,
        "cell_dec_min_deg": dec_deg-S_DEC_DEG/2,
        "cell_dec_max_deg": dec_deg+S_DEC_DEG/2,
        "n_ra_full_circle": n_ra,
        "ra_pitch_deg": pitch,
        "local_parallel_pitch_deg": pitch*math.cos(math.radians(dec_deg)),
        "great_circle_center_separation_deg": separation,
        "boundary_projection": projected,
        "nominal_models": models,
    }


def summarize(rows: list[dict]) -> dict:
    summary = {}
    for name in ("nominal_1.15_arcsec_per_pix", "nominal_1.16_arcsec_per_pix"):
        worst = min(rows, key=lambda row: row["nominal_models"][name]["minimum_any_edge_margin_deg"])
        m = worst["nominal_models"][name]
        summary[name] = {
            "row_count": len(rows),
            "worst_row_index": worst["row_index"],
            "worst_center_dec_deg": worst["center_dec_deg"],
            "minimum_any_edge_margin_deg": m["minimum_any_edge_margin_deg"],
            "minimum_any_edge_margin_arcmin": m["minimum_any_edge_margin_deg"]*60,
            "minimum_any_edge_margin_pixels": m["minimum_any_edge_margin_deg"]*3600/m["scale_arcsec_per_pix"],
            "all_nominal_uncropped_cells_covered": all(
                row["nominal_models"][name]["nominal_uncropped_cell_covered"] for row in rows),
            "rows_passing_0p15_deg_per_edge_crop": sum(
                row["nominal_models"][name]["cell_covered_after_0p15_deg_per_edge_crop"] for row in rows),
        }
    return summary


def main() -> None:
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mp", type=Path, default=base.parent/"sharp_source/targeted_observing/data/mp_footprints.json")
    parser.add_argument("--output", type=Path, default=base/"geometry_audit.json")
    args = parser.parse_args()
    j_min, j_max = math.ceil(DEC_MIN/S_DEC_DEG), math.floor(DEC_MAX/S_DEC_DEG)
    rows = [audit_row(j*S_DEC_DEG, j) for j in range(j_min, j_max+1)]
    cover_min = math.ceil((DEC_MIN-S_DEC_DEG/2)/S_DEC_DEG)
    cover_max = math.floor((DEC_MAX+S_DEC_DEG/2)/S_DEC_DEG)
    cover_rows = [audit_row(j*S_DEC_DEG, j) for j in range(cover_min, cover_max+1)]
    data = {
        "audit_kind": "independent_nominal_TAN_geometry_only",
        "measured_image_WCS_used": False,
        "production_source_modified": False,
        "numpy_version": np.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "mp_source_audit": audit_mp(args.mp),
        "grid_definition": {
            "ra_anchor_deg": 0.0, "dec_anchor_deg": 0.0,
            "nominal_east_west_pitch_deg": S_X_DEG,
            "constant_dec_pitch_deg": S_DEC_DEG,
            "full_ring_formula": "N_j=ceil(360*cos(dec_j)/2.67264); ra_jk=k*360/N_j",
            "center_dec_formula": "dec_j=j*2.67496",
            "requested_center_dec_range_deg": [DEC_MIN, DEC_MAX],
            "center_in_range_row_indices": [j_min, j_max],
            "cell_intersects_requested_range_row_indices": [cover_min, cover_max],
            "orientation": "local north-up PA=0; detector x may be mirrored",
            "boundary_samples_per_edge": EDGE_SAMPLES,
            "detector_pixel_shape_xy": [NX, NY],
        },
        "summary_center_in_range": summarize(rows),
        "rows_center_in_range": rows,
        "summary_cell_intersects_range": summarize(cover_rows),
        "rows_cell_intersects_range": cover_rows,
        "hypothetical_rows_not_part_of_frozen_grid": [audit_row(54.26), audit_row(60.0)],
        "limitations": [
            "Nominal pixel scales 1.15 and 1.16 arcsec/pixel are alternative assumptions; neither is a measured WCS for the observing setup.",
            "No TPV/SIP distortion, bad pixels, chip/readout masks, pointing scatter, PA errors, vignetting or template common-area losses are included.",
            "Positive nominal margins are necessary geometric room, not a proven operational safety allowance.",
            "0.15 deg is tested as a single-sided crop on each of four detector edges; it cannot be adopted without reducing pitch or independently validating real footprints.",
            "Full-ring row counts describe the mother grid; no assumption is made that every field around a full RA ring is selected for TJN.",
            "Rows_center_in_range has centers inside [-26,57] deg; rows_cell_intersects_range also includes boundary rows needed to cover that full declination interval.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({"output": str(args.output.resolve()),
                      "summary_center_in_range": data["summary_center_in_range"],
                      "boundary_row_indices": [cover_min, cover_max]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
