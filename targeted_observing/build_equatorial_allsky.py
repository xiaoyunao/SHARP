"""Build the unfiltered legacy equatorial grid and compare its geometry to the local FITS."""

import argparse
import hashlib
import json
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np

from .extend_mp_allsky import plot_catalog

WIDTH = 9216 * 1.16 / 3600
HEIGHT = 9232 * 1.16 / 3600
OVERLAP = 0.10
ECLIPTIC_LIMIT = 20
XINGLONG_LATITUDE = 40.393
MIN_ALTITUDE = 18.0
GEOMETRY_COLUMNS = ["center_ra", "center_dec"] + [
    f"corner_{axis}_{i}" for axis in ("ra", "dec") for i in range(1, 5)
]


def field_ids(count):
    return np.array([f"{i:04d}" for i in range(1, count + 1)], dtype="U4")


def build_grid():
    step = HEIGHT * (1 - OVERLAP)
    north_center = 90 - HEIGHT / 2
    dec_centers = np.arange(-north_center, north_center + 1e-6, step)
    legacy_bands = len(dec_centers)
    if dec_centers[-1] < north_center - 1e-9:
        dec_centers = np.append(dec_centers, north_center)

    rows = []
    for band, dec in enumerate(dec_centers):
        cos_dec = np.cos(np.deg2rad(dec))
        ra_step = WIDTH * (1 - OVERLAP) / cos_dec
        n_ra = int(np.ceil(360 / ra_step))
        ras = (np.arange(n_ra) + 0.5) * (360 / n_ra)
        half_width = WIDTH / (2 * cos_dec)
        bottom, top = dec - HEIGHT / 2, dec + HEIGHT / 2
        for index, ra in enumerate(ras):
            left, right = (ra - half_width) % 360, (ra + half_width) % 360
            rows.append((band, index, ra, dec, left, right, right, left, bottom, bottom, top, top))

    table = Table(rows=rows, names=["dec_band", "ra_band", *GEOMETRY_COLUMNS])
    table.add_column(field_ids(len(table)), index=0, name="field_id")
    centers = SkyCoord(table["center_ra"], table["center_dec"], frame="icrs", unit="deg")
    table["ecliptic_lat_deg"] = centers.barycentrictrueecliptic.lat.deg
    table.meta["EDGEGEOM"] = "RADEC"
    table.meta["FOVXDEG"] = WIDTH
    table.meta["FOVYDEG"] = HEIGHT
    table.meta["OVERLAP"] = OVERLAP
    table.meta["NBANDLEG"] = legacy_bands
    return table


def check_coverage(table):
    """Check strip overlap, including the periodic RA seam and both polar endpoints."""
    edges = []
    minimum_ra_overlap = np.inf
    for group in table.group_by("dec_band").groups:
        lower, upper = group["corner_dec_1"], group["corner_dec_3"]
        assert np.all(lower == lower[0]) and np.all(upper == upper[0])
        ra = np.asarray(group["center_ra"])
        west = (ra - group["corner_ra_1"]) % 360
        east = (group["corner_ra_2"] - ra) % 360
        gap = np.diff(np.r_[ra, ra[0] + 360])
        overlap = east + np.roll(west, -1) - gap
        projected_overlap = overlap * np.cos(np.deg2rad(group["center_dec"][0]))
        minimum_ra_overlap = min(minimum_ra_overlap, projected_overlap.min())
        edges.append((lower[0], upper[0]))
    bottom, top = np.asarray(edges).T
    minimum_dec_overlap = np.min(top[:-1] - bottom[1:])
    np.testing.assert_allclose([bottom[0], top[-1]], [-90, 90], atol=1e-10, rtol=0)
    assert minimum_ra_overlap > 0 and minimum_dec_overlap > 0
    return float(minimum_ra_overlap), float(minimum_dec_overlap)


def compare_legacy(table, legacy):
    in_band = np.abs(table["ecliptic_lat_deg"]) <= ECLIPTIC_LIMIT
    subset = table[in_band]
    for column in ("dec_band", "ra_band"):
        np.testing.assert_array_equal(subset[column], legacy[column])

    max_error = {}
    for column in GEOMETRY_COLUMNS:
        delta = np.asarray(subset[column]) - np.asarray(legacy[column])
        if "_ra" in column:
            delta = (delta + 180) % 360 - 180
        np.testing.assert_allclose(delta, 0, atol=1e-10, rtol=0)
        max_error[column] = float(np.max(np.abs(delta)))

    table.meta["NOLD"] = int(in_band.sum())
    return in_band, max_error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy", type=Path, help="Local equatorial_footprints.fits")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    table = build_grid()
    legacy = Table.read(args.legacy)
    in_band, max_error = compare_legacy(table, legacy)
    ra_overlap, dec_overlap = check_coverage(table)

    args.output.mkdir(parents=True, exist_ok=True)
    fits_path = args.output / "equatorial_footprints_allsky.fits"
    table.write(fits_path, overwrite=True)
    result = Table.read(fits_path)
    for column in table.colnames:
        np.testing.assert_array_equal(result[column], table[column])
    check_coverage(result)

    culmination = 90 - np.abs(result["center_dec"] - XINGLONG_LATITUDE)
    observable = culmination >= MIN_ALTITUDE
    retained = result[observable]
    retained["field_id"] = field_ids(len(retained))
    retained.meta["NOLD"] = int(np.count_nonzero(in_band & observable))
    retained.meta["SITELAT"] = XINGLONG_LATITUDE
    retained.meta["MINALT"] = MIN_ALTITUDE
    retained.meta["NREMOVED"] = int(np.count_nonzero(~observable))
    visible_path = args.output / "equatorial_footprints_xinglong18.fits"
    retained.write(visible_path, overwrite=True)
    readback = Table.read(visible_path)
    for column in retained.colnames:
        np.testing.assert_array_equal(readback[column], retained[column])

    report = {
        "legacy_file_sha256": hashlib.sha256(args.legacy.read_bytes()).hexdigest(),
        "catalog_sha256": hashlib.sha256(fits_path.read_bytes()).hexdigest(),
        "fields": len(result),
        "xinglong_latitude_deg": XINGLONG_LATITUDE,
        "minimum_culmination_altitude_deg": MIN_ALTITUDE,
        "xinglong_fields": len(retained),
        "excluded_fields": len(result) - len(retained),
        "southern_culmination_declination_limit_deg": XINGLONG_LATITUDE - (90 - MIN_ALTITUDE),
        "lowest_retained_center_dec_deg": float(retained["center_dec"].min()),
        "xinglong_catalog_sha256": hashlib.sha256(visible_path.read_bytes()).hexdigest(),
        "dec_bands": len(np.unique(result["dec_band"])),
        "northern_polar_band_fields": int(np.count_nonzero(result["dec_band"] == result.meta["NBANDLEG"])),
        "legacy_fields_matched": int(in_band.sum()),
        "legacy_fields_total": len(legacy),
        "legacy_filter_max_abs_ecliptic_lat_deg": ECLIPTIC_LIMIT,
        "legacy_geometry_max_absolute_error_deg": max_error,
        "minimum_adjacent_ra_overlap_at_band_center_deg": ra_overlap,
        "minimum_adjacent_dec_overlap_deg": dec_overlap,
        "coverage_model": "constant-declination edges and meridian edges",
    }
    (args.output / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    plot_catalog(result, args.output / "equatorial_footprints_allsky.png", original_mask=in_band)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
