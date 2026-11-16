"""Extend the production MP grid offline and certify footprint coverage."""

import argparse
import hashlib
import json
from pathlib import Path

import astropy.units as u
import numpy as np
from astroplan import Observer
from astropy.coordinates import GeocentricTrueEcliptic, SkyCoord
from astropy.table import Table, vstack
from astropy.time import Time
from astropy.utils import iers
from scipy.spatial import cKDTree

SOURCE_SHA = "13f394877ce26e50e5cee5053a4ad845f412be4520421e46931b3e157d9e0225"
WIDTH, HEIGHT = np.array([9216, 9224]) * 1.16 / 3600
OVERLAP = 0.10
FIRST_BETA = -20 + HEIGHT / 2
BETA_STEP = HEIGHT * (1 - OVERLAP)
LAST_LEGACY_RING = int((20 - HEIGHT / 2 - FIRST_BETA) / BETA_STEP)
POLAR_DEC = 75
ECLIPTIC = GeocentricTrueEcliptic(equinox="J2000")
SEARCH_RADIUS = 2 * np.sin(np.deg2rad(3.2) / 2)


def grid(rings):
    lon, lat = [], []
    for ring in rings:
        beta = FIRST_BETA + ring * BETA_STEP
        step = WIDTH * (1 - OVERLAP) / np.cos(np.deg2rad(beta))
        lambdas = np.arange((ring % 2) * step / 2, 360.001, step)
        lon.extend(lambdas % 360)
        lat.extend(np.full(len(lambdas), beta))
    return SkyCoord(lon=lon * u.deg, lat=lat * u.deg, frame=ECLIPTIC).icrs


def unit_vectors(ra, dec):
    lon, lat = np.deg2rad([ra, dec])
    return np.stack((np.cos(lon) * np.cos(lat), np.sin(lon) * np.cos(lat), np.sin(lat)), axis=-1)


def footprint_table(centers, first_id):
    ra, dec = centers.ra.deg, centers.dec.deg
    spans = WIDTH / (2 * np.cos(centers.dec.rad))
    table = Table()
    table["field_id"] = [f"{i:04d}" for i in range(first_id, first_id + len(ra))]
    table["center_ra"] = ra
    table["center_dec"] = dec
    polar = np.abs(dec) >= POLAR_DEC
    for corner, (east, north) in enumerate([(-1, -1), (1, -1), (1, 1), (-1, 1)], 1):
        corner_ra = (ra + east * spans) % 360
        corner_dec = dec + north * HEIGHT / 2
        # Invert the gnomonic projection where the RA/cos(Dec) approximation breaks down.
        x, y = np.tan(np.deg2rad([east * WIDTH / 2, north * HEIGHT / 2]))
        polar_centers = centers[polar]
        exact = polar_centers.directional_offset_by(
            np.arctan2(x, y) * u.rad, np.arctan(np.hypot(x, y)) * u.rad
        )
        corner_ra[polar], corner_dec[polar] = exact.ra.deg, exact.dec.deg
        table[f"corner_ra_{corner}"] = corner_ra
        table[f"corner_dec_{corner}"] = corner_dec
    return table


class Coverage:
    """Signed angular margin to the edge of the best nearby footprint."""

    def __init__(self, table):
        self.ra, self.dec = np.asarray(table["center_ra"]), np.asarray(table["center_dec"])
        self.center = unit_vectors(self.ra, self.dec)
        self.tree = cKDTree(self.center)
        corners = np.stack([
            unit_vectors(table[f"corner_ra_{i}"], table[f"corner_dec_{i}"])
            for i in range(1, 5)
        ], axis=1)
        normals = np.cross(corners, np.roll(corners, -1, axis=1))
        self.normals = normals / np.linalg.norm(normals, axis=-1, keepdims=True)
        self.polar = np.abs(self.dec) >= POLAR_DEC
        self.half_ra = (table["corner_ra_2"] - self.ra) % 360

    def margin(self, ra, dec, positions, fields):
        # The ordinary footprints have two constant-Dec edges and two meridians.
        dra = (ra[positions] - self.ra[fields] + 180) % 360 - 180
        ddec = dec[positions] - self.dec[fields]
        side = np.sin(np.deg2rad(self.half_ra[fields] - np.abs(dra)))
        meridian_margin = np.rad2deg(np.arcsin(side * np.cos(np.deg2rad(dec[positions]))))
        margin = np.minimum(meridian_margin, HEIGHT / 2 - np.abs(ddec))
        polar = self.polar[fields]
        # The polar TAN rectangles are intersections of four great-circle half-spaces.
        vectors = unit_vectors(ra[positions[polar]], dec[positions[polar]])
        dot = np.einsum("nki,ni->nk", self.normals[fields[polar]], vectors)
        margin[polar] = np.rad2deg(np.arcsin(np.clip(dot.min(axis=1), -1, 1)))
        return margin

    def best_margin(self, ra, dec):
        distances, indices = self.tree.query(
            unit_vectors(ra, dec), k=24, distance_upper_bound=SEARCH_RADIUS
        )
        if np.any(np.isfinite(distances[:, -1])):
            raise RuntimeError("Coverage search found too many nearby fields")
        positions, neighbors = np.nonzero(np.isfinite(distances))
        margin = np.full(len(ra), -np.inf)
        np.maximum.at(margin, positions, self.margin(ra, dec, positions, indices[positions, neighbors]))
        return margin


def uncovered_cells(table):
    coverage = Coverage(table)
    ra, dec = np.meshgrid(np.arange(0.5, 360, 1), np.arange(-89.5, 90, 1))
    ra, dec = ra.ravel(), dec.ravel()
    half = 0.5
    for level in range(6):
        margin = coverage.best_margin(ra, dec)
        # Travel along the cell's center parallel, then along a meridian.
        # No point in the cell can be farther from its center than this bound.
        radius = half * (1 + np.cos(np.deg2rad(dec)))
        uncertain = margin < radius + 1e-9
        ra, dec, margin = ra[uncertain], dec[uncertain], margin[uncertain]
        if not len(ra) or level == 5:
            return ra, dec, margin, half * 2
        half /= 2
        ra = (ra[:, None] + np.array([-half, half, -half, half])).ravel()
        dec = (dec[:, None] + np.array([-half, -half, half, half])).ravel()


def gap_centers(ra, dec, margin, cell_size):
    """Place supplemental fields where they cover the most residual cells."""
    vectors = unit_vectors(ra, dec)
    tree = cKDTree(vectors)
    covered = np.zeros(len(ra), dtype=bool)
    chosen = []
    offsets = np.deg2rad(np.array([-1.2, 0, 1.2]))
    east, north = np.meshgrid(offsets, offsets)
    east, north = east.ravel(), north.ravel()
    separation = np.arctan(np.hypot(np.tan(east), np.tan(north))) * u.rad
    angle = np.arctan2(np.tan(east), np.tan(north)) * u.rad
    for i in np.argsort(margin, kind="stable"):
        if covered[i]:
            continue
        seed = SkyCoord(ra[i], dec[i], unit="deg")
        candidates = seed.directional_offset_by(angle, separation)
        model = Coverage(footprint_table(candidates, 1))
        seed_margin = model.margin(ra, dec, np.full(len(candidates), i), np.arange(len(candidates)))
        usable = np.flatnonzero(seed_margin >= cell_size + 1e-9)
        neighbors = tree.query_ball_point(model.center[usable], SEARCH_RADIUS)
        best_cells = np.array([], dtype=int)
        best_field = None
        for field, nearby in zip(usable, neighbors):
            cells = np.asarray(nearby, dtype=int)
            cells = cells[~covered[cells]]
            margins = model.margin(ra, dec, cells, np.full(len(cells), field))
            cells = cells[margins >= cell_size + 1e-9]
            if len(cells) > len(best_cells):
                best_field, best_cells = field, cells
        if best_field is None:
            raise RuntimeError("Could not fill a residual sky cell")
        chosen.append((candidates.dec.deg[best_field], candidates.ra.deg[best_field]))
        covered[best_cells] = True
    coords = np.array(sorted(chosen)).reshape(-1, 2)
    return SkyCoord(ra=coords[:, 1] * u.deg, dec=coords[:, 0] * u.deg)


def visible_days(centers, earth_orientation):
    observer = Observer(latitude=40.394 * u.deg, longitude=117.575 * u.deg, elevation=960 * u.m)
    days = Time("2025-01-01") + np.arange(365) * u.day
    with iers.earth_orientation_table.set(earth_orientation):
        dusk = observer.twilight_evening_astronomical(days, which="next")
        dawn = observer.twilight_morning_astronomical(days + 0.5 * u.day, which="next")
        length = dawn - dusk
        times = dusk[:, None] + length[:, None] * np.linspace(0, 1, 50)
        # Cache the common, time-dependent transform and bound the broadcast size.
        frame = observer.altaz(times)
        result = []
        for start in range(0, len(centers), 64):
            sky = centers[start : start + 64, None, None].transform_to(frame)
            above = np.count_nonzero(sky.alt.deg > 30, axis=-1)
            result.extend(np.count_nonzero(above / 50 * length.to_value(u.hour) >= 1, axis=-1))
    return np.asarray(result, dtype=np.int64)


def plot_catalog(table, original, grid_count, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={"projection": "mollweide"})
    groups = [
        (slice(original, original + grid_count), "#078889", 3, "Extended grid"),
        (slice(original + grid_count, None), "#e78124", 9, "Gap fills"),
        (slice(0, original), "#354c9a", 3, "Original MP fields"),
    ]
    for rows, color, size, label in groups:
        ra = -((np.asarray(table["center_ra"][rows]) + 180) % 360 - 180)
        dec = table["center_dec"][rows]
        ax.scatter(np.deg2rad(ra), np.deg2rad(dec), color=color, s=size, linewidths=0, label=label)
    ax.grid(alpha=0.3)
    ticks = ax.get_xticks()
    ax.set_xticks(ticks, labels=[f"{round(-x) % 360}°" for x in np.rad2deg(ticks)])
    ax.set_title("SHARP all-sky field centers (ICRS)", pad=18)
    ax.set_xlabel("Right ascension (increases to the left)")
    ax.set_ylabel("Declination")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], loc="lower center", ncol=3, frameon=False)
    fig.subplots_adjust(bottom=0.15)
    fig.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Original production MP footprint FITS")
    parser.add_argument("output", type=Path, help="Directory for the offline catalog and coverage audit")
    parser.add_argument("--iers", type=Path, required=True, help="IERS-A finals2000A file covering 2025")
    args = parser.parse_args()
    earth_orientation = iers.IERS_A.open(str(args.iers))

    source_sha = hashlib.sha256(args.source.read_bytes()).hexdigest()
    if source_sha != SOURCE_SHA:
        raise ValueError("Input must be the verified production footprint table")
    legacy = Table.read(args.source)
    original = len(legacy)
    reconstructed = footprint_table(grid(range(LAST_LEGACY_RING + 1)), 1)
    np.testing.assert_array_equal(legacy["field_id"], reconstructed["field_id"])
    for name in reconstructed.colnames[1:]:
        np.testing.assert_allclose(legacy[name], reconstructed[name], rtol=0, atol=1e-10)

    low = int(np.ceil((-90 - FIRST_BETA) / BETA_STEP))
    high = int(np.floor((90 - FIRST_BETA) / BETA_STEP))
    rings = [*range(low, 0), *range(LAST_LEGACY_RING + 1, high + 1)]
    extended = grid(rings)
    polar_caps = SkyCoord(lon=[0, 0] * u.deg, lat=[-90, 90] * u.deg, frame=ECLIPTIC).icrs
    grid_centers = SkyCoord(
        ra=np.concatenate([extended.ra.deg, polar_caps.ra.deg]) * u.deg,
        dec=np.concatenate([extended.dec.deg, polar_caps.dec.deg]) * u.deg,
    )
    grid_table = footprint_table(grid_centers, original + 1)
    geometry_columns = reconstructed.colnames
    gaps_ra, gaps_dec, margin, cell_size = uncovered_cells(vstack([legacy[geometry_columns], grid_table]))
    gaps = gap_centers(gaps_ra, gaps_dec, margin, cell_size)
    gap_table = footprint_table(gaps, original + len(grid_table) + 1)
    added = vstack([grid_table, gap_table], join_type="exact")
    print(f"Fields: {original} original, {len(grid_table)} extension, {len(gap_table)} gap fills", flush=True)

    if original + len(added) > 9999:
        raise RuntimeError("Catalog exceeds the production four-digit field ID format")
    added_coords = SkyCoord(added["center_ra"], added["center_dec"], unit="deg")
    added["visible_day"] = visible_days(added_coords, earth_orientation)
    table = vstack([legacy, added], join_type="exact")
    table.meta["SRC_SHA"] = source_sha
    table.meta["NOLD"] = original
    table.meta["NGRID"] = len(grid_table)
    table.meta["NGAP"] = len(gap_table)
    table.meta["VISYEAR"] = 2025

    args.output.mkdir(parents=True, exist_ok=True)
    fits_path = args.output / "survey_fov_footprints_allsky_with_visibility.fits"
    table.write(fits_path, overwrite=True)
    result = Table.read(fits_path)
    for name in legacy.colnames:
        np.testing.assert_array_equal(result[name][:original], legacy[name])
    np.testing.assert_array_equal(result["field_id"], [f"{i:04d}" for i in range(1, len(result) + 1)])
    min_distance, _ = cKDTree(unit_vectors(result["center_ra"], result["center_dec"])).query(
        unit_vectors(result["center_ra"], result["center_dec"]), k=2
    )
    if np.any(min_distance[:, 1] < 1e-8):
        raise RuntimeError("Duplicate field centers")
    ra, _, _, resolution = uncovered_cells(result)
    if len(ra):
        raise RuntimeError(f"{len(ra)} sky cells could not be certified")

    legacy_visibility = visible_days(
        SkyCoord(legacy["center_ra"], legacy["center_dec"], unit="deg"), earth_orientation
    )
    visibility_diff = legacy_visibility - legacy["visible_day"]
    report = {
        "source_sha256": source_sha,
        "iers_sha256": hashlib.sha256(args.iers.read_bytes()).hexdigest(),
        "catalog_sha256": hashlib.sha256(fits_path.read_bytes()).hexdigest(),
        "original_fields": original,
        "extended_ring_fields": len(extended),
        "ecliptic_pole_fields": len(polar_caps),
        "gap_fill_fields": len(gaps),
        "total_fields": len(result),
        "new_first_id": str(added["field_id"][0]),
        "new_last_id": str(added["field_id"][-1]),
        "finest_coverage_cell_deg": resolution,
        "uncertified_cells": len(ra),
        "pre_fill_uncertified_cells": len(gaps_ra),
        "pre_fill_uncovered_cell_centers": int(np.count_nonzero(margin < 0)),
        "new_fields_with_visible_days": int(np.count_nonzero(added["visible_day"])),
        "all_fields_with_visible_days": int(np.count_nonzero(result["visible_day"])),
        "legacy_visibility_recalculation_differing_fields": int(np.count_nonzero(visibility_diff)),
        "legacy_visibility_recalculation_max_difference_days": int(np.abs(visibility_diff).max()),
    }
    (args.output / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    plot_catalog(result, original, len(grid_table), args.output / "mp_allsky.png")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
