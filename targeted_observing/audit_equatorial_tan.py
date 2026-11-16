"""Map gaps of the equatorial survey grid with a representative camera TAN footprint."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
from scipy.spatial import cKDTree

from .build_equatorial_allsky import MIN_ALTITUDE, XINGLONG_LATITUDE
from .extend_mp_allsky import unit_vectors

FINEST_CELL = 1 / 256


def camera_geometry(evidence):
    spans = []
    roll = []
    for row in json.loads(evidence.read_text()):
        wcs = WCS(row["full_wcs_header"])
        ny, nx = row["detector_shape"]
        center = (np.array([nx, ny]) + 1) / 2
        offsets = np.array([[0, -ny / 2], [0, ny / 2], [-nx / 2, 0], [nx / 2, 0]])
        sky = SkyCoord(*wcs.all_pix2world(center + offsets, 1).T, unit="deg")
        spans.append([sky[0].separation(sky[1]).deg, sky[2].separation(sky[3]).deg])
        # Pixel +y is approximately east and pixel +x approximately north.
        east, north = wcs.pixel_scale_matrix.T[::-1]
        roll.append(np.mean([
            np.rad2deg(np.arctan2(east[0], east[1])) - 90,
            np.rad2deg(np.arctan2(north[0], north[1])),
        ]))
    return np.asarray(spans), np.asarray(roll)


class TanCoverage:
    def __init__(self, table, width, height, roll_deg):
        ra, dec = table["center_ra"], table["center_dec"]
        center = unit_vectors(ra, dec)
        east = unit_vectors(ra + 90, np.zeros(len(ra)))
        north = np.cross(center, east)
        angle = np.deg2rad(roll_deg)
        east, north = np.cos(angle) * east - np.sin(angle) * north, np.sin(angle) * east + np.cos(angle) * north
        x, y = np.tan(np.deg2rad([width / 2, height / 2]))
        self.normals = np.stack(
            [x * center + east, x * center - east, y * center + north, y * center - north], axis=1
        )
        self.normals /= np.linalg.norm(self.normals, axis=2, keepdims=True)
        self.tree = cKDTree(center)
        # Any field intersecting a 1-degree cell must be within this distance of its center.
        max_angle = np.arctan(np.hypot(x, y)) + np.deg2rad(np.sqrt(2) / 2)
        self.search_radius = 2 * np.sin(max_angle / 2) + 1e-10

    def margin(self, ra, dec):
        vectors = unit_vectors(ra, dec)
        distances, fields = self.tree.query(vectors, k=16, distance_upper_bound=self.search_radius)
        if np.any(np.isfinite(distances[:, -1])):
            raise RuntimeError("TAN coverage search reached the neighbor limit")
        pos, col = np.nonzero(np.isfinite(distances))
        dot = np.einsum("nki,ni->nk", self.normals[fields[pos, col]], vectors[pos])
        margins = np.full(len(ra), -np.inf)
        np.maximum.at(margins, pos, np.rad2deg(np.arcsin(np.clip(dot.min(axis=1), -1, 1))))
        return margins


def gap_cells(coverage):
    ra, dec = np.meshgrid(np.arange(0.5, 360, 1), np.arange(-89.5, 90, 1))
    ra, dec, half = ra.ravel(), dec.ravel(), 0.5
    gaps = []
    boundary_area = 0.0
    while len(ra):
        margins = coverage.margin(ra, dec)
        # Integrate the metric bound along the straight path in RA/Dec coordinates.
        cos_lat = np.cos(np.deg2rad(np.maximum(0, np.abs(dec) - half)))
        radius = half * np.sqrt(1 + cos_lat ** 2) + 1e-9
        missing = margins < -radius
        uncertain = np.abs(margins) <= radius
        if half * 2 <= FINEST_CELL:
            boundary_area = float(cell_areas(dec[uncertain], half).sum())
            missing |= uncertain & (margins < 0)
        gaps.append(np.column_stack([
            ra[missing], dec[missing], np.full(missing.sum(), half),
            (margins[missing] < -radius[missing]).astype(float),
        ]))
        if half * 2 <= FINEST_CELL:
            break
        half /= 2
        ra = (ra[uncertain, None] + np.array([-half, half, -half, half])).ravel()
        dec = (dec[uncertain, None] + np.array([-half, -half, half, half])).ravel()
    return np.vstack(gaps), boundary_area


def cell_areas(dec, half):
    return np.rad2deg(1.0) * (2 * half) * (
        np.sin(np.deg2rad(dec + half)) - np.sin(np.deg2rad(dec - half))
    )


def merge_cells_for_plot(gaps):
    """Merge neighboring rectangles within each two-degree longitude block."""
    ra, dec, half = gaps[:, :3].T
    order = np.lexsort((ra, half, dec))
    ra, dec, half = ra[order], dec[order], half[order]
    breaks = np.flatnonzero(
        (np.diff(dec) != 0) | (np.diff(half) != 0)
        | (np.diff(ra) != half[:-1] + half[1:]) | (np.diff(ra // 2) != 0)
    ) + 1
    first, last = np.r_[0, breaks], np.r_[breaks - 1, len(ra) - 1]
    left, right = ra[first] - half[first], ra[last] + half[last]
    return (left + right) / 2, dec[first], (right - left) / 2, half[first]


def plot_gaps(gaps, poles, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, Patch

    red = "#bf251e"
    fig = plt.figure(figsize=(12, 10), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=(1.65, 1))
    sky = fig.add_subplot(grid[0, :], projection="mollweide")
    ra, dec, half_ra, half_dec = merge_cells_for_plot(gaps)
    lon = 180 - ra
    dx = half_ra[:, None] * np.array([-1, 1, 1, -1])
    dy = half_dec[:, None] * np.array([-1, -1, 1, 1])
    sky.add_collection(PolyCollection(
        np.deg2rad(np.stack((lon[:, None] + dx, dec[:, None] + dy), axis=-1)),
        facecolors=red, edgecolors="none", antialiaseds=False,
    ))
    cutoff = XINGLONG_LATITUDE - (90 - MIN_ALTITUDE)
    line = np.linspace(-np.pi, np.pi, 721)
    sky.plot(line, np.full_like(line, np.deg2rad(cutoff)), c="#536979", lw=1,
             linestyle="--", zorder=2)
    sky.grid(alpha=0.3)
    ticks = sky.get_xticks()
    sky.set_xticks(ticks, labels=[f"{round(180 - x) % 360}°" for x in np.rad2deg(ticks)])
    sky.set_title("TAN coverage gaps in the original 5,808-field grid", fontsize=17, pad=14)
    sky.set_xlabel("ICRS right ascension (increases to the left)")
    sky.legend(handles=[
        Patch(color=red, label="Uncovered sky"),
        Line2D([], [], c="#536979", lw=1, ls="--", label="Xinglong 18° culmination limit"),
    ], loc="lower center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False)

    for col, (sign, title) in enumerate([(1, "North celestial pole"), (-1, "South celestial pole")]):
        ax = fig.add_subplot(grid[1, col], aspect="equal")
        visible = sign * dec + half_dec > 80
        angle = np.deg2rad(ra[visible, None] - dx[visible])
        radius = 90 - sign * (dec[visible, None] + dy[visible])
        polygons = np.stack((radius * np.sin(angle), sign * radius * np.cos(angle)), axis=-1)
        ax.add_collection(PolyCollection(polygons, facecolors=red, edgecolors="none", antialiaseds=False))
        for r in range(2, 11, 2):
            ax.add_patch(Circle((0, 0), r, fc="none", ec="#b4bdc3", lw=0.6))
            ax.text(0.2, r - 0.2, f"{sign * (90 - r)}°", fontsize=8, color="#52626e", va="top")
        for angle_deg in [0, 90, 180, 270]:
            angle_rad = np.deg2rad(angle_deg)
            x, y = np.sin(angle_rad), sign * np.cos(angle_rad)
            ax.plot([0, x * 10], [0, y * 10], color="#b4bdc3", lw=0.55, zorder=0)
            ax.text(x * 10.85, y * 10.85, f"{angle_deg}°", fontsize=8,
                    ha="center", va="center", color="#52626e")
        if poles[(sign + 1) // 2] < 0:
            ax.plot(0, 0, marker="o", mfc="none", mec=red, ms=4.5, mew=0.8, zorder=3)
            ax.annotate("Uncovered pole", xy=(0, 0), xytext=(6, 7),
                        textcoords="offset points", color=red, fontsize=8)
        ax.set(xlim=(-12, 12), ylim=(-12, 12), title=title)
        ax.set_axis_off()

    fig.savefig(output, dpi=250, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("catalog", type=Path, help="Unfiltered equatorial footprint FITS")
    parser.add_argument("evidence", type=Path, help="Archived measured TPV WCS headers")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    table = Table.read(args.catalog)
    spans, rolls = camera_geometry(args.evidence)
    width, height = np.median(spans, axis=0)
    roll = np.median(rolls)
    coverage = TanCoverage(table, width, height, roll)
    gaps, boundary_area = gap_cells(coverage)
    ra, dec, half, certified = gaps.T
    areas = cell_areas(dec, half)
    cutoff = XINGLONG_LATITUDE - (90 - MIN_ALTITUDE)
    poles = coverage.margin(np.array([0, 0]), np.array([-90, 90]))
    report = {
        "catalog_sha256": hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
        "wcs_evidence_sha256": hashlib.sha256(args.evidence.read_bytes()).hexdigest(),
        "wcs_samples": len(spans),
        "wcs_span_range_east_north_deg": [spans.min(axis=0).tolist(), spans.max(axis=0).tolist()],
        "tan_width_east_deg": float(width),
        "tan_height_north_deg": float(height),
        "clockwise_roll_deg": float(roll),
        "pointings": len(table),
        "finest_cell_deg": FINEST_CELL,
        "gap_area_estimate_sq_deg": float(areas.sum()),
        "gap_area_lower_bound_sq_deg": float(areas[certified == 1].sum()),
        "gap_area_upper_bound_sq_deg": float(areas[certified == 1].sum() + boundary_area),
        "gap_area_above_xinglong_culmination_limit_sq_deg": float(areas[dec >= cutoff].sum()),
        "gap_declination_range_deg": [float(dec.min()), float(dec.max())],
        "gap_closest_to_equator_deg": float(dec[np.argmin(np.abs(dec))]),
        "gap_closest_to_equator_by_hemisphere_deg": {
            "south": float(dec[dec < 0].max()), "north": float(dec[dec > 0].min()),
        },
        "south_pole_signed_margin_arcmin": float(poles[0] * 60),
        "north_pole_signed_margin_arcmin": float(poles[1] * 60),
        "gap_area_by_hemisphere_sq_deg": {
            "south": float(areas[dec < 0].sum()), "north": float(areas[dec > 0].sum()),
        },
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "tan_gap_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    plot_gaps(gaps, poles, args.output / "equatorial_tan_gaps.png")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
