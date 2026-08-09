#!/usr/bin/env python3
"""Analyze the frozen unknown-link population and twilight context."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import (
    AltAz,
    EarthLocation,
    GeocentricMeanEcliptic,
    GCRS,
    SkyCoord,
    get_sun,
)
from astropy.time import Time
from astropy.utils import iers


iers.conf.auto_download = False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_site(config: dict) -> EarthLocation:
    # This is an explicitly provisional analysis coordinate, not a claim about
    # the final surveyed station location.
    site = config["site_configuration_as_executed"]["scheduler"]
    return EarthLocation.from_geodetic(
        site["longitude_deg_east"] * u.deg,
        site["latitude_deg_north"] * u.deg,
        site["height_m"] * u.m,
    )


def finite_median(values: pd.Series) -> float:
    array = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    array = array[np.isfinite(array)]
    return float(np.median(array)) if array.size else float("nan")


def twilight_boundaries(night: str, location: EarthLocation) -> tuple[float, float]:
    start = Time(f"{night[:4]}-{night[4:6]}-{night[6:]}T00:00:00", scale="utc")
    grid_mjd = start.mjd + np.linspace(0, 1, 1441)
    grid = Time(grid_mjd, format="mjd", scale="utc")
    altitude = get_sun(grid).transform_to(AltAz(obstime=grid, location=location)).alt.deg
    value = altitude + 18.0
    crossings = []
    for index in np.where(np.signbit(value[:-1]) != np.signbit(value[1:]))[0]:
        x0, x1 = grid_mjd[index], grid_mjd[index + 1]
        y0, y1 = value[index], value[index + 1]
        root = x0 - y0 * (x1 - x0) / (y1 - y0)
        direction = "dawn" if value[index] < value[index + 1] else "dusk"
        crossings.append((direction, float(root)))
    dusk = next((mjd for kind, mjd in crossings if kind == "dusk"), np.nan)
    dawn = next((mjd for kind, mjd in crossings if kind == "dawn" and mjd > dusk), np.nan)
    return dusk, dawn


def add_twilight_metrics(frame: pd.DataFrame, location: EarthLocation) -> pd.DataFrame:
    output = frame.copy()
    boundaries = {night: twilight_boundaries(night, location) for night in sorted(output["night"].unique())}
    output["astronomical_dusk_mjd"] = output["night"].map(lambda night: boundaries[night][0])
    output["astronomical_dawn_mjd"] = output["night"].map(lambda night: boundaries[night][1])
    output["hours_after_astronomical_dusk"] = (
        output["median_mjd"] - output["astronomical_dusk_mjd"]
    ) * 24.0
    output["hours_before_astronomical_dawn"] = (
        output["astronomical_dawn_mjd"] - output["median_mjd"]
    ) * 24.0
    dusk_distance = output["hours_after_astronomical_dusk"].abs()
    dawn_distance = output["hours_before_astronomical_dawn"].abs()
    output["nearest_astronomical_twilight"] = np.where(
        dusk_distance <= dawn_distance, "dusk", "dawn"
    )
    output["nearest_twilight_abs_hours"] = np.minimum(dusk_distance, dawn_distance)
    output["nearest_twilight_signed_hours"] = np.where(
        output["nearest_astronomical_twilight"].eq("dusk"),
        output["hours_after_astronomical_dusk"],
        -output["hours_before_astronomical_dawn"],
    )
    return output


def add_sky_context(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    valid = (
        np.isfinite(output["median_ra_deg"])
        & np.isfinite(output["median_dec_deg"])
        & np.isfinite(output["median_mjd"])
    )
    output["ecliptic_lon_j2000_deg"] = np.nan
    output["ecliptic_lat_j2000_deg"] = np.nan
    output["solar_elongation_deg"] = np.nan
    if valid.any():
        coordinates = SkyCoord(
            ra=output.loc[valid, "median_ra_deg"].to_numpy() * u.deg,
            dec=output.loc[valid, "median_dec_deg"].to_numpy() * u.deg,
            frame="icrs",
        )
        ecliptic = coordinates.transform_to(
            GeocentricMeanEcliptic(equinox=Time("J2000"))
        )
        times = Time(output.loc[valid, "median_mjd"].to_numpy(), format="mjd", scale="utc")
        output.loc[valid, "ecliptic_lon_j2000_deg"] = ecliptic.lon.deg
        output.loc[valid, "ecliptic_lat_j2000_deg"] = ecliptic.lat.deg
        apparent_coordinates = coordinates.transform_to(GCRS(obstime=times))
        output.loc[valid, "solar_elongation_deg"] = apparent_coordinates.separation(
            get_sun(times)
        ).deg
    return output


def parse_l2_midpoints(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path, dtype={"night": "string"})
    frame["night"] = frame["night"].str.zfill(8)
    strict = frame["strict_standard_catalog"].astype(str).str.lower().isin({"true", "1"})
    frame = frame.loc[strict].copy()
    mjd = pd.to_numeric(frame.get("mjd"), errors="coerce")
    exposure_raw = pd.to_numeric(frame.get("exposure_s"), errors="coerce")
    frame["exposure_fallback_used"] = exposure_raw.isna()
    exposure = exposure_raw.fillna(30.0)
    frame["median_mjd"] = mjd + exposure / (2.0 * 86400.0)
    frame["median_ra_deg"] = pd.to_numeric(frame.get("center_ra_deg"), errors="coerce")
    frame["median_dec_deg"] = pd.to_numeric(frame.get("center_dec_deg"), errors="coerce")
    frame["exposure_s"] = exposure
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--products-dir", type=Path, required=True)
    parser.add_argument("--review-snapshot-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--l2-manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    location = parse_site(config)

    products = args.products_dir.resolve()
    catalog = pd.read_parquet(products / "unknown_catalog.parquet")
    review_detections = pd.read_parquet(products / "unknown_review_detections.parquet")
    catalog["night"] = catalog["night"].astype(str).str.zfill(8)
    review_detections["night"] = review_detections["night"].astype(str).str.zfill(8)
    catalog["trk_sub"] = catalog["trk_sub"].astype(str).str.strip()
    review_detections["trk_sub"] = review_detections["trk_sub"].astype(str).str.strip()

    review_detections["detection_key"] = (
        review_detections["night"].astype(str)
        + "::"
        + review_detections["image_name"].astype(str)
        + "::"
        + review_detections["obj_id"].astype(str)
    )
    review_detections["linkage_detection_key"] = (
        review_detections["night"].astype(str)
        + "::"
        + review_detections["linkage_id"].astype(str)
        + "::"
        + review_detections["detection_key"]
    )

    grouped = (
        review_detections.groupby(["night", "trk_sub", "linkage_id"], as_index=False)
        .agg(
            review_member_rows=("detection_key", "size"),
            unique_detection_n=("detection_key", "nunique"),
            first_mjd=("mjd", "min"),
            median_mjd=("mjd", finite_median),
            last_mjd=("mjd", "max"),
            median_ra_deg=("ra_win_deg", finite_median),
            median_dec_deg=("dec_win_deg", finite_median),
            median_mag_aper4=("mag_aper4", finite_median),
            median_mag_psf=("mag_psf", finite_median),
            median_x_px=("x_win_px", finite_median),
            median_y_px=("y_win_px", finite_median),
        )
    )
    link = catalog.merge(
        grouped,
        on=["night", "trk_sub", "linkage_id"],
        how="left",
        validate="one_to_one",
    )
    link["time_span_hours"] = (link["last_mjd"] - link["first_mjd"]) * 24.0
    link["speed_arcsec_per_hour"] = link["lin_speed_arcsec_per_day"] / 24.0

    review_status = pd.read_csv(
        args.review_snapshot_dir / "review_and_mpc_status.csv",
        dtype={"origin_night": "string", "trk_sub": "string"},
    )
    review_status["origin_night"] = review_status["origin_night"].str.zfill(8)
    review_status["trk_sub"] = review_status["trk_sub"].str.strip()
    status = review_status.rename(columns={"origin_night": "night"})[
        ["night", "trk_sub", "linkage_id", "final_paper_status"]
    ]
    link = link.merge(status, on=["night", "trk_sub", "linkage_id"], how="left", validate="one_to_one")
    link["initial_human_selected"] = link["final_paper_status"].notna()
    link["posthoc_retained"] = link["final_paper_status"].eq("retained_after_posthoc_audit")
    review_detections = review_detections.merge(
        status,
        on=["night", "trk_sub", "linkage_id"],
        how="left",
        validate="many_to_one",
    )

    link = add_sky_context(link)
    link = add_twilight_metrics(link, location)
    link.to_csv(output / "unknown_all_links.csv", index=False)
    review_detections.to_csv(output / "unknown_all_link_memberships.csv", index=False)
    review_detections.drop_duplicates("detection_key").to_csv(
        output / "unknown_unique_detections.csv", index=False
    )

    multiplicity = (
        review_detections.groupby("detection_key", as_index=False)
        .agg(
            linkage_membership_n=("linkage_id", "size"),
            linkage_n=("linkage_id", "nunique"),
            nights=("night", lambda values: ";".join(sorted(set(values)))),
            image_name=("image_name", "first"),
            obj_id=("obj_id", "first"),
        )
    )
    multiplicity.to_csv(output / "unknown_detection_multiplicity.csv", index=False)

    retained_link = link.loc[link["posthoc_retained"]].copy()
    retained_detection = review_detections.loc[
        review_detections["final_paper_status"].eq("retained_after_posthoc_audit")
    ].copy()
    retained_link.to_csv(output / "unknown_high_confidence_links_recomputed.csv", index=False)
    retained_detection.to_csv(output / "unknown_high_confidence_detections_recomputed.csv", index=False)

    exposure_context = pd.DataFrame()
    if args.l2_manifest:
        exposure_context = parse_l2_midpoints(args.l2_manifest)
        exposure_context = add_sky_context(exposure_context)
        exposure_context = add_twilight_metrics(exposure_context, location)
        exposure_context.to_csv(output / "exposure_twilight_context.csv", index=False)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "site_status": "provisional_scheduler_known_coordinates_pending_author_confirmation",
        "iers_status": "bundled_table_auto_download_disabled; out-of-range UT1 uses Astropy degraded-accuracy fallback",
        "site": config["site_configuration_as_executed"]["scheduler"],
        "catalog_linkages": int(len(link)),
        "catalog_nights": int(link["night"].nunique()),
        "fit_ok_linkages": int(link["fit_ok"].fillna(False).sum()),
        "is_good_linkages": int(link["is_good"].fillna(False).sum()),
        "linkage_detection_memberships": int(len(review_detections)),
        "globally_unique_detection_keys": int(review_detections["detection_key"].nunique()),
        "duplicate_memberships": int(len(review_detections) - review_detections["detection_key"].nunique()),
        "detection_keys_shared_by_multiple_links": int((multiplicity["linkage_n"] > 1).sum()),
        "initial_human_selected_linkages": int(link["initial_human_selected"].sum()),
        "posthoc_retained_linkages": int(link["posthoc_retained"].sum()),
        "posthoc_retained_memberships": int(len(retained_detection)),
        "posthoc_retained_unique_detection_keys": int(retained_detection["detection_key"].nunique()),
        "exposure_denominator_rows": int(len(exposure_context)),
        "exposure_midpoint_fallback_30s_rows": int(
            exposure_context.get("exposure_fallback_used", pd.Series(dtype=bool)).fillna(False).sum()
        ),
        "input_hashes": {
            "unknown_catalog.parquet": sha256(products / "unknown_catalog.parquet"),
            "unknown_review_detections.parquet": sha256(products / "unknown_review_detections.parquet"),
        },
    }
    if args.l2_manifest and args.l2_manifest.exists():
        summary["input_hashes"]["l2_manifest.csv"] = sha256(args.l2_manifest)
    (output / "unknown_population_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
