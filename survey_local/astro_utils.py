from __future__ import annotations

from datetime import date as dt_date, datetime, timezone
import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, CIRS, EarthLocation, SkyCoord, get_body, get_sun
from astropy.coordinates.errors import NonRotationTransformationWarning
from astropy.time import Time
from astropy.utils import iers
from astroplan import moon

iers.conf.auto_download = False
iers.conf.auto_max_age = None
iers.conf.iers_degraded_accuracy = "warn"
warnings.filterwarnings("ignore", category=NonRotationTransformationWarning)


def build_site(lat_deg: float, lon_deg: float, height_m: float) -> EarthLocation:
    return EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg, height=height_m * u.m)


def local_night_window_utc(night: dt_date, site: EarthLocation, sun_altitude_deg: float = -12.0) -> tuple[Time, Time]:
    utc_start = datetime.combine(night, datetime.min.time(), tzinfo=timezone.utc)
    grid = Time(utc_start) + np.arange(0, 24 * 60 + 1) * u.min
    sun_alt = get_sun(grid).transform_to(AltAz(obstime=grid, location=site)).alt.deg
    mask = sun_alt < sun_altitude_deg
    idx = np.where(mask)[0]
    if idx.size == 0:
        raise RuntimeError(f"No night window found for {night.isoformat()}")
    return grid[idx[0]], grid[idx[-1]]


def altitude_ok(coords: SkyCoord, obstime: Time, site: EarthLocation, min_altitude_deg: float) -> np.ndarray:
    altaz = coords.transform_to(AltAz(obstime=obstime, location=site))
    return altaz.alt.deg >= min_altitude_deg


def apparent_coords(coords: SkyCoord, obstime: Time, site: EarthLocation) -> SkyCoord:
    frame = CIRS(obstime=obstime, location=site)
    return coords.transform_to(frame)


def moon_separation_limit_deg(illumination: float, new_deg: float, half_deg: float, full_deg: float) -> float:
    a = 2.0 * full_deg + 2.0 * new_deg - 4.0 * half_deg
    b = full_deg - new_deg - a
    return a * illumination * illumination + b * illumination + new_deg


def moon_clearance_ok(
    coords: SkyCoord,
    obstime: Time,
    site: EarthLocation,
    new_deg: float,
    half_deg: float,
    full_deg: float,
    margin_deg: float = 0.0,
) -> np.ndarray:
    moon_pos = get_body("moon", obstime, location=site).transform_to(CIRS(obstime=obstime, location=site))
    field_pos = apparent_coords(coords, obstime, site)
    illum = float(moon.moon_illumination(obstime))
    limit_deg = moon_separation_limit_deg(illum, new_deg, half_deg, full_deg)
    sep_deg = field_pos.separation(moon_pos).deg
    return sep_deg >= (limit_deg + margin_deg)


def stable_visibility_mask(
    coords: SkyCoord,
    start_time: Time,
    site: EarthLocation,
    min_altitude_deg: float,
    visibility_window_minutes: int,
    min_visible_fraction: float,
    moon_sep_args: tuple[float, float, float, float],
) -> np.ndarray:
    effective_window = max(int(visibility_window_minutes), 0)
    offsets = np.arange(0, effective_window + 1, 10)
    if offsets.size == 0:
        offsets = np.array([0], dtype=int)
    visible_at_start = altitude_ok(coords, start_time, site, min_altitude_deg)
    visible_at_start &= moon_clearance_ok(coords, start_time, site, *moon_sep_args)
    visible_counts = np.zeros(len(coords), dtype=np.int16)
    for minute in offsets:
        t = start_time + minute * u.min
        visible_now = altitude_ok(coords, t, site, min_altitude_deg)
        visible_now &= moon_clearance_ok(coords, t, site, *moon_sep_args)
        visible_counts += visible_now.astype(np.int16)
    visible_fraction = visible_counts / float(offsets.size)
    return visible_at_start & (visible_fraction >= min_visible_fraction)
