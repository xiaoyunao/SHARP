"""Airless site geometry and spherical footprints; no production imports."""
from __future__ import annotations
import hashlib
from pathlib import Path
from datetime import datetime, timedelta
from zoneinfo import ZoneInfo
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time
from astropy.utils import iers
from scipy.optimize import brentq


def setup_iers(path):
    path = Path(path)
    table = iers.IERS_A.open(str(path))
    iers.conf.auto_download = False
    iers.conf.iers_degraded_accuracy = 'error'
    iers.earth_orientation_table.set(table)
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'mjd_min': float(table['MJD'][0].value), 'mjd_max': float(table['MJD'][-1].value),
            'purpose': 'offline planning; forecast EOP allowed, not an execution certification'}


def site(config):
    s = config['site']
    return EarthLocation.from_geodetic(s['lon_deg']*u.deg, s['lat_deg']*u.deg, s['height_m']*u.m)


def coord(ra, dec):
    return SkyCoord(ra=np.asarray(ra)*u.deg, dec=np.asarray(dec)*u.deg, frame='icrs')


def offset(center, x, y):
    return SkyCoord(lon=np.asarray(x)*u.deg, lat=np.asarray(y)*u.deg,
                    frame=center.skyoffset_frame()).icrs


def footprint(ra, dec, width=3):
    # All boundary vertices, edge midpoints and center. PA=0, spherical offset frame.
    h = width/2
    x = [-h, 0, h, h, h, 0, -h, -h, 0]
    y = [-h, -h, -h, 0, h, h, h, 0, 0]
    return offset(coord(ra, dec), x, y)


def in_field(ra, dec, center_ra, center_dec, width=3, margin=.15):
    ra, dec, r0, d0 = np.deg2rad(ra), np.deg2rad(dec), np.deg2rad(center_ra), np.deg2rad(center_dec)
    dl = ra-r0
    x = np.cos(dec)*np.cos(dl)*np.cos(d0)+np.sin(dec)*np.sin(d0)
    y = np.cos(dec)*np.sin(dl)
    z = -np.cos(dec)*np.cos(dl)*np.sin(d0)+np.sin(dec)*np.cos(d0)
    lon = np.rad2deg(np.arctan2(y,x)); lat = np.rad2deg(np.arcsin(np.clip(z,-1,1)))
    return (np.abs(lon) <= width/2-margin) & (np.abs(lat) <= width/2-margin)


def isot(unix):
    return Time(unix, format='unix', precision=3).isot+'Z'


def unix(value):
    return float(Time(value).unix)


def local_time(value):
    return datetime.fromtimestamp(value, ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')


class Sky:
    def __init__(self, config):
        self.config = config
        self.location = site(config)

    def environment(self, seconds):
        t = Time(np.atleast_1d(seconds), format='unix')
        table = iers.earth_orientation_table.get()
        if np.any(t.mjd < table['MJD'][0].value) or np.any(t.mjd > table['MJD'][-1].value):
            raise ValueError('Date outside supplied IERS table; obtain a covering table')
        frame = AltAz(obstime=t, location=self.location, pressure=0*u.hPa)
        sun = get_sun(t).transform_to(frame)
        moon_body = get_body('moon', t, self.location)
        moon = moon_body.transform_to(frame)
        # Phase triangle; not merely elongation/2 approximation.
        sun_body = get_body('sun', t, self.location)
        elong = sun_body.separation(moon_body)
        phase = np.arctan2(sun_body.distance*np.sin(elong),
                           moon_body.distance-sun_body.distance*np.cos(elong))
        fraction = (1+np.cos(phase.to_value(u.rad)))/2
        return frame, sun, moon, fraction

    def metrics(self, ra, dec, seconds):
        frame, sun, moon, fraction = self.environment(seconds)
        cs = coord(np.atleast_1d(ra), np.atleast_1d(dec))
        aa = cs[:, None].transform_to(frame[None, :])
        return {'alt': aa.alt.deg, 'sun_alt': sun.alt.deg,
                'elong': aa.separation(sun[None, :]).deg,
                'moon_sep': aa.separation(moon[None, :]).deg,
                'moon_alt': moon.alt.deg, 'moon_fraction': fraction}

    def night(self, date):
        noon = datetime.fromisoformat(date).replace(hour=12, tzinfo=ZoneInfo('Asia/Shanghai'))
        a = noon.timestamp()
        seconds = a+np.arange(0, 24*3600+1, 300)
        _, sun, _, _ = self.environment(seconds)
        alt = sun.alt.deg
        roots = {}
        for level in (-12, -18):
            found = []
            for i in np.where((alt[:-1]-level)*(alt[1:]-level)<0)[0]:
                def f(s):
                    t=Time(s,format='unix')
                    return get_sun(t).transform_to(AltAz(obstime=t,location=self.location,pressure=0*u.hPa)).alt.deg-level
                found.append(brentq(f, seconds[i], seconds[i+1], xtol=.01))
            if len(found)!=2:
                raise ValueError(f'No normal dusk/dawn crossings for {date}, level {level}')
            roots[str(level)] = found
        w = self.config['twilight']['window_s']
        ev, am = roots['-12']
        return {'date':date, 'evening':[ev,ev+w], 'morning':[am-w,am],
                'dark':[max(roots['-18'][0], ev+w),min(roots['-18'][1],am-w)]}

    def check_field(self, field, start, end, mode, extended=False, step=30):
        c=self.config
        seconds=np.linspace(start,end,max(2,int(np.ceil((end-start)/step))+1))
        # Add a conservative bound for between-grid altitude change and unsampled edge arcs.
        # Earth's rotation <=0.0042 deg/s; half-step bound plus 0.01 deg numerical allowance.
        guard=.0043*step/2+.01
        fp=footprint(field['ra_deg'],field['dec_deg'],c['fov_deg'])
        m=self.metrics(fp.ra.deg,fp.dec.deg,seconds)
        moon_limit=np.where(m['moon_fraction']>=c['bright_moon_fraction'],c['bright_moon_sep_deg'],c['moon_sep_deg'])
        # Circumscribed angular cap bounds EVERY footprint point, including edge interiors.
        radius=float(fp[-1].separation(fp).deg.max())
        min_alt=float(m['alt'][-1].min()-radius-guard)
        reasons=[]
        if min_alt<=c['min_alt_deg'] or m['alt'][-1].min()<c['center_alt_deg']+guard: reasons.append('altitude')
        moon_up=m['moon_alt']> -guard
        if np.any(m['moon_sep'][-1,moon_up]-radius < moon_limit[moon_up]+guard): reasons.append('moon')
        sun_limit=c['trojan']['sun_alt_deg'] if mode=='trojan' else c['twilight']['sun_alt_deg']
        if np.any(m['sun_alt']>sun_limit+0.002): reasons.append('sun')
        if mode=='twilight':
            lo=c['twilight']['elongation_deg'][0]
            hi=c['twilight']['extended_elongation_deg'][1] if extended else c['twilight']['elongation_deg'][1]
            if m['elong'][-1].min()<lo or m['elong'][-1].max()>hi: reasons.append('elongation')
        return {'ok':not reasons,'reasons':reasons,'min_footprint_alt_bound_deg':min_alt,
                'min_center_alt_deg':float(m['alt'][-1].min()),'min_moon_sep_deg':float(m['moon_sep'].min()),
                'elongation_min_deg':float(m['elong'][-1].min()),'elongation_max_deg':float(m['elong'][-1].max())}
