"""Readiness of a completed observing night's raw images and L2 catalogs."""
from datetime import datetime, timedelta
import json
from pathlib import Path

from .server_common import TZ, runtime


def science_stem(path):
    name = path.name
    for suffix in ('.fits.gz', '.fits'):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    else:
        return None
    lower = name.lower()
    if lower == 'obj_start' or lower.startswith('obj_start_'):
        return None
    if any(token in lower for token in ('focus', 'bias', 'dark', 'flat', 'calib')):
        return None
    return name.removesuffix('_cat')


def night_finished(night, config, now):
    published = runtime(config)/'published'/f'{night}.json'
    if published.exists():
        metadata = json.loads(published.read_text())
        plan = json.loads((Path(metadata['revision'])/'plan.json').read_text())
        return now >= datetime.fromisoformat(plan['end_local']).astimezone(TZ)
    # No v2 plan: wait for morning civil twilight, rather than treating a
    # quiet gap between exposures during the night as completion.
    morning = datetime.strptime(night, '%Y%m%d').date()+timedelta(days=1)
    if now.date() != morning:
        return now.date() > morning
    import astropy.units as u
    from astropy.coordinates import AltAz, get_sun
    from astropy.time import Time
    from .geometry import setup_iers, site
    from .server_common import BASE
    defaults = json.loads((BASE/'defaults.json').read_text())
    setup_iers(BASE/'cache'/config['iers_file'])
    frame = AltAz(obstime=Time(now), location=site(defaults), pressure=0*u.hPa)
    return bool(get_sun(Time(now)).transform_to(frame).alt.deg >= -6)


def readiness(night, config, now=None):
    now = now or datetime.now(TZ)
    result = dict(night=night, ready=False)
    if not night_finished(night, config, now):
        return {**result, 'status': 'waiting_for_observing_end'}
    l2 = Path(config['processed_root'])/night/'L2'
    raw = Path(config.get('raw_root', '/raw1'))/night/'OBJ'
    catalogs = [p for p in l2.iterdir() if p.is_file() and science_stem(p)] if l2.is_dir() else []
    images = [p for p in raw.iterdir() if p.is_file() and science_stem(p)] if raw.is_dir() else []
    if not catalogs:
        return {**result, 'status': 'no_l2'}
    if not images:
        return {**result, 'status': 'waiting_for_raw_inventory'}
    by_stem = {}
    for p in catalogs:
        by_stem.setdefault(science_stem(p), []).append(p)
    missing, stale = [], []
    for p in images:
        matches = by_stem.get(science_stem(p), [])
        if not matches:
            missing.append(p.name)
        elif max(q.stat().st_ctime_ns for q in matches) < p.stat().st_ctime_ns:
            stale.append(p.name)
    latest = max(max(p.stat().st_mtime, p.stat().st_ctime) for p in images+catalogs)
    result.update(raw_files=len(images), l2_files=len(catalogs), missing_l2=missing, stale_l2=stale,
                  last_change_local=datetime.fromtimestamp(latest, TZ).isoformat())
    if missing or stale:
        return {**result, 'status': 'waiting_for_l2_completion'}
    if now.timestamp()-latest < config.get('l2_stable_s', 120):
        return {**result, 'status': 'waiting_for_stable_input'}
    return {**result, 'ready': True, 'status': 'ready'}


def recent_nights(config, now):
    start = config.get('l2_watch_start', '20261118')
    return [night for days in range(1, config.get('l2_watch_lookback_days', 3)+1)
            if (night := (now-timedelta(days=days)).strftime('%Y%m%d')) >= start]
