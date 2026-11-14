"""Audit a static controller list in sky hour angle, without claiming encoder safety."""
import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import HADec
from astropy.time import Time

from .catalog import BASE, write_csv
from .full_night_test import read_lines, replay, READOUT_S, FILTER_CHANGE_S
from .geometry import coord, footprint, local_time, setup_iers
from .offline_twilight import make_sky, load_grid
from .controller_fields import is_fixed_name, canonical_name
from .fixed_fields import FixedFields, REGISTRY_PATH


def signed_ra_step(ra_from, ra_to):
    return (ra_to-ra_from+180) % 360-180


def audit(source, start, output, ha_reference_hours=5):
    source, output = Path(source), Path(output)
    raw = source.read_bytes()
    commands = read_lines(raw.decode(), 'ALL')
    if not commands:
        raise ValueError('Empty controller list')
    if any(c['line'].split()[4] not in ('g','r','i','z') for c in commands):
        raise ValueError('Controller filter tokens must be lowercase')
    setup_iers(BASE/'cache/iers_20261104.all')
    sky = make_sky()
    start_time = datetime.fromisoformat(start)
    if start_time.tzinfo is None:
        raise ValueError('Start must have an explicit timezone')
    start_unix = start_time.timestamp()
    rows = replay(commands, start_unix)
    grid, _ = load_grid()
    mp = {f['field_id']: f for f in grid}

    def ha(ra, dec, seconds):
        t = Time(np.asarray(seconds), format='unix')
        return coord(ra, dec).transform_to(HADec(
            obstime=t, location=sky.location, pressure=0*u.hPa)).ha.wrap_at(180*u.deg).hour

    # A <=5-second sample includes exposure, readout and endpoint checks.
    seconds, indices = [], []
    for i, r in enumerate(rows):
        tt = np.linspace(r['start_unix'], r['release_unix'],
                         int(np.ceil((r['release_unix']-r['start_unix'])/5))+1)
        seconds.extend(tt)
        indices.extend([i]*len(tt))
    seconds, indices = np.asarray(seconds), np.asarray(indices)
    ra = np.asarray([r['ra_deg'] for r in rows])[indices]
    dec = np.asarray([r['dec_deg'] for r in rows])[indices]
    hours = ha(ra, dec, seconds)
    summaries = []
    for name in dict.fromkeys(r['field_id'] for r in rows):
        ids = [i for i,r in enumerate(rows) if r['field_id']==name]
        mask = np.isin(indices, ids)
        first = rows[ids[0]]
        m = sky.metrics(first['ra_deg'], first['dec_deg'], seconds[mask])
        if is_fixed_name(name):
            registered=FixedFields.read(REGISTRY_PATH).by_id[canonical_name(name)]
            if coord(first['ra_deg'],first['dec_deg']).separation(
                    coord(registered['ra_deg'],registered['dec_deg'])).arcsec>.03:
                raise ValueError('Fixed field centre changed: '+name)
            fp=footprint(first['ra_deg'],first['dec_deg'])
            radius=float(fp[-1].separation(fp).deg.max())
            alt_limit,sun_limit=30,-18
        elif name.startswith('MP_'):
            radius = mp[name.replace('MP_','MP')]['radius_deg']
            alt_limit, sun_limit = 25, -10
        else:
            fp = footprint(first['ra_deg'],first['dec_deg'])
            radius = float(fp[-1].separation(fp).deg.max())
            alt_limit, sun_limit = 30, -18
        # This audit assumes each target label identifies one fixed sky position.
        if any(abs(rows[i]['ra_deg']-first['ra_deg'])>1e-7 or
               abs(rows[i]['dec_deg']-first['dec_deg'])>1e-7 for i in ids):
            raise ValueError('Same label has different coordinates: '+name)
        moon_bad = (m['moon_alt']>-.021) & (m['moon_sep'][0]-radius <
            np.where(m['moon_fraction']>=.75,40,30)+.021)
        record = dict(field=name, frames=len(ids),
            first_exposure_local=local_time(first['start_unix']),
            last_readout_local=local_time(rows[ids[-1]]['release_unix']),
            ha_min_h=float(hours[mask].min()), ha_max_h=float(hours[mask].max()),
            center_alt_min_deg=float(m['alt'].min()), center_alt_max_deg=float(m['alt'].max()),
            footprint_alt_lower_bound_deg=float(m['alt'].min()-radius-.021),
            sun_alt_max_deg=float(m['sun_alt'].max()),
            sky_geometry_ok=bool(m['alt'].min()-radius-.021>alt_limit and
                m['sun_alt'].max()<=sun_limit and not np.any(moon_bad)))
        summaries.append(record)
    transitions = []
    for a,b in zip(rows,rows[1:]):
        if not b['slew_s']:
            continue
        raw_step = b['ra_deg']-a['ra_deg']
        step = signed_ra_step(a['ra_deg'],b['ra_deg'])
        ts = [b['overhead_start_unix'], b['start_unix']]
        # Envelope of old/new sightlines while slewing: conditional on a monotonic
        # short-way move in the same HA branch. NOT a measured mechanical path.
        bounds = np.r_[ha(a['ra_deg'],a['dec_deg'],ts),ha(b['ra_deg'],b['dec_deg'],ts)]
        transitions.append(dict(command=b['command_index'], departure_local=local_time(ts[0]),
            from_field=a['field_id'], to_field=b['field_id'], raw_ra_step_deg=raw_step,
            short_ra_step_deg=step, crosses_ra_zero=bool(abs(raw_step-step)>180),
            modeled_slew_s=b['slew_s'], ha_envelope_min_h=float(bounds.min()),
            ha_envelope_max_h=float(bounds.max())))
    min_ha=min([float(hours.min())]+[x['ha_envelope_min_h'] for x in transitions])
    max_ha=max([float(hours.max())]+[x['ha_envelope_max_h'] for x in transitions])
    result=dict(input=str(source.resolve()),input_sha256=hashlib.sha256(raw).hexdigest(),
        start_local=local_time(start_unix),end_local=local_time(rows[-1]['release_unix']),
        controller_lines=len(commands),exposures=len(rows),reference_abs_ha_h=ha_reference_hours,
        tracking_ha_min_h=float(hours.min()),tracking_ha_max_h=float(hours.max()),
        conditional_slew_and_tracking_ha_min_h=min_ha,conditional_slew_and_tracking_ha_max_h=max_ha,
        within_reference_ha=bool(min_ha>-ha_reference_hours and max_ha<ha_reference_hours),
        sky_geometry_ok=all(x['sky_geometry_ok'] for x in summaries),
        mechanical_safety_certified=False,
        limitations=['HA reference is customary, not a supplied encoder hard limit.',
            'Initial mount position, pier/pointing state, cable turns and actual slew branch are unknown.',
            'Slew HA envelope assumes monotonic shortest-way motion; no physical trajectory or collision model.',
            f'Nominal {READOUT_S:g}s readout, {FILTER_CHANGE_S:g}s filter change and empirical slew; no live timing/stop feedback.'],
        fields=summaries,large_or_wrap_transitions=[x for x in transitions
            if abs(x['short_ra_step_deg'])>30 or x['crosses_ra_zero']])
    output.mkdir(parents=True,exist_ok=True)
    (output/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    write_csv(output/'transitions.csv',transitions)
    write_csv(output/'exposures.csv',[dict(r,start_local=local_time(r['start_unix']),
        readout_end_local=local_time(r['release_unix'])) for r in rows])
    write_csv(output/'fields.csv',summaries)
    assert source.read_bytes()==raw, 'Source was modified during audit'
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input',required=True)
    p.add_argument('--start',required=True,help='ISO date-time including UTC offset')
    p.add_argument('--output',required=True)
    a=p.parse_args()
    print(json.dumps(audit(a.input,a.start,a.output),indent=2))
