"""Fixed MP lists for a disconnected telescope: no feedback, waits or publishing.

The background prior controls a precomputed exposure list, not a saturation
guarantee. Original footprints, solar geometry and the complete cyclic route
are checked before a list is exported.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from functools import lru_cache
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, get_sun
from astropy.time import Time
from scipy.optimize import brentq
from .catalog import BASE, write_csv
from .geometry import Sky, coord, isot, setup_iers


def load_config(path=None):
    c = json.loads(Path(path or BASE/'offline_defaults.json').read_text())
    if not 3 <= c['minimum_visits'] <= c['visits'] or c['visits'] != 5:
        raise ValueError('One scheme: five visits, minimum three')
    if c['vmin_arcsec_hour'] != 3 or c['vmax_arcsec_hour'] <= 63:
        raise ValueError('Preserve the original lower speed limit; expand upper only')
    if not 0 < c['cadence_min_s'] <= c['cadence_target_s'] <= c['cadence_max_s']:
        raise ValueError('Invalid cadence bounds')
    if c['min_footprint_alt_deg'] < 25 or c['sun_start_deg'] > -10:
        raise ValueError('Outside approved geometry')
    return c


def load_grid():
    path = BASE/'data/mp_footprints.json'
    source = json.loads(path.read_text())
    fields = []
    for row in source['rows']:
        corners = coord([row[f'corner_ra_{k}'] for k in range(1,5)],
                        [row[f'corner_dec_{k}'] for k in range(1,5)])
        center = coord(row['center_ra'], row['center_dec'])
        fields.append({'field_id':f"MP{int(row['field_id']):04d}",
                       'ra_deg':row['center_ra'], 'dec_deg':row['center_dec'],
                       'radius_deg':float(center.separation(corners).deg.max()),
                       'corners_ra_deg':corners.ra.deg.tolist(),
                       'corners_dec_deg':corners.dec.deg.tolist()})
    if len(fields) != 1855 or len({f['field_id'] for f in fields}) != len(fields):
        raise ValueError('Wrong original MP grid')
    return fields, {'original_fits_sha256':source['source_sha256'],
                    'frozen_json_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def make_sky():
    return Sky(json.loads((BASE/'defaults.json').read_text()))


def night_windows(sky, day, c):
    noon = datetime.fromisoformat(day).replace(hour=12, tzinfo=ZoneInfo('Asia/Shanghai')).timestamp()
    grid = noon + np.arange(0, 86401, 300)
    _, sun, _, _ = sky.environment(grid)
    roots = {}
    for level in (c['sun_start_deg'], c['dark_sun_deg']):
        ix = np.flatnonzero((sun.alt.deg[:-1]-level)*(sun.alt.deg[1:]-level)<0)
        def equation(t):
            tt = Time(t, format='unix')
            return float(get_sun(tt).transform_to(AltAz(obstime=tt,location=sky.location,
                                                      pressure=0*u.hPa)).alt.deg)-level
        values = [brentq(equation, grid[i], grid[i+1], xtol=.001) for i in ix]
        if len(values)!=2:
            raise ValueError('No ordinary pair of solar crossings')
        roots[str(level)] = values
    ev, am = roots[str(c['sun_start_deg'])]
    return {'EV':[ev,ev+c['window_s']], 'AM':[am-c['window_s'],am],
            'dark_crossings':roots[str(c['dark_sun_deg'])]}


def slew(a, b):
    if a['field_id']==b['field_id']:
        return 0.
    dra = abs((b['ra_deg']-a['ra_deg']+180)%360-180)
    ddec = abs(b['dec_deg']-a['dec_deg'])
    return float(9.11+max(11.95*np.sqrt(dra),12.30*np.sqrt(ddec)))


def sky_factor(sun_alt, altitude, elongation, moon_fraction, moon_sep, moon_alt, c):
    """Relative to dark zenith; geometry factors are uncalibrated heuristics."""
    h = np.asarray(sun_alt)
    x = np.clip(-h-5, 0, 10)
    mu = 11.40+1.567*x-.064*x*x
    mu = np.where(h < -15, 20.67+(c['dark_sky_R_mag_arcsec2']-20.67)*
                  np.clip((-h-15)/3,0,1), mu)
    excess = np.maximum(0, 10**(.4*(c['dark_sky_R_mag_arcsec2']-mu))-1)
    excess = np.where(h <= -18, 0., excess)
    airmass = 1/np.sin(np.deg2rad(np.maximum(altitude, 5)))
    directional = 1+c['near_sun_enhancement']*np.exp(-np.asarray(elongation)/35)
    lunar = np.where(np.asarray(moon_alt)>0, 8*np.asarray(moon_fraction)*
                     np.exp(-np.asarray(moon_sep)/40), 0)
    return airmass*(1+excess*directional+lunar)


@lru_cache(maxsize=32)
def anchor_background(dark_sky, enhancement):
    return float(sky_factor(-10,25,40,0,180,-10,
                 {'dark_sky_R_mag_arcsec2':dark_sky,'near_sun_enhancement':enhancement}))


def exposure(sun_alt, altitude, elongation, moon_fraction, moon_sep, moon_alt, c):
    if sun_alt <= c['dark_sun_deg']:
        return c['max_exposure_s']
    anchor = anchor_background(c['dark_sky_R_mag_arcsec2'],c['near_sun_enhancement'])
    value = c['anchor_exposure_s']*anchor/sky_factor(
        sun_alt,altitude,elongation,moon_fraction,moon_sep,moon_alt,c)
    quantum = c['exposure_quantum_s']
    return int(np.clip(np.floor(value/quantum)*quantum,
                       c['min_exposure_s'],c['max_exposure_s']))


class WindowGeometry:
    """Fast planning interpolation, followed by independent direct verification."""
    def __init__(self, sky, fields, bounds, c):
        self.fields=fields;self.bounds=bounds;self.c=c
        self.times=np.arange(bounds[0],bounds[1]+c['late_tolerance_s']+61,60.)
        self.m=sky.metrics([f['ra_deg'] for f in fields],[f['dec_deg'] for f in fields],self.times)

    def at(self, i, t):
        if t < self.times[0]-.001 or t > self.times[-1]+.001:
            return None
        out={k:float(np.interp(t,self.times,self.m[k][i])) for k in ('alt','elong','moon_sep')}
        for k in ('sun_alt','moon_fraction','moon_alt'):
            out[k]=float(np.interp(t,self.times,self.m[k]))
        return out

    def valid(self, i, a, b):
        if a<self.bounds[0]-.001 or b>self.bounds[1]+.001:
            return False
        times=np.r_[a,self.times[(self.times>a)&(self.times<b)],b]
        f=self.fields[i];rad=f['radius_deg'];s=self.m
        alt=np.interp(times,self.times,s['alt'][i])
        moon=np.interp(times,self.times,s['moon_sep'][i])
        moon_alt=np.interp(times,self.times,s['moon_alt'])
        phase=np.interp(times,self.times,s['moon_fraction'])
        elong=np.interp(times,self.times,s['elong'][i])
        sun=np.interp(times,self.times,s['sun_alt'])
        # 0.15deg also guards interpolation and between-sample rotation.
        return bool(np.all(alt-rad>.15+self.c['min_footprint_alt_deg']) and
                    np.all((moon_alt<-.15)|(moon-rad>np.where(phase>=.75,40,30)+.15)) and
                    np.all(elong<=self.c['max_elongation_deg']) and
                    np.all(sun<=self.c['sun_start_deg']+.003))


def route_order(ids, fields):
    """Compare RA order and nearest-neighbour + 2-opt with a closed return."""
    n=len(ids);pts=[fields[i] for i in ids]
    cost=np.array([[slew(a,b) for b in pts] for a in pts])
    ra=np.array([p['ra_deg'] for p in pts]);order=np.argsort(ra).tolist()
    gaps=np.diff(np.r_[ra[order],ra[order[0]]+360]);pivot=(int(np.argmax(gaps))+1)%n
    order=order[pivot:]+order[:pivot]
    def total(path):return sum(cost[path[k-1],path[k]] for k in range(n))
    baseline=total(order);best=order
    for start in range(n):
        path=[start];left=set(range(n))-{start}
        while left:
            nxt=min(left,key=lambda j:(cost[path[-1],j],j));path.append(nxt);left.remove(nxt)
        changed=True
        while changed:
            changed=False
            for i in range(1,n-1):
                for j in range(i+1,n):
                    # Symmetric cost: reversal changes only the two boundary edges.
                    a,b,d,e=path[i-1],path[i],path[j],path[(j+1)%n]
                    delta=cost[a,d]+cost[b,e]-cost[a,b]-cost[d,e]
                    if delta < -1e-8:
                        path=path[:i]+path[i:j+1][::-1]+path[j+1:];changed=True
        if total(path)<total(best):best=path
    return [ids[i] for i in best],baseline,total(best)


def candidate_routes(g, c):
    start,end=g.bounds
    mid=g.at(0,(start+end)/2)  # shared solar state, not used for selection
    del mid
    alt=g.m['alt'];elong=g.m['elong']
    radius=np.array([f['radius_deg'] for f in g.fields])
    moon_ok=(g.m['moon_alt'][0]<-.15)|(g.m['moon_sep'][:,0]-radius>
              (40 if g.m['moon_fraction'][0]>=.75 else 30)+.15)
    available=np.flatnonzero((alt[:,0]-radius>c['min_footprint_alt_deg']+.15)&moon_ok&
                             (elong.min(axis=1)<c['max_elongation_deg']))
    available=sorted(available,key=lambda i:float(elong[i].mean()))[:180]
    if not available:return []
    points=coord([g.fields[i]['ra_deg'] for i in available],[g.fields[i]['dec_deg'] for i in available])
    distance=points[:,None].separation(points[None,:]).deg
    dec=np.array([g.fields[i]['dec_deg'] for i in available])
    routes=[];seen=set()
    nseed=min(c['seed_candidates'],len(available))
    seeds=sorted(set(list(range(nseed//2))+np.linspace(0,len(available)-1,nseed).astype(int).tolist()))
    for j in seeds:
        # Both near-Dec strips and compact near-elongation sets; same objective.
        for penalty in (0.,1.5):
            ranking=np.argsort(distance[j]+penalty*abs(dec-dec[j]))
            for n in c['field_counts']:
                ids=[available[k] for k in ranking[:n]]
                key=tuple(sorted(ids))
                if len(ids)!=n or key in seen:continue
                seen.add(key)
                ordered,baseline,optimized=route_order(ids,g.fields)
                pathpts=coord([g.fields[i]['ra_deg'] for i in ordered],
                              [g.fields[i]['dec_deg'] for i in ordered])
                if np.max(pathpts.separation(pathpts[np.roll(np.arange(n),1)]).deg)>c['max_step_deg']:
                    continue
                score=float(np.mean(elong[ids]))
                routes.append((score,ordered,baseline,optimized))
    # Keep diverse sizes, rather than pruning all longer-cadence groups.
    selected=[]
    for n in c['field_counts']:
        options=sorted((r for r in routes if len(r[1])==n),key=lambda x:(x[0],x[3]))
        if options:
            keep=set(range(min(c['route_candidates'],len(options))))
            keep.update(np.linspace(0,len(options)-1,c['route_candidates']).astype(int))
            selected.extend(options[k] for k in sorted(keep))
    return selected


def simulate_route(g, ids, cycles, c):
    rows=[];t=g.bounds[0];last=None
    for visit in range(1,cycles+1):
        for i in ids:
            f=g.fields[i];turn=c['initial_slew_s'] if last is None else slew(g.fields[last],f)
            change=c['filter_change_s'] if last is None else 0
            begin=t+turn+change
            if begin+c['max_exposure_s']+c['readout_s']+c['late_tolerance_s']>g.bounds[1]:return None
            # Fixed dose uses the brighter geometry over an early/late timing envelope.
            early=max(g.bounds[0],begin-c['initial_early_tolerance_s']-c['early_per_frame_s']*len(rows))
            late=begin+c['late_tolerance_s']+c['max_exposure_s']
            options=[g.at(i,x) for x in (early,begin,late)]
            if any(x is None for x in options):return None
            ex=min(exposure(x['sun_alt'],x['alt'],x['elong'],x['moon_fraction'],
                            x['moon_sep'],x['moon_alt'],c) for x in options)
            finish=begin+ex
            if not g.valid(i,early,finish+c['late_tolerance_s']):return None
            m=g.at(i,begin+ex/2)
            rows.append({'field_id':f['field_id'],'ra_deg':f['ra_deg'],'dec_deg':f['dec_deg'],
                         'field_index':int(i),'visit':visit,'filter':'r','exposure_s':ex,
                         'overhead_start_unix':t,'start_unix':begin,'mid_unix':begin+ex/2,
                         'end_unix':finish,'release_unix':finish+c['readout_s'],
                         'slew_s':turn,'readout_s':c['readout_s'],'filter_change_s':change,
                         'sun_alt_deg':m['sun_alt'],'alt_deg':m['alt'],'elongation_deg':m['elong'],
                         'sky_factor':float(sky_factor(m['sun_alt'],m['alt'],m['elong'],
                                     m['moon_fraction'],m['moon_sep'],m['moon_alt'],c))})
            t=finish+c['readout_s'];last=i
    intervals=[np.diff([r['mid_unix'] for r in rows if r['field_index']==i]) for i in ids]
    gaps=np.concatenate(intervals)
    if gaps.min()<c['cadence_min_s'] or gaps.max()>c['cadence_max_s']:return None
    return rows


def plan_window(sky, fields, bounds, side, c):
    full=WindowGeometry(sky,fields,bounds,c)
    alternatives=[]
    for delay in c['start_offsets_s']:
        g=WindowGeometry.__new__(WindowGeometry)
        g.fields=fields;g.c=c;g.bounds=[bounds[0]+delay,bounds[1]]
        keep=full.times>=g.bounds[0]
        g.times=full.times[keep];g.m={k:v[...,keep] for k,v in full.m.items()}
        alternatives.append((g,candidate_routes(g,c)))
    for cycles in range(c['visits'],c['minimum_visits']-1,-1):
        feasible=[]
        for g,routes in alternatives:
            for _,ids,baseline,optimized in routes:
                # Starts/reversal change which low fields are reached first and last.
                for order in (ids,ids[::-1]):
                    for rotation in (0,len(ids)//2):
                        path=order[rotation:]+order[:rotation]
                        rows=simulate_route(g,path,cycles,c)
                        if rows:
                            gaps=np.concatenate([np.diff([r['mid_unix'] for r in rows if r['field_index']==i]) for i in ids])
                            score=(np.mean([r['elongation_deg'] for r in rows]),
                                   abs(np.median(gaps)-c['cadence_target_s']),optimized)
                            feasible.append((score,rows,path,baseline,optimized,g.bounds))
        if feasible:
            _,rows,path,baseline,optimized,dispatch=min(feasible,key=lambda x:x[0])
            for row in rows:row['side']=side
            return {'side':side,'window_unix':dispatch,'reserved_window_unix':bounds,
                    'rows':rows,'fields':[fields[i] for i in path],
                    'visits':cycles,'ra_order_cycle_slew_s':baseline,'optimized_cycle_slew_s':optimized,
                    'reason':'five visits preferred; minimum three if no longer feasible route',
                    'candidate_routes':sum(len(r) for _,r in alternatives)}
    return {'side':side,'window_unix':bounds,'rows':[],'fields':[],'visits':0,
            'reason':'No MP group meets cadence, full-footprint geometry and timing envelope',
            'candidate_routes':sum(len(r) for _,r in alternatives)}


def verify_window(sky, block, c):
    """Direct ephemerides on <=5s grid; independent of planning interpolation."""
    errors=[];rows=block['rows'];lookup={f['field_id']:f for f in block['fields']}
    samples=[];indices=[];fields=list(lookup.values());idmap={f['field_id']:i for i,f in enumerate(fields)}
    last=None;t=block['window_unix'][0]
    for k,r in enumerate(rows):
        f=lookup[r['field_id']]
        overhead=(c['initial_slew_s']+c['filter_change_s']) if last is None else slew(last,f)
        if abs(r['start_unix']-t-overhead)>.005:errors.append('sequence timing')
        if abs(r['mid_unix']-(r['start_unix']+r['exposure_s']/2))>.001:errors.append('midpoint')
        if r['release_unix']+c['late_tolerance_s']>block['window_unix'][1]+.001:errors.append('deadline')
        a=max(block['window_unix'][0],r['start_unix']-c['initial_early_tolerance_s']-k*c['early_per_frame_s'])
        b=r['end_unix']+c['late_tolerance_s']
        tt=np.linspace(a,b,max(2,int(np.ceil((b-a)/5))+1))
        samples.extend(tt);indices.extend([idmap[r['field_id']]]*len(tt))
        last=f;t=r['release_unix']
    if samples:
        m=sky.metrics([f['ra_deg'] for f in fields],[f['dec_deg'] for f in fields],samples)
        ix=np.asarray(indices);jj=np.arange(len(samples));radius=np.array([f['radius_deg'] for f in fields])[ix]
        if np.any(m['alt'][ix,jj]-radius<=c['min_footprint_alt_deg']+.011):errors.append('footprint altitude')
        if np.any((m['moon_alt']>-.011)&(m['moon_sep'][ix,jj]-radius<np.where(m['moon_fraction']>=.75,40,30)+.011)):errors.append('moon')
        if np.any(m['sun_alt']>c['sun_start_deg']+.003):errors.append('sun')
        if np.any(m['elong'][ix,jj]>c['max_elongation_deg']):errors.append('elongation')
        for f in fields:
            rr=[r for r in rows if r['field_id']==f['field_id']]
            gaps=np.diff([r['mid_unix'] for r in rr])
            if not c['minimum_visits']<=len(rr)<=c['visits']:errors.append('visits')
            if np.any(gaps<c['cadence_min_s']) or np.any(gaps>c['cadence_max_s']):errors.append('cadence')
    return {'passed':not errors,'errors':sorted(set(errors)),'samples':len(samples),
            'scope':'nominal and configured timing envelope, original footprint enclosing cap'}


def make_plan(day, c):
    eop=setup_iers(BASE/'cache'/c['iers_file']);sky=make_sky();fields,provenance=load_grid()
    windows=night_windows(sky,day,c)
    blocks=[plan_window(sky,fields,windows[side],side,c) for side in ('EV','AM')]
    for block in blocks:
        block['verification']=verify_window(sky,block,c)
        if not block['verification']['passed']:raise ValueError(block['verification'])
    return {'date':day,'config':c,'grid':provenance,'iers':eop,'blocks':blocks,
            'status':'OFFLINE_FIXED_LIST_LOCAL_TEST','dark_crossings':windows['dark_crossings'],
            'unallocated_interval_unix':[windows['EV'][1],windows['AM'][0]],
            'assumptions':c['assumptions']}


def write_plan(plan, directory):
    path=Path(directory);path.mkdir(parents=True,exist_ok=True)
    (path/'plan.json').write_text(json.dumps(plan,indent=2,allow_nan=False)+'\n')
    c=plan['config'];token=c['controller_filter_token']
    if token is not None and (not token or any(x.isspace() for x in token)):
        raise ValueError('Invalid controller filter token')
    manifest=[]
    for block in plan['blocks']:
        side=block['side'];rows=block['rows'];lines=[]
        for r in rows:
            p=coord(r['ra_deg'],r['dec_deg'])
            lines.append('\t'.join(['MP_'+r['field_id'][2:],p.ra.to_string(unit=u.hour,sep=':',precision=3,pad=True),
                p.dec.to_string(unit=u.deg,sep=':',precision=2,pad=True,alwayssign=True),
                '0',token or 'UNCONFIRMED_R','1',str(r['exposure_s']),'0']))
        file=path/(side+('.txt' if token else '.draft.txt'))
        file.write_text('\n'.join(lines)+('\n' if lines else ''))
        if rows:write_csv(path/(side+'_timing.csv'),[{**r,'start_utc':isot(r['start_unix']),
                    'mid_utc':isot(r['mid_unix']),'end_utc':isot(r['end_unix'])} for r in rows])
        else:(path/(side+'_timing.csv')).write_text('field_id,start_utc,exposure_s\n')
        start,end=block['window_unix']
        local=lambda t:datetime.fromtimestamp(t,ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')
        manifest.append({'side':side,'list':file.name,'start_local':local(start),
                         'latest_start_local':local(start+c['late_tolerance_s']),
                         'hard_end_local':local(end),'exposures':len(rows),'reason':block['reason'],
                         'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
    (path/'START_TIMES.json').write_text(json.dumps({'date':plan['date'],'lists':manifest,
        'read_before_use':'Sequential fixed lists; no waits and no feedback. Start-time allowance and accumulated slew/readout delays SHARE the configured late envelope. Beyond it stop/replan; TXT cannot enforce this. Morning is next local day. No server deployment.'},indent=2)+'\n')


def main():
    p=argparse.ArgumentParser();p.add_argument('--date',required=True);p.add_argument('--config');p.add_argument('--output',required=True)
    a=p.parse_args();plan=make_plan(a.date,load_config(a.config));write_plan(plan,a.output)
    print(json.dumps({b['side']:{'fields':len(b['fields']),'visits':b['visits'],'frames':len(b['rows']),
                                   'verified':b['verification']['passed']} for b in plan['blocks']}))


if __name__=='__main__':main()
