"""Block scheduler. Pure local products; no cron, publishing, SSH or legacy followup."""
from __future__ import annotations
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, GeocentricMeanEcliptic, get_sun
from astropy.time import Time
from .geometry import Sky, coord, offset, isot, unix, local_time, footprint, in_field
from .catalog import interpolate


def field(ra,dec,prefix='TW',version='20261103'):
    key=f'{float(ra)%360:.6f},{float(dec):.6f},{version}'
    return {'field_id':prefix+'_'+hashlib.sha256(key.encode()).hexdigest()[:10],
            'version':version,'ra_deg':round(float(ra)%360,6),'dec_deg':round(float(dec),6),
            'pa_deg':0.0,'core_targets':[]}


def slew(previous, current, config):
    if previous is None: return config['far_slew_s']
    if previous['field_id']==current['field_id']: return 0.0
    sep=coord(previous['ra_deg'],previous['dec_deg']).separation(coord(current['ra_deg'],current['dec_deg'])).deg
    budget=config['adjacent_slew_s'] if sep<=config['adjacent_max_deg'] else config['far_slew_s']
    if config.get('use_header_slew_floor',True):
        # Existing 30-night normal-branch calibration; add P90 residual margin.
        dra=abs((previous['ra_deg']-current['ra_deg']+180)%360-180)
        ddec=abs(previous['dec_deg']-current['dec_deg'])
        empirical=9.24145+max(11.92277*np.sqrt(dra),12.22922*np.sqrt(ddec))+2.0
        budget=max(budget,float(empirical))
    return budget


def visits(fields,start,end,exposure,cycles,config,block_id,mode):
    rows=[];t=start;previous=None
    for visit in range(1,cycles+1):
        cycle=[]
        for f in fields:
            turn=slew(previous,f,config)
            change=config['filter_change_s'] if previous is None else 0.0
            begin=t+turn+change
            finish=begin+exposure
            if finish+config['readout_s']>end: return rows
            cycle.append({'block_id':block_id,'mode':mode,'field_id':f['field_id'],
                          'field_version':f['version'],'ra_deg':f['ra_deg'],'dec_deg':f['dec_deg'],
                          'filter':config['filter'],'visit':visit,'exposure_s':exposure,
                          'slew_s':turn,'filter_change_s':change,'readout_s':config['readout_s'],
                          'overhead_start_utc':isot(t),'start_utc':isot(begin),'mid_utc':isot((begin+finish)/2),
                          'end_utc':isot(finish),'release_utc':isot(finish+config['readout_s'])})
            t=finish+config['readout_s'];previous=f
        rows.extend(cycle)  # Never retain an incomplete final cycle.
    return rows


def free_intervals(window,reservations):
    a,b=window
    if b<=a:return []
    spans=[(a,b)]
    for r in sorted(reservations,key=lambda r:r['start_utc']):
        x,y=unix(r['start_utc']),unix(r['end_utc'])
        if y<=x:raise ValueError('Invalid reservation interval')
        out=[]
        for p,q in spans:
            if y<=p or x>=q:out.append((p,q))
            else:
                if p<x:out.append((p,x))
                if y<q:out.append((y,q))
        spans=out
    return spans


def block_record(name,fields,rows,start,end,quality,extended=False):
    actual_end=unix(rows[-1]['release_utc'])
    return {'block_id':name,'mode':rows[0]['mode'],'fields':fields,'rows':rows,
            'not_before_utc':isot(start),'latest_start_utc':isot(start),
            'hard_deadline_utc':isot(end),'release_utc':isot(actual_end),
            'not_before_local':local_time(start),'hard_deadline_local':local_time(end),
            'execution':'timed per-exposure dispatch; early completion must wait; late start requires revalidation',
            'extended_search':extended,'quality':quality}


def check_rows(sky,fields,rows,mode,extended=False,tracks=None):
    # Whole per-field span is stricter than checking only its individual exposures.
    quality={}
    for f in fields:
        selected=[r for r in rows if r['field_id']==f['field_id']]
        if not selected:return None
        a,b=unix(selected[0]['start_utc']),unix(selected[-1]['end_utc'])
        q=sky.check_field(f,a,b,mode,extended)
        if not q['ok']:return None
        if tracks is not None:
            for k in f['core_targets']:
                if k not in tracks:raise ValueError(f'Missing track {k}')
                ra,dec=interpolate(tracks[k],np.linspace(a,b,17))
                if not np.all(in_field(ra,dec,f['ra_deg'],f['dec_deg'],sky.config['fov_deg'],sky.config['edge_margin_deg'])): return None
        quality[f['field_id']]=q
    return quality


def twilight_candidates(sky,start,end,side,version):
    c=sky.config;tc=c['twilight'];mid=(start+end)/2
    t=Time(mid,format='unix');ef=GeocentricMeanEcliptic(equinox=Time('J2000'),obstime=t)
    sun=get_sun(t).transform_to(ef)
    sign=1 if side=='EV' else -1
    lon,lat=np.meshgrid(sun.lon.deg+sign*np.arange(43,90,3),np.arange(-22,23,3))
    centers=SkyCoord(lon=lon.ravel()*u.deg,lat=lat.ravel()*u.deg,frame=ef).icrs
    # Compact perimeter traversal: every edge, including wrap, ~2.8 deg.
    xs,ys=[-1.4,1.4,1.4,-1.4],[-1.4,-1.4,1.4,1.4]
    groups=[]
    for center in centers:
        pts=offset(center,xs,ys)
        groups.append([field(r,d,version=version) for r,d in zip(pts.ra.deg,pts.dec.deg)])
    ras=np.array([[f['ra_deg'] for f in g] for g in groups]);decs=np.array([[f['dec_deg'] for f in g] for g in groups])
    m=sky.metrics(ras.ravel(),decs.ravel(),[start+70,(start+end)/2,end-15])
    alt=m['alt'].reshape(-1,4,3);elong=m['elong'].reshape(-1,4,3)
    moon=m['moon_sep'].reshape(-1,4,3)
    moonlim=np.where(m['moon_fraction']>=c['bright_moon_fraction'],c['bright_moon_sep_deg'],c['moon_sep_deg'])
    good=(alt.min(axis=(1,2))>=max(c['center_alt_deg'],c['min_alt_deg']+2.2))
    good&=(elong.min(axis=(1,2))>=tc['elongation_deg'][0])&(elong.max(axis=(1,2))<=tc['extended_elongation_deg'][1])
    good&=np.all((moon>=moonlim+2.2)|(m['moon_alt']<0),axis=(1,2))
    beta=coord(ras.ravel(),decs.ravel()).transform_to(ef).lat.deg.reshape(-1,4)
    good&=np.max(abs(beta),axis=1)<=tc['beta_max_deg']
    result=[]
    for n in np.where(good)[0]:
        el=float(elong[n].mean());beta_abs=float(np.abs(beta[n]).mean())
        extended=bool(elong[n].max()>tc['elongation_deg'][1])
        pref=max(tc['preferred_elongation_deg'][0]-el,0,el-tc['preferred_elongation_deg'][1])
        gal_lat=abs(float(centers[n].galactic.b.deg))
        crowd_penalty=.7*max(0,15-gal_lat)
        score=1000-10*pref-beta_abs+.1*float(alt[n].mean())-500*extended-crowd_penalty
        result.append((score,el,extended,groups[n]))
    return sorted(result,key=lambda r:r[0],reverse=True)


def twilight(sky,window,side,previous=None):
    c=sky.config;tc=c['twilight'];a,b=window
    # Fixed half-windows enforce dawn A not-before even if B finishes early.
    names=['A','B'] if side=='EV' else ['B','A']
    slots=[(a,(a+b)/2),((a+b)/2,b)];blocks=[];used=[]
    for label,(start,end) in zip(names,slots):
        name=f'{side}_{label}'
        if any(max(start,unix(r['start_utc']))<min(end,unix(r['end_utc'])) for r in c['reservations']):continue
        candidates=twilight_candidates(sky,start,end,side,c['version'])
        if previous:
            old=next((x for x in previous if x['block_id']==name),None)
            if old:
                # Try retaining three anchors before replacing an entire group.
                # Order each mixed group by its complete cyclic slew budget.
                replacements=[]
                pool={f['field_id']:f for _,_,_,group in candidates[:40] for f in group}
                for anchors in itertools.combinations(old['fields'],3):
                    if not all(sky.check_field(f,start+70,end-15,'twilight',True)['ok'] for f in anchors):continue
                    for f in pool.values():
                        sep=coord(f['ra_deg'],f['dec_deg']).separation(coord([x['ra_deg'] for x in anchors],[x['dec_deg'] for x in anchors])).deg
                        if sep.min()<2.6 or sep.max()>8:continue
                        group=list(anchors)+[f]
                        routes=list(itertools.permutations(group))
                        route=min(routes,key=lambda rr:sum(slew(rr[i-1],rr[i],c) for i in range(4)))
                        # Extended label is based on actual new time geometry, not prior-night status.
                        em=sky.metrics([x['ra_deg'] for x in route],[x['dec_deg'] for x in route],[start+70,end-15])['elong']
                        ext=bool(em.max()>tc['elongation_deg'][1])
                        replacements.append((1e8,None,ext,list(route)))
                candidates=replacements+candidates
                candidates.insert(0,(1e9,None,old['extended_search'],old['fields']))
        for _,el,extended,fs in candidates:
            if used:
                sep=coord([f['ra_deg'] for f in used],[f['dec_deg'] for f in used])[:,None].separation(coord([f['ra_deg'] for f in fs],[f['dec_deg'] for f in fs])[None,:]).deg
                if sep.min()<2.6:continue
                # Preserve A/B inner/outer ordering, including previous-night reuse.
                nowel=float(sky.metrics([f['ra_deg'] for f in fs],[f['dec_deg'] for f in fs],[(start+end)/2])['elong'].mean())
                oldel=float(sky.metrics([f['ra_deg'] for f in used],[f['dec_deg'] for f in used],[(start+end)/2])['elong'].mean())
                if (side=='EV' and nowel<oldel) or (side=='AM' and nowel>oldel):continue
            rows=visits(fs,start,end,tc['exposure_s'],tc['visits'],c,name,'twilight')
            if len(rows)!=tc['visits']*tc['fields_per_group']:continue
            q=check_rows(sky,fs,rows,'twilight',extended)
            if q:
                block=block_record(name,fs,rows,start,end,q,extended)
                oldids={f['field_id'] for x in (previous or []) if x['block_id']==name for f in x['fields']}
                block['retained_anchor_count']=sum(f['field_id'] in oldids for f in fs)
                blocks.append(block);used.extend(fs);break
    return blocks


def trojan(sky,window,fields,tracks,shift_minutes=0):
    c=sky.config;tc=c['trojan']
    intervals=free_intervals(window,c['reservations'])
    durations=sorted({min(tc['block_s'],b-a) for a,b in intervals if b-a>=3600} | {d for d in (5400.,3600.) if d<=tc['block_s']},reverse=True)
    for duration in durations:
        options=[]
        fs=fields if duration>=5400 else fields[:1]
        for a,b in intervals:
            if b-a<duration:continue
            starts=np.arange(a,b-duration+.01,300)
            if not len(starts):continue
            times=np.column_stack((starts+70,starts+duration-tc['handover_s'])).ravel()
            m=sky.metrics([f['ra_deg'] for f in fs],[f['dec_deg'] for f in fs],times)
            alt=m['alt'].reshape(len(fs),-1,2)
            moon=m['moon_sep'].reshape(len(fs),-1,2)
            moonlim=np.where(m['moon_fraction']>=c['bright_moon_fraction'],c['bright_moon_sep_deg'],c['moon_sep_deg']).reshape(-1,2)
            moonup=m['moon_alt'].reshape(-1,2)>0
            for i,start in enumerate(starts):
                if alt[:,i].min()<33:continue
                if np.any((moon[:,i]<moonlim[i]+2.2)&moonup[i]):continue
                score=float(alt[:,i].mean())
                options.append((score,start,start+duration))
        if not options:continue
        best=max(options,key=lambda x:x[0]);desired=best[1]+shift_minutes*60
        # Deliberate 30--60 minute offset when feasible, rather than a negligible soft jitter.
        options.sort(key=lambda x:(abs(x[1]-desired),-x[0]))
        for _,start,end in options:
            rows=visits(fs,start,end-tc['handover_s'],tc['exposure_s'],tc['max_visits_per_field'],c,'TJ','trojan')
            if not rows:continue
            q=check_rows(sky,fs,rows,'trojan',tracks=tracks)
            if q:return block_record('TJ',fs,rows,start,end,q)
    return None


def make_plan(config,date,campaign=None,previous=None,shift_minutes=0):
    from .config import validate_config
    validate_config(config)
    sky=Sky(config);night=sky.night(date)
    blocks=twilight(sky,night['evening'],'EV',previous)
    reasons=[]
    if campaign:
        if not campaign['start_date']<=date<campaign['stop_date']:
            reasons.append('TJ: date outside frozen campaign; design a new version')
        else:
            tj=trojan(sky,night['dark'],campaign['fields'],campaign['tracks'],shift_minutes)
            if tj:blocks.append(tj)
            else:reasons.append('TJ: no continuous eligible >=60min interval with core targets')
    else:reasons.append('TJ: no frozen campaign supplied')
    blocks+=twilight(sky,night['morning'],'AM',previous)
    present={b['block_id'] for b in blocks}
    for name in ['EV_A','EV_B','AM_B','AM_A']:
        if name not in present:reasons.append(name+': no complete eligible group; time released')
    return {'schema_version':1,'status':'LOCAL_SIMULATION_NOT_FOR_TELESCOPE','date':date,'config':config,
            'night':night,'blocks':blocks,'skips':reasons,'observations_completed':False}
