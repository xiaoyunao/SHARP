"""HEALPix footprint unions and conditional moving-source injection recovery."""
import argparse
import json
from collections import Counter,defaultdict
from functools import lru_cache
import numpy as np
import healpy as hp
from .catalog import BASE,write_csv
from .geometry import in_field,coord,offset,unix
from .annual_sources import dates
from .annual_weather import ensemble


@lru_cache(maxsize=16384)
def pixels(ra,dec,nside):
    ids=hp.query_disc(nside,hp.ang2vec(ra,dec,lonlat=True),np.deg2rad(2.123),inclusive=True)
    r,d=hp.pix2ang(nside,ids,lonlat=True)
    return ids[in_field(r,d,ra,dec,3,0)]


def window_pixels(blocks,nside):
    all_pixels=[];counts=[];first=[];last=[]
    for b in blocks:
        for f in b['fields']:
            ts=[unix(r['mid_utc']) for r in b['rows'] if r['field_id']==f['field_id']]
            ids=pixels(f['ra_deg'],f['dec_deg'],nside)
            all_pixels.extend(ids);counts.extend([len(ts)]*len(ids));first.extend([min(ts)]*len(ids));last.extend([max(ts)]*len(ids))
    if not all_pixels:return np.array([],dtype=int),np.array([],dtype=int)
    ids,inv=np.unique(all_pixels,return_inverse=True);n=np.zeros(len(ids),dtype=int)
    lo=np.full(len(ids),np.inf);hi=np.full(len(ids),-np.inf)
    np.add.at(n,inv,counts);np.minimum.at(lo,inv,first);np.maximum.at(hi,inv,last)
    return ids,ids[(n>=4)&(hi-lo>=720)]


def union_area(plans,nside):
    area=hp.nside2pixarea(nside,degrees=True);rows=[];by_night=[];qualified=[]
    for p in plans:
        daily=[];qday=[]
        for side in ('EV','AM'):
            blocks=[b for b in p['blocks'] if b['block_id'].startswith(side)]
            pix,q=window_pixels(blocks,nside);daily.extend(pix);qday.extend(q)
            exposures=sum(len(b['rows']) for b in blocks)
            rows.append({'date':p['date'],'side':side,'nside':nside,'pixel_area_deg2':area,
                'exposure_area_deg2':exposures*9,'unique_area_deg2':len(pix)*area,'qualified_short_arc_area_deg2':len(q)*area})
        by_night.append(np.unique(daily).astype(int));qualified.append(np.unique(qday).astype(int))
    flat=np.concatenate(by_night);unique,counts=np.unique(flat,return_counts=True)
    summary={'nside':nside,'unique_year_deg2':len(unique)*area,'area_on_ge2_nights_deg2':int(np.sum(counts>=2))*area,
             'area_on_ge3_nights_deg2':int(np.sum(counts>=3))*area,'exposure_area_deg2':sum(r['exposure_area_deg2'] for r in rows),
             'qualified_short_arc_unique_year_deg2':len(np.unique(np.concatenate(qualified)))*area,
             'method':'HEALPix pixel-center inclusion in spherical PA=0 3x3 footprints; finite-resolution estimate'}
    for gap in (1,3):
        summary[f'qualified_area_day_plus_{gap}_sum_deg2']=sum(len(np.intersect1d(qualified[i],qualified[i+gap])) for i in range(len(plans)-gap))*area
    return rows,summary


def move(ra0,dec0,bearing,speed,dt_days):
    r=np.deg2rad(ra0)[:,None];d=np.deg2rad(dec0)[:,None];b=bearing[:,None]
    theta=np.deg2rad(speed[:,None]*dt_days)
    dd=np.arcsin(np.clip(np.sin(d)*np.cos(theta)+np.cos(d)*np.sin(theta)*np.cos(b),-1,1))
    rr=r+np.arctan2(np.sin(b)*np.sin(theta)*np.cos(d),np.cos(theta)-np.sin(d)*np.sin(dd))
    return np.rad2deg(rr)%360,np.rad2deg(dd)


def recover(blocks,initial,weather,night_index,side,limiting_mag,rng):
    """One declared hypothetical detection curve, with common realized weather.

    Linear great-circle sky motion is an injection model, not orbital prediction.
    The seeds sample field area, not a debiased or normalized NEO population.
    """
    n=len(initial['ra']);count=np.zeros(n,int);lo=np.full(n,np.inf);hi=np.full(n,-np.inf)
    realizations=initial['realization'];slot=0 if side=='EV' else 2
    for block in blocks:
        rows=block['rows'];ts=np.array([unix(r['mid_utc']) for r in rows])
        w=lambda key:weather[key][realizations,night_index,slot]
        shifts=w('delay_s')[:,None]+np.arange(len(rows))[None,:]*(w('extra_readout_s')+w('extra_slew_s'))[:,None]
        actual=ts[None,:]+shifts
        ra,dec=move(initial['ra'],initial['dec'],initial['bearing'],initial['speed'],(actual-initial['epoch'])/86400)
        inside=in_field(ra,dec,np.array([r['ra_deg'] for r in rows])[None,:],np.array([r['dec_deg'] for r in rows])[None,:],3,.15)
        start=unix(block['not_before_utc']);end=unix(block['hard_deadline_utc'])
        fraction=(actual-start)/(end-start)
        clear=w('observable')&~w('reserved')
        cloud=(fraction>=w('cloud_start_fraction')[:,None])&(fraction<w('cloud_start_fraction')[:,None]+w('cloud_duration_fraction')[:,None])
        safe=(actual+30<end)&~cloud&clear[:,None]
        # A conservative motion-in-time margin avoids declaring shifted edge frames safe.
        margins=np.array([block['quality'][r['field_id']]['min_footprint_alt_bound_deg']-30 for r in rows])
        safe&=shifts*.0043<margins[None,:]
        trail=initial['speed'][:,None]*3600/86400*30
        loss=1.25*np.log10(1+(trail/(2*w('seeing_arcsec')[:,None]))**2)
        # Approximate twilight penalty from window position; explicitly scenario-based.
        penalty=.8*(1-fraction if side=='EV' else fraction)
        depth=limiting_mag-loss-penalty
        probability=1/(1+np.exp((initial['mag'][:,None]-depth)/.25))
        detected=inside&safe&(rng.random(inside.shape)<probability)
        count+=detected.sum(axis=1)
        lo=np.minimum(lo,np.min(np.where(detected,actual,np.inf),axis=1))
        hi=np.maximum(hi,np.max(np.where(detected,actual,-np.inf),axis=1))
    return (count>=4)&(hi-lo>=720)


def injections(plans,spec,limiting_mag=20.,per_cell=12):
    nights=[p['date'] for p in plans];weather=ensemble(nights,spec['weather_realizations'],spec['seed']);rng=np.random.default_rng(spec['seed']+71)
    output=[];speeds=np.array([.1,.3,.7,1.5,5.]);magnitudes=np.array([17.,18.,19.,20.,21.])
    speed,mag=np.meshgrid(speeds,magnitudes);speed=np.repeat(speed.ravel(),per_cell);mag=np.repeat(mag.ravel(),per_cell)
    for index,plan in enumerate(plans):
        for side in ('EV','AM'):
            blocks=[b for b in plan['blocks'] if b['block_id'].startswith(side)]
            fields=[f for b in blocks for f in b['fields']]
            if not fields:continue
            n=len(speed);choices=rng.integers(0,len(fields),n);ras=[];decs=[]
            for j in choices:
                f=fields[j];p=offset(coord(f['ra_deg'],f['dec_deg']),rng.uniform(-1.5,1.5),rng.uniform(-1.5,1.5));ras.append(p.ra.deg);decs.append(p.dec.deg)
            epoch=min(unix(r['mid_utc']) for b in blocks for r in b['rows'])
            initial={'ra':np.array(ras),'dec':np.array(decs),'bearing':rng.uniform(0,2*np.pi,n),'speed':speed,'mag':mag,
                     'epoch':epoch,'realization':np.arange(n)%spec['weather_realizations']}
            first=recover(blocks,initial,weather,index,side,limiting_mag,rng)
            later={}
            for gap in (1,3):
                if index+gap>=len(plans):later[gap]=None;continue
                p=plans[index+gap];result=np.zeros(n,bool)
                for s in ('EV','AM'):
                    bs=[b for b in p['blocks'] if b['block_id'].startswith(s)]
                    result|=recover(bs,initial,weather,index+gap,s,limiting_mag,rng)
                later[gap]=result
            for i in range(n):
                output.append({'date':plan['date'],'side':side,'speed_deg_day':float(speed[i]),'magnitude':float(mag[i]),
                    'weather_realization':int(initial['realization'][i]),'tracklet':bool(first[i]),
                    'recovered_plus1':None if later[1] is None else bool(first[i]&later[1][i]),
                    'recovered_plus3':None if later[3] is None else bool(first[i]&later[3][i]),
                    'within_existing_rate_range':bool(3<=speed[i]*150<=63),
                    'limiting_r30_dark_scenario':limiting_mag})
        if (index+1)%50==0:print(f'Moving-source injection {index+1}/{len(plans)} nights',flush=True)
    return output


def main():
    p=argparse.ArgumentParser();p.add_argument('--visits',type=int,choices=(5,10),required=True);p.add_argument('--nside',type=int,default=512)
    p.add_argument('--days',type=int);p.add_argument('--skip-injections',action='store_true');args=p.parse_args()
    spec=json.loads((BASE/'annual_defaults.json').read_text());nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    if args.days:nights=nights[:args.days]
    root=BASE/'products/annual'/f'twilight_{args.visits}';plans=[json.loads((root/f'{d}.json').read_text()) for d in nights]
    rows,summary=union_area(plans,args.nside);summary['night_count']=len(plans)
    write_csv(root/f'coverage_nside{args.nside}.csv',rows)
    (root/f'coverage_nside{args.nside}.json').write_text(json.dumps(summary,indent=2)+'\n')
    if not args.skip_injections:
        rows=injections(plans,spec);write_csv(root/'moving_injections.csv',rows)
        (root/'injection_assumptions.json').write_text(json.dumps({'population':'uniform seeds within scheduled footprints, not normalized NEO population',
            'motion':'great-circle constant angular velocity over 3 days; no acceleration/orbit ambiguity',
            'depth':'r30 dark 20 mag scenario; logistic width .25mag; .8mag twilight-position penalty',
            'linking':'truth-associated recovery upper bound; no false detections or orbit-link ambiguity model',
            'weather':'shared 200-realization ensemble; seeds balanced across realizations; nominal field choices',
            'limitations':'Nominal anchors are not rescheduled here after bad weather. Report this separately from a causal execution queue.'},indent=2)+'\n')


if __name__=='__main__':main()
