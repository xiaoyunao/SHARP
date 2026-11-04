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
from .annual_timing import prepare,mask,row_times


@lru_cache(maxsize=16384)
def pixels(ra,dec,nside,margin=0):
    ids=hp.query_disc(nside,hp.ang2vec(ra,dec,lonlat=True),np.deg2rad(2.123),inclusive=True)
    r,d=hp.pix2ang(nside,ids,lonlat=True)
    return ids[in_field(r,d,ra,dec,3,margin)]


def window_pixels(blocks,nside):
    all_pixels=[];counts=[];first=[];last=[];full_pixels=[]
    for b in blocks:
        for f in b['fields']:
            ts=[unix(r['mid_utc']) for r in b['rows'] if r['field_id']==f['field_id']]
            full_pixels.extend(pixels(f['ra_deg'],f['dec_deg'],nside))
            ids=pixels(f['ra_deg'],f['dec_deg'],nside,.15)
            all_pixels.extend(ids);counts.extend([len(ts)]*len(ids));first.extend([min(ts)]*len(ids));last.extend([max(ts)]*len(ids))
    if not all_pixels:return np.array([],dtype=int),np.array([],dtype=int)
    ids,inv=np.unique(all_pixels,return_inverse=True);n=np.zeros(len(ids),dtype=int)
    lo=np.full(len(ids),np.inf);hi=np.full(len(ids),-np.inf)
    np.add.at(n,inv,counts);np.minimum.at(lo,inv,first);np.maximum.at(hi,inv,last)
    return np.unique(full_pixels),ids[(n>=4)&(hi-lo>=720)]


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
        rows=block['rows'];ts=row_times(block)
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
        safe&=mask(block,[r['field_id'] for r in rows],actual)
        trail=initial['speed'][:,None]*3600/86400*30
        loss=1.25*np.log10(1+(trail/(2*w('seeing_arcsec')[:,None]))**2)
        # Approximate twilight penalty from window position; explicitly scenario-based.
        second_half=(block['block_id']=='EV_B' or block['block_id']=='AM_A') and len(rows)==20
        window_start=start-(1800 if second_half else 0)
        window_fraction=np.clip((actual-window_start)/3600,0,1)
        penalty=.8*(1-window_fraction if side=='EV' else window_fraction)
        moon_loss=np.column_stack([np.interp(actual[:,j],block['timing_grid'],block['timing_moon_loss'][r['field_id']]) for j,r in enumerate(rows)])
        altitude=np.column_stack([np.interp(actual[:,j],block['timing_grid'],block['timing_alt'][r['field_id']]) for j,r in enumerate(rows)])
        extinction=.15*(1/np.sin(np.deg2rad(altitude))-1)
        seeing_loss=2.5*np.log10(w('seeing_arcsec')[:,None]/2)
        depth=limiting_mag-loss-penalty-moon_loss-extinction-seeing_loss
        probability=1/(1+np.exp((initial['mag'][:,None]-depth)/.25))
        detected=inside&safe&(rng.random(inside.shape)<probability)
        count+=detected.sum(axis=1)
        lo=np.minimum(lo,np.min(np.where(detected,actual,np.inf),axis=1))
        hi=np.maximum(hi,np.max(np.where(detected,actual,-np.inf),axis=1))
    return (count>=4)&(hi-lo>=720)


def injections(plans,spec,limiting_mag=20.,per_cell=12,causal=None):
    nights=[p['date'] for p in plans];weather=ensemble(nights,spec['weather_realizations'],spec['seed']);rng=np.random.default_rng(spec['seed']+71)
    output=[];speeds=np.array([.1,.3,.7,1.5,5.]);magnitudes=np.array([17.,18.,19.,20.,21.])
    speed,mag=np.meshgrid(speeds,magnitudes);speed=np.repeat(speed.ravel(),per_cell);mag=np.repeat(mag.ravel(),per_cell)
    def selected(index,side,realization):
        if causal is None:return [b for b in plans[index]['blocks'] if b['block_id'].startswith(side)]
        menus,choices,sky=causal;menu=menus[index];slot=0 if side=='EV' else 2
        indices=choices[realization,index,slot]
        names=[side+'_A'] if menu['config']['twilight']['visits']==10 else [side+'_'+x for x in (['A','B'] if side=='EV' else ['B','A'])]
        blocks=[menu['options'][name][int(i)] for name,i in zip(names,indices) if i>=0]
        for block in blocks:prepare(block,sky)
        return blocks
    for index,plan in enumerate(plans):
        for side in ('EV','AM'):
            cohort=(index*2+(side=='AM'))%spec['weather_realizations']
            blocks=selected(index,side,cohort)
            fields=[f for b in blocks for f in b['fields']]
            if not fields:continue
            n=len(speed);choices=rng.integers(0,len(fields),n);xy=rng.uniform(-1.5,1.5,(n,2))
            p=offset(coord([fields[j]['ra_deg'] for j in choices],[fields[j]['dec_deg'] for j in choices]),xy[:,0],xy[:,1])
            ras=p.ra.deg;decs=p.dec.deg
            epoch=min(unix(r['mid_utc']) for b in blocks for r in b['rows'])
            initial={'ra':np.array(ras),'dec':np.array(decs),'bearing':rng.uniform(0,2*np.pi,n),'speed':speed,'mag':mag,
                     'epoch':epoch,'realization':np.full(n,cohort) if causal is not None else np.arange(n)%spec['weather_realizations']}
            first=recover(blocks,initial,weather,index,side,limiting_mag,rng)
            later={}
            for gap in (1,3):
                if index+gap>=len(plans):later[gap]=None;continue
                p=plans[index+gap];result=np.zeros(n,bool)
                for s in ('EV','AM'):
                    bs=selected(index+gap,s,cohort)
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
    p=argparse.ArgumentParser();p.add_argument('--visits',type=int,choices=(5,10),required=True);p.add_argument('--nside',type=int,default=1024)
    p.add_argument('--days',type=int);p.add_argument('--skip-injections',action='store_true');p.add_argument('--causal',action='store_true')
    p.add_argument('--skip-area',action='store_true');p.add_argument('--limiting-mag',type=float,default=20.);args=p.parse_args()
    spec=json.loads((BASE/'annual_defaults.json').read_text());nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    if args.days:nights=nights[:args.days]
    root=BASE/'products/annual'/f'twilight_{args.visits}';plans=[json.loads((root/f'{d}.json').read_text()) for d in nights]
    source_root=root
    if args.days:
        root=root/f'smoke_{args.days}days';root.mkdir(exist_ok=True)
    if not args.skip_area:
        rows,summary=union_area(plans,args.nside);summary['night_count']=len(plans)
        write_csv(root/f'coverage_nside{args.nside}.csv',rows)
        (root/f'coverage_nside{args.nside}.json').write_text(json.dumps(summary,indent=2)+'\n')
    if not args.skip_injections:
        from .geometry import Sky,setup_iers
        setup_iers(BASE/'cache'/spec['iers_file'])
        for plan in plans:
            sky=Sky(plan['config'])
            for block in plan['blocks']:prepare(block,sky)
        causal=None
        if args.causal:
            menus=[json.loads((source_root/'menus'/f'{d}.json').read_text()) for d in nights]
            choices={tuple(row[:3]):tuple(row[3:]) for row in np.load(source_root/'causal_choices.npz')['choices']}
            causal=(menus,choices,Sky(plans[0]['config']))
        rows=injections(plans,spec,limiting_mag=args.limiting_mag,causal=causal)
        tag=f'{"causal_" if args.causal else ""}r{args.limiting_mag:g}'
        write_csv(root/f'moving_injections_{tag}.csv',rows)
        (root/f'injection_assumptions_{tag}.json').write_text(json.dumps({'population':'uniform seeds within scheduled footprints, not normalized NEO population',
            'motion':'great-circle constant angular velocity over 3 days; no acceleration/orbit ambiguity',
            'depth':f'r30 dark {args.limiting_mag:g} mag scenario; logistic width .25mag; .8mag twilight-position penalty; explicit airmass, Moon, seeing and trailing losses',
            'linking':'truth-associated recovery upper bound; no false detections or orbit-link ambiguity model',
            'weather':'shared 200-realization ensemble; one rotating realization per window cohort' if args.causal else 'shared 200-realization ensemble; seeds balanced across realizations',
            'field_policy':'causal weather-feedback anchor queue' if args.causal else 'nominal anchors, not rescheduled after weather',
            'conditioning':'initial window has a scheduled group; closed initial windows are counted separately in causal execution ledger',
            'limitations':'Truth-associated recovery is an upper bound on real linking; end-of-year follow-ups are right-censored.'},indent=2)+'\n')


if __name__=='__main__':main()
