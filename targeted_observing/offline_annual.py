"""365 fixed night lists and 200 weather replays; schedules never read weather."""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from .catalog import BASE, write_csv
from .annual_sources import dates
from .annual_weather import ensemble
from .geometry import setup_iers, coord
from .offline_twilight import load_config, make_plan, write_plan, make_sky, sky_factor


def day_job(job):
    day,c,root=job
    plan=make_plan(day,c)
    write_plan(plan,Path(root)/day)
    return day,sum(len(b['rows']) for b in plan['blocks'])


def audit_nominal(plans):
    rows=[];exposures=[]
    for p in plans:
        for b in p['blocks']:
            rr=b['rows'];gaps=[];spans=[]
            for f in b['fields']:
                times=[r['mid_unix'] for r in rr if r['field_id']==f['field_id']]
                gaps.extend(np.diff(times));spans.append(times[-1]-times[0])
            spent=sum(r['exposure_s']+r['readout_s']+r['slew_s']+r['filter_change_s'] for r in rr)
            if rr and abs(spent-(rr[-1]['release_unix']-b['window_unix'][0]))>.01:
                raise ValueError('Time ledger not closed')
            rows.append({'date':p['date'],'side':b['side'],'fields':len(b['fields']),
                         'exposures':len(rr),'visits':b['visits'],'exposure_s':sum(r['exposure_s'] for r in rr),
                         'overhead_s':spent-sum(r['exposure_s'] for r in rr),'unused_s':3600-spent,
                         'cadence_median_s':float(np.median(gaps)) if gaps else None,
                         'cadence_min_s':float(min(gaps)) if gaps else None,
                         'cadence_max_s':float(max(gaps)) if gaps else None,
                         'span_median_s':float(np.median(spans)) if spans else None,
                         'elongation_median_deg':float(np.median([r['elongation_deg'] for r in rr])) if rr else None,
                         'ra_cycle_slew_s':b.get('ra_order_cycle_slew_s'),
                         'optimized_cycle_slew_s':b.get('optimized_cycle_slew_s'),
                         'reason':b['reason'],'verified':b['verification']['passed']})
            exposures.extend({'date':p['date'],**r} for r in rr)
    return rows,exposures


def replay(plans,c):
    """Fixed list timing, no reactive abort/reordering or exposure changes.

    Out-of-envelope and unsafe frames are counted, not silently rescued by
    a fictitious realtime controller. Toy moving detections use truth IDs;
    they do not measure false association or orbital recovery.
    """
    nights=[p['date'] for p in plans];R=c['weather_realizations']
    weather=ensemble(nights,R,c['seed']);rng=np.random.default_rng(c['seed']+91)
    counts=np.zeros(R,int);unsafe=np.zeros(R,int);envelope=np.zeros(R,int)
    metrics=[];recovery={}
    setup_iers(BASE/'cache'/c['iers_file']);sky=make_sky()
    for day,p in enumerate(plans):
        for b in p['blocks']:
            rows=b['rows']
            if not rows:continue
            slot=0 if b['side']=='EV' else 2
            w=lambda key:weather[key][:,day,slot]
            n=len(rows);ordinal=np.arange(n)
            shifts=w('delay_s')[:,None]+ordinal[None,:]*(w('extra_readout_s')+w('extra_slew_s'))[:,None]
            actual=np.array([r['mid_unix'] for r in rows])[None,:]+shifts
            half=np.array([r['exposure_s']/2 for r in rows])[None,:]
            fields=b['fields'];lookup={f['field_id']:i for i,f in enumerate(fields)}
            ix=np.array([lookup[r['field_id']] for r in rows])
            grid=np.arange(b['window_unix'][0],float((actual+half).max())+16,15)
            geometry=sky.metrics([f['ra_deg'] for f in fields],[f['dec_deg'] for f in fields],grid)
            m={}
            for key in ('alt','elong','moon_sep'):
                m[key]=np.column_stack([np.interp(actual[:,j],grid,geometry[key][ix[j]]) for j in range(n)])
            for key in ('sun_alt','moon_fraction','moon_alt'):
                m[key]=np.interp(actual,grid,geometry[key])
            radius=np.array([fields[i]['radius_deg'] for i in ix])[None,:]
            motion_guard=(half+7.5)*.0043+.01
            safe=(m['alt']-radius-motion_guard>c['min_footprint_alt_deg'])
            safe&=(m['sun_alt']+motion_guard<=c['sun_start_deg']+.003)
            safe&=(m['moon_alt']<-.15)|(m['moon_sep']-radius-motion_guard>np.where(m['moon_fraction']>=.75,40,30))
            safe&=actual+half+c['readout_s']<=b['window_unix'][1]
            enabled=w('observable')&~w('reserved')
            cloud_start=b.get('reserved_window_unix',b['window_unix'])[0]+3600*w('cloud_start_fraction')
            cloud_end=cloud_start+3600*w('cloud_duration_fraction')
            clear=(actual+half<cloud_start[:,None])|(actual-half>cloud_end[:,None])
            valid=safe&clear&enabled[:,None]
            unsafe+=np.sum(~safe&enabled[:,None],axis=1)
            envelope+=np.sum((shifts>c['late_tolerance_s'])&enabled[:,None],axis=1)
            bg=sky_factor(m['sun_alt'],m['alt'],m['elong'],m['moon_fraction'],m['moon_sep'],m['moon_alt'],c)
            x=1/np.sin(np.deg2rad(np.maximum(m['alt'],5)))
            # r30 dark-zenith depth 20 is an assumption, not white-light 20.8.
            depth=20+1.25*np.log10(half*2/30)-1.25*np.log10(bg)-.15*(x-1)-2.5*np.log10(w('seeing_arcsec')[:,None]/2)
            for f in fields:
                jj=np.flatnonzero(np.array([r['field_id']==f['field_id'] for r in rows]))
                passed=valid[:,jj].sum(axis=1)>=c['minimum_visits']
                counts+=passed
                # Moving seeds uniform in the ORIGINAL RA/Dec footprint at first nominal exposure.
                dr=(np.asarray(f['corners_ra_deg'])-f['ra_deg']+180)%360-180
                dmin,dmax=min(f['corners_dec_deg']),max(f['corners_dec_deg'])
                ra0=f['ra_deg']+rng.uniform(dr.min(),dr.max(),R)
                dec0=np.rad2deg(np.arcsin(rng.uniform(np.sin(np.deg2rad(dmin)),np.sin(np.deg2rad(dmax)),R)))
                bearing=rng.uniform(0,2*np.pi,R)
                dt=(actual[:,jj]-rows[jj[0]]['mid_unix'])/3600
                for speed in (15,63,150,300,750):
                    theta=np.deg2rad(speed*dt/3600)
                    r0=np.deg2rad(ra0)[:,None];d0=np.deg2rad(dec0)[:,None];bb=bearing[:,None]
                    dec=np.arcsin(np.clip(np.sin(d0)*np.cos(theta)+np.cos(d0)*np.sin(theta)*np.cos(bb),-1,1))
                    ra=r0+np.arctan2(np.sin(bb)*np.sin(theta)*np.cos(d0),np.cos(theta)-np.sin(d0)*np.sin(dec))
                    delta=(np.rad2deg(ra)-f['ra_deg']+180)%360-180
                    inside=(delta>=dr.min())&(delta<=dr.max())&(np.rad2deg(dec)>=dmin)&(np.rad2deg(dec)<=dmax)
                    trail=speed*(half[:,jj]*2)/3600
                    loss=1.25*np.log10(1+(trail/w('seeing_arcsec')[:,None])**2)
                    for mag in (16,17,18,19,20):
                        probability=.9/(1+np.exp(np.clip((mag-depth[:,jj]+loss)/.3,-50,50)))
                        detected=(rng.random(probability.shape)<probability)&valid[:,jj]&inside
                        recovered=detected.sum(axis=1)>=c['minimum_visits']
                        key=(b['side'],speed,mag)
                        val=recovery.setdefault(key,[0,0,0])
                        val[0]+=int(recovered.sum());val[1]+=R;val[2]+=int(enabled.sum())
            metrics.append({'date':p['date'],'side':b['side'],
                            'scheduled_frames':n,'realizations':R,
                            'safe_clear_frames':int(valid.sum()),'unsafe_frames_if_list_runs_through':int((~safe&enabled[:,None]).sum()),
                            'timing_envelope_exceeded_frames':int(((shifts>c['late_tolerance_s'])&enabled[:,None]).sum())})
    summary={'realizations':R,'seed':c['seed'],'field_windows_with_at_least_three_valid_frames_P10_P50_P90':np.percentile(counts,[10,50,90]).tolist(),
             'unsafe_frames_if_list_runs_through_P10_P50_P90':np.percentile(unsafe,[10,50,90]).tolist(),
             'frames_outside_timing_envelope_P10_P50_P90':np.percentile(envelope,[10,50,90]).tolist(),
             'scope':'fixed-list replay with hypothetical weather and depth; no reactive controller, no NEO discovery rate',
             'depth_assumption':'r30=20 at zenith/dark, 2arcsec seeing; uncalibrated, background/trailing loss included'}
    injection=[{'side':k[0],'speed_arcsec_hour':k[1],'r_mag':k[2],'recovered_at_least_three':v[0],
                'injected':v[1],'observable_not_reserved':v[2],'fraction_all_weather':v[0]/v[1],
                'fraction_observable_not_reserved':v[0]/v[2] if v[2] else 0,
                'association':'truth-known injection; not real linking'} for k,v in sorted(recovery.items())]
    return summary,metrics,injection


def main():
    p=argparse.ArgumentParser();p.add_argument('--config');p.add_argument('--output',required=True)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--days',type=int)
    p.add_argument('--reuse',action='store_true');p.add_argument('--resume',action='store_true')
    a=p.parse_args();c=load_config(a.config)
    root=Path(a.output);root.mkdir(parents=True,exist_ok=True)
    nights=dates(c['start_date'],c['stop_date_exclusive']);nights=nights[:a.days] if a.days else nights
    jobs=[(d,c,str(root/'nights')) for d in nights]
    if a.resume:
        pending=[]
        for job in jobs:
            file=root/'nights'/job[0]/'plan.json'
            try:
                old=json.loads(file.read_text())
                valid=(old['config']==c and all(b['verification']['passed'] for b in old['blocks']) and
                       old['grid']['frozen_json_sha256']==hashlib.sha256((BASE/'data/mp_footprints.json').read_bytes()).hexdigest() and
                       old['iers']['sha256']==hashlib.sha256((BASE/'cache'/c['iers_file']).read_bytes()).hexdigest())
            except (OSError,ValueError,KeyError):valid=False
            if not valid:pending.append(job)
            else:write_plan(old,file.parent)  # Repair an interrupted export after plan.json was written.
        print(f'Reusing {len(jobs)-len(pending)} verified fixed plans',flush=True)
        jobs=pending
    if not a.reuse:
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            for i,f in enumerate(as_completed([pool.submit(day_job,j) for j in jobs]),1):
                day,n=f.result()
                if i%10==0 or i==len(jobs):print(f'Nominal {i}/{len(jobs)}: {day}, {n} frames',flush=True)
    plans=[json.loads((root/'nights'/d/'plan.json').read_text()) for d in nights]
    if any(p['config']!=c for p in plans):raise ValueError('Stale fixed-plan configuration')
    ledger,exposures=audit_nominal(plans);write_csv(root/'time_ledger.csv',ledger)
    write_csv(root/'exposures.csv',exposures)
    print('Nominal complete; replaying fixed lists under 200 weather realizations',flush=True)
    weather,metrics,injection=replay(plans,c)
    write_csv(root/'weather_windows.csv',metrics);write_csv(root/'conditional_detections.csv',injection)
    summary={'nights':len(nights),'windows':len(ledger),'zero_windows':sum(r['exposures']==0 for r in ledger),
             'exposures':len(exposures),'distinct_MP_fields':len({r['field_id'] for r in exposures}),
             'all_verified':all(r['verified'] for r in ledger),'weather':weather,
             'exposure_counts':{str(t):sum(r['exposure_s']==t for r in exposures) for t in range(5,31,5)},
             'source_hashes':{name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in
                              ('offline_twilight.py','offline_annual.py','offline_defaults.json','data/mp_footprints.json')},
             'assumptions':c['assumptions']}
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
