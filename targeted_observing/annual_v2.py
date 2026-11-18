"""Isolated annual replay of the deployed v2 with causal simulated histories.

Run only in a disposable working directory under a Slurm allocation. This entry
never publishes telescope lists, changes cron, submits MPC records, or reads a
future weather draw when selecting fields. Products are planning experiments.
"""
import argparse
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import pickle
import sys
import time

import numpy as np
from astropy.time import Time

from targeted_observing import server_scheduler as scheduler
from targeted_observing import offline_twilight as twilight
from targeted_observing import shared_fields as shared_module
from targeted_observing.annual_sources import dates
from targeted_observing.annual_weather import ensemble, random_error
from targeted_observing import annual_weather
from targeted_observing.annual_periods import injected, search
from targeted_observing.catalog import BASE, write_csv
from targeted_observing.fixed_fields import tan_xy
from targeted_observing.geometry import setup_iers, coord
from targeted_observing.server_common import load_config, write_json
from targeted_observing.server_history import sampling
from targeted_observing.shared_fields import SharedFields, cloud_labels

sys.path.insert(0,str(Path(__file__).resolve().parent))
from annual_v2_fast import simulate_route, GridSky

CTX=None


class Context:
    def __init__(self,root):
        self.root=Path(root);self.config=load_config();self.config['base']=str(self.root/'isolated')
        self.config['processed_root']=str(self.root/'processed')
        self.weather_prior_path=self.root.parent/'data/weather_prior.json'
        if not self.weather_prior_path.exists():self.weather_prior_path=BASE/'data/weather_prior.json'
        annual_weather.BASE=self.weather_prior_path.parent.parent
        self.catalog=SharedFields();self.grid=self.catalog.twilight_grid();self.clouds=cloud_labels()
        self.tracks=scheduler.load_tracks();self.ids=list(self.tracks)
        self.arrays={id(t):dict(jd=np.asarray(t['datetime_jd']),ra=np.unwrap(np.deg2rad(t['RA'])),dec=np.asarray(t['DEC']),
                              V=np.asarray(t['V']),delta=np.asarray(t['delta'])) for t in self.tracks.values()}
        self.evidence=json.loads((BASE/'data/tjn_queue_evidence_20261118.json').read_text())['objects']
        self.spec=json.loads((BASE/'annual_defaults.json').read_text())
        self.iers=setup_iers(BASE/'cache'/self.config['iers_file']);self.sky=twilight.make_sky()
        physical_sources=[Path(scheduler.__file__),Path(twilight.__file__),Path(__file__),Path(__file__).with_name('annual_v2_fast.py'),BASE/'offline_defaults.json',BASE/'server_defaults.json',BASE/'annual_defaults.json',self.weather_prior_path]
        self.model_hash=hashlib.sha256(b''.join(p.read_bytes() for p in physical_sources)).hexdigest()
        self.policy_hash=hashlib.sha256(Path(twilight.__file__).read_bytes()+(BASE/'offline_defaults.json').read_bytes()+Path(__file__).with_name('annual_v2_fast.py').read_bytes()).hexdigest()
        self.stars=scheduler.load_stars();self.original_geometry=scheduler.Geometry
        self.original_make_plan=scheduler.make_plan;self.original_validation=scheduler.final_validation
        self.active=None;self.history={}
        # Read immutable data once per process, not once per simulated night.
        scheduler.load_tracks=lambda:self.tracks
        scheduler.SharedFields=lambda:self.catalog
        scheduler.cloud_labels=lambda:self.clouds
        scheduler.load_stars=lambda:self.stars
        scheduler.setup_iers=lambda path:self.iers
        twilight.setup_iers=lambda path:self.iers
        scheduler.states=lambda config,periods=None:self.history
        scheduler.make_sky=lambda:self.sky
        twilight.simulate_route=simulate_route
        scheduler.interpolate=self.interpolate
        shared_module.interpolate=self.interpolate

    def interpolate(self,track,seconds):
        a=self.arrays[id(track)];q=np.atleast_1d(seconds)/86400+2440587.5
        if q.min()<a['jd'][0] or q.max()>a['jd'][-1]:raise ValueError('Ephemeris outside frozen coverage')
        return np.rad2deg(np.interp(q,a['jd'],a['ra']))%360,np.interp(q,a['jd'],a['dec'])

    def load_day(self,day):
        with (self.root/'geometry'/f'{day}.pkl').open('rb') as stream:plan,shared,yso=pickle.load(stream)
        values=[]
        for attributes in (shared,yso):
            item=self.original_geometry.__new__(self.original_geometry);item.__dict__.update(attributes);values.append(item)
        self.active=(plan,*values)
        scheduler.make_plan=lambda *a,**k:self.active[0]
        scheduler.Geometry=lambda sky,fields,*a:self.active[2] if len(fields)==2 else self.active[1]
        cached=GridSky(self.sky,self.active[1],self.active[2])
        p=self.active[0];begin=p['blocks'][0]['reserved_window_unix'][0];end=p['blocks'][1]['reserved_window_unix'][1]
        duration=2*self.config['yso_frames_per_band']*(self.config['yso_exposure_s']+self.config['readout_s'])+self.config['filter_change_s']
        lo=max(begin+3600,p['dark_crossings'][0]);hi=min(end-3600,p['dark_crossings'][1])
        mids=np.arange(lo+224,hi-duration+1,300)+duration/2
        cached.scalar_exact={}
        if len(mids):
            fields=list(scheduler.YSO_FIELDS.values());exact=self.sky.metrics([f['ra_deg'] for f in fields],[f['dec_deg'] for f in fields],mids)
            for i,f in enumerate(fields):
                for j,t in enumerate(mids):
                    cached.scalar_exact[(f['ra_deg'],f['dec_deg'],float(t))]={k:np.asarray([[v[i,j]]]) if v.ndim==2 else np.asarray([v[j]]) for k,v in exact.items()}
        scheduler.make_sky=lambda:cached
        scheduler.final_validation=lambda rows,sky,catalog,config,tracks:self.original_validation(rows,cached,catalog,config,tracks)
        return cached


def prepare_day(day):
    c=CTX;start=time.time();path=c.root/'geometry'/f'{day}.pkl'
    marker=c.root/'geometry_versions'/f'{day}.json'
    if path.exists() and marker.exists() and json.loads(marker.read_text()).get('policy_hash')==c.policy_hash:
        return dict(date=day,status='cached')
    plan=c.original_make_plan(day,twilight.load_config(),grid=c.grid,allow_empty=True,shorter_fallback=True)
    if path.exists():
        with path.open('rb') as stream:_,shared,yso=pickle.load(stream)
    else:
        begin=plan['blocks'][0]['reserved_window_unix'][0];end=plan['blocks'][1]['reserved_window_unix'][1]
        radius=c.grid[0][0]['radius_deg']
        shared=c.original_geometry(c.sky,list(c.catalog.by_id.values()),begin,end+1800,radius).__dict__
        yso=c.original_geometry(c.sky,[dict(field_id=k,**v) for k,v in scheduler.YSO_FIELDS.items()],begin,end+1800,radius).__dict__
    path.parent.mkdir(parents=True,exist_ok=True);temporary=path.with_suffix('.tmp')
    with temporary.open('wb') as stream:pickle.dump((plan,shared,yso),stream,protocol=5)
    temporary.replace(path)
    result=dict(date=day,seconds=time.time()-start,twilight_frames=sum(len(b['rows']) for b in plan['blocks']),
                shorter_segments=sum(b.get('shorter_fallback',False) for b in plan['blocks']),policy_hash=c.policy_hash)
    write_json(marker,result)
    return result


def curve_assumption(key,evidence):
    entry=evidence.get(key,{});periods=entry.get('periods_h',[])
    rng=np.random.default_rng(int(hashlib.sha256(('v2-period:'+key).encode()).hexdigest()[:8],16))
    period=float(rng.choice(periods)) if periods else float(np.exp(rng.uniform(np.log(2),np.log(1000))))
    known=bool(periods and entry.get('role') in ('control','known_other'))
    return dict(period_h=period,external_prior=known,role=entry.get('role','literature_unchecked'),
                source='frozen literature alternatives' if periods else 'illustrative log-uniform 2-1000h',
                shape=str(rng.choice(['asymmetric','symmetric','complex'],p=[.6,.2,.2])),phase=float(rng.uniform(0,2*np.pi)))


def frame_geometry(plan,grid,offset):
    rows=plan['rows'];q=np.array([r['mid_unix'] for r in rows])+offset
    ra=np.array([r['ra_deg'] for r in rows]);dec=np.array([r['dec_deg'] for r in rows])
    metrics={key:np.zeros(len(rows)) for key in ('alt','moon_sep','sun_alt','moon_alt','moon_fraction')}
    valid=np.ones(len(rows),bool);radius=CTX.grid[0][0]['radius_deg']
    for j,row in enumerate(rows):
        g,i=grid.lookup[(row['ra_deg'],row['dec_deg'])]
        duration=row['exposure_s']/2+row['readout_s']
        check=np.array([q[j]-row['exposure_s']/2,q[j],q[j]+duration])
        if check.min()<g.times[0] or check.max()>g.times[-1]:valid[j]=False;continue
        a=np.interp(check,g.times,g.metrics['alt'][i]);moon=np.interp(check,g.times,g.metrics['moon_sep'][i])
        up=np.interp(check,g.times,g.metrics['moon_alt']);ph=np.interp(check,g.times,g.metrics['moon_fraction'])
        sun=np.interp(check,g.times,g.metrics['sun_alt']);ha=np.interp(check,g.times,g.ha[i])
        altitude=25 if row['mode'] in ('EV','AM','FOLLOW') else 30
        sunlim=-10 if row['mode'] in ('EV','AM','FOLLOW') else -18
        valid[j]=bool(np.all(a-radius>altitude+.01) and np.all((up<-.01)|(moon-radius>np.where(ph>=.75,40,30)+.01))
                      and np.all(sun<=sunlim+.003) and np.all(abs(ha)<5))
        for key in ('alt','moon_sep'):metrics[key][j]=np.interp(q[j],g.times,g.metrics[key][i])
        for key in ('sun_alt','moon_alt','moon_fraction'):metrics[key][j]=np.interp(q[j],g.times,g.metrics[key])
    return q,metrics,valid,ra,dec


def observe(plan,grid,scenario,weather,realization,ni,rng):
    """Realize photons only after the day's complete queue has been chosen."""
    c=CTX;rows=plan['rows'];n=len(rows);nominal=scenario=='clear_baseline'
    photo=c.spec['photometry_scenarios']['baseline' if nominal else scenario]
    slot=np.array([0 if r['mode']=='EV' else 2 if r['mode']=='AM' else 1 for r in rows])
    w=lambda key:weather[key][realization,ni]
    if nominal:
        offset=np.zeros(n);usable=np.ones(n,bool);seeing=np.full(n,2.);zp=0.
    else:
        # A static controller does not reset the clock at program transitions.
        turns=np.cumsum([r['slew_s']>0 for r in rows]);offset=w('delay_s')[0]+np.arange(n)*w('extra_readout_s')[0]+turns*w('extra_slew_s')[0]
        usable=w('observable')[slot]&~w('reserved')[slot]
        nominal_time=np.array([r['mid_unix'] for r in rows])
        for s in (0,1,2):
            hit=np.flatnonzero(slot==s)
            if not len(hit):continue
            fraction=(nominal_time[hit]-nominal_time[hit].min())/max(1,np.ptp(nominal_time[hit]))
            cloud=(fraction>=w('cloud_start_fraction')[s])&(fraction<w('cloud_start_fraction')[s]+w('cloud_duration_fraction')[s])
            usable[hit]&=~cloud
        seeing=w('seeing_arcsec')[slot];zp=float(weather['zeropoint_standard_normal'][realization,ni]*photo['night_zeropoint'])
    q,m,geometry,ra0,dec0=frame_geometry(plan,grid,offset);usable&=geometry
    # All scientific programs contribute TJN detections, including MP and YSO.
    jd=q/86400+2440587.5;ras=[];decs=[];mags=[];dist=[]
    for key in c.ids:
        a=c.arrays[id(c.tracks[key])]
        ras.append(np.rad2deg(np.interp(jd,a['jd'],a['ra']))%360);decs.append(np.interp(jd,a['jd'],a['dec']))
        mags.append(np.interp(jd,a['jd'],a['V']));dist.append(np.interp(jd,a['jd'],a['delta']))
    ras=np.asarray(ras);decs=np.asarray(decs);mags=np.asarray(mags);dist=np.asarray(dist)
    x,y=tan_xy(ras,decs,ra0[None,:],dec0[None,:]);finite=np.isfinite(x)&np.isfinite(y)
    x=np.where(finite,x,0.);y=np.where(finite,y,0.);angle=np.deg2rad(c.catalog.manifest['clockwise_roll_deg'])
    xx=np.cos(angle)*x-np.sin(angle)*y;yy=np.sin(angle)*x+np.cos(angle)*y
    hx,hy=np.rad2deg(np.tan(np.deg2rad([c.catalog.manifest['tan_width_east_deg']/2,c.catalog.manifest['tan_height_north_deg']/2])))
    margin=c.catalog.manifest['quality_margin_deg']+c.catalog.manifest['quality_guard_deg']
    inside=finite&(abs(xx)<=hx-margin)&(abs(yy)<=hy-margin)&usable[None,:]
    # Outside candidate images skip the expensive error/probability calculation.
    obj,frame=np.nonzero(inside)
    if not len(obj):return usable,{},0,0
    ex=np.array([r['exposure_s'] for r in rows]);airmass=1/np.sin(np.deg2rad(np.maximum(5,m['alt'][frame])))
    phase=np.where(m['moon_alt'][frame]>0,m['moon_fraction'][frame],0.)
    errors=random_error(mags[obj,frame],airmass,phase,m['moon_sep'][frame],seeing[frame],ex[frame],photo)
    if nominal:detected=errors<=.217147
    else:
        probability=1/(1+np.exp(np.clip((errors-.217147)/.035,-100,100)))
        detected=(rng.random(len(obj))<probability)&(rng.random(len(obj))<photo['usable_fraction'])
    samples=defaultdict(list);r_count=0
    for o,j,error,d in zip(obj[detected],frame[detected],errors[detected],dist[obj[detected],frame[detected]]):
        band=rows[j]['filter'];good=bool(band=='r' and error<=.05 and (nominal or abs(zp)<=.02))
        samples[c.ids[o]].append(dict(mid=float(q[j]),night=plan['night'],band=band,error=float(error),good=good,
                                     zeropoint=zp,emission=float(q[j]-d*499.004783836),mode=rows[j]['mode']))
        r_count+=band=='r'
    return usable,samples,int(detected.sum()),int(r_count)


def interpolation_audit(grid,day):
    c=CTX;g=c.active[1];rng=np.random.default_rng(int(day.replace('-','')))
    fields=[g.fields[i] for i in rng.choice(len(g.fields),16,replace=False)]
    ts=np.linspace(g.times[2],g.times[-3],33)+7.123
    ra=[f['ra_deg'] for f in fields];dec=[f['dec_deg'] for f in fields]
    exact=c.sky.metrics(ra,dec,ts);approx=grid.metrics(ra,dec,ts)
    errors={k:float(np.max(abs(exact[k]-approx[k]))) for k in ('alt','moon_sep','sun_alt','moon_alt','moon_fraction')}
    all_errors=dict(errors)
    errors['alt']=float(np.max(abs(exact['alt']-approx['alt'])[exact['alt']<80]))
    errors['moon_sep']=float(np.max(abs(exact['moon_sep']-approx['moon_sep'])[exact['moon_sep']>20]))
    if max(errors[k] for k in ('alt','moon_sep','sun_alt','moon_alt'))>.003:
        raise ValueError(dict(date=day,interpolation_error=errors))
    return dict(date=day,max_error_away_from_zenith_or_moon=errors,max_error_all_samples=all_errors)


def simulate_case(job):
    scenario,realization,nights=job;c=CTX;name=scenario if scenario=='clear_baseline' else f'{scenario}_r{realization:02d}'
    root=c.root/'cases'/name;root.mkdir(parents=True,exist_ok=True)
    result_path=root/'summary.json'
    if result_path.exists() and json.loads(result_path.read_text()).get('model_hash')==c.model_hash:return json.loads(result_path.read_text())
    weather=ensemble(nights,max(32,realization+1),20261118)
    rng=np.random.default_rng(int(hashlib.sha256(name.encode()).hexdigest()[:8],16))
    c.history={};history=defaultdict(list);field_daily=np.zeros((len(nights),len(c.catalog.ids)),np.int32)
    planned_daily=np.zeros_like(field_daily);yso_names=list(scheduler.YSO_FIELDS)
    yso_daily=np.zeros((len(nights),len(yso_names)),np.int32);yso_planned=np.zeros_like(yso_daily);nightly=[];detected=Counter();r_detected=Counter();mode_count=Counter();audits=[];samples_saved={}
    checkpoint=root/'checkpoint.pkl';next_index=0
    if checkpoint.exists():
        with checkpoint.open('rb') as stream:saved=pickle.load(stream)
        if saved.get('model_hash')==c.model_hash:
            next_index=saved['next_index'];history=saved['history'];c.history=saved['states']
            for key in ('field_daily','planned_daily','yso_daily','yso_planned','nightly','detected','r_detected','mode_count','audits','samples_saved'):
                if key=='field_daily':field_daily=saved[key]
                elif key=='planned_daily':planned_daily=saved[key]
                elif key=='yso_daily':yso_daily=saved[key]
                elif key=='yso_planned':yso_planned=saved[key]
                elif key=='nightly':nightly=saved[key]
                elif key=='detected':detected=saved[key]
                elif key=='r_detected':r_detected=saved[key]
                elif key=='mode_count':mode_count=saved[key]
                elif key=='audits':audits=saved[key]
                elif key=='samples_saved':samples_saved=saved[key]
            rng.bit_generator.state=saved['rng_state']
            print(json.dumps(dict(case=name,resumed_from_night=next_index)),flush=True)
    for ni in range(next_index,len(nights)):
        day=nights[ni];grid=c.load_day(day)
        try:plan=scheduler.build(day,c.config,[])
        except Exception as exc:
            write_json(root/'error.json',dict(date=day,error=str(exc),last_checkpoint=ni-ni%7))
            raise ValueError(f'{name} {day}: {exc}') from exc
        if scenario=='clear_baseline' and ni%30==0:
            independent=c.original_validation(plan['rows'],c.sky,c.catalog,c.config,c.tracks)
            if not independent['passed']:raise ValueError(independent)
            audits.append(interpolation_audit(grid,day))
            write_json(root/f'audit_plan_{day}.json',plan)
        usable,samples,count,r_count=observe(plan,grid,scenario,weather,realization,ni,rng)
        for j,row in enumerate(plan['rows']):
            suffix=row['field_id'].rsplit('_',1)[-1]
            if suffix in c.catalog.by_id:
                k=int(suffix)-1;planned_daily[ni,k]+=1;field_daily[ni,k]+=int(usable[j])
            if row['field_id'] in yso_names:
                k=yso_names.index(row['field_id']);yso_planned[ni,k]+=1;yso_daily[ni,k]+=int(usable[j])
            mode_count[row['mode']]+=int(usable[j])
        for key,new in samples.items():
            detected[key]+=len(new);r_detected[key]+=sum(row['band']=='r' for row in new)
            history[key].extend(new)
            if not c.history.get(key,{}).get('sampling_ready'):
                state=sampling(history[key],c.config,c.evidence.get(key,{}).get('periods_h',[]))
                c.history[key]=state
                if state['sampling_ready']:
                    # Freeze a reproducible completion certificate, not a forecast count.
                    good=[r for r in history[key] if r['good'] and r['band']=='r']
                    samples_saved[key]=good.copy()
        nightly.append(dict(date=day,planned_frames=len(plan['rows']),usable_frames=int(usable.sum()),
                            tjn_measurements=count,tjn_r_measurements=r_count,objects_detected=len(samples),
                            cumulative_detected=len(detected),sampling_ready=len(samples_saved),
                            yso=plan['yso'],mp_sampling_failures=sum(not v['sampling_ok'] for v in plan['verification']['mp_sampling'])))
        if (ni+1)%7==0 or ni+1==len(nights):
            saved=dict(model_hash=c.model_hash,next_index=ni+1,history=history,states=c.history,field_daily=field_daily,
                planned_daily=planned_daily,yso_daily=yso_daily,yso_planned=yso_planned,nightly=nightly,detected=detected,
                r_detected=r_detected,mode_count=mode_count,audits=audits,samples_saved=samples_saved,rng_state=rng.bit_generator.state)
            temporary=checkpoint.with_suffix('.tmp')
            with temporary.open('wb') as stream:pickle.dump(saved,stream,protocol=5)
            temporary.replace(checkpoint)
        if (ni+1)%30==0 or ni+1==len(nights):print(json.dumps(dict(case=name,night=ni+1,detected=len(detected),sampling_ready=len(samples_saved))),flush=True)
    np.savez_compressed(root/'field_counts.npz',actual=field_daily,planned=planned_daily,dates=nights,field_ids=c.catalog.ids,yso_ids=yso_names,yso_actual=yso_daily,yso_planned=yso_planned)
    write_csv(root/'nightly.csv',nightly);write_json(root/'sampling_certificates.json',c.history)
    ids=[];offsets=[0];times=[];errors=[];zeropoints=[];assumptions=[]
    for key,rows in samples_saved.items():
        ids.append(key);times.extend(r['emission'] for r in rows);errors.extend(r['error'] for r in rows);zeropoints.extend(r['zeropoint'] for r in rows)
        offsets.append(len(times));assumptions.append(curve_assumption(key,c.evidence))
    np.savez_compressed(root/'curves.npz',object_id=ids,offsets=offsets,emission=times,error=errors,zeropoint=zeropoints)
    write_json(root/'period_assumptions.json',dict(zip(ids,assumptions)))
    object_rows=[dict(object_id=k,cloud=c.clouds[k],role=c.evidence.get(k,{}).get('role'),measurements=detected[k],r_measurements=r_detected[k],
                      **{x:c.history.get(k,{}).get(x) for x in ('good_points','good_nights','baseline_days','sampling_ready')}) for k in sorted(detected)]
    write_csv(root/'objects.csv',object_rows)
    result=dict(case=name,model_hash=c.model_hash,scenario=scenario,realization=realization,nights=len(nights),detected_objects=len(detected),
                detected_r_objects=sum(v>0 for v in r_detected.values()),tjn_measurements=sum(detected.values()),tjn_r_measurements=sum(r_detected.values()),
                sampling_ready_objects=len(samples_saved),planned_frames=sum(r['planned_frames'] for r in nightly),
                usable_frames=sum(r['usable_frames'] for r in nightly),shared_fields_visited=int(np.count_nonzero(field_daily.sum(axis=0))),
                mode_usable=dict(mode_count),interpolation_audits=audits)
    write_json(result_path,result)
    checkpoint.unlink(missing_ok=True)
    (root/'error.json').unlink(missing_ok=True)
    return result


def period_case(name):
    c=CTX;root=c.root/'cases'/name;curves=dict(np.load(root/'curves.npz'));assumptions=json.loads((root/'period_assumptions.json').read_text())
    records=[];counts={str(a):dict(blind=0,known_prior=0,blind_false_unique=0,blind_weak_or_unchecked=0) for a in (.1,.3,.6)}
    for i,key in enumerate(curves['object_id']):
        lo,hi=curves['offsets'][i:i+2];ts=curves['emission'][lo:hi];t=(ts-ts.min())/3600
        err=curves['error'][lo:hi];zp=curves['zeropoint'][lo:hi];a=assumptions[str(key)];period=a['period_h']
        for amplitude in (.1,.3,.6):
            rng=np.random.default_rng(int(hashlib.sha256(f'{name}:{key}:{amplitude}'.encode()).hexdigest()[:8],16))
            y=injected(t,period,amplitude,a['phase'],a['shape'])+zp+rng.normal(0,err)
            for mode in ('blind','known_prior'):
                if mode=='known_prior' and not a['external_prior']:continue
                fit=search(t,y,err,None if mode=='blind' else period) if 2<=period<=1000 else None
                recovered=bool(fit and abs(fit['period_h']/period-1)<.01)
                unique=bool(fit and fit['unique_in_searched_candidates']);success=recovered and unique
                counts[str(amplitude)][mode]+=success
                if mode=='blind':
                    counts[str(amplitude)]['blind_false_unique']+=unique and not recovered
                    counts[str(amplitude)]['blind_weak_or_unchecked']+=success and a['role'] not in ('control','known_other')
                records.append(dict(object_id=str(key),amplitude=amplitude,mode=mode,**a,points=len(t),
                                    recovered=recovered,unique=unique,success=success,best_period_h=fit['period_h'] if fit else None))
    write_csv(root/'period_recovery.csv',records);write_json(root/'period_summary.json',counts)
    print(json.dumps(dict(period_case=name,qualified_objects=len(curves['object_id']),counts=counts)),flush=True)
    return dict(case=name,counts=counts)


def run(root,workers=16,day_limit=None,prepare_only=False):
    global CTX
    CTX=Context(root);nights=dates('2026-11-01','2027-11-01')
    if day_limit:nights=nights[:day_limit]
    provenance=dict(model_hash=CTX.model_hash,twilight_policy='full hour preferred; shorter >=3 visits fallback; skip unavailable; 4.95h planning HA guard',start=nights[0],stop=nights[-1],nights=len(nights),candidates=len(CTX.ids),config=load_config(),
                    weather=json.loads(CTX.weather_prior_path.read_text()),photometry=CTX.spec['photometry_scenarios'],
                    zero_seed_history=True,followups=[],engine='deployed shared-grid server_scheduler.build',
                    verification='native twilight direct verification; native final validator using exact-Astropy 60s grids; 13 dates direct full validation and interpolation audits',
                    limitations=['1063 preselected frozen candidates, not complete TJN population','r=V approximation; non-r detections also assume same magnitude',
                        'weather and photometry are scenarios, not fitted site performance','simulated sampling completion is not an observed period',
                        'unknown periods log-uniform 2-1000h, fixed per object; shape mixture .6 asymmetric/.2 symmetric/.2 complex',
                        'amplitude sensitivity .1/.3/.6mag; uniqueness heuristic delta chi2>=9, not calibrated posterior',
                        'no future weather used by queue; every case updates its own realized detections causally'])
    write_json(Path(root)/'provenance.json',provenance)
    ctx=multiprocessing.get_context('fork')
    with ProcessPoolExecutor(max_workers=workers,mp_context=ctx) as pool:
        futures={pool.submit(prepare_day,d):d for d in nights}
        for i,f in enumerate(as_completed(futures),1):
            value=f.result()
            if i%10==0 or i==len(nights):print(json.dumps(dict(prepared=i,total=len(nights),last=value)),flush=True)
    if prepare_only:
        write_json(Path(root)/'prepare_complete.json',dict(nights=len(nights),policy_hash=CTX.policy_hash))
        return
    jobs=[('clear_baseline',0,nights)]+[(s,i,nights) for s,n in [('baseline',16),('optimistic',8),('conservative',8)] for i in range(n)]
    if day_limit:jobs=[('clear_baseline',0,nights),('baseline',0,nights)]
    with ProcessPoolExecutor(max_workers=workers,mp_context=ctx) as pool:
        results=[f.result() for f in as_completed([pool.submit(simulate_case,j) for j in jobs])]
    write_json(Path(root)/'case_summaries.json',results)
    with ProcessPoolExecutor(max_workers=workers,mp_context=ctx) as pool:
        period_results=[f.result() for f in as_completed([pool.submit(period_case,r['case']) for r in results])]
    write_json(Path(root)/'period_case_summaries.json',period_results)
    print(json.dumps(dict(complete=True,cases=len(results),nights=len(nights))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--workers',type=int,default=16);p.add_argument('--days',type=int);p.add_argument('--prepare-only',action='store_true')
    a=p.parse_args();run(a.root,a.workers,a.days,a.prepare_only)
