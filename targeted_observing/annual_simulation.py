"""Causal annual Trojan queue, shared weather and point-level quality scenarios."""
import argparse
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from .catalog import BASE,write_csv
from .annual_sources import dates
from .annual_trojan import load_tracks
from .annual_weather import ensemble,random_error,lightcurve_quality
from .geometry import in_field
from .science import load_periods,assess


def transition(a,b,c):
    if a is None:return c['far_slew_s']
    if a['field_id']==b['field_id']:return 0.
    ra1,ra2=np.deg2rad([a['ra_deg'],b['ra_deg']]);d1,d2=np.deg2rad([a['dec_deg'],b['dec_deg']])
    sep=np.rad2deg(np.arccos(np.clip(np.sin(d1)*np.sin(d2)+np.cos(d1)*np.cos(d2)*np.cos(ra1-ra2),-1,1)))
    nominal=c['adjacent_slew_s'] if sep<=c['adjacent_max_deg'] else c['far_slew_s']
    if not c.get('use_header_slew_floor',True):return nominal
    dra=abs((a['ra_deg']-b['ra_deg']+180)%360-180);ddec=abs(a['dec_deg']-b['dec_deg'])
    return max(nominal,9.24145+max(11.92277*np.sqrt(dra),12.22922*np.sqrt(ddec))+2)


def numeric_visits(fields,dwell,c,delay=0.,extra_readout=0.,extra_slew=0.):
    """Relative midpoints, field indices and actual overheads; complete cycles only."""
    end=c['trojan']['block_s']-c['trojan']['handover_s'];t=delay;previous=None;rows=[]
    route=[j for j in range(len(fields)) for _ in range(dwell)]
    if not route:return np.empty((0,5))
    for _ in range(100):
        cycle=[]
        for j in route:
            slew=transition(previous,fields[j],c)
            if slew:slew+=extra_slew
            change=c['filter_change_s'] if previous is None else 0
            begin=t+slew+change;release=begin+c['trojan']['exposure_s']+c['readout_s']+extra_readout
            if release>end:return np.array(rows).reshape(-1,5)
            cycle.append([j,begin+c['trojan']['exposure_s']/2,release,slew,change])
            t=release;previous=fields[j]
        rows.extend(cycle)
    raise RuntimeError('Unexpected unbounded visit loop')


def update_quality(state,thresholds):
    """Additional noisy companion frames cannot erase an accepted lightcurve.

    Keep all raw simulated points, but preserve the qualifying prefix and its
    certificate. Real calibration revisions would require a separate review.
    """
    if state.get('quality',{}).get('quality_lightcurve'):return
    state['quality']=lightcurve_quality(state['times'],state['night_indices'],state['errors'],
                                       state['zeropoints'],state['valid'],thresholds)
    if state['quality']['quality_lightcurve']:state['completion_index']=len(state['times'])


def pick(menu,data,values,count,history,night_index,completion=None):
    fields=menu['fields'];valid=data['valid']
    if not valid.any():return None
    scores=np.zeros(len(fields))
    def value(k):
        state=history.get(k);q=state.get('quality',{}) if state else {}
        if q.get('quality_lightcurve'):return 0.
        if completion is None:return values.get(k,0)*(1.7 if state and state['times'] else 1.)
        remaining,last,expected_points=completion
        first=q.get('first_good_night')
        if first is None:first=state['night_indices'][0] if state and state['times'] else night_index
        if last.get(k,-1)-first<21:return 0.
        need=max(1.,10-q.get('good_nights',0),(100-q.get('good_points',0))/expected_points)
        feasible=min(1.,remaining.get(k,0)/need)
        continuity=1+min(1.,q.get('good_nights',0)/10)
        # Once point and night requirements are met, let new targets use the
        # waiting interval; return after the 21-day baseline can actually be met.
        if q.get('good_points',0)>=100 and q.get('good_nights',0)>=10 and night_index-first<21:
            continuity*=.05
        return values.get(k,0)*feasible*continuity
    active=np.flatnonzero(np.any(valid,axis=1))
    object_values={k:value(k) for j in active for k in fields[j]['core_targets']}
    for j in active:
        f=fields[j]
        for k in f['core_targets']:
            scores[j]+=object_values.get(k,0)
    available=np.where(valid,scores[:,None],0)
    top=np.sort(available,axis=0)[-count:].sum(axis=0)
    order=np.argsort(-top,kind='stable')
    # Prefer changing hour angle among near-equal utility starts.
    if len(order):
        near=order[top[order]>=top[order[0]]*.98]
        if len(near):order=np.r_[np.roll(near,night_index%3*3),order[~np.isin(order,near)]]
    for start_index in order:
        if top[start_index]<=0:break
        options=list(np.flatnonzero(valid[:,start_index]&(scores>0)))
        options=sorted(options,key=lambda j:-scores[j])[:max(32,count*8)]
        chosen=[];used=set()
        for _ in range(count):
            if not options:break
            def gain(j):return sum(object_values.get(k,0) for k in fields[j]['core_targets'] if k not in used)
            j=max(options,key=gain)
            if gain(j)<=0:break
            chosen.append(j);used.update(fields[j]['core_targets']);options.remove(j)
        if chosen:return int(start_index),chosen
    return None


def simulate(job):
    field_count,dwell,scenario_name,spec,day_limit,realization_limit,smoke=job
    c=json.loads((BASE/'defaults.json').read_text());scenario=dict(spec['photometry_scenarios'][scenario_name])
    nominal=spec.get('nominal_clear_weather',False)
    if nominal:scenario['usable_fraction']=1.
    duration=spec.get('experimental_duration_s',7200)
    c['trojan']['block_s']=duration
    tracks=load_tracks(60 if smoke else None)
    # Unix-seconds arrays avoid millions of repeated Time constructions; no leap
    # second occurs inside the frozen UTC table used by this planning experiment.
    track_arrays={k:{'t':(np.array(v['datetime_jd'])-2440587.5)*86400,'ra':np.unwrap(np.deg2rad(v['RA'])),
                     'dec':np.array(v['DEC']),'V':np.array(v['V']),'delta':np.array(v['delta'])} for k,v in tracks.items()}
    nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    if day_limit:nights=nights[:day_limit]
    nreal=1 if nominal else realization_limit or spec['weather_realizations'];weather=ensemble(nights,nreal,spec['seed'])
    if nominal:
        weather['observable'].fill(True);weather['reserved'].fill(False);weather['seeing_arcsec'].fill(2.)
        for key in ['cloud_duration_fraction','delay_s','extra_readout_s','extra_slew_s','zeropoint_standard_normal']:weather[key].fill(0)
    vroot=BASE/'products/annual'/('visibility_smoke' if smoke else 'visibility')
    menus=[json.loads((vroot/f'{day}.json').read_text()) for day in nights]
    arrays=[dict(np.load(vroot/f'{day}.npz')) for day in nights]
    if duration!=7200:
        for menu,data in zip(menus,arrays):
            grid=data['grid'];starts=grid[grid+duration<=menu['night']['dark'][1]+.01]
            valid=np.zeros((len(menu['fields']),len(starts)),bool)
            if len(starts) and 'alt' in data:
                guard=.0043*150+.01;altlim=max(c['center_alt_deg'],c['min_alt_deg']+2.122)+guard
                moonlim=np.where(data['moon_fraction']>=c['bright_moon_fraction']-.001,c['bright_moon_sep_deg'],c['moon_sep_deg'])+2.122+guard
                good=(data['alt']>=altlim)&((data['moon_sep']>=moonlim[None,:])|(data['moon_alt'][None,:]<-guard))
                width=int(duration/300)+1
                for j in range(len(starts)):valid[:,j]=np.all(good[:,j:j+width],axis=1)
                for j,f in enumerate(menu['fields']):
                    if f['very_bright_star_risk']:valid[j]=False
            data['starts']=starts;data['valid']=valid
    night_values=[]
    strength=spec.get('science_priority_strength',.25)
    science={k:assess(v) for k,v in load_periods(BASE/'cache').items()}
    for menu in menus:
        when=sum(menu['night']['dark'])/2;keys=list(track_arrays)
        magnitudes=np.array([np.interp(when,track_arrays[k]['t'],track_arrays[k]['V']) for k in keys])
        expected=random_error(magnitudes,1.3,0,90,2,120,scenario)
        weights=np.array([science.get(k,{}).get('weight',.1) for k in keys])
        night_values.append(dict(zip(keys,((1-strength+strength*weights)/(1+(expected/.03)**4)).tolist())))
    completion_forecast=None
    if strength:
        ids=list(track_arrays);index={k:j for j,k in enumerate(ids)};possible=np.zeros((len(nights),len(ids)),bool)
        for ni,(menu,data) in enumerate(zip(menus,arrays)):
            for j in np.flatnonzero(np.any(data['valid'],axis=1)):
                for k in menu['fields'][j]['core_targets']:possible[ni,index[k]]=True
        last={k:int(np.flatnonzero(possible[:,j])[-1]) if possible[:,j].any() else -1 for j,k in enumerate(ids)}
        prior=json.loads((BASE/'data/weather_prior.json').read_text())
        monthly=np.mean(list(prior['observable_percent_by_year'].values()),axis=0)/100
        availability=np.array([monthly[int(day[5:7])-1] for day in nights])*.85*scenario['usable_fraction']
        remaining=np.cumsum((possible*availability[:,None])[::-1],axis=0)[::-1]
        expected_points=max(1,int((duration-370)/(field_count*(135+60/dwell)))*scenario['usable_fraction'])
        completion_forecast=[(dict(zip(ids,row)),last,expected_points) for row in remaining]
    name=f'queue_{field_count}fields_dwell{dwell}_{scenario_name}'+('_smoke' if smoke else '')
    if duration!=7200:name+=f'_budget{duration}'
    if strength:name+=f'_science{int(strength*100):02d}'
    if strength and spec.get('queue_policy_version',1)>1:name+=f'_v{spec["queue_policy_version"]}'
    if nominal:name='nominal_'+name
    root=BASE/'products/annual'/name;root.mkdir(parents=True,exist_ok=True)
    nightly=[];objects=[];samples={};actual_rows=[];totals=[]
    curve_times=[];curve_errors=[];curve_zeropoints=[];curve_offsets=[0];curve_ids=[];curve_realizations=[]
    for realization in range(nreal):
        history={};cumulative=0
        for ni,(day,menu,data) in enumerate(zip(nights,menus,arrays)):
            w=lambda key:weather[key][realization,ni,1]
            record={'realization':realization,'date':day,'fields':0,'exposures':0,'allocated_s':0,'exposure_s':0,
                    'readout_s':0.,'slew_s':0.,'filter_s':0.,'delay_s':0.,'handover_s':0.,'slack_s':0.,
                    'quality_objects_cumulative':sum(s.get('quality',{}).get('quality_lightcurve',False) for s in history.values()),'reason':''}
            if not w('observable'):record['reason']='weather_closed';nightly.append(record);continue
            if w('reserved'):record['reason']='external_allocation_cancelled';nightly.append(record);continue
            chosen=pick(menu,data,night_values[ni],field_count,history,ni,completion_forecast[ni] if completion_forecast else None)
            if chosen is None:record['reason']='no_eligible_unfinished_targets';nightly.append(record);continue
            slot,indices=chosen;fields=[menu['fields'][j] for j in indices];start=float(data['starts'][slot])
            rows=numeric_visits(fields,dwell,c,float(w('delay_s')),float(w('extra_readout_s')),float(w('extra_slew_s')))
            if not len(rows):record['reason']='no_complete_cycle_after_delay';nightly.append(record);continue
            record.update(fields=len(fields),exposures=len(rows),allocated_s=duration,exposure_s=120*len(rows),
                readout_s=(c['readout_s']+float(w('extra_readout_s')))*len(rows),slew_s=float(rows[:,3].sum()),
                filter_s=float(rows[:,4].sum()),delay_s=float(w('delay_s')),handover_s=c['trojan']['handover_s'])
            record['slack_s']=duration-sum(record[k] for k in ['exposure_s','readout_s','slew_s','filter_s','delay_s','handover_s'])
            rng=np.random.default_rng(int(w('point_seed')));touched=set()
            zp=float(weather['zeropoint_standard_normal'][realization,ni])*scenario['night_zeropoint']
            for j,f in enumerate(fields):
                rr=rows[rows[:,0]==j];seconds=start+rr[:,1];grid=data['grid'];fi=indices[j]
                altitude=np.interp(seconds,grid,data['alt'][fi]);airmass=1/np.sin(np.deg2rad(altitude))
                moonsep=np.interp(seconds,grid,data['moon_sep'][fi]);moonfrac=np.interp(seconds,grid,data['moon_fraction'])
                moonfrac=np.where(np.interp(seconds,grid,data['moon_alt'])>0,moonfrac,0)
                cloudy=(rr[:,1]/duration>=w('cloud_start_fraction'))&(rr[:,1]/duration<w('cloud_start_fraction')+w('cloud_duration_fraction'))
                for k in f['core_targets']:
                    t=track_arrays[k];ra=np.rad2deg(np.interp(seconds,t['t'],t['ra']))%360;dec=np.interp(seconds,t['t'],t['dec']);mag=np.interp(seconds,t['t'],t['V'])
                    inside=in_field(ra,dec,f['ra_deg'],f['dec_deg'],3,.15)
                    errors=random_error(mag,airmass,moonfrac,moonsep,w('seeing_arcsec'),120,scenario)
                    # Bad-pixel/cosmic/cloud losses are scenarios; close bright-star
                    # field exclusions have already been applied in the visibility menu.
                    valid=inside&~cloudy&(rng.random(len(rr))<scenario['usable_fraction'])
                    state=history.setdefault(k,{'times':[],'night_indices':[],'errors':[],'zeropoints':[],'valid':[],'magnitudes':[]})
                    state['times'].extend(seconds.tolist());state['night_indices'].extend([ni]*len(rr));state['errors'].extend(errors.tolist())
                    state['zeropoints'].extend([zp]*len(rr));state['valid'].extend(valid.tolist());state['magnitudes'].extend(mag.tolist());touched.add(k)
                if realization==0:
                    for row in rr:actual_rows.append({'date':day,'field_id':f['field_id'],'ra_deg':f['ra_deg'],'dec_deg':f['dec_deg'],
                        'start_unix':start+row[1]-60,'mid_unix':start+row[1],'end_unix':start+row[1]+60,
                        'release_unix':start+row[2],'exposure_s':120,'readout_s':c['readout_s']+float(w('extra_readout_s')),
                        'slew_s':row[3],'filter_change_s':row[4],
                        'block_start_unix':start,'block_deadline_unix':start+duration})
            for k in touched:
                update_quality(history[k],spec['quality'])
            record['quality_objects_cumulative']=sum(s.get('quality',{}).get('quality_lightcurve',False) for s in history.values())
            nightly.append(record)
        qualified=0
        for k,s in history.items():
            q=s['quality'];qualified+=q['quality_lightcurve'];objects.append({'realization':realization,'object_id':k,
                'median_V_proxy':float(np.median(s['magnitudes'])),
                'raw_points_total':len(s['times']),'completion_index':s.get('completion_index'),**q})
            if spec.get('save_qualified_curves') and q['quality_lightcurve']:
                good=np.array(s['valid'])&(np.array(s['errors'])<=.05)&(abs(np.array(s['zeropoints']))<=.02)
                seconds=np.array(s['times'])[good];tr=track_arrays[k]
                emission=seconds-np.interp(seconds,tr['t'],tr['delta'])*499.004783836
                curve_times.extend(emission.tolist());curve_errors.extend(np.array(s['errors'])[good].tolist())
                curve_zeropoints.extend(np.array(s['zeropoints'])[good].tolist())
                curve_offsets.append(len(curve_times));curve_ids.append(k);curve_realizations.append(realization)
        totals.append(qualified)
        if realization==0:samples=history
        if (realization+1)%20==0:print(f'{name}: {realization+1}/{nreal} weather realizations',flush=True)
    actual_rows.sort(key=lambda r:r['start_unix']);counters=defaultdict(int)
    for row in actual_rows:
        counters[row['date']]+=1;row['exposure_id']=f"{row['date']}_TJ_{counters[row['date']]:04d}"
    write_csv(root/'nightly_execution_ledger.csv',nightly);write_csv(root/'object_yield.csv',objects);write_csv(root/'execution_realization0.csv',actual_rows)
    (root/'samples_realization0.json').write_text(json.dumps(samples)+'\n')
    (root/'science_roles.json').write_text(json.dumps({k:science.get(k,assess(None)) for k in tracks},indent=2)+'\n')
    if spec.get('save_qualified_curves'):
        np.savez_compressed(root/'qualified_curves.npz',emission_unix=np.array(curve_times),error=np.array(curve_errors),
            zeropoint=np.array(curve_zeropoints),offsets=np.array(curve_offsets),object_id=np.array(curve_ids),realization=np.array(curve_realizations))
    result={'strategy':name,'nights':len(nights),'realizations':nreal,'quality_count_P10_P50_P90':np.percentile(totals,[10,50,90]).tolist(),
            'probability_ge101':float(np.mean(np.array(totals)>=101)),'counts':totals,
            'scope':'causal target completion; observed current-night closures; no future weather in selection',
            'scenario':scenario,'smoke':smoke,'period_recovery_completed':False,
            'moon_policy':menus[0]['moon_policy'],'qualified_curves_saved':bool(spec.get('save_qualified_curves')),
            'allocation_s':duration,'allocation_scope':'standard <=2h mode' if duration==7200 else 'additional-time proposal experiment, not dispatchable by standard mode'}
    result['science_priority_strength']=strength
    result['queue_policy_version']=spec.get('queue_policy_version',1)
    result['nominal_clear_weather']=nominal
    result['selection_policy']='photometry * (1-strength+strength*literature_weight) * expected completion availability; past outcomes only' if strength else 'photometry and observed continuity; throughput comparison'
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--days',type=int)
    p.add_argument('--realizations',type=int);p.add_argument('--fields',type=int,nargs='+');p.add_argument('--scenario',nargs='+')
    p.add_argument('--dwells',type=int,nargs='+',choices=[1,2])
    p.add_argument('--duration-s',type=int,choices=[7200,10800,14400],default=7200)
    p.add_argument('--science-strength',type=float)
    p.add_argument('--nominal',action='store_true')
    p.add_argument('--smoke',action='store_true');p.add_argument('--save-curves',action='store_true')
    args=p.parse_args();spec=json.loads((BASE/'annual_defaults.json').read_text());spec['save_qualified_curves']=args.save_curves
    spec['experimental_duration_s']=args.duration_s
    spec['nominal_clear_weather']=args.nominal
    if args.science_strength is not None:spec['science_priority_strength']=args.science_strength
    if not 0<=spec['science_priority_strength']<=1:raise ValueError('Science strength must lie in [0,1]')
    jobs=[(n,d,s,spec,args.days,args.realizations,args.smoke) for n in args.fields or spec['field_counts'] for d in args.dwells or spec['dwells']
          for s in args.scenario or list(spec['photometry_scenarios'])]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(simulate,jobs):print({k:v for k,v in result.items() if k!='counts'},flush=True)


if __name__=='__main__':main()
