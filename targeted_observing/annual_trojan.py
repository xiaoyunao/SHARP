"""Rolling multi-field annual Trojan plans using topocentric cached ephemerides.

Science campaigns last 35 days; pointing versions last seven days. Scheduling
uses accumulated nominal observations, never a future weather realization.
"""
import argparse
import copy
import hashlib
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import date, timedelta
import numpy as np
from astropy.time import Time
from .catalog import BASE, interpolate, write_csv
from .annual_sources import dates
from .geometry import Sky, setup_iers, unix, coord, in_field
from .scheduler import field, trojan
from .validate import validate_plan
from .science import load_periods, assess
from .bright_stars import load as load_stars


def load_tracks(limit=None):
    root=BASE/'products/annual';index=root/'ephemeris_index.json'
    if limit:
        paths=sorted((BASE/'cache/annual').glob('horizons_*_6h.json'))[:limit]
    else:
        manifest=json.loads(index.read_text())
        if manifest['failures']:raise ValueError('Annual ephemeris acquisition is incomplete')
        paths=map(__import__('pathlib').Path,manifest['paths'].values())
    tracks={}
    for path in paths:
        t=json.loads(path.read_text());tracks[t['pdes']]=t
    if not tracks:raise ValueError('No annual ephemerides')
    return tracks


def snapshot(tracks,seconds):
    ids=list(tracks);jd=Time(seconds,format='unix').jd
    ra=[];dec=[];mag=[]
    for k in ids:
        t=tracks[k];r,d=interpolate(t,seconds)
        ra.append(float(r[0]));dec.append(float(d[0]));mag.append(float(np.interp(jd,t['datetime_jd'],t['V'])))
    return ids,np.array(ra),np.array(dec),np.array(mag)


def campaign_groups(tracks,start,stop):
    """Candidate-centred overlapping footprints, deduplicated by member set.

    No imposed two-cloud season calendar. Actual ephemeris and joint visibility
    determine whether a candidate group can be scheduled on each night.
    """
    midpoint=(unix(start)+unix(stop))/2
    ids,ra,dec,mag=snapshot(tracks,midpoint)
    eligible=(mag<=19.2)&(dec>=-26)&np.isfinite(mag)
    groups={}
    for n in np.flatnonzero(eligible):
        members=np.flatnonzero(eligible&in_field(ra,dec,ra[n],dec[n],3,.3))
        keys=tuple(sorted(ids[i] for i in members))
        if keys not in groups:groups[keys]={'anchor':ids[n],'members':list(keys)}
    return list(groups.values())


def version_fields(groups,tracks,start,stop,stars):
    middle=(unix(start)+unix(stop))/2
    ids,ra,dec,mag=snapshot(tracks,middle);idx={k:i for i,k in enumerate(ids)}
    star_coord=coord([s['ra_deg'] for s in stars],[s['dec_deg'] for s in stars]);star_mag=np.array([s['mag'] for s in stars])
    result=[]
    for group in groups:
        i=idx[group['anchor']];f=field(ra[i],dec[i],'AY',start)
        f['campaign_members']=group['members'];f['anchor_object']=group['anchor']
        members=np.flatnonzero((mag<=19.2)&in_field(ra,dec,ra[i],dec[i],3,.15))
        f['core_targets']=[ids[j] for j in members]
        separations=coord(f['ra_deg'],f['dec_deg']).separation(star_coord).deg
        # Conservative scenario exclusion, not a measured detector halo model.
        f['very_bright_star_risk']=bool(np.any((star_mag<2)&(separations<4.122)))
        f['bright_stars_within_cap']=int(np.sum(separations<2.122))
        if f['core_targets']:result.append(f)
    return result


def priority(f,tracks,when,history,science):
    if f['very_bright_star_risk']:return 0.
    jd=Time(when,format='unix').jd;score=0.
    for k in f['core_targets']:
        past=history[k];nights=len(past)
        baseline=max(past)-min(past) if past else 0
        if nights>=14 and baseline>=28:continue
        v=float(np.interp(jd,tracks[k]['datetime_jd'],tracks[k]['V']))
        # Brightness smoothly weights expected useful photometry; no hard claim of r calibration.
        utility=1/(1+10**(.8*(v-18.5)))
        role=science.get(k,{}).get('role','literature_unchecked')
        science_factor=1.2 if role=='alias_review' else .8 if role=='control' else 1.
        continuity=1.3 if past else 1.
        need=.25 if nights>=14 else 1.
        score+=utility*science_factor*continuity*need
    return score


def choose_block(sky,night,fields,tracks,count,history,science,night_index):
    a,b=night['dark'];duration=sky.config['trojan']['block_s']
    if b-a<duration:return None,'insufficient dark time for requested allocation'
    weights=np.array([priority(f,tracks,(a+b)/2,history,science) for f in fields])
    candidates=[fields[i] for i in np.flatnonzero(weights>0)]
    weights=weights[weights>0]
    if not candidates:return None,'no unfinished non-excluded campaign candidates'
    starts=np.arange(a,b-duration+.01,600)
    ts=np.column_stack([starts+70,starts+duration-300]).ravel()
    m=sky.metrics([f['ra_deg'] for f in candidates],[f['dec_deg'] for f in candidates],ts)
    alt=m['alt'].reshape(len(candidates),-1,2);moon=m['moon_sep'].reshape(len(candidates),-1,2)
    lim=np.where(m['moon_fraction']>=.75,40,30).reshape(-1,2)
    valid=(alt.min(axis=2)>33)&np.all((moon>lim[None,:,:]+2.2)|(m['moon_alt'].reshape(-1,2)[None,:,:]<0),axis=2)
    choices=[]
    for j,start in enumerate(starts):
        available=list(np.flatnonzero(valid[:,j]));selected=[];used=set();total=0.
        for _ in range(count):
            if not available:break
            # Overlapping fields never get credit for the same object twice.
            def marginal(i):
                members=set(candidates[i]['core_targets'])
                return weights[i]*len(members-used)/len(members)
            i=max(available,key=marginal);gain=marginal(i)
            if gain<=0:break
            selected.append(candidates[i]);used.update(candidates[i]['core_targets']);total+=gain;available.remove(i)
        if selected:choices.append((total,len(selected),float(alt[:,j].mean()),start,selected))
    # Rotate among equally useful hour angles to avoid locking every night to one phase.
    choices.sort(key=lambda x:(-x[0],-x[1],abs((x[3]-a)/60-(night_index%3)*40-60),-x[2]))
    for _,_,_,start,selected in choices[:8]:
        block=trojan(sky,(start,start+duration),selected,tracks)
        if block:return block,None
    return None,'no shared full-footprint/moon-safe block among campaign candidates'


def run_strategy(job):
    count,dwell,spec,limit,day_limit=job
    tracks=load_tracks(limit);evidence=load_periods(BASE/'cache');science={k:assess(evidence.get(k)) for k in tracks}
    stars=load_stars(epoch=2027.3);config=json.loads((BASE/'defaults.json').read_text());config['trojan']['exposures_per_dwell']=dwell
    provenance=setup_iers(BASE/'cache'/spec['iers_file']);sky=Sky(config)
    nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    if day_limit:nights=nights[:day_limit]
    name=f'trojan_{count}fields_dwell{dwell}';root=BASE/'products/annual'/name
    if limit:root=root.with_name(name+'_smoke')
    root.mkdir(parents=True,exist_ok=True);history=defaultdict(set);ledger=[];all_fields={};campaigns=[]
    for i,day in enumerate(nights):
        if i%spec['campaign_days']==0:
            stop=min(date.fromisoformat(day)+timedelta(days=spec['campaign_days']),date.fromisoformat(spec['stop_date_exclusive']))
            groups=campaign_groups(tracks,day,stop.isoformat());campaign_start=day;campaign_stop=stop.isoformat()
            campaigns.append({'start_date':day,'stop_date_exclusive':campaign_stop,'candidate_groups':len(groups)})
        if (i%spec['campaign_days'])%spec['field_version_days']==0:
            version_start=day;version_stop=min(date.fromisoformat(day)+timedelta(days=spec['field_version_days']),date.fromisoformat(campaign_stop)).isoformat()
            fields=version_fields(groups,tracks,day,version_stop,stars)
            all_fields.update({f['field_id']:f for f in fields})
        night=sky.night(day)
        block,reason=choose_block(sky,night,fields,tracks,count,history,science,i)
        plan={'date':day,'config':config,'blocks':[block] if block else [],'night':night,'skips':[reason] if reason else [],
              'status':'LOCAL_ANNUAL_NOMINAL_NOT_FOR_TELESCOPE','campaign_start':campaign_start,'campaign_stop':campaign_stop,
              'field_version_start':version_start,'field_version_stop':version_stop,'iers':provenance,
              'nominal_feedback':'scheduled nights only; weather simulation must recompute object completion'}
        campaign={'start_date':version_start,'stop_date':version_stop,'fields':fields,'tracks':tracks}
        # Per-night target flags may fall below 80%; accumulate over full pointing version
        # in the final retention audit rather than silently interpreting nightly loss as full-campaign failure.
        audit=validate_plan(plan,campaign)
        geometric_errors=[e for e in audit['errors'] if 'planned-visit retention below threshold' not in e]
        if geometric_errors:raise RuntimeError(json.dumps(geometric_errors))
        plan['validation']=audit
        (root/f'{day}.json').write_text(json.dumps(plan)+'\n')
        rows=block['rows'] if block else []
        by_object=defaultdict(list)
        for row in rows:
            for k,state in row['target_states'].items():
                if state['in_field']:by_object[k].append(unix(row['mid_utc']))
        for k,ts in by_object.items():
            if len(ts)>=5 and max(ts)-min(ts)>=1800:history[k].add(i)
        ledger.append({'date':day,'strategy':name,'exposures':len(rows),'fields':len(block['fields']) if block else 0,
                       'allocated_s':unix(block['hard_deadline_utc'])-unix(block['not_before_utc']) if block else 0,
                       'exposure_s':sum(r['exposure_s'] for r in rows),'readout_s':sum(r['readout_s'] for r in rows),
                       'slew_s':sum(r['slew_s'] for r in rows),'filter_s':sum(r['filter_change_s'] for r in rows),
                       'sampled_objects':len(by_object),'zero_reason':reason or ''})
        if (i+1)%14==0:print(f'{name}: {i+1}/{len(nights)} nights',flush=True)
    # Final core/secondary labels use all actual planned visits within a pointing
    # version. Demotion preserves every exposure and target flag for science audit.
    from .retention import report
    version_plans=defaultdict(list)
    for day in nights:
        p=json.loads((root/f'{day}.json').read_text());version_plans[p['field_version_start']].append(p)
    retention_audit=[]
    for version,plans in version_plans.items():
        fs={f['field_id']:copy.deepcopy(f) for p in plans for b in p['blocks'] for f in b['fields']}
        by_night={p['date']:[r for b in p['blocks'] for r in b['rows']] for p in plans}
        sampling=report(list(fs.values()),tracks,by_night,config)
        for f in fs.values():
            secondary=[k for k in f['core_targets'] if sampling[f['field_id']+':'+k]['geometric_retention']<.8]
            f['secondary_targets']=secondary;f['core_targets']=[k for k in f['core_targets'] if k not in secondary]
            all_fields[f['field_id']]=f
        for key,row in sampling.items():
            retention_audit.append({'version':version,'object_field_key':key,'geometric_retention':row['geometric_retention'],
                                    'role':'core' if row['geometric_retention']>=.8 else 'secondary'})
        for p in plans:
            for b in p['blocks']:b['fields']=[fs[f['field_id']] for f in b['fields']]
            campaign={'start_date':p['field_version_start'],'stop_date':p['field_version_stop'],
                      'fields':list(fs.values()),'tracks':tracks,'planned_visits':by_night}
            p['validation']=validate_plan(p,campaign)
            if not p['validation']['passed']:raise RuntimeError(json.dumps(p['validation']))
            (root/f"{p['date']}.json").write_text(json.dumps(p)+'\n')
    write_csv(root/'retention_audit.csv',retention_audit)
    write_csv(root/'time_ledger.csv',ledger)
    (root/'fields.json').write_text(json.dumps(list(all_fields.values()))+'\n')
    (root/'campaigns.json').write_text(json.dumps(campaigns,indent=2)+'\n')
    (root/'object_nominal_nights.json').write_text(json.dumps({k:sorted(v) for k,v in history.items()})+'\n')
    (root/'science_evidence.json').write_text(json.dumps(science,indent=2)+'\n')
    return {'strategy':name,'nights':len(nights),'exposures':sum(r['exposures'] for r in ledger),'sampled_objects':len(history),'smoke':bool(limit)}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);parser.add_argument('--limit',type=int)
    parser.add_argument('--days',type=int);parser.add_argument('--fields',type=int,nargs='+');args=parser.parse_args()
    spec=json.loads((BASE/'annual_defaults.json').read_text())
    jobs=[(n,d,spec,args.limit,args.days) for n in args.fields or spec['field_counts'] for d in spec['dwells']]
    results=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(run_strategy,jobs):results.append(result);print(result,flush=True)
    suffix='_smoke' if args.limit else ''
    (BASE/f'products/annual/trojan_completeness{suffix}.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
