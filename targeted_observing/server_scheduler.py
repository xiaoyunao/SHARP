"""Daily shared-grid scheduler, explicit website reservations and actual TJN history."""
from collections import defaultdict
from datetime import datetime
import json
import math
from itertools import permutations
from pathlib import Path

import numpy as np
from astropy.time import Time

from .annual_trojan import load_tracks
from .bright_stars import load as load_stars
from .catalog import BASE, interpolate
from .full_night_test import command, movement, YSO_FIELDS
from .geometry import coord, setup_iers, local_time
from .offline_twilight import load_config as twilight_config, make_plan, make_sky, exposure
from .server_common import TZ
from .server_history import states
from .shared_fields import SharedFields, cloud_labels


def subtract(window, reservations):
    free = [tuple(window)]
    for a,b in sorted(reservations):
        result = []
        for x,y in free:
            if a>=y or b<=x: result.append((x,y)); continue
            if x<a: result.append((x,a))
            if b<y: result.append((b,y))
        free=result
    return free


class Geometry:
    """60s lookup for selection; every final exposure is independently rechecked."""
    def __init__(self, sky, fields, start, end, radius):
        self.fields=fields; self.index={f['field_id']:i for i,f in enumerate(fields)}
        self.times=np.arange(start-300,end+361,60.)
        self.metrics=sky.metrics([f['ra_deg'] for f in fields], [f['dec_deg'] for f in fields],self.times)
        lst=Time(self.times,format='unix').sidereal_time('apparent',longitude=sky.location.lon).deg
        self.ha=((lst[None,:]-np.array([f['ra_deg'] for f in fields])[:,None]+180)%360-180)/15
        self.radius=radius

    def check(self, field, start, end, altitude=30, sun=-18, elongation=False):
        if start<self.times[0] or end>self.times[-1]: return False
        i=self.index[field['field_id']]
        t=np.r_[start,self.times[(self.times>start)&(self.times<end)],end]
        m=self.metrics
        a=np.interp(t,self.times,m['alt'][i])-self.radius
        phase=np.interp(t,self.times,m['moon_fraction'])
        moon=np.interp(t,self.times,m['moon_sep'][i])-self.radius
        up=np.interp(t,self.times,m['moon_alt'])
        solar=np.interp(t,self.times,m['sun_alt'])
        good=(a>altitude+.15)&((up<-.15)|(moon>np.where(phase>=.75,40,30)+.15))&(solar<=sun+.001)
        good &= abs(np.interp(t,self.times,self.ha[i]))<4.95
        if elongation: good &= np.interp(t,self.times,m['elong'][i])<=90
        return bool(good.all())


class Timeline:
    def __init__(self,start,config,rows=None):
        self.config=config;self.start=start;self.rows=list(rows or []);self.gaps=[]
        self.protected=len(self.rows)
        self.time=self.rows[-1]['release_unix'] if self.rows else start

    @property
    def previous(self): return self.rows[-1] if self.rows else None

    def overhead(self,field,band):
        previous=self.previous
        turn=90. if previous is None else movement(previous,field)
        change=self.config['filter_change_s'] if previous is None or previous['filter']!=band else 0.
        return turn,change

    def append(self,field,band,seconds,mode,targets=(),candidate=None,request_sha=None):
        turn,change=self.overhead(field,band);begin=self.time+turn+change
        row=dict(field_id=field['field_id'],ra_deg=field['ra_deg'],dec_deg=field['dec_deg'],
                 filter=band,exposure_s=int(seconds),mode=mode,block=mode,
                 target_ids=list(targets),candidate_id=candidate,
                 request_sha256=request_sha,
                 overhead_start_unix=self.time,slew_s=turn,filter_change_s=change,
                 readout_s=self.config['readout_s'],start_unix=begin,mid_unix=begin+seconds/2,
                 end_unix=begin+seconds,release_unix=begin+seconds+self.config['readout_s'])
        self.rows.append(row);self.time=row['release_unix'];return row

    def wait_to(self,time,reason):
        if time>self.time+.001:
            self.gaps.append(dict(start_unix=self.time,end_unix=time,reason=reason))
            self.time=time

    def fill_tail(self, release):
        """Fill a boundary on the current pointing, including a sub-readout tail."""
        while self.rows and release-self.time>=1:
            last=self.rows[-1];slack=int(math.floor(release-self.time))
            if last['mode'] in ('FOLLOW','YSO'):break
            if slack>=5+self.config['readout_s']:
                maximum=30 if last['mode'] in ('EV','AM') else 120
                self.append(last,last['filter'],min(maximum,slack-self.config['readout_s']),last['mode'],last['target_ids'])
            else:
                delta=5+self.config['readout_s']-slack
                if last['exposure_s']-delta<5 or len(self.rows)<=self.protected:
                    break
                last['exposure_s']-=delta;last['mid_unix']-=delta/2
                last['end_unix']-=delta;last['release_unix']-=delta;self.time-=delta
                self.append(last,last['filter'],5,last['mode'],last['target_ids'])


def object_priority(key, history, evidence, when, planned):
    state=history.get(key,{})
    if state.get('sampling_ready'): return 0.
    science=evidence.get(key,{})
    base=.8+.4*science.get('science_weight',.1)
    nights=state.get('good_nights',0); points=state.get('good_points',0)
    base*=1+min(nights/10,1)
    if points>=100 and nights>=10 and state.get('baseline_days',0)<state.get('required_baseline_days',21): base*=.1
    # Planned points only reduce the marginal reward in this plan; they never
    # enter the durable actual-exposure history or count as observations.
    ts=planned.get(key,[])
    if ts:
        span=max(ts)-min(ts)
        if len(ts)>=10 and span>=1800: base*=.2
        if when-ts[-1]<300: base*=.7
    last=state.get('last_mid_unix')
    periods=science.get('periods_h',[])
    if last and periods:
        # Prefer a new phase relative to the last measured epoch, including all
        # frozen alternatives; no invented period for an unmeasured object.
        phase=np.mod((when-last)/(3600*np.asarray(periods)),1)
        novel=[int((when/(3600*p)%1)*20) not in state.get('prior_phase_bins',{}).get(str(p),[]) for p in periods]
        base*=1+.3*float(np.mean(np.abs(np.sin(np.pi*phase))))+.6*float(np.mean(novel))
    return base


def select_fields(catalog,tracks,clouds,geometry,begin,end,history,evidence,planned,config,stars):
    times=np.linspace(begin,end,max(3,int((end-begin)/600)+1))
    members=catalog.memberships(tracks,times,clouds,config['tjn_mag_limit'],retained_fraction=1.)
    grouped={}
    for row in members:
        key=(row['cloud'],row['field_id'])
        f=grouped.setdefault(key,dict(field_id=row['execution_name'],shared_field_id=row['field_id'],
            geometry_id=row['field_id'],ra_deg=row['ra_deg'],dec_deg=row['dec_deg'],targets=[],values={}))
        f['targets'].append(row['object_id'])
        value=object_priority(row['object_id'],history,evidence,(begin+end)/2,planned)
        f['values'][row['object_id']]=value/(1+10**(.8*(row['V_max']-18.5)))
    viable=[]
    for f in grouped.values():
        if not geometry.check({**f,'field_id':f['geometry_id']},begin,end):continue
        if stars is not None and np.any(coord(f['ra_deg'],f['dec_deg']).separation(stars).deg<4.122):continue
        if sum(f['values'].values())>0:viable.append(f)
    chosen=[];used=set()
    for _ in range(config['tjn_field_count']):
        if not viable: break
        best=max(viable,key=lambda f:sum(v for k,v in f['values'].items() if k not in used))
        if not sum(v for k,v in best['values'].items() if k not in used)>0: break
        chosen.append(best);used.update(best['targets']);viable.remove(best)
    if len(chosen)>1:
        chosen=list(min(permutations(chosen),key=lambda fs:sum(movement(fs[i-1],f) for i,f in enumerate(fs))))
    return chosen


def final_validation(rows,sky,catalog,config,tracks):
    errors=[];retention=[];unique={r['field_id']:r for r in rows}
    fields=list(unique.values());indices={f['field_id']:i for i,f in enumerate(fields)}
    radius=float(catalog.twilight_grid()[0][0]['radius_deg']);seconds=[];codes=[];row_codes=[]
    for j,r in enumerate(rows):
        if j and r['overhead_start_unix']<rows[j-1]['release_unix']-.001:errors.append('overlapping exposures')
        if abs(r['mid_unix']-(r['start_unix']+r['exposure_s']/2))>.001:errors.append('invalid midpoint')
        if r['field_id'] not in YSO_FIELDS:
            f=catalog.resolve(r['field_id'])
            if coord(f['ra_deg'],f['dec_deg']).separation(coord(r['ra_deg'],r['dec_deg'])).arcsec>.03:
                errors.append('shared center changed: '+r['field_id'])
        ts=np.linspace(r['start_unix'],r['release_unix'],int(np.ceil((r['release_unix']-r['start_unix'])/5))+1)
        seconds.extend(ts);codes.extend([indices[r['field_id']]]*len(ts));row_codes.extend([j]*len(ts))
    if seconds:
        m=sky.metrics([f['ra_deg'] for f in fields],[f['dec_deg'] for f in fields],seconds)
        ix=np.array(codes);jj=np.arange(len(seconds));rc=np.array(row_codes)
        alt=np.array([25 if r['mode'] in ('EV','AM','FOLLOW') else 30 for r in rows])[rc]
        solar=np.array([-10 if r['mode'] in ('EV','AM','FOLLOW') else -18 for r in rows])[rc]
        if np.any(m['alt'][ix,jj]-radius<=alt+.011):errors.append('full-footprint altitude')
        if np.any(m['sun_alt']>solar+.003):errors.append('Sun altitude')
        if np.any((m['moon_alt']>-.011)&(m['moon_sep'][ix,jj]-radius<np.where(m['moon_fraction']>=.75,40,30)+.011)):
            errors.append('Moon clearance')
        from astropy.coordinates import HADec
        import astropy.units as u
        # This is a geometric HA envelope, not encoder/meridian-flip certification.
        hour=coord([r['ra_deg'] for r in rows],[r['dec_deg'] for r in rows]).transform_to(
            HADec(obstime=Time([r['mid_unix'] for r in rows],format='unix'),location=sky.location,pressure=0*u.hPa)).ha.wrap_at(180*u.deg).hour
        if max(abs(hour))>=5:errors.append('reference tracking hour-angle envelope')
    for fid in {r['field_id'] for r in rows if r['mode']=='TJN'}:
        rr=[r for r in rows if r['field_id']==fid and r['mode']=='TJN']
        for key in {k for r in rr for k in r['target_ids']}:
            ts=[r['mid_unix'] for r in rr if key in r['target_ids']]
            ra,dec=interpolate(tracks[key],ts)
            good=catalog.edge_distance(ra,dec,fid.rsplit('_',1)[1])>=.16
            retention.append(dict(field=fid,object_id=key,fraction=float(good.mean())))
            if good.mean()<.8:errors.append('TJN target retention '+key)
    mp_sampling=[]
    for side in ('EV','AM'):
        for fid in sorted({r['field_id'] for r in rows if r['mode']==side}):
            ts=[r['mid_unix'] for r in rows if r['mode']==side and r['field_id']==fid]
            gaps=np.diff(ts)
            mp_sampling.append(dict(side=side,field=fid,frames=len(ts),
                min_gap_s=float(gaps.min()) if len(gaps) else None,
                max_gap_s=float(gaps.max()) if len(gaps) else None,
                sampling_ok=bool(len(ts)>=3 and np.all((gaps>=300)&(gaps<=720)))))
    return dict(passed=not errors,errors=sorted(set(errors)),target_retention=retention,mp_sampling=mp_sampling,
                mechanical_safety_certified=False,sampling_is_not_period_solution=True)


def build(day,config,requests,now=None,previous_plan=None):
    catalog=SharedFields();c=twilight_config()
    setup_iers(BASE/'cache'/config['iers_file']);sky=make_sky()
    twilight=make_plan(day,c,grid=catalog.twilight_grid(),allow_empty=True,shorter_fallback=True)
    begin=twilight['blocks'][0]['reserved_window_unix'][0];end=twilight['blocks'][1]['reserved_window_unix'][1]
    ev_end=begin+3600;am_start=end-3600
    tjn_end=min(end,twilight['dark_crossings'][1])
    tracks=load_tracks();clouds=cloud_labels()
    evidence=json.loads((BASE/'data/tjn_queue_evidence_20261118.json').read_text())['objects']
    history=states(config,{k:v['periods_h'] for k,v in evidence.items()})
    geometry=Geometry(sky,list(catalog.by_id.values()),begin,end,catalog.twilight_grid()[0][0]['radius_deg'])
    prefix=[]
    if now is not None and now>begin and previous_plan:
        cutoff=now
        active=[r for r in previous_plan['rows'] if r['mode']=='FOLLOW' and r['overhead_start_unix']<=now<r['release_unix']]
        if active:
            current=active[0]
            cutoff=max(r['release_unix'] for r in previous_plan['rows'] if r['candidate_id']==current['candidate_id']
                       and r.get('request_sha256')==current.get('request_sha256'))
        prefix=[r for r in previous_plan['rows'] if (r['overhead_start_unix']<cutoff-.001 if active else r['overhead_start_unix']<=cutoff)]
    cutoff=prefix[-1]['release_unix'] if prefix else begin
    for r in prefix:r.setdefault('request_sha256',None)
    statuses=[];events=[]
    for request in sorted(requests,key=lambda r:(r['start_time'],r['candidate_id'])):
        requested=datetime.fromisoformat(request['start_time']).timestamp()
        f={**catalog.by_id[request['field_id']],'field_id':request['execution_name']}
        duration=request['frames']*(request['exposure_s']+config['readout_s'])
        status=dict(candidate_id=request['candidate_id'],request_sha256=request.get('sha256'),requested_start=request['start_time'],frames=request['frames'],scheduled=False)
        prior=[r for r in prefix if r['candidate_id']==request['candidate_id'] and r.get('request_sha256')==request.get('sha256')]
        if prior:
            status['reason']='preserved prior execution block; no automatic repetition'
        else:
            max_delay=max(0,end-duration-requested)
            tolerance=config['followup_start_tolerance_s']
            offsets=np.r_[np.arange(0,min(tolerance,max_delay)+1,1.),
                          np.arange(tolerance+1,max_delay+1,30.)]
            for shift in offsets:
                start=requested+shift;finish=start+duration
                if now is not None and now>begin and start<max(cutoff+224,now+config['late_replan_lead_s']):continue
                # Reserve worst-case slew from the preceding target, then fill
                # the actual handover with legal exposures where possible.
                a=start-200-config['filter_change_s']
                if a<begin or finish>end or any(a<e['end'] and finish>e['setup'] for e in events):continue
                if not geometry.check({**f,'field_id':request['field_id']},start,finish,
                        config['followup_min_alt_deg'],-10):continue
                events.append(dict(mode='FOLLOW',field=f,start=start,setup=a,end=finish,
                                   request=request,frames=request['frames'],band=request['filter'],seconds=request['exposure_s']))
                status.update(scheduled=True,nominal_start=local_time(start),nominal_offset_s=float(shift),
                              preferred_time_met=bool(shift<=tolerance))
                break
            if not status['scheduled']:status['reason']='requested block conflicts with another follow-up or sky/time constraints'
        statuses.append(status)
    # Preserve the already independently verified full-hour twilight sequences.
    # Tail repeats on the same pointing would violate their 5-minute cadence.
    am_block=twilight['blocks'][1]
    am_fixed=bool(am_block['rows']) and not any(e['mode']=='FOLLOW' and e['setup']<end and e['end']>am_start-224 for e in events)
    if am_fixed:
        first=am_block['rows'][0]
        f={**first,'field_id':'MP_'+first['field_id'][2:]}
        events.append(dict(mode='AM',field=f,start=first['start_unix'],setup=first['start_unix']-224,
            end=am_block['rows'][-1]['release_unix'],sequence=am_block['rows'],band='r'))
    # Reserve the best continuous YSO block in the remaining dark intervals.
    yso_candidates=[]
    remaining={band:max(0,config['yso_frames_per_band']-sum(r['mode']=='YSO' and r['filter']==band for r in prefix)) for band in ('g','r')}
    yso_duration=sum(remaining.values())*(config['yso_exposure_s']+config['readout_s'])+(config['filter_change_s'] if all(remaining.values()) else 0)
    dark=(max(ev_end,twilight['dark_crossings'][0],cutoff),min(am_start,twilight['dark_crossings'][1]))
    yso_geometry=Geometry(sky,[dict(field_id=k,**v) for k,v in YSO_FIELDS.items()],begin,end,geometry.radius)
    for a,b in subtract(dark,[(e['setup'],e['end']) for e in events]):
        for name,values in YSO_FIELDS.items():
            if not yso_duration:continue
            previous_yso={r['field_id'] for r in prefix if r['mode']=='YSO'}
            if previous_yso and name not in previous_yso:continue
            f=dict(field_id=name,**values)
            for start in np.arange(a+224,b-yso_duration+1,300):
                if yso_geometry.check(f,start,start+yso_duration):
                    altitude=float(sky.metrics(f['ra_deg'],f['dec_deg'],[(start+start+yso_duration)/2])['alt'][0,0])
                    score=altitude+(100 if name==config['yso_preferred'] else 0)
                    yso_candidates.append((score,start,f))
    yso_status='no full three-hour YSO block fits after higher-priority followups'
    if not yso_duration:yso_status='completed in preserved execution prefix'
    if yso_candidates:
        _,start,f=max(yso_candidates,key=lambda x:x[0]);events.append(dict(mode='YSO',field=f,start=start,
            setup=start-224,end=start+yso_duration,frames=config['yso_frames_per_band'],remaining=remaining,
            band='g' if remaining['g'] else 'r',seconds=config['yso_exposure_s']))
        yso_status=f['field_id']
    events.sort(key=lambda e:e['setup'])
    # Preserve already released or in-progress rows of a previously published
    # plan during a late-night replan; no restart at the beginning of the night.
    if prefix:events=[e for e in events if e['setup']>=cutoff]
    timeline=Timeline(begin,config,prefix)
    planned=defaultdict(list)
    for r in prefix:
        for key in r['target_ids']:planned[key].append(r['mid_unix'])
    if now is not None and now>timeline.time and now>begin:
        timeline.wait_to(now,'replan begins after current time; past time is not rescheduled')
    ev_fixed=False
    if not prefix and timeline.time==begin and not any(e['mode']=='FOLLOW' and e['setup']<ev_end for e in events):
        if twilight['blocks'][0]['rows']:
            timeline.wait_to(twilight['blocks'][0]['rows'][0]['overhead_start_unix'],'shorter EV starts later')
        for r in twilight['blocks'][0]['rows']:
            timeline.append({**r,'field_id':'MP_'+r['field_id'][2:]},'r',r['exposure_s'],'EV')
        ev_fixed=bool(twilight['blocks'][0]['rows'])
    elif prefix and timeline.time>=ev_end-10:ev_fixed=True
    bright=[r for r in load_stars() if r['mag']<2]
    stars=coord([r['ra_deg'] for r in bright],[r['dec_deg'] for r in bright]) if bright else None
    twilight_routes={b['side']:[{**f,'field_id':'MP_'+f['field_id'][2:]} for f in b['fields']] for b in twilight['blocks']}
    choices=[];current=[];dwell=0;route_index=0;campaign_end=0;campaign_dwell=config['tjn_dwell_frames']

    def fill(deadline,next_event=None):
        nonlocal current,dwell,route_index,campaign_end,campaign_dwell
        def wait_until(target,reason):
            if next_event is not None:
                turn,change=timeline.overhead(next_event['field'],next_event['band'])
                target=min(target,next_event['start']-turn-change)
            if target<=timeline.time+.001:return False
            timeline.wait_to(target,reason)
            return True
        while timeline.time<deadline-.01:
            t=timeline.time
            mode='EV' if t<ev_end else 'AM' if t>=am_start else 'TJN'
            if mode=='EV' and ev_fixed:mode='TJN'
            if mode=='AM' and (am_fixed or not twilight_routes['AM']):mode='TJN'
            if mode=='EV' and not twilight_routes['EV']:mode='TJN'
            if mode in ('EV','AM'):
                available=twilight_routes[mode]
                if not available:
                    if not wait_until(min(deadline,ev_end if mode=='EV' else end),'no safe twilight route'):break
                    continue
                f=available[route_index%len(available)];route_index+=1;targets=[]
                turn,change=timeline.overhead(f,'r');mid=t+turn+change+15
                m=sky.metrics(f['ra_deg'],f['dec_deg'],[mid])
                seconds=exposure(m['sun_alt'][0],m['alt'][0,0],m['elong'][0,0],m['moon_fraction'][0],m['moon_sep'][0,0],m['moon_alt'][0],c)
            else:
                if not current or t>=campaign_end:
                    horizon=min(deadline,tjn_end,t+config['tjn_campaign_s'])
                    if horizon-t<300:
                        horizon=min(deadline,tjn_end)
                    current=select_fields(catalog,tracks,clouds,geometry,t+224,horizon,
                                          history,evidence,planned,config,stars) if horizon>t+224 else []
                    campaign_end=horizon;dwell=0;route_index=0
                    route_s=sum(movement(current[i-1],f) for i,f in enumerate(current)) if len(current)>1 else 0
                    cycle_frame=config['tjn_exposure_s']+config['readout_s']
                    campaign_dwell=(min(config['tjn_dwell_frames'],max(1,int((config['tjn_max_sampling_gap_s']-route_s-cycle_frame)/((len(current)-1)*cycle_frame))))
                                    if len(current)>1 else config['tjn_dwell_frames'])
                    choices.append(dict(start_local=local_time(t),end_local=local_time(horizon),
                                        fields=[f['field_id'] for f in current],dwell_frames=campaign_dwell,
                                        modeled_largest_gap_s=(len(current)-1)*campaign_dwell*cycle_frame+route_s+cycle_frame))
                if not current:
                    next_time=min(deadline,t+300)
                    if not wait_until(next_time,'no unfinished safe TJN field for this interval'):break
                    campaign_end=0;continue
                f=current[(dwell//campaign_dwell)%len(current)];dwell+=1;targets=f['targets']
                seconds=config['tjn_exposure_s']
            turn,change=timeline.overhead(f,'r');limit=min(deadline,ev_end if mode=='EV' else tjn_end if mode=='TJN' else end)
            if next_event is not None:
                # Stop filling when the next event can start at its requested shutter time.
                nturn=movement(f,next_event['field']);nchange=config['filter_change_s'] if next_event['band']!='r' else 0
                limit=min(limit,next_event['start']-nturn-nchange)
            max_seconds=int(math.floor(limit-t-turn-change-config['readout_s']))
            if max_seconds<5:
                if next_event is not None:
                    nturn,nchange=timeline.overhead(next_event['field'],next_event['band'])
                    actual_limit=min(limit,next_event['start']-nturn-nchange)
                else:actual_limit=limit
                if mode=='TJN':timeline.fill_tail(actual_limit)
                if limit < deadline-1 and limit in (ev_end,am_start):
                    if not wait_until(limit,'sub-frame phase-boundary remainder'):break
                    current=[];campaign_end=0;continue
                break
            seconds=min(seconds,max_seconds)
            check_field={**f,'field_id':f.get('geometry_id',f['field_id'].rsplit('_',1)[-1])}
            if not geometry.check(check_field,t+turn+change,t+turn+change+seconds+config['readout_s'],
                    25 if mode in ('EV','AM') else 30,-10 if mode in ('EV','AM') else -18,mode in ('EV','AM')):
                current=[];campaign_end=0
                # Try another route position rather than exporting an unsafe frame.
                if mode=='TJN':
                    if not wait_until(min(limit,t+60),'no safe TJN exposure at boundary'):break
                elif route_index>len(twilight_routes[mode])*3:
                    if not wait_until(limit,'twilight geometry preempted by commitments'):break
                continue
            r=timeline.append(f,'r',seconds,mode,targets)
            for key in targets:planned[key].append(r['mid_unix'])
        return

    for event in events:
        if event['end']<=timeline.time:continue
        # Use the real preceding pointing/filter, not the worst-case setup allowance.
        fill(event['start'],event)
        turn,change=timeline.overhead(event['field'],event['band'])
        wanted=event['start']-turn-change
        if timeline.time<wanted:
            timeline.wait_to(wanted,'fixed start time requires timed handover')
        start=timeline.time+turn+change
        if start>event['start']+1.1:
            raise ValueError(dict(error='Commitment timing exceeded',field=event['field']['field_id'],
                wanted=local_time(event['start']),actual=local_time(start),delay_s=start-event['start']))
        if event['mode']=='AM':
            for r in event['sequence']:
                timeline.append({**r,'field_id':'MP_'+r['field_id'][2:]},'r',r['exposure_s'],'AM')
        else:
            for band in ([event['band']] if event['mode']=='FOLLOW' else ['g','r']):
                count=event['frames'] if event['mode']=='FOLLOW' else event['remaining'][band]
                for _ in range(count):timeline.append(event['field'],band,event['seconds'],event['mode'],
                    candidate=(event.get('request') or {}).get('candidate_id'),request_sha=(event.get('request') or {}).get('sha256'))
        current=[];campaign_end=0
    fill(end)
    verification=final_validation(timeline.rows,sky,catalog,config,tracks)
    if not verification['passed']:raise ValueError(verification)
    for status in statuses:
        rr=[r for r in timeline.rows if r['candidate_id']==status['candidate_id'] and r.get('request_sha256')==status.get('request_sha256')]
        status['scheduled']=bool(rr)
        if rr:
            status.update(actual_start=local_time(rr[0]['start_unix']),
                          offset_s=rr[0]['start_unix']-datetime.fromisoformat(status['requested_start']).timestamp(),
                          scheduled_frames=len(rr))
        elif 'reason' not in status:status['reason']='request overlaps preserved execution prefix'
    return dict(date=day,night=day.replace('-',''),version=config['version'],rows=timeline.rows,
                start_local=local_time(timeline.rows[0]['overhead_start_unix'] if timeline.rows else begin),
                window_start_local=local_time(begin),end_local=local_time(timeline.time),
                grid=catalog.manifest,config=config,followups=statuses,yso=yso_status,
                tjn_selections=choices,gaps=timeline.gaps,verification=verification,
                twilight_preemptions=any(e['mode']=='FOLLOW' and (e['setup']<ev_end or e['end']>am_start) for e in events),
                execution_update_required=bool(now is not None and now>begin and previous_plan),
                protected_prefix_frames=len(prefix),
                twilight_blocks=[dict(side=b['side'],frames=len(b['rows']),shorter_fallback=b.get('shorter_fallback',False),
                                      skipped=not bool(b['rows']),reason=b['reason']) for b in twilight['blocks']],
                requires_timed_segments=any(g['end_unix']-g['start_unix']>1 for g in timeline.gaps),
                execution_note='Static eight-column TXT has no WAIT or live cutoff. Segment start times must be honoured; revision does not restart telescope execution.')
