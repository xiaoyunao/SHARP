"""Feasible historical field menus and weather-feedback twilight anchor queues."""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import date,timedelta
import numpy as np
from .catalog import BASE,write_csv
from .geometry import Sky,setup_iers,unix,coord
from .scheduler import visits,check_rows,block_record
from .annual_sources import dates
from .annual_weather import ensemble
from .annual_timing import prepare,mask as timing_mask,row_times


def identity(block):return tuple(f['field_id'] for f in block['fields'])


def build_day(job):
    day,cycles,spec=job
    root=BASE/'products/annual'/f'twilight_{cycles}';dest=root/'menus';dest.mkdir(exist_ok=True)
    setup_iers(BASE/'cache'/spec['iers_file']);current=json.loads((root/f'{day}.json').read_text());c=current['config'];sky=Sky(c)
    prior=[];d=date.fromisoformat(day)
    for ago in range(15):
        path=root/f'{(d-timedelta(days=ago)).isoformat()}.json'
        if path.exists():prior.extend((ago,b) for b in json.loads(path.read_text())['blocks'])
    result={}
    for side,key in [('EV','evening'),('AM','morning')]:
        a,b=current['night'][key];names=['A'] if cycles==10 else (['A','B'] if side=='EV' else ['B','A'])
        for j,name in enumerate(names):
            blockid=side+'_'+name;lo=a+j*(b-a)/len(names);hi=lo+(b-a)/len(names);seen=set();options=[]
            for ago,old in prior:
                if old['block_id']!=blockid or identity(old) in seen:continue
                seen.add(identity(old));rows=visits(old['fields'],lo,hi,30,cycles,c,blockid,'twilight')
                if len(rows)!=4*cycles:continue
                quality=check_rows(sky,old['fields'],rows,'twilight',True)
                if not quality:continue
                extended=any(q['elongation_max_deg']>75 for q in quality.values())
                block=block_record(blockid,old['fields'],rows,lo,hi,quality,extended);block['nominal_age_days']=ago
                prepare(block,sky);options.append(block)
            result[blockid]=options
    out={'date':day,'config':c,'night':current['night'],'options':result,'lookback_days':14,
         'scope':'only geometrically feasible previously generated groups; no future fields or weather'}
    (dest/f'{day}.json').write_text(json.dumps(out)+'\n');return day


def compatible(a,b,side):
    if a is None or b is None:return True
    f,g=a['fields'],b['fields'];p=coord([x['ra_deg'] for x in f],[x['dec_deg'] for x in f]);q=coord([x['ra_deg'] for x in g],[x['dec_deg'] for x in g])
    if p[:,None].separation(q[None,:]).deg.min()<2.6:return False
    x=np.mean([v['elongation_min_deg'] for v in a['quality'].values()]);y=np.mean([v['elongation_min_deg'] for v in b['quality'].values()])
    return x<=y if side=='EV' else x>=y


def run(cycles,spec):
    root=BASE/'products/annual'/f'twilight_{cycles}';nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    menus=[json.loads((root/'menus'/f'{d}.json').read_text()) for d in nights]
    setup_iers(BASE/'cache'/spec['iers_file'])
    for menu in menus:
        sky=Sky(menu['config'])
        for options in menu['options'].values():
            for block in options:prepare(block,sky)
    weather=ensemble(nights,spec['weather_realizations'],spec['seed']);records=[];selections=[]
    # Compatibility is invariant across weather realizations; evaluate once.
    pair_cache={}
    for ni,menu in enumerate(menus):
        for side in ('EV','AM'):
            names=[side+'_A'] if cycles==10 else [side+'_'+x for x in (['A','B'] if side=='EV' else ['B','A'])]
            opts=[menu['options'][n] for n in names]
            pairs=[(i,-1) for i in range(len(opts[0]))] if len(opts)==1 else [
                (i,j) for i in range(-1,len(opts[0])) for j in range(-1,len(opts[1]))
                if (i>=0 or j>=0) and compatible(opts[0][i] if i>=0 else None,opts[1][j] if j>=0 else None,side)]
            pair_cache[ni,side]=(names,opts,pairs)
    for realization in range(spec['weather_realizations']):
        effective={};previous={};field_nights={}
        def unfinished(ident):return any(len(field_nights.get(fid,set()))<3 for fid in ident)
        for ni,menu in enumerate(menus):
            for side,slot in [('EV',0),('AM',2)]:
                names,opts,pairs=pair_cache[ni,side];w=lambda key:weather[key][realization,ni,slot]
                record={'realization':realization,'date':menu['date'],'side':side,'scheduled_exposures':0,'good_exposures':0,
                        'qualified_fields':0,'retained_anchors':0,'reason':'','anchor_replacement_reason':''}
                if not w('observable') or w('reserved'):
                    record['reason']='weather_closed' if not w('observable') else 'external_allocation_cancelled';records.append(record)
                    selections.append([realization,ni,slot,-1,-1]);continue
                def score(pair):
                    count=retained=novel=0;age=0
                    for gi,option in enumerate(pair[:len(opts)]):
                        if option<0:continue
                        b=opts[gi][option];ident=identity(b);count+=1;age+=b['nominal_age_days']
                        if ident==previous.get(names[gi]) and unfinished(ident):retained+=4
                        if unfinished(ident):novel+=1
                    return (novel,retained,count,-age)
                if not pairs:
                    record['reason']='no_feasible_current_or_historical_group';records.append(record);selections.append([realization,ni,slot,-1,-1]);continue
                pair=max(pairs,key=score);rng=np.random.default_rng(int(w('point_seed')))
                selections.append([realization,ni,slot,*pair])
                for gi,option in enumerate(pair[:len(opts)]):
                    if option<0:continue
                    b=opts[gi][option];ident=identity(b);record['scheduled_exposures']+=len(b['rows'])
                    if ident==previous.get(names[gi]) and unfinished(ident):record['retained_anchors']+=4
                    old=previous.get(names[gi])
                    if old and ident!=old:
                        reason='three effective nights completed' if not unfinished(old) else 'old group infeasible/unavailable or incompatible with other slot'
                        record['anchor_replacement_reason']+=names[gi]+': '+reason+'; '
                    previous[names[gi]]=ident
                    ts=row_times(b);shift=w('delay_s')+np.arange(len(ts))*(w('extra_readout_s')+w('extra_slew_s'))
                    actual=ts+shift;lo=unix(b['not_before_utc']);hi=unix(b['hard_deadline_utc']);fraction=(actual-lo)/(hi-lo)
                    cloud=(fraction>=w('cloud_start_fraction'))&(fraction<w('cloud_start_fraction')+w('cloud_duration_fraction'))
                    safe=timing_mask(b,[r['field_id'] for r in b['rows']],actual)
                    good=(actual+30<hi)&~cloud&safe&(rng.random(len(ts))<.95)
                    record['good_exposures']+=int(good.sum());good_fields=0
                    for f in b['fields']:
                        mask=np.array([r['field_id']==f['field_id'] for r in b['rows']])&good;t=actual[mask]
                        if len(t)>=4 and np.ptp(t)>=720:
                            good_fields+=1;field_nights.setdefault(f['field_id'],set()).add(menu['date'])
                    record['qualified_fields']+=good_fields
                    if good_fields==4:effective.setdefault(ident,set()).add(menu['date'])
                records.append(record)
        if (realization+1)%20==0:print(f'Twilight {cycles}: {realization+1}/200 causal queues',flush=True)
    write_csv(root/'causal_execution_ledger.csv',records)
    np.savez_compressed(root/'causal_choices.npz',choices=np.array(selections,dtype=np.int32))
    totals=[sum(r['qualified_fields'] for r in records if r['realization']==i) for i in range(spec['weather_realizations'])]
    result={'strategy_visits':cycles,'realizations':spec['weather_realizations'],'qualified_field_windows_P10_P50_P90':np.percentile(totals,[10,50,90]).tolist(),
            'window_records':len(records),'policy':'three successful group nights, then prefer a new group; retain only while geometrically legal',
            'limitation':'14-day historical menu; reset if geometry fails; same-field reacquisition is not moving-object recovery'}
    (root/'causal_queue_summary.json').write_text(json.dumps(result,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser();p.add_argument('--visits',type=int,choices=(5,10),required=True);p.add_argument('--workers',type=int,default=3)
    p.add_argument('--menus-only',action='store_true');p.add_argument('--simulate-only',action='store_true');args=p.parse_args()
    spec=json.loads((BASE/'annual_defaults.json').read_text());nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    if not args.simulate_only:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for i,future in enumerate(as_completed([pool.submit(build_day,(d,args.visits,spec)) for d in nights]),1):
                future.result()
                if i%35==0:print(f'Twilight {args.visits}: {i}/365 menus',flush=True)
    if not args.menus_only:run(args.visits,spec)


if __name__=='__main__':main()
