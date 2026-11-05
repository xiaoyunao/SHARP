"""Replay weather realization zero for actual-frame footprint unions.

This is geometric coverage after timing/cloud/bad-frame losses. Magnitude and
motion efficiency are separately measured by annual_coverage injections.
"""
import json
import numpy as np
import pandas as pd
from .catalog import BASE,write_csv
from .annual_sources import dates
from .annual_weather import ensemble
from .annual_timing import mask,row_times,prepare
from .annual_coverage import union_area
from .geometry import unix,isot,Sky,setup_iers


def run(visits,root_name=None):
    root=BASE/'products/annual'/(root_name or f'twilight_{visits}')
    spec=json.loads((BASE/'annual_defaults.json').read_text());nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    setup_iers(BASE/'cache'/spec['iers_file'])
    weather=ensemble(nights,200,spec['seed'])
    choices={tuple(r[:3]):tuple(r[3:]) for r in np.load(root/'causal_choices.npz')['choices'] if r[0]==0}
    expected=pd.read_csv(root/'causal_execution_ledger.csv');expected=expected[expected.realization==0].set_index(['date','side'])
    plans=[]
    for ni,day in enumerate(nights):
        menu=json.loads((root/'menus'/f'{day}.json').read_text());blocks=[];sky=Sky(menu['config'])
        for side,slot in [('EV',0),('AM',2)]:
            w=lambda key:weather[key][0,ni,slot]
            rng=np.random.default_rng(int(w('point_seed')));actual_count=qualified=0
            from .annual_twilight_queue import group_names
            names=group_names(menu,side)
            for name,choice in zip(names,choices[0,ni,slot]):
                if choice<0:continue
                b=menu['options'][name][choice];prepare(b,sky);ts=row_times(b)
                actual=ts+w('delay_s')+np.arange(len(ts))*(w('extra_readout_s')+w('extra_slew_s'))
                lo,hi=unix(b['not_before_utc']),unix(b['hard_deadline_utc']);fraction=(actual-lo)/(hi-lo)
                cloud=(fraction>=w('cloud_start_fraction'))&(fraction<w('cloud_start_fraction')+w('cloud_duration_fraction'))
                good=(actual+30<hi)&~cloud&mask(b,[r['field_id'] for r in b['rows']],actual)&(rng.random(len(ts))<.95)
                rows=[{**r,'mid_utc':isot(t)} for r,t,g in zip(b['rows'],actual,good) if g]
                fields=[f for f in b['fields'] if any(r['field_id']==f['field_id'] for r in rows)]
                actual_count+=len(rows)
                for f in fields:
                    t=[unix(r['mid_utc']) for r in rows if r['field_id']==f['field_id']]
                    qualified+=len(t)>=4 and np.ptp(t)>=720
                if rows:blocks.append({'block_id':name,'rows':rows,'fields':fields})
            assert actual_count==expected.loc[(day,side),'good_exposures']
            assert qualified==expected.loc[(day,side),'qualified_fields']
        plans.append({'date':day,'blocks':blocks})
    rows,summary=union_area(plans,1024);summary.update(realization=0,nights=365,
        scope='causal executed geometry after cloud/timing/bad-frame loss; no magnitude-dependent detection applied')
    write_csv(root/'effective_coverage_realization0.csv',rows)
    (root/'effective_coverage_realization0.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(visits,summary,flush=True)


if __name__=='__main__':
    for n in [5,10]:run(n)
