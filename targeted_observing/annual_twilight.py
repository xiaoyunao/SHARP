"""Annual nominal twilight plans; explicit zero windows and three-night anchors.

Nominal plans assume successful observations. Execution simulations must apply
feedback to AnchorState rather than treating a scheduled visit as observed.
"""
import argparse
import copy
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import date, timedelta
from .catalog import BASE, write_csv
from .annual_sources import dates
from .geometry import Sky, setup_iers, unix
from .scheduler import twilight
from .validate import validate_plan


class AnchorState:
    def __init__(self):self.blocks={};self.effective_nights={}

    def previous(self,side):
        return [block for key,block in self.blocks.items() if key.startswith(side)
                and len(self.effective_nights.get(key,set()))<3]

    def feedback(self,blocks,night,successful_field_ids):
        for block in blocks:
            key=block['block_id'];old=self.blocks.get(key)
            current={f['field_id'] for f in block['fields']}
            old_ids={f['field_id'] for f in old['fields']} if old else set()
            if len(current&old_ids)<3 or len(self.effective_nights.get(key,set()))>=3:
                self.effective_nights[key]=set()
            self.blocks[key]=block
            if current.issubset(successful_field_ids):self.effective_nights.setdefault(key,set()).add(night)


def run_chunk(job):
    nights,cycles,spec=job
    config=json.loads((BASE/'defaults.json').read_text());config['twilight']['visits']=cycles
    config['trojan']['block_s']=0
    provenance=setup_iers(BASE/'cache'/spec['iers_file'])
    state=AnchorState();sky=Sky(config);ledger=[]
    root=BASE/'products/annual'/f'twilight_{cycles}';root.mkdir(parents=True,exist_ok=True)
    for day in nights:
        destination=root/f'{day}.json'
        if destination.exists():
            plan=json.loads(destination.read_text())
            if plan.get('config')!=config:raise ValueError('Stale configuration cache: '+str(destination))
        else:
            night=sky.night(day);blocks=[];zero=[]
            for side,key in [('EV','evening'),('AM','morning')]:
                chosen=twilight(sky,night[key],side,state.previous(side))
                if not chosen:zero.append(side+': no complete group satisfying geometry and chronology')
                blocks.extend(chosen)
            plan={'date':day,'config':config,'blocks':blocks,'night':night,'skips':zero,
                  'status':'LOCAL_ANNUAL_NOMINAL_NOT_FOR_TELESCOPE',
                  'anchor_feedback_assumption':'all planned fields deliver a qualifying short arc; not real observations',
                  'iers':provenance}
            plan['validation']=validate_plan(plan)
            if not plan['validation']['passed']:raise RuntimeError(json.dumps(plan['validation']))
            destination.write_text(json.dumps(plan)+'\n')
        state.feedback(plan['blocks'],day,{f['field_id'] for b in plan['blocks'] for f in b['fields']})
        for side,key in [('EV','evening'),('AM','morning')]:
            selected=[b for b in plan['blocks'] if b['block_id'].startswith(side)]
            rows=[r for b in selected for r in b['rows']]
            exposure=sum(r['exposure_s'] for r in rows);readout=sum(r['readout_s'] for r in rows)
            slew=sum(r['slew_s'] for r in rows);filters=sum(r['filter_change_s'] for r in rows)
            ledger.append({'date':day,'side':side,'strategy':cycles,'window_s':3600,'exposures':len(rows),
                           'exposure_s':exposure,'readout_s':readout,'slew_s':slew,'filter_s':filters,
                           'unused_s':3600-exposure-readout-slew-filters,
                           'zero_reason':'; '.join(plan['skips']) if not rows else '',
                           'validated':plan['validation']['passed']})
    return ledger


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--days',type=int);args=parser.parse_args()
    spec=json.loads((BASE/'annual_defaults.json').read_text());nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    if args.days:nights=nights[:args.days]
    # Each 35-day campaign begins a fresh anchor episode. No future weather is used.
    chunks=[nights[i:i+35] for i in range(0,len(nights),35)]
    jobs=[(chunk,n,spec) for n in (5,10) for chunk in chunks];ledger=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_chunk,job) for job in jobs]
        for i,future in enumerate(as_completed(futures),1):
            ledger.extend(future.result());print(f'Twilight campaign {i}/{len(jobs)} complete',flush=True)
    ledger.sort(key=lambda r:(r['strategy'],r['date'],r['side']))
    root=BASE/'products/annual';write_csv(root/'twilight_time_ledger.csv',ledger)
    (root/'twilight_completeness.json').write_text(json.dumps({'night_count':len(nights),'windows_per_strategy':len(nights)*2,
        'strategies':[5,10],'ledger_rows':len(ledger),'zero_windows':sum(r['exposures']==0 for r in ledger),
        'all_validated':all(r['validated'] for r in ledger),'scope':'nominal clear-weather geometry, not weather-realized yield'},indent=2)+'\n')


if __name__=='__main__':main()
