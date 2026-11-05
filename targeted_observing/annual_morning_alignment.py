"""EV five visits / AM ten visits aligned to the -12 degree morning cutoff.

Uses the identical frozen 14-day candidate menus, then rechecks shifted geometry.
It does not claim an optimal newly tiled morning search or real-object linking.
"""
import argparse
import copy
import json
from concurrent.futures import ProcessPoolExecutor,as_completed
from .catalog import BASE,write_csv
from .geometry import Sky,setup_iers,unix,isot
from .scheduler import check_rows,block_record
from .annual_timing import prepare
from .annual_sources import dates
from .annual_twilight_queue import run
from .validate import validate_plan

NAME='twilight_morning10_aligned'
TIME_KEYS=('overhead_start_utc','start_utc','mid_utc','end_utc','release_utc')


def aligned(block,sky,guard_s=30.):
    lo=unix(block['not_before_utc']);hi=unix(block['hard_deadline_utc'])
    rows=copy.deepcopy(block['rows']);shift=hi-guard_s-unix(rows[-1]['release_utc'])
    if shift<0:return None
    for r in rows:
        for k in TIME_KEYS:r[k]=isot(unix(r[k])+shift)
    # The full hour remains reserved. Explicit idle time precedes the first slew;
    # the validator's overhead chronology still begins at the window opening.
    rows[0]['overhead_start_utc']=isot(lo)
    rows[0]['planned_wait_s']=shift
    q=check_rows(sky,block['fields'],rows,'twilight',True)
    if not q:return None
    out=block_record(block['block_id'],block['fields'],rows,lo,hi,q,any(v['elongation_max_deg']>75 for v in q.values()))
    out.update(nominal_age_days=block['nominal_age_days'],end_alignment_shift_s=shift,end_alignment_guard_s=guard_s)
    prepare(out,sky);return out


def build_day(day):
    spec=json.loads((BASE/'annual_defaults.json').read_text());setup_iers(BASE/'cache'/spec['iers_file'])
    root=BASE/'products/annual';five=json.loads((root/'twilight_5/menus'/f'{day}.json').read_text())
    ten=json.loads((root/'twilight_10/menus'/f'{day}.json').read_text());c=copy.deepcopy(five['config'])
    c['twilight']['visits_by_side']={'EV':5,'AM':10};sky=Sky(c)
    options={k:v for k,v in five['options'].items() if k.startswith('EV')}
    options['AM_A']=[b for old in ten['options']['AM_A'] if (b:=aligned(old,sky)) is not None]
    menu={**five,'config':c,'options':options,'scope':'EV unchanged five visits; AM ten visits end-aligned and revalidated'}
    dest=root/NAME;(dest/'menus').mkdir(parents=True,exist_ok=True)
    (dest/'menus'/f'{day}.json').write_text(json.dumps(menu)+'\n')
    # Preserve the exact old nominal EV groups; AM nominal uses the newest legal group.
    original=json.loads((root/'twilight_5'/f'{day}.json').read_text())
    blocks=[b for b in original['blocks'] if b['block_id'].startswith('EV')]
    if options['AM_A']:blocks.append(min(options['AM_A'],key=lambda b:b['nominal_age_days']))
    plan={**original,'config':c,'blocks':blocks,'skips':[] if options['AM_A'] else ['AM: no end-aligned ten-visit group'],
          'status':'LOCAL_ANNUAL_NOMINAL_NOT_FOR_TELESCOPE'}
    audit=validate_plan(plan)
    if not audit['passed']:raise ValueError(audit)
    plan['validation']=audit;(dest/f'{day}.json').write_text(json.dumps(plan)+'\n')
    return {'date':day,'morning_options_before':len(ten['options']['AM_A']),'morning_options_after':len(options['AM_A']),
            'nominal_morning_exposures':40 if options['AM_A'] else 0,'validated':True}


def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=3);p.add_argument('--simulate-only',action='store_true');a=p.parse_args()
    spec=json.loads((BASE/'annual_defaults.json').read_text());spec['twilight_root']=NAME
    if not a.simulate_only:
        rows=[]
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            for i,f in enumerate(as_completed([pool.submit(build_day,d) for d in dates(spec['start_date'],spec['stop_date_exclusive'])]),1):
                rows.append(f.result())
                if i%50==0:print(f'Morning alignment {i}/365',flush=True)
        write_csv(BASE/'products/annual'/NAME/'alignment_ledger.csv',sorted(rows,key=lambda r:r['date']))
    run(10,spec)


if __name__=='__main__':main()
