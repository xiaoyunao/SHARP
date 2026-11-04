"""Reusable bounded-geometry menus for causal annual queue simulations.

Every allowed two-hour interval is checked on a 300-second grid with a half-grid
rotation guard. A simulated queue can choose fields using its own past outcome;
it does not have to replay a queue that knew future weather.
"""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import date,timedelta
import numpy as np
from .catalog import BASE,write_csv
from .annual_sources import dates
from .annual_trojan import load_tracks,campaign_groups,version_fields
from .geometry import Sky,setup_iers
from .bright_stars import load as load_stars


def visibility_chunk(job):
    nights,spec,limit=job
    tracks=load_tracks(limit);config=json.loads((BASE/'defaults.json').read_text())
    setup_iers(BASE/'cache'/spec['iers_file']);sky=Sky(config);stars=load_stars(epoch=2027.3)
    stop=(date.fromisoformat(nights[-1])+timedelta(days=1)).isoformat()
    groups=campaign_groups(tracks,nights[0],stop)
    root=BASE/'products/annual'/('visibility_smoke' if limit else 'visibility');root.mkdir(parents=True,exist_ok=True)
    result=[]
    for i,day in enumerate(nights):
        if i%spec['field_version_days']==0:
            version_stop=min(date.fromisoformat(day)+timedelta(days=spec['field_version_days']),date.fromisoformat(stop)).isoformat()
            fields=version_fields(groups,tracks,day,version_stop,stars)
        night=sky.night(day);a,b=night['dark'];grid=np.arange(a,b+.01,300.)
        starts=grid[grid+7200<=b+.01]
        valid=np.zeros((len(fields),len(starts)),dtype=bool)
        if fields and len(starts):
            m=sky.metrics([f['ra_deg'] for f in fields],[f['dec_deg'] for f in fields],grid)
            # Entire spherical footprint; 0.645deg between-grid time bound and 0.01deg numeric margin.
            altlim=max(config['center_alt_deg'],config['min_alt_deg']+2.122)+.0043*150+.01
            # Moon separation varies slowly; use the same larger rotation bound conservatively,
            # and the bright-Moon limit at all phases so a threshold crossing cannot slip through.
            moonlim=config['bright_moon_sep_deg']+2.122+.0043*150+.01
            good=(m['alt']>=altlim)&((m['moon_sep']>=moonlim)|(m['moon_alt'][None,:]<-.0043*150-.01))
            for j in range(len(starts)):valid[:,j]=np.all(good[:,j:j+25],axis=1)
            for j,f in enumerate(fields):
                if f['very_bright_star_risk']:valid[j]=False
        np.savez_compressed(root/f'{day}.npz',starts=starts,valid=valid)
        (root/f'{day}.json').write_text(json.dumps({'date':day,'night':night,'fields':fields,
            'campaign_start':nights[0],'campaign_stop':stop,'menu_span_s':7200,
            'grid_step_s':300,'time_guard_deg':.655,
            'moon_policy':'conservative bright-moon limit for all phases',
            'zero_reason':None if valid.any() else 'no full two-hour interval for campaign fields',
            'state':'planning envelope; dispatch still requires individual exposure validation'})+'\n')
        result.append({'date':day,'candidate_fields':len(fields),'eligible_fields':int(np.any(valid,axis=1).sum()),
                       'eligible_field_intervals':int(valid.sum()),'night_dark_s':b-a})
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--days',type=int);p.add_argument('--limit',type=int)
    args=p.parse_args();spec=json.loads((BASE/'annual_defaults.json').read_text());nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    if args.days:nights=nights[:args.days]
    jobs=[(nights[i:i+spec['campaign_days']],spec,args.limit) for i in range(0,len(nights),spec['campaign_days'])]
    results=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i,future in enumerate(as_completed([pool.submit(visibility_chunk,job) for job in jobs]),1):
            results.extend(future.result());print(f'Visibility campaign {i}/{len(jobs)}',flush=True)
    results.sort(key=lambda r:r['date']);suffix='_smoke' if args.limit else ''
    write_csv(BASE/f'products/annual/visibility_ledger{suffix}.csv',results)


if __name__=='__main__':main()
