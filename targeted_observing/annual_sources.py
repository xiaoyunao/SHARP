"""Full-catalog annual screening and resumable, bounded Horizons acquisition."""
import argparse
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, timedelta, datetime, timezone
import numpy as np
from astroquery.jplhorizons import Horizons
from .catalog import BASE, approximate, write_csv


def dates(start, stop):
    a, b = date.fromisoformat(start), date.fromisoformat(stop)
    return [(a+timedelta(days=i)).isoformat() for i in range((b-a).days)]


def screen(spec, output):
    path = BASE/'cache/trojans.json'
    catalog = json.loads(path.read_text())
    records = {str(r['pdes']): r for r in (dict(zip(catalog['fields'], x)) for x in catalog['data'])}
    best = {}; accepted = set()
    nights = dates(spec['start_date'], spec['stop_date_exclusive'])
    samples = sorted(set(nights[::spec['preselection_step_days']]+[nights[-1]]))
    for day in samples:
        for r in approximate(catalog, day, spec['candidate_max_H']):
            k = str(r['pdes'])
            if k not in best or r['V_proxy'] < best[k]['V_proxy']: best[k] = r
            if r['V_proxy'] <= spec['candidate_V_limit'] and r['dec_deg'] >= spec['candidate_min_dec_deg']:
                accepted.add(k)
    audit = []
    for k, r in records.items():
        reason = 'accepted_for_Horizons' if k in accepted else ('missing_H' if r['H'] is None else
                 'H_screen' if float(r['H']) >= spec['candidate_max_H'] else 'weekly_brightness_or_declination_screen')
        b = best.get(k, {})
        audit.append({'object_id': k, 'name': r['full_name'], 'H': r['H'],
                      'minimum_sampled_V_proxy': b.get('V_proxy'), 'cloud_proxy': b.get('cloud_proxy'),
                      'status': reason, 'classification': 'JPL TJN; cloud is instantaneous longitude proxy'})
    write_csv(output/'candidate_audit.csv', audit)
    selected = [best[k] for k in sorted(accepted)]
    (output/'candidates.json').write_text(json.dumps(selected, indent=2)+'\n')
    (output/'screening_manifest.json').write_text(json.dumps({'source': 'https://ssd-api.jpl.nasa.gov/sbdb_query.api',
        'catalog_sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'catalog_rows': len(records),
        'accepted': len(selected), 'sample_dates': samples, 'specification': spec,
        'limitation': 'Conservative two-body screen, not population completeness proof. Missing-H rows remain explicitly excluded.'}, indent=2)+'\n')
    return selected


def fetch_one(k, spec, config):
    stop = (date.fromisoformat(spec['stop_date_exclusive'])+timedelta(days=spec['lookahead_days']+1)).isoformat()
    path = BASE/'cache/annual'/f"horizons_{k}_{spec['start_date']}_{stop}_{spec['horizons_step']}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists(): return str(path)
    s = config['site']; error = None
    for attempt in range(3):
        try:
            query = Horizons(id=k, id_type='smallbody', location={'lon': s['lon_deg'], 'lat': s['lat_deg'], 'elevation': s['height_m']/1000},
                             epochs={'start': spec['start_date'], 'stop': stop, 'step': spec['horizons_step']})
            table = query.ephemerides(quantities='1,3,9,19,20,23,24', cache=False)
            data = {key: np.ma.filled(table[key], np.nan).tolist() for key in
                    ['datetime_jd', 'RA', 'DEC', 'RA_rate', 'DEC_rate', 'V', 'r', 'delta', 'elong', 'alpha']}
            data.update(pdes=k, targetname=str(table['targetname'][0]), site=s,
                        source='JPL Horizons topocentric astrometric ICRF; UTC', query=query.uri,
                        retrieved_utc=datetime.now(timezone.utc).isoformat())
            content=json.dumps(data, allow_nan=False)
            path.write_text(content+'\n')
            return str(path)
        except Exception as exc:
            error=str(exc)
            time.sleep(2*(attempt+1))
    raise RuntimeError(f'{k}: {error}')


def main():
    p=argparse.ArgumentParser();p.add_argument('--screen-only',action='store_true');p.add_argument('--workers',type=int,default=2)
    args=p.parse_args();spec=json.loads((BASE/'annual_defaults.json').read_text());config=json.loads((BASE/'defaults.json').read_text())
    output=BASE/'products/annual';output.mkdir(parents=True,exist_ok=True)
    candidates=screen(spec,output);print(f'Full TJN screen: {len(candidates)} candidates',flush=True)
    if args.screen_only:return
    success={};failures={}
    with ThreadPoolExecutor(max_workers=max(1,min(3,args.workers))) as pool:
        jobs={pool.submit(fetch_one,str(r['pdes']),spec,config):str(r['pdes']) for r in candidates}
        for i,job in enumerate(as_completed(jobs),1):
            k=jobs[job]
            try:success[k]=job.result()
            except Exception as exc:failures[k]=str(exc)
            if i%20==0 or i==len(jobs): print(f'Horizons {i}/{len(jobs)}; failed {len(failures)}',flush=True)
    (output/'ephemeris_index.json').write_text(json.dumps({'paths':success,'failures':failures},indent=2)+'\n')
    if failures:raise SystemExit(f'{len(failures)} ephemerides failed; rerun resumes successful cache')


if __name__=='__main__':main()
