"""Select a frozen 14-day campaign using JPL observer ephemerides and period evidence."""
import argparse
from datetime import date,timedelta
import hashlib
import json
from pathlib import Path
import numpy as np
from .catalog import BASE,download,approximate,horizons_track,field_members,write_csv
from .geometry import setup_iers,Sky,coord,offset
from .scheduler import field,trojan
from .science import load_periods,evaluate


def design(config,start,cache,output,limit=60):
    from .config import validate_config
    validate_config(config)
    cache=Path(cache);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    d=date.fromisoformat(start);stop=(d+timedelta(days=config['trojan']['freeze_days'])).isoformat()
    mid=(d+timedelta(days=7)).isoformat()
    catalog=json.loads(download('https://ssd-api.jpl.nasa.gov/sbdb_query.api',cache/'trojans.json',
        {'sb-class':'TJN','fields':'spkid,pdes,full_name,epoch,a,e,i,om,w,ma,H','full-prec':'true'}).read_text())
    evidence=load_periods(cache)
    records=approximate(catalog,mid)
    records=[r for r in records if r['V_proxy']<=config['trojan']['extended_V_limit']+.3 and r['dec_deg']>-15]
    sky=Sky(config);night=sky.night(start)
    grid=np.linspace(*night['dark'],19)
    m=sky.metrics([r['ra_deg'] for r in records],[r['dec_deg'] for r in records],grid)
    keep=(m['alt']>35).sum(axis=1)>=4
    records=[r for r,k in zip(records,keep) if k]
    cs=coord([r['ra_deg'] for r in records],[r['dec_deg'] for r in records])
    distances=cs[:,None].separation(cs[None,:]).deg
    weights=np.array([evaluate(evidence.get(r['pdes']))[1] if r['V_proxy']<=18.5 else .2 for r in records])
    scores=(distances<4).dot(weights)
    selected=np.argsort(-scores,kind='stable')[:limit]
    shortlist=[records[i] for i in selected]
    tracks={};failures=[]
    for n,r in enumerate(shortlist):
        try:
            t=horizons_track(r['pdes'],start,stop,config,cache)
            tracks[r['pdes']]=t
            print(f'{n+1}/{len(shortlist)} Horizons {r["pdes"]}',flush=True)
        except Exception as e:failures.append({'pdes':r['pdes'],'error':str(e)})
    # Each target's full hourly path must stay in the inset footprint for all 14 days.
    bright={k:t for k,t in tracks.items() if max(t['V'])<=config['trojan']['core_V_limit']}
    candidates=[]
    for k,t in bright.items():
        ra=float(np.rad2deg(np.unwrap(np.deg2rad(t['RA']))).mean())%360;dec=float(np.mean(t['DEC']))
        for dx in [-1,-.5,0,.5,1]:
            for dy in [-1,-.5,0,.5,1]:
                pt=offset(coord(ra,dec),dx,dy)
                f=field(pt.ra.deg,pt.dec.deg,'TJ',start)
                f['core_targets']=field_members(f,bright,config['fov_deg'],config['edge_margin_deg'])
                if not f['core_targets']:continue
                value=sum(evaluate(evidence.get(v))[1] for v in f['core_targets'])
                candidates.append((value,f))
    candidates.sort(key=lambda x:x[0],reverse=True)
    # Keep best representatives per target set; retain shifted alternatives for close pairs.
    unique=[];seen=set()
    for value,f in candidates:
        key=(tuple(f['core_targets']),round(f['ra_deg'],0),round(f['dec_deg'],0))
        if key not in seen:unique.append((value,f));seen.add(key)
    pairs=[]
    for i,(v,f) in enumerate(unique):
        for w,g in unique[i+1:]:
            if set(f['core_targets'])&set(g['core_targets']):continue
            sep=float(coord(f['ra_deg'],f['dec_deg']).separation(coord(g['ra_deg'],g['dec_deg'])).deg)
            if not 2.5<=sep<=config['adjacent_max_deg']:continue
            ids=f['core_targets']+g['core_targets'];statuses=[evaluate(evidence.get(k))[0] for k in ids]
            bonus=1 if 'control_candidate' in statuses else 0
            pairs.append((v+w+bonus,sep,[f,g]))
    pairs.sort(key=lambda x:(-x[0],x[1]))
    chosen=None
    for value,sep,fs in pairs[:40]:
        if trojan(sky,night['dark'],fs,tracks):chosen=fs;break
    if chosen is None:raise RuntimeError('No eligible adjacent pair with 14-day core retention; inspect shortlist/failures')
    memberids=set(k for f in chosen for k in f['core_targets'])
    selected_tracks={k:t for k,t in tracks.items() if k in memberids}
    rows=[]
    for r in shortlist:
        k=r['pdes'];ev=evidence.get(k);status,_=evaluate(ev)
        rows.append({'pdes':k,'name':r['full_name'].strip(),'cloud_proxy':r['cloud_proxy'],
            'V_min':min(tracks[k]['V']) if k in tracks else None,'V_max':max(tracks[k]['V']) if k in tracks else None,
            'magnitude_band':'V (not measured r)','period_status':status,'selected_core':k in memberids,
            'evidence':json.dumps(ev,ensure_ascii=False),'ephemeris_status':'ok' if k in tracks else 'failed'})
    result={'schema_version':1,'start_date':start,'stop_date':stop,'fields':chosen,'tracks':selected_tracks,
            'selection_scope':{'catalog_count':catalog['count'],'preselected_count':len(records),'horizons_count':len(tracks),
            'method':'two-body density/period ranking then hourly topocentric Horizons; heuristic, not exhaustive optimum',
            'cloud':'heliocentric-position proxy only; membership from JPL TJN class'},
            'period_evidence':{k:evidence.get(k) for k in memberids},'failures':failures,
            'config_sha256':hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest(),
            'retention_margin_deg':config['edge_margin_deg'],'photometry':'V proxy only; r performance unmeasured'}
    (output/'campaign.json').write_text(json.dumps(result,indent=2)+'\n')
    write_csv(output/'target_audit.csv',rows)
    print(json.dumps({'fields':chosen,'failures':failures},indent=2))
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--start',default='2026-11-03');p.add_argument('--config',default=str(BASE/'defaults.json'));p.add_argument('--cache',default=str(BASE/'cache'));p.add_argument('--output',required=True);p.add_argument('--limit',type=int,default=60)
    a=p.parse_args();c=json.loads(Path(a.config).read_text());setup_iers(Path(a.cache)/'iers.all');design(c,a.start,a.cache,a.output,a.limit)

if __name__=='__main__':main()
