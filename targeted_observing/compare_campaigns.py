"""Compare frozen baseline and two alternatives from the cached audited target pool."""
import argparse
import csv
from datetime import date,timedelta
import json
from pathlib import Path
import numpy as np
from .catalog import BASE,horizons_track
from .geometry import coord,offset,in_field,Sky,setup_iers
from .scheduler import field,trojan
from .science import load_periods,assess
from .retention import report
from .bright_stars import assess as bright_assess,load as bright_load

def compare(config,start,cache,output):
    from .config import validate_config
    validate_config(config)
    cache=Path(cache);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    stop=(date.fromisoformat(start)+timedelta(days=config['trojan']['freeze_days'])).isoformat()
    tracks={}
    for row in csv.DictReader((BASE/'data/target_audit_20261103.csv').open()):
        horizons_track(row['pdes'],start,stop,config,cache)
    for p in sorted(cache.glob(f'horizons_*_{start}_{stop}_1h.json')):
        t=json.loads(p.read_text());tracks[t['pdes']]=t
    if not tracks:raise ValueError('No cached audited ephemerides for this campaign; fetch a new target pool first')
    evidence=load_periods(cache);roles={k:assess(evidence.get(k)) for k in tracks}
    bright={k:t for k,t in tracks.items() if max(t['V'])<=config['trojan']['core_V_limit']}
    stars=[s for s in bright_load(cache) if s['band']=='V' and s['mag']<2]
    sc=coord([s['ra_deg'] for s in stars],[s['dec_deg'] for s in stars])
    candidates={}
    for k,t in bright.items():
        center=coord(float(np.rad2deg(np.unwrap(np.deg2rad(t['RA']))).mean())%360,float(np.mean(t['DEC'])))
        for dx in [-.7,0,.7]:
            for dy in [-.7,0,.7]:
                p=offset(center,dx,dy);f=field(p.ra.deg,p.dec.deg,'TJ',start+'-review')
                # Loose hourly preselection only; final denominator is scheduled visits.
                members={j:float(np.mean(in_field(v['RA'],v['DEC'],f['ra_deg'],f['dec_deg'],3,.15))) for j,v in bright.items()}
                f['core_targets']=[j for j,v in members.items() if v>=.6]
                penalty=.25 if np.any(sc.separation(p).deg<4.122) else 1.
                f['selection_contributions']={j:penalty*roles[j]['weight']*members[j]*min(1.,10**(-.2*(max(bright[j]['V'])-17))) for j in f['core_targets']}
                value=sum(f['selection_contributions'].values())
                key=tuple(sorted(f['core_targets']))
                if key and (key not in candidates or value>candidates[key][0]):candidates[key]=(value,f)
    ranked=sorted(candidates.values(),key=lambda x:-x[0])[:35];pairs=[]
    for i,(v,f) in enumerate(ranked):
        for w,g in ranked[i+1:]:
            sep=float(coord(f['ra_deg'],f['dec_deg']).separation(coord(g['ra_deg'],g['dec_deg'])).deg)
            if sep<2.5:continue
            ids=set(f['core_targets']+g['core_targets'])
            # Unique target scoring; small explicit cost for distant pairs, not an exclusion.
            value=sum(max(f['selection_contributions'].get(k,0),g['selection_contributions'].get(k,0)) for k in ids)-.01*sep
            pairs.append((value,[f,g]))
    pairs.sort(key=lambda x:-x[0])
    baseline=json.loads((BASE/'data/campaign_20261103.json').read_text())
    proposals=([('baseline',baseline['fields'])] if baseline['start_date']==start else [])+[(f'alternative_{i+1}',fs) for i,(_,fs) in enumerate(pairs)]
    sky=Sky(config);nights={}
    for i in range(config['trojan']['freeze_days']):
        d=(date.fromisoformat(start)+timedelta(days=i)).isoformat();nights[d]=sky.night(d)
    accepted=[];rejected=[];seen=set()
    for label,fs in proposals:
        identity=tuple(sorted(set(k for f in fs for k in f['core_targets'])))
        if identity in seen:continue
        seen.add(identity);visits={}
        for d,n in nights.items():
            block=trojan(sky,n['dark'],fs,tracks)
            if block:visits[d]=block['rows']
        retain=report(fs,tracks,visits,config)
        if any(r['geometric_retention']<config['trojan']['min_core_retention_fraction'] for r in retain.values()) or not visits:
            rejected.append({'label':label,'reason':'planned_visit_retention_or_visibility','retention':retain});continue
        camp={'schema_version':2,'start_date':start,'stop_date':stop,'fields':fs,
            'tracks':{k:tracks[k] for k in identity},'planned_visits':visits,
            'period_evidence':{k:evidence.get(k) for k in identity},'science_review':{k:roles[k] for k in identity},
            'retention':retain,'selection_scope':'cached audited ephemeris pool; heuristic alternatives, not global optimum',
            'candidate_pool_count':len(tracks),'photometric_usability':'unknown until images'}
        risks=bright_assess(camp,visits,cache);camp['bright_star_risks']=risks
        name='baseline' if label=='baseline' else f'alternative_{len(accepted)}'
        dest=output/name;dest.mkdir(exist_ok=True);(dest/'campaign.json').write_text(json.dumps(camp,indent=2)+'\n')
        accepted.append({'label':name,'fields':fs,'core_objects':list(identity),'roles':{k:roles[k]['role'] for k in identity},
            'scheduled_nights':len(visits),'points':sum(len(r) for r in visits.values()),'retention':retain,
            'bright_risk_rows':len(risks['risks']),'campaign_file':str(dest/'campaign.json')})
        print(name,identity,len(visits),'nights',flush=True)
        if len(accepted)==3:break
    summary={'candidates':accepted,'rejections':rejected,'scope':'fixed 14-day cached candidate pool; no weather or measured photometry'}
    (output/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    if len(accepted)<3:raise ValueError('Fewer than three feasible alternatives; see comparison.json')
    return summary

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',default=str(BASE/'defaults.json'));p.add_argument('--start',default='2026-11-03');p.add_argument('--cache',default=str(BASE/'cache'));p.add_argument('--output',default=str(BASE/'products/pro_review/campaigns'));a=p.parse_args()
    setup_iers(Path(a.cache)/'iers.all');compare(json.loads(Path(a.config).read_text()),a.start,a.cache,a.output)
if __name__=='__main__':main()
