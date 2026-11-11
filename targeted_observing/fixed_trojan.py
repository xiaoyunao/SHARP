"""Opt-in fixed-registry TJN smoke planner; leaves every existing entry point alone."""
import json
from collections import defaultdict
from datetime import date,timedelta
import numpy as np
import pandas as pd
from .catalog import BASE,interpolate
from .fixed_fields import FixedFields,edge_distance
from .fixed_field_build import OUTPUT,write_json
from .annual_trojan import load_tracks,choose_block
from .bright_stars import load as load_stars
from .geometry import Sky,setup_iers,coord,unix
from .validate import validate_plan


def plan_night(catalog,tracks,day,high_dec=False):
    config=json.loads((BASE/'defaults.json').read_text());config['trojan']['exposures_per_dwell']=2
    setup_iers(BASE/'cache/iers_20261104.all');sky=Sky(config);night=sky.night(day)
    fields=catalog.at(tracks,sum(night['dark'])/2)
    if high_dec:fields=[f for f in fields if f['dec_deg']>=40]
    stars=load_stars(epoch=2027.3);sc=coord([s['ra_deg'] for s in stars],[s['dec_deg'] for s in stars]);sm=np.array([s['mag'] for s in stars])
    for f in fields:
        sep=coord(f['ra_deg'],f['dec_deg']).separation(sc).deg
        f['very_bright_star_risk']=bool(np.any((sm<2)&(sep<4.122)));f['bright_stars_within_cap']=int(np.sum(sep<2.122))
    block,reason=choose_block(sky,night,fields,tracks,4,defaultdict(set),{},0)
    plan={'date':day,'config':config,'night':night,'blocks':[block] if block else [],
          'status':'LOCAL_FIXED_FIELD_SMOKE_NOT_FOR_TELESCOPE','skips':[reason] if reason else [],
          'selection':'existing two-hour selector; current-night memberships from fixed registry, no moving centres'}
    retained=[]
    if block:
        for f in block['fields']:
            rows=[r for r in block['rows'] if r['field_id']==f['field_id']]
            times=[unix(r['mid_utc']) for r in rows]
            keep=[]
            for k in f['core_targets']:
                ra,dec=interpolate(tracks[k],times);edges=edge_distance(ra,dec,f['ra_deg'],f['dec_deg'],catalog.grid)
                inside=edges>=catalog.grid.margin+catalog.grid.guard
                retained.append({'field_id':f['field_id'],'object_id':k,'fraction':float(inside.mean()),'minimum_edge_deg':float(edges.min())})
                if inside.mean()>=.8:keep.append(k)
                else:f['secondary_targets'].append(k)
                for r,good,e in zip(rows,inside,edges):r['target_states'][k]={'in_field':bool(good),'edge_distance_deg':float(e),'photometric_usable':None}
            f['core_targets']=keep
        if any(not f['core_targets'] for f in block['fields']):raise ValueError('Fixed field lacks retained targets')
    campaign={'start_date':day,'stop_date':(date.fromisoformat(day)+timedelta(days=1)).isoformat(),
              'fields':fields,'tracks':tracks}
    audit=validate_plan(plan,campaign)
    audit['fixed_quality_retention']=retained
    audit['registry_identity_pass']=all(all(f[k]==catalog.by_id[f['field_id']][k] for k in ['ra_deg','dec_deg','pa_deg','grid_cell_id']) for b in plan['blocks'] for f in b['fields'])
    if not audit['passed'] or not audit['registry_identity_pass']:raise ValueError(audit)
    plan['validation']=audit
    return plan


def main():
    catalog=FixedFields.read(OUTPUT/'registry.json');tracks=load_tracks();out=OUTPUT/'representative_nights';out.mkdir(exist_ok=True)
    jobs=[('2026-12-15','L5_winter',False),('2027-04-15','L4_spring',False),
          ('2027-09-15','season_transition',False),('2027-01-15','high_declination',True)]
    summary=[]
    for day,label,high in jobs:
        plan=plan_night(catalog,tracks,day,high);write_json(out/f'{label}_{day}.json',plan)
        rows=[r for b in plan['blocks'] for r in b['rows']];fs=[f for b in plan['blocks'] for f in b['fields']]
        summary.append({'date':day,'case':label,'fields':len(fs),'field_ids':';'.join(f['field_id'] for f in fs),
            'exposures':len(rows),'exposure_s':sum(r['exposure_s'] for r in rows),
            'slew_s':sum(r['slew_s'] for r in rows),'readout_s':sum(r['readout_s'] for r in rows),
            'unique_targets':len({k for f in fs for k in f['core_targets']}),
            'min_dec_deg':min([f['dec_deg'] for f in fs],default=np.nan),
            'max_dec_deg':max([f['dec_deg'] for f in fs],default=np.nan),
            'passed':plan['validation']['passed'],'reason':';'.join(plan['skips'])})
        print(summary[-1],flush=True)
    pd.DataFrame(summary).to_csv(OUTPUT/'representative_nights.csv',index=False)

if __name__=='__main__':main()
