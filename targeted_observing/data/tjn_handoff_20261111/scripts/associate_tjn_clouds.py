"""Read-only source association for SHARP fixed-field planning."""
from pathlib import Path
import json, hashlib, math
import pandas as pd

import argparse
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--inputs',type=Path,default=Path(__file__).resolve().parents[1]/'inputs')
p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'regenerated_inputs')
args=p.parse_args()
SRC=args.inputs
OUT=args.output
OUT.mkdir(parents=True,exist_ok=True)
parent_file=SRC/'full_parent_selection_function.csv.gz'
fields_file=SRC/'B_fields_original.csv'
special_file=SRC/'special_science_candidates.json'
parent=pd.read_csv(parent_file,dtype={'object_id':str})
base_cols=['object_id','name','H','minimum_sampled_V_proxy','cloud_proxy','status','classification']
for key,group in parent.groupby('object_id'):
    assert len(group[base_cols].drop_duplicates())==1, f'Inconsistent source screening: {key}'
parent_unique=parent.drop_duplicates('object_id')[base_cols].copy()
lookup=parent_unique.set_index('object_id').to_dict(orient='index')
fields=pd.read_csv(fields_file,dtype={'field_id':str,'anchor_object':str})
fields=fields.loc[fields['mode']=='trojan'].copy()

def object_id(value):
    if pd.isna(value) or str(value).strip()=='': return None
    value=str(value).strip()
    if value.endswith('.0') and value[:-2].isdigit(): value=value[:-2]
    return value

def members(row):
    ids=[]
    for col in ['core_targets_json','secondary_targets_json']:
        if pd.notna(row[col]): ids.extend(object_id(x) for x in json.loads(row[col]))
    ids.append(object_id(row['anchor_object']))
    return sorted(set(x for x in ids if x is not None))

association=[]
for _,row in fields.iterrows():
    ids=members(row)
    labels=sorted(set(lookup[k]['cloud_proxy'] for k in ids if k in lookup and pd.notna(lookup[k]['cloud_proxy'])))
    missing=[k for k in ids if k not in lookup or pd.isna(lookup[k]['cloud_proxy'])]
    cloud=labels[0] if len(labels)==1 and not missing else ('mixed' if len(labels)>1 else 'unknown')
    if len(labels)==1 and missing: cloud=labels[0]+'_partially_unknown'
    association.append({'cloud_proxy':cloud,'associated_objects_json':json.dumps(ids),
        'associated_object_count':len(ids),'cloud_label_set_json':json.dumps(labels),
        'missing_cloud_objects_json':json.dumps(missing),
        'cloud_assignment_basis':'union(core_targets,secondary_targets,anchor_object) joined to frozen parent cloud_proxy'})
fields=pd.concat([fields.reset_index(drop=True),pd.DataFrame(association)],axis=1)
fields.to_csv(OUT/'tjn_B_field_nights_with_cloud.csv',index=False)
unique_rows=[]
for fid,group in fields.groupby('field_id',sort=True):
    for col in ['ra_deg','dec_deg','pa_deg','field_version','cloud_proxy']:
        assert group[col].nunique(dropna=False)==1,(fid,col)
    row=group.iloc[0].to_dict()
    ids=sorted(set(k for s in group['associated_objects_json'] for k in json.loads(s)))
    row.update(first_night=group['date'].min(),last_night=group['date'].max(),
               field_nights=len(group),associated_objects_json=json.dumps(ids),associated_object_count=len(ids))
    unique_rows.append(row)
unique=pd.DataFrame(unique_rows)
unique.to_csv(OUT/'tjn_B_unique_fields_with_cloud.csv',index=False)
candidates=parent_unique.loc[parent_unique['status']=='accepted_for_Horizons'].copy()
candidates.to_csv(OUT/'tjn_annual_candidates_1063.csv',index=False)
parent_unique.to_csv(OUT/'tjn_parent_screening_16445.csv',index=False)

stats=[]
for cloud,g in fields.groupby('cloud_proxy'):
    ids=sorted(set(k for s in g['associated_objects_json'] for k in json.loads(s)))
    cores=sorted(set(str(k) for s in g['core_targets_json'].dropna() for k in json.loads(s)))
    stats.append({'cloud_proxy':cloud,'field_nights':len(g),'unique_fields':int(g['field_id'].nunique()),
                  'unique_associated_objects':len(ids),'unique_core_objects':len(cores),
                  'ra_min_deg':float(g['ra_deg'].min()),'ra_max_deg':float(g['ra_deg'].max()),
                  'dec_min_deg':float(g['dec_deg'].min()),'dec_max_deg':float(g['dec_deg'].max()),
                  'object_ids':ids})
bands=[]
for cloud,g in unique.groupby('cloud_proxy'):
    for band,h in g.groupby(g['dec_deg'].map(lambda d:3*math.floor(d/3))):
        ras=sorted(h['ra_deg'].unique());segments=[];start=last=ras[0]
        for r in ras[1:]:
            if r-last>10: segments.append([start,last]);start=r
            last=r
        segments.append([start,last])
        bands.append({'cloud_proxy':cloud,'dec_band_min_deg':band,'dec_band_max_deg':band+3,
                      'unique_fields':len(h),'field_nights':int(h['field_nights'].sum()),
                      'ra_min_deg':min(ras),'ra_max_deg':max(ras),
                      'ra_groups_gap_gt10deg_json':json.dumps(segments),
                      'group_note':'Descriptive split at gaps >10 deg in selected centres; not a science coverage boundary'})
pd.DataFrame(bands).to_csv(OUT/'tjn_B_declination_bands_3deg.csv',index=False)
special=json.loads(special_file.read_text())
special_list=special['targets'] if 'targets' in special else special.get('candidates',[])
if not special_list:
    special_list=next(v for v in special.values() if isinstance(v,list) and v and isinstance(v[0],dict) and 'mpc_number' in v[0])
special_rows=[]
for entry in special_list:
    key=str(entry['mpc_number']);p=lookup.get(key,{})
    special_rows.append({'object_id':key,'name':entry.get('name'),
                         'cloud_proxy':p.get('cloud_proxy'),'screen_status':p.get('status'),
                         'in_annual_1063':p.get('status')=='accepted_for_Horizons',
                         'source_parent_name':p.get('name'),'tier':entry.get('tier'),
                         'tags_json':json.dumps(entry.get('tags',[]))})
pd.DataFrame(special_rows).to_csv(OUT/'tjn_special_16_cloud_labels.csv',index=False)
summary={'source_commit':'7747cc1650fdeb255e677e496f50a376a99cb57b',
         'inputs_sha256':{str(p.relative_to(SRC)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [parent_file,fields_file,special_file]},
         'parent_rows_unique':len(parent_unique),'candidate_rows':len(candidates),
         'candidate_cloud_counts':candidates['cloud_proxy'].value_counts(dropna=False).to_dict(),
         'field_nights':len(fields),'unique_fields':len(unique),'cloud_stats':stats,
         'missing_cloud_field_nights':int((fields['missing_cloud_objects_json']!='[]').sum()),
         'mixed_cloud_field_nights':int((fields['cloud_proxy']=='mixed').sum()),
         'tonight_targets':{k:lookup.get(k) for k in ['17365','4722']},
         'special_targets':special_rows,
         'limitations':['Cloud labels are frozen instantaneous heliocentric longitude proxies, not new dynamical classifications.',
                        'Selected B fields are not the full candidate ephemeris envelope.',
                        'Candidates table contains IDs and original screening metadata, no trajectories.',
                        'Three-degree declination bins and ten-degree RA gaps are descriptive reporting choices.']}
(OUT/'tjn_cloud_association_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='cloud_stats'},ensure_ascii=False,indent=2))
print(json.dumps([{k:v for k,v in s.items() if k!='object_ids'} for s in stats],ensure_ascii=False,indent=2))
print(pd.DataFrame(bands).to_string(index=False))
