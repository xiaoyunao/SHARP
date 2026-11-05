"""Fail-closed paired annual review tables from completed simulations."""
import gzip
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from .catalog import BASE,write_csv

ROOT=BASE/'products/annual'
OUT=BASE/'docs/continuation_20261105'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def archive_csv(frame,path):
    with path.open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as stream:
            stream.write(frame.to_csv(index=False).encode())


def ledger(name):
    d=pd.read_csv(ROOT/name/'nightly_execution_ledger.csv')
    if len(d)!=365*200 or d.duplicated(['realization','date']).any():raise ValueError('Incomplete or duplicate ledger '+name)
    cols=['exposure_s','readout_s','slew_s','filter_s','delay_s','handover_s','slack_s']
    if np.max(np.abs(d.allocated_s-d[cols].sum(axis=1)))>1e-5:raise ValueError('Time ledger does not close '+name)
    if (d.allocated_s>7200).any() or (d.slack_s<-.001).any():raise ValueError('Invalid budget '+name)
    return d


def compare_realization0(a,b):
    def read(name):
        d=pd.read_csv(ROOT/name/'execution_realization0.csv')
        return d.groupby('date').agg(start=('block_start_unix','first'),points=('field_id','size'),
            fields=('field_id',lambda x:json.dumps(sorted(set(x)))))
    x=read(a).join(read(b),how='outer',lsuffix='_a',rsuffix='_b')
    both=x.start_a.notna()&x.start_b.notna()
    return {'realization':0,'both_active_nights':int(both.sum()),
            'changed_start_both_active_nights':int((both&(abs(x.start_a-x.start_b)>.01)).sum()),
            'changed_active_state_nights':int((x.start_a.isna()!=x.start_b.isna()).sum()),
            'changed_field_sets_nights':int((x.fields_a.fillna('[]')!=x.fields_b.fillna('[]')).sum()),
            'actual_exposures_before':int(x.points_a.sum()),'actual_exposures_after':int(x.points_b.sum())}


def main():
    OUT.mkdir(exist_ok=True);summaries=[];paired=[];special=[];objects=[];inputs={};configs={}
    names={'A':'continuation_A_baseline','P0':'queue_4fields_dwell2_baseline_science25_v4',
           **{x:'continuation_'+x+'_baseline' for x in 'BCD'}}
    base=json.loads((ROOT/names['A']/'summary.json').read_text());old=json.loads((ROOT/'queue_4fields_dwell2_baseline_science25_v3/summary.json').read_text())
    if base['counts']!=old['counts']:raise ValueError('A no longer reproduces frozen baseline')
    a_ledger=ledger(names['A']).set_index(['realization','date']);a_counts=np.array(base['counts'])
    for label,name in names.items():
        root=ROOT/name;s=json.loads((root/'summary.json').read_text());counts=np.array(s['counts'])
        if len(counts)!=200 or s['nights']!=365:raise ValueError('Incomplete annual summary '+name)
        d=ledger(name);o=pd.read_csv(root/'object_yield.csv',dtype={'object_id':str})
        if o.duplicated(['realization','object_id']).any():raise ValueError('Object duplicate')
        calculated=o[o.quality_lightcurve].groupby('realization').size().reindex(range(200),fill_value=0).to_numpy()
        np.testing.assert_array_equal(counts,calculated)
        total=d.groupby('realization').allocated_s.sum();sp=d.groupby('realization').special_allocation_s.sum() if 'special_allocation_s' in d else total*0
        cap=s.get('special_budget_fraction',0.)
        if (sp>cap*total+1e-5).any():raise ValueError('Special budget exceeded')
        if 'special_allocation_s' in d:
            for _,rows in d.groupby('realization'):
                if (rows.special_allocation_s.cumsum()>cap*rows.allocated_s.cumsum()+1e-5).any():raise ValueError('Prefix quota exceeded')
        summary={'strategy':label,'quality_p10':np.percentile(counts,10),'quality_p50':np.median(counts),'quality_p90':np.percentile(counts,90),
                 'quality_mean':counts.mean(),'ge101_realizations':int((counts>=101).sum()),
                 'paired_delta_mean':float((counts-a_counts).mean()),'paired_delta_p10':np.percentile(counts-a_counts,10),
                 'paired_delta_p50':np.median(counts-a_counts),'paired_delta_p90':np.percentile(counts-a_counts,90),
                 'sampled_objects_mean':o.groupby('realization').size().mean(),
                 'raw_points_mean':o.groupby('realization').raw_points_total.sum().mean(),
                 'exposures_mean':d.groupby('realization').exposures.sum().mean(),
                 'allocated_hours_mean':total.mean()/3600,'special_quota':cap,'special_hours_mean':sp.mean()/3600,
                 'special_actual_fraction_mean':float((sp/total).fillna(0).mean()),
                 **{k:v for k,v in compare_realization0(names['A'],name).items() if k!='realization'}}
        if 'selected_object_ids' in d:
            b=d.set_index(['realization','date']);both=a_ledger.block_start_unix.notna()&b.block_start_unix.notna()
            summary['changed_start_all200_both_active']=int((both&((a_ledger.block_start_unix-b.block_start_unix).abs()>.01)).sum())
            summary['changed_selected_objects_all200']=int((a_ledger.selected_object_ids!=b.selected_object_ids).sum())
            # Compact full decision ledger enables inspection without ignored products.
            archive_csv(d,OUT/f'{label}_nightly_ledger.csv.gz')
        summaries.append(summary)
        paired.extend({'strategy':label,'realization':i,'quality_count':int(v),'delta_vs_A':int(v-a_counts[i])} for i,v in enumerate(counts))
        taskpath=root/'science_task_status.csv'
        if taskpath.exists():
            task=pd.read_csv(taskpath,dtype={'object_id':str})
            archive_csv(task,OUT/f'{label}_science_tasks.csv.gz')
            for k,g in task.groupby('object_id'):
                special.append({'strategy':label,'object_id':k,'name':g.name.iloc[0],
                    'photometric_probability':g.photometric_complete.mean(),'sampling_ready_probability':g.sampling_ready.mean(),
                    'task_nights_median':g.task_nights.median(),'task_baseline_days_median':g.task_baseline_days.median(),
                    'phase_coverage_median':g.min_phase_coverage.median(),'selection_reasons':json.dumps(g.selection_reason.value_counts().to_dict()),
                    'period_solution_complete':None,'science_case_complete':None})
        roles=json.loads((root/'science_roles.json').read_text())
        for k,g in o.groupby('object_id'):
            objects.append({'strategy':label,'object_id':k,'sampled_realizations':len(g),'quality_probability':g.quality_lightcurve.sum()/200,
                            'role':roles.get(k,{}).get('role','unchecked'),'mean_raw_points_when_sampled':g.raw_points_total.mean()})
        for file in ['summary.json','nightly_execution_ledger.csv','object_yield.csv','execution_realization0.csv','science_task_status.csv','run_config.json']:
            path=root/file
            if path.exists():inputs[name+'/'+file]=sha(path)
        cfg=root/'run_config.json'
        if cfg.exists():configs[label]=json.loads(cfg.read_text())
    write_csv(OUT/'paired_strategy_summary.csv',summaries);write_csv(OUT/'paired_realizations.csv',paired)
    write_csv(OUT/'special_task_summary.csv',special);write_csv(OUT/'object_selection_function.csv',objects)
    (OUT/'strategy_configs.json').write_text(json.dumps(configs,indent=2)+'\n')
    optimistic=[]
    for scenario in ['baseline','optimistic']:
        before='queue_4fields_dwell2_'+scenario+'_science25_v3';after=before[:-1]+'4'
        x=json.loads((ROOT/before/'summary.json').read_text());y=json.loads((ROOT/after/'summary.json').read_text())
        delta=np.array(y['counts'])-x['counts']
        optimistic.append({'scenario':scenario,'before':x['quality_count_P10_P50_P90'],'after':y['quality_count_P10_P50_P90'],
            'paired_delta_mean':float(delta.mean()),'paired_delta_p10_p50_p90':np.percentile(delta,[10,50,90]).tolist(),
            **compare_realization0(before,after)})
    (OUT/'P0_paired_comparison.json').write_text(json.dumps(optimistic,indent=2)+'\n')
    source_files=list(BASE.glob('annual_*.py'))+list(BASE.glob('continuation_*.py'))+[BASE/n for n in ['validate.py','geometry.py','scheduler.py','science.py','bright_stars.py']]
    sources={str(p.relative_to(BASE)):sha(p) for p in source_files}
    sources.update({str(p.relative_to(BASE)):sha(p) for p in [BASE/'annual_defaults.json',BASE/'defaults.json',BASE/'data/special_science_candidates.json',BASE/'data/weather_prior.json']})
    (OUT/'comparison_provenance.json').write_text(json.dumps({'baseline_commit':'a5cb9454e56771ebe4211b7a544412056cce4c40',
        'scope':'all strategies 365 nights x 200 paired weather realizations; P0 same menus, BCD additional free centres',
        'delivery_source_sha256':sources,'input_product_sha256':inputs,'frozen_A_counts_reproduced':True,
        'note':'Git commit containing these source hashes identifies the new code; original baseline plot manifest remains a5cb945.'},indent=2)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
