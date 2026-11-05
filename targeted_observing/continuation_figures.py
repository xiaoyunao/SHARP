"""Source-backed paired strategy figures and candidate rejection audit."""
import json
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .catalog import BASE,write_csv
from .continuation_report import ROOT,OUT,archive_csv
from .annual_trojan import load_tracks
from .annual_tasks import load_tasks


def main():
    s=pd.read_csv(OUT/'paired_strategy_summary.csv');x=np.arange(len(s));fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    axes[0].errorbar(x,s.quality_p50,yerr=[s.quality_p50-s.quality_p10,s.quality_p90-s.quality_p50],fmt='o',capsize=5,color='#245c94')
    axes[0].set_xticks(x,s.strategy);axes[0].set_ylabel('Quality-qualified objects, P10 / P50 / P90');axes[0].set_ylim(0,110)
    axes[0].axhline(101,color='#a4472b',ls='--');axes[0].text(.05,103,'Unchanged goal: 101',color='#a4472b')
    axes[0].set_title('Same parent sample, weather seeds and 2 h cap')
    axes[1].bar(x,s.special_hours_mean,color='#6e5c98');axes[1].set_xticks(x,s.strategy)
    axes[1].set_ylabel('Mean annual special-science held hours');axes[1].set_title('Included in the original Trojan allocation')
    for i,r in s.iterrows():axes[1].text(i,r.special_hours_mean+1,f'{r.special_actual_fraction_mean:.1%}',ha='center')
    fig.suptitle('200 paired annual simulations; quality compliance is not a recovered period')
    fig.tight_layout();fig.savefig(OUT/'strategy_comparison.png',dpi=170);plt.close(fig)
    t=pd.read_csv(OUT/'special_task_summary.csv',dtype={'object_id':str});ids=['17365','22180','4902','5511','5012','7543']
    matrix=t[t.object_id.isin(ids)].pivot(index='object_id',columns='strategy',values='sampling_ready_probability').reindex(ids)[list('ABCD')]
    fig,ax=plt.subplots(figsize=(8,4.8));im=ax.imshow(matrix.to_numpy(),vmin=0,vmax=1,cmap='Blues',aspect='auto')
    ax.set_xticks(range(4),list('ABCD'));ax.set_yticks(range(6),ids)
    for i in range(6):
        for j in range(4):ax.text(j,i,f'{matrix.iloc[i,j]:.1%}',ha='center',va='center',color='white' if matrix.iloc[i,j]>.6 else 'black')
    ax.set_title('Special-task SAMPLING readiness only\nPhysical/period solutions remain unassessed')
    fig.colorbar(im,ax=ax,label='Fraction of 200 simulated years');fig.tight_layout();fig.savefig(OUT/'science_sampling_readiness.png',dpi=170);plt.close(fig)
    tasks=load_tasks();tracks=load_tracks();avail={};fields_present={}
    for name in ['visibility','visibility_free']:
        counts=Counter();present=Counter()
        for p in sorted((ROOT/name).glob('202*.json')):
            menu=json.loads(p.read_text());data=np.load(p.with_suffix('.npz'))
            allids={k for f in menu['fields'] for k in f['core_targets']};present.update(allids)
            selected=np.flatnonzero(np.any(data['valid'],axis=1));counts.update({k for j in selected for k in menu['fields'][j]['core_targets']})
        avail[name]=counts;fields_present[name]=present
    audit=[]
    for k,t in tasks.items():
        n=avail['visibility_free'][k]
        reason='outside unchanged annual preselection' if k not in tracks else 'no candidate footprint at V<=19.2' if not fields_present['visibility_free'][k] else 'no whole-footprint/Moon/bright-star-safe two-hour interval' if not n else 'eligible; competition, task precision and budget decide actual sampling'
        audit.append({'object_id':k,'name':t['name'],'original_review_tier':t['tier'],'annual_parent_member':k in tracks,
            'nominal_V_min_proxy':min(tracks[k]['V']) if k in tracks else None,
            'anchor_menu_eligible_nights':avail['visibility'][k],'free_menu_eligible_nights':n,
            'selection_gate_reason':reason,'periods_h':json.dumps(t['periods_h']),'tags':json.dumps(t['tags'])})
    write_csv(OUT/'special_availability_audit.csv',audit)
    # Include zero-selected objects and all original exclusions in the selection function.
    parent=pd.read_csv(ROOT/'candidate_audit.csv',dtype={'object_id':str})
    selection=pd.read_csv(OUT/'object_selection_function.csv',dtype={'object_id':str});frames=[]
    for label in ['A','P0','B','C','D']:
        merged=parent.merge(selection[selection.strategy==label].drop(columns='strategy'),on='object_id',how='left')
        merged['strategy']=label
        for col in ['sampled_realizations','quality_probability']:merged[col]=merged[col].fillna(0)
        frames.append(merged)
    archive_csv(pd.concat(frames,ignore_index=True),OUT/'full_parent_selection_function.csv.gz')


if __name__=='__main__':main()
