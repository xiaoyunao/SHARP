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


from .plot_style import STRATEGIES, ERRORBAR, apply_style, save_pair


def plot_comparisons():
    apply_style()
    s = pd.read_csv(OUT/'paired_strategy_summary.csv')
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.3), layout='constrained')
    for i, row in s.iterrows():
        color = STRATEGIES[row.strategy]
        axes[0].errorbar(i, row.quality_p50,
            yerr=[[row.quality_p50-row.quality_p10], [row.quality_p90-row.quality_p50]],
            fmt='o', color=color, **ERRORBAR)
        axes[1].bar(i, row.special_hours_mean, color=color)
        axes[1].text(i, row.special_hours_mean+1, f'{row.special_actual_fraction_mean:.1%}', ha='center')
    for ax in axes:
        ax.set_xticks(np.arange(len(s)), s.strategy)
        ax.grid(axis='y')
    axes[0].set_ylabel('Qualified objects [count]')
    axes[0].set_title('P10 / P50 / P90')
    axes[0].set_ylim(0, 115)
    axes[0].axhline(101, color='0.45', ls='--', lw=1, label='Goal: 101')
    axes[0].legend()
    axes[1].set_ylabel('Special-science allocation [h/year]')
    axes[1].set_title('Mean and share of allocation')
    axes[1].set_ylim(0, max(s.special_hours_mean)*1.14)
    fig.suptitle('Baseline precision · 200 paired weather realizations')
    save_pair(fig, OUT/'strategy_comparison')
    plt.close(fig)

    t = pd.read_csv(OUT/'special_task_summary.csv', dtype={'object_id':str})
    ids = ['17365','22180','4902','5511','5012','7543','5284','617']
    matrix = t.pivot(index='object_id', columns='strategy', values='sampling_ready_probability').reindex(ids)[list('ABCD')]
    fig, ax = plt.subplots(figsize=(8.4, 6.3), layout='constrained')
    im = ax.imshow(matrix.to_numpy(), vmin=0, vmax=1, cmap='viridis', aspect='auto')
    ax.set_xticks(range(4), list('ABCD'))
    ax.set_yticks(range(len(ids)), ids)
    ax.set_ylabel('MPC number')
    ax.set_xlabel('Strategy')
    for i, values in enumerate(matrix.to_numpy()):
        for j, value in enumerate(values):
            ax.text(j, i, f'{value:.1%}', ha='center', va='center', color='black' if value>.57 else 'white')
    ax.set_title('Sampling readiness · baseline precision')
    fig.colorbar(im, ax=ax, label='Fraction of 200 simulated years')
    save_pair(fig, OUT/'science_sampling_readiness')
    plt.close(fig)


    cadence = pd.read_csv(OUT/'cadence_intervals.csv.gz')
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5), sharey=True, layout='constrained')
    for ax, kind, title in zip(axes, ['within_dwell', 'same_field_first_to_first'],
                               ['Within two-frame dwell', 'Full same-field revisit']):
        for label, group in cadence[cadence.kind == kind].groupby('scope'):
            times = np.sort(group.seconds.to_numpy())
            style = ':' if 'nominal' in label else '-' if label.startswith('B') else '--'
            ax.step(times, np.arange(1,len(times)+1)/len(times), where='post',
                    ls=style, color=STRATEGIES[label[0]], label=label.replace('baseline realization', 'weather'))
        ax.set_xlabel('Midpoint separation [s]')
        ax.set_title(title)
        ax.grid()
    axes[0].set_ylabel('Cumulative fraction')
    axes[1].legend(loc='upper left')
    save_pair(fig, OUT/'cadence_distribution')
    plt.close(fig)


def main():
    plot_comparisons()
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
