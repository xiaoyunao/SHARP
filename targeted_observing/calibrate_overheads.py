"""Three-night calibration from frozen headers and verified delivered commands.

No scheduler settings are changed. FITS timings measure effective overhead,
not separate hardware start/stop times. Mixed slew/filter pairs are not fitted.
"""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr, theilslopes

from .catalog import BASE, write_csv
from .full_night_test import read_lines, movement
from .analyze_full_night_overhead import timestamp

OUT = BASE/'products/overhead_three_nights_20261114'
NIGHTS = ('20261111', '20261112', '20261113')
READOUT_REFERENCE = 18.0


def describe(values):
    a = np.asarray(values, dtype=float)
    if not len(a):
        return {'n': 0}
    return dict(n=len(a), mean=float(a.mean()), median=float(np.median(a)),
                p10=float(np.quantile(a, .1)), p90=float(np.quantile(a, .9)),
                min=float(a.min()), max=float(a.max()), rms=float(np.sqrt(np.mean(a*a))))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def blocks_for(night):
    source = BASE/f'products/full_night_{night}/{night}.txt'
    expanded = []
    for i, c in enumerate(read_lines(source.read_text(), '')):
        expanded.extend([{**c, 'script_line': i+1}] * c['count'])
    if night == '20261111':
        middle = False
        for p in expanded:
            if p['field_id'].startswith('M31'):
                p['block'] = 'M31'; middle = True
            elif p['field_id'] == 'TJ_test': p['block'] = 'TJ'
            elif p['field_id'] == 'YSO_ORI': p['block'] = 'YSO'
            else: p['block'] = 'AM' if middle else 'EV'
    else:
        ledger = list(csv.DictReader((source.parent/'exposures.csv').open()))
        assert len(ledger) == len(expanded)
        for p, e in zip(expanded, ledger):
            assert p['field_id'] == e['field_id'] and p['filter'] == e['filter']
            assert p['exposure_s'] == float(e['exposure_s'])
            p['block'] = e['block']
    return {k: [p for p in expanded if p['block'] == k] for k in dict.fromkeys(p['block'] for p in expanded)}, source


def frames_for(night):
    if night == '20261113':
        source = OUT/'headers_20261113.json'
        snapshot = json.loads(source.read_text())
        mapping = []
        spans = [('EV',1,47,1),('TJN_EARLY',48,68,1),('TJN_EARLY',69,69,80),
                 ('TJN_MAIN',70,164,1),('YSO_g',165,233,1),('YSO_r',234,302,1),
                 ('TJN_LATE',303,312,1),('AM',313,342,1)]
        for block, lo, hi, offset in spans:
            mapping.extend(dict(sequence=k, block=block, planned_index_in_block=k-lo+offset) for k in range(lo,hi+1))
    else:
        source = BASE/f'products/overhead_full_{night}/headers.json'
        snapshot = json.loads(source.read_text())
        mapping = list(csv.DictReader((source.parent/'matched_frames.csv').open()))
    by_seq = {r['sequence']:r for r in snapshot['rows']}
    assert len(by_seq) == len(snapshot['rows'])
    blocks, script = blocks_for(night)
    result = []
    for m in mapping:
        seq = int(m['sequence']); block = m['block']
        idx = int(m.get('planned_index_in_block', m.get('plan_exposure_index')))
        p = blocks['EV' if block == 'EV_REPLAY' else block][idx-1]
        r = by_seq[seq]; h = r['header']
        assert h['OBJECT'] == 'OBJ_'+p['field_id']
        assert h['FILTER'] == p['filter'] and h['EXPREQ'] == p['exposure_s']
        assert abs((h['OBJ_RA']*15-p['ra_deg']+180)%360-180)*3600 < 2
        assert abs(h['OBJ_DEC']-p['dec_deg'])*3600 < 2
        assert abs(timestamp(h['EXPEND'])-timestamp(h['EXPSTA'])-h['EXPOSURE']) < .003
        assert h['READMODE'] == 'Image'
        result.append({**r, 'night':night, 'block':block, 'script_line':p['script_line'], 'planned_index':idx})
    assert len(result) == {'20261111':354,'20261112':323,'20261113':342}[night]
    fresh = json.loads((OUT/'headers_20261113.json').read_text())['inventories'][night]
    assert all(r['name'] in fresh for r in snapshot['rows'])
    old_names = {r['name'] for r in snapshot['rows']+snapshot['excluded']}
    additions = sorted(set(fresh)-old_names)
    assert all(n.startswith('CAL_') for n in additions), additions
    return result, dict(night=night, snapshot=str(source.relative_to(BASE)), snapshot_sha256=sha(source),
                        script=str(script.relative_to(BASE)), script_sha256=sha(script),
                        selected_frames=len(result), excluded_snapshot_files=len(snapshot['excluded']),
                        additional_non_script_files=additions)


def separation(ra1, dec1, ra2, dec2):
    # Haversine is stable at small angular separation. This is not the RA axis displacement.
    r1,d1,r2,d2=np.deg2rad([ra1,dec1,ra2,dec2])
    q=np.sin((d2-d1)/2)**2+np.cos(d1)*np.cos(d2)*np.sin((r2-r1)/2)**2
    return float(np.rad2deg(2*np.arcsin(np.sqrt(np.clip(q,0,1)))))


def assemble():
    pairs=[];excluded=[];provenance=[];frames=[]
    for night in NIGHTS:
        rr, prov = frames_for(night);provenance.append(prov);frames.extend(rr)
        for a,b in zip(rr,rr[1:]):
            x,y=a['header'],b['header'];gap=timestamp(y['EXPSTA'])-timestamp(x['EXPEND'])
            reason = None
            if b['sequence'] != a['sequence']+1: reason = 'TOO' if a['sequence'] in (207,331) else 'focus/off-script'
            if night=='20261112' and b['sequence']==141: reason='evening-list restart'
            if night=='20261113' and b['sequence']==69: reason='unidentified 2h15m interruption (not readout)'
            if reason:
                excluded.append(dict(night=night,from_sequence=a['sequence'],to_sequence=b['sequence'],gap_s=gap,reason=reason));continue
            ra1,ra2=x['OBJ_RA']*15,y['OBJ_RA']*15;dec1,dec2=x['OBJ_DEC'],y['OBJ_DEC']
            dra=abs((ra2-ra1+180)%360-180);ddec=abs(dec2-dec1)
            move=movement(dict(ra_deg=ra1,dec_deg=dec1),dict(ra_deg=ra2,dec_deg=dec2))
            change=x['FILTER']!=y['FILTER']
            same_exp=x['EXPREQ']==y['EXPREQ'];same_command=a['script_line']==b['script_line']
            kind=('slew_filter' if move and change else 'slew' if move else 'filter' if change else
                  'repeat' if same_exp and same_command else 'command_change')
            context='within_block' if a['block']==b['block'] else 'block_boundary'
            route=f"{ra1:.5f},{dec1:.5f}>{ra2:.5f},{dec2:.5f}"
            pairs.append(dict(night=night,from_sequence=a['sequence'],to_sequence=b['sequence'],
                from_field=x['OBJECT'],to_field=y['OBJECT'],from_block=a['block'],to_block=b['block'],
                from_filter=x['FILTER'],to_filter=y['FILTER'],previous_exposure_s=x['EXPREQ'],next_exposure_s=y['EXPREQ'],
                context=context,kind=kind,route=route,evening_replay=b['block']=='EV_REPLAY',
                gap_s=gap,dra_axis_deg=dra,ddec_axis_deg=ddec,sky_separation_deg=separation(ra1,dec1,ra2,dec2),
                max_axis_deg=max(dra,ddec),dominant_axis='RA' if 11.95*np.sqrt(dra)>=12.30*np.sqrt(ddec) else 'Dec',
                theory_slew_s=move,inferred_slew_s=gap-READOUT_REFERENCE if kind=='slew' else None,
                residual_s=gap-READOUT_REFERENCE-move if kind=='slew' else None,
                filter_added_s=gap-READOUT_REFERENCE if kind=='filter' else None,
                actual_next_start_utc=y['EXPSTA']))
        # Elapsed-time closure, including all explicitly excluded gaps.
        elapsed=timestamp(rr[-1]['header']['EXPEND'])-timestamp(rr[0]['header']['EXPSTA'])
        computed=sum(r['header']['EXPOSURE'] for r in rr)+sum(p['gap_s'] for p in pairs if p['night']==night)+sum(p['gap_s'] for p in excluded if p['night']==night)
        assert abs(elapsed-computed)<.02
    assert len(frames)==1019 and len(pairs)==1011 and len(excluded)==5
    return pairs,excluded,provenance,frames


def trend(routes):
    x=np.array([r['max_axis_deg'] for r in routes]);y=np.array([r['residual_median_s'] for r in routes])
    if len(x)<3 or len(set(x))<2:return dict(n=len(x))
    fit=theilslopes(y,x)
    return dict(n=len(x),spearman_rho=float(spearmanr(x,y).statistic),
        absolute_residual_spearman_rho=float(spearmanr(x,np.abs(y)).statistic),
        theil_sen_s_per_deg=float(fit.slope),intercept_s=float(fit.intercept),
        slope_low_s_per_deg=float(fit.low_slope),slope_high_s_per_deg=float(fit.high_slope),
        caution='Descriptive route medians; observational/repeated-route effects and boundary waits prevent a causal distance interpretation.')


def analyze():
    pairs,excluded,provenance,frames=assemble()
    readouts=[p for p in pairs if p['kind']=='repeat']
    filters=[p for p in pairs if p['kind']=='filter']
    slews=[p for p in pairs if p['kind']=='slew']
    routes=[]
    for key in dict.fromkeys(p['route'] for p in slews):
        pp=[p for p in slews if p['route']==key];first=pp[0]
        routes.append({k:first[k] for k in ('route','from_field','to_field','context','dra_axis_deg','ddec_axis_deg','max_axis_deg','sky_separation_deg','theory_slew_s','dominant_axis')}|
            dict(n=len(pp),nights=','.join(sorted({p['night'] for p in pp})),
                 residual_median_s=float(np.median([p['residual_s'] for p in pp])),
                 residual_min_s=min(p['residual_s'] for p in pp),residual_max_s=max(p['residual_s'] for p in pp)))
    bins=[]
    for lo,hi in ((0,3),(3,10),(10,50),(50,181)):
        pp=[p for p in slews if lo<=p['max_axis_deg']<hi]
        bins.append(dict(axis_range_deg=f'[{lo},{hi})',routes=len({p['route'] for p in pp}),
            **describe([p['residual_s'] for p in pp])))
    report=dict(nights=list(NIGHTS),frames=len(frames),pairs=len(pairs),excluded=excluded,
        per_night={n:dict(frames=sum(r['night']==n for r in frames),pairs=sum(p['night']==n for p in pairs),
                         kinds=dict(Counter(p['kind'] for p in pairs if p['night']==n))) for n in NIGHTS},
        readout=describe([p['gap_s'] for p in readouts]),
        readout_per_night={n:describe([p['gap_s'] for p in readouts if p['night']==n]) for n in NIGHTS},
        readout_by_exposure={str(e):describe([p['gap_s'] for p in readouts if p['next_exposure_s']==e]) for e in sorted({p['next_exposure_s'] for p in readouts})},
        readout_long_gaps=[p for p in readouts if p['gap_s']>30],
        filter_added=describe([p['filter_added_s'] for p in filters]),filter_events=filters,
        filter_by_direction={k:describe([p['filter_added_s'] for p in filters if p['from_filter']+'>'+p['to_filter']==k]) for k in dict.fromkeys(p['from_filter']+'>'+p['to_filter'] for p in filters)},
        readout_reference_s=READOUT_REFERENCE,
        representative_filter_added_s=90.0,
        filter_single_value_basis='Six-event median 88.3595 s rounded to 90 s; includes filter-associated control/wait time, not a hardware constant.',
        slew_all=describe([p['residual_s'] for p in slews]),
        slew_within=describe([p['residual_s'] for p in slews if p['context']=='within_block']),
        slew_boundary=describe([p['residual_s'] for p in slews if p['context']=='block_boundary']),
        slew_per_night={n:describe([p['residual_s'] for p in slews if p['night']==n and p['context']=='within_block']) for n in NIGHTS},
        slew_by_dominant_axis={k:describe([p['residual_s'] for p in slews if p['dominant_axis']==k and p['context']=='within_block']) for k in ('RA','Dec')},
        slew_without_replay=describe([p['residual_s'] for p in slews if p['context']=='within_block' and not p['evening_replay']]),
        distance_bins=bins,trend_all_routes=trend(routes),trend_within_routes=trend([r for r in routes if r['context']=='within_block']),
        unique_routes=len(routes),source_provenance=provenance,
        decision='Evaluation only: keep current readout 18s, filter 10s, and slew code unchanged until user requests adjustment.')
    write_csv(OUT/'all_pairs.csv',pairs);write_csv(OUT/'slew_pairs.csv',slews);write_csv(OUT/'routes.csv',routes)
    write_csv(OUT/'filter_events.csv',filters);write_csv(OUT/'excluded_intervals.csv',excluded)
    write_csv(OUT/'matched_frames.csv',[dict(night=r['night'],sequence=r['sequence'],name=r['name'],block=r['block'],script_line=r['script_line'],planned_index=r['planned_index']) for r in frames])
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    plots(pairs,routes,report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('filter_events','source_provenance','readout_long_gaps')},indent=2))


def plots(pairs,routes,report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    colors=['#0072B2','#D55E00','#009E73'];markers=['o','s','^']
    # 1. Readout cluster and every same-position filter change; no filter outlier removal.
    fig,ax=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for i,n in enumerate(NIGHTS):
        v=np.array([p['gap_s'] for p in pairs if p['kind']=='repeat' and p['night']==n])
        ax[0].scatter(np.full(len(v),i)+(np.arange(len(v))%13-6)*.018,v,s=9,alpha=.35,color=colors[i],marker=markers[i])
        ax[0].plot(i,np.median(v),'k_',ms=22,mew=2)
    ax[0].set(xticks=range(3),xticklabels=['Nov 11','Nov 12','Nov 13'],ylabel='Same-command inter-exposure gap (s)',ylim=(16,25),title='Fixed-field baseline: about 18 s')
    ax[0].axhline(18,color='gray',ls='--');ax[0].text(.02,.97,f"{report['readout']['n']} pairs; {len(report['readout_long_gaps'])} gaps >30 s outside zoom",transform=ax[0].transAxes,va='top',fontsize=9)
    ff=[p for p in pairs if p['kind']=='filter'];xx=np.arange(len(ff))
    ax[1].bar(xx,[p['filter_added_s'] for p in ff],color=[colors[NIGHTS.index(p['night'])] for p in ff])
    for i,p in enumerate(ff):ax[1].text(i,p['filter_added_s']+2,f"{p['filter_added_s']:.1f}",ha='center',fontsize=9)
    ax[1].set(xticks=xx,xticklabels=[f"{p['night'][-4:]}\n{p['from_filter']}→{p['to_filter']}" for p in ff],ylabel='Additional gap after subtracting 18 s (s)',ylim=(0,125),title='Filter change: six events, substantial scatter')
    ax[1].axhline(report['filter_added']['median'],color='gray',ls='--',label='Median');ax[1].legend()
    for a in ax:a.grid(axis='y',alpha=.2)
    fig.savefig(OUT/'01_readout_filter.png',dpi=180);plt.close(fig)
    # 2. Full distances + local-range inset panel. Large boundaries remain visible.
    slews=[p for p in pairs if p['kind']=='slew']
    fig,ax=plt.subplots(1,2,figsize=(13,5.5),layout='constrained')
    for i,n in enumerate(NIGHTS):
        for context in ('within_block','block_boundary'):
            pp=[p for p in slews if p['night']==n and p['context']==context]
            style=dict(marker=markers[i],s=25 if context=='within_block' else 75,
                       facecolors=colors[i] if context=='within_block' else 'none',edgecolors=colors[i],alpha=.6)
            ax[0].scatter([p['max_axis_deg'] for p in pp],[p['residual_s'] for p in pp],label=f"{n[-4:]} {'loop' if context=='within_block' else 'boundary'}",**style)
            if context=='within_block':ax[1].scatter([p['max_axis_deg'] for p in pp],[p['residual_s'] for p in pp],**style)
    for r in routes:
        if r['context']=='within_block':ax[1].plot(r['max_axis_deg'],r['residual_median_s'],'k_',ms=14,mew=2)
    for a in ax:
        a.axhline(0,color='black',ls='--',lw=1);a.set(xlabel='Maximum RA/Dec axis displacement (deg)',ylabel='Inferred slew − theory (s)');a.grid(alpha=.2)
    ax[0].set_title('All same-filter slews; open marks = block boundaries');ax[0].legend(fontsize=8,ncol=2)
    ax[1].set_title('Routine loops; black ticks = route medians')
    ax[1].set_xlim(1.4,5.0)
    fig.suptitle('Inferred slew = measured inter-exposure gap − 18 s; no mixed filter/slew events')
    fig.savefig(OUT/'02_slew_residual_distance.png',dpi=180);plt.close(fig)
    # 3. Theory-vs-effective measured time, and axis geometry of each distinct route.
    fig,ax=plt.subplots(1,2,figsize=(13,5.5),layout='constrained')
    for i,n in enumerate(NIGHTS):
        for context in ('within_block','block_boundary'):
            pp=[p for p in slews if p['night']==n and p['context']==context]
            ax[0].scatter([p['theory_slew_s'] for p in pp],[p['inferred_slew_s'] for p in pp],s=30,
                facecolors=colors[i] if context=='within_block' else 'none',edgecolors=colors[i],marker=markers[i],label=f'{n[-4:]} {context}')
    maxval=max(max(p['inferred_slew_s'],p['theory_slew_s']) for p in slews)*1.06
    ax[0].plot([0,maxval],[0,maxval],'k--',lw=1,label='1:1');ax[0].set(xlim=(0,maxval),ylim=(0,maxval),xlabel='Theoretical slew (s)',ylabel='Inferred slew: gap − 18 s (s)',title='Long transfers occur at block boundaries')
    ax[0].legend(fontsize=7,ncol=2)
    rr=[r for r in routes if r['context']=='within_block'];lim=max(abs(r['residual_median_s']) for r in rr)
    sc=ax[1].scatter([r['dra_axis_deg'] for r in rr],[r['ddec_axis_deg'] for r in rr],s=[40+6*r['n'] for r in rr],
        c=[r['residual_median_s'] for r in rr],cmap='coolwarm',norm=TwoSlopeNorm(0,-lim,lim),edgecolors='black',linewidths=.5,alpha=.8)
    z=np.linspace(0,5,100);ax[1].plot(z,(11.95/12.30)**2*z,'k--',lw=1,label='Equal modeled axis time')
    ax[1].set(xlim=(0,5),ylim=(0,5),xlabel='RA axis displacement (deg; no cos Dec)',ylabel='Dec axis displacement (deg)',title='Routine routes; marker area scales with count')
    fig.colorbar(sc,ax=ax[1],label='Route median residual (s)');ax[1].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2)
    fig.savefig(OUT/'03_slew_theory_axes.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    analyze()
