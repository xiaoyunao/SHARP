"""Frozen Nov 12 audit: planned commands, explicit twilight replay, no TOO bridges."""
import csv
import hashlib
import json
from datetime import datetime
from zoneinfo import ZoneInfo

import numpy as np

from .catalog import BASE, write_csv
from .full_night_test import read_lines, movement
from .analyze_full_night_overhead import timestamp


def analyze():
    out = BASE / 'products/overhead_full_20261112'
    source = out / 'headers.json'
    script = BASE / 'products/full_night_20261112/20261112.txt'
    ledger = BASE / 'products/full_night_20261112/exposures.csv'
    data = json.loads(source.read_text())
    rows = data['rows']
    planned = list(csv.DictReader(ledger.open()))
    expanded = [c for c in read_lines(script.read_text(), '') for _ in range(c['count'])]
    assert len(planned) == len(expanded) == 430
    for p, c in zip(planned, expanded):
        assert p['field_id'] == c['field_id'] and p['filter'] == c['filter']
        assert float(p['exposure_s']) == c['exposure_s']
        assert abs(float(p['ra_deg'])-c['ra_deg']) < 1e-7
        assert abs(float(p['dec_deg'])-c['dec_deg']) < 1e-7
    plans = {k: [p for p in planned if p['block']==k] for k in dict.fromkeys(p['block'] for p in planned)}
    # Observed execution, not assumed completion: main cut after 65; YSO shortened.
    spans = [('EV',91,140,0), ('EV_REPLAY',141,190,0), ('TJN_EARLY',191,271,0),
             ('TJN_MAIN',272,331,0), ('TJN_MAIN',392,396,60),
             ('YSO_g',397,427,0), ('YSO_r',428,457,0), ('TJN_LATE',458,466,0), ('AM',467,473,0)]
    matched = []
    for block, lo, hi, offset in spans:
        rr = [r for r in rows if lo <= r['sequence'] <= hi]
        assert len(rr) == hi-lo+1
        pp = plans['EV' if block=='EV_REPLAY' else block]
        for index, r in enumerate(rr,offset):
            h=r['header'];p=pp[index]
            assert h['OBJECT']=='OBJ_'+p['field_id'] and h['FILTER']==p['filter']
            assert h['EXPREQ']==float(p['exposure_s'])
            assert abs((h['OBJ_RA']*15-float(p['ra_deg'])+180)%360-180)*3600 < 2
            assert abs(h['OBJ_DEC']-float(p['dec_deg']))*3600 < 2
            assert abs(timestamp(h['EXPEND'])-timestamp(h['EXPSTA'])-h['EXPOSURE'])<.003
            matched.append({**r,'block':block,'planned_index_in_block':index+1})
    assert len(matched)==len(rows)==323
    assert len({r['sequence'] for r in matched})==323
    excluded=[];pairs=[]
    for a,b in zip(matched,matched[1:]):
        x,y=a['header'],b['header']
        actual=timestamp(y['EXPSTA'])-timestamp(x['EXPEND'])
        assert actual>0
        if b['sequence']!=a['sequence']+1 or b['sequence']==141:
            excluded.append(dict(from_sequence=a['sequence'],to_sequence=b['sequence'],elapsed_s=actual,
                reason='TOO interruption' if b['sequence']==392 else 'Restart of complete evening list'))
            continue
        turn=movement(dict(ra_deg=x['OBJ_RA']*15,dec_deg=x['OBJ_DEC']),
                      dict(ra_deg=y['OBJ_RA']*15,dec_deg=y['OBJ_DEC']))
        change=10. if x['FILTER']!=y['FILTER'] else 0.
        same_group=a['block']==b['block'] or {a['block'],b['block']}=={'YSO_g','YSO_r'}
        category=('cross_block' if not same_group else 'filter_change' if change else
                  b['block'] if b['block'] in ('EV','EV_REPLAY','AM') else
                  'TJ_slew' if turn else 'TJ_repeat' if b['block'].startswith('TJN') else 'YSO_repeat')
        pairs.append(dict(from_sequence=a['sequence'],to_sequence=b['sequence'],category=category,
            from_field=x['OBJECT'],to_field=y['OBJECT'],from_filter=x['FILTER'],to_filter=y['FILTER'],
            previous_exposure_s=x['EXPREQ'],next_exposure_s=y['EXPREQ'],
            actual_overhead_s=actual,model_slew_s=turn,model_total_s=15+turn+change,
            residual_s=actual-15-turn-change))
    assert len(pairs)==320 and len(excluded)==2
    # Independent closure over the three continuous analyzed runs.
    for lo,hi in ((91,140),(141,331),(392,473)):
        rr=[r for r in matched if lo<=r['sequence']<=hi]
        elapsed=timestamp(rr[-1]['header']['EXPEND'])-timestamp(rr[0]['header']['EXPSTA'])
        total=sum(r['header']['EXPOSURE'] for r in rr)+sum(p['actual_overhead_s'] for p in pairs if lo<p['to_sequence']<=hi)
        assert abs(elapsed-total)<.01
    def stats(rr):
        v=np.array([p['actual_overhead_s'] for p in rr])
        return dict(count=len(rr),actual_mean_s=float(v.mean()),actual_median_s=float(np.median(v)),
            min_s=float(v.min()),max_s=float(v.max()),p10_s=float(np.quantile(v,.1)),p90_s=float(np.quantile(v,.9)),
            model_mean_s=float(np.mean([p['model_total_s'] for p in rr])),
            residual_sum_s=float(sum(p['residual_s'] for p in rr)))
    groups={k:stats([p for p in pairs if p['category']==k]) for k in dict.fromkeys(p['category'] for p in pairs)}
    blocks=[]
    for block in dict.fromkeys(r['block'] for r in matched):
        rr=[r for r in matched if r['block']==block]
        local=lambda h,k:datetime.fromtimestamp(timestamp(h[k]),ZoneInfo('Asia/Shanghai')).isoformat()
        blocks.append(dict(block=block,frames=len(rr),first_start_local=local(rr[0]['header'],'EXPSTA'),
            last_end_local=local(rr[-1]['header'],'EXPEND')))
    boundaries=[p for p in pairs if p['category'] in ('filter_change','cross_block')]
    summary=dict(night='20261112',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),ledger_sha256=hashlib.sha256(ledger.read_bytes()).hexdigest(),
        matched_frames=323,replayed_evening_frames=50,normal_pairs=320,excluded_files=data['excluded'],
        excluded_gaps=excluded,blocks=blocks,groups=groups,all_pairs=stats(pairs),boundaries=boundaries,
        largest_normal_gaps=sorted(pairs,key=lambda p:p['residual_s'],reverse=True)[:10],
        limitations=['Includes 50 exact-command evening replay frames, reported separately; restart wait excluded.',
            'No pair crosses TOO exposures; no TOO headers read.',
            'Effective end-to-next-start overhead includes controller and operator waits; not pure hardware time.',
            'FILTER used rather than stale CHANNEL; final readout and initial startup unavailable.',
            'Matching proves agreement with delivered commands; controller-side edited list is unavailable.'])
    write_csv(out/'pairs.csv',pairs);write_csv(out/'excluded_gaps.csv',excluded)
    write_csv(out/'matched_frames.csv',[dict(sequence=r['sequence'],block=r['block'],name=r['name'],
        planned_index_in_block=r['planned_index_in_block']) for r in matched])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(12,10),layout='constrained')
    for k in groups:
        pp=[p for p in pairs if p['category']==k]
        axes[0].scatter([p['to_sequence'] for p in pp],[p['actual_overhead_s'] for p in pp],s=13,label=k)
    axes[0].axvspan(332,391,color='gray',alpha=.15,label='TOO excluded')
    axes[0].set(xlabel='Exposure sequence',ylabel='Measured overhead (s)');axes[0].legend(ncol=4,fontsize=8)
    for k in ('TJ_repeat','YSO_repeat'):
        pp=[p for p in pairs if p['category']==k]
        axes[1].scatter([p['to_sequence'] for p in pp],[p['actual_overhead_s'] for p in pp],s=13,label=k)
    axes[1].axhline(15,color='black',ls='--',label='Model 15 s')
    axes[1].axhline(18,color='gray',ls=':',label='18 s reference')
    axes[1].set(xlabel='Exposure sequence',ylabel='Same-field overhead (s)',ylim=(14,26));axes[1].legend(ncol=4)
    axes[1].set_title('Zoom; longer gaps retained in top panel and totals',fontsize=10)
    xx=np.arange(len(boundaries))
    axes[2].bar(xx-.2,[p['actual_overhead_s'] for p in boundaries],.4,label='Measured')
    axes[2].bar(xx+.2,[p['model_total_s'] for p in boundaries],.4,label='Model')
    axes[2].set_xticks(xx,[f"{p['from_sequence']}→{p['to_sequence']}\n{p['from_filter']}→{p['to_filter']}" for p in boundaries])
    axes[2].set_ylabel('Boundary overhead (s)');axes[2].legend()
    for ax in axes:ax.grid(axis='y',alpha=.2)
    fig.suptitle('2026-11-12 night | 320 gaps | Evening replay labeled; restart and TOO excluded')
    fig.savefig(out/'overhead_full_night.png',dpi=160);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    analyze()
