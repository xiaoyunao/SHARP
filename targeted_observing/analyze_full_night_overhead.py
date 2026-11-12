"""Audit Nov 11 planned observations; never bridge focus/TOO interruptions."""
import hashlib
import json
from datetime import datetime, timezone
from zoneinfo import ZoneInfo

import numpy as np

from .catalog import BASE, write_csv
from .full_night_test import read_lines, movement


def timestamp(value):
    return datetime.fromisoformat(value).replace(tzinfo=timezone.utc).timestamp()


def analyze():
    out = BASE / 'products/overhead_full_20261111'
    source = out / 'headers.json'
    plan_path = BASE / 'products/full_night_20261111/20261111.txt'
    snapshot = json.loads(source.read_text())
    rows = sorted(snapshot['rows'], key=lambda r: r['sequence'])
    plans = {k: [] for k in ('EV', 'M31', 'TJ', 'YSO', 'AM')}
    passed_m31 = False
    for command in read_lines(plan_path.read_text(), ''):
        name = command['field_id']
        if name.startswith('M31'):
            block = 'M31'; passed_m31 = True
        elif name == 'TJ_test': block = 'TJ'
        elif name == 'YSO_ORI': block = 'YSO'
        else:
            assert name.startswith('MP_')
            block = 'AM' if passed_m31 else 'EV'
        plans[block].extend([command] * command['count'])
    ranges = {'EV': [(31, 79)], 'M31': [(94, 197)],
              'TJ': [(198, 207), (238, 260)], 'YSO': [(261, 404)], 'AM': [(405, 428)]}
    matched = []
    for block, spans in ranges.items():
        selected = [r for r in rows if any(lo <= r['sequence'] <= hi for lo, hi in spans)]
        assert len(selected) == sum(hi-lo+1 for lo, hi in spans)
        for index, row in enumerate(selected):
            p = plans[block][index]; h = row['header']
            assert h['OBJECT'] == 'OBJ_' + p['field_id']
            assert h['FILTER'] == p['filter'] and h['EXPREQ'] == p['exposure_s']
            assert abs((h['OBJ_RA']*15-p['ra_deg']+180)%360-180)*3600 < 2
            assert abs(h['OBJ_DEC']-p['dec_deg'])*3600 < 2
            assert abs(timestamp(h['EXPEND'])-timestamp(h['EXPSTA'])-h['EXPOSURE']) < .003
            matched.append({**row, 'block': block, 'plan_exposure_index': index+1})
    matched.sort(key=lambda r:r['sequence'])
    assert len(matched) == 354
    unmatched = [r for r in rows if r['sequence'] not in {m['sequence'] for m in matched}]
    assert [r['sequence'] for r in unmatched] == [80]
    assert unmatched[0]['header']['OBJECT'] == matched[48]['header']['OBJECT']
    assert unmatched[0]['header']['OBJECT'] != 'OBJ_' + plans['EV'][49]['field_id']
    pairs = []; interruptions = []
    for a, b in zip(matched, matched[1:]):
        x, y = a['header'], b['header']
        actual = timestamp(y['EXPSTA'])-timestamp(x['EXPEND'])
        if b['sequence'] != a['sequence']+1:
            interruptions.append(dict(from_sequence=a['sequence'], to_sequence=b['sequence'],
                elapsed_s=actual, reason='Off-script retry/focus' if a['sequence']==79 else 'TOO interruption'))
            continue
        before = dict(ra_deg=x['OBJ_RA']*15, dec_deg=x['OBJ_DEC'])
        after = dict(ra_deg=y['OBJ_RA']*15, dec_deg=y['OBJ_DEC'])
        turn = movement(before, after)
        change = 10. if x['FILTER'] != y['FILTER'] else 0.
        category = ('cross_block' if a['block'] != b['block'] else
                    'filter_change' if change else
                    a['block'] + ('_slew' if turn else '_repeat'))
        pairs.append(dict(from_sequence=a['sequence'], to_sequence=b['sequence'],
            category=category, from_field=x['OBJECT'], to_field=y['OBJECT'],
            from_filter=x['FILTER'], to_filter=y['FILTER'],
            previous_exposure_s=x['EXPREQ'], next_exposure_s=y['EXPREQ'],
            exposure_setting_changed=x['EXPREQ'] != y['EXPREQ'],
            actual_overhead_s=actual, model_slew_s=turn, model_filter_s=change,
            model_total_s=15.+turn+change, residual_s=actual-15.-turn-change))
    assert len(pairs) == 351 and len(interruptions) == 2
    # Independent elapsed-time closure for each uninterrupted run, including final exposure.
    for lo, hi in ((31,79), (94,207), (238,428)):
        run = [r for r in matched if lo <= r['sequence'] <= hi]
        measured = timestamp(run[-1]['header']['EXPEND'])-timestamp(run[0]['header']['EXPSTA'])
        accounted = sum(r['header']['EXPOSURE'] for r in run) + sum(
            p['actual_overhead_s'] for p in pairs if lo < p['to_sequence'] <= hi)
        assert abs(measured-accounted) < .01
    def stats(rr):
        v = np.array([r['actual_overhead_s'] for r in rr])
        return dict(count=len(rr), actual_mean_s=float(v.mean()), actual_median_s=float(np.median(v)),
            actual_min_s=float(v.min()), actual_max_s=float(v.max()),
            actual_p10_s=float(np.quantile(v,.1)), actual_p90_s=float(np.quantile(v,.9)),
            model_mean_s=float(np.mean([r['model_total_s'] for r in rr])),
            residual_sum_s=float(sum(r['residual_s'] for r in rr)))
    groups = {k:stats([r for r in pairs if r['category']==k]) for k in dict.fromkeys(r['category'] for r in pairs)}
    same = [r for r in pairs if r['category'].endswith('_repeat')]
    settings = {str(changed): stats([r for r in same if r['exposure_setting_changed']==changed]) for changed in (False,True)}
    blocks = []
    for block in plans:
        rr = [r for r in matched if r['block']==block]
        local = lambda h,key: datetime.fromtimestamp(timestamp(h[key]),ZoneInfo('Asia/Shanghai')).isoformat()
        blocks.append(dict(block=block, planned_frames=len(plans[block]), matched_frames=len(rr),
            first_start_local=local(rr[0]['header'],'EXPSTA'), last_end_local=local(rr[-1]['header'],'EXPEND')))
    summary = dict(night='20261111', source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        plan_sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest(), matched_frames=len(matched),
        normal_pairs=len(pairs), groups=groups, repeat_by_exposure_setting_change=settings,
        blocks=blocks, interruptions_excluded=interruptions, all_pairs=stats(pairs),
        extra_frame_excluded=dict(sequence=80, reason='Repeated MP_1831, while next planned exposure is MP_1832'),
        other_excluded_files=len(snapshot['excluded']),
        limitations=['End-to-next-start is total effective overhead, not separate readout/slew/filter hardware timing.',
            'Internal waits remain included; causes cannot be inferred from FITS headers.',
            'First startup and final readout are not measured. No pairs bridge excluded exposures.',
            'FILTER used; CHANNEL is stale for non-r exposures.'])
    write_csv(out/'pairs.csv', pairs)
    write_csv(out/'interruptions_excluded.csv',interruptions)
    write_csv(out/'matched_frames.csv',[dict(sequence=r['sequence'],block=r['block'],name=r['name'],
        plan_exposure_index=r['plan_exposure_index']) for r in matched])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3,1,figsize=(12,10),layout='constrained')
    for category in groups:
        pp=[r for r in pairs if r['category']==category]
        axes[0].scatter([r['to_sequence'] for r in pp],[r['actual_overhead_s'] for r in pp],s=12,label=category)
    axes[0].set(xlabel='Exposure sequence',ylabel='Measured overhead (s)')
    axes[0].legend(ncol=4,fontsize=8)
    for left,right in ((80,93),(208,237)):
        axes[0].axvspan(left,right,color='gray',alpha=.15)
    for block in ('M31','TJ','YSO'):
        rr=[r for r in same if r['category']==block+'_repeat']
        axes[1].scatter([r['to_sequence'] for r in rr],[r['actual_overhead_s'] for r in rr],s=12,label=block)
    axes[1].axhline(15,color='black',ls='--',label='Model 15 s')
    axes[1].set(xlabel='Exposure sequence',ylabel='Same-field overhead (s)',ylim=(14,30))
    axes[1].set_title(f"Zoom: {sum(r['actual_overhead_s'] > 30 for r in same)} same-field gaps exceed 30 s (retained in totals)",fontsize=10)
    axes[1].legend(ncol=4,fontsize=9)
    boundaries=[r for r in pairs if r['category'] in ('filter_change','cross_block')]
    xx=np.arange(len(boundaries))
    axes[2].bar(xx-.2,[r['actual_overhead_s'] for r in boundaries],.4,label='Measured')
    axes[2].bar(xx+.2,[r['model_total_s'] for r in boundaries],.4,label='Model')
    axes[2].set_xticks(xx,[f"{r['from_sequence']}→{r['to_sequence']}\n{r['from_filter']}→{r['to_filter']}" for r in boundaries])
    axes[2].set_ylabel('Boundary overhead (s)');axes[2].legend()
    for ax in axes: ax.grid(axis='y',alpha=.2)
    fig.suptitle('2026-11-11 night | 354 planned frames, 351 adjacent pairs | TOO/focus gaps excluded')
    fig.savefig(out/'overhead_full_night.png',dpi=160);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    analyze()
