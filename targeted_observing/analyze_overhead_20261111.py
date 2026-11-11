"""Compare frozen 0030--0079 headers with the existing offline overhead model."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from .catalog import BASE, write_csv
from .offline_twilight import slew


def timestamp(value):
    return datetime.fromisoformat(value).replace(tzinfo=timezone.utc).timestamp()


def analyze():
    out=BASE/'products/overhead_20261111'
    source=out/'headers_0030_0079.json'
    rows=sorted(json.loads(source.read_text())['rows'],key=lambda r:r['sequence'])
    assert [r['sequence'] for r in rows]==list(range(30,80))
    records=[];excluded=[]
    for a,b in zip(rows,rows[1:]):
        x,y=a['header'],b['header']
        observed=timestamp(y['EXPSTA'])-timestamp(x['EXPEND'])
        assert observed>0
        assert abs(timestamp(x['EXPEND'])-timestamp(x['EXPSTA'])-x['EXPOSURE'])<.003
        if any(k not in h for h in (x,y) for k in ('OBJ_RA','OBJ_DEC')):
            excluded.append(dict(from_sequence=a['sequence'],to_sequence=b['sequence'],
                observed_overhead_s=observed,reason='Startup image has no target coordinates; slew/idle split unknown'))
            continue
        before=dict(field_id=x['OBJECT'],ra_deg=x['OBJ_RA']*15,dec_deg=x['OBJ_DEC'])
        after=dict(field_id=y['OBJECT'],ra_deg=y['OBJ_RA']*15,dec_deg=y['OBJ_DEC'])
        assert x['FILTER']==y['FILTER']=='r'
        move=slew(before,after)
        records.append(dict(from_sequence=a['sequence'],to_sequence=b['sequence'],
            from_field=x['OBJECT'],to_field=y['OBJECT'],previous_exposure_s=x['EXPOSURE'],
            actual_overhead_s=observed,model_slew_s=move,model_readout_s=15.,
            model_total_s=15+move,residual_s=observed-15-move,
            dra_deg=abs((after['ra_deg']-before['ra_deg']+180)%360-180),
            ddec_deg=abs(after['dec_deg']-before['dec_deg'])))
    assert len(records)==48 and len(excluded)==1
    actual=np.array([r['actual_overhead_s'] for r in records])
    predicted=np.array([r['model_total_s'] for r in records]);residual=actual-predicted
    groups=[]
    for route in dict.fromkeys((r['from_field'],r['to_field']) for r in records):
        rr=[r for r in records if (r['from_field'],r['to_field'])==route]
        groups.append(dict(from_field=route[0],to_field=route[1],count=len(rr),
            actual_mean_s=float(np.mean([r['actual_overhead_s'] for r in rr])),
            model_mean_s=float(np.mean([r['model_total_s'] for r in rr])),
            residual_mean_s=float(np.mean([r['residual_s'] for r in rr]))))
    observed_elapsed=timestamp(rows[-1]['header']['EXPSTA'])-timestamp(rows[1]['header']['EXPSTA'])
    integration=sum(r['header']['EXPOSURE'] for r in rows[1:-1])
    assert abs(observed_elapsed-integration-actual.sum())<.003
    original=BASE/'products/tonight_full_hour_20261111/plan.json'
    plan=json.loads(original.read_text())['blocks'][0]['rows']
    for r,p in zip(rows[1:],plan):
        assert r['header']['OBJECT'].replace('OBJ_MP_','MP')==p['field_id']
        assert r['header']['EXPREQ']==p['exposure_s']
    summary=dict(header_snapshot_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        frames=50,science_frames=49,normal_pairs=48,startup_pairs=excluded,
        science_first_utc=rows[1]['header']['EXPSTA'],science_last_end_utc=rows[-1]['header']['EXPEND'],
        actual_mean_s=float(actual.mean()),model_mean_s=float(predicted.mean()),
        actual_median_s=float(np.median(actual)),model_median_s=float(np.median(predicted)),
        residual_mean_s=float(residual.mean()),residual_median_s=float(np.median(residual)),
        residual_min_s=float(residual.min()),residual_max_s=float(residual.max()),
        residual_rmse_s=float(np.sqrt(np.mean(residual**2))),
        fractional_mean_residual=float(residual.mean()/predicted.mean()),
        actual_sum_s=float(actual.sum()),model_sum_s=float(predicted.sum()),
        residual_sum_s=float(residual.sum()),
        residual_mean_if_readout_1447_s=float(residual.mean()+.53),
        early_5_to_25s_pairs=int(sum(r['previous_exposure_s']<29 for r in records)),
        last_sample_end_not_readout=True,
        limitations=['End-to-next-start includes readout plus slew/settle/controller overhead; not pure slew.',
          'No filter changes or fixed-field repeat pairs in this sample; cannot separately measure readout/filter change.',
          'Only nine short local routes; no evidence on large block-to-block moves or long 120s repeat cycles.',
          '0030 has no coordinates; startup interval is excluded from model residual.',
          '50th planned MP frame is outside requested sample; do not infer its status.'])
    write_csv(out/'pairs.csv',records);write_csv(out/'routes.csv',groups)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True,layout='constrained')
    seq=[r['to_sequence'] for r in records]
    axes[0].plot(seq,actual,'o-',ms=3,label='Measured: next EXPSTA - previous EXPEND')
    axes[0].plot(seq,predicted,'--',label='Model: readout 15 s + slew')
    axes[0].set_ylabel('Overhead (s)');axes[0].legend(fontsize=9)
    axes[1].plot(seq,residual,'o-',ms=3,color='tab:purple')
    axes[1].axhline(0,color='gray',ls='--');axes[1].set_ylabel('Measured - model (s)')
    axes[1].set_xlabel('Next exposure sequence')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('2026-11-11 | 48 normal transitions | startup excluded')
    fig.savefig(out/'overhead_comparison.png',dpi=160);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    analyze()
