#!/usr/bin/env python3
"""Compare all accepted sources and registered common stars, with spatial checks."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


def read_catalog(path):
    return np.genfromtxt(path,delimiter=',',names=True)


def match(ref,cat):
    xy=lambda c:np.column_stack((c['x'],c['y']))
    r,c=xy(ref),xy(cat)
    bright_r=r[np.argsort(ref['flux'])[-150:]]
    bright_c=c[np.argsort(cat['flux'])[-150:]]
    delta=(bright_c[:,None,:]-bright_r[None,:,:]).reshape(-1,2)
    bins=np.round(delta/4).astype(int)
    un,cnt=np.unique(bins,axis=0,return_counts=True)
    guess=un[np.argmax(cnt)]*4
    nearby=np.linalg.norm(delta-guess,axis=1)<6
    shift=np.median(delta[nearby],axis=0)
    for _ in range(2):
        d,idx=cKDTree(c).query(r+shift)
        good=d<4
        if good.sum()<20:raise ValueError('Too few matched stars')
        shift=np.median(c[idx[good]]-r[good],axis=0)
    back=cKDTree(r+shift).query(c)[1]
    good=(d<3)&(back[idx]==np.arange(len(r)))
    return np.flatnonzero(good),idx[good],shift


def summarize_ratio(values,cells):
    goodcells=[c for c in np.unique(cells) if np.sum(cells==c)>=10]
    if len(goodcells)<8:raise ValueError('Insufficient field coverage')
    groups=[values[cells==c] for c in goodcells]
    # Equal-area cells prevent a crowded corner dominating the decision.
    med=np.array([np.median(g) for g in groups])
    rng=np.random.default_rng(20261111)
    boot=np.median(np.stack([np.median(rng.choice(g,size=(500,len(g)),replace=True),axis=1) for g in groups]),axis=0)
    return dict(ratio=float(np.exp(np.median(med))),ci95=np.exp(np.percentile(boot,[2.5,97.5])).tolist(),
                cells=len(groups),cell_ratios={str(int(c)):float(np.exp(m)) for c,m in zip(goodcells,med)})


def compare(directory, final_through=None):
    directory=Path(directory)
    summaries=sorted([json.loads(p.read_text()) for p in directory.glob('*_summary.json')],key=lambda s:s['utc'])
    if not summaries: raise ValueError('No measured frames')
    if final_through is not None:
        seq=[int(s['name'].rsplit('_',1)[1]) for s in summaries]
        if max(seq)!=int(final_through) or sorted(seq)!=list(range(min(seq),int(final_through)+1)):
            raise ValueError('Final range must be complete and end at --final-through')
    ref=summaries[0]; rc=read_catalog(directory/(ref['name']+'_sources.csv'))
    results=[]; pairs=[]
    for s in summaries:
        cat=read_catalog(directory/(s['name']+'_sources.csv'))
        row=dict(name=s['name'],utc=s['utc'],n_sources=s['n_sources'],d50_arcsec=s['d50_arcsec'],d80_arcsec=s['d80_arcsec'],ellipticity=s['ellipticity_median'])
        try:
            a,b,shift=match(rc,cat)
            for key in ['d50','d80']:
                row[key+'_relative']=summarize_ratio(np.log(cat[key][b]/rc[key][a]),rc['cell'][a])
            row['n_matched']=len(a);row['shift_pixels']=shift.tolist()
            row['status']='measured';pairs.append((s,cat))
        except ValueError as e:row['status']=str(e)
        results.append(row)
    valid=[r for r in results if r['status']=='measured']
    valid.sort(key=lambda r:r['d50_relative']['ratio'])
    winner=valid[0] if valid else None
    # Direct best-vs-each comparison; avoids equating independent error bars.
    if winner:
        best=read_catalog(directory/(winner['name']+'_sources.csv'))
        for row in valid:
            cat=read_catalog(directory/(row['name']+'_sources.csv'))
            a,b,shift=match(best,cat)
            row['vs_best']={key:summarize_ratio(np.log(cat[key][b]/best[key][a]),best['cell'][a]) for key in ['d50','d80']}
            row['vs_best']['n_matched']=len(a)
            strong=(best['snr'][a]>80)&(cat['snr'][b]>80)
            try: row['vs_best']['bright_d50']=summarize_ratio(np.log(cat['d50'][b[strong]]/best['d50'][a[strong]]),best['cell'][a[strong]])
            except ValueError:row['vs_best']['bright_d50']=None
    report=dict(status='PROVISIONAL: acquisition completion not asserted',reference=ref['name'],best=winner['name'] if winner else None,ranking=valid,excluded=[r for r in results if r['status']!='measured'],
                caveat='Bootstrap intervals quantify source sampling only, not time-variable seeing. Focus cannot be isolated from seeing without repeated exposures.')
    if final_through is not None:
        if report['excluded']: raise ValueError('Cannot finalize with unranked frames')
        report['status']=f'FINAL for user-approved exposure range {min(seq):04d}-{int(final_through):04d}'
        report['final_through']=int(final_through)
    (directory/'comparison.json').write_text(json.dumps(report,indent=2))
    lines=['# r-band focus comparison','',report['status'],'','Smaller D50/D80 is better; dimensions in arcsec. D50 is half-light diameter, not FWHM.','',
           '| Rank | Image | Usable | Matched | D50 | D80 | Relative D50 |','|---|---|---:|---:|---:|---:|---:|']
    for i,r in enumerate(valid):lines.append(f"| {i+1} | {r['name']} | {r['n_sources']} | {r['n_matched']} | {r['d50_arcsec']:.3f} | {r['d80_arcsec']:.3f} | {r['d50_relative']['ratio']:.4f} |")
    lines+=['',report['caveat'],'','Ranking: matched-star size ratios, median in each of 16 spatial cells, then median across cells. All usable stars retained in per-frame catalogs.']
    (directory/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    ordered=sorted(valid,key=lambda r:r['utc'])
    for key,ax in zip(['d50','d80'],axes):
        vals=np.array([r[key+'_relative']['ratio'] for r in ordered]);ci=np.array([r[key+'_relative']['ci95'] for r in ordered])
        ax.errorbar(range(len(ordered)),vals,yerr=np.maximum(0,np.vstack((vals-ci[:,0],ci[:,1]-vals))),fmt='o-',capsize=3)
        ax.set_xticks(range(len(ordered)),[r['name'].split('_')[-1] for r in ordered]);ax.set(xlabel='Exposure sequence',ylabel=key.upper()+' relative to first frame',title='Same stars / equal-area spatial cells');ax.grid(alpha=.2)
    fig.suptitle('r-band focus comparison'+(' / final selected range' if final_through is not None else ' (provisional)'));fig.tight_layout();fig.savefig(directory/'comparison.png',dpi=160);plt.close(fig)
    print('\n'.join(lines))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory');p.add_argument('--final-through',type=int);a=p.parse_args();compare(a.directory,a.final_through)
