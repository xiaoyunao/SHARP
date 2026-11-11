#!/usr/bin/env python3
"""Compare observed PSFs on identical stars across passbands; not a focus-offset fit."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from compare_focus import match, summarize_ratio


def compare_passbands(paths, out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    summaries=[json.loads(Path(p).read_text()) for p in paths]
    bands=[s['band'] for s in summaries]
    if len(bands)!=len(set(bands)) or 'unconfirmed' in bands:
        raise ValueError('Require one explicitly assigned frame per distinct passband')
    cats=[np.genfromtxt(Path(p).parent/(s['name']+'_sources.csv'),delimiter=',',names=True) for p,s in zip(paths,summaries)]
    ref=cats[0];maps=[];shifts=[];keep=np.ones(len(ref),bool)
    for cat in cats:
        a,b,shift=match(ref,cat);m=np.full(len(ref),-1,dtype=int);m[a]=b
        keep &=m>=0;maps.append(m);shifts.append(shift.tolist())
    ids=np.flatnonzero(keep);cells=ref['cell'][ids]
    bright=np.ones(len(ids),bool)
    for cat,m in zip(cats,maps):
        rows=cat[m[ids]];bright &=(rows['snr']>100)&(rows['peak']<20000)
    rows=[]
    for p,s,cat,m,shift in zip(paths,summaries,cats,maps,shifts):
        ix=m[ids]
        r=dict(band=s['band'],name=s['name'],exposure=s['exposure'],utc=s['utc'],
               n_sources=s['n_sources'],d50_arcsec=s['d50_arcsec'],d80_arcsec=s['d80_arcsec'],
               ellipticity=s['ellipticity_median'],shift_pixels=shift,
               summary_sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest(),
               catalog_sha256=hashlib.sha256((Path(p).parent/(s['name']+'_sources.csv')).read_bytes()).hexdigest())
        for key in ['d50','d80']:
            r[key+'_relative']=summarize_ratio(np.log(cat[key][ix]/ref[key][ids]),cells)
        try:r['bright_d50_relative']=summarize_ratio(np.log(cat['d50'][ix[bright]]/ref['d50'][ids[bright]]),cells[bright])
        except ValueError as e:r['bright_d50_relative']={'unavailable':str(e)}
        r['matched_cells']={str(c):dict(n=int(np.sum(cells==c)),d50_arcsec=float(np.median(cat['d50'][ix[cells==c]])*s['pixel_scale']),ellipticity=float(np.median(cat['ellipticity'][ix[cells==c]]))) for c in range(16)}
        rows.append(r)
    report=dict(reference=bands[0],n_common=len(ids),n_common_bright=int(bright.sum()),frames=rows,
                caveat='Observed image quality at one nominal setting. Exposure times, epochs, passbands and atmospheric conditions differ; these are not independently measured best-focus positions.')
    (out/'passband_comparison.json').write_text(json.dumps(report,indent=2))
    lines=['# Nominal 5.88: observed passband image quality','',report['caveat'],'',f'Identical stars in all frames: {len(ids)}; bright unsaturated subset: {bright.sum()}.','',
           '| Band | Image | Exposure (s) | Usable | D50 (arcsec) | D80 (arcsec) | Same-star D50 / reference | Same-star D80 / reference |',
           '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:lines.append(f"| {r['band']} | {r['name']} | {r['exposure']:.3f} | {r['n_sources']} | {r['d50_arcsec']:.3f} | {r['d80_arcsec']:.3f} | {r['d50_relative']['ratio']:.4f} | {r['d80_relative']['ratio']:.4f} |")
    lines+=['','Absolute widths above use all usable stars per image; ratios use only the identical stars, with equal weight per spatial cell. Bootstrap intervals in JSON cover source sampling, not seeing or other systematics.']
    (out/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,len(rows),figsize=(3.5*len(rows),4),constrained_layout=True)
    vals=np.array([[r['matched_cells'][str(c)]['d50_arcsec'] for c in range(16)] for r in rows])
    for ax,r,v in zip(np.atleast_1d(axes),rows,vals):
        sc=ax.imshow(v.reshape(4,4),origin='lower',vmin=float(vals.min()),vmax=float(vals.max()),cmap='viridis',interpolation='nearest')
        for c,value in enumerate(v):ax.text(c%4,c//4,f'{value:.2f}',ha='center',va='center',fontsize=10,color='black' if value>np.mean([vals.min(),vals.max()]) else 'white')
        ax.set(title=f"{r['band'].upper()} / {r['name'].rsplit('_',1)[1]} / {r['exposure']:.0f}s",xlabel='Detector x quarter',ylabel='Detector y quarter',xticks=range(4),yticks=range(4))
    fig.colorbar(sc,ax=np.atleast_1d(axes).tolist(),label='Median D50 [arcsec]',shrink=.8)
    fig.suptitle(f'Nominal 5.88 / same {len(ids)} stars in all passbands / no best-focus inference')
    fig.savefig(out/'passband_spatial_comparison.png',dpi=160);plt.close(fig)
    print('\n'.join(lines));return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('summaries',nargs='+');p.add_argument('--out',required=True);a=p.parse_args();compare_passbands(a.summaries,a.out)
