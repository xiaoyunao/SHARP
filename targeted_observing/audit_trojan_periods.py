"""Exploratory differential photometry and alias-aware period fits.

Uses independently extracted catalogs, never imports production code. Fits are
preliminary: sparse nights and detector-position systematics limit inference.
"""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks


def collect(root):
    base=Path(__file__).parent
    tracks={}
    for number in ['17365','4722','2357','51364','32513']:
        path=next(base.glob(f'cache/horizons_{number}_*.json'),None)
        if path is None:path=next(base.glob(f'cache/annual/horizons_{number}_*.json'))
        tracks[number]=json.loads(path.read_text())
    config=dict(nights=['20261111','20261112','20261113'],numbers=list(tracks),tracks=tracks)
    script='CONFIG = '+repr(config)+'\n'+(base/'collect_trojan_photometry.py').read_text()
    result=subprocess.run(['ssh','-p','20093','smtpipeline@www.xinglong-naoc.cn',
                           '/home/smtpipeline/Softwares/miniconda3/bin/python -'],
                          input=script,text=True,capture_output=True)
    root.mkdir(parents=True,exist_ok=True)
    (root/'collection_stderr.txt').write_text(result.stderr)
    result.check_returncode()
    json.loads(result.stdout)
    (root/'catalog_extract.json').write_text(result.stdout)
    (root/'input_manifest.json').write_text(json.dumps(dict(
        script_sha256=hashlib.sha256(script.encode()).hexdigest(),
        extract_sha256=hashlib.sha256(result.stdout.encode()).hexdigest()),indent=2))


def differential(rows, method):
    key, err = 'Mag_'+method, 'MagErr_'+method
    values=np.array([[r[key] if r is not None and r['Flag']==0 and
                      r['Flag_ISO_Num']==0 and 0<r[err]<.05 and r[key]<30
                      else np.nan for r in row['refs']] for row in rows])
    usable=np.mean(np.isfinite(values),axis=0)>.9
    offsets=np.zeros(len(rows))
    for _ in range(5):
        means=np.nanmedian(values-offsets[:,None],axis=0)
        residual=values-means
        offsets=np.nanmedian(residual[:,usable],axis=1)
        scatter=1.4826*np.nanmedian(abs(residual-offsets[:,None]),axis=0)
        usable=(np.mean(np.isfinite(values),axis=0)>.9)&(scatter<.03)
    if usable.sum()<15:
        raise ValueError(f'Only {usable.sum()} stable references')
    ref_res=values-means-offsets[:,None]
    # Independent reference halves diagnose coherent calibration instability.
    ii=np.flatnonzero(usable)
    half=np.nanmedian(residual[:,ii[::2]],axis=1)-np.nanmedian(residual[:,ii[1::2]],axis=1)
    zeroerr=1.253*np.nanstd(ref_res[:,usable],axis=1)/np.sqrt(np.sum(np.isfinite(values[:,usable]),axis=1))
    y=np.array([r['target'][key] for r in rows])-offsets
    error=np.sqrt(np.array([r['target'][err] for r in rows])**2+zeroerr**2+.01**2)
    valid=np.array([r['target']['Flag']==0 and r['target']['Flag_ISO_Num']==0 and
                    0<r['target'][err]<.10 and r['target'][key]<30 and r['exposure_s']>100
                    for r in rows])&np.isfinite(y)&(np.sum(np.isfinite(values[:,usable]),axis=1)>=15)
    info=dict(n_references=int(usable.sum()),median_reference_scatter=float(np.nanmedian(scatter[usable])),
              half_reference_difference_rms=float(np.nanstd(half)),
              max_abs_half_reference_difference=float(np.nanmax(abs(half))),
              reference_ids=ii.tolist())
    return y,error,valid,offsets,info


def fit(t,y,e,night,period,harmonics=2,night_offsets=False):
    cols=[np.ones(len(t))]
    if night_offsets:
        cols.extend((night==n).astype(float) for n in np.unique(night)[1:])
    angle=2*np.pi*t/period
    for k in range(1,harmonics+1):
        cols.extend([np.sin(k*angle),np.cos(k*angle)])
    design=np.column_stack(cols)
    coef=np.linalg.lstsq(design/e[:,None],y/e,rcond=None)[0]
    model=design@coef
    chi2=np.sum(((y-model)/e)**2)
    return float(chi2),model,coef


def search(t,y,e,night,harmonics=2,night_offsets=False):
    frequencies=np.linspace(1/100,1/2,10000)
    periods=1/frequencies
    chi=np.array([fit(t,y,e,night,p,harmonics,night_offsets)[0] for p in periods])
    peaks=find_peaks(-chi)[0]
    peaks=peaks[np.argsort(chi[peaks])]
    candidates=[]
    for i in peaks[:8]:
        opt=minimize_scalar(lambda p:fit(t,y,e,night,p,harmonics,night_offsets)[0],
                            bounds=(periods[i+1],periods[i-1]),method='bounded')
        candidates.append(dict(period_h=float(opt.x),chi2=float(opt.fun),delta_chi2=float(opt.fun-min(chi))))
    best=min(c['chi2'] for c in candidates)
    for i in (0,len(periods)-1):
        candidates.append(dict(period_h=float(periods[i]),chi2=float(chi[i]),boundary=True))
    best=min(c['chi2'] for c in candidates)
    for c in candidates:c['delta_chi2']=c['chi2']-best
    return sorted(candidates,key=lambda c:c['chi2']),periods,chi


def thymbraeus_checks(root):
    data=pd.read_csv(root/'differential_photometry.csv')
    data=data[data.number==17365]
    result=[]
    for method in ['PSF','Aper3','Aper4','Kron']:
        for free in ([False,True] if method=='PSF' else [False]):
            for omit in [None,20261111,20261112,20261113] if method=='PSF' and not free else [None]:
                d=data[(data.method==method)&data.good]
                if omit is not None:d=d[d.night!=omit]
                t,y,e,n=[d[k].to_numpy() for k in ['hours','mag','err','night']]
                # Local sensitivity check, not an independent global period search.
                grid=np.linspace(6.0,6.6,1201)
                chi=np.array([fit(t,y,e,n,p,4,free)[0] for p in grid])
                best=float(grid[np.argmin(chi)])
                _,m,c=fit(t,y,e,n,best,4,free)
                result.append(dict(method=method,night_offsets=free,omitted_night=omit,
                    local_photometric_period_h=best,double_period_h=2*best,
                    chi2=float(chi.min()),n=len(d),rms=float(np.std(y-m))))
    (root/'thymbraeus_sensitivity.json').write_text(json.dumps(result,indent=2))
    d=data[(data.method=='PSF')&data.good]
    t,y,e,n=[d[k].to_numpy() for k in ['hours','mag','err','night']]
    high={f'H{h}':search(t,y,e,n,h,False)[0] for h in [4,6]}
    (root/'higher_harmonics.json').write_text(json.dumps(high,indent=2))
    return result


def main(root):
    data=json.loads((root/'catalog_extract.json').read_text())
    all_rows=[];summaries={};grids={}
    literature={'17365':12.672,'2357':14.4028,'51364':30.4567,'4722':18.4573}
    for number in data['config']['numbers']:
        rows=sorted([r for r in data['measurements'] if r['number']==number],key=lambda r:r['jd'])
        jd=np.array([r['jd'] for r in rows]);night=np.array([r['night'] for r in rows])
        track=data['config']['tracks'][number]
        delta=np.interp(jd,track['datetime_jd'],track['delta'])
        radius=np.interp(jd,track['datetime_jd'],track['r'])
        # UTC emission times; absolute zero is arbitrary for this short-baseline fit.
        emission=jd-delta*499.004783836/86400
        time=(emission-emission.min())*24
        correction=5*np.log10(radius*delta);correction-=np.median(correction)
        summary=dict(number=number,n_detections=len(rows),nights={},methods={},search={},literature_h=literature.get(number))
        for n in np.unique(night):
            take=night==n
            summary['nights'][n]=dict(n=int(take.sum()),span_h=float(np.ptp(jd[take])*24),
                 first_jd=float(jd[take].min()),last_jd=float(jd[take].max()),
                 exposure_h=sum(r['exposure_s'] for r in rows if r['night']==n)/3600)
        for method in ['PSF','Aper3','Aper4','Kron']:
            y,e,valid,zp,info=differential(rows,method)
            y-=correction
            info.update(n_good=int(valid.sum()),median_error=float(np.median(e[valid])),
                        observed_p95_p05=float(np.ptp(np.percentile(y[valid],[5,95]))))
            summary['methods'][method]=info
            for i,r in enumerate(rows):
                all_rows.append(dict(number=number,method=method,night=r['night'],file=r['file'],
                    jd=r['jd'],emission_jd=emission[i],hours=time[i],mag=y[i],err=e[i],
                    good=bool(valid[i]),zero_offset=zp[i],sep_arcsec=r['sep_arcsec'],
                    x=r['target']['X_Win'],y=r['target']['Y_Win'],fwhm=r['target']['FWHM']))
            if len(np.unique(night))<2:continue
            for free in [False,True] if method=='PSF' else [False]:
                for harmonics in [1,2,3] if method=='PSF' else [2]:
                    key=f'{method}_H{harmonics}_'+('night_offsets' if free else 'shared')
                    cand,periods,chi=search(time[valid],y[valid],e[valid],night[valid],harmonics,free)
                    result=dict(candidates=cand,ndata=int(valid.sum()),nparameters=1+2*harmonics+(len(np.unique(night))-1 if free else 0))
                    if number in literature:
                        result['literature_chi2']=fit(time[valid],y[valid],e[valid],night[valid],literature[number],harmonics,free)[0]
                    summary['search'][key]=result
                    grids[number+'_'+key]=np.array([periods,chi])
        summaries[number]=summary
        print(number,json.dumps(dict(nights=summary['nights'],methods=summary['methods'],
              best={k:[round(c['period_h'],4) for c in v['candidates'][:3]] for k,v in summary['search'].items()})),flush=True)
    pd.DataFrame(all_rows).to_csv(root/'differential_photometry.csv',index=False)
    (root/'period_summary.json').write_text(json.dumps(summaries,indent=2))
    np.savez_compressed(root/'period_grids.npz',**grids)
    thymbraeus_checks(root)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=Path('targeted_observing/products/trojan_period_20261114'))
    parser.add_argument('--collect',action='store_true',help='Refresh inputs with read-only SSH extraction')
    args=parser.parse_args()
    if args.collect:collect(args.root)
    main(args.root)
