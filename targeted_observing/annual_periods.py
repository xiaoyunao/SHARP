"""Injection recovery at realized timestamps, with blind and external-prior fits."""
import argparse
import json
from pathlib import Path
import numpy as np
from astropy.timeseries import LombScargle
from scipy.optimize import minimize_scalar
from .catalog import BASE,write_csv


def fit_frequency(t,y,error,frequency):
    angle=2*np.pi*frequency*t
    design=np.column_stack([np.ones(len(t)),np.cos(angle),np.sin(angle),np.cos(2*angle),np.sin(2*angle)])
    coefficients=np.linalg.lstsq(design/error[:,None],y/error,rcond=None)[0]
    residual=(y-design@coefficients)/error
    return float(residual@residual),coefficients


def search(t,y,error,prior_period=None):
    span=np.ptp(t)
    if len(t)<15 or span<=0:return None
    low,high=1/1000,1/2;df=1/(span*6)
    if prior_period is not None:low,high=1/(prior_period*1.05),1/(prior_period*.95)
    frequencies=np.linspace(low,high,max(32,int((high-low)/df)+1))
    power=LombScargle(t,y,error).power(frequencies,method='fast',assume_regular_frequency=True)
    maxima=np.flatnonzero((power[1:-1]>=power[:-2])&(power[1:-1]>=power[2:]))+1
    maxima=np.r_[0,maxima,len(power)-1];indices=maxima[np.argsort(power[maxima])[-12:]]
    candidates=[]
    for frequency in frequencies[indices]:
        for factor in (0.5,1.,2.):
            center=frequency*factor
            if not low<=center<=high:continue
            bounds=(max(low,center-df*1.5),min(high,center+df*1.5))
            result=minimize_scalar(lambda f:fit_frequency(t,y,error,f)[0],bounds=bounds,method='bounded',options={'xatol':1e-10})
            chi,coeff=fit_frequency(t,y,error,result.x)
            candidates.append((chi,float(result.x),coeff))
    candidates.sort(key=lambda r:r[0])
    if not candidates:return None
    best=candidates[0]
    independent=[r for r in candidates[1:] if abs(r[1]-best[1])>1/(span*2)]
    delta=independent[0][0]-best[0] if independent else None
    phase=np.linspace(0,2*np.pi,1000,endpoint=False);coef=best[2]
    model=coef[0]+coef[1]*np.cos(phase)+coef[2]*np.sin(phase)+coef[3]*np.cos(2*phase)+coef[4]*np.sin(2*phase)
    local_delta=min(fit_frequency(t,y,error,best[1]/factor)[0]-best[0] for factor in (.99,1.01))
    return {'period_h':1/best[1],'amplitude_mag':float(np.ptp(model)),
            'delta_chi2_alternative':delta,'delta_chi2_one_percent':float(local_delta),
            'unique_in_searched_candidates':delta is not None and delta>=9 and local_delta>=9,
            'search_scope':'12 strongest one-harmonic peaks plus half/double frequencies, local two-harmonic refinement'}


def injected(t,period,amplitude,phase,shape_name='asymmetric'):
    x=2*np.pi*t/period+phase
    grid=np.linspace(0,2*np.pi,2048,endpoint=False)
    shapes={
        'asymmetric':lambda z:.65*np.cos(2*z)+.35*np.cos(z+.3),
        'symmetric':lambda z:np.cos(2*z),
        'complex':lambda z:.55*np.cos(2*z)+.25*np.cos(z+.4)+.20*np.cos(4*z+.6)}
    shape=shapes[shape_name]
    return amplitude*shape(x)/np.ptp(shape(grid))


def run(root,max_objects=24,phases=3):
    samples=json.loads((root/'samples_realization0.json').read_text())
    eligible=[(k,s) for k,s in samples.items() if s['quality']['good_nights']>=5 and s['quality']['good_points']>=40]
    eligible.sort(key=lambda x:(x[1]['quality']['baseline_days'],x[1]['quality']['median_random_error_mag']))
    if len(eligible)>max_objects:eligible=[eligible[i] for i in np.linspace(0,len(eligible)-1,max_objects,dtype=int)]
    rng=np.random.default_rng(20261104);rows=[]
    for k,s in eligible:
        mask=np.array(s['valid'])&(np.array(s['errors'])<=.05)&(abs(np.array(s['zeropoints']))<=.02)
        t=np.array(s['times'])[mask];err=np.array(s['errors'])[mask];zp=np.array(s['zeropoints'])[mask]
        t=(t-t.min())/3600
        for period in [2.,4.,8.,12.,24.,48.,100.,200.,500.,1000.]:
            for amplitude in [.1,.3,.6]:
                for trial in range(phases*3):
                    shape_name=['asymmetric','symmetric','complex'][trial%3]
                    phase=rng.uniform(0,2*np.pi);signal=injected(t,period,amplitude,phase,shape_name)
                    observed=signal+zp+rng.normal(0,err)
                    for mode in ('blind','external_5percent_prior'):
                        result=search(t,observed,err,None if mode=='blind' else period)
                        recovered=bool(result and abs(result['period_h']/period-1)<.01)
                        rows.append({'object_id':k,'period_injected_h':period,'amplitude_injected_mag':amplitude,'phase_rad':phase,'shape':shape_name,
                            'mode':mode,'good_points':len(t),'baseline_h':float(np.ptp(t)),
                            'period_recovered':recovered,'period_unique_recovered':bool(recovered and result['unique_in_searched_candidates']),
                            'declared_unique':bool(result and result['unique_in_searched_candidates']),
                            'false_unique':bool(result and result['unique_in_searched_candidates'] and not recovered),
                            'amplitude_recovered':bool(result and abs(result['amplitude_mag']-amplitude)<=max(.05,.2*amplitude)),
                            'best_period_h':result['period_h'] if result else None,
                            'best_amplitude_mag':result['amplitude_mag'] if result else None,
                            'delta_chi2_alternative':result['delta_chi2_alternative'] if result else None})
        print(f'{root.name}: injection object {k}',flush=True)
    write_csv(root/'period_injection.csv',rows)
    summary={'objects':len(eligible),'trials':len(rows),'period_range_h':[2,1000],
        'method':'asymmetric, symmetric double-peaked, and four-harmonic shapes; random phases; actual per-point errors plus correlated night zeropoints; no nightly demeaning',
        'uniqueness':'delta chi-square >=9 to separated candidate frequency AND refitted periods displaced by +/-1%; heuristic, not calibrated posterior probability',
        'scope':'representative realized sampling from weather realization 0; not annual reliable-period yield for all 200 realizations',
        'prior':'external ±5% information around true injected period; hypothetical known-period case, never claimed as new period',
        'limitations':'ideal corrected mean magnitude and constant rotation; no shape evolution, timing-system error, or orbital-model error'}
    (root/'period_injection_manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--objects',type=int,default=24);p.add_argument('--phases',type=int,default=3)
    args=p.parse_args();print(run(args.directory,args.objects,args.phases),flush=True)


if __name__=='__main__':main()
