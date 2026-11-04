"""Illustrative period yield across all saved weather-realized quality curves."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .catalog import BASE,write_csv
from .science import load_periods,assess
from .annual_periods import injected,search


def phase_diagnostics(hours,period):
    phases=np.sort(np.mod(np.asarray(hours)/period,1.))
    return {'phase_bin_coverage_20':len(np.unique(np.floor(phases*20).astype(int)))/20,
            'largest_phase_gap':float(np.diff(np.r_[phases,phases[0]+1]).max())}


def object_assumption(k,evidence):
    seed=int(hashlib.sha256(('period-population:'+k).encode()).hexdigest()[:8],16)
    rng=np.random.default_rng(seed);entry=evidence.get(k,{})
    manual=entry.get('manual_review',{});model=entry.get('models2026');lc=entry.get('LCDB') or {}
    prior=None;source='illustrative log-uniform 2–1000 h, not measured population'
    if manual.get('measurements'):
        prior=manual['measurements'][-1]['period_h'];source=manual['citation']
    elif model and model.get('period_h'):
        prior=model['period_h'];source='2026 model table'
    elif lc.get('U','').startswith('3') and lc.get('period_h'):
        prior=lc['period_h'];source='frozen LCDB U>=3-class entry'
    period=float(prior or np.exp(rng.uniform(np.log(2),np.log(1000))))
    return {'period_h':period,'amplitude_mag':float(rng.choice([.1,.3,.6])),
            'phase_rad':float(rng.uniform(0,2*np.pi)),'external_prior_available':prior is not None,
            'shape':str(rng.choice(['asymmetric','symmetric','complex'],p=[.6,.2,.2])),
            'period_source':source,'role':assess(entry)['role']}


def run(root):
    summary=json.loads((root/'summary.json').read_text())
    if not summary.get('qualified_curves_saved'):raise ValueError('Re-run annual_simulation --save-curves')
    curves=dict(np.load(root/'qualified_curves.npz'));evidence=load_periods(BASE/'cache')
    assumptions={k:object_assumption(k,evidence) for k in set(curves['object_id'].tolist())}
    results=[];totals=np.zeros((summary['realizations'],3),int)
    for i,(k,realization) in enumerate(zip(curves['object_id'],curves['realization'])):
        a,b=curves['offsets'][i:i+2];t=curves['emission_unix'][a:b];t=(t-t.min())/3600
        error=curves['error'][a:b];zp=curves['zeropoint'][a:b];assumption=assumptions[str(k)];p=assumption['period_h']
        seed=int(hashlib.sha256(f'period-noise:{k}:{realization}'.encode()).hexdigest()[:8],16)
        rng=np.random.default_rng(seed);y=injected(t,p,assumption['amplitude_mag'],assumption['phase_rad'],assumption['shape'])+zp+rng.normal(0,error)
        for mode in ('blind','external_prior'):
            if mode=='external_prior' and not assumption['external_prior_available']:continue
            fit=search(t,y,error,None if mode=='blind' else p) if 2<=p<=1000 else None
            recovered=bool(fit and abs(fit['period_h']/p-1)<.01)
            unique=bool(recovered and fit['unique_in_searched_candidates'])
            amp=bool(fit and abs(fit['amplitude_mag']-assumption['amplitude_mag'])<=max(.05,.2*assumption['amplitude_mag']))
            totals[realization,0 if mode=='blind' else 1]+=unique
            if mode=='blind':totals[realization,2]+=amp
            results.append({'realization':int(realization),'object_id':str(k),'mode':mode,**assumption,**phase_diagnostics(t,p),
                'period_recovered':recovered,'period_unique_recovered':unique,'amplitude_recovered':amp,
                'declared_unique':bool(fit and fit['unique_in_searched_candidates']),
                'false_unique':bool(fit and fit['unique_in_searched_candidates'] and not recovered),
                'best_period_h':fit['period_h'] if fit else None,'delta_chi2_alternative':fit['delta_chi2_alternative'] if fit else None})
        if (i+1)%500==0:print(f'{root.name}: {i+1}/{len(curves["object_id"])} realized period curves',flush=True)
    write_csv(root/'period_yield.csv',results)
    (root/'period_population_assumptions.json').write_text(json.dumps(assumptions,indent=2)+'\n')
    out={'realizations':summary['realizations'],'quality_curves_tested':len(curves['object_id']),
        'blind_reliable_periods_P10_P50_P90':np.percentile(totals[:,0],[10,50,90]).tolist(),
        'known_prior_reliable_periods_P10_P50_P90':np.percentile(totals[:,1],[10,50,90]).tolist(),
        'amplitudes_P10_P50_P90':np.percentile(totals[:,2],[10,50,90]).tolist(),
        'counts':totals.tolist(),'lighttime':'observer times minus Horizons delta * 499.004783836 s/AU',
        'interpretation':'conditional illustrative injection forecast among quality lightcurves; neither observations nor guaranteed new periods',
        'population':'known catalogue/model periods where available; otherwise fixed log-uniform 2–1000h per object; amplitudes .1/.3/.6mag; illustrative shape mixture asymmetric/symmetric/complex=.6/.2/.2',
        'independence':'external-prior and blind counts overlap and must not be added',
        'uniqueness':'heuristic delta chi-square >=9 to separated candidates AND refitted +/-1% period alternatives; not a calibrated posterior',
        'limitations':'ideal phase/distance correction; fixed shape and spin; limited candidate-frequency search; 200 weather draws do not describe unknown population-prior uncertainty'}
    (root/'period_yield_summary.json').write_text(json.dumps(out,indent=2)+'\n');return out


def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);args=p.parse_args();print(run(args.directory),flush=True)


if __name__=='__main__':main()
