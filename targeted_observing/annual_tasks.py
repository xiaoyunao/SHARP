"""Explicit special-science sampling tasks, never physical-solution certificates."""
import json
import numpy as np
from .catalog import BASE


def load_tasks():
    return {t['mpc_number']:t for t in json.loads((BASE/'data/special_science_candidates.json').read_text())['targets']}


def status(task,state,quality):
    state=state or {};ts=np.asarray(state.get('times',[]));errors=np.asarray(state.get('errors',[]))
    # A tighter floor for low-amplitude work is a requirement, not a calibration claim.
    amp=task.get('amplitude_mag');limit=min(.03,amp/5) if amp else .03
    good=np.asarray(state.get('valid',[]),bool)&(errors<=limit)&(np.abs(state.get('zeropoints',[]))<=quality['max_zeropoint_mag'])
    times=ts[good];night=np.asarray(state.get('night_indices',[]))[good]
    periods=task['periods_h'];baseline=float(np.ptp(times)/86400) if len(times) else 0.
    covers=[len(set(np.floor((times/(p*3600)%1)*20).astype(int)))/20 for p in periods] if len(times) else [0.]*len(periods)
    longest=max(periods);required_baseline=max(quality['baseline_days'],3*longest/24)
    nights=len(set(night));ready=(nights>=14 and len(times)>=100 and baseline>=required_baseline and min(covers)>=.8)
    return {'photometric_complete':bool(state.get('quality',{}).get('quality_lightcurve',False)),
            'period_solution_complete':None,'science_case_complete':None,'sampling_ready':bool(ready),
            'good_task_points':len(times),'task_nights':nights,'task_baseline_days':baseline,
            'required_baseline_days':required_baseline,'min_phase_coverage':min(covers),
            'required_single_point_error_mag':limit,
            'completion_scope':'sampling readiness proxy only; no flux/model residuals or unique period fit'}


def weights(tasks,history,values,quality,night_index):
    result={}
    for k,t in tasks.items():
        if k not in values:continue
        state=history.get(k);s=status(t,state,quality)
        if s['sampling_ready']:continue
        past=state.get('night_indices',[]) if state else []
        # Sparse long-period tasks: do not take a dense block every consecutive night.
        spacing=3 if max(t['periods_h'])>100 else 1
        if past and night_index-max(past)<spacing:continue
        result[k]=values[k]*(2 if t['tier'].startswith('A:') else .5 if t['tier'].startswith('C:') else 1)*(1+min(s['task_nights']/14,1))
    return result


def window_bonus(task,state,when):
    """Novel phase and alias separation relative to observed epochs, not flux prediction."""
    if not state or not state['times']:return 1.
    times=np.asarray(state['times'])[np.asarray(state['valid'],bool)]
    if not len(times):return 1.
    periods=np.asarray(task['periods_h'])*3600
    phases=(times[:,None]/periods)%1
    bins=np.floor(phases*20).astype(int)
    now=np.floor(((when/periods)%1)*20).astype(int)
    novelty=np.mean([now[j] not in bins[:,j] for j in range(len(periods))])
    # Relative epoch removes arbitrary catalogue phase zeros. This proxy cannot
    # predict a particular minimum without a shape/phase model.
    delta=(when-times[0])/periods
    separation=max((abs(np.sin(np.pi*(delta[i]-delta[j]))) for i in range(len(periods)) for j in range(i)),default=0.)
    return 1+.3*novelty+.2*separation
