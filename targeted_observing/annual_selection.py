"""Compare start times using distinct-object greedy combinations.

The greedy set cover is an approximation, not a global scheduling optimum.
Every legal field participates; no unproven top-field pruning is applied.
"""
import json
import numpy as np
from scipy.sparse import csr_matrix
from .catalog import BASE

DEFAULT_QUALITY=json.loads((BASE/'annual_defaults.json').read_text())['quality']


def object_values(fields,values,history,night_index,completion,quality):
    result={}
    for k in {k for f in fields for k in f['core_targets']}:
        state=history.get(k);q=state.get('quality',{}) if state else {}
        if q.get('quality_lightcurve'):
            result[k]=0.;continue
        if completion is None:
            result[k]=values.get(k,0)*(1.7 if state and state['times'] else 1.);continue
        remaining,last,expected_points=completion
        first=q.get('first_good_night')
        if first is None:first=state['night_indices'][0] if state and state['times'] else night_index
        if last.get(k,-1)-first<quality['baseline_days']:
            result[k]=0.;continue
        need=max(1.,quality['nights']-q.get('good_nights',0),
                 (quality['points']-q.get('good_points',0))/max(expected_points,1e-9))
        continuity=1+min(1.,q.get('good_nights',0)/quality['nights'])
        if q.get('good_points',0)>=quality['points'] and q.get('good_nights',0)>=quality['nights'] and night_index-first<quality['baseline_days']:
            continuity*=.05
        result[k]=values.get(k,0)*min(1.,remaining.get(k,0)/need)*continuity
    return result


def combinations(menu,data,weights,count):
    """Vectorized sparse greedy union, independently for each available start."""
    fields=menu['fields'];valid=np.asarray(data['valid'],bool)
    if not valid.size or not valid.any():return [],np.zeros(valid.shape[1])
    if '_selection_membership' not in data:
        ids=sorted({k for f in fields for k in f['core_targets']});index={k:i for i,k in enumerate(ids)}
        rows=[];cols=[]
        for i,f in enumerate(fields):
            for k in set(f['core_targets']):rows.append(i);cols.append(index[k])
        matrix=csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(fields),len(ids)))
        data['_selection_membership']=(ids,matrix)
    ids,matrix=data['_selection_membership'];w=np.array([max(0.,weights.get(k,0.)) for k in ids])
    gains=np.asarray(matrix@w)[:,None]*np.ones((1,valid.shape[1]));gains[~valid]=-np.inf
    used=np.zeros((len(ids),valid.shape[1]),bool);chosen=[[] for _ in range(valid.shape[1])]
    for _ in range(min(count,len(fields))):
        winners=np.argmax(gains,axis=0);active=np.flatnonzero(gains[winners,np.arange(valid.shape[1])]>1e-12)
        if not len(active):break
        newly=np.zeros_like(used)
        for slot in active:
            j=winners[slot];chosen[slot].append(int(j))
            members=matrix.indices[matrix.indptr[j]:matrix.indptr[j+1]]
            newly[members,slot]=~used[members,slot]
        used|=newly
        gains-=matrix@(newly*w[:,None])
        gains[winners[active],active]=-np.inf
    return chosen,w@used


def pick(menu,data,values,count,history,night_index,completion=None,quality=None):
    quality=DEFAULT_QUALITY if quality is None else quality
    weights=object_values(menu['fields'],values,history,night_index,completion,quality)
    chosen,gains=combinations(menu,data,weights,count)
    if not chosen or not np.any(gains>0):return None
    order=np.argsort(-gains,kind='stable')
    near=order[gains[order]>=gains[order[0]]*.98]
    # Rotate only after all starts have comparable distinct-object scores.
    slot=int(np.roll(near,night_index%3*3)[0])
    return slot,chosen[slot]


def pick_route(menu,data,values,count,history,night_index,completion,quality,config,dwell,
               delay=0.,extra_readout=0.,extra_slew=0.,special=None,tasks=None):
    """Re-score each start's greedy union with its best exact small cyclic route.

    All <=4! route permutations are considered. Field-set search remains greedy;
    this is not a proof of an optimal joint field/route solution.
    """
    from itertools import permutations
    from .annual_simulation import numeric_visits,transition
    from .annual_tasks import window_bonus
    weights=special if special is not None else object_values(menu['fields'],values,history,night_index,completion,quality)
    chosen,_=combinations(menu,data,weights,count)
    if not chosen:return None
    cache=data.setdefault('_route_orders',{})
    candidates=[];evaluated={};reference=max(quality['points_per_night'],quality['points']/quality['nights'])
    for slot,indices in enumerate(chosen):
        if not indices:continue
        key=tuple(sorted(indices))
        if key not in cache:
            routes=[]
            for order in permutations(key):
                fs=[menu['fields'][j] for j in order]
                path=sum(transition(fs[j-1],fs[j],config) for j in range(1,len(fs)))
                closure=transition(fs[-1],fs[0],config)
                routes.append((order,path,closure))
            cache[key]=routes
        if key not in evaluated:
            n=len(key);frames=n*dwell;end=config['trojan']['block_s']-config['trojan']['handover_s']
            frame_time=frames*(config['trojan']['exposure_s']+config['readout_s']+extra_readout)
            def route_cost(route):
                order,path,closure=route
                first=delay+config['far_slew_s']+extra_slew+config['filter_change_s']+frame_time+path+(n-1)*extra_slew
                repeat=frame_time+path+closure+(n if n>1 else 0)*extra_slew
                cycles=0 if first>end else 1+int((end-first)/repeat)
                return (-cycles,first+max(0,cycles-1)*repeat)
            order,_,_=min(cache[key],key=route_cost)
            fs=[menu['fields'][j] for j in order]
            rows=numeric_visits(fs,dwell,config,delay,extra_readout,extra_slew)
            credits={}
            for j,f in enumerate(fs):
                rr=rows[rows[:,0]==j]
                if len(rr)<quality['points_per_night'] or np.ptp(rr[:,1])<quality['night_span_s']:continue
                credit=min(1.,len(rr)/reference)
                for k in f['core_targets']:credits[k]=max(credits.get(k,0),credit)
            evaluated[key]=(order,credits)
        order,credits=evaluated[key]
        gain=sum(weights.get(k,0)*v for k,v in credits.items())
        if special is not None and gain>0:
            when=float(data['starts'][slot])+config['trojan']['block_s']/2
            gain=sum(weights.get(k,0)*v*window_bonus(tasks[k],history.get(k),when) for k,v in credits.items() if k in weights)
        if gain>0:candidates.append((gain,slot,list(order)))
    if not candidates:return None
    best=max(v[0] for v in candidates)
    near=[v for v in candidates if v[0]>=best*.98]
    _,slot,order=near[(-night_index%3*3)%len(near)]
    return slot,order
