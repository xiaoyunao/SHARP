"""Independent overlap, swept-trajectory, legacy and TJ_test checks."""
import json
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from astropy.time import Time
from .catalog import BASE,interpolate
from .geometry import coord,offset,setup_iers,Sky
from .fixed_fields import FixedFields,unit_vectors,tan_xy,edge_distance,Grid
from .fixed_field_build import OUTPUT,HANDOFF,sha,write_json
from .annual_trojan import load_tracks


def edges_to_nearby(ra,dec,rows,width=2.944,height=2.944):
    """Vectorized complete regional neighbor search (cap wider than field diagonal)."""
    ar=np.array([r['ra_deg'] for r in rows]);dr=np.array([r['dec_deg'] for r in rows])
    tree=cKDTree(unit_vectors(ar,dr))
    # At these pitches and Dec<=60 the enclosing cap has at most nine centres.
    dist,ix=tree.query(unit_vectors(ra,dec),k=min(12,len(rows)))
    if ix.ndim==1:ix=ix[:,None];dist=dist[:,None]
    if np.any(dist[:,-1]<2*np.sin(np.deg2rad(2.2)/2)) and len(rows)>=12:
        raise ValueError('Neighbor cap not exhausted')
    x,y=tan_xy(np.asarray(ra)[:,None],np.asarray(dec)[:,None],ar[ix],dr[ix])
    return np.minimum(width/2-abs(x),height/2-abs(y)),ix


def run(output=OUTPUT):
    cat=FixedFields.read(output/'registry.json');grid=cat.grid
    data=np.load(output/'annual_samples.npz');ts=data['seconds'];ids=data['ids'];ra=data['ra'];dec=data['dec'];good=data['strict'];env=data['envelope']
    seed=pd.read_csv(HANDOFF/'tjn_fixed_grid_seed_v0.csv').to_dict('records')
    objects=[];sweep_count=0;sweep_min=90.;overlap_counts=defaultdict(int)
    for i,k in enumerate(ids):
        take=good[i]
        edge,_=edges_to_nearby(ra[i,take],dec[i,take],seed,width=9216*1.15/3600,height=9232*1.15/3600)
        seed_covered=(edge.max(axis=1)>=0)
        e,ix=edges_to_nearby(ra[i,take],dec[i,take],cat.fields)
        quality=e>=grid.margin+grid.guard
        if not np.all(quality.any(axis=1)):raise ValueError('Final overlap audit failed')
        indices=np.flatnonzero(take)
        for rr,cc in zip(*np.where(quality)):
            overlap_counts[(int(ix[rr,cc]),int(indices[rr]))]+=1
        objects.append({'object_id':k,'strict_samples':int(take.sum()),'seed_nominal_covered':int(seed_covered.sum()),
                        'seed_nominal_missed':int((~seed_covered).sum()),'final_quality_covered':int(quality.any(axis=1).sum())})
        # Check three interior points in every possibly observable half-hour segment.
        adjacent=(np.diff(ts)<1800.1)&(env[i,:-1]|env[i,1:])
        a=np.flatnonzero(adjacent)
        dra=(ra[i,a+1]-ra[i,a]+180)%360-180
        for f in (.25,.5,.75):
            ar=(ra[i,a]+dra*f)%360;dr=dec[i,a]+(dec[i,a+1]-dec[i,a])*f
            edge,_=edges_to_nearby(ar,dr,cat.fields)
            best=edge.max(axis=1)
            if len(best):sweep_min=min(sweep_min,float(best.min()))
            sweep_count+=len(best)
            if np.any(best<grid.margin):raise ValueError('Swept segment loses quality footprint')
        if (i+1)%300==0:print(f'Independent sweep: {i+1}/{len(ids)}',flush=True)
    ob=pd.DataFrame(objects);ob.to_csv(output/'seed_vs_final_objects.csv',index=False)
    oc=defaultdict(list)
    for (fi,t),n in overlap_counts.items():oc[fi].append(n)
    pd.DataFrame([{'field_id':f['field_id'],'active_samples':len(oc[i]),'median_simultaneous_quality_objects':float(np.median(oc[i])) if oc[i] else 0,
        'max_simultaneous_quality_objects':max(oc[i],default=0),'scope':'all overlapping quality footprints; counts within each sample are unique object IDs'} for i,f in enumerate(cat.fields)]).to_csv(output/'overlap_occupancy.csv',index=False)
    # Independent denser legacy footprint points, different lattice from build.
    old=pd.read_csv(HANDOFF/'inputs/tjn_B_unique_fields_with_cloud.csv');legacy_count=legacy_bad=0
    for r in old.itertuples():
        q=np.linspace(-1.5,1.5,61);x,y=np.meshgrid(q,q);p=offset(coord(r.ra_deg,r.dec_deg),x.ravel(),y.ravel())
        edge,_=edges_to_nearby(p.ra.deg,p.dec.deg,cat.fields)
        legacy_count+=len(edge);legacy_bad+=int(np.sum(edge.max(axis=1)<grid.margin))
    if legacy_bad:raise ValueError('Legacy footprint hole')
    # Geometric adjacency includes corner contacts and wrap; nominal overlap checked densely.
    adjacent=[];centres=unit_vectors([f['ra_deg'] for f in cat.fields],[f['dec_deg'] for f in cat.fields])
    for i,j in sorted(cKDTree(centres).query_pairs(2*np.sin(np.deg2rad(4.3)/2))):
        f,h=cat.fields[i],cat.fields[j]
        from .fixed_fields import boundary
        b=boundary(f['ra_deg'],f['dec_deg'],grid,samples=65)
        overlaps=bool(np.any(edge_distance(b[:,0],b[:,1],h['ra_deg'],h['dec_deg'],grid)>=0))
        if overlaps:adjacent.append({'field_id':f['field_id'],'neighbor_id':h['field_id'],'relation':'sampled nominal boundary overlap'})
    pd.DataFrame(adjacent).to_csv(output/'neighbors.csv',index=False)
    result={'swept_interior_points':sweep_count,'swept_minimum_best_edge_deg':sweep_min,
        'legacy_footprint_points':legacy_count,'legacy_quality_uncovered':legacy_bad,
        'strict_samples':int(ob.strict_samples.sum()),'seed_nominal_covered':int(ob.seed_nominal_covered.sum()),
        'seed_nominal_missed':int(ob.seed_nominal_missed.sum()),'final_quality_covered':int(ob.final_quality_covered.sum()),
        'seed_objects_with_any_nominal_miss':int((ob.seed_nominal_missed>0).sum()),
        'overlap_max_simultaneous_targets':max(overlap_counts.values()),
        'overlap_mean_simultaneous_targets':float(np.mean(list(overlap_counts.values()))),
        'neighbor_pairs':len(adjacent),'coverage_scope':'nominal and quality model; measured usable masks unknown'}
    write_json(output/'independent_audit.json',result)
    print(result,flush=True)
    tonight_and_special(output,cat)


def tonight_and_special(output,cat):
    tracks=load_tracks();config=json.loads((BASE/'defaults.json').read_text());setup_iers(BASE/'cache/iers_20261104.all');sky=Sky(config)
    source=BASE/'products/full_night_20261111/equatorial_audit/exposures.csv'
    rows=pd.read_csv(source);rows=rows[rows.field_id=='TJ_test']
    if len(rows)!=102:raise ValueError('Tonight list changed; re-read controller audit')
    fixed=cat.fields[cat.tree.query(unit_vectors([79.002262],[23.719287]))[1][0]]
    # Nearest spherical centre, not an assertion of simultaneous target retention.
    result=[];points=[]
    for label,times in [('tonight',rows.mid_unix.to_numpy()),('35days',np.arange(Time('2026-11-11T04:00:00').unix,Time('2026-12-16T04:00:00').unix,1800))]:
        for k in ['17365','4722']:
            ra,dec=interpolate(tracks[k],times);mag=np.interp(Time(times,format='unix').jd,tracks[k]['datetime_jd'],tracks[k]['V'])
            e=edge_distance(ra,dec,fixed['ra_deg'],fixed['dec_deg'],cat.grid)
            olde=edge_distance(ra,dec,79.002262,23.719287,cat.grid)
            memberships=cat.locate(ra,dec,quality=True)
            frame,sun,_,_=sky.environment(times);alt=coord(ra,dec).transform_to(frame).alt.deg
            visible=(sun.alt.deg<=-18)&(alt>=32.5)&(mag<=19.2)
            result.append({'period':label,'object_id':k,'new_nearby_field_id':fixed['field_id'],'samples':len(times),
                'new_same_field_quality_fraction':float(np.mean(e>=cat.grid.margin+cat.grid.guard)),
                'old_nominal_model_quality_fraction':float(np.mean(olde>=cat.grid.margin+cat.grid.guard)),
                'new_same_field_minimum_edge_deg':float(e.min()),'visible_samples':int(visible.sum()),
                'visible_same_field_quality_fraction':float(np.mean(e[visible]>=cat.grid.margin+cat.grid.guard)) if visible.any() else None,
                'any_fixed_field_quality_fraction':float(np.mean([bool(m) for m in memberships])),
                'all_fields_json':json.dumps(sorted({f for m in memberships for f in m}))})
            for t,r,d,v,edge,m,vis in zip(times,ra,dec,mag,e,memberships,visible):
                points.append({'period':label,'object_id':k,'time_utc':Time(t,format='unix').isot+'Z','ra_deg':r,'dec_deg':d,'V':v,
                    'new_same_field_edge_deg':edge,'quality_field_ids':';'.join(m),'visible':bool(vis)})
    pd.DataFrame(result).to_csv(output/'TJ_test_retention.csv',index=False)
    pd.DataFrame(points).to_csv(output/'TJ_test_tracks.csv.gz',index=False)
    write_json(output/'TJ_test_provenance.json',{'input_exposures_sha256':sha(source),
        'controller_sha256':sha(BASE/'products/full_night_20261111/20261111.txt'),
        'historical_identity':'TJ_test / AF_148ecdaae5','original_centre':[79.002262,23.719287],
        'nearby_new_field':{k:fixed[k] for k in ['field_id','ra_deg','dec_deg','grid_cell_id']},
        'comparison':'nominal TAN model; historical images keep original geometry; no executable file changes'})
    special=pd.read_csv(HANDOFF/'inputs/tjn_special_16_cloud_labels.csv',dtype={'object_id':str})
    objects=pd.read_csv(output/'objects.csv',dtype={'object_id':str})
    special=special.merge(objects,on='object_id',how='left',validate='one_to_one')
    special['coverage_status']=np.where(special.in_annual_1063,'see annual object audit','outside frozen annual parent; no annual track; weekly brightness/declination rejected')
    special.to_csv(output/'special_targets.csv',index=False)

if __name__=='__main__':run()
