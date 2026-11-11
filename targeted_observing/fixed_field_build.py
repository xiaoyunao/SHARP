"""Reproducible annual fixed-field construction, using only frozen local inputs."""
import argparse
from dataclasses import asdict
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, get_sun
from astropy.table import Table
from astropy.time import Time
from .catalog import BASE
from .geometry import setup_iers, site, coord, offset
from .fixed_fields import Grid, registry, edge_distance, cell_boundary, boundary, REGISTRY_PATH

HANDOFF=BASE/'data/tjn_handoff_20261111'
OUTPUT=BASE/'products/tjn_fixed_20261111'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def load_inputs():
    candidates=pd.read_csv(HANDOFF/'inputs/tjn_annual_candidates_1063.csv',dtype={'object_id':str})
    index=json.loads((BASE/'products/annual/ephemeris_index.json').read_text())
    qa=json.loads((BASE/'products/annual/input_verification.json').read_text())
    manifest=json.loads((HANDOFF/'MANIFEST_SHA256.json').read_text())
    for name,digest in manifest.items():
        if sha(HANDOFF/name)!=digest:raise ValueError('Handoff hash mismatch: '+name)
    if index['failures'] or set(index['paths'])!=set(candidates.object_id):raise ValueError('Incomplete annual sample')
    config=json.loads((BASE/'defaults.json').read_text());tracks={};hashes={}
    for k,p in index['paths'].items():
        path=Path(p)
        if not path.exists():path=BASE/'cache/annual'/path.name
        t=json.loads(path.read_text());hashes[k]=sha(path)
        if hashes[k]!=qa['ephemeris_sha256'][k]:raise ValueError('Ephemeris changed: '+k)
        if t['pdes']!=k or t['site']!=config['site'] or 'UTC' not in t['source'] or 'ICRF' not in t['source']:
            raise ValueError('Ephemeris identity/site/frame mismatch')
        jd=np.array(t['datetime_jd'])
        if jd[0]>Time('2026-11-01').jd or jd[-1]<Time('2027-11-02').jd or not np.allclose(np.diff(jd),.25):
            raise ValueError('Missing six-hour samples')
        if any(len(t[x])!=len(jd) or not np.isfinite(t[x]).all() for x in ['RA','DEC','V','RA_rate','DEC_rate']):
            raise ValueError('Invalid track values')
        tracks[k]=t
    return candidates,tracks,config,{'ephemeris_sha256':hashes,'ephemeris_qa':{k:v for k,v in qa.items() if k!='ephemeris_sha256'},
        'handoff_manifest_sha256':sha(HANDOFF/'MANIFEST_SHA256.json'),
        'candidate_sha256':sha(HANDOFF/'inputs/tjn_annual_candidates_1063.csv'),
        'config_sha256':sha(BASE/'defaults.json'),'annual_index_sha256':sha(BASE/'products/annual/ephemeris_index.json')}


def best_grid_edge(ra,dec,grid):
    j,_=grid.assign(ra,dec);best=np.full(np.shape(ra),-np.inf)
    for dj in (-1,0,1):
        d,n,step=grid.row(j+dj);k=np.floor((ra%360)/step+.5).astype(int)%n
        for dk in (-1,0,1):best=np.maximum(best,edge_distance(ra,dec,((k+dk)%n)*step,d,grid))
    return best


def geometry_rows(grid,js):
    result=[]
    for j in sorted(set(js)):
        a,d=cell_boundary(j,grid);e=edge_distance(a,d,0,j*grid.sy,grid)
        result.append({'grid':grid.version,'row_index':int(j),'dec_deg':float(j*grid.sy),
                       'boundary_samples':len(a),'minimum_edge_deg':float(e.min()),
                       'nominal_pass':bool(e.min()>=0),'quality_guard_pass':bool(e.min()>=grid.margin+grid.guard)})
    return result


def prepare(output):
    output.mkdir(parents=True,exist_ok=True)
    candidates,tracks,config,provenance=load_inputs()
    provenance['iers']=setup_iers(BASE/'cache/iers_20261104.all')
    # UTC+8 observing nights run local noon to noon, not UTC calendar dates.
    start=Time('2026-11-01T04:00:00').unix;stop=Time('2027-11-01T04:00:00').unix
    all_seconds=np.arange(start,stop,1800.)
    frame=AltAz(obstime=Time(all_seconds,format='unix'),location=site(config))
    sun=get_sun(frame.obstime).transform_to(frame).alt.deg
    time_guard=.0043*900+.01
    keep=sun<=config['trojan']['sun_alt_deg']+time_guard
    seconds=all_seconds[keep];sun=sun[keep];frame=frame[keep];jd=frame.obstime.jd
    ids=candidates.object_id.tolist();ras=[];decs=[];mags=[];strict=[];envelopes=[];max_motion=0
    for n,k in enumerate(ids):
        t=tracks[k];ra=np.rad2deg(np.interp(jd,t['datetime_jd'],np.unwrap(np.deg2rad(t['RA']))))%360
        dec=np.interp(jd,t['datetime_jd'],t['DEC']);v=np.interp(jd,t['datetime_jd'],t['V'])
        alt=coord(ra,dec).transform_to(frame).alt.deg
        dv=float(np.max(abs(np.diff(t['V'])))) / 12  # half of a 30-minute interval / six hours
        good=(sun<=-18)&(alt>=config['center_alt_deg'])&(v<=19.2)
        possible=(alt>=config['center_alt_deg']-time_guard)&(v<=19.2+dv)
        ras.append(ra);decs.append(dec);mags.append(v);strict.append(good);envelopes.append(possible)
        max_motion=max(max_motion,float(np.hypot(t['RA_rate'],t['DEC_rate']).max())*.25/3600)
        if (n+1)%200==0:print(f'Annual tracks: {n+1}/{len(ids)}',flush=True)
    interpolation_guard=3/3600
    if max_motion+interpolation_guard>.01:raise ValueError('Increase grid trajectory guard for observed speed')
    provenance['sampling']={'step_s':1800,'visibility_time_guard_deg':time_guard,
        'max_half_step_motion_deg':max_motion,'interpolation_guard_deg':interpolation_guard,
        'original_qa_max_arcsec':provenance['ephemeris_qa']['max_position_residual_arcsec'],
        'selection':'all possible dark/altitude/brightness segments; Moon never permanently clips registry',
        'strict':'sampled Sun<=-18, object alt>=32.5, V<=19.2; not full two-hour block feasibility',
        'quality_core':'V<=18.5 proxy; not calibrated r photometric certification',
        'continuous_coverage':'visibility upper envelope plus swept-position guard; sampled statistics are not exact durations'}
    np.savez_compressed(output/'annual_samples.npz',ids=np.array(ids),seconds=seconds,ra=ras,dec=decs,mag=mags,strict=strict,envelope=envelopes)
    write_json(output/'provenance.json',provenance)
    print('Annual sample cache complete',flush=True)


def build(output,grid=Grid(),previous=None,release='20261111-rc1'):
    data=np.load(output/'annual_samples.npz');ids=data['ids'];ra=data['ra'];dec=data['dec'];mag=data['mag'];good=data['strict'];env=data['envelope'];ts=data['seconds']
    candidates=pd.read_csv(HANDOFF/'inputs/tjn_annual_candidates_1063.csv',dtype={'object_id':str}).set_index('object_id')
    old=pd.read_csv(HANDOFF/'inputs/tjn_B_unique_fields_with_cloud.csv')
    cells={};reasons={};legacy=[]
    def add(j,k,cloud,reason):
        key=(int(j),int(k));cells.setdefault(key,cloud);reasons.setdefault(key,set()).add(reason)
    for i,k in enumerate(ids):
        js,ks=grid.assign(ra[i,env[i]],dec[i,env[i]])
        for j,c in set(zip(js,ks)):add(j,c,candidates.loc[k,'cloud_proxy'],'annual_envelope')
    # Preserve all legacy footprint regions, with a separate traceable selection reason.
    for row in old.itertuples():
        q=np.linspace(-1.5,1.5,51);x,y=np.meshgrid(q,q)
        points=offset(coord(row.ra_deg,row.dec_deg),x.ravel(),y.ravel())
        js,ks=grid.assign(points.ra.deg,points.dec.deg)
        for j,k in sorted(set(zip(js,ks))):
            add(j,k,row.cloud_proxy,'legacy_B_envelope')
            legacy.append({'legacy_field_id':row.field_id,'grid_cell_id':grid.key(j,k)})
    doc=registry(cells,grid,previous=previous,release=release);by_key={f['grid_cell_id']:f for f in doc['fields']};lookup={}
    for f in doc['fields']:
        key=(f['row_index'],f['column_index']);f['selection_basis']=sorted(reasons.get(key,{'retained_from_previous_release'}));lookup[key]=f
    rows_audit=geometry_rows(grid,[f['row_index'] for f in doc['fields']])
    if not all(r['quality_guard_pass'] for r in rows_audit):raise ValueError('Ownership cells fail swept quality coverage')
    proposal=Grid(sx=2.67264,sy=2.67496,version='MP-pitch-proposal')
    proposal_audit=geometry_rows(proposal,range(-11,25))
    pd.DataFrame(rows_audit+proposal_audit).to_csv(output/'row_geometry_audit.csv',index=False)
    write_json(output/'registry.json',doc)
    flat=[{**{k:v for k,v in f.items() if not isinstance(v,(list,dict))},
           'selection_basis':';'.join(f['selection_basis']),
           'nominal_boundary_ref':'nominal_boundaries.json#'+f['field_id']} for f in doc['fields']]
    df=pd.DataFrame(flat);df.to_csv(output/'registry.csv',index=False,float_format='%.12f')
    fitsdf=df.fillna('UNMEASURED');Table.from_pandas(fitsdf).write(output/'registry.fits',overwrite=True)
    write_json(output/'nominal_boundaries.json',{f['field_id']:boundary(f['ra_deg'],f['dec_deg'],grid).tolist() for f in doc['fields']})
    for cloud in ('L4','L5'):df[df.planning_cloud==cloud].to_csv(output/f'TJ{cloud}.csv',index=False,float_format='%.12f')
    for r in legacy:r['field_id']=by_key[r['grid_cell_id']]['field_id']
    pd.DataFrame(legacy).to_csv(output/'legacy_to_fixed.csv',index=False)
    daily=[];objects=[];intervals=[];field_members={};counts={}
    iso_times=Time(ts,format='unix').isot
    dates=np.array([str(np.datetime64('2026-11-01')+np.timedelta64(int((t-Time('2026-11-01T04:00:00').unix)//86400),'D')) for t in ts])
    # Point table is compressed; interval and daily tables are convenient inspectable summaries.
    with gzip.open(output/'object_samples.csv.gz','wt') as stream:
        stream.write('object_id,time_utc,object_cloud,field_id,ra_deg,dec_deg,V,in_nominal,in_quality_model,in_usable,edge_distance_deg,target_role,strict_visible\n')
        for i,k in enumerate(ids):
            inds=np.flatnonzero(env[i]);js,ks=grid.assign(ra[i,inds],dec[i,inds]);cloud=candidates.loc[k,'cloud_proxy']
            fs=[lookup[(int(j),int(c))] for j,c in zip(js,ks)];fids=np.array([f['field_id'] for f in fs]);ar=np.array([f['ra_deg'] for f in fs]);dr=np.array([f['dec_deg'] for f in fs])
            edges=edge_distance(ra[i,inds],dec[i,inds],ar,dr,grid)
            if np.any(edges<grid.margin+grid.guard):raise ValueError('Trajectory falls outside quality envelope')
            oldedge=best_grid_edge(ra[i,inds],dec[i,inds],proposal)
            sampled=good[i,inds];core=sampled&(mag[i,inds]<=18.5)
            for idx,fid,e in zip(inds,fids,edges):
                role='core_proxy' if mag[i,idx]<=18.5 else 'secondary'
                stream.write(f'{k},{iso_times[idx]}Z,{cloud},{fid},{ra[i,idx]:.9f},{dec[i,idx]:.9f},{mag[i,idx]:.5f},True,True,,{e:.9f},{role},{bool(good[i,idx])}\n')
                if good[i,idx]:
                    field_members.setdefault(fid,set()).add(k)
                    counts[(fid,int(idx))]=counts.get((fid,int(idx)),0)+1
            visible_days=[];core_days=[];primary=[]
            for day in np.unique(dates[inds]):
                take=dates[inds]==day;s=sampled&take;c=core&take
                values=fids[s]
                if len(values):
                    visible_days.append(day);primary.append(str(values[len(values)//2]))
                if c.any():core_days.append(day)
                daily.append({'object_id':k,'date':day,'cloud':cloud,'envelope_samples':int(take.sum()),
                    'visible_samples':int(s.sum()),'core_samples':int(c.sum()),'covered_visible_samples':int(s.sum()),
                    'field_ids_json':json.dumps(sorted(set(fids[take]))),'min_quality_edge_deg':float(edges[take].min()),
                    'recommended_pitch_quality_loss_samples':int(np.sum(s&(oldedge<grid.margin+grid.guard)))})
            # Consecutive sampling runs are NOT a promise of continuously observable intervals.
            if len(inds):
                breaks=np.r_[0,np.flatnonzero((np.diff(ts[inds])>1800.1)|(fids[1:]!=fids[:-1]))+1,len(inds)]
                for a,b in zip(breaks[:-1],breaks[1:]):intervals.append({'object_id':k,'object_cloud':cloud,'field_id':str(fids[a]),
                    'first_sample_utc':iso_times[inds[a]]+'Z','last_sample_utc':iso_times[inds[b-1]]+'Z',
                    'samples':int(b-a),'min_edge_deg':float(edges[a:b].min()),'usable_measured':None})
            def longest(days):
                if not days:return 0
                ints=np.array(days,dtype='datetime64[D]').astype(int);splits=np.r_[0,np.where(np.diff(ints)!=1)[0]+1,len(ints)]
                return int(np.diff(splits).max())
            objects.append({'object_id':k,'cloud':cloud,'visible_samples':int(sampled.sum()),'core_samples':int(core.sum()),
                'covered_visible_samples':int(sampled.sum()),'visible_days':len(visible_days),'core_days':len(core_days),
                'longest_consecutive_visible_days':longest(visible_days),'longest_consecutive_core_days':longest(core_days),
                'fixed_fields_used':len(set(fids[sampled])),'primary_day_field_changes':int(np.sum(np.array(primary)[1:]!=np.array(primary)[:-1])),
                'envelope_samples':len(inds),'minimum_edge_deg':float(edges.min()) if len(edges) else None,
                'proposed_pitch_quality_loss_samples':int(np.sum(sampled&(oldedge<grid.margin+grid.guard))),
                'uncovered_reason':'' if sampled.any() else ('no_V19.2_dark_alt32.5_sample; see conservative envelope' if len(inds) else 'no_possible_brightness_dark_altitude_segment')})
    pd.DataFrame(daily).to_csv(output/'object_days.csv.gz',index=False)
    pd.DataFrame(intervals).to_csv(output/'object_field_intervals.csv.gz',index=False)
    ob=pd.DataFrame(objects);ob.to_csv(output/'objects.csv',index=False)
    occupancy=[]
    occupancy_values={}
    for (fid,t),value in counts.items():occupancy_values.setdefault(fid,[]).append(value)
    for f in doc['fields']:
        fid=f['field_id'];vals=occupancy_values.get(fid,[])
        occupancy.append({'field_id':fid,'planning_cloud':f['planning_cloud'],
            'unique_visible_objects':len(field_members.get(fid,set())), 'active_samples':len(vals),
            'median_simultaneous_assigned_objects':float(np.median(vals)) if vals else 0,
            'max_simultaneous_assigned_objects':max(vals,default=0),
            'recommended_role':'core_active' if len(vals)>=20 else 'candidate_boundary',
            'assignment':'unique ownership, overlap companions available through query API'})
    pd.DataFrame(occupancy).to_csv(output/'field_occupancy.csv',index=False)
    summary={'status':doc['freeze_status'],'counts':df.groupby('planning_cloud').size().to_dict(),
        'objects':len(ob),'with_strict_visible_samples':int((ob.visible_samples>0).sum()),
        'with_core_samples':int((ob.core_samples>0).sum()),'visible_samples':int(ob.visible_samples.sum()),
        'covered_visible_samples':int(ob.covered_visible_samples.sum()),'conservative_envelope_samples':int(ob.envelope_samples.sum()),
        'proposed_pitch_quality_loss_samples':int(ob.proposed_pitch_quality_loss_samples.sum()),
        'grid':asdict(grid),'nominal_min_cell_edge_deg':min(r['minimum_edge_deg'] for r in rows_audit),
        'geometry_certification':'nominal cell union + transported measured WCS checks; actual masks and repeated pointing still unknown',
        'quality_yield':'not recomputed; legacy B lightcurve/period yields do not apply',
        'r_template_count':len(df),'griz_template_count':4*len(df),
        'mean_active_simultaneous_assigned_objects':float(np.mean(list(counts.values()))),
        'median_object_field_changes_visible_sample':float(ob.loc[ob.visible_samples>0,'primary_day_field_changes'].median())}
    write_json(output/'summary.json',summary);print(json.dumps(summary,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--reuse-samples',action='store_true')
    p.add_argument('--previous-registry',type=Path,default=REGISTRY_PATH)
    p.add_argument('--release',default='20261111-rc1');a=p.parse_args()
    if not a.reuse_samples:prepare(a.output)
    else:
        _,_,_,fresh=load_inputs()
        prior=json.loads((a.output/'provenance.json').read_text())
        for key in fresh:
            if fresh[key]!=prior[key]:raise ValueError('Cached sample inputs changed: '+key)
    previous=json.loads(a.previous_registry.read_text()) if a.previous_registry.exists() else None
    build(a.output,previous=previous,release=a.release)

if __name__=='__main__':main()
