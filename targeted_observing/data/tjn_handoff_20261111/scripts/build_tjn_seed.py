"""Build a REVIEW-ONLY fixed TJN grid from legacy B footprint envelopes.

This is a geometry/data handoff, not a scheduler. It neither reads nor writes
production or a telescope. Full annual candidate ephemerides and measured WCS
are needed before freezing the registry. Dependencies: numpy, pandas, matplotlib.
"""
from pathlib import Path
import argparse
import json
import math
import hashlib
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

SX = 2.67264
SY = 2.67496
W = 9216 * 1.15 / 3600
H = 9232 * 1.15 / 3600


def row_parameters(j):
    dec = j * SY
    n = int(np.ceil(360 * np.cos(np.deg2rad(dec)) / SX))
    return dec, n, 360. / n


def inv_local(ra0, dec0, x, y, projection='tan'):
    """Return ICRS-like directions; no aberration or epoch transformation."""
    a,d=np.deg2rad([ra0,dec0])
    p=np.array([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
    east=np.array([-np.sin(a),np.cos(a),0.])
    north=np.array([-np.sin(d)*np.cos(a),-np.sin(d)*np.sin(a),np.cos(d)])
    x,y=np.broadcast_arrays(np.asarray(x),np.asarray(y))
    if projection=='tan':
        vec=p[:,None]+east[:,None]*x.ravel()+north[:,None]*y.ravel()
    else:
        lon,lat=np.deg2rad(x.ravel()),np.deg2rad(y.ravel())
        vec=(p[:,None]*(np.cos(lat)*np.cos(lon))+
             east[:,None]*(np.cos(lat)*np.sin(lon))+north[:,None]*np.sin(lat))
    vec/=np.linalg.norm(vec,axis=0)
    return (np.rad2deg(np.arctan2(vec[1],vec[0]))%360).reshape(x.shape),np.rad2deg(np.arcsin(vec[2])).reshape(x.shape)


def boundary(ra,dec,w=W,h=H,projection='tan',samples=25):
    if projection=='tan':
        # Pixel scale is linear in the TAN projection plane, not atan(x).
        wx,hy=np.deg2rad([w/2,h/2])
    else:
        wx,hy=w/2,h/2
    q=np.linspace(-1,1,samples)
    x=np.r_[q*wx,np.full(samples,wx),-q*wx,np.full(samples,-wx)]
    y=np.r_[np.full(samples,-hy),q*hy,np.full(samples,hy),-q*hy]
    return inv_local(ra,dec,x,y,projection)


def sky_to_tan(ra,dec,ra0,dec0):
    da=np.deg2rad((np.asarray(ra)-ra0+180)%360-180)
    d,d0=np.deg2rad(dec),np.deg2rad(dec0)
    den=np.sin(d0)*np.sin(d)+np.cos(d0)*np.cos(d)*np.cos(da)
    x=np.cos(d)*np.sin(da)/den
    y=(np.cos(d0)*np.sin(d)-np.sin(d0)*np.cos(d)*np.cos(da))/den
    return x,y


def assign_cells(ra,dec):
    ra,dec=np.broadcast_arrays(np.asarray(ra),np.asarray(dec))
    js=np.floor(dec/SY+.5).astype(int)
    keys=[]
    for a,j in zip(ra.ravel(),js.ravel()):
        _,n,step=row_parameters(int(j))
        k=int(np.floor((a%360)/step+.5))%n
        keys.append((int(j),k))
    return keys


def main(source,out):
    out.mkdir(parents=True,exist_ok=True)
    old=pd.read_csv(source/'tjn_B_unique_fields_with_cloud.csv')
    groups={}
    evidence=[]
    # Conservative row/column envelope intersection. Envelope corners can add
    # extra seed cells; this intentionally does not optimize template count.
    for _,f in old.iterrows():
        ra0,dec0=float(f.ra_deg),float(f.dec_deg)
        br,bd=boundary(ra0,dec0,3.,3.,'skyoffset',81)
        delta=(br-ra0+180)%360-180
        lo,hi=ra0+delta.min(),ra0+delta.max()
        dmin,dmax=float(bd.min()),float(bd.max())
        j0=math.ceil((dmin-SY/2)/SY-1e-10)
        j1=math.floor((dmax+SY/2)/SY+1e-10)
        for j in range(j0,j1+1):
            dec,n,step=row_parameters(j)
            k0=math.ceil(lo/step-.5-1e-10)
            k1=math.floor(hi/step+.5+1e-10)
            for rawk in range(k0,k1+1):
                k=rawk%n;key=(f.cloud_proxy,j,k)
                info=groups.setdefault(key,{'fields':set(),'objects':set(),'first':f.first_night,'last':f.last_night,'field_nights':0})
                info['fields'].add(f.field_id)
                info['objects'].update(json.loads(f.associated_objects_json))
                info['first']=min(info['first'],f.first_night)
                info['last']=max(info['last'],f.last_night)
                info['field_nights']+=int(f.field_nights)
                evidence.append({'cloud_proxy':f.cloud_proxy,'row_index':j,'column_index':k,'legacy_field_id':f.field_id,'basis':'legacy nominal 3-degree footprint bounding-envelope intersects grid ownership cell'})
    rows=[];idmap={}
    for cloud in ('L4','L5'):
        for number,key in enumerate(sorted(k for k in groups if k[0]==cloud),1):
            _,j,k=key;dec,n,step=row_parameters(j);ra=k*step
            info=groups[key];fid=f'TJ{cloud}_{number:03d}';idmap[key]=fid
            corners_ra,corners_dec=inv_local(ra,dec,np.deg2rad(W/2)*np.array([-1,1,1,-1]),np.deg2rad(H/2)*np.array([-1,-1,1,1]))
            rows.append({'field_id':fid,'grid_cell_id':f'EQv0_R{j:+04d}_C{k:03d}','planning_cloud':cloud,
                'row_index':j,'column_index':k,'columns_in_full_row':n,'ra_deg':ra,'dec_deg':dec,'pa_deg':0.,
                'ra_step_deg':step,'dec_step_deg':SY,'frame':'ICRS','geometry_model':'ideal_TAN_9216x9232_1.15arcsec',
                'nominal_tan_plane_width_deg':W,'nominal_tan_plane_height_deg':H,'catalog_status':'PROPOSAL_NOT_FROZEN',
                'scope':'legacy_B_footprint_seed_only','source_legacy_fields':len(info['fields']),
                'legacy_field_night_links':info['field_nights'],'legacy_associated_object_ids_json':json.dumps(sorted(info['objects']),separators=(',',':')),
                'source_first_night':info['first'],'source_last_night':info['last'],
                'corner_ra_deg_json':json.dumps(np.round(corners_ra,10).tolist(),separators=(',',':')),
                'corner_dec_deg_json':json.dumps(np.round(corners_dec,10).tolist(),separators=(',',':')),
                'template_status':'unassessed','note':'legacy members are provenance only; not proven occupants of this fixed field'})
    grid=pd.DataFrame(rows)
    grid.to_csv(out/'tjn_fixed_grid_seed_v0.csv',index=False,float_format='%.10f')
    for cloud in ('L4','L5'):
        grid[grid.planning_cloud==cloud].to_csv(out/f'TJ{cloud}_seed_v0.csv',index=False,float_format='%.10f')
    for e in evidence:e['proposed_field_id']=idmap[(e['cloud_proxy'],e['row_index'],e['column_index'])]
    pd.DataFrame(evidence).to_csv(out/'legacy_to_seed_cells.csv',index=False)
    lookup={(r.planning_cloud,int(r.row_index),int(r.column_index)):r for r in grid.itertuples()}
    missing=0;tested=0;minx=1e9;miny=1e9
    # Independent dense legacy footprint coverage against ownership cell + TAN.
    # This is nominal footprint coverage, not candidate-track coverage.
    for _,f in old.iterrows():
        q=np.linspace(-1.5,1.5,31);xx,yy=np.meshgrid(q,q)
        rr,dd=inv_local(float(f.ra_deg),float(f.dec_deg),xx,yy,'skyoffset')
        keys=assign_cells(rr,dd)
        for a,d,(j,k) in zip(rr.ravel(),dd.ravel(),keys):
            tested+=1;r=lookup.get((f.cloud_proxy,j,k))
            if r is None:missing+=1;continue
            x,y=sky_to_tan(a,d,r.ra_deg,r.dec_deg)
            dx=np.deg2rad(W/2)-abs(float(x));dy=np.deg2rad(H/2)-abs(float(y))
            minx=min(minx,np.rad2deg(dx)*60);miny=min(miny,np.rad2deg(dy)*60)
            if dx < -1e-12 or dy < -1e-12:missing+=1
    # Current TJ_test center maps to one ownership cell; target positions were
    # not present in the exported dataset, so do not claim target containment.
    j,k=assign_cells([79.002262],[23.719287])[0]
    tonight=lookup[('L5',j,k)]
    summary={'status':'PROPOSAL_NOT_FROZEN','planning_start':'2026-11-01','planning_stop_exclusive':'2027-11-01',
      'source_commit':'9734af0a04253bf3e7f2caccc68c0a87f48bfe06','latest_reviewed_commit':'7747cc1650fdeb255e677e496f50a376a99cb57b',
      'source_B_generation_commit':'4605b372fde47692ddcd66ed809905b07caf05bb',
      'input_sha256':hashlib.sha256((source/'tjn_B_unique_fields_with_cloud.csv').read_bytes()).hexdigest(),
      'grid':{'dec_anchor_deg':0.,'ra_anchor_deg':0.,'nominal_dec_step_deg':SY,'nominal_east_west_pitch_deg':SX,
        'row_columns':'ceil(360*cos(dec_row)/east_west_pitch)','ra_step':'360/row_columns','ra_centers':'k*ra_step',
        'pa_deg':0.,'nominal_detector_arcsec_per_pixel':1.15,'mp_reference_arcsec_per_pixel':1.16},
      'counts':{c:int(sum(grid.planning_cloud==c)) for c in ('L4','L5')},'legacy_fields':len(old),
      'nominal_coverage_audit':{'legacy_footprint_test_points':tested,'uncovered_points':missing,
        'minimum_sampled_edge_margin_x_arcmin':minx,'minimum_sampled_edge_margin_y_arcmin':miny,
        'scope':'31x31 samples per legacy nominal spherical-offset footprint, ideal 1.15-arcsec TAN new footprints; no real WCS or target ephemerides'},
      'tonight_test':{'legacy_field_id':'AF_148ecdaae5','legacy_name':'TJ_test','legacy_ra_deg':79.002262,'legacy_dec_deg':23.719287,
        'proposed_nearby_field_id':tonight.field_id,'proposed_ra_deg':tonight.ra_deg,'proposed_dec_deg':tonight.dec_deg,
        'target_ids':['17365','4722'],'target_containment_verified':False},
      'limitations':['Seed only: covers legacy B footprint envelopes, not all 1063 annual candidate tracks.',
        'cloud_proxy is inherited instantaneous heliocentric longitude proxy; no new resonance classification.',
        'Proposed field IDs are not frozen and must not replace tonight executable list.',
        'Actual usable WCS, pointing repeatability, PA and template masks were not available.',
        'Field coverage is not number of successful lightcurves or recovered periods.']}
    summary['regions']={c:{'ra_center_min_deg':float(g.ra_deg.min()),'ra_center_max_deg':float(g.ra_deg.max()),
                        'dec_center_min_deg':float(g.dec_deg.min()),'dec_center_max_deg':float(g.dec_deg.max())}
                        for c,g in grid.groupby('planning_cloud')}
    (out/'seed_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    draw(grid,old,out,tonight)
    print(json.dumps(summary,indent=2))


def draw(grid,old,out,tonight):
    plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans','font.size':12,
      'axes.labelsize':16,'axes.titlesize':17,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':10,
      'figure.dpi':140,'savefig.dpi':240,'figure.facecolor':'white','axes.facecolor':'white','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(16,8.5),layout='constrained')
    colors={'L4':'#1f77b4','L5':'#ff8c00'}
    for ax,cloud in zip(axes,['L4','L5']):
        g=grid[grid.planning_cloud==cloud];b=old[old.cloud_proxy==cloud]
        segments=[np.column_stack(boundary(r.ra_deg,r.dec_deg,samples=10)) for r in g.itertuples()]
        ax.add_collection(LineCollection(segments,colors=colors[cloud],linewidths=.60,alpha=.65))
        ax.scatter(g.ra_deg,g.dec_deg,s=6,color=colors[cloud],alpha=.65,zorder=3)
        ax.scatter(b.ra_deg,b.dec_deg,s=12,marker='x',color='.2',linewidths=.7,zorder=4)
        ax.set_xlim(g.ra_deg.max()+4,g.ra_deg.min()-4)
        ax.set_ylim(min(-19,g.dec_deg.min()-3),max(58,g.dec_deg.max()+3))
        ax.set_xlabel('Right ascension (deg)');ax.set_ylabel('Declination (deg)')
        ax.set_title(f'{cloud}: {len(g)} proposed seed fields / {len(b)} legacy pointings')
        ax.grid(alpha=.15);ax.tick_params(which='both',direction='in',top=True,right=True)
        for d in g.dec_deg.unique(): ax.axhline(d,color=colors[cloud],lw=.3,alpha=.15,zorder=0)
        if cloud=='L5':
            ax.scatter([79.002262],[23.719287],marker='*',s=100,color='#2ca02c',edgecolors='white',linewidths=.5,zorder=6)
            ax.annotate('TJ_test',(79.002262,23.719287),xytext=(68,28),ha='center',fontsize=11,
                        arrowprops={'arrowstyle':'-','lw':.8,'color':'#2ca02c'},color='#207120')
    axes[0].legend(handles=[Line2D([],[],color=colors['L4'],lw=1,label='Proposed fixed footprints (ideal TAN)'),
            Line2D([],[],color='.2',marker='x',ls='',label='Legacy B centers (7-day pointing versions)')],loc='upper right')
    fig.suptitle('TJN fixed-field seed atlas  |  2026 Nov - 2027 Oct',fontsize=19)
    fig.supxlabel('Review draft from legacy B coverage; complete annual candidate tracks and measured WCS are required before freezing IDs.',fontsize=11)
    fig.savefig(out/'01_tjn_fixed_seed_atlas.png',dpi=240,bbox_inches='tight')
    fig.savefig(out/'01_tjn_fixed_seed_atlas.pdf',bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,8),layout='constrained')
    g=grid[(grid.planning_cloud=='L5')&(abs(grid.ra_deg-79)<7)&(abs(grid.dec_deg-24)<7)]
    for r in g.itertuples():
        a,d=boundary(r.ra_deg,r.dec_deg);ax.plot(a,d,color=colors['L5'],lw=1)
        ax.plot(r.ra_deg,r.dec_deg,'+',color=colors['L5'])
        if 74 < r.ra_deg < 84 and 19 < r.dec_deg < 29:
            ax.text(r.ra_deg,r.dec_deg+.35,r.field_id,ha='center',va='bottom',fontsize=9,color='#975000',clip_on=True)
    a,d=boundary(79.002262,23.719287,3,3,'skyoffset');ax.plot(a,d,color='#2ca02c',lw=1.6,ls='--',label='Legacy TJ_test nominal footprint')
    ax.scatter(79.002262,23.719287,color='#2ca02c',marker='*',s=130,zorder=5,label='Legacy center 79.002262, +23.719287')
    ax.set(xlim=(85,73),ylim=(18,30),xlabel='Right ascension (deg)',ylabel='Declination (deg)',title='Fixed rows around the current L5 test field')
    ax.grid(alpha=.15);ax.tick_params(which='both',direction='in',top=True,right=True)
    ax.legend(loc='lower left',fontsize=10)
    fig.supxlabel('Orange field IDs are proposed; asteroid containment needs the cached 17365 / 4722 ephemerides.',fontsize=10)
    fig.savefig(out/'02_tj_test_grid_comparison.png',dpi=240,bbox_inches='tight')
    fig.savefig(out/'02_tj_test_grid_comparison.pdf',bbox_inches='tight');plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();main(args.source,args.output)
