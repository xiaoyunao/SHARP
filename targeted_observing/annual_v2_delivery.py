"""Validate annual replay products and render a fixed-log-scale coverage movie."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
from collections import defaultdict

import numpy as np
from astropy.table import Table
from PIL import Image

from .catalog import write_csv
from .fixed_fields import from_tan
from .server_common import write_json
from .shared_fields import SharedFields
from .full_night_test import YSO_FIELDS


def read_csv(path):
    with Path(path).open() as stream:return list(csv.DictReader(stream))


def summarize(root,output):
    root=Path(root);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    summaries=json.loads((root/'case_summaries.json').read_text());periods=json.loads((root/'period_case_summaries.json').read_text())
    provenance=json.loads((root/'provenance.json').read_text())
    assert len(summaries)==33 and len({r['case'] for r in summaries})==33
    assert {r['model_hash'] for r in summaries}=={provenance['model_hash']}
    period_by={r['case']:r['counts'] for r in periods};rows=[];checks=[];field_totals=defaultdict(list)
    for s in summaries:
        case=root/'cases'/s['case'];data=dict(np.load(case/'field_counts.npz'));nightly=read_csv(case/'nightly.csv')
        assert len(nightly)==365 and s['nights']==365
        assert 0<=s['detected_r_objects']<=s['detected_objects']<=1063
        assert s['tjn_r_measurements']<=s['tjn_measurements']
        assert len(set(r['date'] for r in nightly))==365
        assert nightly[0]['date']=='2026-11-01' and nightly[-1]['date']=='2027-10-31'
        assert data['actual'].shape==data['planned'].shape==(365,4397)
        assert data['yso_actual'].shape==(365,2)
        assert np.all(data['actual']>=0) and np.all(data['actual']<=data['planned'])
        assert int(data['actual'].sum()+data['yso_actual'].sum())==s['usable_frames']
        assert int(data['planned'].sum()+data['yso_planned'].sum())==s['planned_frames']
        field_totals[s['scenario']].append(np.r_[data['actual'].sum(axis=0),data['yso_actual'].sum(axis=0)])
        record={**{k:s[k] for k in ('case','scenario','realization','detected_objects','detected_r_objects','tjn_measurements','tjn_r_measurements',
                                     'sampling_ready_objects','planned_frames','usable_frames','shared_fields_visited')},
                'mp_sampling_failures':sum(int(r['mp_sampling_failures']) for r in nightly)}
        for amplitude,count in period_by[s['case']].items():
            for key,value in count.items():
                assert value<=s['sampling_ready_objects']
                record[f'amplitude_{amplitude}_{key}']=value
        rows.append(record);checks.append(dict(case=s['case'],night_rows=365,exposure_count_balance=True))
    write_csv(output/'annual_case_summary.csv',rows)
    by_scenario={}
    for scenario in ('clear_baseline','baseline','optimistic','conservative'):
        group=[r for r in rows if r['scenario']==scenario]
        by_scenario[scenario]=dict(realizations=len(group))
        for key in rows[0]:
            if key in ('case','scenario','realization'):continue
            by_scenario[scenario][key]=np.percentile([r[key] for r in group],[10,50,90]).tolist()
    summary=dict(start_date='2026-11-01',stop_date='2027-10-31',nights=365,candidates=1063,shared_fields=4397,
                 scenarios=by_scenario,period_range_h=[2,1000],amplitudes_mag=[.1,.3,.6],
                 period_counts_overlap=True,weather_percentiles='P10/P50/P90 across scenario realizations, not confidence bounds',
                 movie='clear_baseline actual exposures versus baseline realization 0; cumulative counts, one frame per night',
                 limitation='Conditional planning forecast for the frozen candidate set, weather, noise and injected spin/shape assumptions; no guaranteed new periods.')
    write_json(output/'annual_summary.json',summary);write_json(output/'validation.json',dict(passed=True,cases=checks))
    for name in ('provenance.json','case_summaries.json','period_case_summaries.json'):
        shutil.copy2(root/name,output/name)
    for name in ('clear_baseline','baseline_r00'):
        target=output/name;target.mkdir(exist_ok=True)
        for filename in ('field_counts.npz','nightly.csv','objects.csv','period_recovery.csv','period_summary.json','period_assumptions.json'):
            p=root/'cases'/name/filename
            if p.exists():shutil.copy2(p,target/filename)
    for case in (s['case'] for s in summaries):
        target=output/'audit/cases'/case;target.mkdir(parents=True,exist_ok=True)
        for filename in ('field_counts.npz','nightly.csv','objects.csv','summary.json','curves.npz',
                         'sampling_certificates.json','period_recovery.csv','period_summary.json','period_assumptions.json'):
            source=root/'cases'/case/filename
            if source.exists():shutil.copy2(source,target/filename)
    cat=SharedFields();cases={name:dict(np.load(root/'cases'/name/'field_counts.npz')) for name in ('clear_baseline','baseline_r00')}
    fields=[]
    for i,fid in enumerate(cat.ids):
        row=dict(field_id=fid,ra_deg=cat.by_id[fid]['ra_deg'],dec_deg=cat.by_id[fid]['dec_deg'])
        for name,arrays in cases.items():row[name+'_exposures']=int(arrays['actual'][:,i].sum())
        fields.append(row)
    for i,(name,values) in enumerate(YSO_FIELDS.items()):
        row=dict(field_id=name,**values)
        for case,arrays in cases.items():row[case+'_exposures']=int(arrays['yso_actual'][:,i].sum())
        fields.append(row)
    baseline=np.asarray(field_totals['baseline'])
    percentiles=np.percentile(baseline,[10,50,90],axis=0)
    for i,row in enumerate(fields):
        for j,label in enumerate(('P10','P50','P90')):row['baseline_'+label]=float(percentiles[j,i])
        row['baseline_fraction_visited']=float(np.mean(baseline[:,i]>0))
    write_csv(output/'field_exposure_counts.csv',fields)
    return summary


def clip(poly,bound,keep_greater):
    result=[]
    for a,b in zip(poly,poly[1:]+poly[:1]):
        inside_a=a[0]>=bound if keep_greater else a[0]<=bound
        inside_b=b[0]>=bound if keep_greater else b[0]<=bound
        if inside_a:result.append(a)
        if inside_a!=inside_b:
            f=(bound-a[0])/(b[0]-a[0]);result.append((bound,a[1]+f*(b[1]-a[1])))
    return result


def polygons():
    catalog=SharedFields();result=[];indices=[]
    for i,row in enumerate(catalog.table):
        ra=np.asarray([row[f'corner_ra_{j}'] for j in range(1,5)],float)
        dec=np.asarray([row[f'corner_dec_{j}'] for j in range(1,5)],float)
        result.append((ra,dec));indices.append(i)
    for i,(name,values) in enumerate(YSO_FIELDS.items()):
        ra,dec=from_tan(values['ra_deg'],values['dec_deg'],[-1.48,1.48,1.48,-1.48],[-1.48,-1.48,1.48,1.48])
        result.append((ra,dec));indices.append(len(catalog.ids)+i)
    pieces=[];mapping=[]
    for (ra,dec),index in zip(result,indices):
        x=np.unwrap(-np.deg2rad(ra));x-=2*np.pi*np.round(x.mean()/(2*np.pi))
        y=np.deg2rad(dec)
        for shift in (-2*np.pi,0,2*np.pi):
            poly=list(zip(x+shift,y));poly=clip(clip(poly,-np.pi,True),np.pi,False) if poly else []
            if len(poly)>=3 and np.ptp(np.asarray(poly)[:,0])>1e-8:
                pieces.append(np.asarray(poly));mapping.append(index)
    return pieces,np.asarray(mapping)


def render(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import LogNorm
    from matplotlib.transforms import IdentityTransform
    from matplotlib.cm import ScalarMappable
    output=Path(output);names=['clear_baseline','baseline_r00'];counts=[];dates=None
    for name in names:
        data=dict(np.load(output/name/'field_counts.npz'))
        counts.append(np.cumsum(np.column_stack([data['actual'],data['yso_actual']]),axis=0));dates=data['dates']
    maximum=max(int(v.max()) for v in counts);norm=LogNorm(vmin=1,vmax=max(2,maximum))
    cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#f1f2f4')
    fig=plt.figure(figsize=(14,7),dpi=90,facecolor='white');pieces,mapping=polygons();collections=[];labels=[]
    for i,title in enumerate(('Clear nights / baseline photometry','Weather + baseline noise / realization 0')):
        ax=fig.add_axes([.025+.5*i,.30,.45,.46],projection='aitoff');ax.set_facecolor('#fafafa')
        ax.grid(alpha=.28,linewidth=.5);ax.set_axisbelow(False);ax.tick_params(labelsize=7,pad=2)
        ticks=np.arange(-150,180,30);ax.set_xticks(np.deg2rad(ticks));ax.set_xticklabels([str(int((-v)%360)) for v in ticks])
        ax.set_title(title,fontsize=12,pad=20)
        collection=PolyCollection(pieces,array=np.ma.masked_all(len(mapping)),cmap=cmap,norm=norm,
                                  edgecolors='#dadde2',linewidths=.10,closed=True,zorder=2)
        ax.add_collection(collection);collections.append((ax,collection))
        labels.append(fig.text(.25+.5*i,.25,'',ha='center',fontsize=10,color='#333333'))
    fig.text(.5,.96,'SHARP v2 · Annual cumulative field coverage',ha='center',fontsize=19,weight='bold')
    stamp=fig.text(.5,.88,'',ha='center',fontsize=13)
    colorax=fig.add_axes([.28,.16,.44,.025]);bar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=colorax,orientation='horizontal')
    ticks=[v for v in (1,3,10,30,100,300,1000,3000,10000,30000,100000) if v<=maximum]
    bar.set_ticks(ticks);bar.set_ticklabels([f'{v:,}' for v in ticks]);bar.ax.tick_params(labelsize=9)
    bar.set_label('Cumulative exposure count N · logarithmic color scale',fontsize=10,labelpad=7)
    fig.text(.5,.055,'Grey = no exposure · MP and TJN share the same field IDs · 4,397 common fields + 2 independent YSO fields',ha='center',fontsize=9,color='#555555')
    fig.text(.5,.025,f'{dates[0]} – {dates[-1]} · simulation; no assumed follow-up requests · weather panel is one illustrative realization',ha='center',fontsize=9,color='#555555')
    fig.canvas.draw()
    # Cache the non-affine Aitoff geometry once; only face colors change per night.
    for ax,collection in collections:
        vertices=[ax.transData.transform_path(path).vertices for path in collection.get_paths()]
        collection.set_verts(vertices,closed=True);collection.set_transform(IdentityTransform())
    frames=[];palette=None
    for day in range(len(dates)):
        stamp.set_text(f'{dates[day]}  |  Night {day+1:03d} / {len(dates)}')
        for i,((ax,collection),values,label) in enumerate(zip(collections,counts,labels)):
            total=values[day];collection.set_array(np.ma.masked_less(total[mapping],1))
            label.set_text(f'Shared fields reached: {np.count_nonzero(total[:4397]):,} / 4,397   ·   Exposures: {int(total.sum()):,}\nYSO Ori: {int(total[-2]):,}   ·   YSO Cyg: {int(total[-1]):,}')
        fig.canvas.draw();rgb=np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()
        image=Image.fromarray(rgb)
        if palette is None:palette=image.quantize(colors=256,method=Image.Quantize.MEDIANCUT)
        frames.append(image.quantize(palette=palette,dither=Image.Dither.NONE))
        if day in (0,len(dates)//2,len(dates)-1):image.save(output/f'coverage_{day+1:03d}.png')
    path=output/'annual_coverage_logn.gif'
    frames[0].save(path,save_all=True,append_images=frames[1:],duration=100,loop=0,optimize=False)
    with Image.open(path) as probe:
        assert probe.n_frames==len(dates)
        info=dict(frames=probe.n_frames,width=probe.width,height=probe.height,duration_ms=100,loop=probe.info.get('loop'),log_vmin=1,log_vmax=maximum)
    write_json(output/'animation_validation.json',info)
    plt.close(fig)
    return info


def main():
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);p.add_argument('--no-render',action='store_true');args=p.parse_args()
    summary=summarize(args.root,args.output)
    if not args.no_render:print(render(args.output))
    manifest={str(p.relative_to(args.output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in args.output.rglob('*')
              if p.is_file() and p.name!='SHA256SUMS.json'}
    write_json(args.output/'SHA256SUMS.json',manifest)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
