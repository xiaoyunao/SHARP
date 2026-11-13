"""Compose offline twilight + registered Trojan fields + the fixed Orion YSO test."""
import argparse
import hashlib
import json
import zipfile
from datetime import datetime, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
from astropy.time import Time

from .catalog import BASE, interpolate, write_csv
from .controller_fields import controller_name, canonical_name, is_fixed_name
from .fixed_fields import FixedFields, REGISTRY_PATH, edge_distance
from .full_night_test import read_lines, replay, command, movement, YSO_FIELDS, READOUT_S, YSO_FRAMES_PER_BAND
from .geometry import setup_iers, local_time, coord
from .offline_twilight import make_sky
from .audit_controller_equatorial import audit
from .annual_trojan import load_tracks
from .bright_stars import load as load_stars
from .science import load_periods, assess

TZ=ZoneInfo('Asia/Shanghai')


def make_command(field,band,seconds,block,count=1):
    return read_lines(command(field,band,count,seconds),block)[0]


def rows_for(commands,start):
    return replay(commands,start)


def fill_to_handover(commands,start,fields,next_command,next_shutter,block,dwell=1):
    """Keep cycling until the next fixed shutter time; integer tail, no WAIT token."""
    if not isinstance(dwell,int) or isinstance(dwell,bool) or dwell<1:
        raise ValueError('Dwell must be a positive integer')
    index=0
    while True:
        trial=make_command(fields[(index//dwell)%len(fields)],'r',120,block)
        rr=rows_for(commands+[trial,next_command],start)
        if next(r['start_unix'] for r in rr if r['command_index']==len(commands)+2)>next_shutter:
            break
        commands.append(trial);index+=1
    if not index:
        raise ValueError('No full Trojan exposure fits before handover')
    # Finish on the last visited fixed centre. A short final exposure fills the
    # residual budget without inventing a controller wait or moving the field.
    while True:
        last=rows_for(commands,start)[-1]
        overhead=movement(last,next_command)+(10 if last['filter']!=next_command['filter'] else 0)
        slack=int(np.floor(next_shutter-last['release_unix']-overhead))
        if slack<=0:
            break
        field={k:last[k] for k in ('field_id','ra_deg','dec_deg')}
        if slack>=5+READOUT_S:
            commands.append(make_command(field,'r',min(120,int(slack-READOUT_S)),block))
        else:
            # Split the previous exposure to spend a residual smaller than one
            # 5s exposure + the configured readout, while keeping both exposures >=5s.
            previous=commands[-1]
            if previous['count']!=1 or previous['exposure_s']-(5+READOUT_S-slack)<5:
                raise ValueError('Cannot fill integer handover without an invalid exposure')
            previous=make_command(field,'r',int(previous['exposure_s']-(5+READOUT_S-slack)),block)
            commands[-1]=previous
            commands.append(make_command(field,'r',5,block))
    check=rows_for(commands+[next_command],start)
    first=next(r['start_unix'] for r in check if r['command_index']==len(commands)+1)
    if not 0<=next_shutter-first<1.01:
        raise ValueError('Handover drift')


def compress_dwells(commands):
    """Merge consecutive identical exposures into the controller repeat column."""
    result=[]
    keys=('block','field_id','ra_deg','dec_deg','filter','exposure_s')
    for c in commands:
        if result and all(c[k]==result[-1][k] for k in keys):
            result[-1]['count']+=c['count']
        else:
            result.append(dict(c))
    return result


def plot(rows,sky,output,day):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as md
    from .fixed_fields import boundary
    names=list(dict.fromkeys(r['field_id'] for r in rows if is_fixed_name(r['field_id'])))
    fig,axes=plt.subplots(2,1,figsize=(13,8),layout='constrained')
    t=np.linspace(rows[0]['start_unix'],rows[-1]['release_unix'],400)
    dt=lambda x:[datetime.fromtimestamp(float(v),TZ) for v in x]
    for name in names+['YSO_ORI']:
        f=next(r for r in rows if r['field_id']==name)
        y=sky.metrics(f['ra_deg'],f['dec_deg'],t)['alt'][0]
        line,=axes[0].plot(dt(t),y,lw=1,label=name)
        rr=[r for r in rows if r['field_id']==name]
        tt=[r['mid_unix'] for r in rr]
        aa=sky.metrics(f['ra_deg'],f['dec_deg'],tt)['alt'][0]
        axes[0].scatter(dt(tt),aa,s=5,color=line.get_color())
    axes[0].axhline(32.15,color='gray',ls='--',label='30 deg whole-field bound (approx.)')
    axes[0].set(ylim=(0,90),ylabel='Field centre altitude (deg)')
    axes[0].legend(ncol=3,fontsize=8)
    blocks=list(dict.fromkeys(r['block'] for r in rows))
    for i,b in enumerate(blocks):
        rr=[r for r in rows if r['block']==b]
        a,z=rr[0]['overhead_start_unix'],rr[-1]['release_unix']
        axes[1].barh(i,(z-a)/86400,left=md.date2num(dt([a]))[0])
        axes[1].text(md.date2num(dt([(a+z)/2]))[0],i,str(len(rr)),ha='center',va='center')
    axes[1].set(yticks=range(len(blocks)),yticklabels=blocks,xlabel='Beijing time; morning is next day')
    for ax in axes:
        ax.xaxis.set_major_formatter(md.DateFormatter('%H:%M',tz=TZ));ax.grid(alpha=.2)
    fig.suptitle(f'SHARP {day} | nominal schedule | numbers = frames')
    fig.savefig(output/'schedule.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,6),layout='constrained')
    reg=FixedFields.read(REGISTRY_PATH)
    targets=json.loads((output/'targets.json').read_text())
    colors=['tab:blue','tab:orange','tab:green']
    for name,color in zip(names,colors):
        ax=axes[0] if name.endswith('_606') else axes[1]
        f=reg.by_id[canonical_name(name)];p=boundary(f['ra_deg'],f['dec_deg'],reg.grid)
        ax.plot(p[:,0],p[:,1],color=color,label=name)
        for r in targets:
            if r['field_id']!=name:continue
            ax.plot(r['ra_mid_deg'],r['dec_mid_deg'],'k.')
            ax.annotate(r['object_id'],(r['ra_mid_deg'],r['dec_mid_deg']),xytext=(5,6),textcoords='offset points',fontsize=9)
    for ax,title in zip(axes,['Early evening field','Later two-field loop']):
        ax.set(xlabel='ICRS RA (deg)',ylabel='ICRS Dec (deg)',title=title)
        ax.invert_xaxis();ax.grid(alpha=.2);ax.margins(.15);ax.legend(fontsize=9,loc='upper right')
    fig.suptitle("Fixed registry footprints and tonight's retained targets")
    fig.savefig(output/'fields.png',dpi=160);plt.close(fig)


def build(day,twilight_dir,output):
    output.mkdir(parents=True,exist_ok=True)
    plan=json.loads((twilight_dir/'plan.json').read_text());assert plan['date']==day
    if plan['config']['readout_s'] != READOUT_S:
        raise ValueError('Regenerate twilight with the current readout setting')
    setup_iers(BASE/'cache'/plan['config']['iers_file']);sky=make_sky()
    registry=FixedFields.read(REGISTRY_PATH)
    identities=['TJL5_606','TJL5_400','TJL5_436']
    fields=[{**registry.by_id[k],'field_id':controller_name(k),'canonical_field_id':k} for k in identities]
    yso={'field_id':'YSO_ORI',**YSO_FIELDS['YSO_ORI']}
    ev=read_lines((twilight_dir/'EV.txt').read_text(),'EV')
    am=read_lines((twilight_dir/'AM.txt').read_text(),'AM')
    first=plan['blocks'][0]['window_unix'][0];commands=list(ev)
    # Three-hour YSO allocation includes readout; the sole g/r change adds 10s.
    yg=make_command(yso,'g',60,'YSO_g',YSO_FRAMES_PER_BAND);yr=make_command(yso,'r',60,'YSO_r',YSO_FRAMES_PER_BAND)
    local_day=datetime.fromisoformat(day).replace(tzinfo=TZ)
    main_first=(local_day+timedelta(hours=21,minutes=52)).timestamp()
    yso_first=(local_day+timedelta(days=1,hours=1,minutes=40)).timestamp()
    main_next=make_command(fields[1],'r',120,'TJN_MAIN')
    fill_to_handover(commands,first,fields[:1],main_next,main_first,'TJN_EARLY')
    fill_to_handover(commands,first,fields[1:],yg,yso_first,'TJN_MAIN',dwell=5)
    commands.extend([yg,yr])
    am_first=next(b for b in plan['blocks'] if b['side']=='AM')['rows'][0]['start_unix']
    fill_to_handover(commands,first,fields[1:],am[0],am_first,'TJN_LATE',dwell=5)
    commands.extend(am)
    rows=rows_for(commands,first)
    exported=compress_dwells(commands)
    text=''.join(command({k:c[k] for k in ('field_id','ra_deg','dec_deg')},c['filter'],c['count'],c['exposure_s']) for c in exported)
    destination=output/(day.replace('-','')+'.txt');destination.write_text(text)
    check=rows_for(read_lines(text,'ALL'),first)
    assert len(check)==len(rows) and all(abs(a['start_unix']-b['start_unix'])<.001 for a,b in zip(rows,check))
    for r,exported_row in zip(rows,check):
        r['command_index']=exported_row['command_index']
        r['repeat']=exported_row['repeat']
    for side in ('EV','AM'):
        actual=[r for r in rows if r['block']==side];ref=next(b for b in plan['blocks'] if b['side']==side)['rows']
        assert len(actual)==len(ref)
        assert all(abs(a['start_unix']-b['start_unix'])<1.1 and a['exposure_s']==b['exposure_s'] for a,b in zip(actual,ref))
    verification=audit(destination,datetime.fromtimestamp(first,TZ).isoformat(),output/'audit')
    assert verification['within_reference_ha'] and verification['sky_geometry_ok']
    tracks=load_tracks();periods=load_periods(BASE/'cache');targets=[]
    for f in fields:
        rr=[r for r in rows if r['field_id']==f['field_id']];times=np.array([r['mid_unix'] for r in rr]);jd=Time(times,format='unix').jd
        for key,tr in tracks.items():
            assert jd.min()>=tr['datetime_jd'][0] and jd.max()<=tr['datetime_jd'][-1]
            mag=np.interp(jd,tr['datetime_jd'],tr['V'])
            if mag.min()>19.2:continue
            ra,dec=interpolate(tr,times);edge=edge_distance(ra,dec,f['ra_deg'],f['dec_deg'],registry.grid)
            good=(edge>=registry.grid.margin+registry.grid.guard)&(mag<=19.2)
            if good.mean()<.8:continue
            targets.append(dict(field_id=f['field_id'],canonical_field_id=f['canonical_field_id'],object_id=key,
                name=tr['targetname'],V_min=float(mag.min()),V_max=float(mag.max()),
                ra_mid_deg=float(ra[len(ra)//2]),dec_mid_deg=float(dec[len(dec)//2]),
                planned_visits=len(rr),retained_fraction=float(good.mean()),minimum_edge_deg=float(edge.min()),
                period_evidence=assess(periods.get(key))))
        if not any(t['field_id']==f['field_id'] and t['V_max']<=18.5 for t in targets):
            raise ValueError('No retained core target for '+f['field_id'])
    (output/'targets.json').write_text(json.dumps(targets,indent=2)+'\n')
    stars=load_stars();sc=coord([s['ra_deg'] for s in stars],[s['dec_deg'] for s in stars]);risks=[]
    for f in fields:
        sep=coord(f['ra_deg'],f['dec_deg']).separation(sc).deg
        for star,d in zip(stars,sep):
            if d<4.122:
                risks.append(dict(field_id=f['field_id'],star=star,separation_deg=float(d),
                    in_nominal=bool(edge_distance(star['ra_deg'],star['dec_deg'],f['ra_deg'],f['dec_deg'],registry.grid)>=0)))
    (output/'bright_stars.json').write_text(json.dumps(risks,indent=2)+'\n')
    mapping=[dict(execution_name=f['field_id'],canonical_field_id=f['canonical_field_id'],
        ra_deg=f['ra_deg'],dec_deg=f['dec_deg'],pa_deg=f['pa_deg'],grid_cell_id=f['grid_cell_id'],registry_release=registry.document['release']) for f in fields]
    (output/'field_mapping.json').write_text(json.dumps(mapping,indent=2)+'\n')
    write_csv(output/'exposures.csv',[dict(r,start_local=local_time(r['start_unix']),end_local=local_time(r['end_unix']),readout_end_local=local_time(r['release_unix'])) for r in rows])
    blocks=[]
    for block in dict.fromkeys(r['block'] for r in rows):
        rr=[r for r in rows if r['block']==block]
        blocks.append(dict(block=block,frames=len(rr),first_exposure=local_time(rr[0]['start_unix']),
            readout_end=local_time(rr[-1]['release_unix']),integration_s=sum(r['exposure_s'] for r in rr),
            exposure_counts={str(e):sum(r['exposure_s']==e for r in rr) for e in sorted({r['exposure_s'] for r in rr})}))
    summary=dict(date=day,start_local=local_time(first),end_local=local_time(rows[-1]['release_unix']),
        controller_lines=len(exported),exposures=len(rows),blocks=blocks,fixed_fields=mapping,
        trojan_exposures_per_dwell=5,readout_s=READOUT_S,yso_frames_per_band=YSO_FRAMES_PER_BAND,
        core_targets=[{k:v for k,v in t.items() if k!='period_evidence'} for t in targets if t['V_max']<=18.5],
        registry_sha256=hashlib.sha256(REGISTRY_PATH.read_bytes()).hexdigest(),
        twilight_plan_sha256=hashlib.sha256((twilight_dir/'plan.json').read_bytes()).hexdigest(),
        geometry_checks_passed=True,mechanical_safety_certified=False,formal_instrument_freeze=False,
        AM_shutter_shift_s=next(r['start_unix'] for r in rows if r['block']=='AM')-am_first,
        assumptions=['Static sequential list, no WAIT or automatic clock gates.',
          f'Readout{READOUT_S:g}s, filter10s, empirical slew; same position has no extra slew.',
          'Main/late two-field loops take five consecutive r120 exposures per visit; final visit may be partial.',
          'Short r exposures at block ends align fixed handovers; main Trojan exposures120s.',
          'MP_ is an execution alias only; registry IDs/centres/PA are unchanged.',
          'Actual masks, halo contamination and physical mount path require observing feedback.'])
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    write_csv(output/'blocks.csv',blocks)
    plot(rows,sky,output,day)
    (output/'READ_ME.txt').write_text(f'SHARP {day}\nExecute {destination.name} once starting {summary["start_local"]} Beijing time.\n'
        'EV -> early fixed Trojan -> two fixed Trojan fields -> YSO_ORI g/r -> late Trojan -> AM.\n'
        f'YSO: {YSO_FRAMES_PER_BAND}x60s each band; {READOUT_S:g}s readout; one extra10s filter change.\n'
        'TJL5_400/436: five consecutive exposures per field before switching; final visits can be partial.\n'
        'Fixed Trojan execution names start MP_TJL5_; canonical IDs/centres stay unchanged.\n'
        'All filter tokens lowercase. No M31, no WAIT, no automatic deadline stop.\n'
        'Observe actual times; preset lists cannot compensate arbitrary delays.\n')
    files=[p for p in sorted(output.rglob('*')) if p.is_file() and p.suffix!='.zip' and p.name!='manifest.json']
    manifest={str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(output/f'SHARP_{day.replace("-","")}_fixed_fields.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in files+[output/'manifest.json']:z.write(p,str(p.relative_to(output)))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--date',required=True);p.add_argument('--twilight',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();build(a.date,Path(a.twilight),Path(a.output))
