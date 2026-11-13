"""Compose a local full-night test while preserving the approved twilight lists."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile
from datetime import datetime
from zoneinfo import ZoneInfo

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from .catalog import BASE, write_csv, interpolate
from .geometry import coord, footprint, in_field, local_time, setup_iers
from .offline_twilight import make_sky, slew, load_grid
from .bright_stars import load as load_bright_stars

TJ = {'field_id':'TJN_17365', 'ra_deg':79.002262, 'dec_deg':23.719287}
YSO_FIELDS = json.loads((BASE/'data/fixed_yso_fields.json').read_text())
YSO = {'field_id':'YSO_ORI', **YSO_FIELDS['YSO_ORI']}
TZ = ZoneInfo('Asia/Shanghai')
READOUT_S = float(json.loads((BASE/'offline_defaults.json').read_text())['readout_s'])
YSO_FRAMES_PER_BAND = int(5400 // (60 + READOUT_S))


def lowercase_filters(text):
    # Change ONLY column 5, preserving object labels, coordinates, counts and whitespace.
    return re.sub(r'^(\S+\s+\S+\s+\S+\s+\S+\s+)([GRIZgriz])(?=\s)',
                  lambda m:m.group(1)+m.group(2).lower(),text,flags=re.MULTILINE)


def read_lines(text, block):
    result=[]
    for line in text.splitlines():
        if not line.strip():continue
        a=line.split()
        if len(a)!=8:raise ValueError('Expected eight controller columns')
        p=SkyCoord(a[1],a[2],unit=(u.hourangle,u.deg))
        count=int(a[5]);exposure=int(a[6])
        if count<1 or exposure<1 or a[4].lower() not in ('g','r','i','z'):raise ValueError('Invalid exposure')
        result.append({'field_id':a[0],'ra_deg':float(p.ra.deg),'dec_deg':float(p.dec.deg),
            'filter':a[4].lower(),'count':count,'exposure_s':exposure,'block':block,'line':line})
    return result


def command(field, band, count, exposure):
    p=coord(field['ra_deg'],field['dec_deg'])
    return '\t'.join([field['field_id'],p.ra.to_string(unit=u.hour,sep=':',precision=3,pad=True),
        p.dec.to_string(unit=u.deg,sep=':',precision=2,pad=True,alwayssign=True),
        '0',band.lower(),str(count),str(exposure),'0'])+'\n'


def movement(a,b):
    # Identical coordinates need no mount movement even if target/filter names differ.
    if abs((a['ra_deg']-b['ra_deg']+180)%360-180)<1e-7 and abs(a['dec_deg']-b['dec_deg'])<1e-7:return 0.
    return slew({**a,'field_id':'previous'},{**b,'field_id':'next'})


def replay(commands, start, previous=None, band=None):
    rows=[];t=start
    for line_index,c in enumerate(commands):
        for repeat in range(c['count']):
            turn=90. if previous is None else movement(previous,c)
            change=10. if band is None or band!=c['filter'] else 0.
            begin=t+turn+change;end=begin+c['exposure_s']
            rows.append({k:c[k] for k in ('block','field_id','ra_deg','dec_deg','filter','exposure_s') } |
                {'command_index':line_index+1,'repeat':repeat+1,'overhead_start_unix':t,'slew_s':turn,
                 'filter_change_s':change,'readout_s':READOUT_S,'start_unix':begin,'mid_unix':begin+c['exposure_s']/2,
                 'end_unix':end,'release_unix':end+READOUT_S})
            t=end+READOUT_S;previous=c;band=c['filter']
    return rows


def validate(rows, twilight, sky):
    grid,_=load_grid();mp={f['field_id']:f for f in grid};summary=[]
    for block in ('EV','M31','TROJAN','YSO_G','YSO_R','AM'):
        rr=[r for r in rows if r['block']==block]
        if not rr:raise ValueError('Empty block '+block)
        keys=list(dict.fromkeys(r['field_id'] for r in rr));fields=[next(r for r in rr if r['field_id']==k) for k in keys]
        index={k:i for i,k in enumerate(keys)};times=[];ix=[]
        for r in rr:
            tt=np.linspace(r['start_unix'],r['end_unix'],max(2,int(np.ceil(r['exposure_s']/5))+1))
            times.extend(tt);ix.extend([index[r['field_id']]]*len(tt))
        m=sky.metrics([f['ra_deg'] for f in fields],[f['dec_deg'] for f in fields],times)
        radii=[]
        for f in fields:
            if block in ('EV','AM'):radii.append(mp[f['field_id'].replace('MP_','MP')]['radius_deg'])
            else:
                fp=footprint(f['ra_deg'],f['dec_deg']);radii.append(float(fp[-1].separation(fp).deg.max()))
        ix=np.asarray(ix);jj=np.arange(len(times));radius=np.asarray(radii)[ix]
        altitude=m['alt'][ix,jj]-radius
        min_alt=25 if block in ('EV','AM') else 30
        sun_limit=-10 if block in ('EV','AM') else -18
        if altitude.min()<=min_alt+.011:raise ValueError((block,'altitude',float(altitude.min())))
        if m['sun_alt'].max()>sun_limit+.003:raise ValueError((block,'Sun',float(m['sun_alt'].max())))
        moon_up=m['moon_alt']>-.011
        if np.any(moon_up & (m['moon_sep'][ix,jj]-radius < np.where(m['moon_fraction']>=.75,40,30)+.011)):
            raise ValueError((block,'Moon'))
        if block in ('EV','AM'):
            ref=next(b for b in twilight['blocks'] if b['side']==block)
            if rr[0]['start_unix']<ref['reserved_window_unix'][0] or rr[-1]['release_unix']>ref['reserved_window_unix'][1]:
                raise ValueError('Twilight boundary')
            assert len(rr)==len(ref['rows'])
            assert all(abs(r['mid_unix']-old['mid_unix'])<.5 and r['exposure_s']==old['exposure_s']
                       for r,old in zip(rr,ref['rows']))
        summary.append({'block':block,'frames':len(rr),'filter':sorted({r['filter'] for r in rr}),
            'overhead_start_local':local_time(rr[0]['overhead_start_unix']),
            'first_exposure_local':local_time(rr[0]['start_unix']),
            'readout_complete_local':local_time(rr[-1]['release_unix']),
            'elapsed_s':rr[-1]['release_unix']-rr[0]['overhead_start_unix'],
            'exposure_s':sum(r['exposure_s'] for r in rr),
            'center_alt_min_deg':float(m['alt'][ix,jj].min()),'center_alt_max_deg':float(m['alt'][ix,jj].max()),
            'footprint_alt_min_deg':float(altitude.min()),'moon_sep_min_deg':float(m['moon_sep'][ix,jj].min()),
            'sun_alt_max_deg':float(m['sun_alt'].max()),'verification':'passed: exposure-wide <=5-second samples + footprint cap'})
    for a,b in zip(rows,rows[1:]):
        if abs(a['release_unix']-b['overhead_start_unix'])>.001:raise ValueError('Timeline gap')
    return summary


def science(rows, root):
    rr=[r for r in rows if r['block']=='TROJAN'];times=[r['mid_unix'] for r in rr];targets=[]
    for number in ('17365','4722'):
        path=BASE/'cache'/f'horizons_{number}_2026-11-03_2026-11-17_1h.json'
        track=json.loads(path.read_text());ra,dec=interpolate(track,times)
        retained=in_field(ra,dec,TJ['ra_deg'],TJ['dec_deg'],3,.15)
        if not np.all(retained):raise ValueError('Trojan not retained: '+number)
        v=np.interp(Time(times,format='unix').jd,track['datetime_jd'],track['V'])
        targets.append({'object':number,'name':track['targetname'],'V_min':float(v.min()),'V_max':float(v.max()),
            'ra_mid_deg':float(np.median(ra)),'dec_mid_deg':float(np.median(dec)),
            'all_exposures_inside_with_015deg_margin':True,'source':track['source'],
            'ephemeris_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    stars=[]
    for s in load_bright_stars():
        sep=float(coord(s['ra_deg'],s['dec_deg']).separation(coord(TJ['ra_deg'],TJ['dec_deg'])).deg)
        if sep<4.122:stars.append({**s,'separation_center_deg':sep,
            'inside_footprint':bool(in_field(s['ra_deg'],s['dec_deg'],TJ['ra_deg'],TJ['dec_deg'],3,0))})
    result={'field':TJ,'targets':targets,'bright_stars':sorted(stars,key=lambda x:x['mag']),
        'purpose':'Representative high-amplitude Trojan photometry test; one night is not a new period measurement.',
        'paper':'https://doi.org/10.1051/0004-6361/202347158',
        'brightness_note':'Horizons V is a prediction, not measured R; external bright-star halos remain image-dependent.'}
    (root/'trojan_selection.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def visibility_report(rows, sky, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from .plot_style import apply_style,save_pair
    times=np.linspace(rows[0]['overhead_start_unix'],rows[-1]['release_unix'],800)
    objects=[('M31',10.684708333333333,41.26875),('Trojan: Thymbraeus field',TJ['ra_deg'],TJ['dec_deg']),
             ('YSO 84 / -4.8',84.,-4.8),('YSO 314 / 44.1',314.,44.1)]
    m=sky.metrics([o[1] for o in objects],[o[2] for o in objects],times)
    apply_style();fig,axes=plt.subplots(2,1,figsize=(14,8),layout='constrained',height_ratios=[2,1])
    dt=lambda tt:[datetime.fromtimestamp(t,TZ) for t in tt]
    for i,(name,ra,dec) in enumerate(objects):axes[0].plot(dt(times),m['alt'][i],label=name,lw=1.8)
    axes[0].axhline(30,color='red',ls='--',lw=1,label='30 deg field-centre reference')
    axes[0].set(ylabel='Field centre altitude (deg)',ylim=(-10,92),title='2026-11-11 night | computed at Xinglong')
    axes[0].legend(fontsize=10,ncol=2)
    names=['EV','M31','TROJAN','YSO_G','YSO_R','AM'];colors=['#1f77b4','#9467bd','#2ca02c','#17becf','#d62728','#ff8c00']
    for i,(name,color) in enumerate(zip(names,colors)):
        rr=[r for r in rows if r['block']==name]
        a,b=rr[0]['overhead_start_unix'],rr[-1]['release_unix']
        axes[1].barh(i,(b-a)/86400,left=mdates.date2num(dt([a]))[0],color=color,height=.6)
        axes[1].text(mdates.date2num(dt([(a+b)/2]))[0],i,f"{len(rr)} frames",ha='center',va='center',fontsize=9)
        axes[0].axvspan(*dt([a,b]),color=color,alpha=.045)
    axes[1].set(yticks=range(6),yticklabels=names,xlabel='Beijing time UTC+8 (morning is 12 Nov)',title='One sequential file, no WAIT instructions')
    for ax in axes:
        ax.set_xlim(*dt([times[0],times[-1]]));ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M',tz=TZ));ax.grid(alpha=.15)
    save_pair(fig,output/'full_night_overview','Fixed full-night sequence and YSO candidate visibility; nominal predictions')
    plt.close(fig)
    comparison=[]
    m31_end=max(r['release_unix'] for r in rows if r['block']=='M31')
    for name,ra,dec in objects[2:]:
        at_end=sky.metrics([ra],[dec],[m31_end])['alt'][0,0]
        # Whole square bound, independent of selected schedule, at five-second steps.
        t=np.arange(m31_end,rows[-1]['release_unix'],5)
        alt=sky.metrics([ra],[dec],t)['alt'][0]-2.122
        valid=t[alt>30.011]
        comparison.append({'field':name,'center_alt_at_M31_end_deg':float(at_end),
            'after_M31_footprint_above30_first':local_time(valid[0]) if len(valid) else None,
            'after_M31_footprint_above30_last':local_time(valid[-1]) if len(valid) else None})
    return comparison


def build(twilight_dir,m31_file,output):
    output.mkdir(parents=True,exist_ok=True)
    twilight=json.loads((twilight_dir/'plan.json').read_text())
    setup_iers(BASE/'cache'/twilight['config']['iers_file']);sky=make_sky()
    ev=lowercase_filters((twilight_dir/'EV.txt').read_text())
    am=lowercase_filters((twilight_dir/'AM.txt').read_text())
    m31=lowercase_filters(m31_file.read_text())
    if twilight['config']['readout_s'] != READOUT_S:
        raise ValueError('Regenerate twilight with the current readout setting')
    first=twilight['blocks'][0]['window_unix'][0]
    initial=read_lines(ev,'EV')+read_lines(m31,'M31')
    prefix=replay(initial,first);last=prefix[-1]
    tj=read_lines(command(TJ,'R',1,120),'TROJAN')[0]
    yg=read_lines(command(YSO,'G',YSO_FRAMES_PER_BAND,60),'YSO_G')[0]
    yr=read_lines(command(YSO,'R',YSO_FRAMES_PER_BAND,60),'YSO_R')[0]
    am_commands=read_lines(am,'AM')
    ref_am=next(b for b in twilight['blocks'] if b['side']=='AM')
    # Preserve the accepted AM shutter times; its real long handover starts earlier.
    desired_yso_end=ref_am['rows'][0]['start_unix']-movement(yr,am_commands[0])
    yso_first=desired_yso_end-(2*YSO_FRAMES_PER_BAND*(60+READOUT_S)+10)
    tj_end=yso_first-movement(tj,yg)-10
    tj_first=last['release_unix']+movement(last,tj)+10
    duration=tj_end-tj_first;n=int(duration//(120+READOUT_S));remaining=duration-n*(120+READOUT_S)
    final_exposure=int(np.floor(remaining-READOUT_S))
    if n<1 or not 5<=final_exposure<=120:raise ValueError('No valid integer Trojan tail exposure')
    tj_text=command(TJ,'R',n,120)+command(TJ,'R',1,final_exposure)
    yso_text=command(YSO,'G',YSO_FRAMES_PER_BAND,60)+command(YSO,'R',YSO_FRAMES_PER_BAND,60)
    commands=initial+read_lines(tj_text,'TROJAN')+[yg,yr]+am_commands
    rows=replay(commands,first);summary=validate(rows,twilight,sky)
    # Preserve all M31/twilight columns except the user-corrected lowercase filter token.
    full=''.join(text if text.endswith('\n') else text+'\n' for text in (ev,m31,tj_text,yso_text,am))
    file=output/'FULL_NIGHT_20261111.txt';file.write_text(full)
    check=read_lines(full,'ALL');expanded=replay(check,first)
    assert len(expanded)==len(rows)
    assert all(abs(a['start_unix']-b['start_unix'])<.001 for a,b in zip(rows,expanded))
    assert m31 in full and full.startswith(ev) and full.endswith(am)
    old_copy=output/'02_M31_ORIGINAL.txt'
    if old_copy.exists():old_copy.unlink()  # Obsolete generated uppercase controller copy.
    for name,text in [('01_EV.txt',ev),('02_M31.txt',m31),('03_TROJAN.txt',tj_text),('04_YSO_GR.txt',yso_text),('05_AM.txt',am)]:
        (output/name).write_text(text)
    write_csv(output/'exposures.csv',[{**r,'start_local':local_time(r['start_unix']),
        'mid_local':local_time(r['mid_unix']),'end_local':local_time(r['end_unix'])} for r in rows])
    write_csv(output/'block_summary.csv',summary)
    selected=science(rows,output);comparison=visibility_report(rows,sky,output)
    result={'date':'2026-11-11','start_local':local_time(first),'finish_local':local_time(rows[-1]['release_unix']),
        'controller_lines':len(check),'expanded_exposures':len(rows),'all_geometry_checks_passed':True,
        'blocks':summary,'yso':{'field':YSO,'G_frames':YSO_FRAMES_PER_BAND,'R_frames':YSO_FRAMES_PER_BAND,'each_exposure_s':60,
        'each_band_exposure_plus_readout_s':YSO_FRAMES_PER_BAND*(60+READOUT_S),'between_band_filter_change_s':10,'order':'G then R near transit'},
        'trojan':{'field':TJ,'R120_frames':n,'final_R_exposure_s':final_exposure,'targets':selected['targets']},
        'YSO_candidate_comparison':comparison,
        'AM_first_exposure_shift_s':next(r['start_unix'] for r in rows if r['block']=='AM')-ref_am['rows'][0]['start_unix'],
        'approved_twilight_preserved_except_filter_case':True,'controller_filter_case':'lowercase',
        'yso_field_registry':YSO_FIELDS,'M31_preserved_except_filter_case':True,'M31_input_sha256':hashlib.sha256(m31_file.read_bytes()).hexdigest(),
        'twilight_plan_sha256':hashlib.sha256((twilight_dir/'plan.json').read_bytes()).hexdigest(),
        'assumptions':['Same fixed position: no extra slew between repeats or filter-only changes.',
            f'Nominal readout={READOUT_S:g}s, filter change=10s, measured two-axis slew formula.',
            'No realtime clock gate or automatic cutoff in controller TXT.',
            'Start the combined list once at the stated evening boundary; do not execute the separate copies again.',
            'No server operation or annual rerun.']}
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    readme=(f"SHARP full-night test 2026-11-11\nStart FULL_NIGHT_20261111.txt at {result['start_local']} once.\n"
        "Sequence: EV -> M31 griz -> TJN_17365 -> YSO_ORI g then r -> AM. Filter tokens are lowercase.\n"
        "M31 coordinates/counts/exposures are unchanged; only filter tokens have been lowercased.\n"
        "The numbered files are reference copies, NOT additional jobs after the combined list.\n"
        f"YSO: {YSO_FRAMES_PER_BAND}x60s per band, readout {READOUT_S:g}s; filter switch adds 10s.\n"
        f"Trojan: {n}x120s + one {final_exposure}s bridging exposure in R, same fixed field.\n"
        "AM long slew begins just before the original boundary to preserve accepted shutter times.\n"
        "All times are nominal; static TXT has no WAIT or automatic deadline enforcement.\n")
    (output/'READ_ME.txt').write_text(readme)
    files=[p for p in sorted(output.iterdir()) if p.is_file() and p.suffix!='.zip' and p.name!='manifest.json']
    manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(output/'SHARP_FULL_NIGHT_20261111.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in files+[output/'manifest.json']:z.write(p,p.name)
    print(json.dumps(result,indent=2))


def main():
    p=argparse.ArgumentParser();p.add_argument('--twilight',required=True);p.add_argument('--m31',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();build(Path(a.twilight),Path(a.m31),Path(a.output))


if __name__=='__main__':main()
