"""Auditable products and guarded controller draft export."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import astropy.units as u
from .geometry import coord,unix,local_time,footprint
from .catalog import write_csv


def summarize(plan):
    result=[]
    for b in plan['blocks']:
        for f in b['fields']:
            rr=[r for r in b['rows'] if r['field_id']==f['field_id']]
            mids=np.array([unix(r['mid_utc']) for r in rr])
            result.append({'block':b['block_id'],'field_id':f['field_id'],'ra_deg':f['ra_deg'],'dec_deg':f['dec_deg'],
                           'points':len(rr),'cadence_median_s':float(np.median(np.diff(mids))) if len(mids)>1 else None,
                           'cadence_min_s':float(np.min(np.diff(mids))) if len(mids)>1 else None,
                           'cadence_max_s':float(np.max(np.diff(mids))) if len(mids)>1 else None,
                           'cadence_intervals_s':json.dumps(np.diff(mids).round(3).tolist()),
                           'arc_span_s':float(mids[-1]-mids[0]),'core_targets':','.join(f['core_targets']),
                           'start_local':local_time(unix(rr[0]['start_utc'])),'end_local':local_time(unix(rr[-1]['end_utc']))})
    return result


def write_products(plan,validation,path):
    p=Path(path);p.mkdir(parents=True,exist_ok=True)
    for name,data in [('plan.json',plan),('validation.json',validation)]:
        (p/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    rows=[r for b in plan['blocks'] for r in b['rows']]
    if rows:
        write_csv(p/'exposures.csv',rows);write_csv(p/'fields.csv',summarize(plan))
    else:
        (p/'exposures.csv').write_text('');(p/'fields.csv').write_text('')
    # Reusing an output directory with a disabled/skipped block must not retain its old list.
    present={b['block_id'] for b in plan['blocks']}
    for name in {'EV_A','EV_B','TJ','AM_B','AM_A'}-present:
        (p/(name+'.csv')).unlink(missing_ok=True)
    manifest={'date':plan['date'],'status':plan['status'],'control_export':'disabled: verified filter tokens and timing integration required',
              'plan_sha256':hashlib.sha256((p/'plan.json').read_bytes()).hexdigest(),
              'blocks':[{k:v for k,v in b.items() if k not in ['rows','quality']} for b in plan['blocks']],
              'skips':plan['skips'],'free_intervals':plan['free_intervals'],'free_intervals_scope':plan['free_intervals_scope'],'validation_passed':validation['passed']}
    for b in plan['blocks']:write_csv(p/(b['block_id']+'.csv'),b['rows'])
    (p/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


def export_controller_drafts(plan,path,token_evidence=None):
    """Only block drafts; cannot certify hardware execution or imply TXT enforces times."""
    mapping=plan['config']['controller_filter_tokens']
    filters={r['filter'] for b in plan['blocks'] for r in b['rows']}
    if not token_evidence or any(not mapping.get(f) for f in filters):
        raise ValueError('Verified filter token mapping and evidence required; no fallback to W')
    if any(any(c.isspace() for c in mapping[f]) for f in filters):raise ValueError('Invalid filter token')
    p=Path(path);p.mkdir(parents=True,exist_ok=True)
    for b in plan['blocks']:
        lines=[]
        for r in b['rows']:
            c=coord(r['ra_deg'],r['dec_deg'])
            lines.append('\t'.join([r['field_id'],c.ra.to_string(unit=u.hour,sep=':',precision=3,pad=True),
                c.dec.to_string(unit=u.deg,sep=':',precision=2,alwayssign=True,pad=True),'0',mapping[r['filter']],'1',str(int(r['exposure_s'])),'0']))
        (p/(b['block_id']+'.draft.txt')).write_text('\n'.join(lines)+'\n')
    (p/'NOT_EXECUTABLE.json').write_text(json.dumps({'evidence':token_evidence,
        'required':'Timed dispatcher must honor every row start_utc; TXT has no waits. Revalidate at actual times; abort if deadline/geometry cannot be met.'},indent=2)+'\n')


def plot_plan(plan,campaign,path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as md
    from datetime import datetime,timezone
    from zoneinfo import ZoneInfo
    p=Path(path)
    fig=plt.figure(figsize=(13,8),constrained_layout=True)
    grid=fig.add_gridspec(2,3,height_ratios=[1,1.4])
    timeline=fig.add_subplot(grid[0,:]);maps={k:fig.add_subplot(grid[1,i]) for i,k in enumerate(['EV','TJ','AM'])}
    colors={'EV_A':'#df8f44','EV_B':'#bd6330','TJ':'#4a72b0','AM_B':'#279b91','AM_A':'#61b8aa'}
    for n,b in enumerate(plan['blocks']):
        color=colors[b['block_id']]
        for r in b['rows']:
            a=md.date2num(datetime.fromtimestamp(unix(r['start_utc']),timezone.utc));dur=r['exposure_s']/86400
            timeline.broken_barh([(a,dur)],(n-.3,.6),facecolors=color)
        key=b['block_id'][:2];ax=maps[key]
        for j,f in enumerate(b['fields']):
            fp=footprint(f['ra_deg'],f['dec_deg']);ras=fp.ra.deg[:8];decs=fp.dec.deg[:8]
            ras=f['ra_deg']+(ras-f['ra_deg']+180)%360-180
            ax.plot(np.r_[ras,ras[0]],np.r_[decs,decs[0]],color=color,lw=1.3)
            label=('T' if key=='TJ' else b['block_id'][-1])+str(j+1)
            ax.annotate(label,(f['ra_deg'],f['dec_deg']+1.3),fontsize=8,ha='center',color=color)
            if campaign:
                for k in f['core_targets']:
                    t=campaign['tracks'][k];ax.plot(t['RA'],t['DEC'],color='black',lw=.8)
                    ax.scatter([t['RA'][0]],[t['DEC'][0]],s=10,color='black')
                    ax.annotate(k,(t['RA'][0],t['DEC'][0]),xytext=(5,3),textcoords='offset points',fontsize=8)
    timeline.set_yticks(range(len(plan['blocks'])),[b['block_id'] for b in plan['blocks']])
    timeline.xaxis.set_major_locator(md.HourLocator(interval=2,tz=ZoneInfo('Asia/Shanghai')))
    timeline.xaxis.set_major_formatter(md.DateFormatter('%H:%M',tz=ZoneInfo('Asia/Shanghai')))
    timeline.set_xlabel('Beijing time (morning belongs to the next day)');timeline.set_title(plan['date']+' | LOCAL SIMULATION | bars = planned shutter-open intervals')
    titles={'EV':'Evening fields','TJ':'Trojan fields / 14-day tracks','AM':'Morning fields'}
    for key,ax in maps.items():
        ax.set(xlabel='ICRF RA (deg)',ylabel='Dec (deg)',title=titles[key]);ax.invert_xaxis();ax.grid(alpha=.2)
        if not ax.lines:ax.text(.5,.5,'No eligible block',transform=ax.transAxes,ha='center')
        else:ax.margins(.15)
    fig.savefig(p/'plan.png',dpi=150);fig.savefig(p/'plan.pdf');plt.close(fig)
