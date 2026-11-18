"""Versioned, locked publication at the existing list and plot paths."""
import csv
from datetime import datetime, timedelta
import hashlib
import io
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
from astropy.table import Table, vstack

from .catalog import BASE
from .fixed_fields import from_tan
from .full_night_test import command, read_lines, replay
from .geometry import local_time, isot
from .server_common import atomic_write, write_json, runtime, utc_now, TZ, observation_start
from .shared_fields import SharedFields


def controller_text(rows):
    result=[]
    for r in rows:
        if result and all(r[k]==result[-1][k] for k in ('field_id','filter','exposure_s','mode','candidate_id')):
            result[-1]['count']+=1
        else:result.append({**r,'count':1})
    return ''.join(command(r,r['filter'],r['count'],r['exposure_s']) for r in result)


def legacy_rows(rows):
    result=[];cycle=0;previous=None
    for row in rows:
        block=(row['mode'],row['candidate_id'])
        if block!=previous:cycle+=1;previous=block
        shared=row['field_id'].rsplit('_',1)[-1]
        result.append({**row,'ra':row['ra_deg'],'dec':row['dec_deg'],'exptime_s':row['exposure_s'],
            'start_utc':isot(row['start_unix']).removesuffix('Z'),
            'mid_utc':isot(row['mid_unix']).removesuffix('Z'),
            'end_utc':isot(row['end_unix']).removesuffix('Z'),
            'cycle_id':cycle,'mode':'followup' if row['mode']=='FOLLOW' else row['mode'].lower(),
            'followup_selected_field_id':shared if shared.isdigit() else row['field_id']})
    return result


def plots(plan,output):
    # Keep the established summary/cycle/GIF filenames and manifest schema used
    # by the dashboard. Add YSO's independent footprints for this plot only.
    from survey.visualize_nightly import generate_plan_visualizations
    from .geometry import setup_iers
    setup_iers(BASE/'cache'/plan['config']['iers_file'])
    catalog=SharedFields();table=catalog.table.copy();table['field_id']=table['field_id'].astype('U32')
    for name in {r['field_id'] for r in plan['rows'] if r['mode']=='YSO'}:
        r=next(r for r in plan['rows'] if r['field_id']==name);extra=table[:1].copy()
        extra['field_id']=[name];extra['center_ra']=[r['ra_deg']];extra['center_dec']=[r['dec_deg']]
        ra,dec=from_tan(r['ra_deg'],r['dec_deg'],[-1.48,1.48,1.48,-1.48],[-1.48,-1.48,1.48,1.48])
        for i in range(4):extra[f'corner_ra_{i+1}']=[ra[i]];extra[f'corner_dec_{i+1}']=[dec[i]]
        table=vstack([table,extra])
    history=Table({'field_id':np.asarray(table['field_id']),'exposure_count':np.zeros(len(table),dtype=int)})
    result=generate_plan_visualizations(legacy_rows(plan['rows']),history,table,Path(output),plan['night'])
    return result


def publish(plan,config,reason):
    night=plan['night'];root=Path(config['base'])
    # Every generation is a reviewable immutable revision before live aliases move.
    revision=runtime(config)/'plans'/night/(utc_now().replace(':','').replace('+','_')+'_'+reason)
    revision.mkdir(parents=True,exist_ok=False)
    write_json(revision/'plan.json',plan)
    text=controller_text(plan['rows']);atomic_write(revision/(night+'_plan.txt'),text)
    rows=legacy_rows(plan['rows']);write_json(revision/(night+'_plan.json'),rows)
    stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=list(rows[0])) if rows else None
    if writer:
        writer.writeheader();writer.writerows(rows)
    atomic_write(revision/'exposures.csv',stream.getvalue())
    # A static controller list cannot express waits. Export exact timed segments
    # whenever a safe, legal filler exposure cannot bridge a commitment.
    segments=[];group=[]
    for row in plan['rows']:
        if group and row['overhead_start_unix']-group[-1]['release_unix']>1:
            segments.append(group);group=[]
        group.append(row)
    if group:segments.append(group)
    starts=[]
    for i,group in enumerate(segments,1):
        filename=f'{night}_segment{i:02d}.txt'
        atomic_write(revision/filename,controller_text(group))
        starts.append(dict(list=filename,start_local=local_time(group[0]['overhead_start_unix']),
                           first_shutter_local=local_time(group[0]['start_unix']),frames=len(group)))
        parsed=read_lines(controller_text(group),'CHECK')
        previous=plan['rows'][plan['rows'].index(group[0])-1] if group[0] is not plan['rows'][0] else None
        reproduced=replay(parsed,group[0]['overhead_start_unix'],previous,
                          previous['filter'] if previous else None)
        if len(reproduced)!=len(group):raise ValueError('Controller frame count changed')
        if max(abs(a['start_unix']-b['start_unix']) for a,b in zip(reproduced,group))>10:
            raise ValueError('Static controller timeline differs by more than 10 seconds')
    write_json(revision/'START_TIMES.json',dict(segments=starts,requires_timed_segments=len(segments)>1,
        note='The controller TXT contains no WAIT. Start each segment at its stated time; do not append a new revision to an executing old plan.'))
    pending=plan['rows'][plan.get('protected_prefix_frames',0):]
    if plan.get('execution_update_required'):
        atomic_write(revision/(night+'_remaining.txt'),controller_text(pending))
        write_json(revision/'RESUME.json',dict(protected_prefix_frames=plan.get('protected_prefix_frames',0),
            remaining_list=night+'_remaining.txt',remaining_frames=len(pending),
            first_overhead_local=local_time(pending[0]['overhead_start_unix']) if pending else None,
            note='Use only the unexecuted suffix; honour segment start times in START_TIMES.json.'))
    visual=plots(plan,revision/'plots')
    # Rendering can outlast scheduling. Recheck immediately before changing any
    # public file, leaving a short budget for copying aliases and the marker.
    starts_at=[datetime.fromisoformat(plan['start_local'])]
    previous_start=observation_start(config,night)
    if previous_start:starts_at.append(previous_start)
    if datetime.now(TZ)+timedelta(seconds=30) >= min(starts_at):
        raise ValueError('Execution start too close/reached during rendering; publication blocked')
    liveplots=root/'survey/plots'/night
    liveplots.mkdir(parents=True,exist_ok=True)
    plotmeta={}
    for key,value in visual.items():
        if isinstance(value,list):plotmeta[key]=[str(liveplots/Path(v).name) for v in value]
        elif isinstance(value,str) and value:plotmeta[key]=str(liveplots/Path(value).name)
        else:plotmeta[key]=value
    for path in (revision/'plots').iterdir():
        if path.is_file():atomic_write(liveplots/path.name,path.read_bytes())
    plotmeta['night_tag']=night;write_json(liveplots/(night+'_plots.json'),plotmeta)
    expected={p.name for p in (revision/'plots').iterdir()}
    for path in liveplots.glob(night+'_cycle*.png'):
        if path.name not in expected:
            saved=revision/'previous_extra_plots';saved.mkdir(exist_ok=True)
            path.replace(saved/path.name)
    for directory in (root/'survey/runtime/plans',root/'script'):
        atomic_write(directory/(night+'_plan.txt'),text)
        write_json(directory/(night+'_plan.json'),rows)
        shutil.copy2(revision/'START_TIMES.json',directory/(night+'_START_TIMES.json'))
        for item in starts:atomic_write(directory/item['list'],(revision/item['list']).read_bytes())
        if plan.get('execution_update_required'):
            atomic_write(directory/(night+'_remaining.txt'),(revision/(night+'_remaining.txt')).read_bytes())
            atomic_write(directory/(night+'_RESUME.json'),(revision/'RESUME.json').read_bytes())
    atomic_write(root/'script/current_plan.txt',text)
    metadata=dict(night=night,revision=str(revision),reason=reason,published_at=utc_now(),
                  plan_sha256=hashlib.sha256(text.encode()).hexdigest(),exposures=len(rows),
                  requires_timed_segments=len(segments)>1,followups=plan['followups'],
                  execution_update_required=plan.get('execution_update_required',False),
                  verification=plan['verification'],version=config['version'])
    write_json(runtime(config)/'published'/f'{night}.json',metadata)
    for row in plan['followups']:write_json(runtime(config)/'followup_status'/night/(row['candidate_id']+'.json'),{**row,'revision':str(revision)})
    return metadata
