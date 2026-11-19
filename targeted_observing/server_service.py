"""Noon scheduling, durable follow-up events and independent L2 processing ticks."""
import argparse
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path
import subprocess
import traceback

from .server_common import (TZ, load_config, runtime, write_json, lock, utc_now,
                            observing_night, followup_cutoff, observation_start)
from .server_followup import queue_snapshot, route, submit


def plan(night,config,reason,finalize=False):
    from .server_scheduler import build
    from .server_publish import publish
    with lock(runtime(config)/'locks/references.lock',shared=True),lock(runtime(config)/'locks'/f'plan_{night}.lock'):
        now=datetime.now(TZ)
        start=observation_start(config,night)
        if now >= followup_cutoff(config,night):
            # Only drain commands accepted before the cutoff, while the telescope
            # has not started. No history-triggered or in-night replanning.
            if not finalize or start is None or now >= start:
                raise ValueError('Night is frozen: follow-up acceptance/execution has closed')
        with lock(runtime(config)/'queue.lock'):
            requests,errors=queue_snapshot(config,night)
        day=datetime.strptime(night,'%Y%m%d').date().isoformat()
        options={}
        if start is not None and now<start:
            metadata=json.loads((runtime(config)/'published'/f'{night}.json').read_text())
            previous=json.loads((Path(metadata['revision'])/'plan.json').read_text())
            window=previous.get('window_start_local')
            if window and now>datetime.fromisoformat(window):
                # A shorter EV can start after its reserved window begins. The
                # telescope is still idle: rebuild only future commands, with
                # generation/transfer lead time, and no execution prefix.
                options['now']=(now+timedelta(seconds=config['late_replan_lead_s'])).timestamp()
        value=build(day,config,requests,**options)
        if datetime.now(TZ) >= datetime.fromisoformat(value['start_local']):
            raise ValueError('Execution start reached while planning; publication blocked')
        value['queue_errors']=errors
        result=publish(value,config,reason)
        write_json(runtime(config)/'last_plan_result.json',result)
        return result


def enqueue_processing(night,config):
    with lock(runtime(config)/'locks'/f'enqueue_{night}.lock'):
        return _enqueue_processing(night,config)


def _enqueue_processing(night,config):
    from .server_l2_watch import readiness
    ready=readiness(night,config)
    write_json(runtime(config)/'l2_watch'/f'{night}.json',dict(**ready,checked_at=utc_now()))
    if not ready['ready']:return ready
    l2=Path(config['processed_root'])/night/'L2'
    if not l2.is_dir():return dict(status='no_l2',night=night)
    files=[p for p in l2.iterdir() if p.is_file() and p.name.endswith(('.fits','.fits.gz'))]
    if not files:return dict(status='no_science_files',night=night)
    now=datetime.now(TZ).timestamp()
    qpath=runtime(config)/'processing_queue'/f'{night}.json'
    previous=json.loads(qpath.read_text()) if qpath.exists() else {}
    signature=hashlib.sha256(str([(p.name,p.stat().st_size,p.stat().st_mtime_ns,p.stat().st_ctime_ns) for p in sorted(files)]).encode()).hexdigest()
    process_path=runtime(config)/'processing'/f'{night}.json'
    process_state=json.loads(process_path.read_text()) if process_path.exists() else {}
    if previous.get('signature')==signature:
        if process_state.get('status') in ('complete','no_valid_science_inputs'):return previous
        age=now-datetime.fromisoformat(previous['queued_at']).timestamp()
        if age<900:return previous
        if previous.get('attempts',1)>=3:
            return {**previous,'status':'needs_attention_after_three_attempts'}
    name='sharp_v2_'+night
    active=subprocess.check_output(['squeue','-h','-u','smtpipeline','-n',name,'-o','%A'],text=True).strip()
    if active:return dict(status='already_running',job=active,night=night)
    logs=runtime(config)/'logs';logs.mkdir(parents=True,exist_ok=True)
    command=['sbatch','--parsable','--job-name',name,'--cpus-per-task',str(config['workers']),
             '--mem','32G','--output',str(logs/(night+'_process_%j.log')),
             str(Path(config['base'])/'targeted_observing/server_process_job.sh'),night]
    if previous.get('signature')==signature and process_state.get('search_complete') and process_state.get('known_complete'):
        command.append('--resume-review')
    job=subprocess.check_output(command,text=True).strip().split(';')[0]
    result=dict(status='queued',night=night,job_id=job,signature=signature,queued_at=utc_now(),
                attempts=previous.get('attempts',0)+1 if previous.get('signature')==signature else 1)
    write_json(qpath,result);return result


def discover_events(night,config):
    directory=runtime(config)/'events'/night;directory.mkdir(parents=True,exist_ok=True)
    known=set()
    for p in directory.glob('*.json'):
        if p.name.endswith('.result.json'):continue
        try:known.add(json.loads(p.read_text())['request_sha256'])
        except Exception:continue
    requests,errors=queue_snapshot(config,night)
    for row in requests:
        if row['sha256'] in known:continue
        # Direct SFTP requests use server inode arrival/change time, never a
        # client-supplied submitted_at. The CLI has already created its event.
        source=runtime(config)/'followups'/night/(row['candidate_id']+'.json')
        arrival=datetime.fromtimestamp(source.stat().st_ctime,TZ)
        try:
            resolved=route(row,config,arrival)
            if resolved['night'] != night:
                result=submit(row,config,arrival)
                with lock(runtime(config)/'queue.lock'):
                    write_json(runtime(config)/'deferred_followups'/night/source.name,
                               dict(original=json.loads(source.read_text()),resolved=result))
                    source.unlink()
                continue
        except ValueError as exc:
            errors.append(dict(file=str(source),error=str(exc)))
            continue
        stamp=datetime.now(TZ).strftime('%Y%m%dT%H%M%S%f')
        write_json(directory/(stamp+'_'+row['candidate_id']+'.json'),dict(
            candidate_id=row['candidate_id'],request_sha256=row['sha256'],detected_at=utc_now()))
    write_json(runtime(config)/'followup_queue_errors'/f'{night}.json',errors)
    return [p for p in sorted(directory.glob('*.json')) if not p.name.endswith('.result.json')
            and not p.with_suffix('.result.json').exists()]


def daily_run(config,night=None):
    now=datetime.now(TZ);night=night or now.strftime('%Y%m%d')
    from .server_history import ingest
    ingest(config)
    result=plan(night,config,'noon')
    # A noon generation includes every pre-noon candidate in one complete snapshot.
    for event in discover_events(night,config):
        included={r.get('request_sha256') for r in result['followups']}
        if json.loads(event.read_text())['request_sha256'] in included:
            write_json(event.with_suffix('.result.json'),dict(status='included_in_noon_run',revision=result['revision']))
    return dict(plan=result)


def daily(config,night=None):
    night=night or datetime.now(TZ).strftime('%Y%m%d')
    with lock(runtime(config)/'locks'/f'daily_{night}.lock'):
        published=runtime(config)/'published'/f'{night}.json'
        if published.exists():
            return dict(plan=json.loads(published.read_text()))
        return daily_run(config,night)


def tick(config):
    try:
        with lock(runtime(config)/'locks/tick.lock',blocking=False):
            now=datetime.now(TZ);night=now.strftime('%Y%m%d')
            from .server_l2_watch import recent_nights, night_finished
            # Processing is independent of noon scheduling and its failures.
            for finished in recent_nights(config,now):
                try:
                    if not night_finished(finished,config,now):
                        continue
                    enqueue_processing(finished,config)
                except Exception as exc:
                    write_json(runtime(config)/'l2_watch_errors'/f'{finished}.json',
                               dict(error=str(exc),last_attempt=utc_now()))
                    traceback.print_exc()
            from .script_retention import cleanup
            cleanup(config,now,once_per_day=True)
            if (now.hour>=12 and now<followup_cutoff(config,night)
                    and not (runtime(config)/'published'/f'{night}.json').exists()):
                daily(config,night)
            # After midnight, discover late raw deliveries to yesterday as well.
            active=observing_night(now)
            if active!=night:discover_events(active,config)
            events=discover_events(night,config)
            for event in events if now.hour>=12 else []:
                try:
                    command=json.loads(event.read_text())
                    closed=now>=followup_cutoff(config,night)
                    if closed:
                        published=runtime(config)/'published'/f'{night}.json'
                        metadata=json.loads(published.read_text()) if published.exists() else {}
                        if any(r.get('request_sha256')==command['request_sha256'] for r in metadata.get('followups',[])):
                            write_json(event.with_suffix('.result.json'),dict(status='already_in_published_plan',revision=metadata['revision']))
                            continue
                        if command.get('reason')=='new_actual_known_detections':
                            write_json(event.with_suffix('.result.json'),dict(status='history_saved_for_next_noon'))
                            continue
                        start=observation_start(config,night)
                        if start is None or now>=start:
                            source=runtime(config)/'followups'/night/(command['candidate_id']+'.json')
                            result=submit(json.loads(source.read_text()),config,now) if source.exists() else {}
                            write_json(event.with_suffix('.result.json'),dict(status='deferred_before_execution_freeze',resolved=result))
                            continue
                    result=plan(night,config,'event_'+event.stem,finalize=closed)
                    write_json(event.with_suffix('.result.json'),dict(status='replanned',revision=result['revision']))
                except Exception as exc:
                    write_json(runtime(config)/'event_errors'/event.name,dict(error=str(exc),last_attempt=utc_now()))
                    traceback.print_exc()
                    # Leave the event pending for automatic retry; never drop it.
                    break
    except BlockingIOError:
        return dict(status='another_tick_running')
    return dict(status='tick_complete')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['daily','plan','tick','process'])
    p.add_argument('--night');p.add_argument('--config');a=p.parse_args();c=load_config(a.config)
    if a.action=='daily':result=daily(c,a.night)
    elif a.action=='plan':result=plan(a.night or datetime.now(TZ).strftime('%Y%m%d'),c,'manual')
    elif a.action=='tick':result=tick(c)
    else:
        from .server_pipeline import process
        result=process(a.night,c)
    print(json.dumps(result,indent=2),flush=True)
