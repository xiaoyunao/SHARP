"""One-candidate JSON submission contract; daily queue, no legacy candidate mining."""
import argparse
from datetime import datetime, time, timedelta
import hashlib
import json
from pathlib import Path
import re

from .server_common import TZ, load_config, runtime, write_json, lock, observing_night, followup_cutoff
from .shared_fields import SharedFields


def validate(document, catalog=None):
    if not isinstance(document, dict):
        raise ValueError('Request must be a JSON object')
    required = {'candidate_id', 'night', 'field_id', 'start_time', 'frames', 'exposure_s', 'filter'}
    if not required <= document.keys():
        raise ValueError('Missing fields: '+', '.join(sorted(required-document.keys())))
    result = {key: document[key] for key in required}
    if not isinstance(result['candidate_id'], str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', result['candidate_id']):
        raise ValueError('candidate_id must be a safe 1-80 character identifier')
    night = result['night']
    if not isinstance(night, str) or not re.fullmatch(r'[0-9]{8}', night):
        raise ValueError('night must be YYYYMMDD text, identifying the evening date')
    day = datetime.strptime(night, '%Y%m%d').date()
    field_id = result['field_id']
    if not isinstance(field_id, str) or not re.fullmatch(r'[0-9]{4}', field_id):
        raise ValueError('field_id must be four-character text including leading zeroes')
    if field_id not in (catalog or SharedFields()).by_id:
        raise ValueError('Unknown shared field_id')
    start = datetime.fromisoformat(result['start_time'])
    if start.tzinfo is None:
        raise ValueError('start_time must include an explicit timezone, normally +08:00')
    local = start.astimezone(TZ)
    # A request belongs to its observing night, not the calendar date after midnight.
    noon = datetime.combine(day, time(12), TZ)
    if not 0 <= (local-noon).total_seconds() < 86400:
        raise ValueError('start_time does not belong to the requested observing night')
    for key, limit in (('frames', 10000), ('exposure_s', 3600)):
        if type(result[key]) is not int or not 1 <= result[key] <= limit:
            raise ValueError(key+' must be a positive integer within '+str(limit))
    if result['filter'] not in ('g', 'r', 'i', 'z'):
        raise ValueError('filter must be lowercase g/r/i/z')
    result['start_time'] = local.isoformat()
    result['execution_name'] = 'FOLLOW_'+field_id
    return result


def queue_snapshot(config, night):
    requests, errors = [], []
    for path in sorted((runtime(config)/'followups'/night).glob('*.json')):
        try:
            request = validate(json.loads(path.read_text()))
            if request['night'] != night or request['candidate_id'] != path.stem:
                raise ValueError('Directory/name does not match request identity')
            request['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            requests.append(request)
        except Exception as exc:
            errors.append(dict(file=str(path), error=str(exc)))
    return requests, errors


def route(document, config, now):
    """Resolve a newly arrived command to an evening still accepting requests."""
    now = now.astimezone(TZ)
    result = validate(document)
    today = now.strftime('%Y%m%d')
    earliest = (now+timedelta(days=1)).strftime('%Y%m%d') if now >= followup_cutoff(config, today) else today
    if result['night'] not in (today, observing_night(now), earliest):
        raise ValueError('Command must target the current night or the next accepting evening')
    original_night, original_start = result['night'], result['start_time']
    if result['night'] < earliest:
        days = (datetime.strptime(earliest, '%Y%m%d')-datetime.strptime(result['night'], '%Y%m%d')).days
        result['night'] = earliest
        result['start_time'] = (datetime.fromisoformat(result['start_time'])+timedelta(days=days)).isoformat()
        result['deferred_from'] = dict(night=original_night, start_time=original_start,
                                       reason='acceptance_closed_for_requested_night')
    if datetime.fromisoformat(result['start_time']).timestamp() <= now.timestamp():
        raise ValueError('Requested exposure time is already past')
    return result


def submit(document, config, now=None):
    now = (now or datetime.now(TZ)).astimezone(TZ)
    result = route(document, config, now)
    night = result['night']
    path = runtime(config)/'followups'/night/(result['candidate_id']+'.json')
    with lock(runtime(config)/'queue.lock'):
        result['submitted_at'] = now.isoformat()
        if path.exists():
            old = json.loads(path.read_text())
            revisions = runtime(config)/'followup_revisions'/night/result['candidate_id']
            write_json(revisions/(hashlib.sha256(path.read_bytes()).hexdigest()+'.json'), old)
        write_json(path, result)
        # Every new/updated command gets its own event; a locked planner takes
        # a complete queue snapshot for each event, never only the last candidate.
        stamp = now.strftime('%Y%m%dT%H%M%S%f')
        event = runtime(config)/'events'/night/(stamp+'_'+result['candidate_id']+'.json')
        write_json(event, dict(candidate_id=result['candidate_id'], submitted_at=now.isoformat(),
                               request_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return dict(accepted=True, request=str(path), event=str(event),
                night=night, start_time=result['start_time'], deferred='deferred_from' in result,
                acceptance_cutoff=followup_cutoff(config,night).isoformat(),
                scheduling='noon run' if night != now.strftime('%Y%m%d') or now.hour < 12 else 'replan queued',
                status_path=str(runtime(config)/'followup_status'/night/(result['candidate_id']+'.json')))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('request', type=Path)
    p.add_argument('--config')
    args = p.parse_args()
    try:
        print(json.dumps(submit(json.loads(args.request.read_text()), load_config(args.config)), indent=2))
    except Exception as exc:
        print(json.dumps(dict(accepted=False,error=str(exc))))
        raise SystemExit(2)
