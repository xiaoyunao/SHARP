"""Shared paths, atomic writes and advisory locks for the deployed v2."""
from contextlib import contextmanager
from datetime import datetime, timedelta, time
import fcntl
import json
import os
from pathlib import Path
import tempfile
from zoneinfo import ZoneInfo

TZ = ZoneInfo('Asia/Shanghai')
BASE = Path(__file__).resolve().parent


def load_config(path=None):
    value = json.loads((BASE/'server_defaults.json').read_text())
    if path is None:
        candidate=Path(value['base'])/'v2/config.json'
        path=os.environ.get('SMT_V2_CONFIG') or (candidate if candidate.is_file() else None)
    if path:
        value.update(json.loads(Path(path).read_text()))
    if value.get('unknown_auto_submit'):
        raise ValueError('Unknown submission is a website-only operation in v2')
    return value


def runtime(config):
    return Path(config['base'])/'v2/runtime'


def observing_night(when):
    when = when.astimezone(TZ)
    if when.hour < 12:
        when -= timedelta(days=1)
    return when.strftime('%Y%m%d')


def observation_start(config, night):
    """Earliest published execution start, including initial setup."""
    path = runtime(config)/'published'/f'{night}.json'
    if not path.exists():
        return None
    metadata = json.loads(path.read_text())
    plan = json.loads((Path(metadata['revision'])/'plan.json').read_text())
    values = [datetime.fromisoformat(plan[key]).astimezone(TZ)
              for key in ('start_local', 'window_start_local') if plan.get(key)]
    return min(values) if values else None


def followup_cutoff(config, night):
    day = datetime.strptime(night, '%Y%m%d').date()
    deadline = datetime.combine(day, time.fromisoformat(config.get('followup_cutoff_local', '17:00')), TZ)
    start = observation_start(config, night)
    return min(deadline, start) if start else deadline


def atomic_write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = value if isinstance(value, bytes) else value.encode()
    fd, temporary = tempfile.mkstemp(prefix='.'+path.name+'.', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as stream:
            os.fchmod(stream.fileno(),(path.stat().st_mode & 0o777) if path.exists() else 0o644)
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_json(path, value):
    atomic_write(path, json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


@contextmanager
def lock(path, blocking=True):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
        yield


def utc_now():
    return datetime.now(TZ).isoformat(timespec='seconds')
