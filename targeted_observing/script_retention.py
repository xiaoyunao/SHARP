"""Prune only recognized dated controller products in the synchronization folder."""
from datetime import datetime, timedelta
from pathlib import Path
import re

from .server_common import TZ, runtime, write_json, lock


PRODUCT = re.compile(r'(\d{8})_(?:plan\.(?:txt|json)|plan_with_times\.txt|START_TIMES\.json|segment\d+\.txt|remaining\.txt|RESUME\.json)')
KEEP = re.compile(r'\d{8}_plan(?:_with_times)?\.txt')


def cleanup(config, now=None, once_per_day=False):
    now = (now or datetime.now(TZ)).astimezone(TZ)
    today = now.date()
    days = int(config.get('script_retention_nights', 7))
    if days < 1:
        raise ValueError('script_retention_nights must be positive')
    first = datetime.strptime(config.get('script_retention_start', '20261118'), '%Y%m%d').date()
    cutoff = max(first, today-timedelta(days=days-1)).strftime('%Y%m%d')
    root = Path(config['base'])/'script'
    record = runtime(config)/'script_cleanup'/(today.strftime('%Y%m%d')+'.json')
    with lock(runtime(config)/'locks/script_cleanup.lock'):
        if once_per_day and record.exists():
            return dict(status='already_checked_today')
        removed = []
        if root.is_dir():
            for path in sorted(root.iterdir()):
                match = PRODUCT.fullmatch(path.name)
                if not match or path.is_symlink() or not path.is_file():
                    continue
                try:
                    datetime.strptime(match[1], '%Y%m%d')
                except ValueError:
                    continue
                if match[1] < cutoff or not KEEP.fullmatch(path.name):
                    path.unlink()
                    removed.append(path.name)
        # Preserve prior deletion evidence if publication invokes cleanup again.
        if record.exists():
            import json
            removed = sorted(set(removed+json.loads(record.read_text()).get('removed', [])))
        result = dict(checked_at=now.isoformat(),directory=str(root),cutoff_night=cutoff,
                      retention_nights=days,removed=removed,current_preserved=True)
        write_json(record,result)
        return result
