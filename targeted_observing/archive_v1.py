"""Explicit, reversible archive of the live v1 code trees and cron configuration."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from datetime import datetime


def archive(base):
    base = Path(base).absolute()
    destination = base/'v1'/datetime.now().strftime('%Y%m%d_%H%M%S')
    destination.mkdir(parents=True, exist_ok=False)
    old = subprocess.check_output(['crontab', '-l'], text=True)
    (destination/'crontab.before').write_text(old)
    disabled = []
    for line in old.splitlines():
        if not line.lstrip().startswith('#') and any(str(prefix/entry) in line
                for prefix in (base,base.resolve()) for entry in (
                'heliolincrr/run_daily_pipeline.sh','heliolincrr/run_review_submit_backlog_watch.sh',
                'survey/run_daily.sh','known_asteroid/run_daily.sh')):
            line = '# V1_DISABLED ' + line
        disabled.append(line)
    stopped = '\n'.join(disabled)+'\n'
    (destination/'crontab.stopped').write_text(stopped)
    subprocess.run(['crontab', '-'], input=stopped, text=True, check=True)
    # Uncompressed tar preserves all original code, dependencies and in-tree
    # runtime/results without rewriting the active scientific data directories.
    target = destination/'production_trees.tar'
    subprocess.run(['tar', '-cpf', str(target)+'.partial', 'survey', 'known_asteroid',
                    'heliolincrr'], cwd=base, check=True)
    Path(str(target)+'.partial').rename(target)
    count = 0
    import tarfile
    with tarfile.open(target) as stream:
        for _ in stream:
            count += 1
    digest = hashlib.file_digest(target.open('rb'), 'sha256').hexdigest()
    (destination/'SHA256SUMS').write_text(digest+'  production_trees.tar\n')
    metadata = dict(archive=str(destination), entries=count, bytes=target.stat().st_size,
                    sha256=digest, original_base=str(base), v1_automatic_entries_disabled=True,
                    shared_external_data=[str(base/'data'), '/processed1', '/raw1'],
                    status='V1_ARCHIVED_STOPPED',
                    note='External data are shared, not duplicated. Original cron is preserved; unrelated cron remains active.')
    (destination/'manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')
    (destination/'STOPPED').write_text('Archived v1. Do not schedule these scripts alongside v2.\n')
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', required=True)
    archive(parser.parse_args().base)
