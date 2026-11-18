"""Install a validated staging tree while preserving public scientific output paths."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3


def install(source,base,archive):
    source=Path(source).absolute();base=Path(base).absolute();archive=Path(archive).absolute()
    expected=json.loads((source/'server_expected_sha.json').read_text())
    changed=[name for name,digest in expected.items() if not (base/name).is_file()
             or hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest]
    if changed:raise RuntimeError('Live source changed since inventory: '+', '.join(changed))
    if not (archive/'STOPPED').exists():raise RuntimeError('A completed stopped v1 archive is required')
    if not (base/'data/MPC_SUBMISSION_HOLD').exists():raise RuntimeError('MPC hold must remain active')
    saved=archive/'external_state';saved.mkdir(exist_ok=True)
    if not (saved/'script').exists():shutil.copytree(base/'script',saved/'script')
    for p in [base/'data/MPC_SUBMISSION_HOLD',base/'data/heliolincrr/trksub_history.jsonl',
              *list((base/'data/heliolincrr').glob('*state*.json'))]:
        if p.is_file():shutil.copy2(p,saved/p.name)
    source_files=[]
    for package in ('survey','known_asteroid','heliolincrr','targeted_observing'):
        for p in (source/package).iterdir():
            if p.is_file() and p.suffix in ('.py','.sh','.json','.xsd','.md'):
                target=base/package/p.name;target.parent.mkdir(exist_ok=True)
                temporary=target.with_name('.'+target.name+'.v2-install')
                shutil.copy2(p,temporary);temporary.replace(target)
                source_files.append(target)
    target=base/'targeted_observing'
    for name in ('data','cache','unknown_pipeline'):
        shutil.copytree(source/'targeted_observing'/name,target/name,dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    p=target/'products/annual/ephemeris_index.json';p.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source/'targeted_observing/products/annual/ephemeris_index.json',p)
    live=base/'v2';live.mkdir(exist_ok=True)
    config=live/'config.json'
    if not config.exists():shutil.copy2(target/'server_defaults.json',config)
    history=source/'v2/runtime/history.sqlite';destination=live/'runtime/history.sqlite'
    if history.exists():
        destination.parent.mkdir(parents=True,exist_ok=True)
        if destination.exists():raise RuntimeError('Live history exists: merge explicitly instead of overwriting it')
        with sqlite3.connect('file:'+str(history)+'?mode=ro',uri=True) as src:
            with sqlite3.connect(destination) as dst:src.backup(dst)
    for action,module in [('run_daily','server_service daily'),('submit_followup','server_followup')]:
        # The module string is a fixed installer constant, never user-controlled shell text.
        text='#!/usr/bin/env bash\nset -euo pipefail\ncd '+str(base)+'\nexport PYTHONPATH='+str(base)+'${PYTHONPATH:+:${PYTHONPATH}}\nexec /home/smtpipeline/Softwares/miniconda3/bin/python -u -m targeted_observing.'+module+' "$@"\n'
        p=live/(action+'.sh');p.write_text(text);p.chmod(0o755)
    manifest=dict(version='SHARP-v2-20261118',v1_archive=str(archive),public_paths_preserved=True,
        code_sha256={str(p.relative_to(base)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
        mpc_submission='held',automatic_schedule='not installed yet; publish and validate first')
    (live/'deployment.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(installed=True,files=len(source_files),base=str(base),archive=str(archive))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True);p.add_argument('--base',required=True);p.add_argument('--archive',required=True)
    a=p.parse_args();install(a.source,a.base,a.archive)
