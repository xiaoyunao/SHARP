"""Enable only the validated v2 noon/tick entries, retaining unrelated cron jobs."""
import argparse
from datetime import datetime
import json
from pathlib import Path
import shlex
import subprocess
from zoneinfo import ZoneInfo


def activate(base):
    base=Path(base).absolute();v2=base/'v2';now=datetime.now(ZoneInfo('Asia/Shanghai'))
    publication=json.loads((v2/'runtime/published'/(now.strftime('%Y%m%d')+'.json')).read_text())
    if not publication['verification']['passed']:raise RuntimeError('Validate a current-night plan before enabling cron')
    if not (base/'data/MPC_SUBMISSION_HOLD').exists():raise RuntimeError('MPC hold is required')
    old=subprocess.check_output(['crontab','-l'],text=True)
    backup=v2/'crontab.before_activation'
    if not backup.exists():backup.write_text(old)
    lines=[];inside=False
    for line in old.splitlines():
        if line in ('# BEGIN SHARP_V2','# BEGIN SHARP_REFERENCE_UPDATE'):inside=True;continue
        if line in ('# END SHARP_V2','# END SHARP_REFERENCE_UPDATE'):inside=False;continue
        if not inside:lines.append(line)
    logs=v2/'runtime/logs';logs.mkdir(parents=True,exist_ok=True)
    def command(action,log):
        args=['/home/smtpipeline/Softwares/miniconda3/bin/python','-u','-m','targeted_observing.server_service',action]
        return 'cd '+shlex.quote(str(base))+' && '+shlex.join(args)+' >> '+shlex.quote(str(logs/log))+' 2>&1'
    lines+=['','# BEGIN SHARP_V2','0 12 * * * '+command('daily','cron_noon.log'),
            '* * * * * '+command('tick','cron_tick.log'),
            '@reboot sleep 60 && '+command('tick','reboot_recovery.log'),'# END SHARP_V2']
    reference_log=v2/'runtime/reference_updates/cron.log';reference_log.parent.mkdir(parents=True,exist_ok=True)
    ref_args=['/home/smtpipeline/Softwares/miniconda3/bin/python','-u','-m',
              'targeted_observing.reference_update','--enqueue']
    lines+=['','# BEGIN SHARP_REFERENCE_UPDATE','20 10 * * * cd '+shlex.quote(str(base))+' && '+
            shlex.join(ref_args)+' >> '+shlex.quote(str(reference_log))+' 2>&1','# END SHARP_REFERENCE_UPDATE']
    value='\n'.join(lines)+'\n';subprocess.run(['crontab','-'],input=value,text=True,check=True)
    (v2/'crontab.active').write_text(value)
    path=v2/'deployment.json';metadata=json.loads(path.read_text())
    metadata.update(automatic_schedule='12:00 Asia/Shanghai scheduler; per-minute L2/event/recovery tick; 10:20 monthly reference due/retry check; reboot recovery',
                    activated_at=now.isoformat())
    path.write_text(json.dumps(metadata,indent=2)+'\n')
    print(value)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--base',required=True)
    activate(p.parse_args().base)
