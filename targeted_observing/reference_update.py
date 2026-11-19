"""Validate and atomically refresh Lowell ASTORB and IERS once per ISO week."""
import argparse
from datetime import datetime, timedelta
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.request

from .server_common import TZ, load_config, runtime, write_json, lock, atomic_write, utc_now

ASTORB_URL = 'https://ftp.lowell.edu/pub/elgb/astorb.dat.gz'
IERS_URLS = ('https://datacenter.iers.org/data/9/finals2000A.all',
             'https://maia.usno.navy.mil/ser7/finals2000A.all')


def sha(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024), b''):value.update(chunk)
    return value.hexdigest()


def week_key(when=None):
    when=when or datetime.now(TZ)
    year,week,_=when.astimezone(TZ).isocalendar()
    return f'{year}-W{week:02d}'


def due(config, week=None):
    week=week or week_key()
    path=runtime(config)/'reference_updates/status.json'
    previous=json.loads(path.read_text()) if path.exists() else {}
    completed=previous.get('completed_week')
    if not completed and previous.get('completed_at'):
        # Migrate a successful monthly receipt without downloading again.
        completed=week_key(datetime.fromisoformat(previous['completed_at']))
    return completed != week


def download(urls, path):
    errors=[]
    for url in urls:
        try:
            request=urllib.request.Request(url,headers={'User-Agent':'SHARP-reference-maintenance/2'})
            with urllib.request.urlopen(request,timeout=60) as response, path.open('wb') as stream:
                shutil.copyfileobj(response,stream)
                expected=response.headers.get('Content-Length')
                if expected and path.stat().st_size!=int(expected):raise ValueError('Incomplete download')
                modified=response.headers.get('Last-Modified')
            return dict(url=url,last_modified=modified,downloaded_at=utc_now(),sha256=sha(path))
        except Exception as exc:
            path.unlink(missing_ok=True);errors.append(f'{url}: {exc}')
    raise RuntimeError('; '.join(errors))


def validate_astorb(path):
    if path.stat().st_size<100_000_000:raise ValueError('ASTORB unexpectedly small')
    # Validate with the same parser as production, without its optional SQL cache.
    from aleph.Query import Query
    query=Query(service='Lowell',filename=str(path),open_db=False)
    count=len(query.asts)
    if count<1_000_000:raise ValueError('ASTORB object count unexpectedly small')
    import numpy as np
    epochs=np.array([float(t.mjd) for t in query.ref_epochs.values()])
    if not np.isfinite(epochs).all():raise ValueError('Invalid ASTORB epochs')
    return dict(objects=count,reference_epoch_mjd_min=float(epochs.min()),reference_epoch_mjd_max=float(epochs.max()))


def validate_iers(path):
    from astropy.utils import iers
    from astropy.time import Time
    table=iers.IERS_A.open(str(path))
    now=Time(datetime.now(TZ));last=float(table['MJD'][-1].value)
    if float(table['MJD'][0].value)>now.mjd or last<now.mjd+40:
        raise ValueError('IERS must cover today and at least the next forty days')
    age=float(now.mjd-table.meta['predictive_mjd'])
    if age>35:raise ValueError('IERS predictions too old')
    return dict(last_mjd=last,prediction_age_days=age)


def prepare_bundle(folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    receipt={'astorb.dat.gz':download([ASTORB_URL],folder/'astorb.dat.gz'),
             'finals2000A.all':download(IERS_URLS,folder/'finals2000A.all')}
    write_json(folder/'download.json',receipt)
    return receipt


def update(config, bundle=None, force=False):
    root=runtime(config)/'reference_updates';root.mkdir(parents=True,exist_ok=True)
    week=week_key();state=root/'status.json'
    prepared=root/('bundle_'+week)
    if bundle is None and (prepared/'download.json').exists():bundle=prepared
    if not force and not due(config,week):return dict(status='already_updated_this_week',week=week)
    with lock(runtime(config)/'locks/reference_update.lock',blocking=False):
        if not force and not due(config,week):return dict(status='already_updated_this_week',week=week)
        previous=json.loads(state.read_text()) if state.exists() else {}
        result=dict(previous,status='updating',attempt_week=week,started_at=utc_now())
        result.pop('attempt_month',None);result.pop('completed_month',None)
        write_json(state,result)
        try:
            # All downloads and parsing precede activation. A failure leaves the
            # production files and last successful week untouched.
            with tempfile.TemporaryDirectory(prefix='.staging_',dir=root) as folder:
                stage=Path(folder)
                if bundle:
                    for name in ('astorb.dat.gz','finals2000A.all'):
                        shutil.copy2(Path(bundle)/name,stage/name)
                    receipt=json.loads((Path(bundle)/'download.json').read_text())
                    for name in ('astorb.dat.gz','finals2000A.all'):
                        if sha(stage/name)!=receipt[name]['sha256']:raise ValueError('Relay SHA mismatch: '+name)
                else:
                    receipt={'astorb.dat.gz':download([ASTORB_URL],stage/'astorb.dat.gz'),
                             'finals2000A.all':download(IERS_URLS,stage/'finals2000A.all')}
                with gzip.open(stage/'astorb.dat.gz','rb') as source,(stage/'astorb.dat').open('wb') as out:
                    shutil.copyfileobj(source,out)
                validation=dict(astorb=validate_astorb(stage/'astorb.dat'),iers=validate_iers(stage/'finals2000A.all'))
                # Avoid an old, already-running process reading mixed versions.
                with lock(runtime(config)/'locks/references.lock'):
                    active=subprocess.check_output(['squeue','-h','-u','smtpipeline','-o','%j'],text=True).splitlines()
                    if any(name.startswith('sharp_v2_') for name in active):
                        raise RuntimeError('Processing is active; defer reference activation to the next retry')
                    destinations={Path(config['base'])/'known_asteroid/astorb.dat':stage/'astorb.dat',
                                  Path(__file__).parent/'cache/iers_current.all':stage/'finals2000A.all'}
                    backup=root/'previous';backup.mkdir(exist_ok=True)
                    for target,source in destinations.items():
                        if target.exists():shutil.copy2(target,backup/target.name)
                    try:
                        for target,source in destinations.items():
                            # Staging is on the same filesystem as production.
                            source.replace(target)
                        from astropy.utils import iers
                        from astropy.utils.data import import_file_to_cache
                        import_file_to_cache(iers.IERS_A_URL,str(Path(__file__).parent/'cache/iers_current.all'),replace=True)
                        config_path=Path(config['base'])/'v2/config.json'
                        current=json.loads(config_path.read_text());current['iers_file']='iers_current.all'
                        write_json(config_path,current)
                    except Exception:
                        for target in destinations:
                            if (backup/target.name).exists():shutil.copy2(backup/target.name,target)
                            else:target.unlink(missing_ok=True)
                        raise
                    result.pop('error',None);result.pop('local_relay_available',None)
                    result.update(status='complete',completed_week=week,completed_at=utc_now(),
                                  downloads=receipt,validation=validation,
                                  installed={str(p):sha(p) for p in destinations})
                    write_json(root/(week+'_receipt.json'),result);write_json(state,result)
            if bundle==prepared:
                shutil.rmtree(prepared)
            return result
        except Exception as exc:
            result.update(status='retry_required',error=str(exc),last_attempt=utc_now(),local_relay_available=True)
            write_json(state,result)
            raise


def relay(config, remote='smtpipeline@www.xinglong-naoc.cn', port=20093):
    """Download on the local computer and let the server validate/activate."""
    import shlex
    ssh=['ssh','-p',str(port),remote]
    py=config['python'];base=config['base']
    check='from targeted_observing.reference_update import due;from targeted_observing.server_common import load_config;print(int(due(load_config())))'
    output=subprocess.check_output(ssh+['cd '+shlex.quote(base)+' && '+shlex.join([py,'-c',check])],text=True)
    if output.strip()=='0':return dict(status='server_already_updated')
    with tempfile.TemporaryDirectory(prefix='sharp_reference_relay_') as folder:
        folder=Path(folder)
        receipt={'astorb.dat.gz':download([ASTORB_URL],folder/'astorb.dat.gz'),
                 'finals2000A.all':download(IERS_URLS,folder/'finals2000A.all')}
        write_json(folder/'download.json',receipt)
        create='import tempfile;print(tempfile.mkdtemp(prefix=".relay_",dir='+repr(str(Path(base)/'v2/runtime/reference_updates'))+'))'
        remote_dir=subprocess.check_output(ssh+[shlex.join([py,'-c',create])],text=True).strip()
        try:
            for file in folder.iterdir():
                subprocess.run(['scp','-P',str(port),str(file),remote+':'+remote_dir+'/'],check=True)
            command=[py,'-u','-m','targeted_observing.reference_update','--bundle',remote_dir]
            subprocess.run(ssh+['cd '+shlex.quote(base)+' && '+shlex.join(command)],check=True)
        finally:
            clean='from pathlib import Path;import shutil;p=Path('+repr(remote_dir)+');assert p.name.startswith(".relay_");shutil.rmtree(p)'
            subprocess.run(ssh+[shlex.join([py,'-c',clean])],check=True)
    return dict(status='relayed_and_installed')


def enqueue(config):
    import shlex
    if not due(config):return dict(status='already_updated_this_week')
    jobs=subprocess.check_output(['squeue','-h','-u','smtpipeline','-o','%A|%j'],text=True).splitlines()
    dependencies=[]
    for row in jobs:
        job,name=row.split('|',1)
        if name=='sharp_ref_update':return dict(status='already_queued',job_id=job)
        if name.startswith('sharp_v2_'):dependencies.append(job)
    root=runtime(config)/'reference_updates';root.mkdir(parents=True,exist_ok=True)
    args=['sbatch','--parsable','--job-name','sharp_ref_update','--cpus-per-task','1','--mem','8G',
          '--time','02:00:00','--chdir',config['base'],'--output',str(root/'update_%j.log')]
    if dependencies:args+=['--dependency','afterany:'+':'.join(dependencies)]
    args+=['--wrap',shlex.join([config['python'],'-u','-m','targeted_observing.reference_update'])]
    job=subprocess.check_output(args,text=True).strip().split(';')[0]
    result=dict(status='queued',job_id=job,after_processing_jobs=dependencies,queued_at=utc_now())
    write_json(root/'queue.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config');p.add_argument('--bundle',type=Path)
    p.add_argument('--force',action='store_true');p.add_argument('--local-relay',action='store_true')
    p.add_argument('--enqueue',action='store_true')
    p.add_argument('--prepare-bundle',type=Path)
    a=p.parse_args();c=load_config(a.config)
    print(json.dumps(prepare_bundle(a.prepare_bundle) if a.prepare_bundle else enqueue(c) if a.enqueue else relay(c) if a.local_relay else update(c,a.bundle,a.force),indent=2))
