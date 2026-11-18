"""Only single-night known extraction and unknown search, sharing one input manifest."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

from astropy.table import Table

from .science_manifest import build as build_manifest, load as load_manifest
from .server_common import load_config, runtime, lock, write_json, utc_now, TZ, atomic_write
from .server_history import ingest


def run(command, config, cwd=None, acceptable=()):
    env=os.environ.copy();env['PYTHONPATH']=str(config['base'])+os.pathsep+env.get('PYTHONPATH','')
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    print('[RUN]', ' '.join(map(str,command)),flush=True)
    result=subprocess.run(list(map(str,command)),cwd=cwd,env=env)
    if result.returncode and result.returncode not in acceptable:
        raise subprocess.CalledProcessError(result.returncode,command)
    return result.returncode


def known_file(task):
    row,night,config,parts=task;root=Path(config['base'])/'known_asteroid'
    run([config['python'],root/'match_single_night.py',night,'--root',config['processed_root'],
         '--outdir',parts.parent,'--parts-dir',parts,'--file',row['file'],
         '--file-regex',r'.*\.(fits|fits\.gz)$','--astorb',root/'astorb.dat',
         '--sep-arcsec',1,'--mask-sep-arcsec',1.5,'--njobs',1,
         '--log',parts/(row['exposure_uid']+'.log')],config)
    return row['file']


def snapshot_previous(night,config,l4,work):
    destination=runtime(config)/'result_revisions'/night/datetime.now(TZ).strftime('%Y%m%dT%H%M%S%f')
    destination.mkdir(parents=True)
    for path in l4.iterdir():
        if path.is_file() and path.name.startswith(night+'_'):
            shutil.copy2(path,destination/path.name)
    # Move only derived work and old package artifacts, never raw/L1/L2.
    for folder in ('mask_gaia','tracklets_linreproj','rr_links','analysis'):
        p=work/folder
        if p.exists():shutil.move(str(p),destination/folder)
    for name in ('known_asteroid_parts',):
        if (l4/name).exists():shutil.move(str(l4/name),destination/name)
    for prefix,label in ((Path(config['base'])/'heliolincrr/review_packages','review_package'),
                         (Path(config['base'])/'heliolincrr/plots','plots'),
                         (Path(config['base'])/'known_asteroid/plots','known_plots')):
        if (prefix/night).exists():shutil.move(str(prefix/night),destination/label)
    return destination


def process(night,config):
    base=Path(config['base']);l4=Path(config['processed_root'])/night/'L4'
    work=base/'data/heliolincrr'/night
    state=runtime(config)/'processing'/f'{night}.json'
    l4.mkdir(parents=True,exist_ok=True);work.mkdir(parents=True,exist_ok=True)
    with lock(runtime(config)/'locks'/f'processing_{night}.lock'):
        manifest_path=runtime(config)/'manifests'/f'{night}_science.json'
        manifest=build_manifest(night,config['processed_root'],manifest_path)
        atomic_write(l4/f'{night}_file_manifest.txt',manifest_path.with_suffix('.txt').read_bytes())
        write_json(l4/f'{night}_science_manifest.json',manifest)
        if not manifest['accepted']:
            result=dict(night=night,status='no_valid_science_inputs',known_complete=False,unknown_complete=False,
                        rejected=manifest['rejected'],updated_at=utc_now(),previous_results_kept_for_audit=True)
            old=l4/f'{night}_known_asteroid_status.json'
            if old.exists():
                saved=runtime(config)/'result_revisions'/night/('invalid_input_status_'+datetime.now(TZ).strftime('%Y%m%dT%H%M%S%f')+'.json')
                saved.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(old,saved)
            write_json(old,result)
            write_json(l4/f'{night}_v2_pipeline_status.json',result)
            write_json(state,result)
            # Do not overwrite a previous valid scientific result with a false
            # empty-complete state when every input is stale or malformed.
            return result
        previous=json.loads(state.read_text()) if state.exists() else {}
        code_files=[Path(__file__),Path(__file__).with_name('science_manifest.py'),
                    Path(__file__).parent/'unknown_pipeline/exposure_timing.py',
                    Path(__file__).with_name('catalog_names.py')]
        code_files+=sorted((base/'known_asteroid').glob('*.py'))+sorted((base/'heliolincrr').glob('*.py'))
        code_hash=hashlib.sha256(b''.join(p.read_bytes() for p in code_files)).hexdigest()
        identity=manifest['input_signature']+code_hash
        if previous.get('identity')==identity and previous.get('status')=='complete':return previous
        result=dict(night=night,status='running',identity=identity,input_signature=manifest['input_signature'],
                    inputs=len(manifest['accepted']),known_complete=False,unknown_complete=False,started_at=utc_now())
        write_json(state,result)
        try:
            result['previous_results']=str(snapshot_previous(night,config,l4,work))
            write_json(l4/f'{night}_known_asteroid_status.json',dict(night=night,known_complete=False,
                status='running',pipeline_version=config['version'],input_signature=manifest['input_signature']))
            parts=l4/'known_asteroid_parts';parts.mkdir(exist_ok=True)
            # Old automatic pipelines are retired: successful per-file tasks are
            # the only denominator, including valid files with zero asteroids.
            with ThreadPoolExecutor(max_workers=config['workers']) as pool:
                completed=list(pool.map(known_file,[(r,night,config,parts) for r in manifest['accepted']]))
            load_manifest(manifest_path)  # Fail if input changed during extraction.
            run([config['python'],base/'known_asteroid/merge_night_parts.py',night,
                 '--parts-dir',parts,'--outdir',l4],config)
            known_status_path=l4/f'{night}_known_asteroid_status.json'
            known_status=json.loads(known_status_path.read_text())
            known_status.update(pipeline_version=config['version'],input_signature=manifest['input_signature'],
                                input_files=completed,known_complete=True,submission_status='held_contacts_pending')
            write_json(known_status_path,known_status)
            result['known_complete']=True;write_json(state,result)
            ingest(config,[night])
            tonight=datetime.now(TZ).strftime('%Y%m%d')
            if night>=config['history_start']:
                event=runtime(config)/'events'/tonight/('history_'+night+'_'+manifest['input_signature'][:16]+'.json')
                write_json(event,dict(reason='new_actual_known_detections',night=night,
                    request_sha256='history:'+manifest['input_signature'],updated_at=utc_now()))
            # Generate the existing known PSV without invoking a submission API.
            matched=l4/f'{night}_matched_asteroids.fits'
            if matched.exists() and len(Table.read(matched)):
                run([config['python'],base/'known_asteroid/export_ades.py','--fits',matched,
                     '--out',l4/f'{night}_matched_asteroids_ades.psv','--obs-code','327'],config)
            if config.get('known_auto_submit'):
                if not config.get('submission_contacts_confirmed'):
                    raise RuntimeError('Known automatic submission requires confirmed contacts')
                if not (base/'data/MPC_SUBMISSION_HOLD').exists():
                    run(['bash',base/'known_asteroid/slurm_merge_submit.sh',night,'true'],config)
            unknown=base/'heliolincrr';py=config['unknown_python'];workers=config['workers']
            # All arguments and paths are explicit; there is no 15-night stage,
            # old followup state, unknown submit export or watcher invocation.
            commands=[
                [py,unknown/'mask_gaia.py',night,'--root',config['processed_root'],
                 '--out-root',base/'data/heliolincrr','--science-manifest',manifest_path,
                 '--gaia-dir','/pipeline/ref/healpix','--gaia-cone-deg','5','--match-arcsec','1.5',
                 '--mag-psf-max','21','--require-flag0','--max-residual-rows','3000','--nproc',workers],
                [py,unknown/'make_tracklet_linreproj.py',night,'--root',base/'data/heliolincrr',
                 '--input-subdir','mask_gaia','--outdir',work/'tracklets_linreproj',
                 '--science-manifest',manifest_path,'--mag-max','0','--vmin','3','--vmax','750',
                 '--dmag-max','1','--r-static','2','--min-repeat','2','--erode-pix','30',
                 '--fallback-margin','.1','--nproc',workers,'--edge-pix','500',
                 '--skip-common-area','--max-tracklets-per-group','100000'],
                [py,unknown/'merge_tracklets_night.py',night,'--root',base/'data/heliolincrr',
                 '--tracklet-dir',work/'tracklets_linreproj','--out',work/'tracklets_linreproj'/f'tracklets_{night}_ALL.fits',
                 '--overwrite','--write-index','--allow-empty'],
                [py,unknown/'run_linear_links_from_tracklets.py','--infile',work/'tracklets_linreproj'/f'tracklets_{night}_ALL.fits',
                 '--outdir',work/'rr_links','--speed-thresh-arcsec-per-hour','5','--direction-thresh-deg','10','--require-shared-endpoint'],
                [py,unknown/'orbit_confirm_links.py','--rr-dir',work/'rr_links','--tracklets',work/'tracklets_linreproj'/f'tracklets_{night}_ALL.fits',
                 '--cores',workers,'--log-every','500'],
                [py,unknown/'summarize_single_night.py',night,'--processed-root',config['processed_root'],
                 '--root-out',base/'data/heliolincrr','--rr-subdir','rr_links',
                 '--summary-json',work/'analysis'/f'{night}_single_night_summary.json',
                 '--unknown-json',l4/f'{night}_unknown_links.json','--unknown-fits',l4/f'{night}_unknown_links.fits','--require-mask15'],
            ]
            for command in commands:
                load_manifest(manifest_path)
                run(command,config,cwd=unknown)
            with lock(runtime(config)/'locks/trksub.lock'):
                run([py,unknown/'assign_unknown_trksub.py',night,'--mode','single-night',
                     '--catalog',l4/f'{night}_unknown_links.json','--fits',l4/f'{night}_unknown_links.fits',
                     '--summary-json',work/'analysis'/f'{night}_single_night_summary.json',
                     '--history',base/'data/heliolincrr/trksub_history.jsonl'],config,cwd=unknown)
            run([py,unknown/'unknown_morphology.py',night,'--processed-root',config['processed_root'],
                 '--workers',min(workers,4)],config,cwd=unknown)
            limited=run([py,unknown/'check_unknown_review_count.py',night,'--processed-root',config['processed_root'],
                         '--maximum',config['max_unknown_review']],config,cwd=unknown,acceptable=(20,))==20
            if not limited:
                run([py,unknown/'plot_unknown_links.py',night,'--processed-root',config['processed_root'],
                     '--root-out',base/'data/heliolincrr','--plot-root',unknown/'plots'],config,cwd=unknown)
                run([py,unknown/'package_unknown_review.py',night,'--processed-root',config['processed_root'],
                     '--root-out',base/'data/heliolincrr','--plot-root',unknown/'plots',
                     '--out-root',unknown/'review_packages','--make-tar'],config,cwd=unknown)
            load_manifest(manifest_path)
            if matched.exists() and len(Table.read(matched)):
                with lock(runtime(config)/'locks/known_history.lock'):
                    run([config['python'],base/'known_asteroid/update_all_matched_history.py',
                         '--processed-root',config['processed_root'],'--night',night,
                         '--history-dir',base/'known_asteroid/runtime/history'],config)
                run([config['python'],base/'known_asteroid/plot_known_asteroids.py',night,
                     '--processed-root',config['processed_root'],'--plot-root',base/'known_asteroid/plots',
                     '--science-manifest',manifest_path,
                     '--all-matched-path',base/'known_asteroid/runtime/history/all_matched_asteroids.fits'],config)
            trk_log=work/'tracklets_linreproj/make_tracklet.log'
            gate_lines=[line for line in trk_log.read_text().splitlines() if 'exceeds max_tracklets_per_group=' in line] if trk_log.exists() else []
            result.update(status='complete',unknown_complete=True,review_package_complete=not limited,
                          dense_groups_excluded=gate_lines,
                          review_limited=limited,unknown_submission='website_only',completed_at=utc_now())
            write_json(l4/f'{night}_v2_pipeline_status.json',result);write_json(state,result)
            return result
        except Exception as exc:
            result.update(status='failed',error=str(exc),updated_at=utc_now())
            write_json(state,result)
            raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--night',required=True);p.add_argument('--config')
    a=p.parse_args();print(json.dumps(process(a.night,load_config(a.config)),indent=2))
