#!/usr/bin/env python3
"""Run explicit images in an ephemeral SSH directory; copy and verify, then remove it."""
import argparse
import hashlib
import json
import shlex
import subprocess
from pathlib import Path

HOST='smtpipeline@www.xinglong-naoc.cn'
PYTHON='/home/smtpipeline/Softwares/miniconda3/bin/python'
HERE=Path(__file__).resolve().parent


def ssh(command,**kwargs):
    return subprocess.run(['ssh','-p','20093','-o','BatchMode=yes',HOST,command],check=True,**kwargs)


def run(images,destination,final_through=None,band="r",allow_header_mismatch=False):
    destination=Path(destination).resolve();destination.mkdir(parents=True,exist_ok=True)
    temp=ssh('mktemp -d /tmp/smt_focus_XXXXXXXX',capture_output=True,text=True).stdout.strip()
    if not temp.startswith('/tmp/smt_focus_') or '/' in temp[len('/tmp/'):]:
        raise RuntimeError('Unexpected temporary directory: '+repr(temp))
    record=dict(remote_directory=temp,images=images,band=band,allow_header_mismatch=allow_header_mismatch,final_through=final_through,cleaned=False)
    (destination/'last_remote_run.json').write_text(json.dumps(record,indent=2))
    try:
        subprocess.run(['scp','-P','20093',str(HERE/'measure_focus.py'),str(HERE/'compare_focus.py'),HOST+':'+temp+'/'],check=True)
        ssh('mkdir '+shlex.quote(temp+'/results'))
        previous=sorted(destination.glob('*_sources.csv'))+sorted(destination.glob('*_summary.json'))
        if previous:
            subprocess.run(['scp','-P','20093',*[str(p) for p in previous],HOST+':'+temp+'/results/'],check=True)
        # Do not let bytecode, plotting caches, or temporary files escape this directory.
        command=' '.join(['env', 'MPLCONFIGDIR='+shlex.quote(temp+'/mpl'), 'TMPDIR='+shlex.quote(temp),
                          'OPENBLAS_NUM_THREADS=1', shlex.quote(PYTHON),'-B',shlex.quote(temp+'/measure_focus.py'),
                          *[shlex.quote(p) for p in images],'--out',shlex.quote(temp+'/results'),'--band',shlex.quote(band),*(['--allow-header-mismatch'] if allow_header_mismatch else [])])
        if images: ssh(command)
        ssh(' '.join(['env','MPLCONFIGDIR='+shlex.quote(temp+'/mpl'),'TMPDIR='+shlex.quote(temp),'OPENBLAS_NUM_THREADS=1',shlex.quote(PYTHON),'-B',shlex.quote(temp+'/compare_focus.py'),shlex.quote(temp+'/results'),*(['--final-through',str(int(final_through))] if final_through is not None else [])]))
        # Capture a manifest before copying; verify every returned artifact locally.
        manifest_code="from pathlib import Path; import hashlib,json; p=Path("+repr(temp+'/results')+"); print(json.dumps({x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in p.iterdir() if x.is_file()}))"
        hashes=json.loads(ssh(shlex.quote(PYTHON)+' -B -c '+shlex.quote(manifest_code),capture_output=True,text=True).stdout)
        subprocess.run(['scp','-P','20093','-r',HOST+':'+temp+'/results/.',str(destination)+'/'],check=True)
        for name,digest in hashes.items():
            if hashlib.sha256((destination/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Transfer verification failed: '+name)
        record['verified_sha256']=hashes
    finally:
        # Delete only the mktemp directory owned by this invocation, never science data.
        code='import shutil; from pathlib import Path; p=Path('+repr(temp)+'); shutil.rmtree(p); assert not p.exists()'
        try:
            ssh(shlex.quote(PYTHON)+' -B -c '+shlex.quote(code));record['cleaned']=True
        finally:(destination/'last_remote_run.json').write_text(json.dumps(record,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('images',nargs='*');p.add_argument('--out',required=True);p.add_argument('--final-through',type=int)
    p.add_argument('--band',choices=['g','r','i','z','unconfirmed'],default='r');p.add_argument('--allow-header-mismatch',action='store_true')
    a=p.parse_args();run(a.images,a.out,a.final_through,a.band,a.allow_header_mismatch)
