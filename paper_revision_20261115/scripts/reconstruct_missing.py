#!/usr/bin/env python3
"""Server-only offline recovery into the paper directory; frozen methods are unmodified."""
import argparse,json,subprocess,os,hashlib
from pathlib import Path
import pandas as pd

BASE=Path('/pipeline/xiaoyunao/paper_revision_20261115')
PY='/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python'
KNOWN_PY='/home/smtpipeline/Softwares/miniconda3/bin/python'

def run(args):
 print('RUN',' '.join(map(str,args)),flush=True)
 subprocess.run(list(map(str,args)),check=True,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})

def main():
 a=argparse.ArgumentParser();a.add_argument('part',choices=['links','known']);x=a.parse_args()
 methods=BASE/'frozen_methods'
 ref=Path('/pipeline/xiaoyunao/v2/runtime/reference_updates/previous/astorb.dat')
 assert hashlib.sha256(ref.read_bytes()).hexdigest()=='ea22c7062e7349d9add34fb64835916d66705ad60a717cc7c86a0c47f5aa3040', 'Frozen ASTORB must be restored before offline execution'
 if x.part=='known':
  for f in ['OBJ_MP_1526_0047_cat.fits.gz','OBJ_MP_1787_0048_cat.fits.gz']:
   run([KNOWN_PY,methods/'match_single_night.py','20260605','--root','/processed1','--file',f,'--outdir',BASE/'offline_known','--parts-dir',BASE/'offline_known/parts','--astorb','/pipeline/xiaoyunao/v2/runtime/reference_updates/previous/astorb.dat','--mask-sep-arcsec','1.5','--njobs','2','--exptime-key','EXPOSURE'])
  return
 night='20260414';root=BASE/'offline_links';inputroot=BASE/'offline_inputs';catdir=inputroot/night/'L2';catdir.mkdir(parents=True,exist_ok=True)
 paths=json.loads((BASE/'reconstruction_catalogs.json').read_text())
 for name in paths:
  p=catdir/name
  if not p.exists():p.symlink_to(Path('/processed1')/night/'L2'/name)
 track=root/night/'tracklets_linreproj';track.mkdir(parents=True,exist_ok=True);allfile=track/f'tracklets_{night}_ALL.fits';rr=root/night/'rr_links'
 run([PY,methods/'mask_gaia.py',night,'--root',inputroot,'--out-root',root,'--gaia-dir','/pipeline/ref/healpix','--gaia-cone-deg','5.0','--match-arcsec','1.5','--pm-zp','2016.0','--hdu-header','1','--mag-psf-max','21.0','--require-flag0','--nproc','4'])
 run([PY,methods/'make_tracklet_linreproj.py',night,'--root',root,'--input-subdir','mask_gaia','--outdir',track,'--mag-max','0','--vmin','3','--vmax','63','--dmag-max','1','--r-static','2','--min-repeat','2','--erode-pix','30','--fallback-margin','0.1','--hdu','1','--nproc','4','--edge-pix','500','--skip-common-area','--max-tracklets-per-group','100000'])
 run([PY,methods/'merge_tracklets_night.py',night,'--root',root,'--tracklet-dir',track,'--out',allfile,'--overwrite','--write-index','--allow-empty'])
 run([PY,methods/'run_linear_links_from_tracklets.py','--infile',allfile,'--outdir',rr,'--speed-thresh-arcsec-per-hour','5','--direction-thresh-deg','10','--require-shared-endpoint'])
 run([PY,methods/'orbit_confirm_links.py','--rr-dir',rr,'--tracklets',allfile,'--outdir',rr/'orbit_confirm','--cores','4','--log-every','500'])

if __name__=='__main__':main()
