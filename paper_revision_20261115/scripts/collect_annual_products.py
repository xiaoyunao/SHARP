#!/usr/bin/env python3
"""Collect existing science columns and stage/review provenance into a paper-only directory."""
import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from datetime import date,timedelta,datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.table import Table

ROOT=Path('/pipeline/xiaoyunao')
KNOWN={'name':'asteroid_name','number':'asteroid_number','source_file':'source_file','epoch':'epoch_mjd',
 'ra':'pred_ra_deg','dec':'pred_dec_deg','mag':'pred_v_mag','RA_Win':'obs_ra_win_deg','DEC_Win':'obs_dec_win_deg',
 'X_Win':'x_win_px','Y_Win':'y_win_px','objID':'obj_id','Mag_Aper4':'mag_aper4','MagErr_Aper4':'magerr_aper4',
 'Flux_Aper4':'flux_aper4','FluxErr_Aper4':'fluxerr_aper4','Flag':'flag','delta':'pred_geocentric_distance_au'}


def record(p,kind,night):
 st=p.stat();return dict(night=night,kind=kind,path=str(p),bytes=st.st_size,mtime_utc=datetime.fromtimestamp(st.st_mtime,timezone.utc).isoformat(),sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def read(p,columns=None):
 t=Table.read(p)
 if columns:t=t[[c for c in columns if c in t.colnames]]
 d=t.to_pandas()
 for c in d.select_dtypes('object'):
  d[c]=d[c].map(lambda x:x.decode() if isinstance(x,bytes) else x)
 return d


def night(n):
 out={};records=[];status=dict(night=n);base=Path('/processed1')/n/'L4'
 for kind,suffix in [('known_matched','matched_asteroids'),('known_mask15','matched_asteroids_mask15')]:
  p=base/f'{n}_{suffix}.fits';status[kind+'_exists']=p.exists()
  if n>'20260715' and p.exists():
   d=read(p,list(KNOWN)).rename(columns=KNOWN);d['night']=n
   d['asteroid_number']=pd.to_numeric(d.asteroid_number,errors='coerce').astype('Int64').astype('string')
   d['obj_id']=pd.to_numeric(d.obj_id).astype('Int64').astype(str)
   out[kind]=d;records.append(record(p,kind,n))
 for kind,columns in [('orbit_links',['linkage_id','fit_ok','is_good']),
                       ('orbit_obs_residuals',['linkage_id','obs_key','used'])]:
  p=ROOT/'data/heliolincrr'/n/'rr_links/orbit_confirm'/f'{kind}.fits';status[kind+'_exists']=p.exists()
  if p.exists():
   try:
    d=read(p,columns);d['night']=n;out[kind]=d;status[kind+'_rows']=len(d);status[kind+'_readable']=True;records.append(record(p,kind,n))
   except Exception as exc:status[kind+'_error']=str(exc)
 p=base/f'{n}_unknown_links.fits';status['automatic_exists']=p.exists()
 if p.exists():
  d=read(p,['trk_sub','linkage_id','n_obs','fit_ok','is_good','image_names','objids','mjds','ras_deg','decs_deg','rms_arcsec','max_arcsec','lin_rms_arcsec','lin_speed_arcsec_per_day'])
  d['night']=n;out['unknown_catalog']=d;status['automatic_rows']=len(d);records.append(record(p,'unknown_catalog',n))
 reviews=ROOT/'heliolincrr/review_packages'/n
 for kind,suffix in [('review','unknown_review'),('submit','submit')]:
  p=reviews/(f'{n}_{suffix}.csv' if kind=='review' else f'{n}_submit.csv')
  status[kind+'_exists']=p.exists()
  if p.exists():
   try:d=pd.read_csv(p,dtype=str);d['night']=n;out[kind]=d;status[kind+'_rows']=len(d);records.append(record(p,kind,n))
   except Exception as exc:status[kind+'_error']=str(exc)
 # Status snapshots can establish that a saved all-zero review is completed.
 for name in [f'{n}_known_asteroid_status.json',f'{n}_unknown_links_submit_stats.json']:
  p=base/name
  if p.exists():
   try:status[name]=json.loads(p.read_text())
   except Exception:pass
 status['unknown_reply_exists']=(base/f'{n}_unknown_mpc_reply.txt').exists()
 summary=ROOT/'data/heliolincrr'/n/'analysis'/f'{n}_single_night_summary.json'
 if summary.exists():
  try:status['summary']=json.loads(summary.read_text())
  except Exception:pass
 return out,status,records


def main():
 a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);a.add_argument('--workers',type=int,default=4);x=a.parse_args();x.output.mkdir(parents=True,exist_ok=True)
 days=[(date(2025,11,15)+timedelta(days=i)).strftime('%Y%m%d') for i in range(366)]
 chunks={};states=[];provenance=[]
 with ProcessPoolExecutor(x.workers) as pool:
  for d,state,rec in pool.map(night,days):
   for k,frame in d.items():chunks.setdefault(k,[]).append(frame)
   states.append(state);provenance.extend(rec)
 for k,parts in chunks.items():
  pd.concat(parts,ignore_index=True).to_csv(x.output/f'{k}_incremental.csv.gz',index=False)
  print(k,sum(map(len,parts)),flush=True)
 (x.output/'annual_archive_states.json').write_text(json.dumps(states,indent=2,default=str))
 pd.DataFrame(provenance).to_csv(x.output/'annual_product_provenance.csv.gz',index=False)
 print('DONE',len(states),len(provenance),flush=True)

if __name__=='__main__':main()
