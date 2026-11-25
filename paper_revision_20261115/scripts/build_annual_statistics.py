#!/usr/bin/env python3
"""Annual science summaries on the single selected exposure set; existing products only."""
import json,re
from pathlib import Path
import astropy.units as u
import healpy as hp
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord,BarycentricMeanEcliptic
from astropy.time import Time
from exposures import ROOT,INPUT,TABLES,DATA
from known_performance import (describe,with_residuals,exposure_pairs,geometry,quality_tables,histogram)
from replay_frozen_recovery import object_key,normalize_obj_id,wilson,binned_rows,MEASURED_MAG_BINS,RATE_BINS

OLD=ROOT.parent/'paper_revision_20261015'


def attach(d,led):
 d=d.copy();d['night']=d.night.astype(str);d['obj_id']=d.obj_id.map(normalize_obj_id)
 d['file_name']=d.source_file.astype(str).str.replace(r'_cat\.fits(?:\.gz)?$','.fits',regex=True)
 cols=['night','file_name','noon_night','quality','quality_exclusion','exposure_key','mid_utc','raw_start_utc','l2_available','start_utc']
 d=d.merge(led[cols],on=['night','file_name'],how='left',validate='many_to_one')
 d['object_key']=object_key(d);d['quality']=d.quality.eq(True)
 return with_residuals(d)


def retention(mask,ledger,states):
 d=mask.copy();d['obs_key']=d.source_file+':'+d.obj_id.map(normalize_obj_id)
 mapping=d[['night','obs_key','object_key']].drop_duplicates();n=mapping.groupby(['night','obs_key']).object_key.transform('nunique')
 ambiguous=int(mapping.loc[n.gt(1),['night','obs_key']].drop_duplicates().shape[0]);mapping=mapping.loc[n.eq(1)]
 oldlinks=pd.read_parquet(OLD/'inputs/orbit_links.parquet');oldobs=pd.read_parquet(OLD/'inputs/orbit_obs_residuals.parquet')
 def stages(name,old):
  new=pd.read_parquet(INPUT/f'{name}_incremental.parquet');old.night=old.night.astype(str);new.night=new.night.astype(str)
  parts=[old,new.loc[~new.night.isin(old.night.unique())]]
  recovered=INPUT/f'{name}_reconstructed.parquet'
  if recovered.exists():parts.append(pd.read_parquet(recovered))
  return pd.concat(parts,ignore_index=True)
 links=stages('orbit_links',oldlinks);obs=stages('orbit_obs_residuals',oldobs)
 obs=obs.drop_duplicates(['night','linkage_id','obs_key']);links=links.drop_duplicates(['night','linkage_id'])
 obs=obs.merge(mapping,on=['night','obs_key'],how='left',validate='many_to_one')
 classified=obs.groupby(['night','linkage_id']).agg(n=('obs_key','nunique'),mapped=('object_key','count'),ids=('object_key','nunique')).reset_index()
 pure=classified.loc[classified.n.ge(3)&classified.mapped.eq(classified.n)&classified.ids.eq(1)]
 members=obs.merge(pure[['night','linkage_id']],on=['night','linkage_id'],validate='many_to_one')
 support=d.loc[d.quality,['night','obs_key','object_key','noon_night','exposure_key']].drop_duplicates()
 members=members.merge(support,on=['night','obs_key','object_key'],how='inner',validate='many_to_one')
 supported=members.groupby(['night','linkage_id','noon_night','object_key']).exposure_key.nunique().reset_index(name='n')
 supported=supported.loc[supported.n.ge(3)].merge(links,on=['night','linkage_id'],validate='many_to_one')
 survived=supported.groupby(['noon_night','object_key']).agg(strict_linked=('linkage_id',lambda x:True),fit_ok_survived=('fit_ok','any'),is_good_survived=('is_good','any')).reset_index()
 good=d.loc[d.quality].sort_values(['epoch_mjd','radial_arcsec']).drop_duplicates(['noon_night','object_key','exposure_key'])
 g=good.groupby(['noon_night','object_key'])
 obj=g.agg(detected_1p5arcsec_n=('exposure_key','size'),first_mjd=('epoch_mjd','min'),last_mjd=('epoch_mjd','max'),
  first_ra=('pred_ra_deg','first'),last_ra=('pred_ra_deg','last'),first_dec=('pred_dec_deg','first'),last_dec=('pred_dec_deg','last'),
  median_measured_mag=('mag_aper4','median'),directory_nights=('night',lambda x:'|'.join(sorted(set(x))))).reset_index()
 hours=(obj.last_mjd-obj.first_mjd)*24;dx=((obj.last_ra-obj.first_ra+180)%360-180)*np.cos(np.deg2rad((obj.first_dec+obj.last_dec)/2))
 obj['first_last_rate_arcsec_per_hour']=np.hypot(dx,obj.last_dec-obj.first_dec)*3600/hours.where(hours.gt(0))
 obj['eligible_ge3_and_rate']=obj.detected_1p5arcsec_n.ge(3)&obj.first_last_rate_arcsec_per_hour.between(3,63)
 state={s['night']: bool(s.get('orbit_links_readable')) and bool(s.get('orbit_obs_residuals_readable')) for s in states}
 # The frozen products themselves remain evaluation evidence if the server copy has since gone.
 for n in set(oldlinks.night)&set(oldobs.night):state[n]=True
 if (INPUT/'orbit_links_reconstructed.parquet').exists():
  for n in pd.read_parquet(INPUT/'orbit_links_reconstructed.parquet').night.astype(str).unique():state[n]=True
 obj['evaluated']=obj.directory_nights.map(lambda x:all(state.get(n,False) for n in x.split('|')))
 obj=obj.merge(survived,on=['noon_night','object_key'],how='left',validate='one_to_one')
 for c in ['strict_linked','fit_ok_survived','is_good_survived']:obj[c]=obj[c].eq(True)
 obj.to_parquet(TABLES/'annual_object_nights.parquet',index=False,compression='zstd')
 eligible=obj.loc[obj.eligible_ge3_and_rate];e=eligible.loc[eligible.evaluated].copy();e['night']=e.noon_night
 e.to_parquet(INPUT/'blinded_known_recovery_by_object.parquet',index=False,compression='zstd')
 rows=[]
 for col,edges,label in [('median_measured_mag',MEASURED_MAG_BINS,'measured_aperture_magnitude'),('first_last_rate_arcsec_per_hour',RATE_BINS,'first_last_rate_arcsec_per_hour')]:rows.extend(binned_rows(e,col,edges,label))
 pd.DataFrame(rows).to_csv(DATA/'fig06_retention_frozen_bins.csv',index=False)
 coverage=[]
 for n in sorted(ledger.loc[ledger.quality,'noon_night'].unique()):
  a=eligible.loc[eligible.noon_night.eq(n)];v=a.loc[a.evaluated]
  archive_exists=all(state.get(z,False) for z in ledger.loc[ledger.quality&ledger.noon_night.eq(n),'night'].unique())
  coverage.append(dict(night=n,science_exposures=int((ledger.quality&ledger.noon_night.eq(n)).sum()),eligible=len(a),evaluated=len(v),linked=int(v.strict_linked.sum()),status=('reconstructed' if n=='20260413' else 'archived') if archive_exists and len(v)==len(a) else ('no_eligible_object_nights' if len(a)==0 else 'not_evaluated')))
 pd.DataFrame(coverage).to_csv(TABLES/'retention_processing_coverage.csv',index=False)
 k=int(e.strict_linked.sum());n=len(e)
 return dict(eligible_all=len(eligible),evaluated=n,not_evaluated=len(eligible)-n,linked=k,fraction=k/n,wilson95=wilson(k,n),
   stage_counts={c:int(e[c].sum()) for c in ['strict_linked','fit_ok_survived','is_good_survived']},ambiguous_detection_keys=ambiguous,
   nights=sum(c['status'] in ['archived','reconstructed'] for c in coverage),eligible_nights=e.noon_night.nunique(),exposures=sum(c['science_exposures'] for c in coverage if c['status'] in ['archived','reconstructed']))


def orbit_class(a,e):
 q=a*(1-e)
 return np.select([q<1.3,q<1.666,a<2.06,a<3.28,a<5.05,a<5.4,a>=5.4],
  ['NEO (q < 1.3 AU)','Mars-crosser','Inner belt / Hungaria','Main belt','Outer belt / Hilda','Jupiter Trojan region','Distant (a ≥ 5.4 au)'],default='Other')


def objects(q):
 orbit_path=INPUT/'frozen_orbit_elements.csv.gz'
 if not orbit_path.exists():orbit_path=INPUT/'annual_used_orbit_elements.csv.gz'
 orbit=pd.read_csv(orbit_path,dtype={'num':str},low_memory=False)
 orbit['num']=orbit.num.str.strip();orbit['normname']=orbit.name.fillna('').str.replace(r'\s+','',regex=True).str.lower()
 o=q.groupby('object_key').agg(asteroid_number=('asteroid_number','first'),asteroid_name=('asteroid_name','first'),detections=('epoch_mjd','size'),nights=('noon_night','nunique'),first=('epoch_mjd','min'),last=('epoch_mjd','max')).reset_index()
 o['span_days']=o['last']-o['first'];o['normname']=o.asteroid_name.fillna('').str.replace(r'\s+','',regex=True).str.lower()
 bynum=orbit.loc[orbit.num.notna()&orbit.num.ne('')].drop_duplicates('num').set_index('num');byn=orbit.drop_duplicates('normname').set_index('normname')
 for source,target in [('a','a_au'),('e','ecc'),('incl','inc_deg'),('epoch','catalog_epoch')]:
  o[target]=o.asteroid_number.map(bynum[source]).combine_first(o.normname.map(byn[source]))
 o['orbit_class']=orbit_class(o.a_au,o.ecc);o['orbit_key']=o.object_key
 o.to_csv(DATA/'fig07_objects_orbits_revisits.csv.gz',index=False)
 used=orbit.loc[orbit.num.isin(o.asteroid_number)|orbit.normname.isin(o.normname)]
 used.to_csv(INPUT/'annual_used_orbit_elements.csv.gz',index=False)
 summary=o.groupby('orbit_class').agg(unique_objects=('object_key','size'),detections=('detections','sum')).reset_index();summary['cohort']='quality'
 summary.rename(columns={'orbit_class':'class'}).to_csv(TABLES/'orbital_classes.csv',index=False)
 assert o.a_au.notna().all(),o.loc[o.a_au.isna(),'object_key'].head().tolist()
 return o


def main():
 ledger=pd.read_csv(TABLES/'exposure_ledger.csv.gz',dtype={'night':str,'noon_night':str},low_memory=False)
 states=json.loads((INPUT/'annual_archive_states.json').read_text())
 ds={}
 for name in ['known_matched','known_mask15']:
  old=pd.read_parquet(OLD/'inputs'/f'{name}.parquet');extra=pd.read_parquet(INPUT/f'{name}_incremental.parquet')
  all=pd.concat([old,extra],ignore_index=True)
  d=attach(all,ledger);d.loc[~d.quality].groupby(['night','quality_exclusion'],dropna=False).size().reset_index(name='associations').to_csv(TABLES/f'{name}_excluded.csv',index=False)
  ds[name]=d
  # Keep compact raw columns for portable replay; ledger supplies the shared selection.
  all.to_parquet(INPUT/f'{name}.parquet',index=False,compression='zstd')
 q=ds['known_matched'].loc[ds['known_matched'].quality].copy();m=ds['known_mask15']
 print('Known annual primary',len(q),'identities',q.object_key.nunique(),flush=True)
 obj=objects(q);pairs,audit=exposure_pairs(q)
 geom=geometry(q);histogram(geom.elongation_deg,np.linspace(0,180,61),'figA1_elongation.csv')
 for key in ['predicted','measured']:histogram(pairs[f'{key}_rate_arcsec_h'],np.r_[0,np.geomspace(.1,10000,76),np.inf],f'figA1_{key}_rate.csv')
 quality_tables(ds['known_matched'],m)
 sky=SkyCoord(q.obs_ra_win_deg.to_numpy()*u.deg,q.obs_dec_win_deg.to_numpy()*u.deg).transform_to(BarycentricMeanEcliptic(equinox=Time('J2000')))
 pix=hp.ang2pix(128,sky.lon.deg,sky.lat.deg,lonlat=True)
 pd.Series(pix).value_counts().sort_index().rename_axis('ring_pixel').reset_index(name='known_detections').to_csv(DATA/'fig03_known_ecliptic_nside128.csv',index=False)
 rr=retention(m,ledger,states)
 result=dict(known=dict(qualified_detected=len(q),qualified_identified=q.object_key.nunique(),qualified_detected_nights=q.noon_night.nunique(),
  revisits={c:describe(obj[c]) for c in ['detections','nights','span_days']}),retention=rr,motion_pair_audit=audit,
  predicted_rate_arcsec_h=describe(pairs.predicted_rate_arcsec_h),measured_rate_arcsec_h=describe(pairs.measured_rate_arcsec_h),
  corrected_elong_quality_deg=describe(geom.elongation_deg),astrometry={c:describe(q[c]) for c in ['dra_arcsec','ddec_arcsec','radial_arcsec']},
  l2_quality_header_metrics={c:describe(ledger.loc[ledger.quality,c].where(ledger.loc[ledger.quality,c].gt(0)) if c=='limiting_mag' else ledger.loc[ledger.quality,c]) for c in ['limiting_mag','seeing_arcsec','ccd_zeropoint','ccd_zeropoint_rms']})
 (TABLES/'performance_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
