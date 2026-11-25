#!/usr/bin/env python3
"""Extend the established exposure selection to a complete physical-night year."""
import json
from pathlib import Path
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.utils import iers
from exposures import (ROOT,INPUT,TABLES,DATA,CORNERS,iso_time,mjd_utc,noon_night,
                       aperture_altitude,cadence_pointings,run_coverage)

OLD=ROOT.parent/'paper_revision_20261015'
START=pd.Timestamp('2025-11-15T04:00:00Z');END=pd.Timestamp('2026-11-15T04:00:00Z')


def main():
 old=pd.read_csv(OLD/'tables/exposure_ledger.csv.gz',dtype={'night':str,'noon_night':str,'field_id':str},low_memory=False)
 old['noon_night']=old.noon_night.str.replace(r'\.0$','',regex=True)
 old['time_end_utc']=pd.to_datetime(old.start_utc,utc=True,format='mixed')+pd.to_timedelta(old.actual_exposure_s,unit='s')
 old['provenance']='1e278897 commissioning ledger, original selections preserved'
 s=pd.read_csv(INPUT/'incremental_headers.csv.gz',dtype={'directory_night':str},low_memory=False)
 a=s.loc[s.kind.eq('raw')].copy().rename(columns={'directory_night':'night','source_file':'file_name','filter':'raw_filter','bytes':'current_raw_bytes'})
 w=s.loc[s.kind.eq('l2_wcs')].copy();w['file_name']=w.source_file.str.replace(r'_cat\.fits(?:\.gz)?$','.fits',regex=True)
 w['l2_file']=w.source_file;w['wcs_verified']=w.error.isna()
 cols=['night','file_name','l2_file','wcs_verified','bytes','mtime_ns',*[c for c in w if c.startswith('l2_') and c not in ['l2_file','l2_present']]]
 w=w.rename(columns={'directory_night':'night'});w=w[[c for c in cols if c in w]].drop_duplicates(['night','file_name'])
 a=a.drop(columns=[c for c in a if c.startswith('l2_') and c!='l2_present'])
 a=a.merge(w,on=['night','file_name'],how='left',suffixes=('','_header'),validate='one_to_one')
 a['l2_available']=a.l2_file.notna();a['raw_start_utc']=iso_time(a.obs_date);a['l2_start_utc']=mjd_utc(a.l2_mjd)
 a['start_utc']=a.raw_start_utc.combine_first(a.l2_start_utc)
 a['actual_exposure_s']=pd.to_numeric(a.exposure,errors='coerce').combine_first(pd.to_numeric(a.l2_exposure,errors='coerce'))
 a['mid_utc']=a.start_utc+pd.to_timedelta(a.actual_exposure_s/2,unit='s');a['time_end_utc']=a.start_utc+pd.to_timedelta(a.actual_exposure_s,unit='s')
 a['noon_night']=noon_night(a.mid_utc);a['exposure_key']=a.file_name+'@'+a.start_utc.dt.floor('ms').astype('string')
 a['field_id']=a.file_name.str.extract(r'^OBJ_MP_(\d{4})_\d{4}\.fits$',expand=False)
 a['strict_name']=a.field_id.notna()
 a['filter_reconciled']=a.raw_filter.combine_first(a.l2_filter).astype('string').str.upper()
 a['observation_mode']=np.select([~a.strict_name,a.filter_reconciled.ne('W').fillna(True),~a.actual_exposure_s.between(29,31)],
  ['other_target_or_followup','other_filter','other_exposure_time'],default='standard_MP')
 a['l2_header_date_error_s']=(iso_time(a.l2_obs_date)-a.l2_start_utc).dt.total_seconds()
 a['raw_l2_time_error_s']=(a.raw_start_utc-a.l2_start_utc).dt.total_seconds()
 # Use actual WCS polygon centre (unit-vector average avoids RA=0 wrap).
 corners=a[[c for c in a if c in CORNERS]].to_numpy(float).reshape(-1,4,2)
 coords=SkyCoord(corners[:,:,0]*u.deg,corners[:,:,1]*u.deg)
 vec=coords.cartesian.xyz.value.mean(axis=2)
 c=SkyCoord(x=vec[0],y=vec[1],z=vec[2],representation_type='cartesian').spherical
 a['center_ra_deg']=c.lon.deg;a['center_dec_deg']=c.lat.deg
 for i in range(4):
  dist=coords[:,i].separation(coords[:,(i+1)%4]).deg
  a[f'side_{i}']=dist
 a['wcs_max_side_deg']=a[[f'side_{i}' for i in range(4)]].max(axis=1)
 with iers.earth_orientation_table.set(iers.IERS_A.open(str(INPUT/'paper_iers.all'))):aperture_altitude(a)
 a['quality_exclusion']=np.select([a.mid_utc.isna(),~a.observation_mode.eq('standard_MP'),~a.l2_available,
   a.wcs_verified.ne(True),a.field_altitude_deg.lt(0),~a.wcs_max_side_deg.between(2.80,3.10),
   a.raw_l2_time_error_s.abs().gt(1),a.l2_header_date_error_s.abs().gt(1)],
  ['missing_exposure_time','nonstandard_observing_mode','missing_l2','invalid_l2_wcs','field_below_horizon_at_header_time',
   'abnormal_WCS_detector_perimeter','raw_L2_start_disagreement_gt_1s','L2_ISO_MJD_disagreement_gt_1s'],default='qualified')
 a['quality']=a.quality_exclusion.eq('qualified')
 for source,target in [('mlim','limiting_mag'),('seeing','seeing_arcsec'),('fwhm','fwhm_header'),('ccdzp','ccd_zeropoint'),('ccdzprms','ccd_zeropoint_rms'),('sky','sky_adu'),('skyrms','sky_rms_adu')]:a[target]=pd.to_numeric(a['l2_'+source],errors='coerce')
 a['path']='/processed1/'+a.night+'/L2/'+a.l2_file.fillna('');a['provenance']='incremental existing FITS headers; no re-reduction'
 all=pd.concat([old,a],ignore_index=True)
 for col in ['mid_utc','start_utc','raw_start_utc']:all[col]=pd.to_datetime(all[col],format='mixed',utc=True)
 all['in_year']=all.mid_utc.between(START,END,inclusive='left')
 all['standard_mode']=all.observation_mode.eq('standard_MP') & all.filter_reconciled.eq('W')
 all.loc[~all.in_year,'quality_exclusion']='outside_annual_interval';all.quality=all.quality.eq(True)&all.in_year&all.standard_mode
 # Preserved commissioning duplicate selection wins; new duplicates prefer the physical directory.
 all['canonical_directory_copy']=all.night.eq(all.noon_night)
 good=all.loc[all.quality].sort_values(['canonical_directory_copy','l2_available'],ascending=False)
 duplicate=good.loc[good.exposure_key.duplicated()].index
 all.loc[duplicate,'quality']=False;all.loc[duplicate,'quality_exclusion']='noncanonical_duplicate'
 assert all.loc[all.quality,'exposure_key'].is_unique
 all.to_csv(TABLES/'exposure_ledger.csv.gz',index=False)
 all.loc[~all.quality].groupby(['night','quality_exclusion']).size().reset_index(name='exposures').to_csv(TABLES/'exposure_exclusions.csv',index=False)
 days=pd.date_range('2025-11-15','2026-11-14',freq='D').strftime('%Y%m%d')
 year=all.loc[all.in_year & all.standard_mode].sort_values(['quality','canonical_directory_copy'],ascending=False).drop_duplicates('exposure_key')
 act=year.groupby('noon_night').agg(raw=('file_name','size'),l2=('l2_available','sum'),qualified=('quality','sum')).reindex(days,fill_value=0).rename_axis('night').reset_index()
 act['month']=act.night.str[:6];act.to_csv(DATA/'fig04_nightly_activity.csv',index=False)
 months=pd.period_range('2025-11','2026-11',freq='M').strftime('%Y%m')
 monthly=act.groupby('month')[['raw','l2','qualified']].sum().reindex(months,fill_value=0)
 monthly['observing_nights']=act.assign(hasobs=act.qualified.gt(0)).groupby('month').hasobs.sum()
 monthly['hours']=all.loc[all.quality].assign(month=all.loc[all.quality,'noon_night'].str[:6]).groupby('month').actual_exposure_s.sum().div(3600)
 monthly['hours']=monthly.hours.fillna(0);monthly['archive_status']='existing raw directory inventory; no July16-Aug25 observations found'
 monthly.to_csv(TABLES/'annual_monthly_counts.csv',index_label='month');monthly.reset_index(names='month').to_csv(DATA/'fig04_monthly_activity.csv',index=False)
 print('MODES',all.loc[all.in_year].groupby(['observation_mode','filter_reconciled']).size().to_dict(),flush=True)
 print('QUALITY',all.loc[all.quality].shape,all.loc[all.quality].noon_night.nunique(),flush=True)
 footprint=Table.read(INPUT/'survey_fov_footprints_with_visibility.fits').to_pandas();footprint.field_id=footprint.field_id.map(lambda x:x.decode() if isinstance(x,bytes) else str(x))
 # The original coverage helper also writes activity tables; restore the annual physical-night tables afterwards.
 saved_night=(DATA/'fig04_nightly_activity.csv').read_bytes();saved_month=(DATA/'fig04_monthly_activity.csv').read_bytes()
 _,cad,intervals=run_coverage(all,all.loc[all.standard_mode & all.field_id.notna()],footprint)
 (DATA/'fig04_nightly_activity.csv').write_bytes(saved_night);(DATA/'fig04_monthly_activity.csv').write_bytes(saved_month)
 print('Coverage and cadence complete',len(cad),flush=True)

if __name__=='__main__':main()
