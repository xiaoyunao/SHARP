#!/usr/bin/env python3
"""Inspect per-month header validity and per-night known-product coverage."""
import json
import pandas as pd
from exposures import ROOT,INPUT,TABLES
from known_performance import describe
from build_annual_statistics import attach


def main():
 l=pd.read_csv(TABLES/'exposure_ledger.csv.gz',dtype={'night':str,'noon_night':str},low_memory=False)
 q=l.loc[l.quality].copy();q['month']=q.noon_night.str[:6];rows=[]
 for month,g in q.groupby('month'):
  for metric in ['limiting_mag','seeing_arcsec','ccd_zeropoint','ccd_zeropoint_rms']:
   values=g[metric].where(g[metric].ge(0) if metric=='ccd_zeropoint_rms' else g[metric].gt(0))
   rows.append(dict(month=month,metric=metric,science_exposures=len(g),unavailable=len(g)-values.notna().sum(),**describe(values)))
 pd.DataFrame(rows).to_csv(TABLES/'monthly_header_performance.csv',index=False)
 d=attach(pd.read_parquet(INPUT/'known_matched.parquet'),l);d=d.loc[d.quality]
 known=d.groupby('noon_night').agg(primary_associations=('object_key','size'),unique_identities=('object_key','nunique'))
 frame=q.groupby('noon_night').agg(science_exposures=('exposure_key','size'),directory_nights=('night',lambda x:'|'.join(sorted(set(x)))))
 frame=frame.join(known).fillna(dict(primary_associations=0,unique_identities=0));frame['primary_associations']=frame.primary_associations.astype(int);frame['unique_identities']=frame.unique_identities.astype(int)
 frame['known_state']='archived';frame.loc[frame.directory_nights.eq('20260605'),'known_state']='reconstructed_zero_matches'
 frame['source']='frozen commissioning tables or pre-endpoint nightly products; recovered June zero-match receipt'
 frame.to_csv(TABLES/'known_processing_coverage.csv',index_label='physical_night')
 assert frame.primary_associations.sum()==len(d)
 print('Monthly header rows',len(rows),'known physical nights',len(frame),'known associations',len(d))

if __name__=='__main__':main()
