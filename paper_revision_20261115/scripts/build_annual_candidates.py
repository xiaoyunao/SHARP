#!/usr/bin/env python3
"""Merge existing image-review decisions, preserving unreviewed and excluded outcomes."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.time import Time
from exposures import ROOT,INPUT,TABLES
from replay_frozen_recovery import normalize_obj_id
OLD=ROOT.parent/'paper_revision_20261015'


def main():
 led=pd.read_csv(TABLES/'exposure_ledger.csv.gz',dtype={'night':str,'noon_night':str},low_memory=False)
 u=pd.read_parquet(INPUT/'unknown_catalog_incremental.parquet');u.night=u.night.astype(str)
 old=pd.read_parquet(OLD/'inputs/unknown_catalog.parquet');old.night=old.night.astype(str)
 u=pd.concat([old,u.loc[u.night.gt('20260715')]],ignore_index=True)
 review=pd.read_parquet(INPUT/'review_incremental.parquet');review['review_positive']=pd.to_numeric(review.is_real,errors='coerce').eq(1)
 review=review.rename(columns={'tracklet_id':'trk_sub'}).drop_duplicates(['night','trk_sub'])
 u=u.merge(review[['night','trk_sub','review_positive']],on=['night','trk_sub'],how='left',validate='one_to_one')
 frozen=pd.read_csv(OLD/'tables/retained_linkages.csv',dtype={'night':str})
 status=pd.read_csv(OLD/'inputs/review_and_mpc_status.csv',dtype={'origin_night':str})
 retained_keys=set(zip(frozen.night,frozen.trk_sub))
 u['prior_retained']=[(n,t) in retained_keys for n,t in zip(u.night,u.trk_sub)]
 u['positive_in_current_review']=u.review_positive.eq(True)
 u['retained_review_decision']=u.prior_retained | (u.night.gt('20260715') & u.positive_in_current_review)
 actions=json.loads((INPUT/'candidate_manual_actions.json').read_text())
 rejected={key for a in actions['manual_actions'] if a.get('action')=='abandon_false_positive_sources' for key in a['source_ids']}
 u['later_manual_false_positive']=[n+':'+t in rejected for n,t in zip(u.night,u.trk_sub)]
 u['retained_review_decision'] &= ~u.later_manual_false_positive
 states=json.loads((INPUT/'annual_archive_states.json').read_text())
 # A bare default zero is not treated as a completed image review.
 completed={s['night'] for s in states if s.get('unknown_reply_exists') or any(k.endswith('_unknown_links_submit_stats.json') for k in s)}
 u['review_state']=np.select([u.prior_retained,u.later_manual_false_positive,u.positive_in_current_review,u.night.isin(completed)&u.review_positive.notna()],
  ['retained_prior_posthoc_audit','manual_false_positive','positive_image_review','completed_review_zero'],default='review_completion_not_evidenced')
 # Original source snapshot encodes reviewed positive/rejected decisions separately from default zeros.
 u.loc[u.night.le('20260715')&u.review_positive.notna()&u.night.isin(completed)&~u.prior_retained,'review_state']='completed_historical_review'
 raw=u.assign(image=u.image_names.str.split(';'),obj=u.objids.str.split(';')).explode(['image','obj'])
 raw['file_name']=raw.image.str.replace(r'\.fits(?:\.gz)?$','.fits',regex=True)
 raw=raw.merge(led[['night','file_name','exposure_key','quality','noon_night','mid_utc']],on=['night','file_name'],how='left',validate='many_to_one')
 raw['quality']=raw.quality.eq(True);keys=['night','linkage_id','trk_sub']
 membership=raw.groupby(keys).agg(members=('image','size'),qualified=('quality','sum'),distinct_exposures=('exposure_key','nunique'),physical_nights=('noon_night','nunique'),physical_night=('noon_night','first')).reset_index()
 membership['science_link']=membership.members.eq(membership.qualified)&membership.distinct_exposures.ge(3)&membership.physical_nights.eq(1)
 u=u.merge(membership,on=keys,validate='one_to_one')
 u['retained']=u.retained_review_decision & u.science_link
 selected=u.loc[u.retained].copy()
 members=raw.merge(selected[keys],on=keys,validate='many_to_one')
 members['midpoint_mjd']=Time(pd.to_datetime(members.mid_utc,format='mixed',utc=True).to_numpy(dtype='datetime64[ns]'),scale='utc').mjd
 spans=members.groupby(keys).midpoint_mjd.agg(lambda x:(x.max()-x.min())*1440).reset_index(name='arc_min_midpoint')
 selected=selected.merge(spans,on=keys,validate='one_to_one')
 known_identity=frozen[['night','trk_sub','is_c2025_y1','frozen_identity_status']]
 selected=selected.merge(known_identity,on=['night','trk_sub'],how='left',validate='one_to_one')
 selected['is_c2025_y1']=selected.is_c2025_y1.eq(True)
 selected['frozen_identity_status']=selected.frozen_identity_status.fillna('no_established_identity_in_available_record')
 selected.to_csv(TABLES/'retained_linkages.csv',index=False);members.to_csv(TABLES/'retained_measurements_midpoint.csv.gz',index=False)
 u.to_parquet(TABLES/'annual_candidate_states.parquet',index=False,compression='zstd')
 unreviewed=u.loc[u.science_link & u.review_state.eq('review_completion_not_evidenced')].copy()
 unreviewed.to_csv(TABLES/'candidate_review_queue.csv.gz',index=False)
 coverage=u.groupby('night').agg(archived_automatic=('linkage_id','size'),science_automatic=('science_link','sum'),positive_reviews=('positive_in_current_review','sum'),retained=('retained','sum')).reset_index()
 coverage.to_csv(TABLES/'candidate_processing_coverage.csv',index=False)
 covered=u.loc[u.science_link];reviewed=covered.loc[~covered.review_state.eq('review_completion_not_evidenced')]
 report=dict(archived_automatic=len(u),automatic=len(covered),search_nights=covered.physical_night.nunique(),
  reviewed=len(reviewed),review_nights=reviewed.physical_night.nunique(),retained=len(selected),measurements=len(members),candidate_nights=selected.physical_night.nunique(),
  comet_links=int(selected.is_c2025_y1.sum()),unidentified=int((~selected.is_c2025_y1).sum()),
  prior_retained_in_annual=int(selected.prior_retained.sum()),new_retained=int((~selected.prior_retained).sum()),
  review_completion_not_evidenced=len(unreviewed),unreviewed_nights=unreviewed.physical_night.nunique(),
  outside_science_membership=int((~u.science_link).sum()),manual_false_positive_excluded=int((u.later_manual_false_positive&u.positive_in_current_review).sum()))
 (TABLES/'candidate_summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
