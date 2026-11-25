#!/usr/bin/env python3
"""Single-source annual paper numbers, synchronized figure copies, and delivery checks."""
import hashlib,json,re,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from pypdf import PdfReader
from PIL import Image
from exposures import ROOT,INPUT,TABLES,DATA

STAT=['fig03_coverage_and_known_sky','fig04_activity_and_cadence','fig05_measurement_quality','fig06_conditional_retention','fig07_orbits_and_revisits','figA1_motion_and_solar_elongation','figA2_candidate_short_arc_diagnostic']


def main():
 p=json.loads((TABLES/'performance_summary.json').read_text());c=json.loads((TABLES/'candidate_summary.json').read_text());plot=json.loads((ROOT/'figure_templates/summary.json').read_text())
 l=pd.read_csv(TABLES/'exposure_ledger.csv.gz',dtype={'night':str,'noon_night':str},low_memory=False);q=l.loc[l.quality]
 cad=pd.read_csv(TABLES/'field_night_cadence.csv');adj=pd.read_csv(TABLES/'adjacent_visit_intervals.csv.gz');g=pd.read_csv(TABLES/'geometric_coverage.csv').set_index('nside').loc[1024]
 cls=pd.read_csv(TABLES/'orbital_classes.csv').set_index('class');r=p['retention'];k=p['known'];dep=p['l2_quality_header_metrics']['limiting_mag'];seeing=p['l2_quality_header_metrics']['seeing_arcsec'];astro=p['astrometry']
 ac=adj.adjacent_min.dropna().quantile([.16,.5,.84]);span=cad.loc[cad.exposures.ge(3),'span_min'].quantile([.16,.5,.84])
 v={};src={}
 def add(name,value,source):v['yr'+name]=value;src['yr'+name]=source
 for n,x in [('Nexp',len(q)),('Nnight',q.noon_night.nunique()),('Hours',q.actual_exposure_s.sum()/3600)]:add(n,x,'tables/exposure_ledger.csv.gz:quality')
 for n,x in [('Area',g.geometric_once_deg2),('AreaThree',g.geometric_ge3_same_night_deg2)]:add(n,x,'tables/geometric_coverage.csv:nside=1024')
 for n,x in [('NGroups',len(cad)),('NOne',int(cad.exposures.eq(1).sum())),('NTwo',int(cad.exposures.eq(2).sum())),('NThree',int(cad.exposures.ge(3).sum())),('FracThree',cad.exposures.ge(3).mean()*100),('NAdjacent',len(adj.adjacent_min.dropna())),('Cadence',ac.loc[.5]),('CadenceLow',ac.loc[.16]),('CadenceHigh',ac.loc[.84]),('Span',span.loc[.5]),('SpanLow',span.loc[.16]),('SpanHigh',span.loc[.84])]:add(n,x,'tables/field_night_cadence.csv; tables/adjacent_visit_intervals.csv.gz')
 for label,d in [('Depth',dep),('FWHM',seeing)]:
  for suffix,key in [('', 'median'),('Low','p25'),('High','p75')]:add(label+suffix,d[key],'tables/performance_summary.json:l2_quality_header_metrics')
 for n,x in [('NAssoc',k['qualified_detected']),('NKnown',k['qualified_identified']),('NKnownNights',k['qualified_detected_nights']),('NMainBelt',int(cls.loc['Main belt','unique_objects'])),('MainBeltFrac',cls.loc['Main belt','unique_objects']/k['qualified_identified']*100),('NDistant',int(cls.loc['Distant (a ≥ 5.4 au)','unique_objects']))]:add(n,x,'tables/performance_summary.json:known; tables/orbital_classes.csv')
 for n,x in [('OC',astro['radial_arcsec']['median']),('RAOffset',astro['dra_arcsec']['median']*1000),('DecOffset',astro['ddec_arcsec']['median']*1000),('RAScatter',astro['dra_arcsec']['sigma_iqr']),('DecScatter',astro['ddec_arcsec']['sigma_iqr'])]:add(n,x,'tables/performance_summary.json:astrometry')
 for n,x in [('NRetNights',r['nights']),('NRetExposures',r['exposures']),('NEligible',r['evaluated']),('NLinked',r['linked']),('Retention',r['fraction']*100),('RetentionLow',r['wilson95'][0]*100),('RetentionHigh',r['wilson95'][1]*100)]:add(n,x,'tables/performance_summary.json:retention; tables/annual_object_nights.parquet')
 for n,x in [('MedAssoc',k['revisits']['detections']['median']),('MedNights',k['revisits']['nights']['median']),('MaxNights',k['revisits']['nights']['max']),('MedBaseline',k['revisits']['span_days']['median']),('Elongation',p['corrected_elong_quality_deg']['median']),('ElongationMin',p['corrected_elong_quality_deg']['min']),('NRatePairs',p['measured_rate_arcsec_h']['n']),('MeasuredRate',p['measured_rate_arcsec_h']['median']),('PredictedRate',p['predicted_rate_arcsec_h']['median'])]:add(n,x,'tables/performance_summary.json')
 for n,key in [('NAuto','automatic'),('NSearchNights','search_nights'),('NReviewed','reviewed'),('NReviewNights','review_nights'),('NCandidates','retained'),('NCandidateMeasurements','measurements'),('NCandidateNights','candidate_nights'),('NCometLinks','comet_links'),('NUnidentified','unidentified')]:add(n,c[key],'tables/candidate_summary.json; tables/annual_candidate_states.parquet')
 # Optional identity-level motion classifications are not asserted for the annual sample.
 for n in ['NDistantMotion','NDistantStatic','NDistantUnresolved']:add(n,None,'not used in annual prose; prior local audit remains historical evidence')
 add('NDepth',dep['n'],'tables/performance_summary.json:l2_quality_header_metrics.limiting_mag.n')
 add('NDepthUnavailable',len(q)-dep['n'],'tables/performance_summary.json:depth_invalid_exposures')
 fields=pd.read_csv(ROOT/'annual_results_fields.csv');rows=[];numbers=['% Generated from annual_results.json; display only, never calculation inputs.']
 for f in fields.itertuples():
  value=v[f.macro]
  if value is None:display='not evaluated'
  elif f.display_format=='integer':display=f'{int(value):,}'
  else:display=f'{value:.2f}'
  rows.append(dict(macro=f.macro,quantity=f.quantity,annual_value=value,display=display,source_field=src[f.macro]))
  numbers.append('\\newcommand{\\'+f.macro+'}{'+display+'}')
 (ROOT/'manuscript/annual_numbers.tex').write_text('\n'.join(numbers)+'\n')
 pd.DataFrame(rows).to_csv(TABLES/'annual_number_source_map.csv',index=False)
 result=dict(interval=dict(start_utc_inclusive='2025-11-15T04:00:00Z',end_utc_exclusive='2026-11-15T04:00:00Z',timezone='Asia/Shanghai',night_boundary_local='12:00',calendar_nights=365),
  macros=v,performance=p,candidates=c,plot_selection=plot,scientific_methods_commit='d2f0057bd4cdd4c7a3c6c2431f6449a72c327284',prior_input_commit='1e278897745dc9d879cb2b0e2da0a6687cb843f1',
  catalog=dict(name='Lowell ASTORB',sha256='ea22c7062e7349d9add34fb64835916d66705ad60a717cc7c86a0c47f5aa3040'),optional_annual_distant_motion_classification='not undertaken; no identity-level claims in main paper')
 (ROOT/'annual_results.json').write_text(json.dumps(result,indent=2)+'\n')
 for stem in ['fig05_measurement_quality','fig06_conditional_retention']:
  for ext in ['pdf','png']:shutil.copy2(ROOT/'figure_templates'/f'{stem}.{ext}',ROOT/'figures'/f'{stem}.{ext}')
 for stem in ['fig01_exposure','fig02_workflow','fig08_candidate']:
  # Examples and workflow are the author-approved R3 versions.
  shutil.copy2(ROOT/'manuscript/figures'/f'{stem}.pdf',ROOT/'figures'/f'{stem}.pdf')
 for stem in STAT:
  for ext in ['pdf']:
   pth=ROOT/'figures'/f'{stem}.{ext}';assert pth.exists();shutil.copy2(pth,ROOT/'manuscript/figures/annual'/pth.name)
 assert len(q)==q.exposure_key.nunique();assert q.noon_night.str.fullmatch(r'\d{8}').all();assert v['yrNLinked']<=v['yrNEligible'];assert r['not_evaluated']==0
 month=pd.read_csv(TABLES/'annual_monthly_counts.csv');assert len(month)==13 and int(month.qualified.sum())==len(q)
 assert int(month.observing_nights.sum())==v['yrNnight']
 assert len(PdfReader(ROOT/'figures/fig06_conditional_retention.pdf').pages)==1
 (ROOT/'manuscript/annual_values.tex').write_text('% Generated annual values; see ../annual_results.json.\n\\newif\\ifannualresultsready\n\\annualresultsreadytrue\n\\input{annual_numbers}\n')
 checks=dict(exposure_keys_unique=True,qualified_exposures=len(q),physical_nights=v['yrNnight'],monthly_labels=month.month.tolist(),monthly_counts_close=True,retention_numerator=v['yrNLinked'],retention_denominator=v['yrNEligible'],retention_stage_counts_equal=len(set(r['stage_counts'].values()))==1,unreviewed_or_completion_not_evidenced=c['review_completion_not_evidenced'],prior_retained_preserved=c['prior_retained_in_annual'],seven_annual_figures_present=True)
 (TABLES/'annual_validation.json').write_text(json.dumps(checks,indent=2)+'\n')
 print(json.dumps(v,indent=2))

if __name__=='__main__':main()
