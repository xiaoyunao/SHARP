#!/usr/bin/env python3
"""Final annual set, figure, source and manuscript consistency checks."""
import hashlib,json,logging,re
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader
from exposures import ROOT,INPUT,TABLES,DATA
logging.getLogger('pypdf').setLevel(logging.ERROR)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
 a=json.loads((ROOT/'annual_results.json').read_text());v=a['macros'];q=pd.read_csv(TABLES/'exposure_ledger.csv.gz',dtype={'night':str,'noon_night':str},low_memory=False);q=q.loc[q.quality]
 assert q.exposure_key.is_unique and len(q)==v['yrNexp'] and q.noon_night.nunique()==v['yrNnight']
 times=pd.to_datetime(q.mid_utc,format='mixed',utc=True);assert times.ge(pd.Timestamp(a['interval']['start_utc_inclusive'])).all() and times.lt(pd.Timestamp(a['interval']['end_utc_exclusive'])).all()
 month=pd.read_csv(TABLES/'annual_monthly_counts.csv');assert len(month)==13 and month.qualified.sum()==len(q)
 nights=pd.read_csv(DATA/'fig04_nightly_activity.csv');assert len(nights)==365 and nights.qualified.sum()==len(q)
 assert pd.read_csv(DATA/'fig03_known_ecliptic_nside128.csv').known_detections.sum()==v['yrNAssoc']
 cls=pd.read_csv(TABLES/'orbital_classes.csv');assert cls.unique_objects.sum()==v['yrNKnown'] and cls.detections.sum()==v['yrNAssoc']
 o=pd.read_parquet(TABLES/'annual_object_nights.parquet');e=o.loc[o.eligible_ge3_and_rate];assert e.evaluated.all() and len(e)==v['yrNEligible'] and e.strict_linked.sum()==v['yrNLinked']
 assert e.strict_linked.eq(e.fit_ok_survived).all() and e.strict_linked.eq(e.is_good_survived).all()
 with np.load(ROOT/'figure_templates/data/density_grids.npz') as g:
  assert g['ret_n'].sum()==a['plot_selection']['fig06']['magnitude_eligible_n'] and g['ret_k'].sum()==a['plot_selection']['fig06']['magnitude_retained_n']
  assert (g['ret_k']<=g['ret_n']).all()
  assert pd.read_csv(ROOT/'figure_templates/data/fig05_a_model_depth_counts.csv')['count'].sum()==v['yrNDepth']
  assert v['yrNDepth']+v['yrNDepthUnavailable']==v['yrNexp']
  assert g['photo_n'].sum()==a['plot_selection']['fig05']['photometry_displayed_n']
  assert g['astro_n'].sum()==a['plot_selection']['fig05']['astrometry_displayed_n']
 cand=pd.read_csv(TABLES/'retained_linkages.csv');assert len(cand)==v['yrNCandidates'] and cand.science_link.all() and cand.prior_retained.sum()==58
 countfiles=[]
 for p in sorted((ROOT/'figures').glob('*.pdf')):
  assert len(PdfReader(p).pages)==1
  image=p.with_suffix('.png');assert image.exists()
  with Image.open(image) as im:assert im.width>1000 and im.height>500
  target=ROOT/'manuscript/figures'/('annual' if p.stem not in ['fig01_exposure','fig02_workflow','fig08_candidate'] else '')/p.name
  assert sha(target)==sha(p)
  countfiles.append(dict(file=p.name,pdf_sha256=sha(p),png_sha256=sha(image)))
 assert len(countfiles)==10
 fonts=[]
 for stem in ['fig02_workflow','fig05_measurement_quality','fig06_conditional_retention']:
  file=ROOT/'figures'/f'{stem}.pdf';reader=PdfReader(file);entries=[]
  for ref in reader.pages[0]['/Resources']['/Font'].get_object().values():
   obj=ref.get_object();desc=obj.get('/DescendantFonts',[obj])[0].get_object();fd=desc.get('/FontDescriptor',{}).get_object();embedded=any(k in fd for k in ['/FontFile','/FontFile2','/FontFile3']);assert embedded;entries.append(dict(name=str(obj.get('/BaseFont')),embedded=embedded))
  fonts.append(dict(file=file.name,bytes=file.stat().st_size,sha256=sha(file),fonts=entries))
 (ROOT/'figure_templates/pdf_validation.json').write_text(json.dumps(fonts,indent=2)+'\n')
 captions=[]
 for fn in ['03_moving_object_analysis.tex','04_observing_performance.tex']:
  source=(ROOT/'manuscript/sections'/fn).read_text()
  for body in re.findall(r'\\begin\{figure\}.*?\\end\{figure\}',source,re.S):
   if any(x in body for x in ['fig:workflow','fig:quality','fig:retention']):captions.append(body)
 (ROOT/'figure_templates/captions.md').write_text('# Annual captions (LaTeX; values from annual_numbers.tex)\n\n'+'\n\n'.join('```latex\n'+x+'\n```' for x in captions)+'\n')
 doc=PdfReader(ROOT/'manuscript/manuscript.pdf');text='\n'.join(p.extract_text() for p in doc.pages)
 for x in ['53,398','797,887','80,591','74.08','114','354']:assert x in text,x
 for x in ['Annual measurements awaiting calculation','Blue fields and reserved figure','39,590','74.14','DATA AND SOFTWARE AVAILABILITY']:assert x not in text,x
 assert (ROOT/'manuscript/sections/06_transient_sources.tex').read_text().strip()=='\\section{Transient Sources}\n\\label{sec:transients}'
 assert '\\annualresultsreadytrue' in (ROOT/'manuscript/annual_values.tex').read_text()
 numbers=(ROOT/'manuscript/annual_numbers.tex').read_text();assert not re.search(r'\d+\.\d{3,}',numbers)
 refs=pd.read_csv(INPUT/'annual_product_provenance.csv.gz');assert pd.to_datetime(refs.mtime_utc,utc=True).lt(pd.Timestamp(a['interval']['end_utc_exclusive'])).all()
 for rec in json.loads((ROOT/'frozen_methods/source_manifest.json').read_text()):assert sha(ROOT/'frozen_methods'/Path(rec['source']).name)==rec['sha256']
 previews=pd.read_csv(ROOT/'review_packet/review_preview_index.csv');assert len(previews)==a['candidates']['review_completion_not_evidenced'] and previews.preview_state.eq('existing_gif').all()
 for f in previews.sheet.unique():assert (ROOT/'review_packet'/f).exists()
 final=dict(all_checks_pass=True,science_exposures=len(q),physical_nights=v['yrNnight'],primary_associations=v['yrNAssoc'],eligible_object_nights=len(e),strict_linked=int(e.strict_linked.sum()),post_cutoff_archived_product_count=0,manuscript_pages=len(doc.pages),publication_figures=countfiles,review_previews=len(previews),review_contact_sheets=previews.sheet.nunique(),transient_title_only=True,two_decimal_statistics=True)
 (TABLES/'annual_delivery_validation.json').write_text(json.dumps(final,indent=2)+'\n')
 # Local TeX build logs and large catalog/kernel dependencies are not scientific deliverables.
 ignored={'de432s.bsp','frozen_orbit_elements.csv.gz','freeze_manifest.json'}
 manifest={str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and p.name not in ignored and '__pycache__' not in p.parts and p.suffix not in ['.log','.blg','.aux','.out'] and not p.name.endswith('_incremental.csv.gz') and not any(x.startswith('.') for x in p.relative_to(ROOT).parts)}
 (ROOT/'freeze_manifest.json').write_text(json.dumps(dict(method_commit=a['scientific_methods_commit'],prior_input_commit=a['prior_input_commit'],interval=a['interval'],files_sha256=manifest),indent=2)+'\n')
 print(f'Passed annual delivery: {len(q)} exposures, {len(e)} eligible nights, {len(countfiles)} figures, {len(doc.pages)} manuscript pages; {len(previews)} review previews')

if __name__=='__main__':main()
