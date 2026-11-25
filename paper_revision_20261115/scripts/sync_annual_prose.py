#!/usr/bin/env python3
"""Apply result-specific selections to the supplied editorial prose, preserving its narrative."""
from pathlib import Path
import json
import pandas as pd
r=Path(__file__).resolve().parents[1]
c=json.loads((r/'tables/candidate_summary.json').read_text());s=json.loads((r/'figure_templates/summary.json').read_text())
p=r/'manuscript/sections/05_initial_samples.tex';t=(r/'editorial_source/sections'/p.name).read_text();a=t.index('For the sparsely sampled distant region,');b=t.index('The median identity',a)
t=t[:a]+'''The distant region contains \\yrNDistant\\ positional-association
identities. The orbital distribution retains these sparse associations;
for such slowly moving sources, a positional coincidence alone need not
establish motion. Detailed local checks of selected earlier associations
are provided separately, rather than treating all distant-region points
as independently confirmed moving objects.

'''+t[b:]
t=t.replace('The independent search produces \\yrNAuto\\ automatic links', 'Restricting archived independent-search links to the science exposures\nyields \\yrNAuto\\ automatic links')
t=t.replace('Image inspection covers\n\\yrNReviewed\\ links from \\yrNReviewNights\\ nights and retains', 'Records with explicit image-review or completion evidence cover\n\\yrNReviewed\\ links from \\yrNReviewNights\\ nights and retain')
t=t.replace('(Table~\\ref{tab:candidates}). Each retained link contains',f'''(Table~\\ref{{tab:candidates}}). A further {c['review_completion_not_evidenced']:,} links from
{c['unreviewed_nights']} nights have no explicit review-completion evidence and are listed
separately; they are not classified as rejected. The retained count
combines the previously audited subset with {c['new_retained']} additional positive
image-review records that satisfy the same exposure selection and excludes
subsequent explicit false-positive decisions. Each retained link contains''')
t=t.replace('The automatic and image-inspected samples have\nthe observing coverage given in the text.', 'The automatic and image-inspected samples have\nthe observing coverage given in the text. The retained sample is the\nreviewed subset, and the full-year review is not yet complete.')
t=t.replace('the retained\nfirst-year sample', 'the retained\nreviewed subset from the first year');p.write_text(t)
p=r/'manuscript/sections/04_observing_performance.tex';t=(r/'editorial_source/sections'/p.name).read_text();t=t.replace('six field--night sequences and compare', 'six field--night sequences selected between 2026 February and May and compare')
t=t.replace('using bins with at least 50 stellar pairs.', 'selected between 2026 February and May, using bins with at least 50 stellar pairs.')
t=t.replace('encircled energy at FWHM/2. Its median value is', f'encircled energy at FWHM/2. There are {s["fig05"]["depth_exposures"]:,} valid model-depth estimates; {s["fig05"]["depth_unavailable_exposures"]} header placeholders are treated as unavailable for this statistic. Its median value is')
t=t.replace('Cells are 0.1~mag by 0.04~dex.', f'''The displayed sample contains
{s['fig05']['photometry_displayed_n']:,} measurements with usable magnitudes and uncertainties
within the plotted ranges; the same selection defines the orange curves.
Cells are 0.1~mag by 0.04~dex.''')
t=t.replace('line and band mark the median and interquartile range.', f'line and band mark the median and interquartile range of {s["fig05"]["depth_exposures"]:,} valid estimates. The {s["fig05"]["depth_unavailable_exposures"]} unavailable header values are excluded only from this metric.')
t=t.replace('The magnitude-dependent panels use measurements with usable\n\\gaper\\ in the displayed range.', f'''The radial panel contains {s['fig05']['astrometry_displayed_n']:,} associations;
{s['fig05']['astrometry_magnitude_outside_n']:,} magnitude values lie outside its displayed range.
The nightly panel uses all \\yrNAssoc\\ associations.''')
t=t.replace('Both maps use unsmoothed counts.', f'''The magnitude-dependent panels contain
{s['fig06']['magnitude_eligible_n']:,} eligible object--nights and
{s['fig06']['magnitude_retained_n']:,} linked cases. The
{s['fig06']['magnitude_excluded_n']} excluded magnitudes include
{s['fig06']['magnitude_excluded_retained_n']} linked cases. Both maps use unsmoothed counts.''');p.write_text(t)
p=r/'manuscript/sections/07_discussion_and_conclusions.tex';t=(r/'editorial_source/sections'/p.name).read_text();t=t.replace('retained first-year candidates','retained candidates in the reviewed first-year subset');p.write_text(t)
p=r/'manuscript/manuscript.tex';t=(r/'editorial_source/manuscript.tex').read_text().replace('% R2 methods and figure designs are retained. Annual measurements are pending.','% R2 methods and figure designs retained; annual values generated from annual_results.json.');p.write_text(t)
l=pd.read_csv(r/'tables/exposure_ledger.csv.gz',dtype={'night':str,'noon_night':str},low_memory=False)
x=l.loc[l.in_year].copy();x['month']=x.noon_night.str[:6];x.sort_values(['quality','canonical_directory_copy'],ascending=False,inplace=True);x=x.drop_duplicates('exposure_key')
mode=x.groupby(['month','observation_mode','filter_reconciled']).agg(exposures=('file_name','size'),hours=('actual_exposure_s',lambda v:v.sum()/3600)).reset_index();mode.to_csv(r/'tables/observing_mode_counts.csv',index=False)
m=pd.read_csv(r/'tables/annual_monthly_counts.csv',dtype={'month':str});exc=x.loc[x.standard_mode & ~x.quality].groupby('month').quality_exclusion.agg(lambda z:json.dumps(z.value_counts().to_dict(),sort_keys=True));m['standard_exclusion_reasons']=m.month.map(exc).fillna('{}');m['excluded_standard_exposures']=m.raw-m.qualified;m.to_csv(r/'tables/annual_monthly_counts.csv',index=False)
