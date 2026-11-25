#!/usr/bin/env python3
"""Publish/check R2 figure copies and update the R2 source map (no science rerun)."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
STEMS = ['fig02_workflow', 'fig05_measurement_quality', 'fig06_conditional_retention']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT, help='Paper support directory')
    args = parser.parse_args()
    root = args.source.resolve()
    revision = root/'figure_revision'
    summary = json.loads((revision/'summary.json').read_text())
    for f, record in summary['source_files'].items():
        assert sha(root/f) == record['sha256'], f
    for stem in STEMS:
        for ext in ['pdf', 'png']:
            shutil.copy2(revision/f'{stem}.{ext}', root/'figures'/f'{stem}.{ext}')
        shutil.copy2(revision/f'{stem}.pdf', root/'manuscript/figures'/f'{stem}.pdf')
    with np.load(revision/'data/density_grids.npz') as a:
        assert a['photo_n'].sum() == 500092
        assert a['astro_n'].sum() == 531416
        assert a['ret_n'].sum() == 102315 and a['ret_k'].sum() == 75897
        assert (a['ret_k'] <= a['ret_n']).all()
    bins = pd.read_csv(revision/'data/fig06_frozen_1d_bins.csv')
    for metric, n, k in [('measured_aperture_magnitude', 102315, 75897),
                         ('first_last_rate_arcsec_per_hour', 102347, 75899)]:
        rows = bins.loc[bins.metric.eq(metric)]
        assert rows.eligible_object_n.sum() == n and rows.survived_object_n.sum() == k
    identity = pd.read_csv(root/'tables/manuscript_figure_identity.csv')
    if 'scientific_source_commit' not in identity:
        identity['scientific_source_commit'] = identity.pinned_commit
    if 'generation' not in identity:
        identity['generation'] = 'archived figure; original scientific input unchanged'
    figures = []
    for idx, row in identity.iterrows():
        pdf = root/'figures'/row.pdf
        assert sha(pdf) == sha(root/'manuscript/figures'/row.pdf)
        assert len(PdfReader(pdf).pages) == 1
        identity.loc[idx, 'sha256'] = sha(pdf)
        if pdf.stem in STEMS:
            identity.loc[idx, 'pinned_commit'] = 'R2 derived figure; see generation and file SHA-256'
            identity.loc[idx, 'scientific_source_commit'] = summary['source_commit']
            identity.loc[idx, 'generation'] = 'figure_revision/redraw_figures.py --prepare --source <paper support directory>'
        figures.append({'file': row.pdf, 'sha256': sha(pdf)})
    identity.to_csv(root/'tables/manuscript_figure_identity.csv', index=False)
    mapping = pd.read_csv(root/'tables/manuscript_number_source_map.csv')
    mapping['manuscript_location'] = mapping.manuscript_location.str.replace('Sections 2, 4, 7', 'Sections 2, 4, 8', regex=False).str.replace('Sections 4, 7', 'Sections 4, 8', regex=False).str.replace('Sections 4, 5, 7', 'Sections 4, 5, 8', regex=False).str.replace('Appendix B', 'Appendix A.2', regex=False)
    mapping.loc[mapping.quantity.eq('Conditional connection retention'), 'manuscript_location'] = 'Section 4.3; Table 2; Fig. 6 (historical cohort)'
    mapping.loc[mapping.quantity.eq('Median radial O−C'), 'manuscript_location'] = 'Abstract; Section 4.2; Table 2; Fig. 5c'
    mapping.loc[mapping.quantity.eq('Historical radial O−C comparison'), 'manuscript_location'] = 'Supporting historical comparison only'
    for name in ['Median signed RA*cos(Dec) / Dec', 'Robust signed RA*cos(Dec) / Dec scatter']:
        mapping.loc[mapping.quantity.eq(name), 'manuscript_location'] = 'Section 4.2 (pooled statistics); Fig. 5d shows nightly medians'
    changes = {
        'Fig. 2': ('figure_revision/data/fig02_nodes.csv; figure_revision/data/fig02_edges.csv',
                   'known and independent motion paths; known-object veto follows short-arc consistency'),
        'Fig. 5': ('figure_revision/summary.json:fig05; figure_revision/data/fig05_*',
                   'a model depth; b 500092 valid-error records; c 531416 displayed radial residuals; d 531612 signed residuals over 128 nights'),
        'Fig. 6': ('figure_revision/summary.json:fig06; figure_revision/data/fig06_*',
                   'historical overall 75899/102347; magnitude 75897/102315; same-grid ratios; N<30 masked; 95% Wilson'),
    }
    for loc, (source, criteria) in changes.items():
        select = mapping.manuscript_location.eq(loc)
        mapping.loc[select, 'source_file_and_key'] = source
        mapping.loc[select, 'denominator_or_filter'] = criteria
        mapping.loc[select, 'generation_script'] = 'figure_revision/redraw_figures.py'
    extra = []

    def row(loc, name, unit, value, grain, criteria, cohort, source, generator):
        extra.append(dict(manuscript_location=loc, quantity=name, units=unit, value=value,
                          statistical_unit=grain, denominator_or_filter=criteria, cohort=cohort,
                          source_file_and_key=source, generation_script=generator))
    gen = 'scripts/finalize_revision_r2.py'
    method = 'inputs/calibration_method_note_20261124_r2.md; inputs/calibration_chain_evidence_r2.json'
    row('Section 2.1; Table 1', 'R2 detector and readout', 'pixels; micrometre; channels',
        '9216x9232; 10; 16', 'CCD', 'author-confirmed, 2026-11-24', 'initial campaign', 'AUTHOR_TODO.md', gen)
    row('Section 2.1', 'R2 Gaia reference and WCS propagation', 'release; Julian year',
        'EDR3; 2016.0 -> UTC exposure start Julian year', 'reference source',
        'pm_ra_cosdec and pm_dec; no parallax/RV; nonfinite result falls back', 'three existing frozen exposures', method,
        'scripts/check_calibration_chain_r2.py')
    row('Section 2.1', 'R2 photometric matching coordinate distinction', 'implementation',
        'W gaia_mags uses original ra/dec; SCAMP uses propagated X_WORLD/Y_WORLD', 'reference source',
        'actual consuming functions traced', 'current upstream version checked at author request', method,
        'scripts/check_calibration_chain_r2.py')
    row('Section 2.1; References', 'R2 in-preparation references', 'publication status',
        'N. Li et al. 2026; Y.-A. Xiao et al. 2026: in preparation', 'reference',
        'author-confirmed 2026-11-24; no invented DOI', 'R2 manuscript', 'manuscript/references.bib:Li2026InPrep,Xiao2026InPrep', gen)
    row('Abstract; Section 4.3; Table 2', 'R2 qualified-exposure retention', 'object-nights; percent',
        '75467 / 101784 = 74.144266', 'physical night / identity',
        'qualified exposures; historical processing-state policy retained', '126 physical nights',
        'tables/retention_reconciliation_summary.json; tables/retention_cohort_reconciliation.csv', 'scripts/reconcile_historical_cohorts.py')
    row('Section 4.3', 'R2 readable non-quarantined stage subset retention', 'object-nights; percent',
        '73787 / 99426 = 74.212983', 'physical night / identity',
        'further restrict to readable stage archive and no quarantine', '122 physical nights',
        'tables/retention_reconciliation_summary.json', 'scripts/reconcile_historical_cohorts.py')
    row('Sections 2.2 and 5.2', 'R2 retained candidate exposure validation', 'links; measurements',
        '58 / 179; all qualified', 'link / member exposure', 'at least 3 distinct exposures in same physical night',
        'retained historical candidate cohort', 'tables/linkage_cohort_reconciliation.csv; tables/linkage_cohort_members.csv',
        'scripts/reconcile_historical_cohorts.py')
    row('Section 5.1; Fig. 7 caption', 'R2 distant association vs identity classification', 'associations; identities',
        '79/9/31 associations; 7/7/3 identities; identity-group totals 107/9/3; 28 unresolved associations in supported identities',
        'association vs catalog identity', 'motion-supported / static-coincident / insufficient; do not mix grains',
        '119 associations to 17 identities', 'tables/distant_association_local_checks.csv; tables/rare_object_local_checks.csv',
        'scripts/review_distant_associations.py')
    row('Fig. 5b', 'R2 aperture-error density selection', 'measurements', 500092, 'primary association',
        'quality; flag=0; positive magerr_aper4, flux_aper4 and fluxerr_aper4', 'qualified exposures',
        'figure_revision/summary.json:fig05; figure_revision/data/fig05_b_counts.csv.gz', 'figure_revision/redraw_figures.py --prepare')
    row('Fig. 5c–d', 'R2 radial display and nightly signed counts', 'measurements; nights',
        '531416 displayed / 531612 total; 196 magnitude exclusions; 128 nights', 'association / physical night',
        'c: 9.5<=g_aper<22; d: all qualified primary associations', 'qualified exposures',
        'figure_revision/summary.json:fig05; figure_revision/data/fig05_c_outside_magnitude_range.csv', 'figure_revision/redraw_figures.py --prepare')
    row('Fig. 6a–c', 'R2 magnitude-panel retention denominator', 'object-nights',
        '75897 / 102315; excluded 2 / 32; overall and rate 75899 / 102347', 'historical archive-night / identity',
        '10<=median g_aper<20; retention grid N>=30; same-grid unsmoothed numerator/denominator',
        'historical eligible sample', 'figure_revision/summary.json:fig06; figure_revision/data/fig06_magnitude_excluded_32.csv',
        'figure_revision/redraw_figures.py --prepare')
    mapping = mapping.loc[~mapping.quantity.str.startswith('R2 ')].copy()
    pd.concat([mapping, pd.DataFrame(extra)], ignore_index=True).to_csv(root/'tables/manuscript_number_source_map.csv', index=False)
    sections = root/'manuscript/sections'
    subsection_counts = {fn: (sections/fn).read_text().count('\\subsection{') for fn in [
        '02_survey_and_data.tex', '03_moving_object_analysis.tex', '04_observing_performance.tex',
        '05_initial_samples.tex', '06_discussion.tex']}
    assert list(subsection_counts.values()) == [2, 0, 3, 2, 0]
    transient = (sections/'06_transient_sources.tex').read_text().strip()
    assert transient == '\\section{Transient Sources}\n\\label{sec:transients}'
    evidence = json.loads((root/'inputs/calibration_chain_evidence_r2.json').read_text())
    assert len(evidence['samples']) == 3
    assert all(s['tile_ref_epochs'] == [2016.0] for s in evidence['samples'])
    assert len(figures) == 10
    artifact = {'input_commit': summary['source_commit'], 'source_file_hashes_verified': True,
                'three_figures_prepared_from_frozen_inputs': True, 'fig05': summary['fig05'],
                'fig06': summary['fig06'], 'figure_pairs': figures,
                'subsection_counts': subsection_counts, 'transient_title_only': True,
                'calibration_exposures': 3, 'reference_examples': 18,
                'publication_pdf_pages': len(PdfReader(root/'manuscript/manuscript.pdf').pages)}
    (root/'tables/revision_r2_validation.json').write_text(json.dumps(artifact, indent=2)+'\n')
    print(f'Validated 10 figures, 3 R2 plotting selections, 3 calibration exposures, and {len(mapping)+len(extra)} source-map entries')


if __name__ == '__main__':
    main()
