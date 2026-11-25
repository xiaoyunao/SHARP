# SHARP Figures 2, 5, and 6: publication revision

## Scope and source

The source is `paper_revision_20261015` at SHARP commit
`2f5f71d5a074c4cfd712b0e3fbef240c40745e2a`. The scientific cutoff remains
2026-10-15. This directory only re-expresses frozen products; it does not run
the production pipeline, search for new objects, fit orbits, or repeat the
previous science validation exercises.

The source scripts were read to establish the selections and definitions:
`known_performance.py`, `replay_frozen_recovery.py`, `stellar_repeatability.py`,
and `exposures.py`. The plotting script does not import or execute them.
SHA-256 hashes of the actual input files are in `summary.json`.

## Deliverable files

- `fig02_workflow.pdf`, `fig02_workflow.png`, `fig02_workflow.svg`:
  compact common-reduction column and two analysis branches, preserving the
  actual search -> known-object veto -> image-review sequence.
- `fig05_measurement_quality.pdf`, `fig05_measurement_quality.png`:
  model depth; raw aperture-error counts and frozen trends; raw radial
  astrometric counts and frozen trends; nightly signed offsets with count support.
- `fig06_conditional_retention.pdf`, `fig06_conditional_retention.png`:
  historical eligible-object-night density; same-grid conditional retention;
  frozen one-dimensional fractions with Wilson intervals and denominator counts.
- `redraw_figures.py`: reproducible preparation and drawing code.
- `data/`: compact NPZ/CSV plotting data, original frozen trend rows, workflow
  nodes/edges, and explicit lists of the magnitude exclusions.
- `summary.json`: exact sample counts, cell widths, definitions, and file hashes.
- `captions.md`: English replacement captions and synchronization table.
- `requirements.txt`: portable Python dependencies.
- `pdf_validation.json`: page geometry, embedded-font and vector-output checks.

All PDFs are 6.6 inches wide. Figure 2 is 2.88 inches tall; Figures 5 and 6 are
5.85 and 6.15 inches tall. Base labels and annotations are 8--9.2 pt at this size.
The PNGs are 320 dpi. Text, lines, cells, and color bars are vector objects in
the PDFs; the text uses embedded TrueType fonts. Statistical figures use the
manuscript's serif font preference (Times New Roman, Liberation Serif, DejaVu
Serif) with STIX math, while the workflow uses DejaVu Sans.

## Redraw from the included compact data

```sh
python -m pip install numpy pandas matplotlib
python redraw_figures.py
```

This mode does not require the original Parquet files or `pyarrow`. The output
defaults to the script directory, making the folder portable. No absolute
workspace path or local `_runtime` directory is used by the script.

## Prepare again from the same frozen source bundle

```sh
python -m pip install -r requirements.txt
python redraw_figures.py --prepare --source /path/to/paper_revision_20261015
```

This optional mode reads the existing frozen files and regenerates only compact
plotting arrays/tables. It joins the primary associations to the frozen exposure
ledger, computes spherical O-C residuals from the already archived coordinates,
and checks that bin counts and medians reproduce the frozen summary tables.
It never makes catalog queries or writes into the source bundle.

## Figure-specific data definitions

**Figure 5.** `quality == True` in the frozen exposure ledger selects 39,590
exposures and 531,612 primary 1-arcsec associations. Aperture-error data additionally
require `flag == 0`, `magerr_aper4 > 0`, `flux_aper4 > 0`, and `fluxerr_aper4 > 0`,
giving 500,092 records. The error density uses 0.1-mag by 0.04-dex cells and includes
all selected records. The radial density uses 0.1-mag by 0.02-arcsec cells; its
display interval, 9.5 <= g_aper < 22, contains 531,416 records. The omitted 196
records are listed in `fig05_c_outside_magnitude_range.csv`: 195 have magnitude
99 and one has magnitude 23.1787414551. The nightly offsets/counts use all 531,612
records, grouped by the ledger's local-noon night, for 128 nights. Gray count
bars are sample support, not errors on the nightly median.

**Figure 6.** Use the archived `eligible_ge3_and_rate` boolean directly; the
102,347 historical archive-night/object identities and 75,899 `strict_linked`
cases are unchanged. No exposure-quality/noon-night regrouping is substituted.
The 32 out-of-range magnitude records, including two retained cases, are listed
in `fig06_magnitude_excluded_32.csv`. They remain in the rate/overall samples.
The magnitude sample has N=102,315 and k=75,897. The 2D grid has 0.25-mag by
3-arcsec/h cells; retention cells with N<30 are masked. The count map is not
masked at N=30. The one-dimensional bins and Wilson intervals are the original
frozen values. Their lower histograms show N per bin, not N divided by bin width.

## Instructions for subsequent Codex edits

1. Start from `redraw_figures.py` and the included compact tables.
2. Retain the fixed source commit, historical sample definitions, `g_aper`
   notation, and the distinct overall versus magnitude denominators.
3. Keep the independent moving-source sequence before the known-object veto
   and image review. Do not add transient detection, new production capabilities,
   later observations, or later code changes.
4. Do not turn binary retention flags into a continuous scatter/KDE density.
   Counts refer to actual object-night positions in magnitude/rate space;
   retention is an unsmoothed cell ratio. Do not fill low-count masked cells.
5. At 6.6-inch output width, keep base fonts >=8 pt, moderate line weights,
   and unobstructed arrows, data, tick labels, and color bars.
6. After any drawing change, rerender all modified PDFs/PNGs and inspect them.
   Update the matching captions if a plotted definition changes.

Exclude `_runtime/`, `__pycache__/`, and temporary QA renderings from handoff
archives. They are local tooling rather than figure deliverables.


## Repository integration (2026-11-25)

Run from `paper_revision_20261015/`:

```sh
python figure_revision/redraw_figures.py --prepare --source .
python scripts/finalize_revision_r2.py --source .
cd manuscript
tectonic --keep-logs manuscript.tex
```

The first command rebuilds compact data and all three figures from the existing frozen inputs (completed in this delivery). Omit `--prepare --source .` to use the distributed compact tables. The finalizer synchronizes figure copies and updates current hashes, R2 number sources and validation; it does not repeat the earlier five science audits. Run it again after compiling to refresh the PDF page count. `scripts/make_figures.py` also invokes the R2 plotter for Figures 2, 5 and 6; its legacy layout functions are retained only as historical source. PDF bytes may vary with font availability and creation metadata; the exact generated hashes are recorded in `tables/manuscript_figure_identity.csv`.
