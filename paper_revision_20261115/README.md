# SHARP first-year paper support

The deliverable is the first-year manuscript and its science support, covering UTC exposure midpoints **[2025-11-15 04:00, 2026-11-15 04:00)**. Physical nights use Beijing local noon; labels are 20251115–20261114. The statistical sample selects the established standard 30-second white-light MP observations. Other filters, exposure times and target modes are separately tabulated, including November 11–13. No production software, raw/L1/L2 files, cron, submission state or calibration zero points were changed.

- [Compiled manuscript and LaTeX](manuscript/) — 16 pages, eight main and two appendix figures; three supplied authors; title-only Transient Sources; combined Discussion and Conclusions; no Data and Software Availability section.
- [Annual results](annual_results.json) — single full-precision number source; [display/source map](tables/annual_number_source_map.csv); generated `manuscript/annual_numbers.tex` and `annual_values.tex`.
- [Exposure ledger](tables/exposure_ledger.csv.gz), [13-month table](tables/annual_monthly_counts.csv), [mode counts](tables/observing_mode_counts.csv), [exclusions](tables/exposure_exclusions.csv).
- [Performance](tables/performance_summary.json), [monthly header validity/distributions](tables/monthly_header_performance.csv), [known processing coverage](tables/known_processing_coverage.csv), [object-night ledger](tables/annual_object_nights.parquet), [retention coverage](tables/retention_processing_coverage.csv), [candidate coverage](tables/candidate_processing_coverage.csv).
- [Final figures](figures/) — PDF/PNG; seven statistics PDFs are synchronized into `manuscript/figures/annual/`. Figures 1/8 use the supplied R3 display precision; Figure 2 retains the approved workflow. Plot cells/trends are in `figure_data/` and `figure_templates/data/`.
- [Exploratory annual rotation-period work](period_analysis/RESULTS_20261128.md) — independent offline port, 520-target search and 16-target review queue; no confirmed new periods; does not change the manuscript statistics.
- [Handoff](ANNUAL_HANDOFF.md) — cohort changes, missing-product recovery, interpretive points and the remaining review-completion evidence.
- [Review queue](tables/candidate_review_queue.csv.gz) and [contact sheets/index](review_packet/) — 2,286 science links without explicit review-completion evidence. These are not rejected links. Three temporal views of each preserved GIF are included in 48 sheets; all source GIFs were available. A bare default zero is not assumed to prove completed review.

## Reproduce the delivered results locally

Use Python and dependencies in `requirements.txt`. Existing commissioning scientific inputs remain in the adjacent `paper_revision_20261015/` directory pinned by the previous commit. They are reused, not overwritten. The complete compact known tables are here; the fixed ASTORB join contains only the elements needed by these identities, not the large catalog dependency. The DE432s kernel is intentionally untracked: obtain the same JPL kernel used by the earlier delivery, SHA-256 `363f32e14f5255359ac32c4d38080cf28ab55564a5e16696a75f63394b666e9b`, into `inputs/de432s.bsp`.

From this directory:

```sh
python scripts/build_annual_ledger.py
python scripts/build_annual_statistics.py
python scripts/build_annual_candidates.py
python scripts/summarize_annual_coverage.py
python figure_templates/redraw_figures.py --prepare --source .
python scripts/make_figures.py
python scripts/finalize_annual.py
python scripts/sync_annual_prose.py
cd manuscript
tectonic --keep-logs manuscript.tex
# or latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript.tex
```

`exposures.py` and `known_performance.py` preserve the earlier helper implementations and are imported by the annual wrappers; their historical standalone entry points are not the annual entry points. `editorial_source/` stores the supplied prose before result-specific scope updates. `figure_templates/` now contains the annual plotting data, not the old example counts. The locked design and original reference figure hashes are documented in `FIGURE_STYLE_LOCK_20261125.md`.

## Server collection and offline recovery

The source collection scripts read existing FITS headers and selected catalog/stage columns. Each collected science file has path, size, mtime and SHA-256 in `inputs/annual_product_provenance.csv.gz`; none of those source products has a post-November-15 mtime. Exposure membership still uses exposure time, not file mtime. The standard EasyTier SSH connection is documented in the project README.

`frozen_methods/` contains byte-preserved source at `d2f0057`, with hashes in `source_manifest.json`. `scripts/reconstruct_missing.py` uses these only in `/pipeline/xiaoyunao/paper_revision_20261115/`: it recovered the missing pre-subtraction links for 105 already-selected April exposures (48 raw reconstructed links, 44 eligible known-object nights recovered), and evaluated two June exposures with zero primary/auxiliary matches. The full ASTORB reference was read from the existing preserved file whose SHA matches the previous scientific input (`ea22c706…3040`); it was not downloaded or updated. No whole-year production rerun occurred.

The three linkage-stage counts agree for all 151,222 eligible annual object-nights. The final ratio is 112,029/151,222 = 74.08%; no eligible object-night remains unevaluated. The two June exposures have no eligible three-exposure object-night; they do not enter that denominator.
