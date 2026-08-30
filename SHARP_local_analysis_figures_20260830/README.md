# SHARP local analysis and figure bundle — 2026-08-30

This directory is the manuscript-independent handoff requested in `CODEX_LOCAL_ANALYSIS_AND_FIGURES_20260830.md`. It contains only data verification, statistics, scientific plots, plotting data, audit metadata, scripts, and QA products. It does not contain or modify manuscript, LaTeX, production-code, MPC, or JPL submission products.

## Executive results

- Frozen baseline reconciliation: 13/13 requested counts pass, including 534,780 detection-weighted known-object associations, 58,482 unique catalog identities, 4,762 automated single-night linkages, 67 visually retained candidates, 9 later artifacts, and 58 final retained linkages comprising 179 measurements on 34 nights.
- Exact production ASTORB SHA-256: `ea22c7062e7349d9add34fb64835916d66705ad60a717cc7c86a0c47f5aa3040`; catalog orbit epoch: 2026-03-01. The join is one row per known-object identity. See `tables/known_orbit_join_audit.csv`.
- The frozen matched table's distance and solar-elongation columns were entirely null. Figures 04b/04c therefore recompute elongation from Astropy solar positions and distances from the frozen ASTORB two-body elements at each observation epoch. Full derived values are included in Parquet form.
- PHA classification is deliberately reported as unavailable because this frozen ASTORB file does not supply an authoritative PHA flag. NEO is a transparent derived class using `q = a(1-e) < 1.3 AU`.
- Figure 06 uses three real calibrated L1 exposures, a common WCS grid, one moving retained linkage, and the same fixed unsaturated reference star. Registration passes at a maximum 0.445 arcsec reference-star offset against a 0.8 arcsec QA threshold.
- Time-domain Figures 07 and 08 were not fabricated because no explicit calibrated-source manifest was found; see `MISSING_INPUTS.md`.

## Important current-code discrepancy

The current repository is not uniformly 960 m. The scheduler and known-object matcher use 960 m, but `heliolincrr/orbit_confirm_links.py` currently defines `MPC327_ALT_M = 868.221`. This analysis preserved that current orbit-confirmation value for exact reproducibility and recorded it in `audit/code_scope.json`. This is a code/manuscript reconciliation item, not a change made by this bundle.

## Data grains

- Detection-weighted: sky, motion, distance, photometry, and residual plots.
- Unique-object: orbital-element and revisit products, one joined orbit per identity.
- Linkage-level: candidate RMS, arc length, rate, direction, and distance-grid diagnostics.
- Measurement-level: representative trajectories, residual vectors, and Figure 06 astrometry.

## Layout

- `audit/`: code/catalog scope and baseline reconciliation.
- `tables/`: manuscript-ready statistical tables.
- `figure_data/`: complete binned or row-level plotting inputs.
- `figures/`: vector PDF and 300-dpi PNG figures.
- `scripts/`: reproducible generation and validation programs.
- `qa/`: contact sheet and machine-readable figure QA.
- `logs/`: run summaries and validation logs.

## Reproduction

Run `scripts/run_all.sh` from the repository root after providing the exact ASTORB and four calibrated L1 files listed in the script. The general analysis uses `/opt/anaconda3/bin/python3`; solution-family recomputation uses Python 3.10 with `poliastro==0.17.0`. Raw calibrated images and the 407-MB ASTORB catalog are intentionally not included in the ZIP; their identities and checksums are preserved in the audit.

All plotted thresholds remain nominal pipeline diagnostics. In particular, the 1 arcsec line is not a completeness or contamination claim, and the earlier 4.02% association estimate is not reinterpreted here.
