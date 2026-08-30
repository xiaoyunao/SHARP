# Handoff to GPT Pro

Please use this ZIP only as the verified local analysis/figure supplement for the SHARP PASP structural revision. It does not replace the manuscript source and contains no manuscript edits.

Start with `README.md`, then use:

- `audit/baseline_reconciliation.csv` for the reconciled headline counts and their grains;
- `audit/code_scope.json` for the exact current-code scope and the confirmed altitude discrepancy;
- `audit/catalog_scope_c2025y1.json` and `tables/known_orbit_join_audit.csv` for the orbital-catalog scope and join quality;
- `tables/candidate_linkage_diagnostics.csv` for the 58 retained linkage-level diagnostics;
- `figures/*.pdf` as publication vector masters and `figures/*.png` as previews;
- `GOTTA_SOURCE_MAP.md` for the exact GOTTA-source mapping used in the complete redraw;
- each matching file in `figure_data/` for plotted values and reproducibility;
- `qa/figure_qa.json` and `qa/figure_contact_sheet.png` for QA;
- `MISSING_INPUTS.md` for why no time-domain Figures 07/08 were generated.

Important factual constraint: current code is not uniformly at 960 m. Scheduler/known matching use 960 m, while the current short-arc orbit-confirmation source still defines 868.221 m. This bundle preserves and audits that current-code discrepancy instead of silently normalizing it.

Other constraints to retain: known-sky products are detection-weighted; orbit/revisit products are unique-object weighted; candidate diagnostics are linkage-level; representative trajectories are measurement-level. The 1 arcsec line is a nominal pipeline threshold, not a completeness/contamination result. PHA is unavailable from the frozen ASTORB source; NEO is transparently derived as `q < 1.3 AU`.

Figure-specific editorial constraints: all photometric labels are calibrated aperture magnitude (`g_aper`), never instrumental magnitude; Figure 05c has been removed; Figure 06 is exactly the three author-selected GIF frames in a 1x3 row with no added material.
