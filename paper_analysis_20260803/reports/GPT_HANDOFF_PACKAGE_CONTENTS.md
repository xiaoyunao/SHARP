# GPT handoff package contents

The upload bundle is intentionally compact. It contains the complete original
LaTeX draft plus the publication-facing products needed for the next manuscript
revision, but excludes multi-gigabyte row-level extracts and live server logs.

## Top level

- `START_HERE.md`: the prompt and non-negotiable evidence boundaries for GPT.
- `manuscript/`: the original PASP draft, expanded from
  `SHARP_PASP_draft.zip`.
- `analysis_support/`: frozen analysis products for manuscript revision.
- `MANIFEST_SHA256.txt`: SHA-256 for every other file in the bundle.

## Manuscript directory

The original `.tex`, `.bib`, `.md`, build files, preview PDF, and planning files
are preserved. The PASP entry point is `manuscript/manuscript.tex`, while
`manuscript/manuscript_preview.tex` is the preview entry point.
`manuscript/figures/` additionally contains all twelve final
figures as both vector PDF and 300 dpi PNG. The exact PDF names match the
original `figure_manifest.csv`.

## Analysis support

- `reports/`: final summary, all 19 manuscript/code conflicts, author-input
  checklist, P0/P1 status, validation report, and the GPT handoff prompt.
- `tables/`: Tables 1--5 as CSV and TeX, plus hashes and the table summary.
- `figure_data/`: compact, paper-facing CSV/JSON inputs used by the figures;
  the 33 MB exposure-level Fig. 11 denominator is omitted because the plotted
  result and retained-link table are already included.
- `evidence_summaries/`: compact nightly and aggregate results for known,
  unknown, scheduler, operations, site sensitivity, JPL, review, coverage,
  random-shift, and time-correction analyses.
- `qa/`: 12-figure contact sheet, figure QA JSON, and independent snapshot
  validation JSON.
- `notebook/`: executed lightweight audit notebook.

## Deliberately excluded

- multi-gigabyte per-detection CSV/Parquet/FITS products;
- raw images, L2 catalogs, production logs, and server temporary files;
- figure-generation source code, because the figures are frozen and the next
  task is manuscript revision rather than re-analysis;
- authoritative MPC information that has not been supplied.

The full 5.7 GB frozen local analysis remains at
`/Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/` if a later
numerical audit needs the omitted row-level evidence.
