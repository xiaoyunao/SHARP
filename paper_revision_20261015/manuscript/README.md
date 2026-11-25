# SHARP manuscript — 2026-11-24 R2 integrated 2026-11-25

This is the GPT Pro R2 LaTeX manuscript integrated into the paper branch, with the current upstream calibration call-chain correction and regenerated Figures 2, 5 and 6. Scientific inputs remain fixed at `2f5f71d5a074c4cfd712b0e3fbef240c40745e2a`; new dates refer to analysis and manuscript preparation.

Build here with `tectonic --keep-logs manuscript.tex` (verified with Tectonic 0.16.9), or use a normal TeX Live installation with `latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript.tex`. The AASTeX 7.0.1 class and bibliography style are supplied with the source manuscript. Compilation byproducts are ignored; the compiled PDF is a deliberate paper deliverable.

The current PDF has 17 pages, eight main figures and two appendix figures. `sections/06_transient_sources.tex` intentionally contains only the section heading and label. Author metadata and the colleague's transient content will be incorporated when supplied. See [the R2 handoff](../REVISION_R2_HANDOFF.md), [figure reproduction](../figure_revision/README.md) and [calibration method audit](../inputs/calibration_method_note_20261124_r2.md).
