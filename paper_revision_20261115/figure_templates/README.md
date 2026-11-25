# Annual Figures 2, 5 and 6

The approved R2 layouts are retained, using the annual exposure selection and object-night ledger in the adjacent parent directory. This folder now holds annual compact data and figures. Its source method/style reference is commit `1e278897`; the science interval is UTC midpoint [2025-11-15T04:00, 2026-11-15T04:00).

From the parent paper directory, run:

```sh
python figure_templates/redraw_figures.py --prepare --source .
```

Without `--prepare --source .`, the script redraws only from the included compact data. Source file hashes and exact selected/displayed counts are in `summary.json`. Figure 5 density and trend lines share the same plotting selection. Figure 6 uses the single annual evaluated object-night cohort, 0.25-mag by 3-arcsec/hour cells, N<30 ratio masking and original 95% Wilson bins. Its magnitude grid extends to 21 and its one-dimensional y range to 0–100% to include the new faint, low-count bin and full intervals; the layout, font hierarchy, colors and line styles remain unchanged.

Figure 2 is unchanged in scientific content. The manuscript and final figures directory use the exact author-approved R3 PDF. The plotter also creates a scientifically identical workflow when reproducing the full set.

Final PDF/PNG copies are synchronized with `scripts/finalize_annual.py`. The scripts never invoke production code or update a scientific catalog. `../FIGURE_STYLE_LOCK_20261125.md` documents the accepted design. Captions are in `captions.md` and the compiled manuscript; their numerical macro values come from `../annual_results.json`.
