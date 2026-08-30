# Exact GOTTA plotting-source map

This redraw was made from the plotting source of `/Users/yunaoxiao/Desktop/gotta_asteroid_1`, not from visual imitation. Checksum-preserved source copies are bundled in `scripts/gotta_reference_source/`; their SHA-256 values are recorded in `audit/gotta_style_provenance.json`.

| SHARP output | Exact source anchor | Retained plotting grammar |
|---|---|---|
| Figure 01 | `make_fig1_review.py` | `ZScaleInterval(contrast=0.25, krej=2.5)`, 61-pixel median-filter background subtraction, shared median/MAD cutout limits, gapped crosshairs, one full-field plus 2x2 cutouts, no N/E labels |
| Figure 02 | `paper_analysis_20260803/scripts/make_fig01_architecture.py`, the earlier full SHARP architecture requested by the author; typography/color grammar follows GOTTA `redraw_text_only.py` | Full production topology, code-path validation, operational feedback loop, rounded workflow nodes and routed arrows |
| Figure 03 | `make_ecliptic_healpix.py` | NSIDE=64 Mollweide HEALPix density map, `rainbow` colormap, `LogNorm(1,100)`, orange celestial equator, embedded horizontal colorbar |
| Figure 04 | `make_orbit_revision.py` | Exact eight-class order, marker shapes, colors, alpha/size rules, log semimajor-axis and inclination axes, and `q=1.3 AU` curve |
| Figures 04b--04e | `generate_paper_products.py` | Times New Roman 30-point base style, 1.4-point axes, stepfilled histograms, `viridis` local-density scatter, appended GOTTA colorbars, red binned medians and 16--84% bands |
| Figure 04f | temporal-sampling panel in `generate_paper_products.py` | Ordinary circular scatter markers, `viridis`/`LogNorm` count encoding, appended GOTTA colorbar, stepfilled log histogram |
| Figures 05 and 05b | `generate_paper_products.py` style helpers | Ordinary circular linkage markers, GOTTA `viridis` colorbar, white/orange track markers, blue orbit track, red residual vectors, unclipped labels and axes |
| Figure 06 | Author-specified `20260210_00000P2_link0035.gif` | The three supplied frames are placed exactly in one row; no additional labels, arrows, directions, registration, or other content |

Photometry is labeled as calibrated aperture magnitude, `$g_{\rm aper}$`. The rejected term “instrumental magnitude” does not occur in the regenerated plots or linkage table. Figure 05c is intentionally absent.
