# Replacement captions

These captions describe only the figures and frozen data in this directory.
Scientific source commit: `2f5f71d5a074c4cfd712b0e3fbef240c40745e2a`.

## Figure 2

Processing sequence of the frozen SHARP survey workflow. Repeated white-light imaging is reduced to calibrated L1 images and L2 source catalogs. The upper branch uses ALEPH predictions for the primary 1-arcsec known-object associations and their measurement products. The independent moving-source branch rejects Gaia and stationary sources, constructs tracklets and linkages, applies short-arc consistency checks, and then removes linkages identified through auxiliary 1.5-arcsec known-object associations. The dashed connection supplies known-object identities at this classification step, after moving-source linking. Image review yields the retained single-night candidate sample. The diagram distinguishes the primary known-object measurement sample from the auxiliary associations used for the known-object veto.

## Figure 5

Measurement quality in the validated exposure sample. (a) Distribution of the catalog 5-sigma model depth for 39,590 exposures, in 0.2-mag bins. The orange line marks the median (20.25 mag), and the shaded interval is the interquartile range (19.76--20.58 mag). (b) Catalog aperture error versus calibrated aperture magnitude for 500,092 primary known-object associations with zero extraction flags and positive magnitude error, flux, and flux error. The background shows measured counts in cells of 0.1 mag by 0.04 dex in aperture error. The solid orange curve and dashed curves in (b,c) show the median and the 16th and 84th percentiles in the frozen 0.5-mag bins with at least 30 measurements. Blue points in (b) show the independent per-epoch stellar repeatability estimate, 0.7413 IQR(Delta m)/sqrt(2), from six selected field-night triplets, retaining bins with at least 50 stellar pairs. (c) Great-circle radial O-C separation versus magnitude, with counts in 0.1-mag by 0.02-arcsec cells. This panel displays 531,416 of the 531,612 primary associations: 195 entries with magnitude 99 and one entry at 23.18 mag lie outside the displayed magnitude range. The residual distribution is conditional on the primary 1-arcsec association selection. The shared logarithmic color bar in (b,c) gives measurements per cell. (d) Nightly medians of the signed Delta-alpha cos(delta) and Delta-delta residuals for all 531,612 primary associations over 128 local-noon observing nights; the lower gray bars give the association count in each night. The historical g_aper symbol denotes the Gaia-referenced white-light aperture scale defined in the text.

## Figure 6

Conditional linkage retention in the historical object-night sample. An eligible object-night contains at least three distinct detections associated with one known object within 1.5 arcsec and a first-to-last predicted rate of 3--63 arcsec per hour. Retention requires at least one frozen linkage containing at least three members, all mapping unambiguously to that identity. (a) Eligible object-night counts in cells of 0.25 mag by 3 arcsec per hour; the color scale is logarithmic. (b) The retained fraction in the same cells, with cells containing fewer than 30 eligible object-nights masked in gray. Neither map uses smoothing. (c,d) The frozen one-dimensional retained fractions versus median measured aperture magnitude and predicted rate. Vertical bars are 95% Wilson intervals, horizontal bars show bin widths, and the lower gray histograms give eligible counts per bin rather than probability densities. The dashed reference is the full historical fraction, 75,899/102,347 = 74.16%. The magnitude-dependent panels (a--c) contain 102,315 eligible and 75,897 retained object-nights; 32 object-nights with unusable catalog magnitudes, including two retained cases, are omitted from these panels but remain in the overall and rate-dependent denominators. This statistic is conditional on archived detections and measures linkage retention, not image-level or L2 detection completeness.

## Exact quantities for manuscript synchronization

| Quantity | Value |
|---|---:|
| Fig. 5(a) exposures | 39,590 |
| Fig. 5(b) selected and displayed associations | 500,092 |
| Fig. 5(c) selected astrometric associations | 531,612 |
| Fig. 5(c) displayed associations, 9.5 <= g_aper < 22 | 531,416 |
| Fig. 5(c) records outside displayed magnitude interval | 196 |
| Fig. 5(d) associations / local-noon nights | 531,612 / 128 |
| Fig. 6 overall / rate denominator and numerator | 102,347 / 75,899 |
| Fig. 6 magnitude denominator and numerator | 102,315 / 75,897 |
| Fig. 6 magnitude exclusions / retained among exclusions | 32 / 2 |
| Fig. 6(b) minimum eligible count per displayed cell | 30 |
| Fig. 6(b) displayed cells / represented object-nights | 290 / 99,353 |

The extra 32 magnitude exclusions and the 196 Figure 5 display exclusions are
provided as explicit record lists in `data/`. They are figure-specific masks;
no historical global denominator, link label, or primary-association selection
has been changed.
