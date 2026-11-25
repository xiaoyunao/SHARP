# Annual captions (LaTeX; values from annual_numbers.tex)

```latex
\begin{figure}[tbp]
\centering
\includegraphics[width=\linewidth]{fig02_workflow.pdf}
\caption{The observing and identification sequence. Repeated white-light
imaging produces calibrated images and source catalogs. The upper branch
associates sources with ALEPH predictions to form the primary
$1\asec$ known-object sample. The moving-source branch constructs
tracklets and links after stationary-source rejection, tests short-arc
consistency, and removes links identified by the auxiliary $1.5\asec$
known-object associations. The dashed connection supplies those
identities after motion linking. Image inspection selects the
single-night candidates. L1 and L2 denote calibrated images and
source catalogs, respectively; ADES is the Astrometric Data Exchange
Standard.}
\label{fig:workflow}
\end{figure}
```

```latex
\begin{figure}[tbp]
\centering
\annualplot{fig05_measurement_quality.pdf}{\linewidth}
\caption{Photometric and astrometric measurements during the first
year. (a) Catalog $5\sigma$ model depth in 0.2-mag bins; the orange
line and band mark the median and interquartile range of 52,452 valid estimates. The 946 unavailable header values are excluded only from this metric.
(b) Aperture-magnitude uncertainty for unflagged primary associations
with positive magnitude uncertainty, flux, and flux uncertainty.
The displayed sample contains
749,572 measurements with usable magnitudes and uncertainties
within the plotted ranges; the same selection defines the orange curves.
Cells are 0.1~mag by 0.04~dex. Orange curves show the median and
16th--84th percentiles in 0.5-mag bins with at least 30 measurements.
Blue points show the stellar repeatability from the six selected
sequences, selected between 2026 February and May, using bins with at least 50 stellar pairs.
(c) Radial O$-$C in 0.1-mag by $0.02\asec$ cells, with the same
trend conventions. The color bar gives the number of measurements
per cell. These associations lie within $1\asec$ of the predictions.
(d) Nightly medians of $\Delta\alpha\cos\delta$ and $\Delta\delta$;
the gray bars show the number of associations per night.
The radial panel contains 796,336 associations;
1,551 magnitude values lie outside its displayed range.
The nightly panel uses all \yrNAssoc\ associations.}
\label{fig:quality}
\end{figure}
```

```latex
\begin{figure}[tbp]
\centering
\annualplot{fig06_conditional_retention.pdf}{\linewidth}
\caption{Conditional linkage retention for first-year known-object
nights. Eligibility requires detections in at least three distinct
science exposures within $1.5\asec$ of one known identity, with
predicted first-to-last angular rate $3$--$63$~\rateunit.
(a) Eligible counts in 0.25-mag by 3~\rateunit\ cells, with a
logarithmic color scale. (b) Linked fractions in the same cells;
cells containing fewer than 30 eligible object--nights are gray.
(c,d) One-dimensional linked fractions with 95\% Wilson intervals
and horizontal bars indicating bin widths. The lower histograms
show eligible counts per bin. The dashed line marks the overall
fraction, \yrRetention\%. Objects with unusable magnitudes remain
in the overall and angular-rate samples and are omitted only from
panels requiring a magnitude. The magnitude-dependent panels contain
151,165 eligible object--nights and
112,026 linked cases. The
57 excluded magnitudes include
3 linked cases. Both maps use unsmoothed counts.}
\label{fig:retention}
\end{figure}
```
