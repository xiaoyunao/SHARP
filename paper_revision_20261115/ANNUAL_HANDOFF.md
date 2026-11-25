# SHARP annual handoff — 2026-11-25

## Final interval and methods

- Exposure midpoint UTC **[2025-11-15T04:00:00Z, 2026-11-15T04:00:00Z)**; Beijing local-noon labels 20251115–20261114 (365 complete possible nights). The final evening of November 15 is outside this interval.
- This replaces the earlier October cutoff for observations. Processing methods remain the paper's pre-November-15 definitions. Raw images, calibrated images/catalogs and production state are unchanged. Current Gaia calibration implementation is inherited from R2, with propagated WCS coordinates and original-coordinate white-light photometric matching; no earlier-backup inquiry or recalibration was repeated.
- Initial input commit: `1e278897745dc9d879cb2b0e2da0a6687cb843f1`. Frozen algorithm: `d2f0057bd4cdd4c7a3c6c2431f6449a72c327284`; executable files are byte-preserved in `frozen_methods/`, with SHA-256 manifest. Lowell ASTORB SHA-256 is `ea22c7062e7349d9add34fb64835916d66705ad60a717cc7c86a0c47f5aa3040`; the existing preserved copy was found and its selected orbit elements exported. The current later catalog was not used.
- The standard sample is 30-second white-light strict MP observations. Other modes/filters are tabulated separately. In particular November 11–13 is not silently treated as white-light standard survey data. Directory labels only locate files; membership and plotting activity use physical nights. Earlier timing/WCS exclusions and proven duplicate rules are retained.

## Annual results and differences from the earlier sample

| Quantity | First-year result | Previous commissioning result |
|---|---:|---:|
| Science exposures / physical nights | 53,398 / 168 | 39,590 / 129 |
| Open shutter | 445.02 h | 329.95 h |
| Once / same-night ≥3 geometric area | 11,025.29 / 10,566.82 deg² | 10,928.09 / 10,470.59 deg² |
| Adjacent cadence median | 31.98 min | 37.65 min |
| ≥3-visit pointing-group fraction | 71.09% | 67.43% |
| Model-depth median / IQR | 20.35 / 19.88–20.67 mag | 20.25 / 19.76–20.58 mag |
| Primary associations / identities | 797,887 / 80,591 | 531,612 / 58,281 |
| Main-belt-region fraction | 95.84% | 95.73% |
| Radial O−C median (1-arcsec selection) | 0.25 arcsec | 0.26 arcsec |
| Conditional linkage ratio | 112,029 / 151,222 = 74.08% | Qualified historical-policy cohort 75,467 / 101,784 = 74.14% |
| Reviewed retained links / measurements | 114 / 354 | 58 / 179 |

The annual single ratio is used throughout Figure 6, abstract, text, tables and conclusions. Earlier 74.16/74.14/74.21 definitions remain in the old directory only. Counts, precision and source fields are generated from `annual_results.json` through `scripts/finalize_annual.py`.

## Product coverage and successful incremental recovery

- Existing archived stage products were reused, including actual products on historical quarantine nights. Quarantine labels were not used to delete a whole night or convert missing products into failures.
- A 105-exposure archive group in directory 20260414 belongs physically to night 20260413. Its missing stages were reconstructed using the original Gaia-mask, tracklet, linear-link and orbit-confirmation code in the paper workspace. This produced 48 pre-subtraction links; 44 eligible known-object nights meet the strict annual membership rule. `inputs/reconstruction_receipt.json` records product hashes.
- Two 20260605 calibrated exposures had no archived known output. The frozen matcher with fixed ASTORB/DE432s was run only on those two inputs. There were 195 and 109 prediction opportunities and **zero** primary or 1.5-arcsec matches. This is a measured empty result, not a missing result. There cannot be an eligible ≥3-exposure known-object night in those two images.
- **All 151,222 eligible object-nights are evaluated.** Motion-stage coverage is 53,396 science exposures on 167 physical nights; two June exposures form the no-eligible exception. Eligible objects occur on 165 nights. `strict_linked`, `fit_ok_survived`, `is_good_survived` all count 112,029; 19 ambiguous detection keys are excluded from identity mapping.
- The annual magnitude panels use 151,165 eligible / 112,026 linked cases; 57 invalid/placeholder magnitudes include three linked cases. Rate and overall counts retain these 57. The 0.25-mag grid extends to 21 mag to include three real object-night medians just fainter than 20; N<30 masking and 95% Wilson definitions are unchanged. The 1D retention y axes extend to 0–100% to show the new low-count confidence intervals in full.
- Monthly depth/seeing/zero-point distributions and missingness are in `tables/monthly_header_performance.csv`. The model-depth statistic uses 52,452 valid estimates; 946 `MLIM=-99` placeholders are unavailable for this metric and remain in the exposure/astrometry cohort. Depth histograms include all valid values through 22.78 mag. Figure 5 uses 749,572 plotted positive-error/flux records; three valid-positive records fall outside the plotting ranges. Density and trend lines share the plotted selection. The radial panel shows 796,336/797,887 associations, excluding 1,551 magnitude-range outliers; nightly offsets use all 797,887 across 167 nights. Six stellar sequences remain the dated February–May example, not an annual stellar error estimate.

## Candidate review and the actual remaining work

- Archived automatic output contains 8,449 links; after annual exposure-member selection, 8,032 remain on 152 physical nights. These numbers have different selections and are both preserved.
- Existing explicit review/completion evidence covers 5,746 science links on 144 physical nights. The retained reviewed subset contains **114 links / 354 measurements / 57 nights**: all previous 58/179, plus 56 additional positive image-review records passing the same exposure selection. The explicit 2026-09-04 false-positive decision removes `00001oI` and `00001oj` despite positive values still present in the earlier review CSV. Manual-action provenance is supplied.
- **2,286 links on 25 nights lack explicit review-completion evidence.** A default saved zero is not interpreted as proof of human rejection. This is a completion-evidence queue; some may already have been inspected without a retained completion marker. Before requesting fresh inspection, compare the queue with any author-held review record. The exact CSV is `tables/candidate_review_queue.csv.gz`; `review_packet/review_preview_index.csv` points to 48 contact sheets and the original GIF paths. All 2,286 source GIFs were available; each card contains three temporal views.
- Six retained links keep the existing numerical C/2025 Y1 association. The remaining 108 lack established identities in the available records; no new identity research was undertaken or invented. Additional identities, if independently documented later, can be merged by night and tracklet ID.
- There are 28 distant-region positional identities in the annual known sample. The old 17/119 local audit is reused as historical support in the prior directory, not presented as annual 79/9/31 or 7/7/3 classifications. The author allowed reducing this optional discussion to overall composition; the three unused annual classification macros are marked not evaluated and make no prose claim.

## Manuscript and points for GPT Pro

The supplied R3 author list, science-led Introduction, merged Discussion and Conclusions, and symbolic calibration coefficients are retained. Transient Sources is title-only; Data and Software Availability remains deleted. Seven annual statistics figures and all used number fields are filled, and `annualresultsready` is true. Publication figures retain approved layouts; Figures 1/8 preserve the supplied two-decimal labels. Full LaTeX and compiled PDF are committed directly.

Actual prose changes are limited to: standard-mode exclusion of the November 11–13 observations; explicit annual plot selections; the reviewed-subset/unfinished-review scope; and removal of unsupported annual distant-identity classifications. GPT Pro can now refine the interpretation of these measured trends:

1. Exposure/identity counts grow strongly, while area changes modestly: the extension principally adds repeat sampling. Revisit median shortens from 37.65 to 31.98 min and the ≥3-visit group fraction rises to 71.09%.
2. Model-depth median is 20.35 mag, with a brighter/fainter condition distribution shown in Figure 5; this remains a noise-model quantity, not completeness.
3. Aggregate conditional retention remains near 74%, with the same low-rate deficit. Sparse faint bins have wide intervals; do not interpret the final three-object bin as a new sharp feature.
4. The known sample remains 95.84% main-belt-region identities. Its median per-object span changes to 5.00 days; the maximum number of nights rises to 33. These are object-level detected sampling statistics, distinct from field cadence.
5. The 114 retained links are a reviewed subset of the first-year search, not a fully reviewed year or 114 independent new objects.

The remaining author/colleague inputs are unchanged: response-curve definition/data, transient section, correspondence/email/funding details and final in-preparation bibliography. The only new practical item is resolving or inspecting the explicitly listed review-completion queue. No additional generic method audit, injection test, orbit refinement, recalibration or production development is requested.
