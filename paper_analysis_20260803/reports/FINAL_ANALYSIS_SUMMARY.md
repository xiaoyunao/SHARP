# Final analysis summary

## Delivery state

All locally executable P0/P1 analysis for the frozen 2026-08-03 paper snapshot
is complete. The manuscript source was not edited.

- Science interval: 2025-11-15 through 2026-07-15, inclusive.
- Tables: 5 CSV and 5 TeX products under `../tables/`.
- Figures: 12 vector PDF and 12 300 dpi PNG products under `../figures/`.
- Reproducibility: executed notebook, 555-artifact SHA256 manifest, and complete
  machine-readable snapshot.
- Validation: 63 PASS, 0 FAIL, 1 SKIP. The sole SKIP is author sign-off on the
  provisional night-quality mask.

## Paper-facing results

| Topic | Frozen result | Interpretation boundary |
|---|---:|---|
| Strict raw science frames | 41,074 over 134 nights | 342.283 h open shutter; 1,430 acquired fields |
| Strict L2 catalogs | 40,399 over 131 nights | Missing-L2 nights remain explicit, not zero-filled |
| Nominal known predictions | 13,310,546 unique keys | V<22.5 plus nominal WCS rectangle; not true detectable completeness |
| Known matches within 1 arcsec | 534,780 (4.017716%) | Nearest-neighbor, threshold-truncated association |
| Known radial residual | median 0.263406 arcsec | Signed medians: -0.017433 arcsec in RA cos Dec, +0.022456 arcsec in Dec |
| Formal unknown catalog | 4,762 linkages | 14,299 linkage-detection memberships; 14,159 global unique detection keys |
| Human states | 68 -> 67 -> 58 | marked real -> submission selected -> post-hoc retained |
| Cross-night screen | 37 candidate components | Linear-motion screen only; not an independent-object count |
| JPL second-pass association | 6 linkage rows with C/2025 Y1 | Numerical candidate association, not MPC attribution or designation |
| Plan compliance | 7,455 / 7,833 = 95.1743% | Plan-active acquired cohort only; not observatory efficiency |
| Operations normal cohort | 2 nights | CPU/RAM and most human/MPC intervals unavailable |

## Observer-location finding

The current production-aligned repository does not use one observer location:

- scheduler and known-object matching: 117.575 deg E, 40.393 deg N, 960 m;
- short-arc orbit confirmation: 117.575 deg E, 40.394239 deg N, 868.221 m;
- L1/L2 headers: contaminated GMG template values and therefore unusable as the
  site authority.

A complete 128-night alternate-location rerun changed neither `fit_ok` nor
`is_good` for any of 87,850 link keys, the formal 4,762, or the retained 58.
This is an observer-location sensitivity result, not proof that the production
code already uses 960 m. Some three-observation short arcs still show very large
element changes; orbital elements must not be used for population inference.

## Manuscript corrections with highest priority

1. Separate `fit_ok` from `is_good`; the formal unknown selector is
   `fit_ok && all_non_asteroid`, while the final 4,762 happen also to be
   `is_good`.
2. Replace “14,299 unique detections” with either 14,299 memberships or 14,159
   unique detection keys.
3. Do not describe historical near-Sun scheduling as minimum-solar-separation
   prioritization: the solar sort was overwritten by RA sorting.
4. Replace “prevents repeat selection” with a down-weighting/reweighting
   description, and disclose that production history stopped updating after
   2026-03-30.
5. Describe static-source filtering as reference-anchored, note the RA=0 link
   wrap defect, and distinguish unknown ADES start times from midpoint times.
6. Use linkage, component, JPL association, and MPC identity/designation as
   separate statistical grains. No discovery count is currently supported.

The evidence and rerun boundaries for all conflicts are documented in
`MANUSCRIPT_CODE_CONFLICTS.md`.

## Inputs still required from the author or collaborators

- surveyed MPC 327 longitude, latitude, height type, and datum;
- signed night-quality policy and anomaly-night inclusion decisions;
- upstream reduction/depth/effective-mask/header documentation;
- authoritative MPC row records for the 58 retained plus 9 rejected links;
- policy for the six C/2025 Y1 candidates and the 15 s unknown-time correction;
- an approved injection grid and upstream image-processing route;
- independent CPU/RAM, human-review, weather, fault, and override records.

The minimum response template is in `AUTHOR_INPUTS_REQUIRED.md`.
