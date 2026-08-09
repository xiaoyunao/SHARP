# Manuscript–Code Conflict Audit

## Scope and evidence boundary

This report is a read-only audit of the draft delivered in
`SHARP_PASP_draft.zip` against the current production-aligned repository in
`/Users/yunaoxiao/Desktop/smt_asteroid` and selected current server products.
It does **not** edit the manuscript.

- Manuscript snapshot date: 2026-08-03.
- Manuscript production commit: `d2f0057`.
- Repository commit inspected: `cfee217` at audit start.
- `git diff d2f0057..cfee217 -- survey known_asteroid heliolincrr` was empty;
  therefore the algorithm findings below apply to both the paper-frozen code
  and the inspected current repository state.
- Server evidence was collected read-only from the production paths defined in
  the project instructions.

Rerun labels used below:

- **Text only**: correct the method wording; existing numerical products do not
  need regeneration.
- **Analysis rebuild**: rebuild paper statistics/tables/figures from existing
  frozen products; no upstream production stages need rerunning.
- **Partial pipeline rerun**: rerun the named stage and all downstream derived
  products.
- **Full unknown rerun**: rerun Gaia/static filtering through review/MPC
  reconciliation because candidate membership can change.
- **User decision required**: do not change production behavior or submit/correct
  MPC data until the author explicitly selects the policy.

## Executive findings

1. The L1/L2 headers sampled on the server contain a wrong GMG observatory
   template and cannot be used to choose the Xinglong site coordinates. A
   completed 128-night `868.221 m / 40.394239 deg` to
   `960 m / 40.393 deg` counterfactual caused zero `fit_ok` or `is_good`
   classification flips, but a small number of short-arc orbital-element
   solutions remained highly unstable and cannot support population inference.
2. The unknown catalog is selected on `fit_ok`, whereas the draft describes
   `is_good`-style thresholds. The 4,762 published rows happen to satisfy both,
   but that equivalence is not guaranteed by the code.
3. The reported 14,299 value is a linkage–detection membership count, not a
   global unique-detection count. The latter is 14,159 for the frozen 124
   catalogs.
4. The near-Sun implementation discards its solar-separation ordering before
   choosing fields.
5. Unknown ADES times are 15 s earlier than the exposure-midpoint convention
   used by the known and orbit branches.
6. The static-source filter and unknown RA-wrap behavior can change candidate
   membership if corrected and therefore define full-rerun boundaries.
7. The 13,311,871 `_all_asteroids` rows are not a unique prediction
   denominator: 1,325 repeated `(night, exposure, object)` rows must be removed,
   leaving 13,310,546 nominal predictions.
8. The 58 retained rows form 37 threshold-dependent linear-motion candidate
   components, not 58 independent sources; six rows in two components are
   independently re-identified by two frozen JPL API passes as C/2025 Y1
   (ATLAS), while their MPC ingest/designation states remain unresolved.
9. Human-review accounting has three distinct gates: 68 rows were marked real
   in the review CSVs, one row (`20260507/00001et`) was withdrawn before
   submission, 67 entered the frozen submission/audit set, and 58 survived the
   post-hoc audit. Calling both 68 and 67 “first review retained” hides a real
   decision transition.

## Conflict 1 — Site coordinates and contaminated L1/L2 metadata

**Severity:** P0 / High

**Decision class:** Author/observatory site confirmation remains required. The
full frozen observer-location sensitivity run is complete; the tested change
requires text clarification but no classification-driven downstream rerun.

### Manuscript evidence

- `sections/02_facility_data.tex:5-8`, especially:
  “The scheduler and known-object matcher use a site altitude of 960 m. A legacy
  constant of 868.221 m remains in the short-arc orbit-confirmation module. The
  difference is negligible...”
- `sections/02_facility_data.tex:15-16` lists
  `117.575 deg E, 40.393 deg N, 960 m`.
- `sections/02_facility_data.tex:27` repeats the legacy-altitude note.
- `sections/appendices.tex:73` requires one surveyed observatory position.

### Code evidence

- Scheduler: `survey/config.py:10-12` uses
  `40.393, 117.575, 960.0`.
- Known matcher defaults: `known_asteroid/match_single_night.py:226-228` and
  observer construction at `:265-268` use the same values.
- Production known wrapper: `known_asteroid/slurm_match_one_file.sh:41-43`
  passes the same values.
- Single-night orbit code: `heliolincrr/orbit_confirm_links.py:58-60,110-114`
  uses `117.5750, 40.394239, 868.221`.
- Experimental 15-night code:
  `heliolincrr/run_rr_from_tracklets_15.py:42-65` uses the same legacy orbit
  values.

### Server/header evidence

Read-only checks of the following L2 files and their corresponding L1 files:

- `/processed1/20260111/L2/OBJ_MP_0028_0205_cat.fits.gz`
- `/processed1/20260529/L2/OBJ_MP_0993_0049_cat.fits.gz`
- `/processed1/20260715/L2/OBJ_MP_1545_0038_cat.fits.gz`

all returned:

```text
OBSERVAT = GMG
OBSLON   = 100.0313
OBSLAT   = 26.6974
OBSALT   = 3227
TELESCOP = SMT
```

Those values are incompatible with Xinglong/MPC 327 and are evidence of an
upstream metadata-template problem. They must not be used to adjudicate 960 m
versus 868.221 m.

### Completed 128-night observer-location sensitivity test

- `scripts/run_orbit_site_sensitivity.py:69-77` imports the production orbit
  module read-only and overrides only its observer globals; `:81-143` reruns
  each frozen candidate night into a non-production output root.
- `snapshot/orbit_site960/run_summary.json:2-14` records 128 candidate nights,
  128 completed, zero errors, and the alternate observer
  `117.575 deg E, 40.393 deg N, 960 m`. One night was already complete and was
  skipped safely; the other 127 ran successfully. Five successful nights had
  zero orbit-link rows, so the link-level comparison contains 123 non-empty
  nights without losing any link keys.
- `scripts/compare_orbit_site_sensitivity.py:166-188` performs an outer join on
  `(night, linkage_id)` and computes classification flips; `:203-228` reports
  the all-link, formal-unknown, and retained-58 cohorts separately.
- `snapshot/orbit_site_comparison/orbit_site_sensitivity_summary.json:2-18`
  shows that all 87,850 baseline keys are present in the alternate run, with no
  missing or new rows. `:19-42` preserves the complete-run metadata and both
  observer definitions. The baseline is the historical production result at
  `117.575 deg E, 40.394239 deg N, 868.221 m`; the alternate changes both height
  and latitude to the scheduler/known values. It is therefore an observer-site
  comparison, not an isolated height derivative.

Classification closure is exact:

| Cohort | Rows | `fit_ok` baseline -> 960 m | `fit_ok` flips | `is_good` baseline -> 960 m | `is_good` flips |
|---|---:|---:|---:|---:|---:|
| All orbit links | 87,850 | 87,462 -> 87,462 | 0 | 87,351 -> 87,351 | 0 |
| Formal unknown catalog | 4,762 | 4,762 -> 4,762 | 0 | 4,762 -> 4,762 | 0 |
| Post-audit retained set | 58 | 58 -> 58 | 0 | 58 -> 58 | 0 |

The all-link result is frozen at
`orbit_site_sensitivity_summary.json:2-12`; the 4,762 and 58 subset results are
at `:49-71`.

### Quantified impact

- Height discrepancy: 91.779 m.
- Latitude discrepancy between scheduler/known and orbit code: 0.001239 deg,
  about 138 m north–south.
- Scheduler and known matching explicitly construct their sites from code/CLI
  constants, so the wrong L1/L2 header coordinates do not alter their current
  statistics.
- For the tested observer change, orbit residual differences do not cross any
  `fit_ok` or `is_good` boundary, and neither the formal 4,762-link catalog nor
  the retained 58-link set changes membership. Existing headline counts and
  classifications therefore do **not** need rerunning merely to substitute the
  tested 960 m observer in this frozen comparison.
- “No classification flips” does **not** make the fitted orbital elements
  scientifically robust. Across all fitted links the 95th-percentile absolute
  changes are small, but extreme absolute changes reach 4,427.289 au in
  semimajor axis, 3.509 in eccentricity, and 42.199 deg in inclination
  (`orbit_site_sensitivity_summary.json:75-121`). Angular-element raw
  differences also approach 180--360 deg (`:83-89,147-153,171-177`), partly
  exposing branch/wrap degeneracies.
- This coexistence is expected from the code: the single-night `is_good` gate
  checks residuals, observation/tracklet counts, and linear-motion RMS, but does
  not bound `a`, `e`, or the angular elements
  (`heliolincrr/orbit_confirm_links.py:920-929`). The output has no formal
  `degenerate` flag; “underconstrained/degenerate” here is the scientific
  interpretation of extreme short-arc element sensitivity, not another hidden
  classification cut.
- The instability reaches the published cohorts even though classifications do
  not flip. Formal-unknown examples include `20251119/270`, whose fitted
  semimajor axis changes from 63.552 au to -2.740 au
  (`orbit_site_sensitivity_by_link.csv:1622`), and `20260410/971`, whose
  inclination changes by -14.160 deg (`:77087`). Within the retained 58,
  `20260601/82` changes by 0.311 au in semimajor axis (`:85957`), while
  `20260420/807` changes by 1.586 deg in inclination and 152.359 deg in the raw
  ascending-node value (`:78242`). These are short-arc solution instabilities,
  not evidence for physical population structure.
- The manuscript may call the observer change negligible **for these frozen
  classification counts only**. It must not call the orbital-element impact
  negligible, and it must not use the fitted `a/e/i/Omega/omega/nu`
  distributions for asteroid-population interpretation.

### Required handling and rerun boundary

1. Obtain the surveyed MPC 327 longitude, latitude, height type, and geodetic
   datum from the observatory owner. The completed alternate run does not make
   the scheduler constants an authoritative survey measurement.
2. Preserve the historical production observer (`868.221 m`) in the
   as-executed method, and report the completed 960 m sensitivity result
   separately. Do not silently rewrite the historical code path.
3. Replace the unqualified “difference is negligible” wording with the bounded
   result: zero classification/membership flips for all 87,850 links, the formal
   4,762, and the retained 58, alongside the orbital-element-instability
   warning.
4. Prohibit population interpretation of these single-night short-arc orbital
   elements. Sky-plane motion and residual distributions may still be reported
   with their stated selection limitations.
5. Correct the upstream GMG header template independently of the asteroid code.
6. No orbit/unknown production rerun is required for the already tested
   observer pair. If the authoritative surveyed site differs from both tested
   definitions, repeat the frozen sensitivity comparison; trigger downstream
   regeneration only if classifications, membership, or paper-consumed
   numerical diagnostics change. Never resubmit to the MPC automatically.

## Conflict 2 — `fit_ok` is the catalog gate, not `is_good`

**Severity:** P0 / High

**Decision class:** Text correction plus analysis rebuild; full unknown rerun if
the production gate changes.

### Manuscript evidence

- `sections/05_unknown_pipeline.tex:51` says a single-night solution is accepted
  at RMS <= 5 arcsec, maximum residual <= 8 arcsec, a 5 arcsec clip, and at
  least 80% retained.
- `sections/05_unknown_pipeline.tex:55,59` describes only accepted orbit fits as
  entering the unknown catalog.
- `sections/05_unknown_pipeline.tex:85-88` presents those values as the frozen
  single-night acceptance criteria.
- `sections/07_results_unknown.tex:5-7` describes the 4,762 candidates as the
  output after that orbit-consistency filter.

### Code evidence

- A successful numerical orbit hypothesis sets `fit_ok=True` in
  `heliolincrr/orbit_confirm_links.py:887-908,942-973`.
- Single-night `is_good` is computed at
  `heliolincrr/orbit_confirm_links.py:920-929` and requires:
  RMS <= configured threshold, maximum residual <= configured maximum,
  `n_used >= 3`, at least two tracklets, and linear RMS <= 10 arcsec.
- The retained-fraction criterion is used only in the multi-night branch at
  `heliolincrr/orbit_confirm_links.py:930-940`.
- CLI defaults are defined at
  `heliolincrr/orbit_confirm_links.py:1080-1086`.
- The unknown catalog gate is explicitly
  `info["fit_ok"] and all_non_asteroid` at
  `heliolincrr/summarize_single_night.py:323-326`; `is_good` is stored but not
  used as the gate.

### Quantified impact

From 128 available nightly summary files:

| Quantity | Count |
|---|---:|
| All links | 87,850 |
| `fit_ok` links | 87,462 |
| `is_good` links | 87,351 |
| all-non-known `fit_ok` links | 7,085 |
| all-non-known `is_good` links | 6,995 |

The observed difference occurred on quarantined night 20251226:

- all-non-known `fit_ok`: 366
- all-non-known `is_good`: 276

For the 124 retained unknown catalogs, all 4,762 rows are both `fit_ok=True`
and `is_good=True`. Thus the current headline candidate list is unchanged by
applying `is_good` post hoc, but the equivalence is accidental and fails on a
known stress night.

### Required handling and rerun boundary

- Methods and Figure 9 must expose two separate stages: numerical `fit_ok` and
  thresholded `is_good`.
- Do not state that 80% retention is a single-night criterion.
- If production continues to select `fit_ok`, this is a **text correction and
  analysis rebuild** only for the current 4,762 rows.
- If production is changed to select `is_good`, rerun every frozen unknown
  night from summarization onward and reconsider quarantine status. If orbit
  fitting or site coordinates also change, rerun the full unknown branch.

## Conflict 3 — 14,299 is not a global unique-detection count

**Severity:** P0 / High for published accounting

**Decision class:** Analysis rebuild; no production rerun.

### Manuscript evidence

- `sections/07_results_unknown.tex:5`: “14,299 unique detection members.”
- `sections/07_results_unknown.tex:21` places the same value in the candidate
  funnel.
- The macro is `paper_numbers.tex:65` (`\NUnknownDet`).

### Code evidence

- `heliolincrr/summarize_single_night.py:342-351` deduplicates endpoint
  detections only within one linkage.
- Different linkages may share an endpoint/detection because the directed link
  graph can branch.

### Quantified impact

Directly aggregating the 124 server unknown FITS catalogs gives:

| Definition | Count |
|---|---:|
| Sum of `n_obs` / linkage–detection memberships | 14,299 |
| Unique `(night, image_name, objID)` detections | 14,159 |
| Repeated memberships across linkages | 140 |

### Required handling and rerun boundary

- Preferred reporting: publish both 14,299 linkage–detection memberships and
  14,159 unique detections, with definitions.
- If only one number is retained, rename 14,299 to “linkage–detection
  memberships”; it must not be called unique.
- Rebuild macros, tables, funnel inputs, and figure annotations from the frozen
  catalogs. No Gaia/tracklet/link/orbit rerun is needed.

## Conflict 4 — Near-Sun selection loses solar-separation priority

**Severity:** P0 / High for the scheduler method claim

**Decision class:** Text correction for historical results; code change affects
future plans only unless the observing history is reinterpreted.

### Manuscript evidence

- `sections/03_scheduler.tex:48`: the scheduler “ranks currently visible fields
  by angular separation from the Sun” and “favor[s] the smallest solar
  elongations.”
- The interpretation propagates to `sections/07_results_unknown.tex:71`,
  `sections/09_discussion.tex:43`, and the scheduler conclusion at
  `sections/10_conclusions.tex:10`.

### Code evidence

- `survey/scheduler.py:186-189` initially sorts `candidate_idx` by separation
  from the Sun.
- `survey/scheduler.py:190` immediately sorts the entire candidate set by RA.
- `survey/scheduler.py:193-203` then chooses a prefix of that RA-sorted array.
  The solar-separation ordering therefore does not determine which fields are
  selected.

### Quantified impact

- The archived planned exposure count and acquired-frame plan-compliance metric
  remain valid descriptions of the generated plans.
- The interpretation “selected smallest available elongations” is unsupported.
- Any attribution of the 58-link timing distribution to near-Sun prioritization
  is weakened.

### Required handling and rerun boundary

- For historical results, describe the actual implementation and compute solar
  elongation directly for archived planned and acquired exposures.
- Figure 4 must use a frozen historical plan and must not illustrate intended
  behavior as though it were executed behavior.
- A code fix changes only future schedules. Historical plans must not be
  regenerated and substituted for the plans actually used.

## Conflict 5 — In-memory history discourages but does not prevent repeat fields

**Severity:** P1 / Medium–High

**Decision class:** Text only plus scheduler-analysis rebuild.

### Manuscript evidence

- `sections/03_scheduler.tex:46`: “This prevents repeated selection of the same
  field later in the same nightly simulation.”

### Code evidence

- Weights depend on history at `survey/scheduler.py:93-120`.
- The in-memory history is incremented at `survey/scheduler.py:243-252`.
- No exclusion mask is added when the next cluster is selected at
  `survey/scheduler.py:313-338`.

### Quantified impact

Among 121 currently archived plan JSON files, 107 contain at least one field
selected in more than one normal cycle, excluding the intended repeats within
one cycle.

### Required handling and rerun boundary

- Replace “prevents” with “discourages by reducing subsequent weight.”
- Scheduler statistics and Figure 4 should retain the actual per-field
  multiplicities.
- No observing or asteroid-processing rerun is required.

## Conflict 6 — Static-source removal is reference-anchored and not exposure-unique

**Severity:** P0 / High for the unknown selection function

**Decision class:** Text correction if preserved; full unknown rerun if fixed.

### Manuscript evidence

- `sections/05_unknown_pipeline.tex:23`: “A source recurring in at least two
  exposures is marked static.”
- `sections/05_unknown_pipeline.tex:78` repeats this as the static-source rule.

### Code evidence

- `heliolincrr/make_tracklet_linreproj.py:356-370` selects the largest nonempty
  catalog as a single reference.
- `:377-393` accumulates all nearest-neighbor hits onto reference rows using
  `bincount`.
- `:395-405` removes sources in every catalog only when they match a reference
  source classified as static.

### Quantified impact

- A source absent from the largest reference catalog but recurring in two other
  exposures is not identified by this algorithm.
- Multiple detections in one comparison exposure can increment one reference
  source more than once, so the counter is not strictly a number of distinct
  exposures.
- The exact number of affected detections/candidates cannot be determined from
  final catalogs alone; it requires replaying the filter.

### Required handling and rerun boundary

- If preserving the algorithm, describe it as a “largest-reference-catalog
  anchored recurrence filter,” not a symmetric two-exposure recurrence test.
- If changing it to a symmetric, per-exposure-unique test, rerun static
  subtraction, tracklets, links, orbit fitting, known subtraction, review/MPC
  reconciliation, and all unknown statistics/figures.

## Conflict 7 — Unknown ADES timestamps use exposure start, not midpoint

**Severity:** P0 / High for submitted astrometry

**Decision class:** Local correction audit and paper-analysis rebuild complete;
user/MPC decision required before generating or transmitting corrected ADES.

### Manuscript evidence

- `sections/05_unknown_pipeline.tex:67-69` describes standards-formatted unknown
  ADES export and production submission but does not disclose a different time
  convention from the known branch.
- `sections/04_known_pipeline.tex:13-17` explicitly defines the ephemeris time as
  the exposure midpoint.

### Code evidence

- Known midpoint: `known_asteroid/match_single_night.py:145-173`.
- Tracklet/orbit midpoint: `heliolincrr/make_tracklet_linreproj.py:162-183`
  adds 15 s to `EXPSTA`, `DATE-OBS`, or `MJD`.
- Unknown ADES reloads the original L2 row at
  `heliolincrr/export_unknown_ades.py:290-305` and assigns `MJD` directly at
  `:350-355,379` without adding half the exposure.
- The analysis-only correction script states its exact-join and no-ADES/no-MPC
  boundary at `scripts/analyze_unknown_ades_time_correction.py:2-10`.

### Quantified impact

- Inspected L2 files have `MJD == DATE-OBS == EXPSTA` and
  `EXPOSURE approximately 30.002 s`; the L2 `MJD` is the exposure start.
- Unknown ADES times are therefore about 15 s early.
- At the pipeline speed limits 3–63 arcsec/hr, 15 s corresponds to an along-track
  displacement of approximately 0.0125–0.2625 arcsec.
- The initial 67 accepted links contain 206 submitted observations; the final
  58-link table contains 179 observations.
- The completed retained-sample audit joins all 179 rows exactly to 164 L2
  files, with no unmatched row, duplicate strict-manifest key, or exposure-time
  fallback (`snapshot/unknown_ades_time_correction/unknown_ades_midpoint_correction_summary.json:21-39,46-52,63-76`).
- The measured midpoint additions are 15.0015 s at the median and 15.003 s at
  maximum (`:53-57`). The corresponding per-link along-track displacement has
  median 0.0617 arcsec, p90 0.1227 arcsec, and maximum 0.2556 arcsec
  (`:15-19`); the production-bound envelope remains 0.0125--0.2626 arcsec.
- The audit records that it generated no ADES/MPC file, contacted no network
  service, and modified no production file (`:64-68`). Figure 11's rebuilt
  frozen products use midpoint epochs; the 15 s issue is therefore no longer a
  pending local twilight-statistics calculation.

### Required handling and rerun boundary

1. Define a single observation-time convention from the upstream schema.
2. Preserve the completed 179-row midpoint-correction table as the
   analysis-only retained-sample audit. No tracklet/link/orbit rerun is needed,
   because those branches already used midpoint epochs.
3. If the author selects local-analysis-only treatment, disclose the historical
   exposure-start ADES behavior and use the frozen corrected analysis products;
   no additional local statistics run is pending.
4. If the author/MPC owner decides that already transmitted rows require
   correction, generate and validate the complete applicable submitted-row
   payload (not merely the retained 179-row subset) under an MPC-approved
   procedure. Do not resubmit automatically.

## Conflict 8 — Unknown link direction is not RA-wrap safe

**Severity:** P0 / High for completeness near RA = 0

**Decision class:** Full unknown rerun if fixed.

### Manuscript evidence

- `sections/05_unknown_pipeline.tex:37-45` describes compatible sky-plane
  motions and a wrapped difference in position angle.
- The manuscript documents an RA-wrap fix for known ephemeris queries at
  `sections/04_known_pipeline.tex:19`, but no corresponding unknown-link
  limitation.

### Code evidence

- `heliolincrr/run_linear_links_from_tracklets.py:128-132` calculates
  `dra_arcsec` from the raw expression `end_ra - start_ra`.
- `:48-50,174-179` wraps only the difference between the two already-computed
  direction angles; it does not wrap the underlying RA displacement.

### Quantified impact

An object crossing 359.999 deg to 0.001 deg is interpreted using an almost
-360 deg RA displacement at the link stage. The number of lost links is not
recoverable from the selected catalog alone and must be measured by rerunning
or a targeted RA-wrap injection/blinded-known test.

### Required handling and rerun boundary

- If not fixed before submission, document it as a selection limitation and
  report recovery by RA-wrap proximity.
- If fixed, rerun linking, orbit confirmation, known subtraction, review/MPC
  reconciliation, and all unknown figures/tables. Tracklet construction need not
  be rerun if its inputs and great-circle pair selection are unchanged.

## Conflict 9 — Daily unknown processing is gated by known reporting, not only mask readiness

**Severity:** P1 / Medium–High for operations and latency

**Decision class:** Text correction or orchestrator change; status/latency rebuild.

### Manuscript evidence

- `sections/08_operations.tex:5` says the wrapper waits for the official 1.5
  arcsec mask or an explicit empty status and then launches the unknown branch.

### Code evidence

- `heliolincrr/run_daily_pipeline.sh:43-68` defines known readiness as either a
  nonempty ADES PSV plus MPC reply, or a completed official 1 arcsec run with
  zero matches.
- `heliolincrr/run_daily_pipeline.sh:101-115` waits for that state before calling
  the unknown wrapper.
- Only then does `heliolincrr/run_daily_unknown.sh:54-94` wait for the 1.5 arcsec
  mask/explicit empty-mask status.

### Quantified impact

A nonzero known night with a valid 1.5 arcsec mask but no production reply is
not eligible for the daily unknown call under the top-level wrapper. This can
inflate latency or create apparently missing unknown nights independently of
scientific mask readiness.

### Required handling and rerun boundary

- Describe both gates exactly, or change the orchestrator so scientific mask
  readiness is independent of MPC submission.
- Existing candidate membership does not need a rerun solely for wording.
- Night-status and latency tables must distinguish mask-ready time, ADES-ready
  time, MPC-reply time, and unknown-start time.

## Conflict 10 — “True WCS footprint” is a nominal detector rectangle

**Severity:** P1 / Medium–High for recovery denominators

**Decision class:** Text correction; new recovery analysis required.

### Manuscript evidence

- `sections/01_introduction.tex:19`: “tests the true WCS footprint.”
- `sections/04_known_pipeline.tex:17,46` refers to detector-footprint selection
  and a future actual-valid-footprint denominator.

### Code evidence

- `known_asteroid/match_single_night.py:25-28` substitutes fixed detector
  dimensions for binary-table HDUs.
- `:300-346` projects the four corners of a 9216 x 9232 rectangle and accepts
  predicted points inside that rectangle.
- No bad-pixel, saturation, missing-region, chip-quality, or effective-search
  mask is applied.

### Quantified impact

- The 534,780 matched detections remain valid as nominal-rectangle associations.
- A predicted-detectable denominator constructed from the same rectangle alone
  will overcount objects in unusable image regions and bias recovery low.

### Required handling and rerun boundary

- Replace “true” with “WCS-projected nominal detector rectangle.”
- Do not rerun the existing association sample solely for wording.
- P1 known recovery must construct an exposure-level effective mask before
  reporting completeness.

## Conflict 11 — Production known matching does not auto-read exposure time

**Severity:** P1 / Low–Medium

**Decision class:** Text only.

### Manuscript evidence

- `sections/04_known_pipeline.tex:17`: “with a 30 s fallback exposure time if the
  header does not provide one.”

### Code evidence

- `known_asteroid/match_single_night.py:163-173` can read a header key only when
  `exptime_key` is explicitly supplied.
- `known_asteroid/slurm_match_one_file.sh:30-32,69-70,87-89` defaults
  `EXPTIME_KEY` to empty and therefore always passes the configured 30 s value.

### Quantified impact

The inspected L2 header uses `EXPOSURE=30.002`; the production 30.0 s midpoint
differs by approximately 0.001 s. This is negligible for the current numbers,
but the implementation description is inaccurate.

### Required handling and rerun boundary

Write: “Production uses a configured 30 s exposure; an optional header-key
override is available.” No statistics need regeneration.

## Conflict 12 — Gaia release and unconditional source cuts are not frozen

**Severity:** P0 / Medium–High for reproducibility

**Decision class:** External provenance required; conditional rerun.

### Manuscript evidence

- `sections/02_facility_data.tex:5-6` describes Gaia DR3/EDR3-era products.
- `sections/05_unknown_pipeline.tex:15` states that Gaia DR3 is used and that
  every source must satisfy `Mag_PSF <= 21` and `Flag == 0`.

### Code evidence

- `heliolincrr/run_single_night.sh:106-117` points to the unversioned directory
  `/pipeline/ref/healpix`.
- `heliolincrr/mask_gaia.py:77-90` applies the magnitude and flag cuts only if
  those columns exist.
- `heliolincrr/mask_gaia.py:142-179` propagates positions from a configured
  epoch, but this does not identify the catalog release.

### Quantified impact

- The server directory contains 12,288 HEALPix files.
- The inspected tile header contains no release/version identifier.
- The impact of conditional cuts is zero only if every frozen L2 catalog has the
  required columns; that schema audit has not yet been documented.

### Required handling and rerun boundary

- Obtain the Gaia source/release provenance and freeze a manifest/hash set.
- Audit all L2 schemas for `Mag_PSF` and `Flag` before making unconditional
  statements.
- If the same files and required columns are confirmed, no rerun is needed.
- A catalog-release change or missing-column discovery requires rebuilding Gaia
  masks and all downstream unknown products.

## Conflict 13 — Follow-up containment uses a bounding circle, not the polygon

**Severity:** P1 / Medium

**Decision class:** Text only unless follow-up is activated.

### Manuscript evidence

- `sections/08_operations.tex:69`: “selects the existing footprint that contains
  the predicted position.”

### Code evidence

- `survey/followup.py:250-260` defines each footprint radius as the maximum
  center-to-corner separation.
- `survey/followup.py:269-277` accepts a target when its center separation is
  within that radius plus 0.05 deg.

### Quantified impact

The selection is a bounding-circle test that includes regions outside the
footprint polygon. The frozen sample has zero follow-up sources, nights, and
exposures, so no published on-sky count changes.

### Required handling and rerun boundary

- Write “within an approximate footprint bounding circle.”
- No frozen-statistics rerun is required.
- Implement polygon containment before claiming validated closed-loop follow-up.

## Conflict 14 — Association identity and ADES-valid identity use different rules

**Severity:** P1 / Medium for accounting

**Decision class:** Analysis/MPC reconciliation rebuild; no rematch required.

### Manuscript evidence

- `sections/04_known_pipeline.tex:30` says the object key uses a permanent number
  or otherwise a cleaned provisional/name field.
- `sections/06_results_survey_known.tex:24` uses that rule for 58,482 distinct
  associated objects.

### Code evidence

- ADES provisional identifiers must match `PROVID_RE` at
  `known_asteroid/export_ades.py:18,39-43`.
- `known_asteroid/export_ades.py:238-275` exports only rows with a valid permanent
  or strict provisional identifier.

### Quantified impact

The 58,482 association identities are not necessarily all ADES-exportable
identities. The draft already reports fewer ADES objects, but the identity
definitions are not explicitly separated in the method/table schema.

### Required handling and rerun boundary

- Define “association identity key” separately from ADES `permID/provID`.
- Recompute reporting-accounting columns from existing matched/PSV products.
- No ephemeris query or catalog rematch is required.

## Conflict 15 — Twilight concentration cannot be attributed to near-Sun mode

**Severity:** P1 / Medium–High for interpretation

**Decision class:** Analysis rebuild and text qualification.

### Manuscript evidence

- `sections/07_results_unknown.tex:71`: the timing concentration is said to be
  consistent with the near-Sun end-of-night mode.
- `sections/09_discussion.tex:43-45` discusses twilight/near-Sun potential.

### Code evidence

- `survey/scheduler.py:186-203` shows that near-Sun selection did not retain
  solar-separation priority; the selected prefix comes from an RA-sorted array.
- The draft itself states that the exposure-area/depth denominator is absent.
- The retained sample contains 38 dusk-nearest and 20 dawn-nearest linkages;
  near-Sun mode is an end-of-night mode and cannot explain the dusk majority by
  itself.

### Quantified impact

- 27/58 are within 1 hr and 37/58 within 2 hr of an astronomical-twilight
  boundary.
- These fractions describe selected candidates, not survey efficiency or mode
  yield.

### Required handling and rerun boundary

- Remove causal language until every exposure has solar altitude, elongation,
  depth/background, plan mode, and searched-area metadata.
- Figure 11 must include the all-exposure denominator; otherwise keep it as
  descriptive appendix/context.
- No production rerun is required, but all twilight statistics must use the
  corrected midpoint-time convention from Conflict 7.

## Conflict 16 — Known matching still has a `Mag_Kron` schema/finite-value gate

**Severity:** P1 / Medium for recovery accounting

**Decision class:** Schema audit first; conditional rematch.

### Manuscript evidence

- `sections/04_known_pipeline.tex:28` says the 99 mag tolerance is effectively
  disabled and the principal association criteria are footprint and angular
  separation.

### Code evidence

- `known_asteroid/match_single_night.py:361-369` rejects all matches in a frame
  if the configured measured-magnitude column is absent.
- `known_asteroid/match_single_night.py:370-379` still requires predicted and
  measured magnitudes to be finite before either the 1 arcsec or 1.5 arcsec
  match is retained.

### Quantified impact

The numerical effect has not yet been counted across the frozen L2 manifest.
It is zero only if `Mag_Kron` exists and is finite for every otherwise matched
source. Therefore the known recovery denominator and failure taxonomy must
include this schema/photometry gate.

### Required handling and rerun boundary

- Qualify the text: the magnitude-difference bound is loose, but availability
  and finiteness of `Mag_Kron` remain required.
- Audit all frozen L2 schemas and per-frame skipped-match logs.
- If the schema is uniform, no rematch is required. If valid geometric matches
  were lost solely through missing/invalid `Mag_Kron` and the policy changes,
  rerun known matching, 1.5 arcsec masks, and all downstream unknown products.

## Conflict 17 — Known recovery must deduplicate prediction rows before forming the denominator

**Severity:** P1 / High for recovery accounting

**Decision class:** Frozen analysis/table/figure rebuild complete; manuscript
consumption and denominator wording remain. No production ephemeris or matching
rerun.

### Manuscript evidence

- `sections/04_known_pipeline.tex:46-52` requires a future detectable-prediction
  denominator and a stage-wise known-object recovery measurement.
- `sections/06_results_survey_known.tex:63-64` reserves Figure 8(c) for recovery
  versus predicted magnitude/rate and says that the denominator has not yet been
  constructed.
- `sections/09_discussion.tex:51-52` says to start “for every predicted object.”
  That wording requires one row per object per exposure, not one row per raw
  ephemeris-query return.

### Local code and frozen-product evidence

- `snapshot/frozen_products/row_counts.json:12-20` records 13,311,871 source
  rows in the 131 available `known_all` products.
- `scripts/analyze_known_population.py:281-292` builds a canonical object key,
  combines it with `source_file`, counts duplicate prediction keys within each
  night, and retains rank zero before any recovery calculation.
- `scripts/analyze_known_population.py:499-506` writes both the raw source-row
  count, deduplicated nominal-prediction count, duplicate count, and recovery
  fraction to the analysis summary.
- The completed summary freezes 13,310,546 nominal predictions, 1,325 duplicate
  keys, 534,780 official matches, and the 4.017716% nominal match fraction at
  `snapshot/derived_known/known_population_summary.json:17-20,32-38`.
- Table 3 already exports the raw, deduplicated, matched, fraction, and duplicate
  rows (`tables/table3_known_recovery_astrometry.csv:2-8`); Figure 8 consumes
  the rebuilt binned recovery products.
- A direct read-only pass over the frozen parquet, using the same canonical
  object-key rule, found 13,310,546 unique `(night, source_file, object_key)`
  predictions and 1,325 duplicate prediction rows. All rows had a usable number
  or name key; none collapsed into the `UNKNOWN` fallback.

### Quantified impact

| Prediction grain | Count |
|---|---:|
| Raw `_all_asteroids` source rows | 13,311,871 |
| Unique `(night, exposure, object)` predictions | 13,310,546 |
| Duplicate prediction rows removed | 1,325 (0.00995%) |

The 534,780 matched rows have no duplicate `(night, exposure, object)` keys
under the same audit. Using the raw source rows would therefore dilute the
overall nominal match fraction from 4.017716% to 4.017317%, a difference of
0.000400 percentage points. The numerical effect is small globally but can be
larger in individual bins or nights and, more importantly, changes the unit of
analysis. A recovery fraction is not auditable if its denominator mixes raw
query returns with unique exposure-level object opportunities.

The 13,310,546 value remains a **nominal footprint denominator**, not a complete
detectable denominator: image masks, missing catalog regions, saturation,
ephemeris uncertainty, and a calibrated limiting-magnitude model are still not
fully represented.

### Required handling and rerun boundary

- Define the primary recovery key as `(night, exposure/source_file,
  canonical_object_key)` and publish both raw and deduplicated counts.
- Use 13,310,546, not 13,311,871, as the denominator before applying any further
  exposure-quality, magnitude, rate, edge, or detectability cuts.
- The known-recovery table, bins, intervals, and Figure 8 have already been
  rebuilt from frozen products. The manuscript and any remaining macros must
  consume those products rather than the raw 13,311,871-row denominator.
- No `aleph` query, known-object rematch, 1.5 arcsec mask rebuild, or unknown
  pipeline rerun is required unless the *upstream prediction membership* itself
  is changed rather than merely deduplicated for analysis.

## Conflict 18 — The 58 retained links are neither 58 independent sources nor an all-unknown sample

**Severity:** P0 / High for the central unknown-object claim

**Decision class:** User decision and MPC reconciliation; analysis rebuild.

### Manuscript evidence

- `frontmatter.tex:16` correctly calls the 58 rows candidate linkages and says
  that cross-night association and final MPC reconciliation are incomplete.
- `sections/07_results_unknown.tex:13-15` reports 58 post-audit linkages and
  explicitly warns that some may be repeated observations or known-object
  recoveries.
- `sections/07_results_unknown.tex:25` labels the 58-row stage “Single-night
  linkages, not confirmed discoveries.”
- `sections/10_conclusions.tex:14-16` repeats the 58-linkage result but leaves
  individual cross-night identity and designation unresolved.
- `AUTHOR_TODO_AND_ANALYSIS_PLAN.md:83-103` asks for a final count of independent
  sky-plane sources. A linear-motion connected-component screen cannot, by
  itself, supply that physical-object count.

### Local linear-motion evidence

- `scripts/analyze_cross_night_repeat_candidates.py:82-87` fixes the primary
  thresholds at maximum two-way residual <= 0.03 deg, speed difference <=
  2 arcsec hr^-1, and direction difference <= 5 deg.
- `scripts/analyze_cross_night_repeat_candidates.py:368-470` evaluates every
  `0 < delta_t <= 7 d` pair by propagating each link's own linear motion forward
  and backward in its local tangent plane.
- `snapshot/derived_unknown/unknown_linear_motion_candidate_group_summary.json:3-12`
  explicitly prohibits interpreting a component as a confirmed independent
  object, orbit, or MPC/JPL identity; components use transitive closure.
- The primary result has 37 components: 25 singletons and 12 non-singletons
  containing 33 links. It is supported by 30 direct edges (25 cross-night and
  five same-night); see the summary at `:1158-1160` and `:1277-1297`.
- Threshold sensitivity is material: the strict test gives 41 components and
  27 links in 10 non-singletons, while the relaxed test gives 35 components and
  36 links in 13 non-singletons
  (`unknown_linear_motion_candidate_group_summary.json:1252-1322`). Thus even
  “37” is a screening result, not an independent-body count.

### Local JPL evidence and authority boundary

- The first frozen all-58 query used JPL's official `sb_ident` API and returned
  six strict numerical candidate associations
  (`snapshot/jpl_identification/jpl_identification_summary.json:2-14`).
- A separately executed targeted pass again queried the official API with a
  0.03 deg half-width and its own saved raw responses
  (`scripts/confirm_jpl_candidates_second_pass.py:35-59,62-93`). Its numerical
  confirmation requires the same object name, <=1 arcsec separation,
  <=1 arcsec hr^-1 speed difference, and <=5 deg direction difference
  (`:95-145`).
- All six second-pass rows uniquely return **C/2025 Y1 (ATLAS)** and pass those
  numerical cuts. Their separations are 0.0981--0.5465 arcsec, speed differences
  0.0123--0.7910 arcsec hr^-1, and direction differences 0.0596--1.6276 deg
  (`snapshot/jpl_identification/second_pass/jpl_second_pass_confirmations.csv:2-7`).
- The six rows are exactly two primary linear-motion components:
  - `LMG027`: 20260309 `000011l/35`, 20260310 `000012H/1026`, and
    20260311 `000013V/1139`;
  - `LMG036`: 20260530 `00001gL/73`, 20260530 `00001gM/75`, and
    20260601 `00001h8/82`.
- The raw-response hashes from both runs are frozen in
  `snapshot/jpl_identification/hashes.sha256` and
  `snapshot/jpl_identification/second_pass/hashes.sha256`.
- This establishes a strong, repeatable **post-hoc JPL ephemeris association**.
  It does **not** establish what the MPC ingested from each SHARP submission,
  whether every submitted row was accepted or attributed, or whether SHARP
  received/owns any designation. The second-pass summary preserves that
  guardrail at
  `snapshot/jpl_identification/second_pass/jpl_second_pass_summary.json:7-12`.

### Quantified impact

- 33/58 links (56.9%) participate in a non-singleton primary component; 25/58
  remain singleton components under the chosen threshold.
- Six of 58 retained links (10.3%) are numerically associated twice with one
  known comet and occupy two of the 37 candidate components.
- If the author classifies those six as post-hoc known-comet recoveries and
  excludes them from an “unattributed at final analysis” subset, that subset is
  52 linkages in 35 screening components. It is still not a discovery or
  independent-object count because the remaining identity states are unresolved
  and the component count is threshold-dependent.

### Required handling and rerun boundary

1. Preserve 58 as the count of post-audit single-night linkages only, with a
   separate post-hoc identity/status column.
2. Decide whether C/2025 Y1 rows remain in descriptive “pipeline output” plots,
   move to a “known-object leakage/recovery” category, or are excluded from an
   explicitly final-unattributed subset. Report both historical-at-run and
   post-hoc identity states if they differ.
3. Obtain the submitter's row-level MPC records for all 58 links, especially the
   six comet-associated links, before populating ingest, attribution,
   identification, rejection, duplicate, or designation fields.
4. Never replace “58 linkages” with “58 independent sky-plane sources,” “58 new
   asteroids,” or a discovery count. Do not promote 37 screening components to
   an object count either.
5. Rebuild the 58-row candidate table, status matrix, unknown results table,
   Figures 3/9/10, abstract/conclusion counts, and machine-readable release if
   the author creates a 52-row final-unattributed subset. No production rerun is
   required for this post-hoc relabeling.
6. Only if the author changes the contemporaneous known-object database or
   matching policy and asks for a counterfactual pipeline result is a known-mask
   and downstream unknown rerun needed; preserve the historical-as-executed
   result separately.

## Conflict 19 — “First review retained 67” collapses a 68 -> 67 -> 58 decision chain

**Severity:** P0 / Medium–High for review provenance and submitted-sample
accounting

**Decision class:** Text correction and frozen analysis/table/figure rebuild;
no production rerun.

### Manuscript evidence

- `frontmatter.tex:16` says the automatic branch produced candidate linkages “of
  which `NHumanPass` were initially approved and submitted,” treating review and
  submission as one gate.
- `sections/05_unknown_pipeline.tex:67-69` correctly says that only final
  `<night>_submit.csv` rows with `is_real=1` reach the reviewed exporter, but
  then calls `NHumanPass` the number approved by “the first human review.”
- `sections/07_results_unknown.tex:11,22` calls `NHumanPass` both the first-review
  retained count and the count entering production submission files.
- `sections/10_conclusions.tex:14` again says first human review retained
  `NHumanPass`. With the current macro value of 67, this describes the
  submission-selected set, not every row initially marked real in the review
  CSVs.

### Frozen evidence

- `scripts/build_night_tables.py:144-170` preserves the production-summary
  fields `review_real` and `submit_real` as separate per-night columns;
  `:204-219` independently constructs the frozen audit-initial and retained
  sets from the row-level review ledger.
- `snapshot/tables/unknown_stage_counts_total.json:2-20` closes the four relevant
  totals: `review_real_n=68`, `submit_real_n=67`, `audit_initial_n=67`, and
  `audit_real_n=58`.
- The only per-night review/submission mismatch is 20260507:
  `snapshot/tables/unknown_stage_counts_by_night.csv:175` has
  `review_real_n=1` and `submit_real_n=0`.
- The frozen taxonomy identifies the exact withdrawn row as
  `20260507/00001et`
  (`snapshot/tables/unknown_false_positive_taxonomy.csv:3`). This withdrawal
  occurs **before** the 67-row audit partition and is not one of the nine
  post-hoc rejects.
- The 67-row frozen ledger partitions exactly into 58
  `retained_after_posthoc_audit` and nine `rejected_posthoc` rows
  (`reports/VALIDATION_REPORT.md:116-119`). Figure 9 now exports the three
  plotted gates and the evidence note at
  `figure_data/fig09_unknown_funnel.csv:16-19`.

### Quantified impact

| Gate | Linkages | Interpretation |
|---|---:|---|
| Review CSV marked real | 68 | Initial human `is_real=1` decisions before submission reconciliation |
| Submission-selected / frozen audit initial set | 67 | One pre-submission withdrawal: `20260507/00001et` |
| Post-hoc retained | 58 | 67 initial audit rows minus nine later artifact rejects |

Thus “67 initially approved” is only defensible when explicitly defined as
**submission-selected** or **frozen audit-initial**, while “68 marked real in
the first review CSVs” is the correct pre-submission review count. Neither 68,
67, nor 58 is an independent-object or discovery count.

### Required handling and rerun boundary

1. Use three explicit labels in prose, tables, macros, and figures:
   `review-marked real (68)`, `submission-selected/frozen audit initial (67)`,
   and `post-hoc retained (58)`.
2. State that `20260507/00001et` was the unique review-real row withdrawn before
   submission. Do not merge it with the nine post-hoc artifact rejects.
3. Keep submission/MPC accounting anchored to the 67-row frozen ledger, while
   review-protocol accounting begins at 68.
4. Rebuild affected manuscript macros, the unknown funnel/table, and review
   prose from the frozen stage/taxonomy products. Figures 3 and 9 and Tables 4
   and 5 already preserve the separated gates in the analysis project; no
   Gaia/tracklet/link/orbit/known-subtraction production rerun is required.
5. Row-level MPC ingest/attribution/designation evidence remains an external
   requirement for the 67 submission-selected rows; this review-count
   reconciliation does not supply it.

## Statements verified against current code

The following draft statements were checked and are consistent with the
paper-frozen/current code, subject to the separate provenance tasks:

- Fixed 35 s slew overhead.
- Official 1 arcsec known association versus 1.5 arcsec unknown-source mask.
- Production `--skip-common-area`, making the optional 30-pixel common-mask
  erosion and 500-pixel endpoint veto inactive.
- Active 300–500 pixel shell filter conditional on `Flag_ISO_Num > 0`.
- Tracklet speed interval 3–63 arcsec/hr.
- Link speed/direction-difference thresholds of 5 arcsec/hr and 10 deg.
- Quarantine when post-known unknown links exceed 200.
- Follow-up had zero on-sky executions in the frozen sample.
- The 58 retained rows must remain “high-confidence single-night linkages,” not
  confirmed discoveries or 58 independent objects.

## Recommended resolution order

1. Obtain author decisions in `AUTHOR_INPUTS_REQUIRED.md`.
2. Freeze the night/file/quality manifests and exact external dependencies.
3. Resolve the unique site and time conventions.
4. Decide whether near-Sun, static filtering, RA wrapping, and the
   `fit_ok/is_good` catalog gate are code changes or documented historical
   behavior.
5. Run the minimum sensitivity tests first; trigger the stated downstream rerun
   boundary only when classifications or memberships change.
6. Rebuild all paper tables and figures from frozen machine-readable products,
   keeping detection, membership, tracklet, linkage, object, exposure, and night
   units separate.
