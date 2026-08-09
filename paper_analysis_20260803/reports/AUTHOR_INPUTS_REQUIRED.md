# Author and External Inputs Required

## Purpose

Most frozen accounting, statistics, and plotting can be completed from the
repository and server products. The items below cannot be resolved safely by
code inspection alone. They require an author decision, observatory/upstream
pipeline information, reviewer records, or external MPC state.

No production algorithm, quality mask, or MPC submission should be changed on
the basis of an assumed answer.

## Current completion boundary (2026-08-09)

The following work is already complete locally and should not be described as
awaiting a run:

- frozen raw/L2/product inventories, file ledgers, hashes, and a provisional
  machine-readable night-quality mask;
- known-object population, signed-residual, random-shift, and nominal-footprint
  recovery analyses;
- the full unknown-stage funnel, including the separated 68 -> 67 -> 58 human
  gates, true-zero versus not-run nights, and detection/linkage units;
- the 58-row audit table, 37-component linear-motion screen, and two independent
  JPL numerical checks of the six C/2025 Y1-associated rows;
- the complete 128-night observer-location sensitivity comparison;
- the local 179-observation unknown-time midpoint-correction analysis;
- the retrospective blinded-known survival proxy, scheduler analysis, and
  machine-evidence operations analysis; and
- Tables 1--5 plus all 12 figures as both PDF and PNG.

The genuine external blockers are author/observatory approval of the site and
quality mask; upstream L1/L2 method and depth/mask metadata; row-level MPC
records and the C/2025 Y1/time-correction policies; authoritative human/MPC,
weather/fault, and resource logs; and an approved catalog/image injection
design and execution route. The assembled snapshot is otherwise technically
closed: 555 hashed artifacts report `complete`, all 12 figures pass the contact
sheet QA, and the independent validator returns 63 PASS, 0 FAIL, and one SKIP.
That sole SKIP is the unsigned author quality-mask decision.

## Immediate decision checklist

The following decisions block the largest number of downstream tasks:

- [ ] Confirm the single surveyed MPC 327 longitude, latitude, height, and
  geodetic datum.
- [ ] Decide whether the paper reports the historical production algorithms as
  executed or first changes and reruns them.
- [ ] Approve the frozen night-quality mask and primary exclusion reason per
  night.
- [ ] Confirm the unknown ADES observation-time correction policy and whether
  any MPC correction is required.
- [ ] Supply/finalize row-level MPC ingest/accept/reject/attribution/designation
  records for the 58 retained links and nine rejected links; a stored reply or
  a JPL positional match is not a substitute.
- [ ] Decide how the six links independently associated twice by the JPL API
  with C/2025 Y1 (ATLAS) are classified in the historical-output sample and in
  any final-unattributed subset.
- [ ] Identify the upstream L1/L2 owner and provide the reduction/version and
  image-quality metadata needed for recovery analyses.

## P0-1 — Frozen data snapshot and quality mask

### Input required from the authors

1. Approve one immutable observation interval and cutoff timestamp.
2. Approve one primary `quality_code` and any secondary flags for every night.
3. Decide the scientific treatment of at least:
   - 20260111 (`WCS_OFFSET`);
   - 20251226, 20260111, 20260528, and 20260611 (high-candidate/quarantine
     nights);
   - 20260209, 20260312, and 20260413 (raw science but no current L2);
   - nights with matched detections but no ADES;
   - nights with ADES but no recorded production reply;
   - eight valid zero-unknown nights versus nights not run.
4. State whether operational/failure plots show all nights while headline
   scientific summaries use only `GOOD` nights.

### Why code cannot decide this

The files establish what exists, but scientific inclusion/exclusion is an
authorial policy. A night must not disappear merely because a downstream product
is absent.

### Current local status

The inventories, hashes, `night_status.csv`, provisional `quality_mask.csv`,
headline denominators, Table 2, and Figures 2, 3, 7, 8, and 9 have already been
generated. The remaining external step is approval of the scientific mask and
its primary reason codes; approval may change the primary-sample aggregation
but does not require rediscovering which files exist.


## P0-2 — Production configuration truth

### A. Unique site coordinates

**Required owner:** observatory/telescope authority.

Provide:

- longitude and sign convention;
- latitude;
- height and whether it is ellipsoidal or orthometric;
- geodetic datum/reference;
- official MPC 327 position if it differs from the local survey value.

The L1/L2 headers cannot be used: sampled files contain the wrong GMG
coordinates `100.0313, 26.6974, 3227 m`.

### B. Historical behavior versus code changes

For each row, select one policy:

| Topic | Keep and document historical behavior | Fix and rerun | Author decision |
|---|---|---|---|
| Orbit site `868.221 m / 40.394239 deg` | Disclose completed counterfactual and element-instability boundary | Adopt the surveyed site; repeat sensitivity and rerun only if classifications/products change | [ ] |
| Unknown catalog gate | Keep `fit_ok`, show `is_good` separately | Select `is_good`, rerun downstream | [ ] |
| Near-Sun sorting | Describe historical RA-prefix behavior | Fix for future plans; do not rewrite history | [ ] |
| Static-source removal | Describe reference-anchored algorithm | Make symmetric/per-exposure unique; full unknown rerun | [ ] |
| RA=0 unknown linking | Document completeness limitation | Fix wrap and rerun linking/downstream | [ ] |
| Common-area/edge policy | Keep `--skip-common-area` and active shell | Enable common mask/erosion/500 px veto; full unknown rerun | [ ] |
| Known `Mag_Kron` gate | Keep and disclose | Change match policy; known+unknown rerun | [ ] |

The observer-location sensitivity calculation is complete for all 128 frozen
candidate nights. It compares the historical orbit observer
`117.575 E, 40.394239 N, 868.221 m` with
`117.575 E, 40.393 N, 960 m`; it is not a height-only experiment. All 87,850
link keys matched, with zero `fit_ok` or `is_good` flips overall, in the formal
4,762, or in the retained 58. This closes the tested classification question,
but not the authoritative surveyed-site decision: a few underconstrained
short-arc orbital elements change drastically and must not be used for
population interpretation.

### C. Photometric/reporting convention

Confirm with the upstream/ADES owners:

- whether “unfiltered,” header `FILTER=W`/`BAND=w`, and ADES `band=G` are the
  intended descriptions of the same system;
- whether `Gaia3E` is valid for both `astCat` and `photCat`;
- whether a color term or explicit instrumental-band label is required;
- whether past submissions used this convention consistently.

### D. External dependency versions

Provide or approve the exact versions/dates of:

- `astorb.dat`;
- Gaia HEALPix source and release;
- JPL ephemeris files;
- `aleph`;
- `poliastro` and the `heliolinc` Python environment;
- default Python environment for survey/known processing;
- MPC validation/production interface behavior during the snapshot.

Local production-file hashes and the observable remote environment inventory
have already been frozen. Author/upstream confirmation is still needed for the
semantic release/date and provenance of evolving catalogs, ephemerides, and
external interfaces that cannot be inferred from a filename or installed
package alone.

## P0-3 — Authorship, facility, and upstream L1/L2 pipeline

### Authorship and acknowledgments

Provide:

- final author order;
- affiliations;
- ORCIDs;
- corresponding author;
- CRediT roles;
- observers, measurers, reviewers, and software/infrastructure contributors;
- funding programs and grant numbers;
- institutional acknowledgments;
- official English telescope and station names.

### Upstream pipeline documentation

The upstream owner must provide or cite:

- bias/overscan and flat-field procedure;
- detector nonlinearity, saturation, bad-pixel, cosmic-ray, and artifact masks;
- astrometric solver, reference catalog, distortion model, and residual metrics;
- photometric calibration/reference and passband convention;
- source-extraction software/configuration;
- definitions of `RA_Win`, `RA_PSF`, uncertainty columns, `Flag`, and
  `Flag_ISO_Num`;
- units for coordinate uncertainties;
- seeing/FWHM, sky background, limiting magnitude, source density, and catalog
  completeness per exposure or night;
- pipeline version/commit and production dates;
- explanation and correction plan for the wrong GMG site metadata in L1/L2
  headers.

### Current local status and remaining boundary

Table 1 and the available-data versions of Figures 8 and 11 have already been
generated. They preserve the nominal-footprint denominator and observable
header context. A true detectable-object denominator, calibrated completeness
claim, image-level injection, and final facility/upstream prose remain blocked
on the metadata and ownership information above.

## P0-4 — The 58 retained links, nine rejects, and MPC/JPL status

### Required reviewer/author records

The frozen human-decision chain is now resolved numerically: 68 rows were
marked real in the review CSVs, `20260507/00001et` was the only row withdrawn
before submission, 67 rows entered the submission-selected/frozen audit set,
and 58 survived the later audit. The withdrawal is separate from the nine
post-hoc rejects.

For every one of the 67 submission-selected/frozen-audit links, provide or
confirm:

- original review label;
- audit label;
- reviewer identity or role;
- audit reason code;
- final paper disposition;
- whether the candidate was submitted;
- whether any correction/retraction was communicated after the nine artifacts
  were identified.

For `20260507/00001et`, separately confirm the withdrawal reason, who made the
decision, whether any observation row had already left the local pipeline, and
whether an MPC-side action exists. Review-protocol accounting starts from 68;
submission/MPC accounting starts from the 67-row frozen ledger.

Recommended reason codes include detector artifact, bright-star structure,
edge artifact, stationary residual, inconsistent morphology, duplicate link,
known-object leakage, and uncertain.

### Required MPC/external state

For every link/submission, obtain:

- submission ID;
- row-level ingest/accept/reject state where available;
- MPC identification or provisional designation;
- known-object identity if attributed;
- duplicate/previous-observation state;
- final query timestamp and source.

Access may require the submitter's MPC records or credentials; a stored HTTP
reply alone is insufficient.

### Cross-night identity decision

The completed linear-motion screen gives 37 threshold-dependent candidate
components. Twelve are non-singletons containing 33 links; the remaining 25
are singletons. Six links in two of those components (`LMG027` and `LMG036`)
were then returned uniquely as C/2025 Y1 (ATLAS) in both the frozen all-sample
JPL query and an independent targeted second pass. These are reproducible
ephemeris associations, not MPC row-level ingest or designation records.

The author must approve:

- positional/motion uncertainty model;
- maximum time separation;
- whether MPC/JPL identity overrides a purely geometric cluster;
- wording for unresolved groups.

The author must also choose one explicit C/2025 Y1 reporting policy:

- [ ] Keep all 58 in historical pipeline-output plots, label the six rows as
  post-hoc known-comet associations, and provide a separately filtered
  final-unattributed subset.
- [ ] Move the six rows to a known-object leakage/recovery category, giving
  52 final-unattributed linkages, while preserving the original 58-row table for
  audit.
- [ ] Use another documented policy supplied by the authors/MPC owner.

For the six rows, supply or confirm the submission ID and row-level MPC state
for each of:

| Night | `trk_sub` | Linkage ID | JPL numerical association |
|---|---|---:|---|
| 20260309 | `000011l` | 35 | C/2025 Y1 (ATLAS) |
| 20260310 | `000012H` | 1026 | C/2025 Y1 (ATLAS) |
| 20260311 | `000013V` | 1139 | C/2025 Y1 (ATLAS) |
| 20260530 | `00001gL` | 73 | C/2025 Y1 (ATLAS) |
| 20260530 | `00001gM` | 75 | C/2025 Y1 (ATLAS) |
| 20260601 | `00001h8` | 82 | C/2025 Y1 (ATLAS) |

Required author-supplied fields remain: submitted yes/no, submission ID,
row-level accepted/rejected/duplicate state, MPC attribution/identification,
designation (if any), timestamp/source of the final check, and the policy for
whether the object was considered unknown at observation/processing time versus
known at final analysis time.

Until complete, report 58 single-night linkages, not 58 independent bodies or
discoveries. Neither the 37 linear-motion components nor the 52 links remaining
after a possible comet reclassification may be presented as a confirmed object
or discovery count.

### Unknown 15 s timestamp correction

The local correction analysis is complete for all 179 retained observation
rows (58 linkages, 164 L2 files): every row joins exactly to an L2 exposure,
the midpoint additions are 15.0015 s at the median and 15.003 s at maximum,
and no production or MPC/ADES file was generated or modified. What remains is
an author/MPC policy decision, not another local calculation.

The author/MPC owner must choose one:

- [ ] Accept the completed local correction table for analysis only and
  document historical submitted times.
- [ ] Ask MPC whether corrected observations should be submitted.
- [ ] Submit corrections using an MPC-approved procedure.

Codex must not choose or execute the latter two without explicit authorization.

## P1-1 — Known-object recovery and astrometric performance

### Author/upstream definitions required

- limiting magnitude or probabilistic detectability model;
- acceptable ephemeris-uncertainty threshold;
- effective detector/search mask, including bad/saturated/missing regions;
- exposure/night quality criteria;
- whether measured magnitude, predicted V, or both define photometric
  detectability;
- treatment of trails and high-rate objects;
- seeing, sky background, Moon, source-density, and edge-distance metadata;
- confidence-interval convention and minimum bin size.

### Completed locally under the nominal-footprint convention

The deduplicated nominal prediction table (13,310,546 rows), 534,780 official
matches, signed `dRA*cos(dec)`/`dDec` residuals, detector/night/magnitude/rate
trends, binned intervals, and deterministic random-shift control have already
been generated. The random-shift control contains 257 sub-arcsecond matches in
352,184 shifted prediction trials. These products support a clearly labeled
nominal-footprint analysis; they do not create a true detectability denominator.

The limiting-magnitude model, effective search mask, ephemeris-uncertainty
policy, and related upstream inputs listed above remain necessary before the
paper can call the denominator “detectable” or interpret the nominal 4.02%
match fraction as survey completeness.

## P1-2 — Full unknown-stage funnel

### Author decision required

Approve the frozen configuration and quality mask before aggregation. In
particular, decide whether funnel stages include both `fit_ok` and `is_good` and
how quarantined/skipped groups are displayed.

### Current local status

The full frozen funnel and Figure 9 are complete. They distinguish detections,
tracklets, linkages, linkage--detection memberships, and globally unique
detections; expose both `fit_ok` and `is_good`; retain quarantined/not-run
states; distinguish 116 non-empty and eight true-zero unknown-catalog nights;
and show 68 review-marked real -> 67 submission-selected -> 58 post-audit
retained, including the unique `20260507/00001et` withdrawal. Author approval
is still required for which quality cohort and gate wording becomes primary.

### Required unit definitions

The author should approve these labels:

- L2 rows/detections;
- post-cut detections;
- post-Gaia/static detections;
- two-detection tracklets;
- shared-endpoint links;
- `fit_ok` links;
- `is_good` links;
- all-non-known links;
- first-review links;
- post-audit links;
- unique detections versus linkage–detection memberships;
- nights.

These definitions are already encoded in the frozen stage-definition table and
figure-data exports. Approval may change labels or the primary cohort, but the
stage computation itself is no longer pending.

## P1-3 — Selection function and injection tests

### Blinded-known experiment decisions

Approve:

- sample size;
- stratification by magnitude, rate, direction, detector radius, density, and
  observing conditions;
- holdout by night to avoid leakage;
- whether known objects are withheld only from final subtraction or from other
  truth-assisted stages;
- success definitions at L2, tracklet, link, orbit, and review stages.

### Completed retrospective proxy and its limit

A frozen retrospective identity-blind proxy has already been computed
conditional on detections that entered the 1.5-arcsec known mask. Among
102,347 object-night cases with at least three detections and rates in the
production interval, 75,899 (74.16%) survive as strict single-object
`fit_ok/is_good` links. This is not an L2 source-detection completeness test,
does not separate all intermediate losses, and does not replace a controlled
injection experiment.

### Catalog-level injection

This can be implemented locally after the author approves the grid and the
frozen code/configuration. At minimum cover 15–21 mag, 3–63 arcsec/hr including
boundary stress tests, all motion directions, RA=0 crossings, detector edges,
source density, and ordinary/pathological nights.

### Image-level injection

This requires the upstream image-processing owner to provide:

- an injection interface or permission to run modified images through L1/L2;
- PSF/trail model;
- gain/noise/saturation characterization;
- safe storage and naming conventions;
- a no-MPC/dry-run execution route.

Without image-level access, the project may report catalog-stage completeness
but not end-to-end source-detection completeness.

## P1-4 — Latency and resource use

### External/manual inputs required

- timezone convention for every log source;
- authoritative L2-ready event definition;
- human-review start/end or final-decision timestamps;
- MPC response versus final-ingest timestamps;
- server reboot/maintenance intervals;
- any logs outside the project tree;
- any independent CPU/RAM accounting source, because Slurm accounting was
  disabled for the available runs and those quantities cannot be reconstructed
  from wall-clock logs.

### Completed from available machine evidence

Known/unknown event evidence, per-night latency rows, and p50/p90/p95 summaries
have already been generated where a normal-daily timing chain exists.
Historical reruns, latest-mtime-only nights, negative intervals, and restart
anomalies are separated rather than pooled into the normal-night percentiles.
Figure 12 explicitly marks CPU and peak RAM unavailable; neither value was
invented.

Automated processing, human review, and MPC delay must remain separate.

## P1-5 — Scheduler performance beyond plan compliance

### External observatory inputs required

- weather closure and cloud logs;
- dome/telescope/instrument fault intervals;
- manual override records;
- telescope-control execution logs and actual slew/readout overheads;
- reason codes for planned but unacquired frames;
- upstream-loss versus non-observation distinction.

### Reporting boundary if inputs are unavailable

The acquired-frame plan-compliance, planned-versus-actual field multiplicity,
start-time residual, revisit-cadence, and exposure-elongation analyses and
Figure 4 are complete. The 95.17% value is defined only on plan-active acquired
nights. Without the listed external logs, it must not be called observatory
efficiency or full plan completion.

## Figure-specific inputs for final acceptance

All 12 figures already exist as PDF and PNG. The table below lists remaining
author/external decisions that can alter interpretation, labeling, or panel
scope; it is not a list of figures still awaiting generation.

| Figure | Author/external input that cannot be inferred safely |
|---|---|
| 1 System architecture | Confirm ownership boundary and which components were exercised on sky. |
| 2 Footprint coverage | Approve final exposure/night mask. |
| 3 Data accounting | Approve all exclusion/status categories and MPC interpretations. |
| 4 Scheduler example | Approve a representative historical night; decide how to label the near-Sun bug. |
| 5 Known residuals | Approve detectable denominator and chance-match protocol. |
| 6 Real/artifact examples | Approve publication of selected cutouts/GIFs and labels. |
| 7 Exposure history | Approve quality-mask display. |
| 8 Known results | Supply upstream seeing/depth/background or accept a reduced panel set. |
| 9 Unknown funnel | Approve `fit_ok`/`is_good`, quarantine, and membership/unique units. |
| 10 High-confidence distributions | Approve final 58/9 sample and MPC status freeze. |
| 11 Twilight context | Supply/approve exposure-level searched-area/depth denominator. |
| 12 Operations | Supply human, MPC, weather/reboot, and missing Slurm timing context. |

## Minimum author response format

To unblock work efficiently, the authors can return one concise record with:

```text
SITE: lon=..., lat=..., height=..., height_type=..., datum=..., source=...
QUALITY_MASK: approved file/version=...
HISTORICAL_ALGORITHM_POLICY: document-as-run | fix-and-rerun
FIT_GATE: fit_ok | is_good
STATIC_FILTER: preserve | fix
RA_WRAP: preserve-and-disclose | fix
COMMON_AREA: skip | enable-and-rerun
UNKNOWN_TIME: local-correction-only | ask-MPC | MPC-approved-correction
MPC_STATUS_SOURCE: ...
C2025Y1_POLICY: keep-58-and-label | move-6-to-known-recovery | other=...
MPC_ROW_RECORDS_58_PLUS_9: supplied path/source=... | unavailable
UPSTREAM_OWNER: ...
UPSTREAM_VERSION/DOC: ...
WEATHER/FAULT_LOGS: available path/source | unavailable
IMAGE_INJECTION: supported | unavailable
```

Any omitted item remains unresolved and must be carried as a limitation rather
than silently assumed.
