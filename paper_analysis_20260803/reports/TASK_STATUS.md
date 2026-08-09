# P0/P1 task status

Status as of 2026-08-09. This tracks the paper-analysis project only; no
manuscript or production-pipeline file is changed by these tasks.

Legend:

- `[x]` complete from frozen/local evidence;
- `[~]` local analysis complete but the scientific interpretation or final
  scope still needs external input;
- `[ ]` genuine author/observatory/upstream/MPC input still required.

## P0

- [x] Freeze raw/L2/known/unknown manifests, product ledgers, hashes,
  `night_status.csv`, and a machine-readable provisional quality mask.
- [ ] Obtain author sign-off on the quality mask and primary exclusion reason
  for each non-standard night.
- [x] Reconcile the production code/configuration and freeze observable
  dependency/environment provenance.
- [x] Complete the 128-night observer-location sensitivity comparison:
  87,850 link keys matched; zero `fit_ok` and `is_good` flips overall, in the
  formal 4,762, and in the retained 58.
- [ ] Obtain the authoritative surveyed MPC 327 position, height type, and
  datum, then choose the historical-as-run versus code-change policy. Repeat
  the sensitivity test only if the surveyed site differs from both tested
  definitions.
- [x] Inventory all upstream L1/L2 metadata that can be established from files
  and code, including the contaminated GMG site headers.
- [ ] Obtain the upstream reduction/mask/depth documentation and owner-approved
  correction for the header template.
- [x] Build the 67-row frozen review ledger, 58-row retained table,
  37-component cross-night linear-motion screen, and two-pass JPL numerical
  reconciliation. The review chain is 68 marked real -> one withdrawal
  (`20260507/00001et`) -> 67 submission-selected -> 58 post-audit retained.
- [ ] Obtain row-level MPC records for the 67 submission-selected rows and
  approve the handling of the six C/2025 Y1-associated links.
- [x] Build the local 179-observation unknown-time midpoint-correction analysis
  with exact L2 joins and no production/MPC mutation.
- [ ] Choose whether the time correction remains analysis-only or requires an
  MPC-approved communication/correction workflow.

## P1

- [x] Complete nominal-footprint known-object recovery, signed astrometric
  residuals, detector/night/magnitude/rate trends, binned intervals, and the
  deterministic random-shift chance-match control.
- [~] A formally “detectable” recovery denominator still needs upstream depth,
  mask, ephemeris-uncertainty, and quality definitions; current outputs are
  explicitly nominal-footprint results.
- [x] Complete the full unknown-stage funnel with nightly distributions,
  detection/tracklet/linkage units, both orbit gates, quarantine, true zero
  versus not-run, and the 68 -> 67 -> 58 human gates.
- [~] Complete the retrospective blinded-known linkage/orbit survival proxy.
  Controlled catalog injection still needs an author-approved grid; image-level
  injection additionally needs an upstream execution interface and dry-run
  route.
- [~] Complete the available machine-evidence latency analysis and Figure 12,
  separating automatic, human-review, and MPC intervals and isolating reruns,
  mtime-only cases, negative intervals, and restart anomalies. Authoritative
  human/MPC timestamps and independent CPU/RAM evidence remain external; CPU
  and peak RAM are currently reported as unavailable.
- [~] Complete scheduler plan/cadence/compliance analysis and Figure 4. Weather,
  fault, override, and control-system logs remain necessary before interpreting
  the result as observatory efficiency or full plan completion.

## Tables and figures

- [x] Generate Tables 1--5 as CSV and TeX with a machine-readable summary and
  hash manifest.
- [x] Generate Figures 1--12 as both PDF and PNG with non-empty,
  positive-dimension machine checks.
- [~] Final scientific labeling, representative-example approval, and editorial
  visual acceptance remain author decisions; figure generation itself is not
  pending.

## Current validation state

A fresh run of `scripts/validate_snapshot.py` against the current config,
inventory, frozen products, review sample, derived known/unknown products,
frozen tables, and all figures returned **63 PASS, 0 FAIL, 1 SKIP**. The sole
SKIP is the genuine external item `quality.author_signoff`; all path, summary,
hash, and unit-metadata checks now pass.

The snapshot assembler reports `complete` with 555 hashed artifacts and zero
blocked, incomplete, or failed components. Figure QA reports 12/12 pass, zero
failed/incomplete figures, and a complete contact sheet. All headline count
checks, formal-unknown `fit_ok/is_good` checks, 58-key consistency checks, core
provenance hashes, and all 12 PDF/PNG pairs pass.

## Genuine external blockers only

1. Surveyed site/datum and historical-versus-change policy.
2. Quality-mask and interpretation sign-off.
3. Upstream L1/L2 reduction, depth, effective-mask, and header-correction
   documentation.
4. Row-level MPC state, C/2025 Y1 handling, and unknown-time correction policy.
5. Human/MPC event logs, weather/fault/control logs, and any independent
   CPU/RAM accounting.
6. Approved controlled-injection design and an image-processing dry-run route.
