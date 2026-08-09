# SHARP PASP paper analysis (2026-08-03 science snapshot)

This directory contains the reproducible data accounting, statistical analyses,
quality checks, tables, and figures requested by the PASP draft's
`AUTHOR_TODO_AND_ANALYSIS_PLAN.md` and `FIGURE_PLAN.md`.

The manuscript itself is deliberately not copied or edited here.

## Scope and source of truth

- Science interval: `2025-11-15` through `2026-07-15` (inclusive).
- Snapshot label: `2026-08-03`; extraction times are recorded separately.
- Production algorithm commit: `d2f0057bd4cdd4c7a3c6c2431f6449a72c327284`.
- Repository commit reviewed: `cfee217611bca4e71b40ea3687a0906a39849997`.
- Server roots: `/raw1`, `/processed1`, and `/pipeline/xiaoyunao`.
- Headline numbers are derived from frozen tables, never hard-coded in plotting
  scripts.

## Layout

- `config/`: frozen analysis configuration and quality-mask rules.
- `scripts/`: collectors, analyses, validation, and figure generators.
- `snapshot/`: frozen manifests, tables, compact source extracts, and hashes.
- `figure_data/`: small, reviewable data tables consumed by figure scripts.
- `figures/`: the twelve requested paper figures in PDF and PNG form.
- `notebooks/`: executed reader-facing analysis notebook.
- `reports/`: task status, source reconciliation, methods, caveats, and QA.
- `qa/`: render previews, contact sheets, and machine-readable validation output.
- `logs/`: reproducibility metadata and command records.

Large FITS extracts, image cutouts, and transient caches remain local to this
directory and are excluded from Git by the nested `.gitignore`.

## Reproduction

Run scripts with the repository's documented Python environments.  The remote
environment, production-file digests, and Gaia-tile inventory are frozen under
`snapshot/provenance/`. Remote collection is read-only with respect to
production data and stages temporary products under a uniquely named `/tmp`
directory.

The final audit entry points are:

```bash
python3 scripts/validate_snapshot.py --tables snapshot/tables
python3 scripts/build_analysis_notebook.py
python3 scripts/make_figure_contact_sheet.py
python3 scripts/assemble_paper_snapshot.py --overwrite-generated
shasum -a 256 -c hashes.sha256
```

Generators refuse to overwrite most existing outputs. Move an old generated
artifact aside before deliberately rebuilding it. The assembled local snapshot
contains 555 hashed artifacts and currently reports `complete`; the independent
validator reports 63 PASS, 0 FAIL, and one intentional SKIP for the still
unsigned author quality-mask decision.

## Evidence boundaries

- Scheduler and known-object code use `117.575 deg E, 40.393 deg N, 960 m`.
  The frozen orbit-confirmation production run used `117.575 deg E,
  40.394239 deg N, 868.221 m`. A full 128-night alternate-location rerun caused
  zero `fit_ok` or `is_good` flips, including the 4,762 formal links and 58
  retained links, but short-arc orbital elements remain unstable and are not a
  population-inference product.
- The 58 retained rows are single-night linkages, not confirmed independent
  objects or discoveries. Linear-motion screening yields 37 candidate
  components; six linkage rows have a repeatable numerical association with
  C/2025 Y1 (ATLAS), pending authoritative MPC reconciliation.
- Human-review states are explicitly separated: 68 marked real, 67 entered the
  frozen submission set, and 58 survived the post-hoc audit.
- The known-object denominator is a nominal WCS-rectangle prediction set, not a
  true photometrically detectable completeness denominator. The blinded-known
  result is a conditional downstream survival proxy, not image-level detection
  completeness.
- Operations percentiles use only two normal-daily nights with explicit event
  chains. CPU/RAM, weather, and device-efficiency claims are unavailable and
  are never zero-filled.

See `reports/MANUSCRIPT_CODE_CONFLICTS.md` for the full code-versus-draft audit
and `reports/AUTHOR_INPUTS_REQUIRED.md` for decisions that cannot be inferred
from code or archived products.
