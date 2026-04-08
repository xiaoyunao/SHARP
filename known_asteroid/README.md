# Known Asteroid Pipeline v2

This directory contains the server-oriented known-object pipeline for
Observatory `327`.

## Scope

- one process handles one `L2/*MP*` file
- concurrency is controlled by `slurm`
- per-night products are written to `/processed1/YYYYMMDD/L4`
- the pipeline can stop after extraction or continue to MPC submission

The current ADES defaults have already passed MPC test ingest and a live known-
object submission workflow.

## Files

- `match_single_night.py`: single-file or single-night known-object matching
- `merge_night_parts.py`: merge per-file FITS parts into night-level FITS
- `export_ades.py`: FITS-to-ADES exporter plus optional MPC submit
- `run_daily.sh`: unattended daily trigger for the previous night
- `plot_known_asteroids.py`: nightly known-object visualization
- `run_visual_daily.sh`: unattended daily trigger for the visualization step
- `update_all_matched_history.py`: accumulate matched-object history products
- `slurm_match_one_file.sh`: one `*MP*` file per slurm task
- `slurm_merge_submit.sh`: merge night-level outputs and optionally submit
- `submit_pipeline_slurm.sh`: submit slurm jobs for one night or a date range
- `cron.example`: daily 09:00 cron entry
- `cron_visual.example`: daily visualization cron entry

## Default server layout

- code directory: `/pipeline/xiaoyunao/known_asteroid`
- ephemeris / catalog files kept in that directory:
  - `de432s.bsp`
  - `astorb.dat`
- data root: `/processed1`
- nightly input: `/processed1/YYYYMMDD/L2`
- nightly output: `/processed1/YYYYMMDD/L4`

## Default submission configuration

- `submitter name`: `Y.-A. Xiao`
- `acknowledgement email (ac2)`: `wsgp2024@163.com`
- `obs code`: `327`
- `obs name`: `Xinglong Station`
- `minimum detections per object per night`: `3`
- `mode`: `CCD`
- `astCat`: `Gaia3E`
- `photCat`: `Gaia3E`
- `band`: `G`
- `err-unit`: `deg`

## Default header layout

The exporter writes:

- `# observatory`
- `# submitter`
- `# observers`
- `# measurers`
- `# telescope`
- `# coinvestigators`
- `# collaborators`
- `# fundingSource`

The telescope block defaults to:

- `name = 60/90cm Schmidt telescope`
- `design = Schmidt`
- `aperture = 0.6`
- `detector = CCD`
- `fRatio = 3`
- `filter = unfiltered`
- `arraySize = 9216 x 9232`
- `pixelScale = 1.15`

The funding source is intentionally truncated to the first sentence only.

## Skip rules

- if `L4/YYYYMMDD_all_asteroids.fits` and
  `L4/YYYYMMDD_matched_asteroids.fits` already exist, extraction is skipped
- if `L4/YYYYMMDD_matched_asteroids_ades.psv` and
  `L4/YYYYMMDD_mpc_reply.txt` already exist, the whole reporting stage is
  skipped
- file-level reruns also skip any `*MP*` file whose part FITS already exist in
  `L4/known_asteroid_parts`

## Slurm runner

Single night, extraction only:

```bash
cd /pipeline/xiaoyunao/known_asteroid
./submit_pipeline_slurm.sh --batch false --submit-mpc false --max-parallel 24 20260318
```

Date range, extraction only:

```bash
cd /pipeline/xiaoyunao/known_asteroid
./submit_pipeline_slurm.sh --batch true --submit-mpc false --max-parallel 24 20260318 20260320
```

Single night with MPC submission:

```bash
cd /pipeline/xiaoyunao/known_asteroid
./submit_pipeline_slurm.sh --batch false --submit-mpc true --max-parallel 24 20260318
```

The slurm flow is:

1. create a manifest of unprocessed `L2/*MP*` files for the night
2. submit an array job with one task per file
3. submit a dependent finalize job to merge night-level FITS
4. optionally export ADES and submit to MPC

`submit-mpc` is still serial even in the slurm flow. It runs as a single
dependent finalize job after extraction is complete, so it does not create
extra parallel pressure on the server.

## Resource control

Recommended defaults on the current server:

- `--max-parallel 24`
- `MATCH_CPUS=1`
- `MATCH_MEM=10G`
- `FINALIZE_CPUS=1`
- `FINALIZE_MEM=10G`

You can override the resource settings at submit time, for example:

```bash
cd /pipeline/xiaoyunao/known_asteroid
MATCH_MEM=8G FINALIZE_MEM=8G ./submit_pipeline_slurm.sh --batch false --submit-mpc false --max-parallel 16 20260318
```

## Unattended operation

Daily trigger at 09:00:

```bash
cd /pipeline/xiaoyunao/known_asteroid
./run_daily.sh
```

This script checks the previous night, verifies that `L2` exists and contains
`*MP*` files, and only then submits the night to slurm.

The current server tree no longer keeps the older sequential helpers
`run_pipeline.sh` / `run_backfill.sh`; the repository now follows that same
server-only layout.

## Expected outputs

For a night `YYYYMMDD`, the pipeline writes under `/processed1/YYYYMMDD/L4`:

- `known_asteroid_parts/*_all_asteroids.fits`
- `known_asteroid_parts/*_matched_asteroids.fits`
- `YYYYMMDD_all_asteroids.fits`
- `YYYYMMDD_matched_asteroids.fits`
- `YYYYMMDD_matched_asteroids_ades.psv`
- `YYYYMMDD_mpc_reply.txt`
- `YYYYMMDD_match.log`
- `YYYYMMDD_file_manifest.txt`

## Notes

- `slurm` handles parallelism; the matcher itself runs with `--njobs 1`
- the pipeline only touches `L2` and `L4`; other nightly subdirectories are
  ignored
- for real MPC submission tracking, keep the returned submission ID and follow
  WAMO / ingest status on the MPC side
