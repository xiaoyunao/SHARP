#!/bin/sh
set -eu

SHARP_REPO="${SHARP_REPO:-/Users/yunaoxiao/Desktop/smt_asteroid}"
SHARP_OUT="${SHARP_OUT:-$SHARP_REPO/SHARP_local_analysis_figures_20260830}"
SHARP_SNAPSHOT="${SHARP_SNAPSHOT:-$SHARP_REPO/paper_analysis_20260803/snapshot}"
SHARP_ASTORB="${SHARP_ASTORB:-$SHARP_REPO/tmp/SHARP_local_sources_20260830/astorb/astorb.dat}"
SHARP_FIG01_L1="${SHARP_FIG01_L1:-$SHARP_REPO/tmp/SHARP_redraw_sources_20260830/OBJ_MP_0940_0273.fits.gz}"
SHARP_FIG06_GIF="${SHARP_FIG06_GIF:-/Users/yunaoxiao/Desktop/submitted_unknown_20251116_20260628/gifs/20260210_00000P2_link0035.gif}"
SHARP_PY="${SHARP_PY:-/opt/anaconda3/bin/python3}"
SHARP_ORBIT_PY="${SHARP_ORBIT_PY:-$SHARP_REPO/tmp/SHARP_poliastro_env_20260830/bin/python}"

"$SHARP_PY" "$SHARP_OUT/scripts/generate_products.py" --output "$SHARP_OUT" --snapshot "$SHARP_SNAPSHOT" --astorb "$SHARP_ASTORB" --repo "$SHARP_REPO" --gotta /Users/yunaoxiao/Desktop/gotta_asteroid_1
"$SHARP_ORBIT_PY" "$SHARP_OUT/scripts/generate_candidate_families.py" --output "$SHARP_OUT" --measurements "$SHARP_SNAPSHOT/frozen_products/unknown_review_detections.parquet" --repo "$SHARP_REPO"
"$SHARP_PY" "$SHARP_REPO/paper_analysis_20260803/scripts/make_fig01_architecture.py" --output-stem "$SHARP_OUT/figures/fig02_moving_object_workflow"
"$SHARP_PY" "$SHARP_OUT/scripts/export_workflow_data.py" --root "$SHARP_OUT" --architecture-script "$SHARP_REPO/paper_analysis_20260803/scripts/make_fig01_architecture.py"
"$SHARP_PY" "$SHARP_OUT/scripts/redraw_image_figures.py" --root "$SHARP_OUT" --fig01-image "$SHARP_FIG01_L1" --known-matched "$SHARP_SNAPSHOT/frozen_products/known_matched.parquet" --fig06-gif "$SHARP_FIG06_GIF"
"$SHARP_PY" "$SHARP_OUT/scripts/redraw_scientific_figures.py" --root "$SHARP_OUT" --snapshot "$SHARP_SNAPSHOT"
"$SHARP_PY" "$SHARP_OUT/scripts/validate_bundle.py" --root "$SHARP_OUT"
