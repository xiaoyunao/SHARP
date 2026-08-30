#!/bin/sh
set -eu

ROOT="/Users/yunaoxiao/Desktop/smt_asteroid"
OUT="$ROOT/SHARP_local_analysis_figures_20260830"
SNAP="$ROOT/paper_analysis_20260803/snapshot"
SRC="$ROOT/tmp/SHARP_local_sources_20260830"
PY="/opt/anaconda3/bin/python3"
ORBIT_PY="$ROOT/tmp/SHARP_poliastro_env_20260830/bin/python"

"$PY" "$OUT/scripts/generate_products.py" --output "$OUT" --snapshot "$SNAP" --astorb "$SRC/astorb/astorb.dat" --repo "$ROOT" --gotta /Users/yunaoxiao/Desktop/gotta_asteroid_1
"$PY" "$OUT/scripts/generate_image_figures.py" --output "$OUT" --fig01-image "$SRC/fig01/OBJ_MP_0940_0273.fits.gz" --known-matched "$SNAP/frozen_products/known_matched.parquet" --fig06-images "$SRC/fig06/OBJ_MP_1602_0203.fits.gz" "$SRC/fig06/OBJ_MP_1602_0229.fits.gz" "$SRC/fig06/OBJ_MP_1602_0255.fits.gz" --fig06-catalogs "$SRC/fig06_l2/OBJ_MP_1602_0203_cat.fits.gz" "$SRC/fig06_l2/OBJ_MP_1602_0229_cat.fits.gz" "$SRC/fig06_l2/OBJ_MP_1602_0255_cat.fits.gz" --measurements "$SNAP/frozen_products/unknown_review_detections.parquet"
"$ORBIT_PY" "$OUT/scripts/generate_candidate_families.py" --output "$OUT" --measurements "$SNAP/frozen_products/unknown_review_detections.parquet" --repo "$ROOT"
"$PY" "$OUT/scripts/validate_bundle.py" --root "$OUT"
