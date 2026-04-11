#!/usr/bin/env bash
set -euo pipefail

# =============================
# Usage
# =============================
# ./run_single_night.sh YYYYMMDD

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 YYYYMMDD"
  exit 1
fi

NIGHT="$1"
CONDA_HOME_DEFAULT="${CONDA_HOME_DEFAULT:-/home/smtpipeline/Softwares/miniconda3}"
HELIOLINC_ENV_NAME="${HELIOLINC_ENV_NAME:-heliolinc}"
HELIOLINC_PYTHON_DEFAULT="${CONDA_HOME_DEFAULT}/envs/${HELIOLINC_ENV_NAME}/bin/python"
PYTHON_BIN="${PYTHON_BIN:-${HELIOLINC_PYTHON_DEFAULT}}"
RUN_W15="${RUN_W15:-0}"
PREP_W15_MISSING_NIGHTS="${PREP_W15_MISSING_NIGHTS:-1}"
FORCE_MASK_GAIA="${FORCE_MASK_GAIA:-0}"
SKIP_PLOTS="${SKIP_PLOTS:-0}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "[fatal] Python executable not found: ${PYTHON_BIN}"
  echo "        Set PYTHON_BIN explicitly or create conda env '${HELIOLINC_ENV_NAME}'."
  exit 1
fi

echo "[info] PYTHON_BIN=${PYTHON_BIN}"

# =============================
# Global paths
# =============================
ROOT_RAW="/processed1"
ROOT_OUT="/pipeline/xiaoyunao/data/heliolincrr"

MASK_GAIA_DIR="${ROOT_OUT}/${NIGHT}/mask_gaia"
TRACKLET_DIR="${ROOT_OUT}/${NIGHT}/tracklets_linreproj"
NIGHTLY_ALL="${TRACKLET_DIR}/tracklets_${NIGHT}_ALL.fits"

RR_NIGHT_DIR="${ROOT_OUT}/${NIGHT}/rr_links"
ANALYSIS_DIR="${ROOT_OUT}/${NIGHT}/analysis"
PLOTS_ROOT="/pipeline/xiaoyunao/heliolincrr/plots"

LINEAR_NIGHT_SPEED=5
LINEAR_NIGHT_DIRECTION=10

night_mask_dir() {
  local night="$1"
  echo "${ROOT_OUT}/${night}/mask_gaia"
}

night_tracklet_dir() {
  local night="$1"
  echo "${ROOT_OUT}/${night}/tracklets_linreproj"
}

nightly_all_path() {
  local night="$1"
  echo "$(night_tracklet_dir "${night}")/tracklets_${night}_ALL.fits"
}

is_valid_table() {
  local fits_path="$1"
  [[ -f "${fits_path}" ]] || return 1
  "${PYTHON_BIN}" - <<PY >/dev/null 2>&1
from astropy.table import Table
Table.read("${fits_path}")
PY
}

count_l2_masked_files() {
  local night="$1"
  local mask_dir
  mask_dir="$(night_mask_dir "${night}")"
  if [[ ! -d "${mask_dir}" ]]; then
    echo 0
    return 0
  fi
  find "${mask_dir}" -maxdepth 1 -name 'OBJ_MP_*_cat.fits.gz' 2>/dev/null | wc -l | tr -d ' '
}

run_mask_gaia_step() {
  local night="$1"
  local mask_dir
  mask_dir="$(night_mask_dir "${night}")"
  mkdir -p "${mask_dir}"

  echo "[run] mask_gaia night=${night}"
  "${PYTHON_BIN}" mask_gaia.py "${night}" \
    --root "${ROOT_RAW}" \
    --gaia-dir "/pipeline/ref/healpix" \
    --out-root "${ROOT_OUT}" \
    --gaia-cone-deg 5.0 \
    --match-arcsec 1.0 \
    --pm-zp 2016.0 \
    --hdu-header 1 \
    --mag-psf-max 21.0 \
    --require-flag0 \
    --nproc 16 \
    --log "${mask_dir}/mask_gaia.log"
}

run_make_tracklets_step() {
  local night="$1"
  local tracklet_dir
  tracklet_dir="$(night_tracklet_dir "${night}")"
  mkdir -p "${tracklet_dir}"

  echo "[run] make_tracklet night=${night}"
  "${PYTHON_BIN}" make_tracklet_linreproj.py "${night}" \
    --root "${ROOT_OUT}" \
    --input-subdir "mask_gaia" \
    --outdir "${tracklet_dir}" \
    --mag-max 0 \
    --vmin 3.0 \
    --vmax 63.0 \
    --dmag-max 1.0 \
    --r-static 2.0 \
    --min-repeat 2 \
    --erode-pix 30 \
    --fallback-margin 0.1 \
    --hdu 1 \
    --nproc 16 \
    --edge-pix 500 \
    --skip-common-area \
    --log "${tracklet_dir}/make_tracklet.log"
}

run_merge_tracklets_step() {
  local night="$1"
  local tracklet_dir nightly_all
  tracklet_dir="$(night_tracklet_dir "${night}")"
  nightly_all="$(nightly_all_path "${night}")"

  echo "[run] merge_tracklets night=${night}"
  "${PYTHON_BIN}" merge_tracklets_night.py "${night}" \
    --root "${ROOT_OUT}" \
    --tracklet-dir "${tracklet_dir}" \
    --out "${nightly_all}" \
    --overwrite \
    --write-index
}

ensure_nightly_all() {
  local night="$1"
  local nightly_all mask_n
  nightly_all="$(nightly_all_path "${night}")"

  if is_valid_table "${nightly_all}"; then
    echo "[info] valid nightly ALL already exists: ${nightly_all}"
    return 0
  fi

  if [[ -f "${nightly_all}" ]]; then
    echo "[warn] invalid nightly ALL detected, rebuilding: ${nightly_all}"
    rm -f "${nightly_all}" "${nightly_all}.inputs.txt"
  fi

  mask_n="$(count_l2_masked_files "${night}")"
  if [[ "${FORCE_MASK_GAIA}" -eq 1 || "${mask_n}" -le 0 ]]; then
    run_mask_gaia_step "${night}"
  else
    echo "[info] reuse existing mask_gaia outputs for ${night} (${mask_n} files)"
  fi

  run_make_tracklets_step "${night}"
  run_merge_tracklets_step "${night}"

  if ! is_valid_table "${nightly_all}"; then
    echo "[fatal] nightly ALL is still invalid after rebuild: ${nightly_all}"
    return 1
  fi
}

if is_valid_table "${NIGHTLY_ALL}"; then
  echo "[info] nightly ALL already exists, skip Step1–3:"
  echo "       ${NIGHTLY_ALL}"
  SKIP_PREP=1
else
  SKIP_PREP=0
fi

if [[ "${SKIP_PREP}" -eq 0 ]]; then
  ensure_nightly_all "${NIGHT}"
else
  echo "[info] Skip Step1–3 (mask / make / merge)"
fi

# ---- hard check: nightly ALL must exist and be non-empty (A)
if [[ ! -f "${NIGHTLY_ALL}" ]]; then
  echo "[fatal] nightly ALL not created: ${NIGHTLY_ALL}"
  exit 1
fi

N_TRK_NIGHT="$("${PYTHON_BIN}" - <<PY
from astropy.table import Table
t = Table.read("${NIGHTLY_ALL}")
print(len(t))
PY
)"

if [[ "${N_TRK_NIGHT}" -le 0 ]]; then
  echo "[fatal] nightly ALL has 0 tracklets: ${NIGHTLY_ALL}"
  exit 1
fi

# =============================
# Step 4: linear linking — single night
# =============================
mkdir -p "${RR_NIGHT_DIR}"
mkdir -p "${ANALYSIS_DIR}"

"${PYTHON_BIN}" run_linear_links_from_tracklets.py \
  --infile "${NIGHTLY_ALL}" \
  --outdir "${RR_NIGHT_DIR}" \
  --speed-thresh-arcsec-per-hour "${LINEAR_NIGHT_SPEED}" \
  --direction-thresh-deg "${LINEAR_NIGHT_DIRECTION}" \
  --require-shared-endpoint

# =============================
# Step 5: orbit confirm
# =============================
rm -rf "${RR_NIGHT_DIR}/orbit_confirm"
"${PYTHON_BIN}" orbit_confirm_links.py \
  --rr-dir "${RR_NIGHT_DIR}" \
  --tracklets "${NIGHTLY_ALL}" \
  --cores 16 \
  --log-every 500

# =============================
# Step 6: summarize
# =============================
"${PYTHON_BIN}" summarize_single_night.py "${NIGHT}" \
  --processed-root "${ROOT_RAW}" \
  --root-out "${ROOT_OUT}" \
  --rr-subdir "rr_links"

# =============================
# Step 7: visualize unknown fit_ok links
# =============================
if [[ "${SKIP_PLOTS}" -eq 1 ]]; then
  echo "[info] skip Step7 plots for ${NIGHT} because SKIP_PLOTS=1"
else
  "${PYTHON_BIN}" plot_unknown_links.py "${NIGHT}" \
    --processed-root "${ROOT_RAW}" \
    --root-out "${ROOT_OUT}" \
    --plot-root "${PLOTS_ROOT}"
fi

echo "[done] single-night pipeline finished for ${NIGHT}"
