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
TRK_SUB_MAP="${TRK_SUB_MAP:-}"
ASSIGN_UNKNOWN_TRKSUB="${ASSIGN_UNKNOWN_TRKSUB:-1}"
EXPORT_UNKNOWN_ADES="${EXPORT_UNKNOWN_ADES:-0}"
UNKNOWN_REVIEW_CSV="${UNKNOWN_REVIEW_CSV:-}"
REQUIRE_UNKNOWN_REVIEW="${REQUIRE_UNKNOWN_REVIEW:-0}"
VALIDATE_UNKNOWN_MPC="${VALIDATE_UNKNOWN_MPC:-0}"
SUBMIT_UNKNOWN_MPC="${SUBMIT_UNKNOWN_MPC:-0}"
MAX_TRACKLETS_PER_GROUP="${MAX_TRACKLETS_PER_GROUP:-100000}"
MAX_UNKNOWN_LINKS_AFTER_KNOWN="${MAX_UNKNOWN_LINKS_AFTER_KNOWN:-200}"
REQUIRE_KNOWN_MASK15="${REQUIRE_KNOWN_MASK15:-1}"
UNKNOWN_LINK_SKIP_CODE=20

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
TRK_SUB_HISTORY="${TRK_SUB_HISTORY:-${ROOT_OUT}/trksub_history.jsonl}"

MASK_GAIA_DIR="${ROOT_OUT}/${NIGHT}/mask_gaia"
TRACKLET_DIR="${ROOT_OUT}/${NIGHT}/tracklets_linreproj"
NIGHTLY_ALL="${TRACKLET_DIR}/tracklets_${NIGHT}_ALL.fits"

RR_NIGHT_DIR="${ROOT_OUT}/${NIGHT}/rr_links"
ANALYSIS_DIR="${ROOT_OUT}/${NIGHT}/analysis"
PLOTS_ROOT="/pipeline/xiaoyunao/heliolincrr/plots"
SUMMARY_JSON="${ANALYSIS_DIR}/${NIGHT}_single_night_summary.json"
UNKNOWN_JSON="${ROOT_RAW}/${NIGHT}/L4/${NIGHT}_unknown_links.json"
UNKNOWN_FITS="${ROOT_RAW}/${NIGHT}/L4/${NIGHT}_unknown_links.fits"

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
    --match-arcsec 1.5 \
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
    --max-tracklets-per-group "${MAX_TRACKLETS_PER_GROUP}" \
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
    --write-index \
    --allow-empty
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
  echo "[info] nightly ALL has 0 tracklets; continue to write empty link/orbit/unknown outputs: ${NIGHTLY_ALL}"
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
SUMMARY_ARGS=(
  summarize_single_night.py "${NIGHT}"
  --processed-root "${ROOT_RAW}" \
  --root-out "${ROOT_OUT}" \
  --rr-subdir "rr_links" \
  --summary-json "${SUMMARY_JSON}" \
  --unknown-json "${UNKNOWN_JSON}" \
  --unknown-fits "${UNKNOWN_FITS}"
)
if [[ -n "${TRK_SUB_MAP}" ]]; then
  SUMMARY_ARGS+=(--trk-sub-map "${TRK_SUB_MAP}")
fi
if [[ "${REQUIRE_KNOWN_MASK15}" -eq 1 ]]; then
  SUMMARY_ARGS+=(--require-mask15)
fi
"${PYTHON_BIN}" "${SUMMARY_ARGS[@]}"

UNKNOWN_LINK_COUNT="$("${PYTHON_BIN}" - <<PY
import json
from pathlib import Path
path = Path("${UNKNOWN_JSON}")
print(len(json.loads(path.read_text(encoding="utf-8"))) if path.exists() else 0)
PY
)"
echo "[info] unknown links after known subtraction=${UNKNOWN_LINK_COUNT} threshold=${MAX_UNKNOWN_LINKS_AFTER_KNOWN}"
if [[ "${MAX_UNKNOWN_LINKS_AFTER_KNOWN}" -gt 0 && "${UNKNOWN_LINK_COUNT}" -gt "${MAX_UNKNOWN_LINKS_AFTER_KNOWN}" ]]; then
  echo "[skip] unknown links ${UNKNOWN_LINK_COUNT} exceeds MAX_UNKNOWN_LINKS_AFTER_KNOWN=${MAX_UNKNOWN_LINKS_AFTER_KNOWN}; skip ${NIGHT}"
  rm -f "${UNKNOWN_JSON}" \
        "${UNKNOWN_FITS}" \
        "${ROOT_RAW}/${NIGHT}/L4/${NIGHT}_unknown_links_ades.psv" \
        "${ROOT_RAW}/${NIGHT}/L4/${NIGHT}_unknown_mpc_reply.txt" \
        "${ROOT_RAW}/${NIGHT}/L4/${NIGHT}_unknown_mpc_validate_reply.txt" \
        "${ROOT_RAW}/${NIGHT}/L4/${NIGHT}_unknown_validate_reply.txt"
  exit "${UNKNOWN_LINK_SKIP_CODE}"
fi

# =============================
# Step 6a: assign globally unique unknown trkSub IDs
# =============================
if [[ "${ASSIGN_UNKNOWN_TRKSUB}" -eq 1 ]]; then
  "${PYTHON_BIN}" assign_unknown_trksub.py "${NIGHT}" \
    --mode "single-night" \
    --catalog "${UNKNOWN_JSON}" \
    --fits "${UNKNOWN_FITS}" \
    --summary-json "${SUMMARY_JSON}" \
    --history "${TRK_SUB_HISTORY}"
else
  echo "[info] skip unknown trkSub assignment for ${NIGHT} because ASSIGN_UNKNOWN_TRKSUB=0"
fi

# =============================
# Step 6b: optional unknown ADES export / submission
# =============================
if [[ "${EXPORT_UNKNOWN_ADES}" -eq 1 ]]; then
  UNKNOWN_ADES_OUT="${ROOT_RAW}/${NIGHT}/L4/${NIGHT}_unknown_links_ades.psv"
  UNKNOWN_ADES_REPLY="${ROOT_RAW}/${NIGHT}/L4/${NIGHT}_unknown_mpc_reply.txt"
  ADES_ARGS=(
    export_unknown_ades.py "${NIGHT}"
    --processed-root "${ROOT_RAW}"
    --catalog "${UNKNOWN_JSON}"
    --out "${UNKNOWN_ADES_OUT}"
    --response-out "${UNKNOWN_ADES_REPLY}"
  )
  if [[ -n "${UNKNOWN_REVIEW_CSV}" ]]; then
    ADES_ARGS+=(--review-csv "${UNKNOWN_REVIEW_CSV}")
  fi
  if [[ "${REQUIRE_UNKNOWN_REVIEW}" -eq 1 ]]; then
    ADES_ARGS+=(--require-review)
  fi
  if [[ "${VALIDATE_UNKNOWN_MPC}" -eq 1 ]]; then
    ADES_ARGS+=(--validate)
  fi
  if [[ "${SUBMIT_UNKNOWN_MPC}" -eq 1 ]]; then
    ADES_ARGS+=(--submit)
  fi
  "${PYTHON_BIN}" "${ADES_ARGS[@]}"
else
  echo "[info] skip unknown ADES export for ${NIGHT} because EXPORT_UNKNOWN_ADES=0"
fi

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
