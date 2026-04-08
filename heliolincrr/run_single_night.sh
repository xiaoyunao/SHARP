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

W15_OUTDIR="${ROOT_OUT}/${NIGHT}/tracklets_windows_15"
W15_NDAYS=15

RR_NIGHT_DIR="${ROOT_OUT}/${NIGHT}/rr_links"
RR_W15_DIR="${ROOT_OUT}/${NIGHT}/rr_links_15"

# =============================
# RR parameters — single night
# =============================
RR_NIGHT_CORES=16
RR_NIGHT_REF_MODE="mid"
RR_NIGHT_REF_DT=0.05
RR_NIGHT_TOL=0.02
RR_NIGHT_MIN_LEN_OBS=3
RR_NIGHT_MIN_NIGHTS=1
RR_NIGHT_KCAP=200
RR_NIGHT_MAXV=200
RR_NIGHT_MININIT=0.01

# =============================
# RR parameters — 15 nights
# =============================
RR_W15_CORES=20
RR_W15_REF_MODE="mid"
RR_W15_REF_DT=0.50
RR_W15_TOL=0.02
RR_W15_MIN_LEN_OBS=4
RR_W15_MIN_NIGHTS=2
RR_W15_KCAP=200
RR_W15_MAXV=200
RR_W15_MININIT=0.02

CACHE_PROP=0
RR_CACHE_FLAG=""
if [[ "${CACHE_PROP}" -eq 1 ]]; then
  RR_CACHE_FLAG="--cache-prop"
fi

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

prepare_previous_nights_for_w15() {
  local end_night="$1"
  local nights
  mapfile -t nights < <("${PYTHON_BIN}" - <<PY
from datetime import datetime, timedelta
end = datetime.strptime("${end_night}", "%Y%m%d")
start = end - timedelta(days=${W15_NDAYS} - 1)
d = start
while d < end:
    print(d.strftime("%Y%m%d"))
    d += timedelta(days=1)
PY
)

  local n
  for n in "${nights[@]}"; do
    if ! is_valid_table "$(nightly_all_path "${n}")"; then
      echo "[info] missing/invalid prior nightly ALL for W15: ${n}"
      ensure_nightly_all "${n}"
    fi
  done
}

# =============================
# Early check: reuse existing nightly ALL
# If exists, skip Step 1–3
# =============================
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
# Step 4a: RR — single night (ALWAYS run if Step3 passed)
# =============================
mkdir -p "${RR_NIGHT_DIR}"

"${PYTHON_BIN}" run_rr_from_tracklets.py \
  --infile "${NIGHTLY_ALL}" \
  --outdir "${RR_NIGHT_DIR}" \
  --cores "${RR_NIGHT_CORES}" \
  --ref-epoch-mode "${RR_NIGHT_REF_MODE}" \
  --ref-dt-days "${RR_NIGHT_REF_DT}" \
  --tol "${RR_NIGHT_TOL}" \
  --min-len-obs "${RR_NIGHT_MIN_LEN_OBS}" \
  --min-nights "${RR_NIGHT_MIN_NIGHTS}" \
  --k-neighbors-cap "${RR_NIGHT_KCAP}" \
  --max-v-kms "${RR_NIGHT_MAXV}" \
  --min-init-earth-au "${RR_NIGHT_MININIT}" \
  ${RR_CACHE_FLAG}

if [[ "${RUN_W15}" -ne 1 ]]; then
  echo "[info] RUN_W15=${RUN_W15}; skip W15 processing"
  exit 0
fi

# =============================
# Step 4b: merge rolling window (15 nights)
#         failure here MUST NOT block single-night RR (B)
# =============================
if [[ "${PREP_W15_MISSING_NIGHTS}" -eq 1 ]]; then
  prepare_previous_nights_for_w15 "${NIGHT}"
fi

mkdir -p "${W15_OUTDIR}"

"${PYTHON_BIN}" merge_tracklets_window15.py "${NIGHT}" \
  --root "${ROOT_OUT}" \
  --ndays "${W15_NDAYS}" \
  --nightly-subpath "tracklets_linreproj" \
  --nightly-name-template "tracklets_{night}_ALL.fits" \
  --outdir "${W15_OUTDIR}" \
  --overwrite \
  --write-index \
  || {
    echo "[skip] W15 RR skipped because window merge failed"
    exit 0
  }

# =============================
# Step 5: find W15 file safely (no ls | tail)
# =============================
shopt -s nullglob
W15_CANDIDATES=(${W15_OUTDIR}/tracklets_*_${NIGHT}_W${W15_NDAYS}.fits)
shopt -u nullglob

if [[ ${#W15_CANDIDATES[@]} -eq 0 ]]; then
  echo "[skip] W15 RR skipped because no window file found"
  exit 0
fi

W15_FILE="${W15_CANDIDATES[-1]}"
echo "[info] W15_FILE=${W15_FILE}"

# =============================
# Step 6: read meta.json, check n_nights_found (B)
# =============================
META_JSON="${W15_FILE}.meta.json"
if [[ ! -f "${META_JSON}" ]]; then
  echo "[skip] W15 RR skipped because meta json missing"
  exit 0
fi

N_FOUND="$("${PYTHON_BIN}" - <<PY
import json
with open("${META_JSON}","r",encoding="utf-8") as f:
    m=json.load(f)
print(int(m.get("n_nights_found",0)))
PY
)"

echo "[info] W15 n_nights_found=${N_FOUND}"

if [[ "${N_FOUND}" -lt 5 ]]; then
  echo "[skip] W15 RR skipped because n_nights_found < 5"
  exit 0
fi

# =============================
# Step 7: RR — 15-night window
# =============================
mkdir -p "${RR_W15_DIR}"

"${PYTHON_BIN}" run_rr_from_tracklets.py \
  --infile "${W15_FILE}" \
  --outdir "${RR_W15_DIR}" \
  --cores "${RR_W15_CORES}" \
  --ref-epoch-mode "${RR_W15_REF_MODE}" \
  --ref-dt-days "${RR_W15_REF_DT}" \
  --tol "${RR_W15_TOL}" \
  --min-len-obs "${RR_W15_MIN_LEN_OBS}" \
  --min-nights "${RR_W15_MIN_NIGHTS}" \
  --k-neighbors-cap "${RR_W15_KCAP}" \
  --max-v-kms "${RR_W15_MAXV}" \
  --min-init-earth-au "${RR_W15_MININIT}" \
  ${RR_CACHE_FLAG}
