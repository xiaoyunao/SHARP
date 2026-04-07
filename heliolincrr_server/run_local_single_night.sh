#!/usr/bin/env bash
set -euo pipefail

# =========================================
# Usage:
#   ./run_local_rr_orbit_one_night.sh YYYYMMDD
# or:
#   ./run_local_rr_orbit_one_night.sh YYYYMMDD /path/to/tracklets_ALL.fits
# Example:
#   ./run_local_rr_orbit_one_night.sh 20260115
# =========================================

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 YYYYMMDD [TRACKLETS_ALL_FITS]"
  exit 1
fi

NIGHT="$1"
CONDA_HOME_DEFAULT="${CONDA_HOME_DEFAULT:-/home/smtpipeline/Softwares/miniconda3}"
HELIOLINC_ENV_NAME="${HELIOLINC_ENV_NAME:-heliolinc}"
HELIOLINC_PYTHON_DEFAULT="${CONDA_HOME_DEFAULT}/envs/${HELIOLINC_ENV_NAME}/bin/python"
PYTHON_BIN="${PYTHON_BIN:-${HELIOLINC_PYTHON_DEFAULT}}"

ROOT_DATA="${HOME}/Desktop"
CODE_DIR="${HOME}/Desktop/smt_asteroid/heliolincrr"

# 默认 tracklets ALL 路径（你按之前流程一般就是这个）
DEFAULT_ALL="${ROOT_DATA}/${NIGHT}/L2/tracklets_linreproj/tracklets_${NIGHT}_ALL.fits"
TRACKLETS_ALL="${2:-$DEFAULT_ALL}"

RR_CORES=4

# RR 输出目录
RR_NIGHT_DIR="${ROOT_DATA}/${NIGHT}/L2/_rr_links/night/${NIGHT}"
ORBIT_OUTDIR="${RR_NIGHT_DIR}/orbit_confirm"

# ---- RR params（单夜常用）----
RR_REF_MODE="mid"
RR_REF_DT=0.05
RR_TOL=0.02
RR_MIN_LEN_OBS=3
RR_MIN_NIGHTS=1
RR_KCAP=200
RR_MAXV=200
RR_MININIT=0.01

# =========================================
# Sanity checks
# =========================================
if [[ ! -d "${CODE_DIR}" ]]; then
  echo "[fatal] CODE_DIR not found: ${CODE_DIR}"
  exit 1
fi

if [[ ! -f "${TRACKLETS_ALL}" ]]; then
  echo "[fatal] tracklets ALL not found: ${TRACKLETS_ALL}"
  echo "        (你可以作为第二个参数显式传入)"
  exit 1
fi

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "[fatal] Python executable not found: ${PYTHON_BIN}"
  echo "        Set PYTHON_BIN explicitly or create conda env '${HELIOLINC_ENV_NAME}'."
  exit 1
fi

cd "${CODE_DIR}"

echo "[info] NIGHT        = ${NIGHT}"
echo "[info] TRACKLETS_ALL= ${TRACKLETS_ALL}"
echo "[info] RR_NIGHT_DIR = ${RR_NIGHT_DIR}"
echo "[info] PYTHON_BIN   = ${PYTHON_BIN}"

# =========================================
# Step 1: RR clustering (single night)
# =========================================
mkdir -p "${RR_NIGHT_DIR}"

"${PYTHON_BIN}" run_rr_from_tracklets.py \
  --infile "${TRACKLETS_ALL}" \
  --outdir "${RR_NIGHT_DIR}" \
  --cores "${RR_CORES}" \
  --ref-epoch-mode "${RR_REF_MODE}" \
  --ref-dt-days "${RR_REF_DT}" \
  --tol "${RR_TOL}" \
  --min-len-obs "${RR_MIN_LEN_OBS}" \
  --min-nights "${RR_MIN_NIGHTS}" \
  --k-neighbors-cap "${RR_KCAP}" \
  --max-v-kms "${RR_MAXV}" \
  --min-init-earth-au "${RR_MININIT}"

# =========================================
# Step 2: orbit confirmation / fitting
# =========================================
mkdir -p "${ORBIT_OUTDIR}"

"${PYTHON_BIN}" orbit_confirm_links.py \
  --rr-dir "${RR_NIGHT_DIR}" \
  --root "${ROOT_DATA}" \
  --outdir "${ORBIT_OUTDIR}" \
  --cores "${RR_CORES}"

echo ""
echo "[done] RR dir   : ${RR_NIGHT_DIR}"
echo "[done] Orbit out: ${ORBIT_OUTDIR}"
