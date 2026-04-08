#!/usr/bin/env bash
#SBATCH -J ka_match
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 YYYYMMDD MANIFEST" >&2
  exit 1
fi

NIGHT="$1"
MANIFEST="$2"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

ROOT_DIR="${ROOT_DIR:-/processed1}"
SCRIPT_DIR="${SCRIPT_DIR:-/pipeline/xiaoyunao/known_asteroid}"
ASTORB="${ASTORB:-${SCRIPT_DIR}/astorb.dat}"
PYTHON="${PYTHON:-/home/smtpipeline/Softwares/miniconda3/bin/python}"

FILE_REGEX="${FILE_REGEX:-.*MP.*\\.(fits|fits\\.gz)$}"
HDU="${HDU:-1}"
NX="${NX:-9216}"
NY="${NY:-9232}"
OBS_DATE_KEY="${OBS_DATE_KEY:-OBS_DATE}"
EXPTIME_KEY="${EXPTIME_KEY:-}"
EXPTIME_SEC="${EXPTIME_SEC:-30.0}"
SEP_ARCSEC="${SEP_ARCSEC:-1.0}"
MAG_LIMIT="${MAG_LIMIT:-22.5}"
MAGDIFF="${MAGDIFF:-99.0}"
CAT_MAG_COL="${CAT_MAG_COL:-Mag_Kron}"
CONFIDENCE_PAD_DEG="${CONFIDENCE_PAD_DEG:-0.5}"
CORNER_PAD_DEG="${CORNER_PAD_DEG:-0.05}"
LON="${LON:-117.575}"
LAT="${LAT:-40.393}"
HEIGHT="${HEIGHT:-960.0}"
VERBOSE="${VERBOSE:-0}"

L4_DIR="${ROOT_DIR}/${NIGHT}/L4"
PARTS_DIR="${L4_DIR}/known_asteroid_parts"
MATCH_LOG="${L4_DIR}/${NIGHT}_match.log"
mkdir -p "${L4_DIR}" "${PARTS_DIR}"

LINE_NO=$((TASK_ID + 1))
FILE="$(sed -n "${LINE_NO}p" "${MANIFEST}")"
if [[ -z "${FILE}" ]]; then
  echo "[FATAL] no manifest entry for task ${TASK_ID}" >&2
  exit 2
fi

CMD=(
  "${PYTHON}" "${SCRIPT_DIR}/match_single_night.py" "${NIGHT}"
  --root "${ROOT_DIR}"
  --outdir "${L4_DIR}"
  --parts-dir "${PARTS_DIR}"
  --file "${FILE}"
  --astorb "${ASTORB}"
  --file-regex "${FILE_REGEX}"
  --hdu "${HDU}"
  --nx "${NX}"
  --ny "${NY}"
  --obs-date-key "${OBS_DATE_KEY}"
  --exptime-sec "${EXPTIME_SEC}"
  --sep-arcsec "${SEP_ARCSEC}"
  --mag-limit "${MAG_LIMIT}"
  --magdiff "${MAGDIFF}"
  --cat-mag-col "${CAT_MAG_COL}"
  --njobs 1
  --confidence-pad-deg "${CONFIDENCE_PAD_DEG}"
  --corner-pad-deg "${CORNER_PAD_DEG}"
  --lon "${LON}"
  --lat "${LAT}"
  --height "${HEIGHT}"
  --log "${MATCH_LOG}"
)

if [[ -n "${EXPTIME_KEY}" ]]; then
  CMD+=(--exptime-key "${EXPTIME_KEY}")
fi
if [[ "${VERBOSE}" == "1" ]]; then
  CMD+=(--verbose)
fi

"${CMD[@]}"
