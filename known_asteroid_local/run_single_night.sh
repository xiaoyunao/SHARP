#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 YYYYMMDD" >&2
  exit 1
fi

NIGHT="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/local_data/processed1}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/local_data/known_asteroid/products}"
ASTORB="${ASTORB:-${PROJECT_ROOT}/resources/known_asteroid/astorb.dat}"
PYTHON="${PYTHON:-python3}"

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
NJOBS="${NJOBS:-5}"
CONFIDENCE_PAD_DEG="${CONFIDENCE_PAD_DEG:-0.5}"
CORNER_PAD_DEG="${CORNER_PAD_DEG:-0.05}"
LON="${LON:-117.575}"
LAT="${LAT:-40.393}"
HEIGHT="${HEIGHT:-960.0}"
VERBOSE="${VERBOSE:-0}"

OBS_CODE="${OBS_CODE:-327}"
OBS_NAME="${OBS_NAME:-Xinglong Station}"
SUBMITTER_NAME="${SUBMITTER_NAME:-Y.-A. Xiao}"
INSTITUTION="${INSTITUTION:-NAOC}"
OBSERVERS="${OBSERVERS:-Xiangnan Guan;Pengfei Liu;Xiaoming Teng;WSGP Team}"
MEASURERS="${MEASURERS:-Niu Li;Yun-Ao Xiao;Jingyi Zhang;Hu Zou;WSGP Team}"
COMMENT="${COMMENT:-}"
TELESCOPE_NAME="${TELESCOPE_NAME:-60/90cm Schmidt telescope}"
APERTURE="${APERTURE:-0.6}"
DESIGN="${DESIGN:-Schmidt}"
DETECTOR="${DETECTOR:-CCD}"
F_RATIO="${F_RATIO:-3}"
FILTER_NAME="${FILTER_NAME:-unfiltered}"
ARRAY_SIZE="${ARRAY_SIZE:-9216 x 9232}"
PIXEL_SCALE="${PIXEL_SCALE:-1.15}"
COINVESTIGATORS="${COINVESTIGATORS:-Hu Zou}"
COLLABORATORS="${COLLABORATORS:-Yun-Ao Xiao;Niu Li;Jingyi Zhang;Wenxiong Li;Zhaobin Chen;Shufei Liu;WSGP Team}"
FUNDING_SOURCE="${FUNDING_SOURCE:-The 60/90-cm Schmidt telescope situated at the Xinglong station is overseen by the research team of the Wide-field Survey and Galaxy Physics (WSGP) at NAOC.}"
ASTCAT="${ASTCAT:-Gaia3E}"
BAND="${BAND:-G}"
PHOTCAT="${PHOTCAT:-Gaia3E}"
ERR_UNIT="${ERR_UNIT:-deg}"
MIN_OBSERVATIONS="${MIN_OBSERVATIONS:-3}"
PROG="${PROG:-}"
INCLUDE_LOGSNR="${INCLUDE_LOGSNR:-0}"

VALIDATE="${VALIDATE:-0}"
SUBMIT="${SUBMIT:-0}"
AC2_EMAIL="${AC2_EMAIL:-wsgp2024@163.com}"
ACK="${ACK:-Known-object ADES submission ${NIGHT}}"
OBJ_TYPE="${OBJ_TYPE:-MBA}"

mkdir -p "${OUTDIR}"

MATCH_LOG="${OUTDIR}/${NIGHT}_match.log"
MATCHED_FITS="${OUTDIR}/${NIGHT}_matched_asteroids.fits"
OUT_PSV="${OUTDIR}/${NIGHT}_known_asteroids_ades.psv"
REPLY_TXT="${OUTDIR}/${NIGHT}_mpc_reply.txt"

MATCH_CMD=(
  "${PYTHON}" "${SCRIPT_DIR}/match_single_night.py" "${NIGHT}"
  --root "${ROOT_DIR}"
  --outdir "${OUTDIR}"
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
  --njobs "${NJOBS}"
  --confidence-pad-deg "${CONFIDENCE_PAD_DEG}"
  --corner-pad-deg "${CORNER_PAD_DEG}"
  --lon "${LON}"
  --lat "${LAT}"
  --height "${HEIGHT}"
  --log "${MATCH_LOG}"
)

if [[ -n "${EXPTIME_KEY}" ]]; then
  MATCH_CMD+=(--exptime-key "${EXPTIME_KEY}")
fi
if [[ "${VERBOSE}" == "1" ]]; then
  MATCH_CMD+=(--verbose)
fi

"${MATCH_CMD[@]}"

if [[ ! -s "${MATCHED_FITS}" ]]; then
  echo "[FATAL] matched FITS not found or empty: ${MATCHED_FITS}" >&2
  exit 2
fi

EXPORT_CMD=(
  "${PYTHON}" "${SCRIPT_DIR}/export_ades.py"
  --fits "${MATCHED_FITS}"
  --out "${OUT_PSV}"
  --obs-code "${OBS_CODE}"
  --obs-name "${OBS_NAME}"
  --submitter-name "${SUBMITTER_NAME}"
  --institution "${INSTITUTION}"
  --observers "${OBSERVERS}"
  --measurers "${MEASURERS}"
  --comment "${COMMENT}"
  --telescope-name "${TELESCOPE_NAME}"
  --aperture "${APERTURE}"
  --design "${DESIGN}"
  --detector "${DETECTOR}"
  --coinvestigators "${COINVESTIGATORS}"
  --collaborators "${COLLABORATORS}"
  --funding-source "${FUNDING_SOURCE}"
  --astcat "${ASTCAT}"
  --band "${BAND}"
  --photcat "${PHOTCAT}"
  --err-unit "${ERR_UNIT}"
  --min-observations "${MIN_OBSERVATIONS}"
  --ac2-email "${AC2_EMAIL}"
  --ack "${ACK}"
  --obj-type "${OBJ_TYPE}"
  --response-out "${REPLY_TXT}"
)

if [[ -n "${F_RATIO}" ]]; then
  EXPORT_CMD+=(--f-ratio "${F_RATIO}")
fi
if [[ -n "${FILTER_NAME}" ]]; then
  EXPORT_CMD+=(--filter-name "${FILTER_NAME}")
fi
if [[ -n "${ARRAY_SIZE}" ]]; then
  EXPORT_CMD+=(--array-size "${ARRAY_SIZE}")
fi
if [[ -n "${PIXEL_SCALE}" ]]; then
  EXPORT_CMD+=(--pixel-scale "${PIXEL_SCALE}")
fi
if [[ -n "${PROG}" ]]; then
  EXPORT_CMD+=(--prog "${PROG}")
fi
if [[ "${INCLUDE_LOGSNR}" == "1" ]]; then
  EXPORT_CMD+=(--include-logsnr)
fi
if [[ "${VALIDATE}" == "1" ]]; then
  EXPORT_CMD+=(--validate)
fi
if [[ "${SUBMIT}" == "1" ]]; then
  EXPORT_CMD+=(--submit)
fi

"${EXPORT_CMD[@]}"

echo "[DONE] night=${NIGHT}"
echo "[DONE] matched_fits=${MATCHED_FITS}"
echo "[DONE] ades_psv=${OUT_PSV}"
