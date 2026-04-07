#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_pipeline.sh --batch false --submit-mpc false YYYYMMDD
  ./run_pipeline.sh --batch true --submit-mpc false START_YYYYMMDD END_YYYYMMDD

Notes:
  - Default ROOT_DIR is /processed1
  - Outputs are written under ROOT_DIR/YYYYMMDD/L4
  - Parallelism is file-level only: one process per *MP* file
EOF
}

if [[ $# -lt 3 ]]; then
  usage >&2
  exit 1
fi

BATCH="false"
SUBMIT_MPC="false"
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch)
      BATCH="$2"
      shift 2
      ;;
    --submit-mpc)
      SUBMIT_MPC="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

set -- "${POSITIONAL[@]}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="${ROOT_DIR:-${PROJECT_ROOT}/local_data/processed1}"
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
AC2_EMAIL="${AC2_EMAIL:-wsgp2024@163.com}"
ACK_PREFIX="${ACK_PREFIX:-Known-object ADES submission}"
OBJ_TYPE="${OBJ_TYPE:-MBA}"

if [[ "${BATCH}" == "true" ]]; then
  if [[ $# -ne 2 ]]; then
    usage >&2
    exit 1
  fi
  START_NIGHT="$1"
  END_NIGHT="$2"
else
  if [[ $# -ne 1 ]]; then
    usage >&2
    exit 1
  fi
  START_NIGHT="$1"
  END_NIGHT="$1"
fi

if [[ ! "${START_NIGHT}" =~ ^[0-9]{8}$ || ! "${END_NIGHT}" =~ ^[0-9]{8}$ ]]; then
  echo "[FATAL] nights must be YYYYMMDD" >&2
  exit 2
fi

NIGHTS=()
while IFS= read -r NIGHT_NAME; do
  NIGHTS+=("${NIGHT_NAME}")
done < <(
  find "${ROOT_DIR}" -maxdepth 1 -mindepth 1 -type d -exec basename {} \; |
  grep -E '^[0-9]{8}$' |
  awk -v s="${START_NIGHT}" -v e="${END_NIGHT}" '$0 >= s && $0 <= e' |
  sort
)

if [[ "${#NIGHTS[@]}" -eq 0 ]]; then
  echo "[INFO] No nights found in range [${START_NIGHT}, ${END_NIGHT}] under ${ROOT_DIR}"
  exit 0
fi

for NIGHT in "${NIGHTS[@]}"; do
  NIGHT_DIR="${ROOT_DIR}/${NIGHT}"
  L2_DIR="${NIGHT_DIR}/L2"
  L4_DIR="${NIGHT_DIR}/L4"
  PARTS_DIR="${L4_DIR}/known_asteroid_parts"
  MATCH_LOG="${L4_DIR}/${NIGHT}_match.log"
  ALL_FITS="${L4_DIR}/${NIGHT}_all_asteroids.fits"
  MATCHED_FITS="${L4_DIR}/${NIGHT}_matched_asteroids.fits"
  OUT_PSV="${L4_DIR}/${NIGHT}_matched_asteroids_ades.psv"
  REPLY_TXT="${L4_DIR}/${NIGHT}_mpc_reply.txt"

  mkdir -p "${L4_DIR}" "${PARTS_DIR}"

  if [[ "${SUBMIT_MPC}" == "true" && -s "${OUT_PSV}" && -s "${REPLY_TXT}" ]]; then
    echo "[SKIP] ${NIGHT} (report products already exist in L4)"
    continue
  fi

  if [[ ! -d "${L2_DIR}" ]]; then
    echo "[SKIP] ${NIGHT} (missing L2 dir)"
    continue
  fi

  if [[ -s "${ALL_FITS}" && -s "${MATCHED_FITS}" ]]; then
    echo "[SKIP] ${NIGHT} extraction (night-level FITS already exist)"
  else
    FILES=()
    while IFS= read -r FILE_NAME; do
      FILES+=("${FILE_NAME}")
    done < <(
      find "${L2_DIR}" -maxdepth 1 -type f -exec basename {} \; |
      grep -E '.*MP.*\.(fits|fits\.gz)$' |
      sort
    )
    if [[ "${#FILES[@]}" -eq 0 ]]; then
      echo "[SKIP] ${NIGHT} (no *MP* files in L2)"
      continue
    fi

    for FILE in "${FILES[@]}"; do
      STEM="${FILE%.fits.gz}"
      STEM="${STEM%.fits}"
      PART_ALL="${PARTS_DIR}/${STEM}_all_asteroids.fits"
      PART_MATCHED="${PARTS_DIR}/${STEM}_matched_asteroids.fits"
      if [[ -s "${PART_ALL}" || -s "${PART_MATCHED}" ]]; then
        echo "[SKIP] ${NIGHT} ${FILE} (part outputs already exist)"
        continue
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
    done

    "${PYTHON}" "${SCRIPT_DIR}/merge_night_parts.py" \
      "${NIGHT}" \
      --parts-dir "${PARTS_DIR}" \
      --outdir "${L4_DIR}"
  fi

  if [[ "${SUBMIT_MPC}" != "true" ]]; then
    echo "[DONE] ${NIGHT} extraction finished (submit disabled)"
    continue
  fi

  if [[ ! -s "${MATCHED_FITS}" ]]; then
    echo "[SKIP] ${NIGHT} submit (missing matched FITS)"
    continue
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
    --f-ratio "${F_RATIO}"
    --filter-name "${FILTER_NAME}"
    --array-size "${ARRAY_SIZE}"
    --pixel-scale "${PIXEL_SCALE}"
    --coinvestigators "${COINVESTIGATORS}"
    --collaborators "${COLLABORATORS}"
    --funding-source "${FUNDING_SOURCE}"
    --astcat "${ASTCAT}"
    --band "${BAND}"
    --photcat "${PHOTCAT}"
    --err-unit "${ERR_UNIT}"
    --min-observations "${MIN_OBSERVATIONS}"
    --ac2-email "${AC2_EMAIL}"
    --ack "${ACK_PREFIX}"
    --obj-type "${OBJ_TYPE}"
    --response-out "${REPLY_TXT}"
    --submit
  )
  if [[ -n "${PROG}" ]]; then
    EXPORT_CMD+=(--prog "${PROG}")
  fi
  if [[ "${INCLUDE_LOGSNR}" == "1" ]]; then
    EXPORT_CMD+=(--include-logsnr)
  fi
  "${EXPORT_CMD[@]}"
  echo "[DONE] ${NIGHT} submit finished"
done
