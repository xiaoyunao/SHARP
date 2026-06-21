#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_one_night_direct.sh YYYYMMDD [submit_mpc=false]

Runs one observing night without submitting a Slurm array.  Existing per-file
parts are reused only when all required outputs are newer than MIN_PART_MTIME.
This is intended for controlled recovery/rematch runs after code fixes.
EOF
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage >&2
  exit 1
fi

NIGHT="$1"
SUBMIT_MPC="${2:-false}"

ROOT_DIR="${ROOT_DIR:-/processed1}"
SCRIPT_DIR="${SCRIPT_DIR:-/pipeline/xiaoyunao/known_asteroid}"
PYTHON="${PYTHON:-/home/smtpipeline/Softwares/miniconda3/bin/python}"
MAX_PARALLEL="${MAX_PARALLEL:-16}"

FILE_REGEX="${FILE_REGEX:-.*MP.*\\.(fits|fits\\.gz)$}"
HDU="${HDU:-1}"
NX="${NX:-9216}"
NY="${NY:-9232}"
OBS_DATE_KEY="${OBS_DATE_KEY:-OBS_DATE}"
OBS_TIMEZONE="${OBS_TIMEZONE:-Asia/Shanghai}"
NIGHT_ROLLOVER_HOUR="${NIGHT_ROLLOVER_HOUR:-12}"
EXPTIME_KEY="${EXPTIME_KEY:-}"
EXPTIME_SEC="${EXPTIME_SEC:-30.0}"
SEP_ARCSEC="${SEP_ARCSEC:-1.0}"
MASK_SEP_ARCSEC="${MASK_SEP_ARCSEC:-1.5}"
MASK_MATCHED_SUFFIX="${MASK_MATCHED_SUFFIX:-_mask15}"
MAG_LIMIT="${MAG_LIMIT:-22.5}"
MAGDIFF="${MAGDIFF:-99.0}"
CAT_MAG_COL="${CAT_MAG_COL:-Mag_Kron}"
CONFIDENCE_PAD_DEG="${CONFIDENCE_PAD_DEG:-0.5}"
CORNER_PAD_DEG="${CORNER_PAD_DEG:-0.05}"
LON="${LON:-117.575}"
LAT="${LAT:-40.393}"
HEIGHT="${HEIGHT:-960.0}"
VERBOSE="${VERBOSE:-0}"
ASTORB="${ASTORB:-${SCRIPT_DIR}/astorb.dat}"
FORCE_ALL_FILES="${FORCE_ALL_FILES:-0}"
MIN_PART_MTIME="${MIN_PART_MTIME:-2026-06-19 12:25:00}"

if [[ ! "${MAX_PARALLEL}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[FATAL] MAX_PARALLEL must be a positive integer" >&2
  exit 2
fi

L2_DIR="${ROOT_DIR}/${NIGHT}/L2"
L4_DIR="${ROOT_DIR}/${NIGHT}/L4"
PARTS_DIR="${L4_DIR}/known_asteroid_parts"
MANIFEST="${L4_DIR}/${NIGHT}_file_manifest.txt"
RUN_MANIFEST="${L4_DIR}/${NIGHT}_file_manifest_direct.txt"
MATCH_LOG="${L4_DIR}/${NIGHT}_match.log"
MARKER="${L4_DIR}/.known_wrapfix_min_part_mtime"

if [[ ! -d "${L2_DIR}" ]]; then
  echo "[SKIP] ${NIGHT} missing L2 dir: ${L2_DIR}"
  exit 0
fi

mkdir -p "${L4_DIR}" "${PARTS_DIR}"
touch -d "${MIN_PART_MTIME}" "${MARKER}"

"${PYTHON}" "${SCRIPT_DIR}/build_file_manifest.py" "${NIGHT}" \
  --l2-dir "${L2_DIR}" \
  --out "${MANIFEST}" \
  --file-regex "${FILE_REGEX}" \
  --hdu "${HDU}" \
  --obs-date-key "${OBS_DATE_KEY}" \
  --timezone "${OBS_TIMEZONE}" \
  --night-rollover-hour "${NIGHT_ROLLOVER_HOUR}"

: > "${RUN_MANIFEST}"
while IFS= read -r FILE; do
  [[ -n "${FILE}" ]] || continue
  STEM="${FILE%.fits.gz}"
  STEM="${STEM%.fits}"
  PART_ALL="${PARTS_DIR}/${STEM}_all_asteroids.fits"
  PART_MATCHED="${PARTS_DIR}/${STEM}_matched_asteroids.fits"
  PART_MASK_MATCHED="${PARTS_DIR}/${STEM}_matched_asteroids${MASK_MATCHED_SUFFIX}.fits"
  if [[ "${FORCE_ALL_FILES}" != "1" \
      && -s "${PART_ALL}" && -s "${PART_MATCHED}" && -s "${PART_MASK_MATCHED}" \
      && "${PART_ALL}" -nt "${MARKER}" \
      && "${PART_MATCHED}" -nt "${MARKER}" \
      && "${PART_MASK_MATCHED}" -nt "${MARKER}" ]]; then
    continue
  fi
  echo "${FILE}" >> "${RUN_MANIFEST}"
done < "${MANIFEST}"

run_one_file() {
  local file="$1"
  local cmd=(
    "${PYTHON}" "${SCRIPT_DIR}/match_single_night.py" "${NIGHT}"
    --root "${ROOT_DIR}"
    --outdir "${L4_DIR}"
    --parts-dir "${PARTS_DIR}"
    --file "${file}"
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
  if [[ -n "${MASK_SEP_ARCSEC}" ]]; then
    cmd+=(--mask-sep-arcsec "${MASK_SEP_ARCSEC}" --mask-matched-suffix "${MASK_MATCHED_SUFFIX}")
  fi
  if [[ -n "${EXPTIME_KEY}" ]]; then
    cmd+=(--exptime-key "${EXPTIME_KEY}")
  fi
  if [[ "${VERBOSE}" == "1" ]]; then
    cmd+=(--verbose)
  fi
  "${cmd[@]}"
}

export NIGHT ROOT_DIR SCRIPT_DIR PYTHON FILE_REGEX HDU NX NY OBS_DATE_KEY EXPTIME_KEY EXPTIME_SEC
export SEP_ARCSEC MASK_SEP_ARCSEC MASK_MATCHED_SUFFIX MAG_LIMIT MAGDIFF CAT_MAG_COL
export CONFIDENCE_PAD_DEG CORNER_PAD_DEG LON LAT HEIGHT VERBOSE ASTORB L4_DIR PARTS_DIR MATCH_LOG
export -f run_one_file

NFILES_TOTAL=$(wc -l < "${MANIFEST}" | tr -d ' ')
NFILES_RUN=$(wc -l < "${RUN_MANIFEST}" | tr -d ' ')
echo "[DIRECT] ${NIGHT} manifest=${NFILES_TOTAL} to_run=${NFILES_RUN} max_parallel=${MAX_PARALLEL} min_part_mtime='${MIN_PART_MTIME}'"

if [[ -s "${RUN_MANIFEST}" ]]; then
  xargs -r -n 1 -P "${MAX_PARALLEL}" bash -c 'run_one_file "$1"' _ < "${RUN_MANIFEST}"
else
  echo "[DIRECT] ${NIGHT} no stale or missing file parts"
fi

FORCE_MERGE=1 FORCE_EXTRACT=1 MASK_MATCHED_SUFFIX="${MASK_MATCHED_SUFFIX}" \
  "${SCRIPT_DIR}/slurm_merge_submit.sh" "${NIGHT}" "${SUBMIT_MPC}"
