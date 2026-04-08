#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./submit_pipeline_slurm.sh --batch false --submit-mpc false [--max-parallel 24] YYYYMMDD
  ./submit_pipeline_slurm.sh --batch true --submit-mpc false [--max-parallel 24] START_YYYYMMDD END_YYYYMMDD

Notes:
  - Extraction uses a slurm array with one task per *MP* file
  - Final merge/submit stays serial in a single dependent job
  - --max-parallel limits how many file-level tasks run at once
EOF
}

if [[ $# -lt 3 ]]; then
  usage >&2
  exit 1
fi

BATCH="false"
SUBMIT_MPC="false"
MAX_PARALLEL="${MAX_PARALLEL:-24}"
MATCH_CPUS="${MATCH_CPUS:-1}"
MATCH_MEM="${MATCH_MEM:-10G}"
FINALIZE_CPUS="${FINALIZE_CPUS:-1}"
FINALIZE_MEM="${FINALIZE_MEM:-10G}"
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch) BATCH="$2"; shift 2 ;;
    --submit-mpc) SUBMIT_MPC="$2"; shift 2 ;;
    --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) POSITIONAL+=("$1"); shift ;;
  esac
done
set -- "${POSITIONAL[@]}"

ROOT_DIR="${ROOT_DIR:-/processed1}"
SCRIPT_DIR="${SCRIPT_DIR:-/pipeline/xiaoyunao/known_asteroid}"

if [[ ! "${MAX_PARALLEL}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[FATAL] --max-parallel must be a positive integer" >&2
  exit 2
fi

if [[ "${BATCH}" == "true" ]]; then
  [[ $# -eq 2 ]] || { usage >&2; exit 1; }
  START_NIGHT="$1"
  END_NIGHT="$2"
else
  [[ $# -eq 1 ]] || { usage >&2; exit 1; }
  START_NIGHT="$1"
  END_NIGHT="$1"
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

for NIGHT in "${NIGHTS[@]}"; do
  L2_DIR="${ROOT_DIR}/${NIGHT}/L2"
  L4_DIR="${ROOT_DIR}/${NIGHT}/L4"
  PARTS_DIR="${L4_DIR}/known_asteroid_parts"
  ALL_FITS="${L4_DIR}/${NIGHT}_all_asteroids.fits"
  MATCHED_FITS="${L4_DIR}/${NIGHT}_matched_asteroids.fits"
  OUT_PSV="${L4_DIR}/${NIGHT}_matched_asteroids_ades.psv"
  REPLY_TXT="${L4_DIR}/${NIGHT}_mpc_reply.txt"
  MANIFEST="${L4_DIR}/${NIGHT}_file_manifest.txt"

  mkdir -p "${L4_DIR}" "${PARTS_DIR}"

  if [[ "${SUBMIT_MPC}" == "true" && -s "${OUT_PSV}" && -s "${REPLY_TXT}" ]]; then
    echo "[SKIP] ${NIGHT} report products already exist"
    continue
  fi

  ARRAY_JOB=""
  if [[ -s "${ALL_FITS}" && -s "${MATCHED_FITS}" ]]; then
    echo "[SKIP] ${NIGHT} extraction already complete"
  else
    : > "${MANIFEST}"
    while IFS= read -r FILE_PATH; do
      FILE="$(basename "${FILE_PATH}")"
      STEM="${FILE%.fits.gz}"
      STEM="${STEM%.fits}"
      PART_ALL="${PARTS_DIR}/${STEM}_all_asteroids.fits"
      PART_MATCHED="${PARTS_DIR}/${STEM}_matched_asteroids.fits"
      if [[ -s "${PART_ALL}" || -s "${PART_MATCHED}" ]]; then
        continue
      fi
      echo "${FILE}" >> "${MANIFEST}"
    done < <(find "${L2_DIR}" -maxdepth 1 -type f | grep -E '.*MP.*\.(fits|fits\.gz)$' | sort)

    if [[ -s "${MANIFEST}" ]]; then
      NFILES=$(wc -l < "${MANIFEST}")
      LAST_INDEX=$((NFILES - 1))
      ARRAY_LIMIT="${MAX_PARALLEL}"
      if (( ARRAY_LIMIT > NFILES )); then
        ARRAY_LIMIT="${NFILES}"
      fi
      ARRAY_JOB=$(sbatch \
        --parsable \
        --cpus-per-task="${MATCH_CPUS}" \
        --mem="${MATCH_MEM}" \
        --array="0-${LAST_INDEX}%${ARRAY_LIMIT}" \
        "${SCRIPT_DIR}/slurm_match_one_file.sh" "${NIGHT}" "${MANIFEST}")
      echo "[SUBMIT] ${NIGHT} array job ${ARRAY_JOB} (${NFILES} files, max_parallel=${ARRAY_LIMIT}, cpu=${MATCH_CPUS}, mem=${MATCH_MEM})"
    else
      echo "[INFO] ${NIGHT} no remaining files to process"
    fi
  fi

  if [[ -n "${ARRAY_JOB}" ]]; then
    MERGE_JOB=$(sbatch \
      --parsable \
      --cpus-per-task="${FINALIZE_CPUS}" \
      --mem="${FINALIZE_MEM}" \
      --dependency="afterok:${ARRAY_JOB}" \
      "${SCRIPT_DIR}/slurm_merge_submit.sh" "${NIGHT}" "${SUBMIT_MPC}")
  else
    MERGE_JOB=$(sbatch \
      --parsable \
      --cpus-per-task="${FINALIZE_CPUS}" \
      --mem="${FINALIZE_MEM}" \
      "${SCRIPT_DIR}/slurm_merge_submit.sh" "${NIGHT}" "${SUBMIT_MPC}")
  fi
  echo "[SUBMIT] ${NIGHT} finalize job ${MERGE_JOB} (cpu=${FINALIZE_CPUS}, mem=${FINALIZE_MEM}, submit=${SUBMIT_MPC})"
done
