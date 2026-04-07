#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_backfill.sh --submit-mpc false START_YYYYMMDD END_YYYYMMDD

Notes:
  - processes one night at a time
  - can run extraction only or extraction + report
  - waits for the finalize job to finish before moving to the next night
EOF
}

BATCH_SUBMIT_MPC="true"
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --submit-mpc)
      BATCH_SUBMIT_MPC="$2"
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

if [[ $# -ne 2 ]]; then
  usage >&2
  exit 1
fi

START_NIGHT="$1"
END_NIGHT="$2"

if [[ ! "${START_NIGHT}" =~ ^[0-9]{8}$ || ! "${END_NIGHT}" =~ ^[0-9]{8}$ ]]; then
  echo "[FATAL] nights must be YYYYMMDD" >&2
  exit 2
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${ROOT_DIR}/.." && pwd)"
PROCESSED_ROOT="${PROCESSED_ROOT:-${PROJECT_ROOT}/local_data/processed1}"
WORKSPACE="${WORKSPACE:-${PROJECT_ROOT}/local_data/known_asteroid/runtime}"
LOG_DIR="${WORKSPACE}/logs"
MAX_PARALLEL="${MAX_PARALLEL:-24}"
POLL_SECONDS="${POLL_SECONDS:-60}"
mkdir -p "${LOG_DIR}"

LOG_PATH="${LOG_DIR}/${START_NIGHT}_${END_NIGHT}_backfill.log"

wait_for_job() {
  local job_id="$1"
  local state=""
  while true; do
    if squeue -h -j "${job_id}" | grep -q .; then
      sleep "${POLL_SECONDS}"
      continue
    fi
    if scontrol show job "${job_id}" >/tmp/known_asteroid_scontrol_${job_id}.txt 2>/dev/null; then
      state="$(awk -F= '/JobState=/{print $2; exit}' /tmp/known_asteroid_scontrol_${job_id}.txt | awk '{print $1}')"
      rm -f "/tmp/known_asteroid_scontrol_${job_id}.txt"
      if [[ -n "${state}" ]]; then
        echo "${state}"
        return 0
      fi
    fi
    echo "COMPLETED"
    return 0
  done
}

{
  echo "[INFO] backfill range ${START_NIGHT} -> ${END_NIGHT}"

  while IFS= read -r NIGHT; do
    NIGHT_DIR="${PROCESSED_ROOT}/${NIGHT}"
    L2_DIR="${NIGHT_DIR}/L2"
    L4_DIR="${NIGHT_DIR}/L4"
    OUT_PSV="${L4_DIR}/${NIGHT}_matched_asteroids_ades.psv"
    REPLY_TXT="${L4_DIR}/${NIGHT}_mpc_reply.txt"

    if [[ ! -d "${NIGHT_DIR}" ]]; then
      echo "[SKIP] ${NIGHT} missing night dir"
      continue
    fi
    if [[ ! -d "${L2_DIR}" ]]; then
      echo "[SKIP] ${NIGHT} missing L2 dir"
      continue
    fi
    MP_COUNT="$(
      find "${L2_DIR}" -maxdepth 1 -type f -exec basename {} \; |
      grep -Ec '.*MP.*\.(fits|fits\.gz)$' || true
    )"
    if [[ "${MP_COUNT}" == "0" ]]; then
      echo "[SKIP] ${NIGHT} no *MP* files"
      continue
    fi
    if [[ "${BATCH_SUBMIT_MPC}" == "true" && -s "${OUT_PSV}" && -s "${REPLY_TXT}" ]]; then
      echo "[SKIP] ${NIGHT} report products already exist"
      continue
    fi

    echo "[RUN] ${NIGHT} submit_mpc=${BATCH_SUBMIT_MPC} max_parallel=${MAX_PARALLEL}"
    OUTPUT="$(ROOT_DIR="${PROCESSED_ROOT}" ASTORB="${PROJECT_ROOT}/resources/known_asteroid/astorb.dat" PYTHON="${PYTHON:-python3}" \
      "${ROOT_DIR}/run_pipeline.sh" \
        --batch false \
        --submit-mpc "${BATCH_SUBMIT_MPC}" \
        "${NIGHT}")"
    printf '%s\n' "${OUTPUT}"

    echo "[DONE] ${NIGHT} sequential local run finished"
  done < <(
    find "${PROCESSED_ROOT}" -maxdepth 1 -mindepth 1 -type d -exec basename {} \; |
    grep -E '^[0-9]{8}$' |
    awk -v s="${START_NIGHT}" -v e="${END_NIGHT}" '$0 >= s && $0 <= e' |
    sort
  )
} 2>&1 | tee -a "${LOG_PATH}"
