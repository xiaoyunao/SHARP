#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE="${WORKSPACE:-${ROOT_DIR}/runtime}"
PROCESSED_ROOT="${PROCESSED_ROOT:-/processed1}"
LOG_DIR="${WORKSPACE}/logs"
MAX_PARALLEL="${MAX_PARALLEL:-24}"
PYTHON="${PYTHON:-/home/smtpipeline/Softwares/miniconda3/bin/python}"
mkdir -p "${LOG_DIR}"

RUN_DATE="${1:-$(TZ=Asia/Shanghai date +%F)}"
TARGET_NIGHT="${TARGET_NIGHT:-$(TZ=Asia/Shanghai date -d "${RUN_DATE} -1 day" +%Y%m%d)}"
LOG_PATH="${LOG_DIR}/${RUN_DATE//-/}_daily.log"

{
  echo "[INFO] run_date=${RUN_DATE} target_night=${TARGET_NIGHT}"

  NIGHT_DIR="${PROCESSED_ROOT}/${TARGET_NIGHT}"
  L2_DIR="${NIGHT_DIR}/L2"
  L4_DIR="${NIGHT_DIR}/L4"
  OUT_PSV="${L4_DIR}/${TARGET_NIGHT}_matched_asteroids_ades.psv"
  REPLY_TXT="${L4_DIR}/${TARGET_NIGHT}_mpc_reply.txt"
  STATUS_JSON="${L4_DIR}/${TARGET_NIGHT}_known_asteroid_status.json"

  if [[ ! -d "${NIGHT_DIR}" ]]; then
    echo "[SKIP] missing night dir: ${NIGHT_DIR}"
    exit 0
  fi
  if [[ ! -d "${L2_DIR}" ]]; then
    echo "[SKIP] missing L2 dir: ${L2_DIR}"
    exit 0
  fi
  MP_COUNT="$(
    find "${L2_DIR}" -maxdepth 1 -type f -exec basename {} \; |
    grep -Ec '.*MP.*\.(fits|fits\.gz)$' || true
  )"
  if [[ "${MP_COUNT}" == "0" ]]; then
    echo "[SKIP] no *MP* files in ${L2_DIR}"
    exit 0
  fi
  if [[ -s "${OUT_PSV}" && -s "${REPLY_TXT}" ]]; then
    echo "[SKIP] report products already exist for ${TARGET_NIGHT}"
    exit 0
  fi
  if [[ -s "${STATUS_JSON}" ]]; then
    if "${PYTHON}" - <<PY
import json, sys
from pathlib import Path
status = json.loads(Path("${STATUS_JSON}").read_text(encoding="utf-8"))
matched = status.get("matched_asteroids", {})
ok = bool(status.get("known_complete")) and int(matched.get("rows") or 0) == 0
sys.exit(0 if ok else 1)
PY
    then
      echo "[SKIP] known extraction complete with zero 1.0 arcsec matched detections for ${TARGET_NIGHT}"
      exit 0
    fi
  fi

  echo "[RUN] submit night ${TARGET_NIGHT} with max_parallel=${MAX_PARALLEL}"
  PYTHON="${PYTHON}" "${ROOT_DIR}/submit_pipeline_slurm.sh" \
    --batch false \
    --submit-mpc true \
    --max-parallel "${MAX_PARALLEL}" \
    "${TARGET_NIGHT}"
} 2>&1 | tee -a "${LOG_PATH}"
