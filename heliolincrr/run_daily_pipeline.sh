#!/usr/bin/env bash
set -euo pipefail

RUN_DATE="${1:-$(TZ=Asia/Shanghai date +%F)}"
BASE_DIR="${BASE_DIR:-/pipeline/xiaoyunao}"
SURVEY_DIR="${SURVEY_DIR:-${BASE_DIR}/survey}"
KNOWN_DIR="${KNOWN_DIR:-${BASE_DIR}/known_asteroid}"
UNKNOWN_DIR="${UNKNOWN_DIR:-${BASE_DIR}/heliolincrr}"
LOG_DIR="${LOG_DIR:-/pipeline/xiaoyunao/data/heliolincrr/daily_logs}"
RECOVERY_LOOKBACK_DAYS="${RECOVERY_LOOKBACK_DAYS:-7}"
ENABLE_RECOVERY_PROCESSING="${ENABLE_RECOVERY_PROCESSING:-1}"
PROCESSED_ROOT="${PROCESSED_ROOT:-/processed1}"
KNOWN_PYTHON="${KNOWN_PYTHON:-/home/smtpipeline/Softwares/miniconda3/bin/python}"
WAIT_KNOWN_REPORT_TIMEOUT_SEC="${WAIT_KNOWN_REPORT_TIMEOUT_SEC:-28800}"
WAIT_KNOWN_REPORT_INTERVAL_SEC="${WAIT_KNOWN_REPORT_INTERVAL_SEC:-120}"

mkdir -p "${LOG_DIR}"
LOG_PATH="${LOG_DIR}/${RUN_DATE//-/}_daily_pipeline.log"

section() {
  printf '\n[%s] %s\n' "$(TZ=Asia/Shanghai date '+%F %T %Z')" "$*"
}

night_for_offset() {
  local offset="$1"
  TZ=Asia/Shanghai date -d "${RUN_DATE} -${offset} day" +%Y%m%d
}

run_date_after_night() {
  local night="$1"
  TZ=Asia/Shanghai date -d "${night:0:4}-${night:4:2}-${night:6:2} +1 day" +%F
}

has_mp_l2() {
  local night="$1"
  local l2_dir="${PROCESSED_ROOT}/${night}/L2"
  [[ -d "${l2_dir}" ]] || return 1
  [[ -n "$(
    find "${l2_dir}" -maxdepth 1 -type f \( -name '*MP*.fits' -o -name '*MP*.fits.gz' \) -print -quit
  )" ]]
}

known_report_ready() {
  local night="$1"
  local l4_dir="${PROCESSED_ROOT}/${night}/L4"
  local status="${l4_dir}/${night}_known_asteroid_status.json"
  local out_psv="${l4_dir}/${night}_matched_asteroids_ades.psv"
  local reply_txt="${l4_dir}/${night}_mpc_reply.txt"

  if [[ -s "${out_psv}" && -s "${reply_txt}" ]]; then
    echo "[INFO] known report products ready for ${night}"
    return 0
  fi
  if [[ -s "${status}" ]]; then
    if "${KNOWN_PYTHON}" - <<PY
import json, sys
from pathlib import Path
status = json.loads(Path("${status}").read_text(encoding="utf-8"))
matched = status.get("matched_asteroids", {})
ok = bool(status.get("known_complete")) and int(matched.get("rows") or 0) == 0
sys.exit(0 if ok else 1)
PY
    then
      echo "[INFO] known complete with zero 1.0 arcsec matched detections for ${night}"
      return 0
    fi
  fi
  return 1
}

wait_for_known_report() {
  local night="$1"
  local waited=0
  if [[ ! -d "${PROCESSED_ROOT}/${night}" ]]; then
    echo "[SKIP] missing night dir while waiting known report: ${PROCESSED_ROOT}/${night}"
    return 1
  fi
  if ! has_mp_l2 "${night}"; then
    echo "[SKIP] no MP L2 input while waiting known report: ${night}"
    return 1
  fi
  while true; do
    if known_report_ready "${night}"; then
      return 0
    fi
    if (( waited >= WAIT_KNOWN_REPORT_TIMEOUT_SEC )); then
      echo "[SKIP] known report not ready after ${WAIT_KNOWN_REPORT_TIMEOUT_SEC}s for ${night}"
      return 1
    fi
    echo "[WAIT] known report not ready for ${night}; sleep ${WAIT_KNOWN_REPORT_INTERVAL_SEC}s"
    sleep "${WAIT_KNOWN_REPORT_INTERVAL_SEC}"
    waited=$((waited + WAIT_KNOWN_REPORT_INTERVAL_SEC))
  done
}

run_processing_for_night() {
  local night="$1"
  local log_date
  log_date="$(run_date_after_night "${night}")"

  section "known asteroid for ${night}"
  if [[ -x "${KNOWN_DIR}/run_daily.sh" ]]; then
    TARGET_NIGHT="${night}" bash "${KNOWN_DIR}/run_daily.sh" "${log_date}"
  else
    echo "[SKIP] missing known daily: ${KNOWN_DIR}/run_daily.sh"
  fi

  if ! wait_for_known_report "${night}"; then
    echo "[SKIP] unknown check package for ${night}: known report not ready"
    return 0
  fi

  section "unknown check package for ${night}"
  if [[ -x "${UNKNOWN_DIR}/run_daily_unknown.sh" ]]; then
    TARGET_NIGHT="${night}" bash "${UNKNOWN_DIR}/run_daily_unknown.sh" "${log_date}"
  else
    echo "[SKIP] missing unknown daily: ${UNKNOWN_DIR}/run_daily_unknown.sh"
  fi
}

main() {
  section "daily pipeline start"
  echo "[INFO] run_date=${RUN_DATE} recovery_lookback_days=${RECOVERY_LOOKBACK_DAYS}"

  section "survey plan"
  if [[ -x "${SURVEY_DIR}/run_daily.sh" ]]; then
    bash "${SURVEY_DIR}/run_daily.sh" "${RUN_DATE}"
  else
    echo "[SKIP] missing survey daily: ${SURVEY_DIR}/run_daily.sh"
  fi

  if [[ "${ENABLE_RECOVERY_PROCESSING}" == "1" ]]; then
    local offset night
    for ((offset=RECOVERY_LOOKBACK_DAYS; offset>=1; offset--)); do
      night="$(night_for_offset "${offset}")"
      run_processing_for_night "${night}"
    done
  else
    run_processing_for_night "$(night_for_offset 1)"
  fi

  section "daily pipeline done"
}

main 2>&1 | tee -a "${LOG_PATH}"
