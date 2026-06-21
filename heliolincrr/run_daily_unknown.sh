#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RUN_DATE="${1:-$(TZ=Asia/Shanghai date +%F)}"
TARGET_NIGHT="${TARGET_NIGHT:-$(TZ=Asia/Shanghai date -d "${RUN_DATE} -1 day" +%Y%m%d)}"

PROCESSED_ROOT="${PROCESSED_ROOT:-/processed1}"
ROOT_OUT="${ROOT_OUT:-/pipeline/xiaoyunao/data/heliolincrr}"
PLOT_ROOT="${PLOT_ROOT:-/pipeline/xiaoyunao/heliolincrr/plots}"
REVIEW_ROOT="${REVIEW_ROOT:-/pipeline/xiaoyunao/heliolincrr/review_packages}"
PYTHON_BIN="${PYTHON_BIN:-/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python}"
LOG_DIR="${LOG_DIR:-${ROOT_OUT}/daily_logs}"
WAIT_KNOWN_TIMEOUT_SEC="${WAIT_KNOWN_TIMEOUT_SEC:-28800}"
WAIT_KNOWN_INTERVAL_SEC="${WAIT_KNOWN_INTERVAL_SEC:-120}"
FORCE_UNKNOWN="${FORCE_UNKNOWN:-0}"
SKIP_PLOTS="${SKIP_PLOTS:-0}"
MAX_UNKNOWN_LINKS_AFTER_KNOWN="${MAX_UNKNOWN_LINKS_AFTER_KNOWN:-200}"
START_SUBMIT_WATCHER="${START_SUBMIT_WATCHER:-1}"
SUBMIT_WATCH_STATE="${SUBMIT_WATCH_STATE:-${ROOT_OUT}/review_submit_state.json}"
SUBMIT_WATCH_LOG="${SUBMIT_WATCH_LOG:-${LOG_DIR}/review_submit_watch.log}"
SUBMIT_WATCH_INTERVAL_SEC="${SUBMIT_WATCH_INTERVAL_SEC:-300}"
SUBMIT_WATCH_VALIDATE="${SUBMIT_WATCH_VALIDATE:-1}"
SUBMIT_WATCH_SUBMIT="${SUBMIT_WATCH_SUBMIT:-1}"
SUBMIT_WATCH_NO_LOGSNR="${SUBMIT_WATCH_NO_LOGSNR:-0}"
SUBMIT_WATCH_RETRY_FAILED="${SUBMIT_WATCH_RETRY_FAILED:-1}"

mkdir -p "${LOG_DIR}"
LOG_PATH="${LOG_DIR}/${RUN_DATE//-/}_unknown_daily.log"

section() {
  printf '\n[%s] %s\n' "$(TZ=Asia/Shanghai date '+%F %T %Z')" "$*"
}

log() {
  printf '%s\n' "$*"
}

has_mp_l2() {
  local night="$1"
  local l2_dir="${PROCESSED_ROOT}/${night}/L2"
  [[ -d "${l2_dir}" ]] || return 1
  [[ -n "$(
    find "${l2_dir}" -maxdepth 1 -type f \( -name '*MP*.fits' -o -name '*MP*.fits.gz' \) -print -quit
  )" ]]
}

known_ready() {
  local night="$1"
  local l4_dir="${PROCESSED_ROOT}/${night}/L4"
  local mask15="${l4_dir}/${night}_matched_asteroids_mask15.fits"
  local status="${l4_dir}/${night}_known_asteroid_status.json"
  if [[ -s "${mask15}" ]]; then
    log "[INFO] known mask15 ready: ${mask15}"
    return 0
  fi
  if [[ -s "${status}" ]]; then
    if "${PYTHON_BIN}" - <<PY
import json, sys
from pathlib import Path
status = json.loads(Path("${status}").read_text(encoding="utf-8"))
mask = status.get("mask_matched_asteroids", {})
ok = bool(status.get("known_complete")) and int(mask.get("rows") or 0) == 0 and not bool(mask.get("exists"))
sys.exit(0 if ok else 1)
PY
    then
      log "[INFO] known complete with empty 1.5 arcsec mask: ${status}"
      return 0
    fi
  fi
  return 1
}

wait_for_known() {
  local night="$1"
  local waited=0
  while true; do
    if known_ready "${night}"; then
      return 0
    fi
    if (( waited >= WAIT_KNOWN_TIMEOUT_SEC )); then
      log "[SKIP] known mask15 not ready after ${WAIT_KNOWN_TIMEOUT_SEC}s for ${night}"
      return 1
    fi
    log "[WAIT] known mask15 not ready for ${night}; sleep ${WAIT_KNOWN_INTERVAL_SEC}s"
    sleep "${WAIT_KNOWN_INTERVAL_SEC}"
    waited=$((waited + WAIT_KNOWN_INTERVAL_SEC))
  done
}

review_unknown_count() {
  local manifest="$1"
  "${PYTHON_BIN}" - <<PY
import json
from pathlib import Path
data = json.loads(Path("${manifest}").read_text(encoding="utf-8"))
print(int(data.get("n_catalog_rows") or data.get("review_full_rows") or 0))
PY
}

start_submit_watcher() {
  local night="$1"
  local watcher_log="${SUBMIT_WATCH_LOG}"
  local cmd=(
    "${PYTHON_BIN}" "${SCRIPT_DIR}/watch_submit_reviews.py"
    --review-root "${REVIEW_ROOT}"
    --processed-root "${PROCESSED_ROOT}"
    --state "${SUBMIT_WATCH_STATE}"
    --start "${night}"
    --end "${night}"
    --only-night "${night}"
    --follow
    --exit-when-complete
    --interval-sec "${SUBMIT_WATCH_INTERVAL_SEC}"
  )
  if [[ "${SUBMIT_WATCH_VALIDATE}" == "1" || "${SUBMIT_WATCH_SUBMIT}" == "1" ]]; then
    cmd+=(--validate)
  fi
  if [[ "${SUBMIT_WATCH_SUBMIT}" == "1" ]]; then
    cmd+=(--submit)
  fi
  if [[ "${SUBMIT_WATCH_NO_LOGSNR}" == "1" ]]; then
    cmd+=(--no-logsnr)
  fi
  if [[ "${SUBMIT_WATCH_RETRY_FAILED}" == "1" ]]; then
    cmd+=(--retry-failed)
  fi

  if pgrep -af "watch_submit_reviews.py .*--only-night ${night}" >/dev/null 2>&1; then
    log "[SKIP] submit watcher already running for ${night}"
    return 0
  fi
  log "[RUN] start submit watcher for ${night}; log=${watcher_log}"
  nohup "${cmd[@]}" >> "${watcher_log}" 2>&1 &
  log "[INFO] submit watcher pid=$!"
}

main() {
  section "unknown daily start"
  log "[INFO] run_date=${RUN_DATE} target_night=${TARGET_NIGHT}"

  local night_dir="${PROCESSED_ROOT}/${TARGET_NIGHT}"
  local l4_dir="${night_dir}/L4"
  local unknown_json="${l4_dir}/${TARGET_NIGHT}_unknown_links.json"
  local review_manifest="${REVIEW_ROOT}/${TARGET_NIGHT}/${TARGET_NIGHT}_unknown_review_manifest.json"

  if [[ ! -d "${night_dir}" ]]; then
    log "[SKIP] missing night dir: ${night_dir}"
    return 0
  fi
  if ! has_mp_l2 "${TARGET_NIGHT}"; then
    log "[SKIP] no MP L2 input for ${TARGET_NIGHT}"
    return 0
  fi
  if [[ "${FORCE_UNKNOWN}" != "1" && -s "${unknown_json}" && -s "${review_manifest}" ]]; then
    log "[SKIP] unknown review package already exists: ${review_manifest}"
    return 0
  fi
  if ! wait_for_known "${TARGET_NIGHT}"; then
    return 0
  fi

  section "run unknown extraction"
  set +e
  (
    cd "${SCRIPT_DIR}"
    REQUIRE_KNOWN_MASK15=1 \
    MAX_UNKNOWN_LINKS_AFTER_KNOWN="${MAX_UNKNOWN_LINKS_AFTER_KNOWN}" \
    SKIP_PLOTS="${SKIP_PLOTS}" \
    PYTHON_BIN="${PYTHON_BIN}" \
    ./run_single_night.sh "${TARGET_NIGHT}"
  )
  local rc=$?
  set -e
  if [[ "${rc}" -eq 20 ]]; then
    log "[SKIP] ${TARGET_NIGHT} unknown count exceeds MAX_UNKNOWN_LINKS_AFTER_KNOWN=${MAX_UNKNOWN_LINKS_AFTER_KNOWN}"
    return 0
  fi
  if [[ "${rc}" -ne 0 ]]; then
    log "[ERROR] run_single_night failed rc=${rc}"
    return "${rc}"
  fi

  section "package review"
  "${PYTHON_BIN}" "${SCRIPT_DIR}/package_unknown_review.py" "${TARGET_NIGHT}" \
    --processed-root "${PROCESSED_ROOT}" \
    --root-out "${ROOT_OUT}" \
    --plot-root "${PLOT_ROOT}" \
    --out-root "${REVIEW_ROOT}" \
    --make-tar

  local n_unknown
  n_unknown="$(review_unknown_count "${review_manifest}")"
  if [[ "${START_SUBMIT_WATCHER}" == "1" && "${n_unknown}" -gt 0 ]]; then
    start_submit_watcher "${TARGET_NIGHT}"
  else
    log "[INFO] submit watcher not started for ${TARGET_NIGHT}: START_SUBMIT_WATCHER=${START_SUBMIT_WATCHER} unknown_count=${n_unknown}"
  fi

  section "unknown daily done"
}

main 2>&1 | tee -a "${LOG_PATH}"
