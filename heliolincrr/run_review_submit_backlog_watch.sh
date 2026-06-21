#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

START_NIGHT="${1:-${START_NIGHT:-20251116}}"
END_NIGHT="${2:-${END_NIGHT:-20260617}}"
PYTHON_BIN="${PYTHON_BIN:-/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python}"
REVIEW_ROOT="${REVIEW_ROOT:-/pipeline/xiaoyunao/heliolincrr/review_packages}"
PROCESSED_ROOT="${PROCESSED_ROOT:-/processed1}"
STATE_PATH="${STATE_PATH:-/pipeline/xiaoyunao/data/heliolincrr/review_submit_backlog_${START_NIGHT}_${END_NIGHT}.json}"
LOG_DIR="${LOG_DIR:-/pipeline/xiaoyunao/data/heliolincrr/daily_logs}"
INTERVAL_SEC="${INTERVAL_SEC:-300}"
VALIDATE="${VALIDATE:-1}"
SUBMIT="${SUBMIT:-1}"
NO_LOGSNR="${NO_LOGSNR:-0}"
RETRY_FAILED="${RETRY_FAILED:-1}"

mkdir -p "${LOG_DIR}"
LOG_PATH="${LOG_PATH:-${LOG_DIR}/review_submit_backlog_${START_NIGHT}_${END_NIGHT}.log}"

if pgrep -af "watch_submit_reviews.py .*--state ${STATE_PATH}" >/dev/null 2>&1; then
  echo "[SKIP] backlog watcher already running for state ${STATE_PATH}"
  pgrep -af "watch_submit_reviews.py .*--state ${STATE_PATH}" || true
  exit 0
fi

CMD=(
  "${PYTHON_BIN}" "${SCRIPT_DIR}/watch_submit_reviews.py"
  --review-root "${REVIEW_ROOT}"
  --processed-root "${PROCESSED_ROOT}"
  --state "${STATE_PATH}"
  --start "${START_NIGHT}"
  --end "${END_NIGHT}"
  --follow
  --exit-when-complete
  --interval-sec "${INTERVAL_SEC}"
)

if [[ "${VALIDATE}" == "1" || "${SUBMIT}" == "1" ]]; then
  CMD+=(--validate)
fi
if [[ "${SUBMIT}" == "1" ]]; then
  CMD+=(--submit)
fi
if [[ "${NO_LOGSNR}" == "1" ]]; then
  CMD+=(--no-logsnr)
fi
if [[ "${RETRY_FAILED}" == "1" ]]; then
  CMD+=(--retry-failed)
fi

echo "[RUN] ${CMD[*]}"
echo "[LOG] ${LOG_PATH}"
nohup "${CMD[@]}" >> "${LOG_PATH}" 2>&1 &
PID=$!
echo "[STARTED] pid=${PID} state=${STATE_PATH}"
