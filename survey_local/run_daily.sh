#!/bin/bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$PACKAGE_DIR/.." && pwd)"
PYTHON_ROOT="${PROJECT_ROOT}"
WORKSPACE="${WORKSPACE:-${PROJECT_ROOT}/local_data/survey/runtime}"
FOOTPRINTS="${FOOTPRINTS:-$PACKAGE_DIR/footprints/equatorial_footprints.fits}"
PROCESSED_ROOT="${PROCESSED_ROOT:-${PROJECT_ROOT}/local_data/processed1}"
PUBLISH_DIR="${PUBLISH_DIR:-${PROJECT_ROOT}/local_data/publish/survey}"
LOG_DIR="${LOG_DIR:-$PACKAGE_DIR/logs}"
mkdir -p "$LOG_DIR"

RUN_DATE="${1:-$(TZ=Asia/Shanghai date +%F)}"
LOG_PATH="$LOG_DIR/${RUN_DATE//-/}.log"

PYTHON_BIN="${PYTHON_BIN:-python3}"

log() {
  printf '[%s] [INFO] %s\n' "$(TZ=Asia/Shanghai date '+%F %T %Z')" "$*"
}

section() {
  printf '[%s] [INFO] ===== %s =====\n' "$(TZ=Asia/Shanghai date '+%F %T %Z')" "$*"
}

exec > >(tee -a "$LOG_PATH") 2>&1

section "cron wrapper start"
log "package_dir=$PACKAGE_DIR"
log "python_root=$PYTHON_ROOT"
log "workspace=$WORKSPACE"
log "footprints=$FOOTPRINTS"
log "processed_root=$PROCESSED_ROOT"
log "publish_dir=${PUBLISH_DIR:-<empty>}"
log "python_bin=$PYTHON_BIN"
log "run_date=$RUN_DATE"

status=0
trap 'status=$?; section "cron wrapper end"; log "exit_code=$status"' EXIT

PYTHONPATH="$PYTHON_ROOT${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON_BIN" -m survey_local.run_daily \
  --date "$RUN_DATE" \
  --root-dir "$PACKAGE_DIR" \
  --workspace "$WORKSPACE" \
  --footprints "$FOOTPRINTS" \
  --processed-root "$PROCESSED_ROOT" \
  --publish-dir "$PUBLISH_DIR"
