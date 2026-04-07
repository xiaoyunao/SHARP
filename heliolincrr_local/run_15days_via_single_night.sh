#!/usr/bin/env bash
set -euo pipefail

# =========================================
# Run run_single_night.sh for 15 consecutive nights
#
# Usage:
#   ./run_15days_via_single_night.sh START_NIGHT
# Example:
#   ./run_15days_via_single_night.sh 20260115
# =========================================

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 START_NIGHT(YYYYMMDD)"
  exit 1
fi

START="$1"
NDAYS=15
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_HOME_DEFAULT="${CONDA_HOME_DEFAULT:-/home/smtpipeline/Softwares/miniconda3}"
HELIOLINC_ENV_NAME="${HELIOLINC_ENV_NAME:-heliolinc}"
HELIOLINC_PYTHON_DEFAULT="${CONDA_HOME_DEFAULT}/envs/${HELIOLINC_ENV_NAME}/bin/python"
PYTHON_BIN="${PYTHON_BIN:-${HELIOLINC_PYTHON_DEFAULT}}"

# ---- path to your existing script (adjust if needed)
RUN_ONE="${PROJECT_ROOT}/heliolincrr_local/run_single_night.sh"

if [[ ! -x "${RUN_ONE}" ]]; then
  echo "[fatal] ${RUN_ONE} not found or not executable"
  exit 2
fi

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "[fatal] Python executable not found: ${PYTHON_BIN}"
  echo "        Set PYTHON_BIN explicitly or create conda env '${HELIOLINC_ENV_NAME}'."
  exit 1
fi

# ---- generate nights list (portable, no GNU date dependency)
NIGHTS=($("${PYTHON_BIN}" - <<PY
from datetime import datetime, timedelta
start = datetime.strptime("${START}", "%Y%m%d")
for i in range(${NDAYS}):
    print((start + timedelta(days=i)).strftime("%Y%m%d"))
PY
))

echo "============================================================"
echo "[info] Will run run_single_night.sh for these nights:"
printf "  %s\n" "${NIGHTS[@]}"
echo "============================================================"

for NIGHT in "${NIGHTS[@]}"; do
  echo "------------------------------------------------------------"
  echo "[run] NIGHT=${NIGHT}"
  echo "------------------------------------------------------------"

  "${RUN_ONE}" "${NIGHT}"

  echo "[done] NIGHT=${NIGHT}"
done

echo "============================================================"
echo "[all done] ${NDAYS} nights finished"
echo "============================================================"
