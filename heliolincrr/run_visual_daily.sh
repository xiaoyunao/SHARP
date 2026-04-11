#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python}"
PROCESSED_ROOT="${PROCESSED_ROOT:-/processed1}"
ROOT_OUT="${ROOT_OUT:-/pipeline/xiaoyunao/data/heliolincrr}"
PLOT_ROOT="${PLOT_ROOT:-${ROOT_DIR}/plots}"

RUN_DATE="${1:-$(TZ=Asia/Shanghai date +%F)}"
TARGET_NIGHT="${TARGET_NIGHT:-$(TZ=Asia/Shanghai date -d "${RUN_DATE} -1 day" +%Y%m%d)}"
RR_DIR="${ROOT_OUT}/${TARGET_NIGHT}/rr_links"
UNKNOWN_JSON="${PROCESSED_ROOT}/${TARGET_NIGHT}/L4/${TARGET_NIGHT}_unknown_links.json"
L1_DIR="${PROCESSED_ROOT}/${TARGET_NIGHT}/L1"
OUT_DIR="${PLOT_ROOT}/${TARGET_NIGHT}"

if [[ ! -d "${L1_DIR}" || ! -s "${UNKNOWN_JSON}" || ! -d "${RR_DIR}" ]]; then
  exit 0
fi

mkdir -p "${OUT_DIR}"

"${PYTHON_BIN}" "${ROOT_DIR}/plot_unknown_links.py" "${TARGET_NIGHT}" \
  --processed-root "${PROCESSED_ROOT}" \
  --root-out "${ROOT_OUT}" \
  --plot-root "${PLOT_ROOT}" \
  --catalog "${UNKNOWN_JSON}"
