#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-/home/smtpipeline/Softwares/miniconda3/bin/python}"
PROCESSED_ROOT="${PROCESSED_ROOT:-/processed1}"
PLOT_ROOT="${PLOT_ROOT:-${ROOT_DIR}/plots}"

RUN_DATE="${1:-$(TZ=Asia/Shanghai date +%F)}"
TARGET_NIGHT="${TARGET_NIGHT:-$(TZ=Asia/Shanghai date -d "${RUN_DATE} -1 day" +%Y%m%d)}"
NIGHT_DIR="${PROCESSED_ROOT}/${TARGET_NIGHT}"
L2_DIR="${NIGHT_DIR}/L2"
MATCHED_FITS="${PROCESSED_ROOT}/${TARGET_NIGHT}/L4/${TARGET_NIGHT}_matched_asteroids.fits"
ALL_FITS="${PROCESSED_ROOT}/${TARGET_NIGHT}/L4/${TARGET_NIGHT}_all_asteroids.fits"
OUT_DIR="${PLOT_ROOT}/${TARGET_NIGHT}"

if [[ ! -d "${NIGHT_DIR}" || ! -d "${L2_DIR}" ]]; then
  exit 0
fi

MP_COUNT="$(
  find "${L2_DIR}" -maxdepth 1 -type f -name 'OBJ_MP_*.fits.gz' | wc -l | tr -d ' '
)"
if [[ "${MP_COUNT}" == "0" ]]; then
  exit 0
fi

if [[ ! -s "${ALL_FITS}" ]]; then
  # known_asteroid has not finished merging yet
  exit 0
fi

if [[ ! -s "${MATCHED_FITS}" ]]; then
  # the night has data, but no matched asteroids to visualize
  exit 0
fi

"${PYTHON}" "${ROOT_DIR}/update_all_matched_history.py" --processed-root "${PROCESSED_ROOT}" --night "${TARGET_NIGHT}"

if [[ -d "${OUT_DIR}" ]] && compgen -G "${OUT_DIR}/*.png" > /dev/null && compgen -G "${OUT_DIR}/*.gif" > /dev/null; then
  exit 0
fi

"${PYTHON}" "${ROOT_DIR}/plot_known_asteroids.py" "${TARGET_NIGHT}" --processed-root "${PROCESSED_ROOT}" --plot-root "${PLOT_ROOT}"
