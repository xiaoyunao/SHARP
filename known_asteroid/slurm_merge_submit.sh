#!/usr/bin/env bash
#SBATCH -J ka_merge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 YYYYMMDD true|false" >&2
  exit 1
fi

NIGHT="$1"
SUBMIT_MPC="$2"

ROOT_DIR="${ROOT_DIR:-/processed1}"
SCRIPT_DIR="${SCRIPT_DIR:-/pipeline/xiaoyunao/known_asteroid}"
PYTHON="${PYTHON:-/home/smtpipeline/Softwares/miniconda3/bin/python}"

OBS_CODE="${OBS_CODE:-327}"
OBS_NAME="${OBS_NAME:-Xinglong Station}"
SUBMITTER_NAME="${SUBMITTER_NAME:-Y.-A. Xiao}"
INSTITUTION="${INSTITUTION:-NAOC}"
OBSERVERS="${OBSERVERS:-Xiangnan Guan;Pengfei Liu;Xiaoming Teng;WSGP Team}"
MEASURERS="${MEASURERS:-Niu Li;Yun-Ao Xiao;Jingyi Zhang;Hu Zou;WSGP Team}"
COMMENT="${COMMENT:-}"
TELESCOPE_NAME="${TELESCOPE_NAME:-60/90cm Schmidt telescope}"
APERTURE="${APERTURE:-0.6}"
DESIGN="${DESIGN:-Schmidt}"
DETECTOR="${DETECTOR:-CCD}"
F_RATIO="${F_RATIO:-3}"
FILTER_NAME="${FILTER_NAME:-unfiltered}"
ARRAY_SIZE="${ARRAY_SIZE:-9216 x 9232}"
PIXEL_SCALE="${PIXEL_SCALE:-1.15}"
COINVESTIGATORS="${COINVESTIGATORS:-Hu Zou}"
COLLABORATORS="${COLLABORATORS:-Yun-Ao Xiao;Niu Li;Jingyi Zhang;Wenxiong Li;Zhaobin Chen;Shufei Liu;WSGP Team}"
FUNDING_SOURCE="${FUNDING_SOURCE:-The 60/90-cm Schmidt telescope situated at the Xinglong station is overseen by the research team of the Wide-field Survey and Galaxy Physics (WSGP) at NAOC.}"
ASTCAT="${ASTCAT:-Gaia3E}"
BAND="${BAND:-G}"
PHOTCAT="${PHOTCAT:-Gaia3E}"
ERR_UNIT="${ERR_UNIT:-deg}"
MIN_OBSERVATIONS="${MIN_OBSERVATIONS:-3}"
PROG="${PROG:-}"
INCLUDE_LOGSNR="${INCLUDE_LOGSNR:-0}"
AC2_EMAIL="${AC2_EMAIL:-wsgp2024@163.com}"
ACK_PREFIX="${ACK_PREFIX:-Known-object ADES submission}"
OBJ_TYPE="${OBJ_TYPE:-MBA}"

L4_DIR="${ROOT_DIR}/${NIGHT}/L4"
PARTS_DIR="${L4_DIR}/known_asteroid_parts"
ALL_FITS="${L4_DIR}/${NIGHT}_all_asteroids.fits"
MATCHED_FITS="${L4_DIR}/${NIGHT}_matched_asteroids.fits"
OUT_PSV="${L4_DIR}/${NIGHT}_matched_asteroids_ades.psv"
REPLY_TXT="${L4_DIR}/${NIGHT}_mpc_reply.txt"

mkdir -p "${L4_DIR}" "${PARTS_DIR}"

if [[ ! -s "${ALL_FITS}" || ! -s "${MATCHED_FITS}" ]]; then
  "${PYTHON}" "${SCRIPT_DIR}/merge_night_parts.py" "${NIGHT}" --parts-dir "${PARTS_DIR}" --outdir "${L4_DIR}"
fi

if [[ "${SUBMIT_MPC}" != "true" ]]; then
  echo "[DONE] ${NIGHT} extraction finished (submit disabled)"
  exit 0
fi

if [[ -s "${OUT_PSV}" && -s "${REPLY_TXT}" ]]; then
  echo "[SKIP] ${NIGHT} report products already exist"
  exit 0
fi

if [[ ! -s "${MATCHED_FITS}" ]]; then
  echo "[FATAL] missing matched FITS: ${MATCHED_FITS}" >&2
  exit 2
fi

CMD=(
  "${PYTHON}" "${SCRIPT_DIR}/export_ades.py"
  --fits "${MATCHED_FITS}"
  --out "${OUT_PSV}"
  --obs-code "${OBS_CODE}"
  --obs-name "${OBS_NAME}"
  --submitter-name "${SUBMITTER_NAME}"
  --institution "${INSTITUTION}"
  --observers "${OBSERVERS}"
  --measurers "${MEASURERS}"
  --comment "${COMMENT}"
  --telescope-name "${TELESCOPE_NAME}"
  --aperture "${APERTURE}"
  --design "${DESIGN}"
  --detector "${DETECTOR}"
  --f-ratio "${F_RATIO}"
  --filter-name "${FILTER_NAME}"
  --array-size "${ARRAY_SIZE}"
  --pixel-scale "${PIXEL_SCALE}"
  --coinvestigators "${COINVESTIGATORS}"
  --collaborators "${COLLABORATORS}"
  --funding-source "${FUNDING_SOURCE}"
  --astcat "${ASTCAT}"
  --band "${BAND}"
  --photcat "${PHOTCAT}"
  --err-unit "${ERR_UNIT}"
  --min-observations "${MIN_OBSERVATIONS}"
  --ac2-email "${AC2_EMAIL}"
  --ack "${ACK_PREFIX}"
  --obj-type "${OBJ_TYPE}"
  --response-out "${REPLY_TXT}"
  --submit
)

if [[ -n "${PROG}" ]]; then
  CMD+=(--prog "${PROG}")
fi
if [[ "${INCLUDE_LOGSNR}" == "1" ]]; then
  CMD+=(--include-logsnr)
fi

"${CMD[@]}"
