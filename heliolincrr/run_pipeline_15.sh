#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 YYYYMMDD"
  exit 1
fi

NIGHT="$1"
CONDA_HOME_DEFAULT="${CONDA_HOME_DEFAULT:-/home/smtpipeline/Softwares/miniconda3}"
HELIOLINC_ENV_NAME="${HELIOLINC_ENV_NAME:-heliolinc}"
HELIOLINC_PYTHON_DEFAULT="${CONDA_HOME_DEFAULT}/envs/${HELIOLINC_ENV_NAME}/bin/python"
PYTHON_BIN="${PYTHON_BIN:-${HELIOLINC_PYTHON_DEFAULT}}"
PREP_W15_MISSING_NIGHTS="${PREP_W15_MISSING_NIGHTS:-1}"

ROOT_OUT="/pipeline/xiaoyunao/data/heliolincrr"
TRACKLET_DIR="${ROOT_OUT}/${NIGHT}/tracklets_linreproj"
NIGHTLY_ALL="${TRACKLET_DIR}/tracklets_${NIGHT}_ALL.fits"
W15_OUTDIR="${ROOT_OUT}/${NIGHT}/tracklets_15"
RR_W15_DIR="${ROOT_OUT}/${NIGHT}/rr_links_15"
W15_NDAYS=15

RR_W15_CORES=20
RR_W15_REF_MODE="mid"
RR_W15_REF_DT=0.50
RR_W15_TOL=0.02
RR_W15_MIN_LEN_OBS=4
RR_W15_MIN_NIGHTS=2
RR_W15_KCAP=200
RR_W15_MAXV=200
RR_W15_MININIT=0.02

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "[fatal] Python executable not found: ${PYTHON_BIN}"
  exit 1
fi

is_valid_table() {
  local fits_path="$1"
  [[ -f "${fits_path}" ]] || return 1
  "${PYTHON_BIN}" - <<PY >/dev/null 2>&1
from astropy.table import Table
Table.read("${fits_path}")
PY
}

prepare_previous_nights_for_15() {
  local end_night="$1"
  local nights
  mapfile -t nights < <("${PYTHON_BIN}" - <<PY
from datetime import datetime, timedelta
end = datetime.strptime("${end_night}", "%Y%m%d")
start = end - timedelta(days=${W15_NDAYS} - 1)
d = start
while d < end:
    print(d.strftime("%Y%m%d"))
    d += timedelta(days=1)
PY
)
  local n
  for n in "${nights[@]}"; do
    if ! is_valid_table "${ROOT_OUT}/${n}/tracklets_linreproj/tracklets_${n}_ALL.fits"; then
      echo "[fatal] missing nightly ALL required for 15-night merge: ${n}"
      exit 1
    fi
  done
}

if ! is_valid_table "${NIGHTLY_ALL}"; then
  echo "[fatal] nightly ALL not found: ${NIGHTLY_ALL}"
  exit 1
fi

if [[ "${PREP_W15_MISSING_NIGHTS}" -eq 1 ]]; then
  prepare_previous_nights_for_15 "${NIGHT}"
fi

mkdir -p "${W15_OUTDIR}" "${RR_W15_DIR}"

"${PYTHON_BIN}" merge_tracklets_15.py "${NIGHT}" \
  --root "${ROOT_OUT}" \
  --ndays "${W15_NDAYS}" \
  --nightly-subpath "tracklets_linreproj" \
  --nightly-name-template "tracklets_{night}_ALL.fits" \
  --outdir "${W15_OUTDIR}" \
  --overwrite \
  --write-index

shopt -s nullglob
W15_CANDIDATES=(${W15_OUTDIR}/tracklets_*_${NIGHT}_W${W15_NDAYS}.fits)
shopt -u nullglob
if [[ ${#W15_CANDIDATES[@]} -eq 0 ]]; then
  echo "[fatal] no merged 15-night tracklets found in ${W15_OUTDIR}"
  exit 1
fi
W15_FILE="${W15_CANDIDATES[-1]}"

"${PYTHON_BIN}" run_rr_from_tracklets_15.py \
  --infile "${W15_FILE}" \
  --outdir "${RR_W15_DIR}" \
  --profile w15 \
  --cores "${RR_W15_CORES}" \
  --ref-epoch-mode "${RR_W15_REF_MODE}" \
  --ref-dt-days "${RR_W15_REF_DT}" \
  --tol "${RR_W15_TOL}" \
  --min-len-obs "${RR_W15_MIN_LEN_OBS}" \
  --min-nights "${RR_W15_MIN_NIGHTS}" \
  --k-neighbors-cap "${RR_W15_KCAP}" \
  --max-v-kms "${RR_W15_MAXV}" \
  --min-init-earth-au "${RR_W15_MININIT}"

"${PYTHON_BIN}" orbit_confirm_links_15.py \
  --rr-dir "${RR_W15_DIR}" \
  --tracklets "${W15_FILE}" \
  --cores 16 \
  --log-every 500

echo "[done] 15-night pipeline finished for ${NIGHT}"
