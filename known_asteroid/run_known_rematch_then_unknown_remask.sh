#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  cat >&2 <<'EOF'
Usage:
  run_known_rematch_then_unknown_remask.sh KNOWN_START KNOWN_END UNKNOWN_START UNKNOWN_END

Environment:
  RUN_ID                  default known_rematch_<timestamp>
  MAX_PARALLEL            default 16
  FORCE_EXTRACT           default 1
  MASK_SEP_ARCSEC         default 1.5
  MASK_MATCHED_SUFFIX     default _mask15
  REMASK_MEM              default 12G
  POLL_SECONDS            default 300
  SBATCH_RETRIES          default 12
  SBATCH_RETRY_SLEEP      default 60
EOF
  exit 1
fi

KNOWN_START="$1"
KNOWN_END="$2"
UNKNOWN_START="$3"
UNKNOWN_END="$4"

RUN_ID="${RUN_ID:-known_rematch_$(date +%Y%m%d_%H%M%S)}"
MAX_PARALLEL="${MAX_PARALLEL:-16}"
FORCE_EXTRACT="${FORCE_EXTRACT:-1}"
MASK_SEP_ARCSEC="${MASK_SEP_ARCSEC:-1.5}"
MASK_MATCHED_SUFFIX="${MASK_MATCHED_SUFFIX:-_mask15}"
REMASK_MEM="${REMASK_MEM:-12G}"
POLL_SECONDS="${POLL_SECONDS:-300}"
SBATCH_RETRIES="${SBATCH_RETRIES:-12}"
SBATCH_RETRY_SLEEP="${SBATCH_RETRY_SLEEP:-60}"

SCRIPT_DIR="${SCRIPT_DIR:-/pipeline/xiaoyunao/known_asteroid}"
HELIOLINCRR_DIR="${HELIOLINCRR_DIR:-/pipeline/xiaoyunao/heliolincrr}"
HELIOLINC_PYTHON="${HELIOLINC_PYTHON:-/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python}"
KNOWN_LOG="${KNOWN_LOG:-/pipeline/xiaoyunao/known_asteroid/runtime/logs/${RUN_ID}.log}"
DRIVER_LOG="${DRIVER_LOG:-/pipeline/xiaoyunao/known_asteroid/runtime/logs/${RUN_ID}_driver.log}"
REMASK_LOG="${REMASK_LOG:-/pipeline/xiaoyunao/data/heliolincrr/batch_logs/${RUN_ID}_unknown_remask.log}"
REMASK_STATUS="${REMASK_STATUS:-/pipeline/xiaoyunao/data/heliolincrr/batch_logs/${RUN_ID}_unknown_remask_status.tsv}"

mkdir -p "$(dirname "${KNOWN_LOG}")" "$(dirname "${REMASK_LOG}")"

submit_sbatch() {
  local attempt=1
  local rc=0
  while true; do
    set +e
    local output
    output="$(sbatch --parsable "$@" 2>&1)"
    rc=$?
    set -e
    if [[ "${rc}" -eq 0 ]]; then
      printf '%s\n' "${output}"
      return 0
    fi
    if (( attempt >= SBATCH_RETRIES )); then
      printf '%s\n' "${output}" >&2
      return "${rc}"
    fi
    printf '[driver] sbatch failed attempt %d/%d: %s\n' "${attempt}" "${SBATCH_RETRIES}" "${output}" >> "${DRIVER_LOG}"
    sleep "${SBATCH_RETRY_SLEEP}"
    attempt=$((attempt + 1))
  done
}

{
  echo "[driver] start $(date -Is)"
  echo "[driver] run_id=${RUN_ID}"
  echo "[driver] known=${KNOWN_START}..${KNOWN_END} unknown=${UNKNOWN_START}..${UNKNOWN_END}"
} >> "${DRIVER_LOG}"

cd "${SCRIPT_DIR}"
FORCE_EXTRACT="${FORCE_EXTRACT}" \
MASK_SEP_ARCSEC="${MASK_SEP_ARCSEC}" \
MASK_MATCHED_SUFFIX="${MASK_MATCHED_SUFFIX}" \
MAX_PARALLEL="${MAX_PARALLEL}" \
./submit_pipeline_slurm.sh --batch true --submit-mpc false --max-parallel "${MAX_PARALLEL}" "${KNOWN_START}" "${KNOWN_END}" >> "${KNOWN_LOG}" 2>&1

DEP_IDS="$(awk '/ finalize job /{print $5}' "${KNOWN_LOG}" | paste -sd: -)"
FINALIZE_JOB_COUNT="$(awk '/ finalize job /{n++} END{print n+0}' "${KNOWN_LOG}")"
if [[ -z "${DEP_IDS}" ]]; then
  echo "[driver] fatal: no finalize jobs found in ${KNOWN_LOG}" >> "${DRIVER_LOG}"
  exit 2
fi

while true; do
  ACTIVE_FINALIZE_COUNT="$(
    awk '/ finalize job /{print $5}' "${KNOWN_LOG}" |
    sort -u |
    comm -12 - <(squeue -h -u "${USER}" -o "%i" | sort -u) |
    wc -l |
    tr -d ' '
  )"
  echo "[driver] active_finalize_count=${ACTIVE_FINALIZE_COUNT} at $(date -Is)" >> "${DRIVER_LOG}"
  if [[ "${ACTIVE_FINALIZE_COUNT}" -eq 0 ]]; then
    break
  fi
  sleep "${POLL_SECONDS}"
done

REMASK_JOB="$(
  submit_sbatch \
    --job-name=unknown_remask \
    --cpus-per-task=1 \
    --mem="${REMASK_MEM}" \
    --output="${REMASK_LOG}" \
    --error="${REMASK_LOG}" \
    --wrap="cd ${HELIOLINCRR_DIR} && ${HELIOLINC_PYTHON} remask_unknown_with_known.py ${UNKNOWN_START} ${UNKNOWN_END} --status-out ${REMASK_STATUS}"
)"

{
  echo "RUN_ID=${RUN_ID}"
  echo "KNOWN_LOG=${KNOWN_LOG}"
  echo "FINALIZE_JOB_COUNT=${FINALIZE_JOB_COUNT}"
  echo "REMASK_JOB=${REMASK_JOB}"
  echo "REMASK_LOG=${REMASK_LOG}"
  echo "REMASK_STATUS=${REMASK_STATUS}"
  echo "[driver] done $(date -Is)"
} >> "${DRIVER_LOG}"
