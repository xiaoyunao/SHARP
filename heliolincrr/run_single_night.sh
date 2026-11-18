#!/usr/bin/env bash
set -euo pipefail
BASE_DIR="${SMT_BASE_DIR:-/pipeline/xiaoyunao}"
cd "${BASE_DIR}"
export PYTHONPATH="${BASE_DIR}${PYTHONPATH:+:${PYTHONPATH}}"
exec /home/smtpipeline/Softwares/miniconda3/bin/python -u -m targeted_observing.server_pipeline --night "${1:?YYYYMMDD required}"
