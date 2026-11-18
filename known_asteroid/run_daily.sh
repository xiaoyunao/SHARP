#!/usr/bin/env bash
set -euo pipefail
BASE_DIR="${SMT_BASE_DIR:-/pipeline/xiaoyunao}"
cd "${BASE_DIR}"
export PYTHONPATH="${BASE_DIR}${PYTHONPATH:+:${PYTHONPATH}}"
RUN_DATE="${1:-$(TZ=Asia/Shanghai date +%F)}"
NIGHT="${TARGET_NIGHT:-$(TZ=Asia/Shanghai date -d "${RUN_DATE} -1 day" +%Y%m%d)}"
exec /home/smtpipeline/Softwares/miniconda3/bin/python -u -m targeted_observing.server_pipeline --night "${NIGHT}"
