# survey

服务器部署版每日调度程序。

## 行为

- 每天上午 `09:00` 运行一次。
- 当天运行日期记为当夜日期，例如 `2026-02-20`。
- footprint 默认路径：`/pipeline/xiaoyunao/survey/footprints/equatorial_footprints.fits`
- 历史曝光文件默认保存在 `workspace/history/exposure_history.fits`
- 如果历史文件不存在：
  - 直接新建空历史文件
  - 不回读 `/processed1`
  - 直接生成当夜观测列表
- 如果历史文件存在：
  - 尝试读取 `/processed1/YYYYMMDD/L2`
  - 只统计文件名里包含 `MP` 的文件
  - 从文件名中提取 `MP` 后面的编号，例如 `OBJ_MP_0023_0140_cat.fits.gz -> MP0023`
  - 每出现一次，对历史曝光次数加 `1`
  - 同一夜重复运行时，不会重复累计已经 ingest 过的文件
  - 如果目录不存在，或者没有 `L2`，或者没有任何 `MP` 文件，则跳过历史更新，直接生成当夜列表
- 如果设置了共享目录，还会额外复制：
  - `YYYYMMDD_plan.txt`
  - `current_plan.txt`

## 输出

- `workspace/plans/YYYYMMDD_plan.json`
- `workspace/plans/YYYYMMDD_plan.txt`
- `workspace/history/exposure_history.fits`
- `workspace/history/l2_ingest_state.json`
- `survey/logs/YYYYMMDD.log`
- `survey/plots/YYYYMMDD/`
  - `YYYYMMDD_cycle01.png`, `YYYYMMDD_cycle02.png`, ...
  - `YYYYMMDD_nightly_cycles.gif`
  - `YYYYMMDD_all_fields.png`
  - `YYYYMMDD_plots.json`

## 手动运行

```bash
python -m survey.run_daily \
  --date 2026-02-20 \
  --root-dir /pipeline/xiaoyunao/survey \
  --workspace /pipeline/xiaoyunao/survey/runtime \
  --footprints /pipeline/xiaoyunao/survey/footprints/equatorial_footprints.fits \
  --processed-root /processed1 \
  --publish-dir /path/to/shared/plan_dir
```

或者：

```bash
/bin/bash /pipeline/xiaoyunao/survey/run_daily.sh 2026-02-20
```

## 定时运行

参考 `survey/cron.example`：

```cron
0 9 * * * cd /pipeline/xiaoyunao && /bin/bash /pipeline/xiaoyunao/survey/run_daily.sh
```

## 部署说明

- 现在 `survey` 目录已经自包含，不再依赖其他 Python 包目录。
- 部署到服务器时，把 `equatorial_footprints.fits` 一起放到 `/pipeline/xiaoyunao/survey/footprints/`。
- 默认运行路径：
  - 包目录：`/pipeline/xiaoyunao/survey`
  - runtime：`/pipeline/xiaoyunao/survey/runtime`
  - processed：`/processed1`
  - 日志：`/pipeline/xiaoyunao/survey/logs/YYYYMMDD.log`
  - 图像：`/pipeline/xiaoyunao/survey/plots/YYYYMMDD/`
