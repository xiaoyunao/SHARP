# heliolincrr

本目录保存服务器当前 `heliolincrr` 代码，并作为仓库中的唯一正式版本。

说明：

- 保留服务器路径和环境假设
- 默认使用 `heliolinc` 环境
- 后续调试和调参都直接基于这一套代码进行
- 单夜正式入口为 `run_single_night.sh`
- 15 夜正式入口为 `run_pipeline_15.sh`
- 单夜统计脚本为 `summarize_single_night.py`
- 未知单夜 link 的可视化已并入 `run_single_night.sh`，底层脚本为 `plot_unknown_links.py`，输出到 `heliolincrr/plots/YYYYMMDD/`
- unknown ADES 导出脚本为 `export_unknown_ades.py`，默认只写 PSV，不自动提交 MPC

## Unknown trkSub and review interface

MPC ADES 的 unknown 临时标识使用 `trkSub`。本项目自动分配的 `trkSub`
采用 7 位 62 进制，满足 MPC 当前要求的 no more than 7 characters：

- 字符表：`0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ`
- `a=10`，`z=35`，`A=36`
- 起始值：`0000001`
- 默认历史文件：`/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl`
- 历史文件旁边会使用 `.lock` 文件加锁，避免单夜和 15 夜流程同时分配时撞号

历史文件一行一个 JSON 记录，包含 `trk_sub`、night、mode、linkage_id、
内部 tracklet ids、原始 image/object ids、MJD、RA/Dec 和轨道拟合摘要，便于从
上报编号追溯到源观测点、tracklet 和 link。

在单夜 summary 阶段可传入编号映射：

```bash
TRK_SUB_MAP=/path/to/trk_sub_map.csv bash run_single_night.sh 20260220
```

映射文件支持 CSV/TSV/JSON。CSV 推荐列：

```text
linkage_id,trk_sub
1,XA000001
3,XA000002
```

也可以用内部 `tracklet_id` 映射：

```text
tracklet_id,trk_sub
T12345,XA000001
```

正常单夜流程默认会在 summary 后自动调用 `assign_unknown_trksub.py`，给还没有
`trk_sub` 的 unknown rows 分配全局唯一编号，并回写：

- `/processed1/<night>/L4/<night>_unknown_links.json`
- `/processed1/<night>/L4/<night>_unknown_links.fits`

如需临时关闭：

```bash
ASSIGN_UNKNOWN_TRKSUB=0 bash run_single_night.sh 20260220
```

人工复核接口预留为两列 CSV：

```text
tracklet_id,is_real
XA000001,1
XA000002,0
```

其中 `tracklet_id` 对应写入 unknown catalog 的 `trk_sub`，`is_real=1`
才允许上报，`0` 会被过滤。

给观测助手打包某晚 unknown GIF，并生成复核 CSV 模板：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python package_unknown_review.py 20260220 --make-tar
```

默认输出：

```text
/pipeline/xiaoyunao/heliolincrr/review_packages/<night>/
  gifs/
    <trkSub>_linkXXXX_<night>.gif
  <night>_unknown_review.csv
  <night>_unknown_review_manifest.json
/pipeline/xiaoyunao/heliolincrr/review_packages/<night>_unknown_review.tar.gz
```

观测助手只需要填写 `<night>_unknown_review.csv` 的 `is_real` 列：

- `1`: 真源，允许进入 unknown ADES
- `0`: 假源，导出和上报时过滤

网页人工 check 的最终提交文件约定为：

```text
/pipeline/xiaoyunao/heliolincrr/review_packages/<night>/<night>_submit.csv
```

格式仍为两列：

```text
tracklet_id,is_real
000001fU,1
000001fV,0
```

`<night>_submit.csv` 表示该夜已经完成最终判定。程序会要求每个 unknown
`tracklet_id` 都有明确的 `0/1`，然后只保留 `is_real=1` 的条目生成 reviewed/masked
产物和 ADES PSV。默认输出：

```text
/processed1/<night>/L4/<night>_unknown_links_submit_masked.json
/processed1/<night>/L4/<night>_unknown_links_submit_masked.fits
/processed1/<night>/L4/<night>_unknown_links_submit_ades.psv
/processed1/<night>/L4/<night>_unknown_links_submit_stats.json
```

单独导出 unknown ADES PSV：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python export_unknown_ades.py 20260220 \
  --processed-root /processed1 \
  --catalog /processed1/20260220/L4/20260220_unknown_links.json \
  --review-csv /path/to/review.csv \
  --require-review \
  --out /processed1/20260220/L4/20260220_unknown_links_ades.psv
```

网页 submit CSV 完成后，先生成 reviewed/masked 产物和 PSV，不 validate、不上报：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python submit_reviewed_unknown.py 20260512
```

`submit_reviewed_unknown.py` 默认在 ADES PSV 中包含 `logSNR` 列，数值来自
`log10(Flux_Aper4 / FluxErr_Aper4)`。如需临时关闭：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python submit_reviewed_unknown.py 20260512 --no-logsnr
```

做 MPC test validation：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python submit_reviewed_unknown.py 20260512 --validate
```

正式上报必须显式打开：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python submit_reviewed_unknown.py 20260512 --validate --submit
```

### Daily automation

每日 09:00 总入口：

```bash
/pipeline/xiaoyunao/heliolincrr/run_daily_pipeline.sh
```

它按顺序运行：

1. `/pipeline/xiaoyunao/survey/run_daily.sh` 只为 `RUN_DATE` 生成当天晚上的观测脚本
2. 按 `RECOVERY_LOOKBACK_DAYS` 回看实际存在数据的夜次，调用 `/pipeline/xiaoyunao/known_asteroid/run_daily.sh`
3. 对同一批夜次调用 `/pipeline/xiaoyunao/heliolincrr/run_daily_unknown.sh`，等待 known 的 1.5 角秒 mask ready 后生成 unknown check 包

`run_daily_unknown.sh` 默认要求 `/processed1/<night>/L4/<night>_matched_asteroids_mask15.fits`。
如果 known 已完成但 1.5 角秒内没有任何 matched detection，则依赖
`<night>_known_asteroid_status.json` 明确记录 empty mask 后继续。

`run_daily_pipeline.sh` 的断电恢复逻辑不补生成连续关机期间已经错过的中间观测脚本；
survey 只针对当前 `RUN_DATE`。known/unknown 数据处理会按回看窗口补实际拍到且尚未完成的夜。

`run_daily_unknown.sh` 只有在 check 包已经生成且 manifest 中 `n_catalog_rows > 0`
时，才为该夜启动 `watch_submit_reviews.py`。unknown 为 `0` 的夜只保留空 check 包，
不挂 submit watcher。

建议 cron：

```cron
0 9 * * * cd /pipeline/xiaoyunao && /bin/bash /pipeline/xiaoyunao/heliolincrr/run_daily_pipeline.sh
@reboot sleep 600 && cd /pipeline/xiaoyunao && /bin/bash /pipeline/xiaoyunao/heliolincrr/run_daily_pipeline.sh
```

`@reboot` 会重新跑当前日期的总入口；survey 只补当天晚上的脚本，known/unknown
按回看窗口补处理数据。

网页生成 `<night>_submit.csv` 后，批量监测并提交。历史补报推荐用 wrapper：

```bash
/pipeline/xiaoyunao/heliolincrr/run_review_submit_backlog_watch.sh 20251116 20260617
```

该 wrapper 会后台运行 `watch_submit_reviews.py --follow --exit-when-complete --validate --submit --retry-failed`，
并把状态写到：

```text
/pipeline/xiaoyunao/data/heliolincrr/review_submit_backlog_20251116_20260617.json
```

底层命令也可手动运行：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python /pipeline/xiaoyunao/heliolincrr/watch_submit_reviews.py \
  --start 20251116 \
  --end 20260617 \
  --validate \
  --submit \
  --follow \
  --exit-when-complete \
  --retry-failed
```

默认正式运行状态记录在：

```text
/pipeline/xiaoyunao/data/heliolincrr/review_submit_state.json
```

watcher 只处理已有 `*_unknown_review_manifest.json` 的夜次；zero unknown check 包会直接记为
`no_observations`。有 unknown 但 submit CSV 全为 `0` 的夜也会记为 `no_observations`，
不会向 MPC 上报。

在 `run_single_night.sh` 里启用导出：

```bash
EXPORT_UNKNOWN_ADES=1 \
TRK_SUB_MAP=/path/to/trk_sub_map.csv \
UNKNOWN_REVIEW_CSV=/path/to/review.csv \
REQUIRE_UNKNOWN_REVIEW=1 \
bash run_single_night.sh 20260220
```

真实提交仍需显式打开：

```bash
SUBMIT_UNKNOWN_MPC=1 EXPORT_UNKNOWN_ADES=1 bash run_single_night.sh 20260220
```
