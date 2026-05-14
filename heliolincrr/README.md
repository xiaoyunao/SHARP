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

MPC ADES 的 unknown 临时标识使用 `trkSub`。当前代码只接受 1-8 位
ASCII 字母、数字、下划线或连字符。

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

人工复核接口预留为两列 CSV：

```text
tracklet_id,is_real
XA000001,1
XA000002,0
```

其中 `tracklet_id` 对应写入 unknown catalog 的 `trk_sub`，`is_real=1`
才允许上报，`0` 会被过滤。

单独导出 unknown ADES PSV：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python export_unknown_ades.py 20260220 \
  --processed-root /processed1 \
  --catalog /processed1/20260220/L4/20260220_unknown_links.json \
  --review-csv /path/to/review.csv \
  --require-review \
  --out /processed1/20260220/L4/20260220_unknown_links_ades.psv
```

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
