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

做 MPC test validation：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python submit_reviewed_unknown.py 20260512 --validate
```

正式上报必须显式打开：

```bash
/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python submit_reviewed_unknown.py 20260512 --validate --submit
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
