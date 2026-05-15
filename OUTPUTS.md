# Runtime Outputs

本文件记录服务器正式流程的主要运行产物位置。仓库代码仍以服务器路径为准。

## Nightly L4 Products

每个观测夜的 L4 目录：

```text
/processed1/<night>/L4/
```

known asteroid 主产物：

```text
<night>_all_asteroids.fits
<night>_matched_asteroids.fits
<night>_file_manifest.txt
<night>_match.log
known_asteroid_parts/
  *_all_asteroids.fits
  *_matched_asteroids.fits
```

如果 known MPC 上报打开，还会生成：

```text
<night>_matched_asteroids_ades.psv
<night>_mpc_reply.txt
```

unknown 主产物：

```text
<night>_unknown_links.json
<night>_unknown_links.fits
<night>_unknown_links_ades.psv
<night>_unknown_mpc_validate_reply.txt
<night>_unknown_mpc_reply.txt
```

其中 `<night>_unknown_mpc_reply.txt` 只有真实 submit 后才会存在。

## Unknown trkSub History

unknown 临时编号使用全局 history，供单夜和未来 15 夜流程共用：

```text
/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl
/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl.lock
```

`trkSub` 规则：

```text
alphabet = 0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ
first    = 00000001
width    = 8
```

history 每行一个 JSON record，记录 `trk_sub`、night、mode、linkage_id、
内部 tracklet ids、image/object ids、MJD、RA/Dec、orbit-fit 摘要和源文件路径。

## heliolincrr Products

单夜 unknown 中间结果：

```text
/pipeline/xiaoyunao/data/heliolincrr/<night>/mask_gaia/
/pipeline/xiaoyunao/data/heliolincrr/<night>/tracklets_linreproj/
/pipeline/xiaoyunao/data/heliolincrr/<night>/rr_links/
/pipeline/xiaoyunao/data/heliolincrr/<night>/rr_links/orbit_confirm/
/pipeline/xiaoyunao/data/heliolincrr/<night>/analysis/
```

summary 文件：

```text
/pipeline/xiaoyunao/data/heliolincrr/<night>/analysis/<night>_single_night_summary.json
/pipeline/xiaoyunao/data/heliolincrr/<night>/analysis/<night>_single_night_summary.txt
```

unknown GIF：

```text
/pipeline/xiaoyunao/heliolincrr/plots/<night>/
  unknown_link_XXXX_<night>.gif
  <night>_unknown_link_summary.json
```

unknown 人工复核包：

```text
/pipeline/xiaoyunao/heliolincrr/review_packages/<night>/
  gifs/
    <trkSub>_linkXXXX_<night>.gif
  <night>_unknown_review.csv
  <night>_unknown_review_full.fits
  <night>_unknown_review_ades.fits
  <night>_unknown_review_manifest.json
/pipeline/xiaoyunao/heliolincrr/review_packages/<night>_unknown_review.tar.gz
```

其中 `<night>_unknown_review.csv` 是人工标注模板，保留 `tracklet_id,is_real`
两列；`*_unknown_review_full.fits` 是逐 detection 复核表，包含
`RA_Win/DEC_Win`、`RA_PSF/DEC_PSF`、`trk_sub`、`linkage_id`、tracklet
来源和轨道摘要；`*_unknown_review_ades.fits` 是复核包内用于核对 ADES 行的表。

## known_asteroid Plots

known asteroid nightly visual products are not written under L4. They are stored in:

```text
/pipeline/xiaoyunao/known_asteroid/plots/<night>/
```

Typical files:

```text
survey_coverage_<night>.png
known_asteroid_allsky_<night>.png
known_asteroid_statistics_<night>.png
known_asteroid_field_yield_<night>.png
telescope_slew_statistics_<night>.png
topXX_<object>_<night>.gif
```

## Survey Scheduler Products

Survey scheduler plans:

```text
/pipeline/xiaoyunao/survey/runtime/plans/<night>_plan.json
```

Survey persistent history:

```text
/pipeline/xiaoyunao/survey/runtime/history/exposure_history.fits
/pipeline/xiaoyunao/survey/runtime/history/l2_ingest_state.json
```

Survey plots and nightly cycle GIF:

```text
/pipeline/xiaoyunao/survey/plots/<night>/
  <night>_all_fields.png
  <night>_cycle01.png
  <night>_cycle02.png
  <night>_cycle03.png
  <night>_cycle04.png
  <night>_cycle05.png
  <night>_nightly_cycles.gif
  <night>_plots.json
```

## Notes

- `known_asteroid/plots/<night>/survey_coverage_<night>.png` 是 known 可视化脚本对 scheduler coverage 的汇总展示。
- `survey/plots/<night>/` 是 survey 自己的调度/循环图。
- unknown 上报默认应保持显式触发；daily 自动化后续不应默认真实 submit。
