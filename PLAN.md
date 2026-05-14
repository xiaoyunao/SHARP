# PLAN

## Current objective

回到 `heliolincrr` 主线，先把单夜 unknown 搜索/提取自动化做稳。`unknown_backfill_continue_20260514_115014.log` 对应的服务器续跑任务已完成；本轮已修复 0-link 空表链路，并给单个 dense group 增加 tracklet 数量上限。

最新覆盖状态：

- `20251115..20260514` 有 MP L2 数据的夜次共 `122` 个
- 当前 `120/122` 个完成 compute 链路
- 未完成：
  - `20251201`: 无重复视场，`make_tracklet` 得到 `n_groups=0`，没有 `tracklets_ALL`
  - `20260414`: 缺 known matched FITS，已有 unknown summary 不可信，需先补 known_asteroid match
- `20260503` 已按新 dense-group 保护重跑，`unknown=95`，GIF/review package 已重建完成

短期目标：

- 自动选择一个可处理夜次
- 调用 `run_single_night.sh` 完成单夜 unknown 搜索
- 稳定产出 `/processed1/<night>/L4/<night>_unknown_links.{fits,json}`
- 给 unknown link 写入符合 MPC ADES `trkSub` 规则的临时编号
- 生成便于人工复核的单夜 summary/候选清单/日志
- 预留人工复核 mask：`tracklet_id,is_real`
- 上报链路默认关闭；导出/提交必须显式打开
- 异常 dense group 默认跳过，避免单个 group 生成数十万 tracklet 并污染 unknown catalog

## Milestones

1. `20260220` 已作为单夜 unknown 主测试夜完整跑通
2. `run_single_night.sh` 已串起 Gaia mask、tracklet、linear link、orbit confirm、summary、unknown GIF
3. `summarize_single_night.py` 已生成 unknown fit_ok catalog，并按 known asteroid match 结果扣除已知小行星
4. 历史主批处理曾推进到 `20260410`，当前 `20251116..20260410` 中有 `98` 个可读 nightly summary
5. 现有 summary 统计：平均 `64.63` 条 unknown fit_ok/night；排除异常夜后平均 `53.17` 条/night
6. 已新增 `--trk-sub-map`、人工复核 CSV 过滤接口和 `export_unknown_ades.py` unknown ADES 导出器
7. 已实现全局 `trkSub` 历史分配：`/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl`
8. `20260220` unknown ADES 已通过 MPC test validation，并完成一次正式 submit
9. 0-link 夜次现在应能生成 empty RR/orbit/summary/unknown JSON/FITS，并让 assign/plot 正常 0 行退出
10. 真实 0-link 夜 `20260128` 和 `20260224` 已重跑成功，均生成 unknown=0 的 summary 和 empty unknown FITS/JSON
11. `make_tracklet_linreproj.py` 已新增 `--max-tracklets-per-group`，默认 `100000`；`20260503` 重跑确认 group `073/075` 被跳过，unknown 从 `1591` 降到 `95`
12. `20260503` 已重新生成 `95/95` 个 unknown GIF，并打包 review package，`n_gifs_missing=0`
13. `20251128`, `20260412`, `20260430` 已补齐空 unknown FITS/JSON

## Outstanding issues

- 服务器续跑任务 `1597630` 已完成；全量产物审计显示 `120/122` 个有 L2 夜次完成 compute 链路
- `20251116..20260410` 仍有 `48` 个日历夜没有 summary，需要区分缺原始数据、失败和未跑
- `20251201` 这类没有重复视场、没有任何 group tracklet 文件的夜次仍会失败，尚未转成成功空结果
- `20260414` 缺 known matched FITS，必须先补 known_asteroid match，再重跑 unknown 扣除
- `20260503` 旧异常运行已写入 `1591` 条 `trkSub` history；本次重跑未清理 history，真实 export/review 前需决定是否清理或过滤旧记录
- 当前单夜自动化还没有“每日选择目标夜 + 防重复 + 日志 + 产物检查”的外层入口
- GIF 可视化很慢，应在自动提取中默认可跳过或限量，避免拖慢主计算
- 15 夜流程尚未接入同一个 `trkSub` history
- 人工复核 GIF 打包和 review CSV 生成器已实现初版，仍需接入 daily wrapper
- `20260220` unknown 正式 submit 已完成，需后续跟进 MPC/WAMO ingest 结果
- unknown 真实提交策略仍需收紧；后续 daily wrapper 默认不应自动 submit

## Validation criteria

- 给定一个存在 `/processed1/<night>/L2` 和 known asteroid match 的夜次，自动入口能完整运行或明确跳过
- 成功夜必须存在，即使 unknown=0 也要写出 schema-only FITS/JSON：
  - `/pipeline/xiaoyunao/data/heliolincrr/<night>/analysis/<night>_single_night_summary.json`
  - `/pipeline/xiaoyunao/data/heliolincrr/<night>/analysis/<night>_single_night_summary.txt`
  - `/processed1/<night>/L4/<night>_unknown_links.json`
- summary 中 `counts.matched_detections_total > 0`，否则不能认为已完成已知小行星扣除
- 若启用 `trkSub`，每个 unknown link 的 `trk_sub` 必须符合 8 位 `[0-9a-zA-Z]`
- `trkSub` 分配必须写入全局 history，并且重复运行同一 catalog 不新增重复记录
- 若启用人工复核，`is_real=0` 的 `trk_sub` 不得进入 ADES PSV
- 自动入口日志需记录 target night、skip/run reason、exit code、核心产物路径和 unknown count

## Next recommended steps

1. 补跑 `20260414` known_asteroid match，确认生成 `/processed1/20260414/L4/20260414_matched_asteroids.fits`
2. 重跑 `20260414` heliolincrr summary/unknown，确保 `matched_detections_total > 0`
3. 单独处理 `20251201` no group/no tracklet files 场景，决定是否也记录为成功空结果
4. 新增一个 `heliolincrr/run_daily_unknown.sh` 或 Python wrapper，负责每日选择目标夜并调用 `run_single_night.sh`
5. 默认 `SKIP_PLOTS=1`，只做提取和 summary；必要时再单独补 GIF
6. 增加产物检查：summary、unknown JSON/FITS、matched count、unknown count、ADES 行数
7. 将 unknown GIF 打包和 review CSV 模板输出接入 daily wrapper
8. 将未来 15 夜 unknown catalog 接入同一个 `assign_unknown_trksub.py`
