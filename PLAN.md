# PLAN

## Current objective

回到 `heliolincrr` 主线，先把单夜 unknown 搜索/提取自动化做稳。当前正在用新 Gaia mask 规则全量重跑 `20251120..20260514`，服务器后台 PID `1652847`。

本轮全量任务：

- 最新日志/状态表用 `ls -t /pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_full_remask_*.log` 查找
- `20251115` 的本轮中间产物已按用户要求清理，本轮不继续处理该晚
- `20251119` orbit fit_ok=`1481/1493`，已按用户要求清理并跳过
- 已撤销按 orbit fit_ok 总数跳过的规则
- `run_single_night.sh` 已新增 `MAX_UNKNOWN_LINKS_AFTER_KNOWN=200` 默认阈值；summarize/扣 known 后 unknown catalog 行数超过阈值才返回 skip code `20` 并清理 unknown/ADES 输出
- 已清理旧 unknown/heliolincrr 产物和 scoped `/tmp` 测试目录
- 保留 `/processed1/<night>/L1/L2` 和 known matched
- 顺序：link 完成 -> 轨道确认 -> summarize/扣 known 并生成 unknown catalog -> 分配 `trkSub` -> 生成 GIF/review package -> 跳过观测助手外部检查 -> 生成 ADES PSV
- 本轮不 validate、不 submit

最新覆盖状态：

- `20251115..20260514` 有 MP L2 数据的夜次共 `122` 个
- 全量重跑前 `120/122` 个完成 compute 链路
- 全量重跑会跳过缺 known matched 的夜次，例如 `20260414`
- `20251201` no-repeat-field / no tracklet 场景已改为可生成 schema-only `tracklets_ALL` 并继续写 empty unknown 结果
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
14. `mask_gaia.py` 已改为同时使用 `RA_Win/DEC_Win` 和 `RA_PSF/DEC_PSF` 匹配 Gaia，`run_single_night.sh` 使用 `--match-arcsec 1.5`
15. `package_unknown_review.py` 已新增 `*_unknown_review_full.fits`，逐 detection 记录 `RA_Win/DEC_Win`, `RA_PSF/DEC_PSF`, `trk_sub`, `linkage_id`, detection 所属 tracklet、整条 link 的 tracklet 列表和轨道摘要
16. `merge_tracklets_night.py --allow-empty` 和 `run_single_night.sh` 已支持 no-group/no-tracklet 夜生成 empty unknown 结果
17. `run_single_night.sh` 已支持 `MAX_UNKNOWN_LINKS_AFTER_KNOWN`，默认 `200`，用于跳过扣 known 后 unknown links 仍过多的异常夜

## Outstanding issues

- 服务器全量重跑任务 `1652847` 正在运行；完成后需要审计状态表和每晚产物
- `20251116..20260410` 仍有 `48` 个日历夜没有 summary，需要区分缺原始数据、失败和未跑
- `20260414` 目录名与 observing-night 归属不一致；用户判断 `0414` 实际应无这些观测，本轮按缺 known matched 跳过
- 本轮全量重跑已清空旧 unknown `trkSub` history；历史编号会随新 catalog 重新分配，且本轮不 submit
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

1. 监控最新 `unknown_full_remask_*.log` 和 status TSV，确认 `20251120` 和后续夜次通过
2. 全量完成后审计 done/fail/skip 数量、unknown count、review full FITS 行数和 ADES 行数
3. 抽查若干 review package tar，确认 GIF、review CSV、`*_unknown_review_full.fits`、`*_unknown_review_ades.fits` 均包含在包内
4. 对失败夜按 status note 分类修复或记录跳过原因
5. 新增一个正式 `heliolincrr/run_daily_unknown.sh` 或 Python wrapper，负责每日选择目标夜并调用 `run_single_night.sh`
6. 默认 `SKIP_PLOTS=1`，只做提取和 summary；必要时再单独补 GIF
7. 增加产物检查：summary、unknown JSON/FITS、matched count、unknown count、review full FITS 行数、ADES 行数
8. 将 unknown GIF 打包和 review CSV 模板输出接入 daily wrapper
9. 将未来 15 夜 unknown catalog 接入同一个 `assign_unknown_trksub.py`
