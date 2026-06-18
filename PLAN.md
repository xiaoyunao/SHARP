# PLAN

## Current objective

回到 `heliolincrr` 主线，把单夜 unknown 搜索/提取自动化做稳，并在
`unknown_full_remask_20260514_171651` 全量重跑完成后做产物审计和 daily wrapper 收口。

最新全量重跑状态：

- 日志：`/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_full_remask_20260514_171651.log`
- 状态表：`/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_full_remask_20260514_171651_status.tsv`
- 范围：`20251119..20260514`
- 已完成时间：`2026-05-15T16:41:11+08:00`
- 状态表共 `177` 行：`done=115`, `skip=62`
- skip 原因：`no_l2=56`, `unknown_links_after_known_gt_200=3`, `no_mp_l2=2`, `missing_known_matched=1`
- 完成夜次累计 unknown link 数：`4764`
- 最后一行：`20260514 skip rc=0 no_l2`
- `2026-06-03` 复查：当前无 unknown 相关后台进程；combined link/detection 表仍与状态表一致

人工复核和上报状态：

- review package manifest/tar 共 `117` 个，包含早期 `20251116` 和 `20251118`
- review CSV 合计 `4848` 行，其中 `52` 行已填写 `is_real`，分布为 `1=4`, `0=48`
- 已填写 review 的夜次：`20260309`, `20260504`, `20260507`, `20260508`, `20260509`, `20260512`
- `export_unknown_ades.py` 已支持 `--review-csv`, `--require-review`, `--validate`, `--submit`
- `export_unknown_ades.py` 已支持网页最终提交文件 `--submit-csv`，默认读取 `/pipeline/xiaoyunao/heliolincrr/review_packages/<night>/<night>_submit.csv`
- `submit_reviewed_unknown.py` 已作为人工 check 完成后的单夜入口，默认生成 masked JSON/FITS、ADES PSV 和 stats，不自动 validate/submit
- `trkSub` 已从旧的 8 位 leading-zero managed base62 迁移为 MPC 要求的 7 位；服务器 history、unknown JSON/FITS/PSV、review/submit CSV、manifest/FITS、GIF 文件名、review tar.gz 和 combined 表均已迁移
- `run_single_night.sh` 默认关闭 unknown ADES export/validate/submit
- crontab 目前只有 known asteroid daily，没有 unknown daily wrapper
- 未发现 `*_unknown_mpc_reply.txt` 或 `*_unknown_validate_reply.txt`，人工筛选后的 unknown 上报尚无完成证据

`20260514` 之后新夜次处理：

- `unknown_after_20260514_20260603_162313` 已确认中断：状态表只有表头，日志停在 `20260528 mask_gaia`
- `unknown_after_20260514_20260617_161428` 已完成
- 日志：`/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_after_20260514_20260617_161428.log`
- 状态表：`/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_after_20260514_20260617_161428_status.tsv`
- 成功夜次：
  - `20260529`: `unknown=29`, `review_full_rows=87`, `ades_rows=87`, `n_gifs_missing=0`
  - `20260530`: `unknown=21`, `review_full_rows=63`, `ades_rows=63`, `n_gifs_missing=0`
  - `20260601`: `unknown=19`, `review_full_rows=57`, `ades_rows=57`, `n_gifs_missing=0`
- 跳过夜次：
  - `20260528`: `unknown_count=1168 > 200`
  - `20260611`: `unknown_count=653 > 200`
  - `20260605`: known-only 补跑无 matched detections，仍缺 matched FITS
- 本轮只生成结果/review package/ADES PSV；未 validate，未 submit

`2026-06-17` 复查补充：

- 当前有 unknown JSON/review manifest 的夜次共 `120` 个
- 其中 positive unknown 夜次 `114` 个，zero unknown 夜次 `6` 个，unknown link 合计 `4917`
- review CSV 共 `120` 个，`4917` 行；已填写 `110` 行，其中 `is_real=1` 为 `15` 条，`is_real=0` 为 `95` 条
- 2025 年子集共有 `37` 个夜次，其中 `35` 个 positive unknown、`2` 个 zero unknown，unknown link/review CSV 行合计 `1948`
- 2025 年 review CSV 尚未开始填写：`filled=0`, `is_real=1` 为 `0`, `is_real=0` 为 `0`
- 2025 年 review manifest、CSV、`*_unknown_review.tar.gz` 和 GIF 均完整，GIF 缺失合计 `0`
- `20260528` 高 unknown 主要来自异常密集 tracklet/link：`tracklets_total=92792`, `links_total=1170`, `unknown_count=1168`，集中在 group `47` 和 `OBJ_MP_1428_0236/0287/0309`
- `20260611` 高 unknown 同样来自密集场/文件组：`tracklets_total=267450`, `links_total=1319`, `unknown_count=653`；同时识别出 `all_same_asteroid=636` 个 known links
- `20260605` 未自动触发是因为缺少应在 `2026-06-06 09:00` 运行的 daily log；手动补跑后两个 MP 文件均无 matched detections，仍不能进入当前 unknown 流程
- 服务器 `2026-06-17 16:00:12` 重启，普通 cron 不会补跑当天 `09:00` 已错过的任务
- 已手动补跑 `2026-06-17` survey daily，生成并发布 `20260617_plan`
- 已手动补跑 `2026-06-17` known daily；target night `20260616` 因 `/processed1/20260616` 缺失而 skip
- 用户已确认 recovery/monitor 暂缓，待 unknown 人工 check 与上报链路稳定后再实现；当前重点是人工 check
- `20260512` 已存在网页生成的 `20260512_submit.csv`；dry run 成功生成 `1` 个 masked unknown link 和 `3` 行 ADES obsData，未 validate，未 submit
- `20260103` 的 6 条 unknown link 均被人工判为真源；该夜 `tracklets_total=718`, `links_total=79`, `fit_ok all_non_asteroid=6`，5 条非已知候选已在 orbit 阶段因 `max_v` 被剔除，可作为 submit CSV/validate 的正例夜次
- `2026-06-18` 迁移后 `trkSub` history 为 `4917` 条、全部 7 位，范围 `0000001..00001hj`；postcheck dry-run 无 8 位残留

短期目标：

- 推进 unknown 人工 check，先按年份/夜次分批生成 `<night>_submit.csv`
- 明确人工 check 完成判定：每个正 unknown 夜次存在 `<night>_submit.csv`，且每个 unknown `tracklet_id` 都有明确 `0/1`；zero unknown 夜次可直接视为完成
- 对已完成 check 的夜次用 `submit_reviewed_unknown.py <night>` 重新导出 filtered/masked 产物和 ADES，并先 validate，不自动 submit
- 抽查 review package tar 是否包含 GIF、review CSV、`*_unknown_review_full.fits`、`*_unknown_review_ades.fits`
- 对 skip 夜按原因分类记录，尤其是 `unknown_links_after_known_gt_200` 和 `missing_known_matched`
- 新增正式 daily unknown 入口，负责选择目标夜、防重复运行、日志记录和产物检查
- daily recovery/monitor 入口暂缓，等 unknown 人工 check 与上报链路稳定后再做
- 上报链路默认关闭；导出/提交必须显式打开
- GIF 可视化默认可跳过或限量，避免拖慢主计算

PPT 素材旁支已完成一批图件，统一在服务器
`/pipeline/xiaoyunao/ppt_assets_20260414` 管理，后续只按汇报版式做小幅微调。

## Milestones

1. `20260220` 已作为单夜 unknown 主测试夜完整跑通
2. `run_single_night.sh` 已串起 Gaia mask、tracklet、linear link、orbit confirm、summary、unknown GIF
3. `summarize_single_night.py` 已生成 unknown fit_ok catalog，并按 known asteroid match 结果扣除已知小行星
4. 已新增 `--trk-sub-map`、人工复核 CSV 过滤接口和 `export_unknown_ades.py` unknown ADES 导出器
5. 已实现全局 `trkSub` 历史分配：`/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl`
6. `20260220` unknown ADES 已通过 MPC test validation，并完成一次正式 submit
7. 0-link 夜次现在能生成 empty RR/orbit/summary/unknown JSON/FITS，并让 assign/plot 正常 0 行退出
8. 真实 0-link 夜 `20260128` 和 `20260224` 已重跑成功，均生成 unknown=0 的 summary 和 empty unknown FITS/JSON
9. `make_tracklet_linreproj.py` 已新增 `--max-tracklets-per-group`，默认 `100000`；`20260503` 重跑确认 dense group 被跳过
10. `20260503` 已重新生成 unknown GIF，并打包 review package，`n_gifs_missing=0`
11. `20251128`, `20260412`, `20260430` 已补齐空 unknown FITS/JSON
12. `mask_gaia.py` 已改为同时使用 `RA_Win/DEC_Win` 和 `RA_PSF/DEC_PSF` 匹配 Gaia，`run_single_night.sh` 使用 `--match-arcsec 1.5`
13. `package_unknown_review.py` 已新增 `*_unknown_review_full.fits` 和 `*_unknown_review_ades.fits`
14. `merge_tracklets_night.py --allow-empty` 和 `run_single_night.sh` 已支持 no-group/no-tracklet 夜生成 empty unknown 结果
15. `run_single_night.sh` 已支持 `MAX_UNKNOWN_LINKS_AFTER_KNOWN`，默认 `200`
16. `unknown_full_remask_20260514_171651` 已跑完整个区间，状态表统计为 `done=115`, `skip=62`
17. PPT 素材已新增已知小行星 detection histogram 和 orbit confirm diagnostic card 脚本
18. 已按人工复核包 `*_unknown_review_full.fits` 合并出全量 unknown link 表，桌面副本在 `/Users/island/Desktop/unknown_review_combined_20260516/`

## Outstanding issues

- 全量重跑完成后只做了状态表、combined 表和 manifest required outputs 的初步审计；还未逐夜深审 GIF/ADES/review full 内容
- `20260414` 目录名与 observing-night 归属不一致；本轮因 `missing_known_matched` 跳过
- `unknown_links_after_known_gt_200` 的 3 个 skip 夜需要记录具体 night 并决定是否人工排除或调参重跑
- 当前单夜自动化还没有“每日选择目标夜 + 防重复 + 日志 + 产物检查”的外层入口
- 当前没有“错过 09:00 cron 后自动补跑/报警”的 recovery 机制；服务器重启后需人工确认当天 survey/known/unknown 任务是否已执行
- GIF 可视化很慢，应在自动提取中默认可跳过或限量
- 15 夜流程尚未接入同一个 `trkSub` history
- unknown 真实提交策略仍需收紧；后续 daily wrapper 默认不应自动 submit，人工筛选后需先 validate 再显式 submit
- 人工复核目前只填写少量条目；已筛选的 `is_real=1` 还未重新导出为 filtered unknown ADES
- `20260528` 和 `20260611` unknown 数过高，需要单独复盘
- `20260605` known-only 补跑无 matched detections，是否允许 no matched 夜继续 unknown 需要决定
- PPT 侧 `known_object_detection_histogram` 已可出图，但 summary JSON 是否稳定记录该字段仍需补查

## Validation criteria

- 给定一个存在 `/processed1/<night>/L2` 和 known asteroid match 的夜次，自动入口能完整运行或明确跳过
- 成功夜必须存在，即使 unknown=0 也要写出 schema-only FITS/JSON：
  - `/pipeline/xiaoyunao/data/heliolincrr/<night>/analysis/<night>_single_night_summary.json`
  - `/pipeline/xiaoyunao/data/heliolincrr/<night>/analysis/<night>_single_night_summary.txt`
  - `/processed1/<night>/L4/<night>_unknown_links.json`
- summary 中 `counts.matched_detections_total > 0`，否则不能认为已完成已知小行星扣除
- 若启用 `trkSub`，每个 unknown link 的 `trk_sub` 必须符合 MPC 当前要求：不超过 7 个字符；自动分配器使用固定 7 位 `[0-9a-zA-Z]`
- `trkSub` 分配必须写入全局 history，并且重复运行同一 catalog 不新增重复记录
- 若启用人工复核或网页 submit CSV，`is_real=0` 的 `trk_sub` 不得进入 ADES PSV
- 网页 submit CSV 模式必须在导出前检查所有 unknown `tracklet_id` 都有明确 `0/1`
- 自动入口日志需记录 target night、skip/run reason、exit code、核心产物路径和 unknown count
- PPT 已知小行星图可稳定写出 `known_object_detection_histogram.png`
- PPT orbit confirm 图可稳定写出 `20260220_link234_orbit_fit_diagnostic*.png`
- 全量人工复核包合并表应满足 link-level 行数 `4764`、detection-level 行数 `14325`、逐夜 link 数与状态表无差异

## Next recommended steps

1. 把 2025 年 `35` 个 positive unknown 夜次交给网页筛选；`20251128` 和 `20251201` 为 zero unknown，可直接标记完成
2. 建立人工 check 完成检查脚本或命令：查找 `<night>_submit.csv`，校验每个 unknown `tracklet_id` 都有明确 `0/1`
3. 对已完成 submit CSV 的夜次用 `submit_reviewed_unknown.py <night> --validate` 生成筛选后 PSV 并做 MPC test validation，可先用正例夜 `20260103`
4. 抽查 `20260529`, `20260530`, `20260601` review package 内容并交给网页筛选
5. 对 `20260528`、`20260611` 单独复盘 high unknown 原因，重点查 known subtraction、mask 后源密度、dense group 和观测场区
6. 决定 `20260605` 这种 known-only 无 matched detections 的夜次是否要写 schema-only matched 文件再跑 unknown
7. 新增正式 `heliolincrr/run_daily_unknown.sh` 或 Python wrapper，负责每日选择目标夜并调用 `run_single_night.sh`
8. 默认 `SKIP_PLOTS=1`，只做提取和 summary；必要时再单独补 GIF
9. 增加产物检查：summary、unknown JSON/FITS、matched count、unknown count、review full FITS 行数、ADES 行数
10. 将 unknown GIF 打包和 review CSV 模板输出接入 daily wrapper
11. 待 unknown 人工 check 与上报链路稳定后，再新增 `daily_monitor`/recovery 脚本
12. 将未来 15 夜 unknown catalog 接入同一个 `assign_unknown_trksub.py`
13. 如 PPT 还要继续打磨，再补 detection histogram 的 cumulative 版或 orbit confirm 的 good/rejected 双栏对照图
