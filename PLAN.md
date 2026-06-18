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
- `20260103` 的 6 条旧 unknown link 虽均被人工判为真源，但 JPL 临时检查确认其中至少 `5` 条为已知小行星；该夜应作为 known 提取漏检回归样本，不再作为 unknown submit 正例
- `2026-06-18` 迁移后 `trkSub` history 为 `4917` 条、全部 7 位，范围 `0000001..00001hj`；postcheck dry-run 无 8 位残留
- 早期 MPC unknown 正式 submit 为 `20260220`，submission ID `2026-05-14T03:00:19.564_0000Bx4c`；当时未经过人工 check，提交版本为 `34` 条 link、`102` 行 obsData。当前服务器 remask 后 `20260220` 版本为 `29` 条 link、`87` 行 obsData
- 8 位 trkSub 迁移备份 `/pipeline/xiaoyunao/data/heliolincrr/backups/trksub_8char_20260618_094807` 已按用户确认删除
- `2026-06-18` JPL 临时检查发现人工真源中至少 `6` 条实际为已知小行星；known 提取存在系统性漏检
- known 提取根因已修复：L2 catalog BINTABLE 不再用 `NAXIS1/NAXIS2` 当图像尺寸，跨 RA=0 视场改用 `RA=0` 查询中心；official known matched 保持 `1.0"` 用于 ADES/MPC，上游另写 `*_matched_asteroids_mask15.fits` 用 `1.5"` 给 unknown 扣除
- 修复后临时验证：`3331 Kvistaberg`, `1522 Kokkola`, `2220 Hicks`, `168 Sibylla`, `1795 Woltjer`, `8219 (1996 JL)` 均可进入 known matched；临时验证目录已删除
- 修复前生成的 known matched、unknown links 和 review packages 需要重跑后才能用于人工 check 后上报
- 已启动批量重跑：`RUN_ID=known_rematch_20260618_111935`
  - known forced rematch：中途因 Slurm `sbatch` 临时 `Resource temporarily unavailable` 停在 `20251228`；已给 submit/driver 脚本加 retry 并从 `20251228` 续跑
  - 最新进度口径改为三分类：有 `OBJ_MP` 输入的夜才作为应跑分母；无 `OBJ_MP` 或无 L2 的夜作为不适用/跳过
  - 当前应跑夜次 `128`，跳过/不适用 `22`；普通 known 和 mask15 本轮均已完成 `29` 夜，完成范围到 `20251217` 且包含 `20251221`
  - retry 后 `20251228` 已成功提交 array `190635` / finalize `190636`
  - `20260106` 曾因 `sbatch` stderr warning 污染 array job id，导致 finalize dependency 非法并停住；已修复 stdout/stderr 分离与 job id 解析，取消坏 array `191015`，并从 `20260106` 续跑成功，新 job 为 array `195437` / finalize `195438`
  - 后续 driver 因缺 L2 夜 `/processed1/20260203/L2` 退出；已改为缺 L2 时 skip 并继续，当前从 `20260205` 续跑
  - duplicate cleanup：并发残留导致的重复/坏 job `185697..185706` 已取消
  - unknown remask 尚未开始
  - unknown remask：driver 完成 known 提交并轮询全部 finalize job 离队后提交 `20251116..20260616`，不处理 `20260617` unknown
  - logs: `/pipeline/xiaoyunao/known_asteroid/runtime/logs/known_rematch_20260618_111935.log`, `/pipeline/xiaoyunao/known_asteroid/runtime/logs/known_rematch_20260618_111935_driver.log`
  - remask status: `/pipeline/xiaoyunao/data/heliolincrr/batch_logs/known_rematch_20260618_111935_unknown_remask_status.tsv`

短期目标：

- 暂停基于旧 unknown 产物的上报，等待 `known_rematch_20260618_111935` 的 known/remask 批处理完成
- 批处理完成后先抽查用户已人工检查的 `20260101..20260107` 和追加检查的 `20260110`，确认已知小行星不再进入 unknown
- 推进 unknown 人工 check 时，以修复后重建的 review package 为准
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
- 人工复核目前只填写少量条目；修复 known 提取并重建对应 unknown/review 产物前，已筛选的 `is_real=1` 不应导出为正式 unknown ADES
- 修复前的 known 提取可能在 L2 BINTABLE 尺寸、RA=0 跨界、`1.0"` 半径三方面漏检已知小行星，需要确定重跑范围
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

1. 监控 `known_rematch_20260618_111935`：确认 known submit driver 完成、dependent remask job 已提交并完成
2. 抽查修复后 `20260102/20260103` 中 JPL 命中的已知小行星是否从 unknown 中消失
3. 汇总 remask status TSV：done/skip/error、unknown_count、review_full_rows、n_gifs_missing
4. 以修复后产物重新交给网页筛选；旧 submit CSV 只能作为人工判断参考，不能直接上报
5. 建立人工 check 完成检查脚本或命令：查找 `<night>_submit.csv`，校验每个 unknown `tracklet_id` 都有明确 `0/1`
6. 对已完成 submit CSV 的夜次用 `submit_reviewed_unknown.py <night> --validate` 生成筛选后 PSV 并做 MPC test validation，不能再用旧 `20260103` 作为正例
7. 抽查 `20260529`, `20260530`, `20260601` review package 内容并交给网页筛选
8. 对 `20260528`、`20260611` 单独复盘 high unknown 原因，重点查 known subtraction、mask 后源密度、dense group 和观测场区
9. 决定 `20260605` 这种 known-only 无 matched detections 的夜次是否要写 schema-only matched 文件再跑 unknown
10. 新增正式 `heliolincrr/run_daily_unknown.sh` 或 Python wrapper，负责每日选择目标夜并调用 `run_single_night.sh`
11. 默认 `SKIP_PLOTS=1`，只做提取和 summary；必要时再单独补 GIF
12. 增加产物检查：summary、unknown JSON/FITS、matched count、unknown count、review full FITS 行数、ADES 行数
13. 将 unknown GIF 打包和 review CSV 模板输出接入 daily wrapper
14. 待 unknown 人工 check 与上报链路稳定后，再新增 `daily_monitor`/recovery 脚本
15. 将未来 15 夜 unknown catalog 接入同一个 `assign_unknown_trksub.py`
16. 如 PPT 还要继续打磨，再补 detection histogram 的 cumulative 版或 orbit confirm 的 good/rejected 双栏对照图
