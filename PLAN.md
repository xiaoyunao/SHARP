# PLAN

## Current objective

2026-06-24 follow-up 调度入口已完成并部署：

- 新增 `survey.followup` 状态机和 `survey.apply_followup` 插入器
- follow-up 只处理 `20260624` 及之后新完成网页 check 的 true unknown 源，不回追历史 backlog
- 状态文件：
  - `/pipeline/xiaoyunao/survey/runtime/followup/followup_state.json`
- 调度规则：
  - 每个源需要两个实际满 `5` 帧的 follow-up 夜
  - 少于 `5` 帧的夜记为 failed，后续继续排
  - 从发现夜起满 `10` 天仍未完成则 abandoned
  - 目标名为 `MP_FU_<trkSub>_MP####`，既区别于普通 MP，又能进入 known/unknown 链条和历史曝光统计
  - 预测位置用当前 link detections 的线性 RA/Dec 模型；选择已有 footprint 中覆盖预测位置且中心距最近的 field
- 自动入口：
  - `survey/run_daily.sh` 生成基础计划后默认调用 `survey.apply_followup`
  - `heliolincrr/run_daily_unknown.sh` 在 unknown link 生成或已有 package skip 时调用 `associate_followup_links.py`
  - 单夜 `watch_submit_reviews.py` 由 daily unknown 启动时带 `--enable-followup`，网页 submit 处理完成后立即更新 follow-up 状态和当天 `current_plan.txt`
- 已部署到服务器：
  - `/pipeline/xiaoyunao/survey/followup.py`
  - `/pipeline/xiaoyunao/survey/apply_followup.py`
  - `/pipeline/xiaoyunao/survey/run_daily.sh`
  - `/pipeline/xiaoyunao/heliolincrr/associate_followup_links.py`
  - `/pipeline/xiaoyunao/heliolincrr/watch_submit_reviews.py`
  - `/pipeline/xiaoyunao/heliolincrr/run_daily_unknown.sh`
- 验证：
  - 本地和服务器 `py_compile`/`bash -n` 通过
  - 本地 smoke 合成 true link 后插入 `5` 条 `MP_FU_...` 行
  - 本地 smoke 合成 unknown link 后 `associate_followup_links.py` 返回 `matches=1`
  - 服务器真实计划 dry-run：`active_sources=0`, `inserted_rows=0`，未改写生产计划
- 待验证：
  - 等 `20260624` 之后第一条网页 submit true 源出现后，检查 `followup_state.json` 和 `current_plan.txt`
  - 第二天检查 L2 中 `MP_FU_...` 文件数是否触发 success/failed reconcile

2026-06-24 reviewed unknown 补报完成：

- `20251116..20260617` 历史 reviewed unknown 补报已经全部完成
- state `/pipeline/xiaoyunao/data/heliolincrr/review_submit_backlog_20251116_20260617.json`
- 最终 summary：
  - `review_packages=122`
  - `complete=122`
  - `pending=0`
  - `failed=0`
  - `submitted=35`
  - `no_observations=87`
- 顺序上从 `20251116` 连续完成到 `20260617`
- submit CSV 共 `116` 个、`4743` 行
  - `is_real=1`: `67`
  - `is_real=0`: `4676`
  - blank/invalid: `0`
- watcher 日志尾部已写出 `status=complete`，历史补报 watcher 进程已退出
- 后续只需要处理新 daily unknown check package，不再需要历史区间 backlog watcher 常驻

2026-06-23 断电恢复检查：

- 服务器 `2026-06-23 12:28 CST` 重启；上午断电属实
- `2026-06-23 09:00 CST` cron daily pipeline 在断电前正常跑完
- `2026-06-23 12:38 CST` `@reboot` daily pipeline 自动补跑完成
- `20260623` 观测脚本正常生成并发布，`n_exp=285`
- 昨晚目标夜 `20260622` 没有 processed data：
  - `/processed1/20260622` 不存在
  - known/unknown 均因 missing night dir 正常 skip
- 历史补报 watcher 已由 `@reboot` 自动恢复：
  - 当前 PID `12242`
  - 当前只有一个 watcher 进程，无重复运行
- 当前 `20251116..20260617` reviewed unknown 补报进度：
  - `review_packages=122`
  - `complete=76`
  - `submitted=19`
  - `no_observations=57`
  - `pending=46`
  - `failed=0`
  - submit CSV 共 `70` 个、`3138` 行
  - `is_real=1` 为 `35`，`is_real=0` 为 `3103`
  - 无空标签/非法值

2026-06-22 断电恢复检查：

- 服务器 `2026-06-22 10:20 CST` 重启，说明 `09:00` cron 当时错过
- `@reboot` daily pipeline 已在 `2026-06-22 10:31 CST` 自动补跑完成
- 今日观测脚本正常生成：
  - `/pipeline/xiaoyunao/survey/runtime/plans/20260622_plan.json`
  - `/pipeline/xiaoyunao/survey/runtime/plans/20260622_plan.txt`
  - `n_exp=282`
- 昨晚目标夜 `20260621` 没观测：
  - `/processed1/20260621` 不存在
  - known/unknown 均因 missing night dir 正常 skip
- recovery 回看 `20260615..20260621`：
  - `20260615/16/18/19/21` 缺 processed night dir，正常 skip
  - `20260617/20260620` 已有 known report 和 unknown review package，幂等 skip
- 历史补报 watcher 在断电后未自动恢复，已手动重启并加入 `@reboot`：
  - 当前 PID `201618`
  - crontab 新增 `@reboot sleep 900 && ... run_review_submit_backlog_watch.sh 20251116 20260617`
- `run_daily_unknown.sh` 已修复已有 review package 的恢复场景：
  - 已有正 unknown package 时不重跑 unknown
  - 会确保该夜 submit watcher 被启动
  - zero unknown package 不挂 watcher
- 当前 `20251116..20260617` reviewed unknown 补报进度：
  - `review_packages=122`
  - `complete=68`
  - `submitted=14`
  - `no_observations=54`
  - `pending=54`
  - `failed=0`
  - submit CSV 共 `62` 个，`2794` 行，`is_real=1` 为 `22`、`is_real=0` 为 `2772`，无空标签/非法值

2026-06-22 补充：

- 已生成项目状态英文单页 PPT：
  - `/Users/yunaoxiao/Desktop/smt_asteroid/outputs/smt_asteroid_project_status.pptx`
  - 桌面副本已修订：`/Users/yunaoxiao/Desktop/smt_asteroid_project_status.pptx`
  - 中文版本已按 DESI MWBP 参考 deck 风格重调：`/Users/yunaoxiao/Desktop/smt_asteroid_project_status_cn.pptx`
  - 版式参考 `/Users/yunaoxiao/Desktop/desi_mwbp_et.pptx` 的白底学术汇报风格
  - 内容概括 nightly scheduler、raw images、collaborator preprocessing/calibration/photometry、L1/L2 catalogs、known-object extraction/reporting、unknown search/review
  - 右侧放大嵌入服务器 known asteroid GIF 和本地 checked unknown candidate GIF
  - 底部状态采用高层统计：`354237` 次 known asteroid detections、`58224` 个 known objects、`22` 条 submitted reviewed unknown links
  - 桌面 known GIF 已替换为 `top03_Bakhchisaraj_20260528.gif`，并改为与 unknown GIF 一致的 `3` 帧、`3000 ms` 总时长
  - known/unknown 两个模块底部说明框已恢复上一版文字
- 若继续改 PPT，优先微调：
  - 标题和状态句措辞
  - 右侧示例 GIF 选择
  - 统计数字需按服务器最新产物刷新

2026-06-21 最新状态：

- Daily automation 已切到服务器 crontab：
  - `0 9 * * * cd /pipeline/xiaoyunao && /bin/bash /pipeline/xiaoyunao/heliolincrr/run_daily_pipeline.sh >> /pipeline/xiaoyunao/data/heliolincrr/daily_logs/cron_daily_pipeline.log 2>&1`
  - `@reboot sleep 600 && cd /pipeline/xiaoyunao && /bin/bash /pipeline/xiaoyunao/heliolincrr/run_daily_pipeline.sh >> /pipeline/xiaoyunao/data/heliolincrr/daily_logs/reboot_recovery_daily_pipeline.log 2>&1`
  - 旧的 `survey/run_daily.sh` 和 `known_asteroid/run_daily.sh` cron 已移除，避免 daily 入口不一致
- `2026-06-21 09:00` unknown 未触发的原因：
  - cron 还停留在旧的 survey/known 两条入口，没有安装 `run_daily_pipeline.sh`
  - known extraction 实际已触发并完成 `20260620`
  - known ADES/reply 初始未生成，是因为服务器 `known_asteroid/export_ades.py` 旧版本不支持 `slurm_merge_submit.sh` 传入的 `--dedup-history-root/--current-night` 参数，finalize 静默失败
- 已把正确的 `export_ades.py` 同步到服务器并通过 `py_compile`/`--help` 检查
- `heliolincrr/run_daily_pipeline.sh` 已更新：
  - survey 只跑 `RUN_DATE` 当天，用于生成今晚观测脚本
  - known/unknown 按 `RECOVERY_LOOKBACK_DAYS` 回看数据夜，适合断电后恢复
  - known daily 后等待 official known report 就绪，再运行 unknown
  - official known report ready 条件为 `<night>_matched_asteroids_ades.psv` + `<night>_mpc_reply.txt` 存在，或 status 明确 official matched `0`
- 自动 smoke 已完整跑通 `20260620`：
  - known official ADES/reply 已生成：`/processed1/20260620/L4/20260620_matched_asteroids_ades.psv`, `/processed1/20260620/L4/20260620_mpc_reply.txt`
  - unknown link `19`
  - GIF `19/19`
  - review package：`/pipeline/xiaoyunao/heliolincrr/review_packages/20260620/20260620_unknown_review_manifest.json`
  - review package manifest 中 `n_catalog_rows=19`, `n_gifs_copied=19`, `n_gifs_missing=0`
  - 网页已生成 `20260620_submit.csv`；19 条均为 `is_real=0`
  - 单夜 watcher 已自动写出 submit masked/ADES/stats，并以 `no_observations` 完成，未向 MPC 提交 unknown obsData
- 历史补报 watcher 仍在后台运行：
  - PID `2605614`
  - state `/pipeline/xiaoyunao/data/heliolincrr/review_submit_backlog_20251116_20260617.json`
  - 当前 summary：`complete=59`, `pending=63`, `failed=0`
  - `20260114` 的坏 submit CSV 已删除；网页随后重写了完整全假 submit，watcher 已处理为 `no_observations`
  - 不应与 daily 新夜 watcher 冲突，因为 state 文件和 night range 不同
- `watch_submit_reviews.py` 已避免对确定性坏输入反复 retry：
  - `blank/missing/unknown is_real`、重复 key 等错误会自动删除对应 `<night>_submit.csv`
  - 该夜回到 pending，等网页重新生成完整 submit CSV
  - submit CSV 或 review package mtime/size 变化后仍会重新处理

2026-06-20 状态：

- RA wrap 修复后的 full known/mask15 已完成并通过覆盖检查；unknown 后续应做 mask-only review rebuild，不应默认重画 GIF
- `unknown_remask_after_known_full_20260620_0346` 已停止；该流程错误地默认删除源 GIF 后重画，已完成到 `20251220`
- `remask_unknown_with_known.py` 已改为：
  - 默认保留 `/pipeline/xiaoyunao/heliolincrr/plots/<night>/unknown_link_*_<night>.gif`
  - 默认只清 `review_packages/<night>` 后重打 package
  - 当前 link 源 GIF 全存在时跳过 `plot_unknown_links.py`
  - `--skip-plots --require-existing-gifs` 可做严格 mask-only，缺图则不打不完整 package
- `20251220` smoke run 已验证：未重画 GIF，`unknown_count=20`, `n_gifs_missing=0`
- `20251221` 被旧残留 `plot_unknown_links.py` 清图后部分重画；已用 `plot_unknown_links.py --linkage-id` 只补缺失的 `8` 张 GIF，并重打 package，`n_gifs_missing=0`
- 严格 mask-only run 已完成：
  - command: `remask_unknown_with_known.py 20251222 20260617 --skip-plots --require-existing-gifs`
  - log: `/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_mask_only_strict_after_known_full_20260620_0945.log`
  - status: `/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_mask_only_strict_after_known_full_20260620_0945_status.tsv`
  - `done=91`, `skip=28`, `error=0`
  - done 夜 unknown link 合计 `3379`
  - 所有 done 夜均 `assigned_new=0`, `n_gifs_missing=0`
  - `20260103` 最终 mask 后 `unknown_count=1`
  - high-unknown skip: `20251226=366`, `20260111=378`, `20260528=926`, `20260611=653`
  - missing-input skip 共 `24` 夜，需后续分类确认
- ADES unknown 最终上报：
  - L2 catalog 已确认有 `Flux_Aper4`/`FluxErr_Aper4`
  - `logSNR` 采用 `log10(Flux_Aper4 / FluxErr_Aper4)`
  - `submit_reviewed_unknown.py` 默认透传 `--include-logsnr`
  - 如 MPC validation 临时不接受，可用 `submit_reviewed_unknown.py <night> --no-logsnr` 关闭
- Daily automation:
  - 服务器 known 匹配代码已确认包含 RA wrap 多中心查询；known Slurm 默认 `MASK_SEP_ARCSEC=1.5`
  - `known_asteroid/merge_night_parts.py` 会写 `/processed1/<night>/L4/<night>_known_asteroid_status.json`
  - `heliolincrr/summarize_single_night.py --require-mask15` 不再 fallback 到 1 角秒 matched；缺 mask15 时只允许 status 明确 empty mask
  - `heliolincrr/run_single_night.sh` 默认 `REQUIRE_KNOWN_MASK15=1`
  - `heliolincrr/run_daily_unknown.sh` 等 known mask15/status ready 后运行 unknown 并打 check 包；只有 check 包生成且 `n_catalog_rows > 0` 时才启动该夜 submit watcher
  - `heliolincrr/run_daily_pipeline.sh` 现在 survey 只跑 `RUN_DATE` 当天；known/unknown 按 `RECOVERY_LOOKBACK_DAYS` 回看实际数据夜补处理，适合 `0 9 * * *` 和 `@reboot` recovery
  - `heliolincrr/watch_submit_reviews.py` 扫描 `<night>_submit.csv`，要求已有 review manifest，用 state JSON 去重，并调用 `submit_reviewed_unknown.py` 导出/validate/submit
  - `heliolincrr/run_review_submit_backlog_watch.sh` 已用于历史补报区间 `20251116..20260617`
  - 已做 smoke：`20260103 --require-mask15` 得到 `unknown_count=1`；`20260605` known status 可表达 empty mask；watcher dry-run 可发现已有 submit CSV
  - `20260617` 真实 daily unknown 已完成：unknown link `0`，check package 已生成
  - 历史补报 watcher 正在后台运行：PID `2530505`，state `/pipeline/xiaoyunao/data/heliolincrr/review_submit_backlog_20251116_20260617.json`
  - 历史补报当前 summary：`review_packages=122`, `complete=42`, `pending=80`, `failed=0`
  - `20251116..20251231` 子区间 `37` 个 review package 已全部完成，真源 `20` 条 link / `60` 次探测
  - 该区间真源 GIF 和测光星表已导出到本地 `/Users/yunaoxiao/Desktop/true_unknown_20251116_20251231`
  - 已正式提交 reviewed unknown：`20251121` -> `2026-06-20T05:00:14.546_0000C0eY`; `20251124` -> `2026-06-20T05:02:22.435_0000C0ea`; `20251125` -> `2026-06-20T05:00:53.443_0000C0eZ`
  - cron 已在 `2026-06-21` 改为总入口 `run_daily_pipeline.sh`
- Known anomaly:
  - `20260111` known candidate 数正常（`all_asteroids=106042`），但 matched 极少（official `180`, mask15 `419`）
  - 相邻夜 matched 为几千到一万量级，且 `20260111` 候选星等分布并不更暗
  - 抽样显示 L2 catalog/WCS 相对 known ephemeris 有约 `8..11"` RA 方向系统偏移；Win/PSF 坐标一致
  - 初步判断不是 known 查询漏候选，而是该夜 L2 astrometry/WCS 有系统偏移；修复前不建议把该夜 known/unknown 作为可靠上报结果

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
- known 提取根因继续修复：L2 catalog BINTABLE 不再用 `NAXIS1/NAXIS2` 当图像尺寸；跨 RA=0 视场不能只改用 `RA=0` 查询中心，必须同时查原始 WCS 中心、`RA=0` 和 `RA=359` 后合并去重；official known matched 保持 `1.0"` 用于 ADES/MPC，上游另写 `*_matched_asteroids_mask15.fits` 用 `1.5"` 给 unknown 扣除
- 修复后临时验证：`3331 Kvistaberg`, `1522 Kokkola`, `2220 Hicks`, `168 Sibylla`, `1795 Woltjer`, `8219 (1996 JL)` 均可进入 known matched；临时验证目录已删除
- `2026-06-19` 再次用当前 `20260103` unknown 实时 RA/Dec 做 JPL/Horizons 复查，发现 remask 后新增的 `00001iF/00001iG/00001iH/00001iI` 仍分别是 `2125 Karl-Ontjes`, `2728 Yatskiv`, `3856 Lutskij`, `1989 Tatry`
- 已修复跨 RA=0 多中心查询并烟测 `OBJ_MP_1060_0055_cat.fits.gz`：上述四个对象均进入 all/matched/mask15，匹配到当前 unknown 的对应 `objID`
- `known_rematch_20260618_111935` 仍基于旧 RA=0 单中心逻辑生成，不能作为最终人工 check 输入；已废弃
- 已停止旧 GIF repair，并清空 review package/L4 中的 submit CSV 和 submit 派生产物
- 已启动新全量重跑：`RUN_ID=known_wrapfix_20260619_122524`
  - known/mask15 forced rematch：`20251116..20260617`
  - unknown remask：`20251116..20260617`
  - 参数：`FORCE_EXTRACT=1`, `MAX_PARALLEL=16`, `MASK_SEP_ARCSEC=1.5`, `MASK_MATCHED_SUFFIX=_mask15`
  - logs: `/pipeline/xiaoyunao/known_asteroid/runtime/logs/known_wrapfix_20260619_122524.log`, `/pipeline/xiaoyunao/known_asteroid/runtime/logs/known_wrapfix_20260619_122524_driver.log`
  - remask status: `/pipeline/xiaoyunao/data/heliolincrr/batch_logs/known_wrapfix_20260619_122524_unknown_remask_status.tsv`
- 本轮 remask 已改为每夜打包前默认清理旧错误 artifact：删除 `plots/<night>/unknown_link_*_<night>.gif` 和整个 `review_packages/<night>` 后再重画 GIF、重建 review package
- 批量 driver 已停用：旧 finalize 在 `FORCE_EXTRACT=1` 时会因为夜级 FITS 已存在而跳过 merge，导致 night-level known/mask15 仍是旧结果；已修复 `slurm_merge_submit.sh`，后续用逐夜手动流程推进
- 手动流程：每夜/每批提交 `ka_match` array，确认 array 离队后运行 `FORCE_MERGE=1 FORCE_EXTRACT=1 ./slurm_merge_submit.sh <night> false`，确认夜级 all/matched/mask15 时间戳为新 run，再进入下一夜
- 当前手动 merge 已完成到 `20251206`；`20251207..20260104` 是已提交但仍 active/pending 的 array，`20260105` 之后尚未继续提交
- 修复前生成的 known matched、unknown links 和 review packages 需要重跑后才能用于人工 check 后上报
- 已启动批量重跑：`RUN_ID=known_rematch_20260618_111935`
  - known forced rematch：中途因 Slurm `sbatch` 临时 `Resource temporarily unavailable` 停在 `20251228`；已给 submit/driver 脚本加 retry 并从 `20251228` 续跑
  - 最新进度口径改为三分类：有 `OBJ_MP` 输入的夜才作为应跑分母；无 `OBJ_MP` 或无 L2 的夜作为不适用/跳过
  - 当前应跑夜次 `128`，跳过/不适用 `22`；普通 known 和 mask15 本轮均已完成 `29` 夜，完成范围到 `20251217` 且包含 `20251221`
  - retry 后 `20251228` 已成功提交 array `190635` / finalize `190636`
  - `20260106` 曾因 `sbatch` stderr warning 污染 array job id，导致 finalize dependency 非法并停住；已修复 stdout/stderr 分离与 job id 解析，取消坏 array `191015`，并从 `20260106` 续跑成功，新 job 为 array `195437` / finalize `195438`
  - 后续 driver 因缺 L2 夜 `/processed1/20260203/L2` 退出；已改为缺 L2 时 skip 并继续，当前从 `20260205` 续跑
  - duplicate cleanup：并发残留导致的重复/坏 job `185697..185706` 已取消
  - `2026-06-19` known rematch 后的 unknown remask 已完成；旧人工 `<night>_submit.csv` 全部删除，review package 下当前 submit CSV 数量为 `0`，所有人工 check 需要基于新 review package 重做
  - unknown remask：已处理 `20251116..20260616`，不处理 `20260617` unknown
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

1. 明早检查 `cron_daily_pipeline.log`、`<run_date>_daily_pipeline.log` 和 `<run_date>_unknown_daily.log`，确认 09:00 cron 自然触发且不是人工 smoke
2. 继续让历史补报 watcher 等待 `20251116..20260617` 剩余 submit CSV，阶段性检查 `complete/pending/failed`
3. 遇到网页生成的不完整 submit CSV，确认 watcher 删除该文件并让该夜回到 pending
4. 抽查修复后 `20260102/20260103` 中 JPL 命中的已知小行星是否从 unknown 中消失
5. 以修复后产物继续网页筛选；旧 submit CSV 只能作为人工判断参考，不能直接上报
6. 对已完成 submit CSV 的夜次用 watcher/`submit_reviewed_unknown.py` 生成筛选后 PSV 并正式 submit；异常时先看 validate/submit reply
7. 对 `20260528`、`20260611` 单独复盘 high unknown 原因，重点查 known subtraction、mask 后源密度、dense group 和观测场区
8. 复盘 `20260111` L2 astrometry/WCS 约 `8..11"` RA 系统偏移；修复前不要把该夜 known/unknown 当可靠上报
9. 如 PPT 还要继续打磨，再补 detection histogram 的 cumulative 版或 orbit confirm 的 good/rejected 双栏对照图
