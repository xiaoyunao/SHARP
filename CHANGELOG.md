# CHANGELOG

本文件记录仓库中值得保留的版本级行为变化，重点覆盖服务器基线流程、
调试中确认过的默认值变更，以及会影响结果解释的输出变化。

## 2026-09-04

### 修复 review 后即时 follow-up 的运行环境

- `watch_submit_reviews.py` 新增 `--followup-python`，默认使用 survey 运行环境 `/home/smtpipeline/Softwares/miniconda3/bin/python`
- unknown 主流程继续使用 `heliolinc` 环境；只有 `survey.apply_followup` 切换到包含 `astroplan` 的 survey Python
- 20260904 最终计划图已按注入后的 453 行重画；follow-up footprint 使用青绿色独立显示，标题记录 55 条 follow-up exposure

### 防止并发 watcher 重复提交 unknown

- `watch_submit_reviews.py` 现在使用与 state 相邻的专用文件锁；每轮扫描拿锁后重新加载最新状态
- “检查是否完成—执行 MPC 提交—原子写回状态”在同一个跨进程事务内完成，避免多个单夜 watcher 用旧快照相互覆盖
- 生产 state 已从 watcher 日志恢复 20260830、20260901、20260902 的首次有效提交记录，并保留两次重复尝试 ID 供审计

## 2026-09-03

### 统一 known 与 survey 全天图的 RA 方向

- 修正 `plot_known_asteroids.py` 中错误的 `180-RA` Aitoff 变换，改为与 survey 图一致的标准 `[-180,180)` RA wrap
- known 的 survey coverage、全天已知小行星分布和 field yield 三类图统一显示 `210,240,270,300,330,0,30,...,150` 横轴标签
- 该改动只影响绘图坐标和标签，不改变 scheduler 计划、天区编号、匹配结果或上报数据

### 放宽 Gaia 残留坏帧门控

- 将 unknown 链生成前的整帧拒绝阈值从 `>2000` 调整为 `>3000`，减少正常密星场被保守丢弃
- `mask_gaia.py`、单夜入口和 daily 入口的默认值保持一致；仍可用 `MAX_GAIA_RESIDUAL_ROWS` 显式覆盖

## 2026-09-02

### unknown 前 Gaia 残留坏帧门控

- `mask_gaia.py` 新增逐帧 `max_residual_rows` 门控，默认值为 `2000`
  - Gaia 去除后残留数严格大于阈值时，将整帧标记为 `rejected_high_residual`
  - 被拒绝帧不写入 `mask_gaia`，并删除可能存在的同名旧输出
  - `mask_gaia.log` 记录拒绝数量、阈值和逐帧是否写出
- `run_single_night.sh` 与 `run_daily_unknown.sh` 新增并透传 `MAX_GAIA_RESIDUAL_ROWS`
- `FORCE_MASK_GAIA=1` 现在会清理目标夜旧的 `mask_gaia` 和 `tracklets_linreproj` 后完整重建，避免历史异常输出被继续复用
- 生产回放结果：
  - `20260830`：拒绝 `11/331` 帧，unknown `847 -> 41`，复核包含 41 个 GIF/125 条观测
  - `20260901`：拒绝 `20/392` 帧，unknown `1339 -> 56`，复核包含 56 个 GIF/169 条观测
- 两夜复核包均无缺失 GIF，JSON/FITS/manifest/tar 完整性验证通过；本次未自动提交 MPC

## 2026-08-09

### SHARP PASP 冻结分析交付

- 新增 `paper_analysis_20260803/`，集中保存 2026-08-03 论文快照的：
  - raw/L2/known/unknown/plan/log/review manifests 与 provenance hashes
  - P0/P1 统计程序、5 张 paper table、12 张 PDF/PNG figure
  - executed notebook、代码-稿件冲突审计、作者输入清单和 QA 报告
- 完成 128 夜 observer-location sensitivity run：
  - 历史 orbit-confirmation 位置与 960 m scheduler/known 位置之间无 `fit_ok`
    或 `is_good` 分类翻转
  - 保留短弧轨道元素不稳定警告，不把单夜元素用于种群推断
- 冻结 unknown 人工口径为 `68 -> 67 -> 58`，并区分 14,299 个
  linkage-detection memberships 与 14,159 个唯一 detection keys
- 全量验证结果为 63 PASS / 0 FAIL / 1 author-signoff SKIP；12/12 图件通过
  PNG/PDF 渲染检查；555 个冻结制品通过 SHA256 manifest 复验
- 本次交付不修改论文工程或生产 pipeline 行为

## 2026-05-14

### heliolincrr 空结果与异常 group 保护

- 单夜 linear link、orbit confirm、summary 现在能写出 0 行但 schema 完整的 FITS
  - `n_links=0` 不再触发 Astropy 空表无列错误
  - unknown catalog 即使为空也保留 JSON `[]` 和 schema-only FITS
  - unknown plot 对空 catalog 写空 summary，并以 `n_gifs=0` 正常退出
- `make_tracklet_linreproj.py` 新增 `--max-tracklets-per-group`
  - 默认值 `100000`
  - 单个 exposure group 的候选 tracklet 超过阈值时直接跳过该 group，不写 group FITS
  - `run_single_night.sh` 可用 `MAX_TRACKLETS_PER_GROUP` 覆盖默认值
- 服务器验证：
  - 临时 0-link 夜 `20990101` 成功跑完 linear link、orbit confirm、summary、assign、plot 的空结果链路
  - 真实 0-link 夜 `20260128` 和 `20260224` 已成功生成 unknown=0 的 summary 和 empty unknown FITS/JSON
  - 清理并重跑 `20260503` 后，group `073` 和 `075` 被阈值跳过
  - `20260503` 从旧结果 `688008` 条 tracklet、`1591` 条 unknown 降为 `53193` 条 tracklet、`95` 条 unknown

## 2026-04-09

### known_asteroid 上报去重保护

- `submit_pipeline_slurm.sh` 不再无条件吞掉整个 `L2/` 目录
- 新增 `known_asteroid/build_file_manifest.py`
  - 读取每个 `L2/*MP*` 文件头里的 `OBS_DATE` / `DATE-OBS`
  - 按 `Asia/Shanghai` 本地时间、`12:00` 换夜规则重判 observing night
  - 只有真正属于目标夜次的文件才进入 manifest
- `export_ades.py` 新增两层提交前去重：
  - 同一批导出内按 `(object, obsTime)` 保留一条最优观测
  - 对照 `/processed1/*/L4/*_matched_asteroids_ades.psv` 过滤已提交过的观测
- `slurm_merge_submit.sh` 现在默认把 `ROOT_DIR` 和当前 `NIGHT` 传给导出器，启用历史去重
- 背景：
  - 服务器历史提交中存在跨相邻夜次的重复上报
  - 已确认 `20251225` 与 `20251226` 两夜含大量相同 `(object, obsTime)` 观测
  - 源头原因是 `20251226/L2` 混入了 `166` 个实际属于 `20251225` observing night 的文件
  - 最近 `20260331` 到 `20260407` 的提交未发现同类重复
  - 该保护用于在上游夜目录再次重叠时，阻止重复观测继续送到 MPC

### heliolincrr RR 默认参数更新

- 将单夜 RR 的默认 `tol` 从 `0.02` 调整为 `0.03`
- 将 RR 主程序 CLI 默认 `k-neighbors-cap` 调整为 `300`
- 将单夜脚本 `run_single_night.sh` 中的 `RR_NIGHT_KCAP` 调整为 `300`
- 将单夜 RR 的默认 `max-v-kms` 从 `200` 调整为 `30`
- 将 RR 主程序 CLI 默认 `min-init-earth-au` 从 `0.02` 调整为 `0.01`
- 选择依据：
  - 当前目标从“优先压歧义”切换为“优先保召回，把噪声留给下一步轨道拟合”
  - 在 `ref-dt-days=0.05, tol=0.03` 下，`kcap=300` 相比 `kcap=200` 将
    `rr_given_tracklet` 从 `96.31%` 提升到 `97.53%`
  - 代价是 linkage 和成员规模继续增加，且 hits-only 歧义度 `p90` 从 `8`
    升到 `9`
  - 在继续固定 `ref-dt-days=0.05, tol=0.03, k-neighbors-cap=300` 后，
    `max-v-kms=30` 相比旧基线 `200` 将
    `n_links` 从 `12274` 降到 `8697`，`n_member_rows` 从 `39806` 降到 `27204`
  - 同时 `rr_given_tracklet` 从 `2247/2304=97.53%` 提升到
    `2259/2304=98.05%`，hits-only `p90` 从 `9` 降到 `8`
  - 在固定 `max-v-kms=30` 后，`min-init-earth-au=0.02` 与 `0.01`
    在 `20260220` 上结果完全一致，因此将 CLI 默认值也统一到
    单夜脚本当前使用的 `0.01`

## 2026-04-08

### heliolincrr 默认参数更新

- 在边缘壳层过滤基线下，围绕 `20260220` 完成 `vmax` / `vmin` 单夜扫描
- 将 `make_tracklet` 默认值更新为：
  - `vmax = 63.0`
  - `vmin = 3.0`
- 选择依据：
  - `vmax=63` 相比旧默认 `80`，显著提高 purity，同时只小幅降低 completeness
  - 在固定 `vmax=63` 时，`vmin` 从 `0.0` 提高到 `3.0` 的过程中，`completeness_object_fraction` 保持 `0.87391` 不变
  - 当 `vmin=5.0` 时，`completeness_object_fraction` 开始下降到 `0.86957`
- 因此 `vmin=3.0, vmax=63.0` 被认定为当前更优的默认参数组合

### 仓库结构改造

- 删除 `survey_local`、`known_asteroid_local`、`heliolincrr_local`
- 将 `survey_server`、`known_asteroid_server`、`heliolincrr_server`
  分别重命名为：
  - `survey`
  - `known_asteroid`
  - `heliolincrr`
- 仓库从“双版本并存”切换为“只保留服务器对齐版本”
- 更新根目录 `README.md`、`AGENTS.md`、`.gitignore`
  - 默认工作流改为服务器运行
  - `known_asteroid` 大依赖文件改为按服务器布局忽略
  - 删除对本地运行目录的说明

## 2026-04-07

### 项目初始整理

- 建立 `survey_local` / `survey_server`
- 建立 `known_asteroid_local` / `known_asteroid_server`
- 建立 `heliolincrr_local` / `heliolincrr_server`
- 明确本地调试版本和服务器备份版本的分工
- 补充项目级 `README.md`、`AGENTS.md`、`WORKLOG.md`、`PLAN.md`

### heliolincrr make tracklet 改进

- `pair_two_exposures()` 不再只保留第二帧中的最近邻，改为：
  - 使用 `search_around_sky()`
  - 根据两帧 `EXPSTA + 15s` 的实际曝光中点时间差计算允许分离范围
  - 保留所有满足 `vmin / vmax / dmag-max` 的候选配对
- 曝光分组中心从 `CRVAL1 / CRVAL2` 改为优先使用 `CEN_RA / CEN_DEC`
- 增加 `--skip-common-area / --no-skip-common-area`
  - 默认跳过 `common_mask_from_group_reproject()` / fallback
  - 分组后直接进入静止源扣除与两点配对
- `r-static` 默认值固定为 `2.0`
- `dmag-max` 默认值固定为 `1.0`
- 在静止源扣除之前，新增边缘壳层伪源过滤：
  - 仅删除最近边缘距离在 `300` 到 `500` 像素之间
  - 且 `Flag_ISO_Num > 0` 的 detection

### heliolincrr 输出与验证辅助

- 新增 `tracklet_completeness_purity.py`
  - 用于只基于 nightly tracklets 和 `matched_asteroids.fits`
  - 统计已知小行星完备度与 tracklet 纯度
- 修复 `merge_tracklets_night.py`
  - 临时输出文件补上 `.fits` 扩展名
  - 避免 `astropy` 无法自动识别输出格式

### 20260220 单夜基准

- 纯 `make tracklet` 基准目录：
  - `/pipeline/xiaoyunao/data/heliolincrr/20260220/tracklets_linreproj_tracklet_only`
- 加入边缘壳层伪源过滤后的对比目录：
  - `/pipeline/xiaoyunao/data/heliolincrr/20260220/tracklets_linreproj_tracklet_only_edgeiso`
- 主要结果变化：
  - 两点 tracklet 总数：`31925 -> 17205`
  - “两端点都是已知小行星”的纯度：`0.04448 -> 0.08184`
  - 按已知小行星个数计算的完备度：`0.8804 -> 0.8772`
  - 按已知探测点数计算的完备度：`0.9030 -> 0.8996`

### known_asteroid 服务器备份同步

- 同步补回服务器当前存在但本地备份缺失的辅助脚本：
  - `plot_known_asteroids.py`
  - `run_visual_daily.sh`
  - `update_all_matched_history.py`
  - `cron_visual.example`
- 更新 `known_asteroid/README.md`
  - 纠正当前服务器目录与历史顺序脚本不一致的描述
