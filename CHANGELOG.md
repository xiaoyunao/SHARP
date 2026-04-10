# CHANGELOG

本文件记录仓库中值得保留的版本级行为变化，重点覆盖服务器基线流程、
调试中确认过的默认值变更，以及会影响结果解释的输出变化。

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
