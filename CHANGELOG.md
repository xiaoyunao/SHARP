# CHANGELOG

本文件记录仓库中值得保留的版本级行为变化，重点覆盖服务器基线流程、
调试中确认过的默认值变更，以及会影响结果解释的输出变化。

## 2026-04-08

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
