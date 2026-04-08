# PLAN

## Current objective

在单一服务器版本目录结构下，围绕 `heliolincrr` 的 `make tracklet`
环节做单夜调参，并用 `20260220` 持续跟踪：

- 已知小行星 completeness
- tracklet purity
- tracklet 总数

## Milestones

1. 完成仓库目录改造，只保留 `survey` / `known_asteroid` / `heliolincrr`
2. 改正 `pair_two_exposures()` 为全候选配对
3. 允许跳过 common-area，默认打开
4. 用 `CEN_RA/CEN_DEC` 做曝光分组
5. 清理 `20260220` 旧测试目录，仅保留 `mask_gaia`
6. 在 `tracklets_linreproj_tracklet_only` 上生成纯 tracklet baseline 统计
7. 加入边缘壳层 + `Flag_ISO_Num` 过滤并形成新基线
8. 基于已知小行星速度分布继续扫描 `vmin` / `vmax`

## Outstanding issues

- 仍需继续调 `vmin` / `vmax`
- 需要判断是否补充更细的 group / exposure-pair 级诊断
- `known_asteroid/astorb.dat` 和 `known_asteroid/de432s.bsp` 不应进入 git
- 后续代码修改要持续与服务器目录保持一致
- 需要基于 `20260220_known_motion_summary.json` 缩小 `vmin` / `vmax` 扫描范围

## Validation criteria

- `20260220` 能稳定重跑纯 `make tracklet`
- 每次调参后都能产出 completeness / purity / tracklet 总数
- baseline 结果已记录，可和后续参数组合直接比较
- 已知小行星速度分布已整理成可直接参考的 summary / plot

## Next recommended steps

1. 固定 `20260220` 为对比样本，先扫 `vmax`
2. 可优先测试 `vmax=40, 50, 60, 80`
3. 在确定 `vmax` 后，再扫 `vmin=0.0, 0.3, 0.5, 1.0`
4. 以 `tracklets_linreproj_tracklet_only_edgeiso` 作为当前比较基线
5. 若服务器端目录有变更，优先同步到仓库正式目录而不是另建副本
