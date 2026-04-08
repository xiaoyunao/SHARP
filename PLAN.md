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

- `mask_gaia + make tracklet` 的 clean rerun 已完成，后续优化重点应转到 `rr link`
- RR 单夜基线已跑出，但当前已知检测对应的 RR 歧义度仍高，需要压低每个 detection 对应的候选 linkage 数
- RR 新逻辑已在服务器上复跑验证，当前问题变成“歧义显著下降，但已知检测召回也有明显回落”
- 本轮算法微调尝试（mutual-neighbor、compactness gate、全局 prune）都未优于 `422a2dd` 稳定版，现阶段先停止算法层修改
- 需要判断是否补充更细的 group / exposure-pair 级诊断
- `known_asteroid/astorb.dat` 和 `known_asteroid/de432s.bsp` 不应进入 git
- 后续代码修改要持续与服务器目录保持一致
- `known_asteroid` 的 finalize 自动出图链路刚恢复，需要确认明天自动任务是否正常衔接
- 服务器热修如果未先回写仓库，后续同步可能再次覆盖运行期修复

## Validation criteria

- `20260220` 能稳定从空目录重跑 `mask_gaia` 和纯 `make tracklet`
- 每次调参后都能产出 completeness / purity / tracklet 总数
- baseline 结果已记录，可和后续参数组合直接比较
- 已知小行星速度分布已整理成可直接参考的 summary / plot
- RR 调参后，除 link 总数外，还要能比较已知检测进入 RR 的覆盖率与歧义度
- 新 RR 代码需要稳定写出 `rr_summary.json`，`compare_with_known_asteroids.py` 需要稳定写出 FITS 和 JSON summary
- `known_asteroid` 每日 09:00 提交后，finalize 成功结束时能自动触发历史更新和出图

## Next recommended steps

1. 以稳定版 `422a2dd` 的 RR 基线结果 `n_links=8393`、`rr_given_tracklet=2082/2304=90.36%`、`p90(n_rr_links_hits_only)=6` 作为参数扫描起点
2. 围绕 `tol`、`k-neighbors-cap`、`ref-dt-days` 做单维扫描，优先观察是否能降低歧义度而不明显损失已知检测覆盖率
3. 若参数扫描后仍不满意，再重新设计更高效的 candidate / prune 架构，而不是继续在当前 `cluster_one()` 上打补丁
4. 核对 `survey` 和 `known_asteroid` 的 09:00 自动任务是否按新时间运行
