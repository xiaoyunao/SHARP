# PLAN

## Current objective

在单一服务器版本目录结构下，围绕 `heliolincrr` 的单夜 `rr link`
参数做扫描，并用 `20260220` 持续跟踪：

- 已知检测进入 RR 的覆盖率
- 每个 detection 对应的 RR 歧义度
- RR linkage 总数与成员规模

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
- `ref-dt-days=0.1` 已试跑，结果退回到高召回高歧义状态：`n_links=29554`、`rr_given_tracklet=2287/2304=99.26%`、`p90=48`
- 当前服务器默认 `rr_links` 目录中的 `0.05` 统计文件对应 `8342 / 2057 / p90=6`，与工作日志里记录的 `422a2dd` 稳定基线 `8393 / 2082 / p90=6` 不一致，需要避免后续比较时混淆基线
- 当前单夜 RR 参数方向已经切为“优先高召回，让下一步轨道拟合压噪声”，因此 `tol=0.03` 比 `0.02` 更符合当前目标
- 在 `tol=0.03` 下，`k-neighbors-cap=100` 会把 `rr_given_tracklet` 从 `96.31%` 压到 `84.16%`，说明这个值对当前目标过小
- 在 `tol=0.03` 下，`k-neighbors-cap=300` 会把 `rr_given_tracklet` 从 `96.31%` 提高到 `97.53%`，代价是 `p90` 从 `8` 升到 `9`
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

1. 继续保留 `ref-dt-days=0.05`
2. 单夜 RR 当前工作基线切到 `tol=0.03, k-neighbors-cap=300`，对应 `n_links=12274`、`rr_given_tracklet=2247/2304=97.53%`、`median=4`、`p90=9`
3. 后续默认以 `ref-dt-days=0.05, tol=0.03, k-neighbors-cap=300` 继续 RR / 轨道拟合调试
4. 下一步重点转到轨道拟合阶段，验证更高召回带来的额外候选能否被有效筛掉
