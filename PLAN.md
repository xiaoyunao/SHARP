# PLAN

## Current objective

优先处理 `known_asteroid` 的 MPC 重复上报问题，先堵住重复提交，再回到
`heliolincrr` 的 RR 扫描。

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

- 历史提交中确认存在跨夜次重复上报；`20251225` 与 `20251226` 至少共享大量相同 `(object, obsTime)` 观测
- 最近提交 `20260331` 到 `20260407` 未发现同类重复
- 已定位更上游根因：`submit_pipeline_slurm.sh` 旧逻辑按目录全量吞 `L2/*MP*`，不检查 FITS 实际观测时间；`20251226/L2` 中有 `166` 个文件经本地 observing-night 重判后应归 `20251225`
- manifest 源头修复已经写入并同步到服务器，但还需要用完整流程 rerun 验证对历史夜次的实际输出
- `export_ades.py` 的历史去重现在只是最后一道保险，不应替代源头夜次切分修复
- 需要继续追查上游 `/processed1/YYYYMMDD/L2` 的夜次切分或落盘规则，解释为什么会出现同名文件和相同 `DATE-OBS` 同时进入相邻夜目录
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

- `known_asteroid` 导出时能自动过滤：
  - 当前批次内同一 `(object, obsTime)` 的重复行
  - 历史 `*_matched_asteroids_ades.psv` 中已提交过的观测
- `submit_pipeline_slurm.sh` 生成 manifest 时，只会纳入实际属于目标 observing night 的 `L2` 文件
- 对历史问题夜次做 dry-run / rerun 时，输出日志能明确给出被历史去重丢弃的行数
- 对历史问题夜次做 manifest dry-run 时，`wrong_night` 文件数应与原相邻夜次重叠文件数一致
- 最近正常夜次在启用去重后仍能顺利导出，不会把整夜新观测误删
- `20260220` 能稳定从空目录重跑 `mask_gaia` 和纯 `make tracklet`
- 每次调参后都能产出 completeness / purity / tracklet 总数
- baseline 结果已记录，可和后续参数组合直接比较
- 已知小行星速度分布已整理成可直接参考的 summary / plot
- RR 调参后，除 link 总数外，还要能比较已知检测进入 RR 的覆盖率与歧义度
- 新 RR 代码需要稳定写出 `rr_summary.json`，`compare_with_known_asteroids.py` 需要稳定写出 FITS 和 JSON summary
- `known_asteroid` 每日 09:00 提交后，finalize 成功结束时能自动触发历史更新和出图

## Next recommended steps

1. 用已同步到服务器的 manifest 源头修复，对 `20251226` 做完整但不提交的 rerun
2. 核查 rerun 后的 `matched_asteroids.fits` / `ades.psv` 是否不再包含 `20251225` 夜的重复观测
3. 继续排查 `/processed1` 上游夜次目录的切分/落盘规则，解释相邻夜目录为何共享同名文件
4. `known_asteroid` 重复上报风险解除后，再恢复 `heliolincrr` 的 RR / 轨道拟合调试
