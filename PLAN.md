# PLAN

## Current objective

恢复 `heliolincrr` 的 RR link 调参，并以已经确认的单夜高召回基线继续往下做
参数扫描和后续轨道拟合衔接。

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

- `known_asteroid` 重复上报问题已经定位并完成程序侧保护，当前不再作为主线阻塞
- `20260220` 的 RR 单夜新基线已切到 `max-v-kms=30`，当前默认参数结果是 `8697 / 27204 / 2259/2304 / p90=8`
- `max-v-kms=30` 相比 `200` 已明显改善规模和高歧义尾部，下一轮主线改为扫描 `min-init-earth-au`
- `min-init-earth-au=0.02` 与当前 `0.01` 基线结果完全一致，当前默认值统一保持 `0.01`
- RR 新逻辑已在服务器上复跑验证，当前问题变成“歧义显著下降，但已知检测召回也有明显回落”
- 本轮算法微调尝试（mutual-neighbor、compactness gate、全局 prune）都未优于 `422a2dd` 稳定版，现阶段先停止算法层修改
- `ref-dt-days=0.1` 已试跑，结果退回到高召回高歧义状态：`n_links=29554`、`rr_given_tracklet=2287/2304=99.26%`、`p90=48`
- 当前服务器默认 `rr_links` 目录中的 `0.05` 统计文件对应 `8342 / 2057 / p90=6`，与工作日志里记录的 `422a2dd` 稳定基线 `8393 / 2082 / p90=6` 不一致，需要避免后续比较时混淆基线
- 当前单夜 RR 参数方向已经切为“优先高召回，让下一步轨道拟合压噪声”，因此 `tol=0.03` 比 `0.02` 更符合当前目标
- 在 `tol=0.03` 下，`k-neighbors-cap=100` 会把 `rr_given_tracklet` 从 `96.31%` 压到 `84.16%`，说明这个值对当前目标过小
- 在 `tol=0.03` 下，`k-neighbors-cap=300` 会把 `rr_given_tracklet` 从 `96.31%` 提高到 `97.53%`，代价是 `p90` 从 `8` 升到 `9`
- 需要重新确认后续服务器实验是否统一基于仓库当前 `run_rr_from_tracklets.py`，避免再混入别的脚本版本
- 当前服务器 `run_rr_from_tracklets.py` 已重新同步到仓库版本，当前默认单夜 `max-v-kms` 已定为 `30`
- `min-init-earth-au` 的当前有效边界还没找到；`0.01 -> 0.02` 没有产生任何变化，后续若继续扫应直接跳更大的值
- 需要判断是否补充更细的 group / exposure-pair 级诊断
- `known_asteroid/astorb.dat` 和 `known_asteroid/de432s.bsp` 不应进入 git
- 后续代码修改要持续与服务器目录保持一致
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

1. 以新的 `/pipeline/xiaoyunao/data/heliolincrr/20260220/rr_links` `max-v-kms=30` 结果作为唯一单夜基线
2. 若继续扫描 `min-init-earth-au`，直接从比 `0.02` 更大的点开始
3. 对每个 `min-init-earth-au` 结果统一比较 `n_links`、`n_member_rows`、`rr_given_tracklet`、命中 RR 的 `n_rr_links median/p90`、总运行时间
4. 若 `min-init-earth-au` 没有进一步收益，再决定是否继续进轨道拟合验证
