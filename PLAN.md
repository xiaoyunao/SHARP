# PLAN

## Current objective

将 `heliolincrr` 的 RR 单夜基线固定为唯一版本，并继续衔接后续轨道拟合
验证。

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
- `20260220` 的 RR 单夜基线已重新清空旧结果并重跑，唯一默认结果固定为 `8697 / 27204 / 2259/2304 / p90=8`
- 当前服务器 `run_rr_from_tracklets.py` 已同步到仓库版本，新的 `rr_links` 基线不再混入旧脚本结果
- 当前 RR 单夜参数固定为 `ref-dt-days=0.05`、`tol=0.03`、`k-neighbors-cap=300`、`max-v-kms=30`、`min-init-earth-au=0.02`
- 当前单夜默认 hypo 已固定为只扫 `r=[1.3, 1.7, 2.1, 2.5, 2.9, 3.3, 3.7, 4.1]`，并固定 `rdot=0.0`、`rddot=0.0`
- `run_rr_from_tracklets.py` 已分离出两组 profile 默认值：`single-night` 与 `w15`
- 当前单夜 orbit fitting 基准已生成：`fit_ok=146/8697`、`is_good=97/8697`
- 已知小行星在当前单夜 orbit fitting 基准上的命中统计为：`fit_ok_any=314`、`is_good_any=242`
- 当前单夜 orbit fitting 默认物理参数已统一到 RR 单夜默认：`min-init-earth-au=0.02`、`max-v-kms=30`、`hypos=r=[1.3..4.1], rdot=0, rddot=0`
- 用统一后的单夜 orbit fitting 默认值重跑后，新的基准更新为：`fit_ok=137/8697`、`is_good=96/8697`
- 对应的已知小行星命中统计更新为：`fit_ok_any=312`、`is_good_any=242`
- 新增 `orbit_fit_stats.py`，当前单夜最终结果统计已写出 `/pipeline/xiaoyunao/data/heliolincrr/20260220/analysis/20260220_orbit_fit_stats.json`
- 服务器正式路径 `/pipeline/xiaoyunao/heliolincrr/run_rr_from_tracklets.py` 只允许放“已入仓库、已提交”的版本
- 任何 RR 试验改动都不再直接修改服务器正式脚本；临时实验统一使用副本，例如 `/tmp/run_rr_from_tracklets_<tag>.py`
- 每次正式重跑基线前，先比较仓库版和服务器正式脚本的 SHA；若不一致，先同步仓库，再跑正式基线
- 若服务器热修被证明有效，必须先回写仓库并提交，再允许覆盖服务器正式脚本
- RR 新逻辑已在服务器上复跑验证，当前问题变成“歧义显著下降，但已知检测召回也有明显回落”
- 本轮算法微调尝试（mutual-neighbor、compactness gate、全局 prune）都未优于 `422a2dd` 稳定版，现阶段先停止算法层修改
- `ref-dt-days=0.1` 已试跑，结果退回到高召回高歧义状态：`n_links=29554`、`rr_given_tracklet=2287/2304=99.26%`、`p90=48`
- 当前单夜 RR 参数方向已经切为“优先高召回，让下一步轨道拟合压噪声”，因此 `tol=0.03` 比 `0.02` 更符合当前目标
- 在 `tol=0.03` 下，`k-neighbors-cap=100` 会把 `rr_given_tracklet` 从 `96.31%` 压到 `84.16%`，说明这个值对当前目标过小
- 在 `tol=0.03` 下，`k-neighbors-cap=300` 会把 `rr_given_tracklet` 从 `96.31%` 提高到 `97.53%`，代价是 `p90` 从 `8` 升到 `9`
- 当前单夜 RR 参数已经固定，不再继续做这组参数扫描
- 需要判断是否补充更细的 group / exposure-pair 级诊断
- 后续若继续调参，应把单夜与 15 夜完全分开记录和比较，避免默认值再次混淆
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

1. 以新的 `/pipeline/xiaoyunao/data/heliolincrr/20260220/rr_links` 结果作为唯一 RR 单夜基线
2. 以当前 `rr_links/orbit_confirm` 结果作为单夜 orbit fitting 当前基准，开始系统整理其参数含义、默认值和优先调参顺序
3. 若后续继续处理 15 夜 RR，统一基于 `w15` profile 单独维护参数与实验记录，不再借用单夜默认值
4. 若继续测试 RR `hypo` 网格，保持正式脚本不动，统一通过 `--hypos <tmpfile>` 或 `/tmp/run_rr_from_tracklets_<tag>.py` 做隔离实验
5. 以 `20260220_orbit_fit_stats.json` 为基础，优先找出最影响 `fit_ok/is_good` 的分桶和 residual 指标，再决定首轮 orbit fitting 调参顺序
