# PLAN

## Current objective

将 `rr_links_compactcore_softsky` 固化为唯一正式单夜 RR 基线，并继续提升单夜
`orbit_confirm_links.py` 对污染 links 的鲁棒轨道拟合能力。

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

- 当前 orbit fitting 的主要问题先不是阈值调参，而是大量 `fit_ok=False` link 直接停在“未产出轨道”阶段；这些行的 `rms_arcsec=inf`
- `orbit_confirm_links.py` 已补充 `fail_reason` 与 `fail_counts` 字段，可用于区分失败发生在 `geod2heliod`、`min_init_earth_au`、`lambert`、`max_v_kms`、传播还是 outlier clip
- `20260220` 单夜 orbit fitting 失败主因已确认高度集中在 `max_v_kms=30`：整体失败 links 中 `8435/8560` 主失败于 `max_v`，known-hit 失败 links 中 `2707/2719` 主失败于 `max_v`
- known-hit 失败 links 的累计 `fail_counts` 为 `max_v=21619`、`outlier_clip=133`，说明当前单夜 `8` 个 hypo 对大多数已知小行星 link 几乎都在 Lambert 后被速度上限统一卡掉
- 但 `max_v` 失败速度分布显示，这不主要是 `30 km/s` 略偏紧的问题：known-hit `max_v` 失败 links 的 `min_rejected_max_v_kms` 只有极少数贴近阈值，主体分布为 `p10=172 km/s`、`median=567 km/s`
- known-hit 且 `fit_ok=True` 的 links，其 `best_v1_kms` 上界只有 `28.91 km/s`；当前可用解与失败解在速度上明显分成两群
- 只看“全部观测都属于同一个已知目标”的 pure-known links 时，当前 orbit fitting 已经 `76/76` 全部拟合成功
- 因此后续若要提升 known-hit 的 orbit fitting 命中，重点不在继续优化 pure-link 的种子，而在为带污染的 RR links 增加 subset / consensus 式鲁棒拟合
- RR 端已试做 compact-core clustering：`pure_same_object` 从 `76` 提到 `95`，`partial_known_mixed_objects` 从 `813` 降到 `377`
- 但当前 compact-core 版本过于激进，known-tracklet 覆盖从 `2259/2304` 降到 `1996/2304`，`rr_given_tracklet` 从 `98.05%` 降到 `86.63%`
- `all_known_mixed_objects` 几乎未变（`244 -> 243`），说明仅靠连通分量后处理还不足以解决“多个 known 对象被合并”的问题
- compact-core 放宽版已试跑：`pure_same_object=89`、`all_known_mixed_objects=235`、`partial_known_mixed_objects=565`，同时 `rr_given_tracklet` 回升到 `90.06%`
- 这说明 compact-core 路线可调，但当前仍未达到兼顾 purity 与 recall 的可接受折中
- 第三轮继续放宽后，`rr_given_tracklet` 只升到 `91.28%`，但 purity 开始回退到 `pure_same_object=82`、`partial_known_mixed_objects=613`
- 因此单纯继续放宽几何 compact-core 阈值已接近收益递减，当前最好折中点仍是第二轮 `rr_links_compactcore_relaxed`
- motion-consistency 版证明了更物理的约束是有效方向：`pure_same_object=94`、`partial_known_mixed_objects=337`
- 但当前 `motion_limit` 过严，导致 `rr_given_tracklet` 掉到 `84.64%`；现阶段仍不如 `rr_links_compactcore_relaxed` 作为折中点
- soft sky-clipping 版在 `rr_links_compactcore_relaxed` 基础上带来了小幅改进：`rr_given_tracklet=90.63%`，`partial_known_mixed_objects=537`
- 但 `pure_same_object` 仍停在 `89`，`all_known_mixed_objects` 还略有回升到 `238`，因此它只是小幅优于 `rr_links_compactcore_relaxed`，还不是明显的新基线
- `known_asteroid` 重复上报问题已经定位并完成程序侧保护，当前不再作为主线阻塞
- `20260220` 的 RR 单夜基线已重新清空旧结果并重跑，唯一默认结果固定为 `8697 / 27204 / 2259/2304 / p90=8`
- 当前服务器 `run_rr_from_tracklets.py` 已同步到仓库版本，新的 `rr_links` 基线不再混入旧脚本结果
- 当前 RR 单夜参数固定为 `ref-dt-days=0.05`、`tol=0.03`、`k-neighbors-cap=300`、`max-v-kms=30`、`min-init-earth-au=0.02`
- 当前单夜默认 hypo 已固定为只扫 `r=[1.3, 1.7, 2.1, 2.5, 2.9, 3.3, 3.7, 4.1]`，并固定 `rdot=0.0`、`rddot=0.0`
- `run_rr_from_tracklets.py` 已分离出两组 profile 默认值：`single-night` 与 `w15`
- 当前正式单夜 RR 基线已切换为 `rr_links_compactcore_softsky` 升格后的 `rr_links`：`7900 / 21325 / 2088/2304 / rr_given_tracklet=90.63%`
- 当前单夜 orbit fitting 默认物理参数已统一到 RR 单夜默认：`min-init-earth-au=0.02`、`max-v-kms=30`、`hypos=r=[1.3..4.1], rdot=0, rddot=0`
- 单夜 orbit fitting 现已接入 tracklet-subset seed + 全 link inlier 重拟合流程
- 新增 relaxed-seed 变体：`seed_max_v_kms` 与 `final_max_v_kms` 分离，当前测试为 `120 -> 30`
- 新版单夜 orbit 输出已带 tracklet 诊断列：`seed_tracklets`、`inlier_tracklets`、`rejected_tracklets`
- 在新的正式 `rr_links` 上，当前最新 orbit fitting 结果为：`fit_ok=263/7900`、`is_good=110/7900`
- 对应的已知小行星命中统计当前仍为：`fit_ok_any=438`、`is_good_any=272`
- 已知命中 link 级统计当前为：`known_hit_links=2421`、`fit_ok_links=222`、`is_good_links=94`
- 新测试下失败主因仍以 `max_v` 为绝对主导：`7532` 个失败 link 主失败于 `max_v`，`105` 个主失败于 `outlier_clip`
- 这说明“放宽 seed `max_v` + 当前 tracklet 级 residual 聚合规则”并未继续提升 known-hit orbit 命中
- `final_max_v_kms=100` 的隔离实验只带来极弱增益：`fit_ok_any 438 -> 440`、`is_good_any` 仍为 `272`
- 同时更宽的最终速度上限已经引入更高速度尾部，`fit_ok` 样本 `best_v1_kms p99` 升到 `73.34 km/s`
- 因此暂不建议把正式 `final_max_v_kms` 直接切到 `100`
- 新增 `orbit_fit_stats.py`，当前单夜最终结果统计已写出 `/pipeline/xiaoyunao/data/heliolincrr/20260220/analysis/20260220_orbit_fit_stats.json`
- 服务器正式路径 `/pipeline/xiaoyunao/heliolincrr/run_rr_from_tracklets.py` 只允许放“已入仓库、已提交”的版本
- 任何 RR 试验改动都不再直接修改服务器正式脚本；临时实验统一使用副本，例如 `/tmp/run_rr_from_tracklets_<tag>.py`
- 每次正式重跑基线前，先比较仓库版和服务器正式脚本的 SHA；若不一致，先同步仓库，再跑正式基线
- 若服务器热修被证明有效，必须先回写仓库并提交，再允许覆盖服务器正式脚本
- RR 新逻辑已在服务器上复跑验证，当前问题变成“正式 RR purity/歧义有所改善，但已知检测召回仍明显低于高召回老基线”
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

1. 以当前 `/pipeline/xiaoyunao/data/heliolincrr/20260220/rr_links` 作为唯一正式单夜 RR 基线，不再回切其它 RR 试验目录
2. 基于 `seed_tracklets / inlier_tracklets / rejected_tracklets` 诊断，先分析 known-hit links 中哪些污染模式仍未被剔掉
3. 暂不把“继续放宽 seed/final `max_v_kms`”作为主线；若改 orbit fitting，优先调整 tracklet 级评分/筛选规则，而不是再放物理上限
4. 若 robust orbit fitting 继续提升有限，再回到 RR 端评估更细的对象级拆分或软评分裁剪
5. 若后续继续处理 15 夜 RR，统一基于 `w15` profile 单独维护参数与实验记录，不再借用单夜默认值
6. 若继续测试 RR `hypo` 网格，保持正式脚本不动，统一通过 `--hypos <tmpfile>` 或 `/tmp/run_rr_from_tracklets_<tag>.py` 做隔离实验
