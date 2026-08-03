# SMT 黄道面小行星巡天项目：论文讨论用项目与结果总结

> 版本：2026-08-03 服务器与本地材料快照
>
> 用途：团队讨论论文基调、定位、结果口径和后续补充工作。本文不是论文正文，也不是最终数据发布说明。
>
> 核心原则：把“程序已实现”“实际已观测”“已生成上报文件”“收到 MPC submission ID”“后验高可信样本”严格分开。

仓库快照为本地 `main` 的 `d2f0057`。本次逐一比较了 17 个 survey/known/unknown 关键入口和算法文件的 SHA-256，均与服务器正式路径一致。

## 1. 一页结论

本项目已经形成一条在兴隆站 60/90 cm Schmidt 望远镜上实际运行的端到端小行星巡天链路：

1. 根据黄道附近预定义天区、可见性、月距、历史覆盖和近期观测记录，自动生成逐夜观测计划；
2. 接收同事完成预处理、定标和测光后产生的 L1/L2 数据；
3. 对 L2 星表并行计算已知小行星星历、完成天球匹配、生成 ADES，并向 MPC 提交；
4. 在扣除 Gaia 恒星和已知小行星后，构造单夜 tracklet、连接候选链、做轨道一致性确认；
5. 为 unknown 候选分配可追溯的 `trkSub`，生成 GIF 和网页人工复核包，只上报人工标记为真的候选；
6. 已开发 reviewed unknown 的自动 follow-up 排程、观测结果核对和跨夜关联状态机。

截至本次快照，最有说服力的定量结果是：

| 层次 | 当前结果 | 推荐论文表述 |
|---|---:|---|
| 标准巡天原始曝光 | 41,074 帧，134 夜，30 s/帧 | 342.28 h 有效快门时间 |
| 标准巡天天区 | 1,430 个 | HEALPix 去重覆盖约 10,387.43 deg² |
| 全部 MP FITS（含 CAL/Tel 等） | 41,152 帧，136 夜，7.735 TB | 只用于数据工程统计，不作为纯科学曝光数 |
| 当前 L2 可处理数据 | 40,399 个 MP 星表，131 夜 | known/unknown 的实际输入规模 |
| 已知天体匹配 | 534,780 次，58,482 个对象，130 夜 | 1″匹配产物；不是全部都上报 |
| 已知天体 ADES 生成 | 350,706 行，29,498 个对象，127 个文件 | “生成的 ADES”口径 |
| 有 MPC reply 的 known 提交 | 348,568 行，29,488 个对象，122 夜 | 服务器收到 submission ID；不等同最终接收/入库 |
| unknown 自动候选 | 4,762 条单夜 link，14,299 个 detection 成员 | 116 个非空夜，另有 8 个零结果夜 |
| 人工初次通过并提交 | 67 条 link，206 次观测，35 夜 | 已收到 35 个 unknown submission ID |
| 后验再次清洗 | 58 条高可信 link，179 次探测，34 夜 | 论文科学分析应使用这一口径 |
| 自动 follow-up 实际执行 | 0 个源、0 夜、0 帧 | 程序已部署，但尚无真实闭环观测结果 |

最稳妥的论文主线不是“发现了 58 颗新小行星”，而是：

> 建立并实地运行了一套面向广域黄道巡天的自动调度、已知天体天测上报和未知移动源搜索系统；用约 4.1 万帧观测验证了系统的工程可运行性，产出大规模已知小行星天测，并给出一批经过人工复核的单夜未知移动源候选。

原因是当前 58 条高可信记录仍然是“单夜 link”，没有完成统一的跨夜关联、轨道精化、MPC designation/identification 核验和本项目 follow-up。论文中应称为 `high-confidence unknown moving-object linkages/candidates`，暂不称为 58 个 confirmed discoveries。

## 2. 推荐论文定位

### 2.1 推荐主定位

推荐定位为“巡天系统 + 运行表现 + 初步科学产出”的方法与数据论文，贡献按以下顺序展开：

1. **动态黄道巡天调度**：不是固定 RA 表，而是按逐夜天文条件和历史覆盖动态选场；
2. **端到端自动化运行**：计划、known、unknown、人工复核、MPC 上报和断电恢复串成统一 daily workflow；
3. **大规模已知小行星天测**：58,482 个已知对象、534,780 次 1″匹配，展示小口径 Schmidt 的广域天测产能；
4. **单夜未知移动源搜索**：给出 4,762 候选到 58 条高可信 link 的完整筛选漏斗和假源问题；
5. **follow-up 框架**：作为已经开发的闭环能力和下一阶段，而不是当前论文的观测结果。

这个定位的优势是证据最完整。known 部分有大样本统计和 MPC 上报链路；unknown 部分有透明的候选、人工复核和后验清洗过程；follow-up 尚未实际发生，不会被迫过度宣称。

### 2.2 不建议作为当前主定位的方向

- **“58 颗新小行星发现”**：当前没有 58 个独立天体的跨夜去重和确认，也没有逐目标 MPC 最终识别/定名状态。
- **“精确轨道分布研究”**：55/58 条高可信 link 只有 3 次单夜探测，短弧轨道中的 `a/e/i` 有明显退化，甚至出现负半长轴或 `e > 1`；轨道解只适合作为一致性过滤。
- **“follow-up 成果论文”**：自动 follow-up 代码已完成，但 state 仍为空，没有 `MP_FU_...` 真实计划或数据。
- **“完备性/发现效率测量”**：目前没有注入恢复实验、控制样本或按星等/速度计算的 completeness，不能只用 58/4,762 推断巡天完备性。

### 2.3 可讨论的英文标题方向

- *An End-to-End Ecliptic Small-Body Survey with the Xinglong 60/90 cm Schmidt Telescope*
- *Automated Scheduling, Known-Asteroid Astrometry, and Single-Night Moving-Object Searches in a Wide-Field Ecliptic Survey*
- *Operational Performance of a Wide-Field Asteroid Survey at Xinglong Station*

## 3. 系统和数据流总览

```text
黄道面 footprint 库
        ↓
survey 动态 scheduler ──→ 当晚观测脚本/current_plan.txt
        ↓                           ↓
      实际观测 ──→ raw images ──→ 同事的预处理/定标/测光
                                      ↓
                                  L1 / L2 catalogs
                                      ↓
                 ┌────────────────────┴────────────────────┐
                 ↓                                         ↓
        known_asteroid                              heliolincrr unknown
   星历预测→1″匹配→ADES                    Gaia mask→tracklet→link→orbit check
          ↓                    1.5″ known mask             ↓
       MPC submit ─────────────────────────────→ unknown 扣除已知天体
                                                            ↓
                                              GIF/review package/网页人工 check
                                                            ↓
                                             reviewed ADES→MPC submit→follow-up
```

三个正式代码目录及服务器路径为：

| 模块 | 仓库目录 | 服务器目录 | 运行环境 |
|---|---|---|---|
| 巡天计划与随访 | `survey/` | `/pipeline/xiaoyunao/survey` | 默认 Python |
| 已知小行星 | `known_asteroid/` | `/pipeline/xiaoyunao/known_asteroid` | 默认 Python + Slurm |
| 未知小行星 | `heliolincrr/` | `/pipeline/xiaoyunao/heliolincrr` | `heliolinc` conda 环境 |

每日 09:00 的正式总入口为：

```bash
/pipeline/xiaoyunao/heliolincrr/run_daily_pipeline.sh
```

服务器 crontab 还配置了 `@reboot` 延迟 600 s 的恢复入口。总入口先生成当天晚上的 survey plan，再在 recovery lookback 窗口内补跑已有数据夜的 known 和 unknown。它不会凭空补拍关机期间错过的夜晚，只补处理实际存在的数据。

## 4. 黄道面巡天规划程序

### 4.1 天区库

正式 footprint 文件为：

```text
/pipeline/xiaoyunao/survey/footprints/survey_fov_footprints_with_visibility.fits
```

当前天区库包含 1,855 个场。每场四角给出的名义面积为 8.826 deg²，对应约 2.94° × 2.95° 的相机视场。天区中心覆盖：

- RA：0.438°–359.908°；
- Dec：−41.952°–+39.699°；
- 黄道纬度：−18.516°–+16.262°。

场表保存 `field_id`、中心/四角 RA/Dec 和 `visible_day(s)`。程序会补算 `field_index`、赤纬带、RA 带和黄道纬度。

### 4.2 调度约束

核心代码是 `survey/scheduler.py`、`survey/astro_utils.py` 和 `survey/run_daily.py`。默认参数为：

| 参数 | 值 | 含义 |
|---|---:|---|
| 站点 | 117.575° E, 40.393° N, 960 m | 兴隆站，MPC 327 |
| 夜晚边界 | Sun altitude < −12° | scheduler 使用航海曙暮光边界 |
| 单次曝光 | 30 s | 观测曝光时间 |
| 读出开销 | 15 s | 每帧 |
| 转场开销 | 35 s + 大角距附加项 | 相邻天区 |
| 最低高度角 | 30° | 可观测约束 |
| 稳定可见窗口 | 90 min | 每 10 min 采样 |
| 窗口内可见比例 | ≥ 50% | 同时满足高度角和月距 |
| 普通 cluster | 最多 45 场 | 围绕最高权重中心选最近邻场 |
| 普通重复次数 | 3 次 | 为移动目标构造时间基线 |
| 重复起始间隔 | 1,800 s | 约 30 min |
| 单个 repeat block | ≤ 1,800 s | 控制一轮长度 |
| 最近观测冷却 | 3 夜 | 避免连续重复同一批场 |
| 夜末触发 | 剩余 < 60 min | 转为 near-Sun 模式 |

月距阈值随月面照明比例二次插值：新月约 2°、半月约 25°、满月约 40°，另加 0.5° margin。每个候选场的权重由以下因素组成：

```text
1.0 × 年内可见天数稀缺度
+ 1.2 × 当前高度角
+ 2.5 × 历史曝光不足度
+ 1.8 × 最近观测冷却
+ 0.4 × 月距是否通过
```

因此该 scheduler 的科学意图是：优先当前可见、较高、长期可见窗口较窄、历史覆盖较少、最近没有重复观测的场。

### 4.3 夜间排程逻辑

普通模式中，程序先找最高权重场作为中心，再取稳定可见场中的最近 45 场，按 RA 跨 0° 连续排序，完成三轮重复。在夜晚最后不足 60 min 时，转入 `near_sun` 模式：按与太阳角距由小到大选择还能在剩余时间内完成重复的场。

输出包括：

```text
runtime/plans/YYYYMMDD_plan.json
runtime/plans/YYYYMMDD_plan.txt
survey/plots/YYYYMMDD/*.png
survey/plots/YYYYMMDD/*nightly_cycles.gif
/pipeline/xiaoyunao/script/current_plan.txt
```

每个计划行保存场号、cycle、repeat、模式、UTC 开始时间、中心 RA/Dec、曝光时间、月相和月距，适合直接追溯一次观测为何被选中。

### 4.4 调度运行统计

服务器从 2026-03-29 到 2026-08-03 保存了 115 夜计划，共 37,619 个计划曝光：

- `normal`：35,576；
- `near_sun`：2,043；
- 使用过 769 个不同场；
- 计划的纯曝光时间合计 313.49 h。

为了只比较已经有原始数据的区间，取 2026-03-29–2026-07-15：

- 96 夜计划，31,639 个计划曝光；
- 其中 38 夜获得标准巡天原始数据，共 7,833 帧；
- 按“夜次 + field_id + 重复次数”的多重集合比较，7,455/7,833 = 95.17% 的实际曝光可在当夜计划中找到对应项；
- 每个有数据夜的匹配比例中位数为 100%，最低为 78.90%。

这里的 95.17% 衡量“实际拍到的帧有多少符合计划”，不衡量天气造成的计划损失。不能把所有计划曝光与实际曝光直接相除后解释为望远镜效率。

### 4.5 历史覆盖机制的技术说明

当前 `exposure_history.fits` 继承了旧调度阶段的累计覆盖状态，现版 daily runner 又会从 L2 文件增量 ingest 实际 `MP` 文件。现版 `build_plan()` 在单夜内部会临时增加 history，防止同一夜反复选场，但这一份 scheduler copy 不会在生成计划后写回磁盘。因此服务器现有 history 是面向选场优先级的历史/迁移状态，不是可直接解释的纯实际曝光台账。本文所有“实际观测帧数”均重新从 `/raw1` 或 `/processed1/L2` 文件统计，不能直接把 history 中的 `exposure_count` 求和当作实际观测数。

## 5. 已完成观测和数据规模

### 5.1 望远镜与探测器参数

当前 ADES header 中使用的仪器描述为：

- 60/90 cm Schmidt telescope；
- 有效口径 0.6 m；
- f/3；
- CCD，9216 × 9232 pixels；
- pixel scale 1.15″/pixel；
- unfiltered；
- 星表天测/测光参考 `Gaia3E`，上报 band 为 `G`；
- MPC observatory code：327。

### 5.2 原始数据统计

统计区间为 2025-11-15–2026-07-15；截至 2026-08-03 服务器没有更晚的原始观测夜。

`/raw1/YYYYMMDD/OBJ` 中：

- 全部文件：47,101 个，8.8538 TB；其中 FITS 为 47,095 个；
- 文件名含 `MP` 的 FITS：41,152 个，7.7351 TB，分布在 136 夜；
- 严格符合 `OBJ_MP_####_####.fits` 的标准巡天科学帧：41,074 个，分布在 134 夜；
- 标准科学帧的 30 s 快门时间合计 1,232,220 s = 342.28 h = 14.26 d。

41,152 与 41,074 的差额 78 个文件由以下非标准产品构成：

- `CAL_MP_...`：24 帧；
- `Tel_MP_...`：53 帧；
- `OBJ_MP_..._overscan_removed.fits`：1 个派生副本。

另有两个名字含 `MP` 的 `.tab` 文件，它们不是曝光，未计入 41,152 帧。

### 5.3 天区覆盖和访问频次

标准科学曝光覆盖 1,430 个 footprint：

- 单场累计曝光最少 1 次；
- 中位数 28 次；
- 最大 95 次；
- 所有名义 footprint 面积直接相加为约 12,621 deg²，但场间有重叠；
- 以前一次 HEALPix `nside=512` 的重叠去重统计为 **10,387.43 deg²**，论文应采用这一有效覆盖面积；
- 实际使用场中心的 Dec 为 −18.754°–+39.699°；
- 实际使用场中心的黄道纬度为 −18.516°–+16.262°。

### 5.4 L2 可用数据

在 `/processed1/YYYYMMDD/L2` 中，2025-11-15–2026-07-15 有：

- 131 个包含 MP 星表的夜次；
- 40,399 个 MP L2 catalog；
- 14,176 个“夜次–场”组合；
- 每夜 MP L2 数中位数 313，最大 627；
- 每夜不同场数中位数 100，最大 213。

134 个标准原始科学夜中，`20260209`、`20260312`、`20260413` 没有可用 MP L2；其他夜也可能存在部分帧未进入 L2。因此 raw 是观测完成口径，L2 是当前能进入小行星分析的口径。

## 6. 已知小行星分析与上报

### 6.1 程序结构

正式入口为 `known_asteroid/submit_pipeline_slurm.sh`。每个 L2 MP 文件作为一个 Slurm array task，推荐同时最多 24 个任务；所有 task 完成后由依赖的 finalize job 归并并串行导出/上报。

关键程序为：

| 程序 | 作用 |
|---|---|
| `build_file_manifest.py` | 从 FITS header 的观测时间重建夜次输入清单；北京时间中午前归入前一观测夜 |
| `match_single_night.py` | 对一个文件计算已知小行星星历并与 L2 源匹配 |
| `slurm_match_one_file.sh` | 一文件一任务，控制 CPU/内存和重试 |
| `merge_night_parts.py` | 合并 per-file 结果，写夜级 FITS 和 status JSON |
| `export_ades.py` | 质量过滤、去重、输出 ADES PSV、validate/submit |
| `slurm_merge_submit.sh` | finalize、ADES 和 MPC submit |
| `plot_known_asteroids.py` | 巡天覆盖、全天分布、夜间产额和动画/GIF |
| `update_all_matched_history.py` | 累积历史 matched 产品 |
| `audit_duplicate_causes.py` | 审计历史重复观测/重复提交来源 |
| `run_known_rematch_then_unknown_remask.sh` | known 强制重匹配完成后，依赖触发 unknown remask |

### 6.2 匹配方法

`match_single_night.py` 使用 Lowell `astorb.dat` 和 `aleph.Query`：

1. 从每个 L2 catalog 的 WCS 得到四角和中心；
2. 用曝光开始时间加半个曝光时间作为星历 epoch；
3. 使用兴隆站地理位置作为 topocentric observer；
4. 查询预测 `V < 22.5` 的候选；
5. 将预测位置投影回 CCD，只保留落在图像内的目标；
6. 与 L2 `RA_Win/DEC_Win` 最近邻匹配。

RA 跨 0° 的图像会额外以 RA=0° 和 359° 为中心查询并去重，修复了此前 Lowell/aleph 圆形查询在 RA wrap 处漏目标的问题。

正式输出分两种半径：

- `*_matched_asteroids.fits`：1.0″，用于 known 科学统计和 MPC 上报；
- `*_matched_asteroids_mask15.fits`：1.5″，只用于 unknown 阶段更保守地扣除已知源。

### 6.3 ADES 规则与防重复

known ADES 的关键规则包括：

- 每个对象每夜至少 3 次观测；
- 同批数据内按 `(object, obsTime)` 只保留质量最佳的一条；
- 与 `/processed1/*/L4/*_matched_asteroids_ades.psv` 历史比较，过滤已经上报过的 `(object, obsTime)`；
- 先 validate，真实 submit 在单个 finalize job 中串行进行；
- 保存 MPC response 和 submission ID。

增加历史去重是因为旧夜目录曾有跨相邻夜的重复观测，特别是 2025-12-25/26；近期流程已加保护。

### 6.4 已知天体统计结果

当前 1″ matched 产品覆盖 2025-11-15–2026-07-15 的 130 夜：

- 534,780 行 matched detection；
- 按编号或未编号名称去重得到 58,482 个已知对象；
- 其中 56,784 个有永久编号，1,698 个按名称/临时标识统计；
- 每对象探测次数：中位数 4，90% 分位 25，最大 90；
- 35,153 个对象至少 3 次；
- 18,902 个对象至少 10 次；
- 3,992 个对象至少 30 次。

位置匹配残差为预测位置与 `RA_Win/DEC_Win` 的球面小角距：

- 中位数：0.263″；
- 90% 分位：0.628″；
- 95% 分位：0.757″；
- 最大值受 1.0″匹配阈值限制，为 0.99997″。

亮度快照：

- 预测 V 星等中位数 18.80，10%–90% 为 17.28–19.69；
- 有效 `Mag_Aper4` 中位数 18.23，10%–90% 为 16.64–19.04。

ADES 有两种必须区分的统计：

1. **服务器生成文件口径**：127 个 `*_matched_asteroids_ades.psv`，350,706 行，29,498 个对象；
2. **有 known MPC reply/submission ID 的夜次口径**：122 夜，348,568 行，29,488 个对象，其中 `(object, obsTime)` 去重后 346,486 条。

论文如果写“submitted”，优先使用第二组数；如果写“ADES products generated”，使用第一组数。submission ID 只能证明服务器收到提交响应，不自动证明全部记录已经通过 MPC 后续 ingest 或进入最终数据库。

### 6.5 已知天体数据质量问题

- `20260111` 有 347 个 MP L2 输入和 106,042 条预测候选，但 1″ matched 仅 180、1.5″ matched 仅 419。抽样显示 L2/WCS 相对星历在 RA 方向有约 8–11″系统偏移；该夜 known/unknown 不宜作为可靠科学结果。
- 当前 matched 存在但没有 ADES 文件的夜次：`20251116`、`20260220`、`20260414`。
- 当前生成 ADES 但没有记录到正式 known reply 的夜次：`20251115`（0 行）、`20251201`、`20260111`、`20260412`、`20260503`。

正式论文统计前应建立统一的 quality mask，至少排除或单列 `20260111`，并解释上述夜次状态。

## 7. 未知小行星搜索、复核与上报

### 7.1 正式单夜流程

正式入口是 `heliolincrr/run_single_night.sh`。主流程如下：

```text
L2 photometry catalogs
  → Gaia 恒星 mask
  → 同场多曝光构造两点 tracklet
  → 合并全夜 tracklet
  → 共享端点的线性 link
  → heliocentric/Lambert 轨道一致性确认
  → 扣除 1.5″ known asteroid matched detections
  → unknown FITS/JSON
  → trkSub
  → GIF + review package
  → 网页人工 check
  → reviewed ADES + MPC submit
```

### 7.2 Gaia 和基础源过滤

`mask_gaia.py` 对每个 L2 catalog：

- 预选 `Mag_PSF ≤ 21`；
- 要求 `Flag == 0`；
- 从 `/pipeline/ref/healpix` 读取 Gaia 源；
- 将 Gaia proper motion 从 epoch 2016.0 推到曝光时刻；
- 同时检查 window 和 PSF 坐标；
- 在 1.5″内匹配 Gaia 的 L2 源被删除；
- 16 进程并行处理文件。

这一步的目的是先消除绝大多数恒星/静态源，而不是直接检测移动目标。

### 7.3 Tracklet 构造

`make_tracklet_linreproj.py` 的正式默认值：

- 角速度范围 3–63″/h；
- 两次探测星等差 `≤ 1.0 mag`；
- 静态源半径 2.0″；
- 至少两次重复曝光；
- 图像边缘剔除 `edge_pix=500`，另有 30 pixel erosion；
- 单 group 最多 100,000 tracklets，避免异常组占满内存；
- 支持空结果并输出 schema 完整的 FITS。

历史调参以 `20260220` 为基准。当前参数在该夜能保留约 87.4% 的 known-object tracklet completeness，但 tracklet purity 仍低，说明后续 link/orbit/人工复核不可省略。

### 7.4 单夜线性连接

正式生产流程已经从早期 RR 聚类试验转为 `run_linear_links_from_tracklets.py`：

- 两个 tracklet 必须共享一个 detection endpoint；
- 速度差不超过 5″/h；
- 方向差不超过 10°；
- 将兼容 tracklet 组合为单夜 linkage；
- 写出 `links_tracklets.fits`、`linkage_members.fits` 和连接边表。

共享端点约束显著降低了早期 RR linkage 的组合爆炸和同一 detection 进入大量 link 的歧义。

### 7.5 轨道一致性确认

`orbit_confirm_links.py` 使用兴隆站 topocentric 位置、日心观察者位置、Lambert 解和 `poliastro` 做后验确认。单夜 profile 的主要参数为：

- 初始日心距离假设：1.3, 1.7, …, 4.1 au；
- 最小初始地心距离：0.02 au；
- 最大允许初速度：35 km/s；
- 单夜 RMS 阈值：5″；
- 单夜最大残差阈值：8″；
- outlier clip：5″；
- 至少使用 80% 的观测。

它同时计算简单的天球线性运动 RMS、速度和方向。对只有 3 个点的单夜短弧，轨道解主要用于拒绝明显不一致的组合，不应作为精确物理轨道。

### 7.6 known 扣除、异常保护和产物

`summarize_single_night.py --require-mask15` 只接受当前 1.5″ known mask，或明确写明 known 已完成但 mask 为空的 status，不再回退到旧 1″产品。

如果扣除 known 后仍有超过 200 条 unknown link，单夜 wrapper 会以专用状态码跳过 review/export，避免把明显异常夜直接送入人工界面。已经出现过的高候选夜包括：

- `20251226=366`；
- `20260111=378`；
- `20260528=926`；
- `20260611=653`。

unknown 主要产物为：

```text
/processed1/<night>/L4/<night>_unknown_links.{json,fits}
/pipeline/xiaoyunao/heliolincrr/plots/<night>/*.gif
/pipeline/xiaoyunao/heliolincrr/review_packages/<night>/
```

### 7.7 trkSub 与人工复核

`assign_unknown_trksub.py` 使用 7 位 base-62 编号，字符表为 `0-9a-zA-Z`，全局历史文件为：

```text
/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl
```

分配过程加文件锁，并保存 night、linkage、tracklet、原图、objID、MJD、RA/Dec 和拟合摘要，使每个上报标识可以追溯到原始探测。

`package_unknown_review.py` 为每夜生成：

- 每条 link 的多帧 GIF；
- 人工 check CSV；
- 全 detection FITS；
- ADES 对照 FITS；
- manifest 和可选 tar 包。

网页最终写出 `<night>_submit.csv`，每条候选必须明确为 `is_real=0/1`。`watch_submit_reviews.py` 监测文件签名，检查缺标签、非法值和重复 key，只将 `is_real=1` 的 link 传给 `submit_reviewed_unknown.py`。全假夜和零候选夜记为 `no_observations`，不发空提交。

### 7.8 自动候选统计

当前服务器有 124 个 `unknown_links.json` 夜级文件，范围为 2025-11-16–2026-07-15：

- 非空夜 116 个；
- 零结果夜 8 个；
- 总候选 link 4,762 条；
- link 成员 detection 共 14,299 个；
- 候选 link 的观测点数中位数为 3；
- 候选轨道拟合 RMS 中位数 0.614″，90% 分位 1.058″。

网页人工 submit CSV 共 117 夜、4,762 行：

- 初次判为真：67 条；
- 判为假：4,695 条；
- 空白/非法标签：0；
- 初次通过比例：67/4,762 = 1.407%。

这说明当前自动候选集的主要成本仍是高假源率。它不等同于“算法 completeness 很低”，因为没有统计被算法漏掉的真目标；它只说明候选 precision 在人工复核前约为 1% 量级。

### 7.9 unknown 上报口径

历史自动补报区间 `20251116..20260617` 的 watcher 已完成：

- 122 个 review package 全部结束；
- `pending=0`、`failed=0`；
- 35 夜有真源并提交；
- 87 夜为 `no_observations`。

服务器当前保留：

- 67 个 submitted `trkSub`；
- 206 行 unknown ADES 观测；
- 35 个 `*_unknown_mpc_reply.txt` 和 submission ID。

这 67 条是“第一次人工 check 后进入 MPC submit 的单夜 link”，不是 67 个已确认的新天体。

### 7.10 后验高可信样本

提交后再次逐个检查 GIF，确认 2026-05-29/30 的 9 条 link 实为假源，并从本地科学分析 CSV 中删除了对应的 27 次 detection。最终建议用于论文图表的样本为：

- 58 条高可信单夜 link；
- 179 次 detection；
- 34 个观测夜；
- 55 条为 3 个点，2 条为 4 个点，1 条为 6 个点；
- 相对全部 4,762 候选的后验保留率为 1.218%。

对应本地证据文件：

```text
/Users/yunaoxiao/Desktop/submitted_unknown_20251116_20260628/
  submitted_unknown_true_detections_photometry.csv
  analysis_after_false_removal/analysis_summary.json
  analysis_after_false_removal/moon_analysis/moon_analysis_summary.json
```

58 条 link 的稳健统计为：

| 指标 | 结果 |
|---|---:|
| `Mag_Aper4` detection 中位数 | 17.73 mag |
| `Mag_Aper4` 10%–90% | 16.25–18.73 mag |
| `Mag_Aper4` 全范围 | 15.05–19.59 mag |
| 轨道拟合 RMS 中位数 | 0.124″ |
| 轨道拟合 RMS 10%–90% | 0.027″–0.379″ |
| 线性速度中位数 | 360.87″/day = 15.04″/h |
| 线性速度 10%–90% | 221.33–707.01″/day |
| detection 黄道纬度范围 | −14.37°–+14.91° |
| `|β|` 中位数 | 4.83° |

短弧拟合给出的 `a` 中位数约 2.20 au、`e` 中位数约 0.31、`i` 中位数约 22.94°，但样本中仍有负 `a`、`e > 1` 和高退化解；这些数字只可作为诊断，不宜作为论文主要物理分布。

### 7.11 晨昏和月光环境

对 58 条 link 取观测时间中值，并用兴隆站位置、Sun altitude = −18° 定义 astronomical dusk end/dawn start：

- 27/58 = 46.6% 距离最近的晨昏蒙影边界小于 1 h；
- 37/58 = 63.8% 小于 2 h；
- 最近的是 dusk end：38 条；
- 最近的是 dawn start：20 条；
- 最近晨昏边界间隔中位数为 1.21 h。

现有探索性月光指标为：

```text
moon_impact = illuminated_fraction × exp(−moon_separation_deg / 40)
```

其结果为：

- `moon_impact` 中位数 0.275；
- 与最近 twilight 间隔的 Spearman 相关 −0.09，基本不显著；
- 与 link 中值 `Mag_Aper4` 的 Spearman 相关 −0.20，只有弱趋势；
- 58 条中有 36 条在中值时刻月亮位于地平线下。

因此这个指标当前只适合探索性图表。若进入论文，应乘以月亮高度或使用天空背景模型，并对全部真/假候选一起分析，才能讨论月光对通过率的影响。

## 8. Follow-up 程序和当前状态

### 8.1 已实现功能

2026-06-24 已完成并部署：

- `survey/followup.py`：follow-up 状态机；
- `survey/apply_followup.py`：将 active source 插入现有观测计划；
- `heliolincrr/associate_followup_links.py`：把新夜 unknown link 关联回 follow-up 源；
- `watch_submit_reviews.py --enable-followup`：人工 check 后立即触发同日计划更新。

默认策略为：

- 只接收 origin night `≥ 20260624` 的新 reviewed true source，不回追旧样本；
- 对单夜 detections 做 RA wrap 安全的线性 RA/Dec 拟合；
- 将运动外推到计划曝光时间；
- 选择能够覆盖预测位置且中心最近的既有 footprint；
- 目标名为 `MP_FU_<trkSub>_<field_id>`，保证继续进入 MP 的 known/unknown 数据链；
- 每源每夜插入 5 帧；
- 实际 L2 满 5 帧才记为该夜成功；
- 需要 2 个成功夜才完成；
- 不满 5 帧的夜记为 failed，后续继续排；
- 发现后超过 10 天仍未完成则 abandoned；
- 新夜 link 与 active source 的预测位置中值距离 ≤ 30″时关联，并用新 detections 更新运动模型。

### 8.2 当前实际结果

截至 2026-08-03：

```json
{
  "start_night": "20260624",
  "sources": {},
  "ingested_nights": {}
}
```

即：

- active source：0；
- completed source：0；
- abandoned source：0；
- 实际 follow-up 夜：0；
- `MP_FU_...` 实际曝光：0。

原因不是程序未部署，而是当前 67 条提交 link 的最新 origin night 是 2026-06-01，早于 follow-up 起点；后来的 2026-06-20 候选全部人工判假，2026-07-15 为零候选，因此没有源进入新的状态机。

论文中可以把 follow-up 写入系统设计和 future work，但当前不能把它列为观测成果。若论文希望形成真正的 discovery/confirmation 故事，最关键的新工作就是用一个真实新源跑通“两夜 × 5 帧 + 跨夜关联 + MPC 状态”的完整闭环。

## 9. 已开发程序清单

### 9.1 `survey`

| 文件 | 功能 |
|---|---|
| `config.py` | 站点、调度器和服务器路径配置 |
| `astro_utils.py` | 夜晚边界、高度角、月距和稳定可见性 |
| `history.py` | footprint 历史表读写 |
| `history_update.py` | 从 L2 MP 文件增量更新历史 |
| `scheduler.py` | 加权选场、cluster、三访和 near-Sun 排程 |
| `io_utils.py` | plan JSON/TXT 输出 |
| `run_daily.py` / `run_daily.sh` | 每日正式入口和发布 |
| `visualize_nightly.py` | 每 cycle、全夜和 GIF 可视化 |
| `make_yearly_coverage_animation.py` | 年度覆盖演化动画 |
| `followup.py` | source 状态、预测、选场、核对和关联 |
| `apply_followup.py` | 将 follow-up 插入 nightly plan |

### 9.2 `known_asteroid`

| 文件 | 功能 |
|---|---|
| `match_single_night.py` | 星历查询、CCD 内筛选、1″/1.5″匹配 |
| `build_file_manifest.py` | 夜次文件 manifest 和时间归属 |
| `submit_pipeline_slurm.sh` | 单夜/日期段 Slurm 提交 |
| `slurm_match_one_file.sh` | per-file array task |
| `slurm_merge_submit.sh` | finalize、ADES、MPC submit |
| `merge_night_parts.py` | 夜级归并与状态文件 |
| `export_ades.py` | ADES、质量筛选、历史去重、MPC 接口 |
| `run_daily.sh` | unattended daily trigger |
| `plot_known_asteroids.py` | 覆盖、全天、统计、field yield 和 GIF |
| `run_visual_daily.sh` | 每日可视化入口 |
| `update_all_matched_history.py` | known 历史累积 |
| `audit_duplicate_causes.py` | 重复提交审计 |
| `run_known_rematch_then_unknown_remask.sh` | known→unknown 的 Slurm dependency 重算 |

### 9.3 `heliolincrr`

| 文件 | 功能 |
|---|---|
| `mask_gaia.py` | Gaia proper-motion mask 和 L2 预筛 |
| `make_tracklet_linreproj.py` | 两点移动 tracklet 构造 |
| `merge_tracklets_night.py` | 全夜 tracklet 合并 |
| `run_linear_links_from_tracklets.py` | 共享端点的单夜线性 linkage |
| `orbit_confirm_links.py` | Lambert/日心轨道一致性确认 |
| `summarize_single_night.py` | 统一统计、known 分类和 unknown 输出 |
| `assign_unknown_trksub.py` | 全局 7 位 base-62 `trkSub` |
| `plot_unknown_links.py` | 多帧 GIF 和 link summary |
| `package_unknown_review.py` | 人工复核包和 manifest |
| `export_unknown_ades.py` | unknown ADES PSV |
| `submit_reviewed_unknown.py` | submit CSV 过滤、validate/submit、stats |
| `watch_submit_reviews.py` | 自动监测网页 submit 和幂等状态 |
| `run_review_submit_backlog_watch.sh` | 历史 review 包补报 |
| `remask_unknown_with_known.py` | known 更新后的 mask-only 重建 |
| `run_daily_unknown.sh` | known ready 后的 unknown daily 入口 |
| `run_daily_pipeline.sh` | survey→known→unknown 总入口和恢复 |
| `associate_followup_links.py` | follow-up 跨夜关联 |
| `aggregate_nightly_stats.py` | 多夜统计和诊断图 |
| `run_pipeline_15.sh` 等 | 15 夜 link/orbit 实验入口；不属于当前主结果口径 |

## 10. 工程可靠性方面完成的工作

除核心算法外，项目还解决了多项实际运行问题：

- Slurm 任务 ID 解析、瞬时提交失败重试、单文件跳过和 finalize 依赖；
- known 夜级强制归并，避免旧 FITS 因“文件已存在”而保留；
- RA=0° wrap 星历查询修复和全区间 known rematch；
- known 1″上报表与 1.5″ unknown mask 表分离；
- known 更新后 unknown mask-only 重打包，复用已有 GIF，避免昂贵重画；
- 0 tracklet、0 link、0 unknown 夜均能写 schema 完整空产物；
- 异常 dense group 和 unknown > 200 夜的保护；
- review submit CSV 的缺标签、非法值、重复 key 检查和自动等待重生成；
- watcher state 按 manifest/submit 文件签名去重，避免重复上报；
- daily 总入口、`@reboot` 恢复和幂等 skip；
- 断电后恢复 survey plan、known/unknown 数据处理和 backlog watcher；
- known 与 unknown 的 ADES/GIF/manifest/response 路径标准化。

这些不是附属细节，而是“系统论文”定位的重要贡献：项目已经从离线脚本集合演变为能在真实数据、Slurm、人工复核和断电条件下长期运行的 pipeline。

## 11. 当前统计的主要限制

### 11.1 unknown 不是 confirmed discoveries

58 条是高可信单夜移动源 link，不一定对应 58 个独立天体。需要：

1. 按时间和轨道做跨夜关联；
2. 用最新 MPC/JPL catalog 回查这些位置和时间；
3. 查询每个 submission ID 的 MPC/WAMO 最终 ingest 状态；
4. 明确哪些获得 designation、哪些被识别为已知对象、哪些被拒绝；
5. 对 9 条后验假源确认是否需要 MPC 更正或在论文中单列。

### 11.2 unknown 假源率高

初次人工通过率只有 1.407%，后验高可信保留率 1.218%。论文应把它作为 pipeline 现状和改进空间，而不是隐去：

- Gaia/known mask 仍会留下坏像素、衍射/边缘、静态源残差和错误组合；
- 3 点短弧可以被非常简单的线性模型拟合；
- `fit_ok/is_good` 是动力学一致性条件，不是真源分类器；
- 需要图像形态特征、背景质量、PSF/shape、跨夜信息或监督分类器进一步提纯。

### 11.3 没有 completeness

目前有 known tracklet benchmark，但没有对整条 unknown pipeline 做可控注入：

- 不同星等、速度、方向、背景和月相下的检测效率未知；
- 不能从高假源率推断漏检率；
- 不能给出严格 limiting magnitude 或 sky-plane rate completeness。

最有价值的补充是把已知小行星当作 truth set，逐级统计：L2 detection → Gaia mask 后保留 → tracklet → linear link → orbit confirm → known subtraction前的恢复率，并做少量人工注入验证。

### 11.4 观测和处理口径不同

- raw 标准科学帧：41,074；
- L2 MP catalog：40,399；
- known matched 夜：130；
- unknown catalog 夜：124；
- 有最终 submit CSV 的 unknown 夜：117。

论文流程图和表格必须把每一层的 drop/skip 原因写清楚，不能用一个“观测夜数”贯穿所有模块。

### 11.5 站点高度参数需要统一

survey 和 known 默认使用 960 m；`orbit_confirm_links.py` 的 MPC 327 常量为 868.221 m。对当前单夜筛选影响很小，但正式 methods 和可复现版本应统一为一个经过确认的站点位置。

## 12. 建议的论文结构

### 1. Introduction

- 小口径广域望远镜在黄道小天体监测、已知目标补充天测和快速候选筛查中的作用；
- 本项目目标：建立可长期运行的调度—观测—known—unknown—人工复核—随访链路；
- 文章贡献和数据区间。

### 2. Telescope, Survey Footprints, and Observations

- 60/90 cm Schmidt、CCD、pixel scale、曝光；
- 1,855 场 footprint 和黄道纬度覆盖；
- 2025-11-15–2026-07-15 的 134 个科学夜、41,074 帧、10,387 deg²；
- 数据预处理由合作同事完成，本项目从 L1/L2 接手。

### 3. Dynamic Survey Scheduling

- 夜晚和可见性约束；
- 月距、历史覆盖、最近访问、稀缺可见期权重；
- cluster 三访和 near-Sun 模式；
- 计划—实际 95.17% field-count 一致率。

### 4. Known-Asteroid Identification and Reporting

- astorb/aleph topocentric ephemeris；
- WCS footprint 和 RA wrap；
- 1″ official / 1.5″ mask；
- Slurm、ADES 和历史去重；
- 58,482 对象、534,780 匹配、残差分布和有 reply 的提交统计。

### 5. Single-Night Unknown Moving-Object Pipeline

- Gaia mask；
- tracklet 参数；
- shared-endpoint linear link；
- orbit consistency；
- known subtraction；
- trkSub、GIF、网页 check 和 MPC submit。

### 6. Results of the Unknown Search

- 4,762 → 67 → 58 的筛选漏斗；
- 58 条高可信 link 的星等、速度、黄道纬度和晨昏分布；
- 假源类型和后验 9 条剔除；
- 明确它们仍是 single-night candidates。

### 7. Operational Automation and Follow-up Framework

- daily pipeline、Slurm、watcher、reboot recovery；
- follow-up 状态机设计；
- 当前尚无实际 follow-up，作为 future validation。

### 8. Discussion

- 小望远镜的 known astrometry 产能；
- 动态调度的运行价值；
- unknown precision 的瓶颈；
- twilight/near-Sun 候选的意义；
- 下一步 completeness、MPC 状态核对和真实 follow-up。

### 9. Conclusions

- 强调端到端系统和实地运行；
- known 大样本成果；
- 58 条高可信单夜候选；
- 不越过当前证据宣称 confirmed discoveries。

## 13. 建议优先制作的图和表

### 主文图

1. 端到端流程图：scheduler → observation → L1/L2 → known/unknown → review → MPC/follow-up；
2. footprint 与实际曝光频次 Aitoff 图，标黄道和有效 10,387 deg²；
3. 逐夜计划/实际曝光和累计覆盖；
4. known 全天分布或 HEALPix 密度图；
5. known 位置残差、星等和每目标访问次数；
6. unknown 筛选漏斗：4,762 → 67 → 58；
7. 58 条 link 的 RA/Dec + 黄道、星等和速度分布；
8. unknown 距最近 astronomical twilight 的分布；
9. 代表性 known GIF 与高可信 unknown GIF。

### 主文表

1. 望远镜、相机和 scheduler 参数；
2. raw/L2/known/unknown 各层数据规模；
3. known 匹配和 ADES 统计；
4. unknown 候选、人工通过、后验清洗统计；
5. 58 条 link 的汇总或完整在线表。

### 附录或补充材料

- 全部 58 条 link 的 `trkSub`、night、linkage_id、MJD、RA/Dec、星等、速度、拟合残差；
- 9 条后验假源清单；
- 所有夜次的处理状态和质量 mask；
- known/unknown ADES submission ID 对照；
- scheduler 伪代码和全部默认参数；
- 软件版本、commit 和服务器环境。

## 14. 团队现在需要拍板的问题

1. 论文主标题更偏“survey system”还是“known-asteroid astrometry”；推荐前者、known 作为最大科学结果。
2. 观测起点是否统一为 2025-11-15，还是只写当前 scheduler 部署后的 2026-03-29；推荐主体使用完整观测期，同时单独评估新 scheduler 时段。
3. known 的主数字使用 534,780 匹配还是 348,568 有 reply 的 ADES；推荐两者都写，分别叫 matched 和 submitted-with-reply。
4. 58 条候选是否逐条完成最新 MPC/JPL crossmatch 和 submission status 核查；这是发表前必须做的。
5. 9 条已提交但后验判假的源如何处理和表述。
6. 是否等待第一条真实自动 follow-up 闭环再投稿；若想强化 discovery 故事，建议等待；若定位系统论文，可以先写并把 follow-up 放 future work。
7. 是否增加 injection/recovery 或 known-object truth-set completeness；建议至少完成后者。
8. 是否固定一个数据冻结日并保存不可变的 night mask、统计 CSV、图表和软件 tag。

## 15. 统计证据和可复现路径

本总结交叉使用以下来源：

- 仓库 git history、`WORKLOG.md`、`PLAN.md` 和三个模块 README；
- 服务器 `/raw1/YYYYMMDD/OBJ` 原始文件；
- 服务器 `/processed1/YYYYMMDD/L2` 和 `L4`；
- `/pipeline/xiaoyunao/survey/runtime/plans`；
- `/pipeline/xiaoyunao/survey/runtime/followup/followup_state.json`；
- `/pipeline/xiaoyunao/data/heliolincrr/review_submit_backlog_20251116_20260617.json`；
- `/pipeline/xiaoyunao/heliolincrr/review_packages`；
- 本地后验清洗后的 58-link photometry CSV 和 analysis summary；
- 以前项目对话中的 2026-07-17 HEALPix 覆盖面积统计。

正式投稿前建议将这些只读统计固化为一个版本化脚本和一个 `paper_data_snapshot.json`，并给每一张图记录输入文件 hash。当前数字是 2026-08-03 快照；daily scheduler 仍在运行，未来新增数据后数字会变化。
