# 2026-11-05 年度设计收尾

常规展示采用 B，明确优先 Thymbraeus 和 Paeon 时采用 C。基准精度的 200 条配对天气实现中，B 的合格光变数 P10/P50/P90 为 41/47/52，C 为 38/44/48。晨昏仍以五访为覆盖基线，晨间十访后对齐作为短弧增强选项。百颗目标仍是待验证的科学目标。

## 图册与复现

默认 B 图册从 `B_plot_data` 的 365 夜、37,264 条名义曝光生成。排程源提交为 `4605b372fde47692ddcd66ed809905b07caf05bb`。发布版本由 `git log -1 --format=%H -- targeted_observing/docs/continuation_20261105/atlas_B/figure_manifest.json` 查询。旧 A 的 365 夜、37,078 条曝光和图册保存在原历史目录。B 图册按北京时间夜次索引，无 Trojan 的日期显示原排程的跳过理由。逐夜图用点显示曝光中点，用条带区分曝光、读出和转场。午夜后的刻度属于下一民用日。

绘图基准取自本地 Pal 5 New 的 `scripts/step5_stream_model.py::plot_density_track_summary`、文件顶部的 rcParams、同文件的配色常量及保存函数。该文件最近提交为 `c57e689b4fa456355f484c089e558a6f67c49f80`，实读文件 SHA256 为 `7a94072df9a9445166b838fb18a729dc8568ffa3be94ca73d7c83298aebd4e82`。复用 13/17/18 pt 字号体系、140 dpi 预览、300 dpi 输出、浅网格、蓝橙绿和 viridis、空心标记及误差棒设置。字体与数学字体显式固定为原入口沿用的 Matplotlib 默认 DejaVu Sans。SHARP 按全天投影、三联月图和两联逐夜图分别设置画幅，并额外输出矢量 PDF。

[默认 B 图册](atlas_B/index.html) · [旧 A 图册](atlas/index.html) · [B 排程源数据](B_plot_data/manifest.json)

```bash
python3 -m targeted_observing.review_tools.plot_exported_plan targeted_observing/docs/continuation_20261105/B_plot_data --output targeted_observing/products/annual/B_publication_figures --all-nights --workers 5 --iers targeted_observing/cache/iers_20261104.all
python3 -m targeted_observing.review_tools.publish_atlas targeted_observing/products/annual/B_publication_figures targeted_observing/docs/continuation_20261105/B_plot_data targeted_observing/docs/continuation_20261105/atlas_B
python3 -m targeted_observing.final_sampling_stats
python3 -m targeted_observing.final_twilight_bins
python3 -c 'from targeted_observing.continuation_figures import plot_comparisons; plot_comparisons()'
python3 -m targeted_observing.continuation_twilight_report --plots-only
```

本地 `products/annual/B_publication_figures` 有 378 组单张 PDF 和 300 dpi PNG。GitHub 交付两卷矢量 PDF 和供日期索引用的轻量预览。打包需要 Pillow 和 pypdf。五张统计图在本目录各有一对 PDF/PNG。

## 驻场、轮巡与跨夜

| 样本 | 驻场两张间隔 P10/P50/P90 (s) | 同场轮巡首张到下轮首张 P10/P50/P90 (s) |
|---|---:|---:|
| B 名义年 | 135 / 135 / 135 | 660 / 1325 / 1366 |
| B 基准天气实现 0 | 135.06 / 135.34 / 136.12 | 1309 / 1340 / 1366 |
| C 基准天气实现 0 | 135.06 / 135.34 / 136.12 | 665 / 1337 / 1367 |

同场轮巡跨度包含中间场的驻场与转场。表格另列上轮末张到下一轮首张。场夜表保留逐场曝光数及中点跨度。样本范围由保存的完整时序决定，B/C 天气 cadence 只统计实现 0。

实现 0 的 B/C 分别有 41/42 个合格对象。在保存的首次合格前缀内，31/33 个无法在单个固定 35 夜 campaign 独立通过。将前缀在有效夜间隔 ≥35 日处切开，10/9 个在切开后没有任何一段独立达标。该诊断使用排程的固定 campaign 和有效夜间隔衡量跨段依赖。200 条实现另有原始证书统计。合格对象×实现合并后的有效夜中位数均为 13 夜，证书基线中位数 B 为 38 日、C 为 42 日。

[分位数](cadence_quantiles.csv) · [分布图](cadence_distribution.pdf) · [场夜统计](cadence_field_nights.csv) · [逐对象时序](BC_object_timing_realization0.csv) · [逐 campaign 时序](BC_campaign_timing_realization0.csv) · [200 次证书统计](BC_ensemble_certificates.csv)

## 特殊科学任务的缺口

C 中 Thymbraeus 和 Paeon 的任务就绪率为 99.5% 和 87.5%，B 中为 56% 和 35.5%。Thymbraeus 在 B 最常缺有效夜，C 同时补足点数、夜数和相位。Paeon 在 B 有 56% 的年份有效夜不足、30% 点数不足，C 分别降至 7% 和 3.5%。Thessandrus 在 C 的基线中位数约 104 日，四分之三的年份通过基线条件，但所有年份的最弱候选周期相位覆盖仍低于 80%。Cloanthus 和 617 在 C 的相位失败率分别为 100% 和 97%。

C 实现 0 中，5012、5284、7543 分别有 126、50、91 个通过可用标记和零点筛选的点，应用各自任务单点误差限后均剩 0 点。该实现明确显示精度先卡住了采样，200 次任务表也显示这些目标的中位任务好点偏少。失效条件可以重叠，物理解与周期解状态保留未知。

| 角色 | B 期望合格数 | C 期望合格数 |
|---|---:|---:|
| 可靠控制 | 32.165 | 29.550 |
| 已有模型 | 8.140 | 7.530 |
| 弱周期／未知 | 4.445 | 3.950 |
| 明确别名复核 | 1.000 | 1.000 |
| 文献未核查 | 1.040 | 1.000 |

各列由逐对象质量频率求和。文献未核查和弱周期对象仍需研究核验。

[八目标门槛与重叠失败率](BC_focus_science_gates.csv) · [实现 0 精度与各候选周期相位](BC_focus_precision_realization0.csv) · [角色产额](BC_role_yield.csv)

## 晨昏距日角

晨间 45°–55° 箱内，五访有 50 / 12,480 个曝光中心、10 / 2,496 个合格场窗口，球面像素并集约 95.1 deg²。十访后对齐在这一箱三个值均为 0，55°–65° 箱的并集也从约 2,276.2 降至 889.5 deg²。两方案晚窗相同。所有值来自名义年。

每个曝光中心在自身中点分箱。场窗口以最接近窗口中点的实际曝光中心分箱，资格为至少四次曝光且中点跨度至少 720 秒。足迹使用 NSIDE=256 球面像素中心，在每次曝光中点重新计算站心太阳距角。每侧、每策略、每箱取像素并集。像素可以随日期进入不同距角箱，各箱面积不能相加当作独立年覆盖。外侧箱保留足迹越过 90° 的面积。恢复图采用基准暗空深度 20、注入 r=19，初始短弧以注入数为分母，后续恢复以已有合格短弧且未右删失的注入为分母。

[六个距角箱完整统计](twilight_elongation_bins.csv) · [晨间中心几何](morning_geometry.pdf) · [条件恢复](morning_recovery.pdf)

现场阶段验证实际 r30/r120 误差和跨夜零点、污染及饱和、控制器输入、读出与转场。快速目标需另行验证生产管线当前 3–63 arcsec/hour 速度范围与注入范围之间的差距。真实周期与物理模型推断留待观测数据。年度设计在此阶段冻结。
