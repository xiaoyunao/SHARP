# SHARP TJN 固定天区实施要求

更新日期：2026-11-11，北京时间。工作对象：`xiaoyunao/SHARP` 的 `codex/twilight-trojan` 分支。本文件可直接交给 Codex 执行。本轮目标是建立可长期复用模板的 L4/L5 固定场库；这一轮不重构整个 pipeline。

## 1. 本轮目标与代码基线

用户要求：新 pipeline 以后由晨昏、特洛伊、YSO、超新星、后随、主程序六部分组成，并接 known/unknown 等处理模块。当前只做 TJN 固定天区，其余需求进入项目备忘。

已核对的版本：

- `9734af0a04253bf3e7f2caccc68c0a87f48bfe06`：晨昏各完整一小时、离线预设可变曝光。
- 后续 `a99f06fc`：EV→M31→TJN→YSO→AM 整夜测试拼接。
- 后续 `1ab29e94`：固定目标名、控制滤镜改为小写。
- 本次核对的分支 HEAD 为 `7747cc1650fdeb255e677e496f50a376a99cb57b`：用户修改列表后的时角审查。当前控制 token 为小写 `g/r/i/z`；不要被较早文档中的大写 `R` 误导。

开始实现前读取当前 HEAD、`AGENTS.md`、`targeted_observing/README.md`、`targeted_observing/PLAN.md`、`WORKLOG.md` 和 `docs/GOAL.md`，检查后续新增提交。不要回退或覆盖用户今日编辑的执行列表。

实现位于 `targeted_observing/`，延续现有生产隔离。保持 `main`、`survey/`、`known_asteroid/`、`heliolincrr/` 及服务器现状。当前任务授权的是固定场数据、生成器、审计、必要的 TJN 场查询接入；不包含现场部署、服务器操作、整套 scheduler 重构或新晨昏全年重跑。

## 2. 必须改变的概念

旧年度流程是 35 日 science campaign、7 日 pointing version。旧 `annual_trojan.version_fields` 和 `annual_free_fields.supplement` 在版本中点跟随目标重新定心，场 ID 又包含日期或坐标 hash。这只能在短期固定，不符合长期暂现源模板需求。

新设计固定的是静止天球上的场注册表：

1. 同一个 `field_id` 永远对应同一 ICRS 中心、同一规定的相机方向及同一明确几何定义。
2. 每晚只从已存在的场中选活动子集，变化的是访问时间、优先级及目标成员。
3. 小行星跨过场边界后转到另一个固定场；光变仍按稳定 `object_id` 跨场积累。
4. 35 日科学任务管理可保留；7 日可用作重评活动场的周期，但不得再改变指向中心。
5. 每个静态场可以在同一幅图像中服务 TJN 与暂现源分析；并不要求整个静态场库每天全部拍摄。

## 3. 原 MP 的实测表结构核查

直接核对 `targeted_observing/data/mp_footprints.json`，共 1,855 场。其来源记录为服务器 FITS 快照 `survey_fov_footprints_with_visibility.fits`，源 SHA256 为 `13f394877ce26e50e5cee5053a4ad845f412be4520421e46931b3e157d9e0225`，提取时间为 2026-11-09 17:30:02 +08:00。

| 参数 | 从 MP 表算得 |
|---|---:|
| RA 宽度 × cos(中心 Dec) | 2.96960000° |
| Dec 高度 | 2.9721777778° |
| 典型最近邻球面中心距，中位数 | 2.67262012° |
| 采用原宽度的 90% 作东西节距 | 2.67264000° |
| 采用原高度的 90% 作 Dec 行距 | 2.67496000° |

原 MP 中心沿近似黄道纬度带排列，1,855 个中心 Dec 均不同。新 TJN 应继承其视场尺度、约 10% 重叠及固定编号思想，改成赤道坐标的恒 Dec 行；不要直接把黄道网格旋转后当成恒 Dec 行。

原 MP 角点为 RA–Dec 矩形近似，不是实际 TPV WCS。计算角点宽度必须先处理 RA 跨 0°，否则会出现约 −357° 的假宽度。

### 视场口径不能混用

9216×9232、1.15″/pix 的名义 TAN 平面尺寸约为 2.9440°×2.9491111°。MP 表的东西尺寸对应 9216×1.16″，而南北尺寸在 1.16″/pix 下对应 9224 像素，较 9232 少 8 像素；原因尚不明确。不能擅自认定它是裁边或 overscan。

请从现有本地图像/保存的 WCS header 核实有效科学像素范围、像元尺度、TPV 畸变、相机方向和模板实际裁边；保留原 MP 文件。没有可用本地实测 WCS 时先完成名义几何，明确留下现场几何验证项，不宣称实测覆盖已通过。

## 4. 推荐的母网格

采用固定锚点 `RA=0°`、`Dec=0°`，并继承 MP 的中心节距：

```text
sx = 2.67264 deg
sdec = 2.67496 deg
dec_j = j * sdec
N_j = ceil(360 * cos(dec_j) / sx)
ra_step_j = 360 / N_j
ra_jk = k * ra_step_j,  k = 0 ... N_j-1
```

三角函数计算时转换为弧度。这个整数整圈构造比对有限区域逐次累加 RA 更容易保证跨 0° 的闭合；每行实际东西节距略小于或等于目标节距。实际只选择需要的 L4/L5 区域，不是要观测每一整圈。

小角度关系 `ΔRA ≈ sx/cos(Dec)` 可作解释和检查，不应对所有 Dec 都固定同一个 RA 角间隔。所有中心在恒 Dec 行上；不要求不同行的 RA 列严格对齐。首次冻结前允许对整体网格相位作一次有明确目标函数的比较，若改相位必须重新生成并说明收益；冻结后不能每季重新优化锚点。

### 正南正北的含义

- 每个场中心的相机北向固定，采用明确记录的 PA 约定；本草案为局部北向 `PA=0°`。
- 相机像素 x 方向的镜像/宇称另行记录，不能把显示图的左右方向当作硬件 PA。
- 恒 Dec 的是场中心行。真实 TAN/TPV 视场边缘在 RA–Dec 图中可能略弯曲或倾斜，不应强行画成经纬矩形后拿该矩形做覆盖检验。
- 若设备不能设置 PA，先用实测 WCS 识别其固定方向和变化范围。不能只写 `pa_deg=0` 就认为硬件北向已经实现。

### 这里的 RA 投影不得改动 slew 模型

网格的物理宽度需要 cos(Dec)，但赤经机械轴的有效位移采用 wrap 后的 RA 坐标差。当前实测 slew 模型的 RA 项不能因为新网格而乘 cos(Dec)。本轮只报告固定网格造成的转场差异，不改 overhead 系数。

## 5. 覆盖范围必须从全年候选轨迹确定

保持目前比较年度：2026-11-01 至 2027-10-31，停止日期为排他 `2027-11-01`。以后年份沿同一母网格追加场，不移动已有场。

冻结原母表有 16,445 对象。原每周预筛采用 `H < 14`，存在一次 `V_proxy ≤ 19.7` 且 `Dec ≥ −29°` 的采样，得到 1,063 个年度候选：L4 610、L5 453。现有云标签来自相对木星的瞬时日心黄经代理，应保留其来源；不要按视赤经左右给对象重新分云。缺标签、矛盾标签及特殊对象应单列。

本轮已把旧 B 场的 core、secondary、anchor 对象 ID 与原母表关联，结果无缺失标签或混云场。

| 旧 B 数据 | L4 | L5 |
|---|---:|---:|
| 旧唯一场 ID | 199 | 184 |
| 场夜记录 | 619 | 511 |
| 旧场中心 RA 范围 | 179.951860–244.775920° | 35.806600–128.520330° |
| 旧场中心 Dec 范围 | −14.342071–+34.236330° | −9.170280–+54.260960° |

这些是已选场的中心包络，不是两云的完整边界。不得把两云混合成一个 RA 36–245°的大矩形，也不得把中心包络矩形全部视作连续高密度区域。

### 生成最终有效区域的顺序

1. 读取冻结的 1,063 对象候选表和全部年度 Horizons 6 小时缓存，核验索引完整性、站点、时间系及输入 hash；缓存缺失不能用 hash 冒充轨迹。
2. 以全年轨迹建立候选占用区域。保留逐日期可见性，按当前 TJN 的暗夜/高度/亮度条件筛选“曾有可能观测”的片段；不以某一个夜晚、某一次月亮遮挡或旧 B 的胜出结果永久裁掉场。
3. 对快过场或靠边片段加密检查；6 小时采样的间隙不能简单当作没有目标。沿用已有插值与星历 QA，并记录星历误差/采样余量。单位向量或 unwrap 后的 RA 插值必须正确处理跨零。
4. 用实际足迹/安全区域与轨迹片段相交选择母网格单元。将旧 B 已选区域、今晚场及可行特殊科学目标作为明确对照；在边界按真实指向/轨迹余量保留必要缓冲。
5. 报告逐对象、逐日及两云分开的覆盖率、连续可用日数、跨场次数、单场同时目标数、边缘损失和未覆盖理由。既要有全候选样本，也要单列可测核心样本和特殊对象。
6. 给出核心活动场与候补/边界场的推荐标签，二者都来自同一固定几何库。标签可更新，中心不改。

静态场库的天区边界、对象的动态归属、每晚的排程选择应分开保存。不要在静态注册表写入会过期的“固定目标成员”。

## 6. 本包提供的 522 场几何草案

`tjn_fixed_grid_seed_v0.csv` 及两云分表是可复现的讨论起点，状态均为 `PROPOSAL_NOT_FROZEN`：

- L4：241 场，草案编号 `TJL4_001` 至 `TJL4_241`。
- L5：281 场，草案编号 `TJL5_001` 至 `TJL5_281`。
- 区域来自旧 B 的 383 个名义 3°×3°球面局部足迹的包络，与固定网格的分配单元相交。这个保守规则可能纳入额外边角单元；没有优化模板数。
- 同一个旧场会与多个新场对应；`legacy_to_seed_cells.csv` 保存关系。
- 独立以每个旧足迹 31×31 采样点检查，合计 368,063 点，在保守 1.15″/pix 名义 TAN 模型下无未覆盖点。这只说明旧足迹的采样几何覆盖，不是全 1,063 对象轨迹覆盖或真实 CCD 完整性。

**522 不是最终场数，也不是每夜场数。** 使用完整星历和真实 WCS 后可增减草案；正式第一次冻结编号发生在本轮最终验收后。冻结后只追加编号，不重排旧编号。草案里的 `legacy_associated_object_ids_json` 只追溯旧场成员，不保证这些对象真实落入该新场。

### 今晚 TJ_test 的对照

今晚旧中心为 `(79.002262°, +23.719287°)`，旧 B 场 ID 为 `AF_148ecdaae5`，关联对象为 17365 Thymbraeus、4722 Agelaos，均属旧标签 L5。

此中心落入草案 `TJL5_169` 的分配单元，其新中心为 `(79.0243902439°, +24.0746400000°)`。只是一个接近旧场的几何候选，并未验证两个移动对象在新场内的留场。

Codex 必须用今晚逐帧及后续 35 日的现有星历检查这两个目标、边缘余量和邻接切场。不要直接用新中心覆盖今晚已确认的执行文件。今晚已经积累的旧场图像可保留为独立历史几何/模板候选，不要硬改其身份。

## 7. 真正的足迹与安全区域验收

必须区分名义足迹、有效科学像素遮罩、模板公共区域、移动目标质量边缘余量。原 `in_field(width=3, margin=.15)` 中的 0.15°是单侧目标边缘条件，不能默认为这次网格一定可用的统一裁边。

本包 `audit/geometry_audit.json` 使用独立 TAN 公式，逐行密集检查网格分配单元边界，每条边 2,049 个点。在 Dec 中心 −26°…57°的实际 31 行中，两种未裁边名义模型均通过；最差行中心为 +56.17416°，1.15″模型最小余量约 0.08946°，即 5.37′。

若四条边各向内裁掉 0.15°，这 31 行在两种名义模型中都不能保证完整覆盖各自单元。这不意味着原观测无效；它说明旧目标质量留边与新网格完整覆盖是不同条件。

因此：

- 用实际 TPV WCS 的完整 world-to-pixel 变换检查对象是否落在有效像素/模板公共区域。
- 加密边缘像素，检查相邻场、相邻 Dec 行和高 Dec 角落的覆盖。`calc_footprint` 的四角、圆形半径或 RA–Dec 包围盒不能单独作最终验收。
- 裁边、指向重复误差、PA 变化和畸变共同消耗几何余量。参数必须有实测或明确工程假设依据。
- 首次冻结前，若安全足迹之间的间隙影响所需轨迹，应缩小全局节距并整体重新生成、重新验收。首次冻结后保留全部已有中心、编号与几何身份；需要补洞时另追加具有明确新几何键/版本的固定补场。若同一母网格相邻场裁边后本身存在间隙，单纯多选该母网格单元不能消除间隙。
- 不强求坏像素、亮星遮罩等不存在；应把实际不可用区域及损失明确统计。不能把无法检测的像素计入有效覆盖。

## 8. 永久编号与模板契约

第一次正式发行：按云分组、Dec 从南到北、同一行按固定 RA 展开顺序编号，分别从 `TJL4_001` 和 `TJL5_001` 起。三位数是最小宽度，超过 999 时继续 1000，不截断。

冻结后的规则：

- 相同 `field_id` 不移动中心、不改变约定 PA；新增区域只追加号码，停用号码不复用。
- 保留独立 `grid_cell_id`（母网格版本、行、列）和 registry release/version；不要按当前排序重建 field_id。
- 真正更改指向几何应创建新的字段身份或显式 pointing variant，旧数据不能被覆盖为新坐标。
- 将 `object_cloud` 与 `planning_cloud` 分离。固定天区以后可能被另一云经过；云前缀是最初规划分区，不是该场内所有对象永恒的物理分类。
- 相同空间单元跨项目复用时以几何主键去重，避免同一天区因为规划标签不同而建立互不关联的模板。

建议静态注册表至少含：`field_id, grid_cell_id, ra_deg, dec_deg, frame, pa_deg, orientation/parity convention, row_index, column_index, geometry_version, nominal_boundary, usable_mask_ref, created_release, status`。不要因格式化控制列表而降低主表坐标精度。

动态表至少含：`object_id, time/time_interval, field_id, predicted_ra/dec, magnitude_band/value, in_nominal, in_usable, edge_distance, target_role`。重叠场中同对象的多次真实测量可保留；覆盖/产额按对象去重。

模板至少按 `instrument + geometry_version + field_id + filter + template_version` 区分，并保留历次模板来源曝光。控制列表当前为小写 `r`；模板需要真实物理滤镜标识，不将 `g` 与 `r` 混成同一个模板。相机方向、光学畸变或共同区域改变时，明确新模板几何关系。

## 9. 必要输入与本包的数据范围

| 输入 | 位置或本包文件 | 用途 |
|---|---|---|
| 原 MP 1,855 场 | `inputs/mp_footprints.json` | 原尺度、重叠与编号比较 |
| B 原始场夜表 | `inputs/B_fields_original.csv` | 原排程指向来源 |
| B 场夜分云表 | `inputs/tjn_B_field_nights_with_cloud.csv` | 1,130 场夜及云来源 |
| B 唯一旧场 | `inputs/tjn_B_unique_fields_with_cloud.csv` | 383 场几何对照 |
| 年度候选 ID | `inputs/tjn_annual_candidates_1063.csv` | 固定原筛选；不含轨迹 |
| 全母表筛选 | `inputs/tjn_parent_screening_16445.csv` | 保留拒绝对象与理由 |
| 16 特殊目标 | `inputs/tjn_special_16_cloud_labels.csv` 和原科学说明 JSON | 独立科学标签/几何验证 |
| 本地完整年星历 | `targeted_observing/cache/annual/horizons_*_6h.json` | 最终轨迹范围的必要输入，未提交 GitHub |
| 本地完整索引 | `targeted_observing/products/annual/ephemeris_index.json` | 核验1,063对象星历完整性 |
| 本地候选表 | `products/annual/candidates.json`、`candidate_audit.csv` | 与本包候选快照对照 |
| 年星历 hash 清单 | `docs/continuation_20261105/annual_input_manifest.json` | 查验输入版本，不是星历本身 |
| 本地实测 WCS/有效图像区 | 由本地已保存图像/header提供 | 正式 footprint/PA/裁边验收 |

特殊名单的 38615、228155 不在原1,063母样本内，原拒绝理由是每周亮度或赤纬预筛。若要给两者做专门场，应单列新增候选及实际可行性，不能静默改变比较母样本。其余14个在原母样本中；“在候选表中”也不保证能完成整块观测。

## 10. Codex 的实施与交付顺序

1. 核对当前专用分支及输入快照，导入本包数据，输出缺失项清单；不下载新母表覆盖旧快照。
2. 建立独立 `fixed_fields` 生成器/读取器，采用上面的母网格和永久 ID 契约。只接 TJN 场查询层。
3. 用本地全年缓存完成区域选择、对象动态落场、真实/名义足迹分层审计；将草案补成可冻结的第一版。
4. 在独立目录输出注册表 CSV+FITS+JSON、真实边缘/遮罩描述、邻接表和对象落场表。保留旧→新映射与来源 hash。
5. 给出两云场数、核心/候补/边界区域、覆盖差异、典型同时目标数、跨场次数、模板数/负担；逐一列出特殊目标，单列17365/4722。
6. 做有限接入烟测：选择L5冬季、L4春季、季节过渡与高Dec代表夜，用原曝光/overhead/科学条件选固定场，验证field_id重复运行和跨夜稳定。年度几何审计不等于新晨昏全年调度或天气Monte Carlo重跑。
7. 绘图至少包含：全年RA–Dec总览、L4/L5局部图、典型相邻行/高Dec边缘放大、TJ_test旧新对照、目标跨固定场轨迹。使用 `targeted_observing/plot_style.py` 的 Pal 5 风格；图中区分静态场库、当前活动场与目标轨迹。
8. 更新 `targeted_observing/PLAN.md`、`WORKLOG.md`，把本包项目备忘纳入该目录的持续文档；完成有意义的单一里程碑提交并报告commit。此次无需开展另外四项长期任务。

### 验收要点

- 坐标/编号/PA在重复运行、重排优先级、跨夜、跨campaign和追加场后保持稳定。
- 所有中心严格落在规定恒Dec行，RA跨零与高Dec正确；相邻场真实足迹与安全区域的缺口有明确审计。
- 全年每个被定义为可用的对象轨迹片段都能关联固定场，或者有明确未覆盖原因；不是只统计已排程赢家。
- 旧B几何覆盖、新候选轨迹覆盖、质量光变/唯一周期产额分别报告。原B P10/P50/P90=41/47/52不能自动沿用到固定网格。
- 旧0.15°留边语义已明确；真实PA和模板共同区域不被配置值代替。
- `TJ_test`保留历史身份，17365/4722的新场留场有星历证据。
- 已观察图像和今晚执行文件不被改写；其他scheduler模块及known/unknown未被顺手重构。

## 11. 本包复现

在解压目录运行，依赖 Python、NumPy、pandas、matplotlib：

```bash
python scripts/geometry_audit.py --mp inputs/mp_footprints.json --output audit/geometry_audit_rerun.json
python scripts/build_tjn_seed.py --source inputs --output regenerated_seed
python scripts/associate_tjn_clouds.py --inputs inputs --output regenerated_inputs
```

`scripts/associate_tjn_clouds.py` 可用原始B表、完整母表策略记录和特殊目标JSON重新生成分云中间数据，详见脚本参数。上述脚本不导入生产调度器、不连接望远镜、不生成可直接执行的整夜控制列表。

## 来源

- [专用分支9734af0a至7747cc16增量](https://github.com/xiaoyunao/SHARP/compare/9734af0a...7747cc1650fdeb255e677e496f50a376a99cb57b)
- [原MP数据](https://github.com/xiaoyunao/SHARP/blob/9734af0a/targeted_observing/data/mp_footprints.json)
- [旧TJN场构建](https://github.com/xiaoyunao/SHARP/blob/7747cc16/targeted_observing/annual_trojan.py)
- [旧自由中心候选](https://github.com/xiaoyunao/SHARP/blob/7747cc16/targeted_observing/annual_free_fields.py)
- [B原始场表](https://github.com/xiaoyunao/SHARP/blob/9734af0a/targeted_observing/docs/continuation_20261105/B_plot_data/fields.csv)
- [全母表选择函数](https://github.com/xiaoyunao/SHARP/blob/9734af0a/targeted_observing/docs/continuation_20261105/full_parent_selection_function.csv.gz)
- [今晚完整测试](https://github.com/xiaoyunao/SHARP/blob/7747cc16/targeted_observing/docs/FULL_NIGHT_TEST_20261111.md)
- [用户修改列表时角审查](https://github.com/xiaoyunao/SHARP/blob/7747cc16/targeted_observing/docs/EQUATORIAL_AUDIT_20261111.md)
- [Astropy完整WCS变换与footprint接口](https://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html)
- [JPL Horizons参考系与观测者星历定义](https://ssd.jpl.nasa.gov/horizons/manual.html)
