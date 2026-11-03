# SHARP v0：交给 Codex 的天区设计与调度实施说明

日期：2026-11-03。状态：设计规格，尚未修改仓库、生成最终 RA/Dec 天区或在现场执行。配套 YAML 是要接入的新模式规格，不是现有代码已经支持的即插即用配置。

## 1. 本轮目标与默认值

将现有黄道条带巡天调度增加为可选择的新策略：**晨昏两段各最多 1 小时 + 每夜最多 2 小时木星特洛伊监测**。旧模式保留为 legacy。其他项目仅通过外部占用区间接入，本轮不设计 YSO、超新星或 proposal 调度算法。

- 晨昏：r30；每窗口两组、每组 4 场、每场 5 次；邻场模式重访 260 s；晚间 A→B，清晨 B→A。
- Trojan：r120；2 场；120 min 连续块；5 min 接入/交接总余量；重访 310 s、每场预算 22 点。
- 硬件预算：3°×3°；读出 15 s；邻场转场 20 s；换滤片 10 s；较远转场暂 60 s；全都按串行，直至日志证明可重叠。
- 日界、时间戳：以北京当地晚间日期标识一夜，清晨属于次日；计算和文件时间用 UTC，摘要同步输出 Asia/Shanghai。
- 初始 3 个有效夜固定 r，不自动颜色循环；从科学帧顺手估计精度和背景，不以完整 ETC/注入实验为启动前提。

分别完整阅读 `SHARP_Twilight_Plan_v0_20261103.md` 和 `SHARP_Trojan_Plan_v0_20261103.md`，科学约束和回退规则以二者为准；共享参数见 `SHARP_Defaults_v0_20261103.yaml`。

## 2. 代码基线与必须核对的现状

本次只读审查的 main 为 `2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757`，提交日期 2026-09-04。实施前重新取 HEAD，保留用户已有改动。本次读取了 config、astro_utils、scheduler、io_utils、run_daily、shell wrapper、followup、history_update 和 visualize_nightly。

| 位置 | 当前行为 | 新模式要做的事 |
|---|---|---|
| survey/config.py | 全局 30 s、15 s 读出、35 s slew、3 次/1800 s 重访、3 夜冷却；未按项目显式配置滤镜 | 加 program/config 输入；明确每个曝光波段与曝光；按新预算分项目 |
| survey/run_daily.py | 直接建立 SchedulerConfig、StripScheduler | 增加明确 `--strategy` 和 `--config`；新建分块模式，不只改旧常数 |
| survey/scheduler.py | normal 先填，末尾剩余不足 60 min 才 near_sun；无晚间独立块 | 先保护晚间/晨间，再分配 Trojan；每帧受硬截止约束 |
| survey/scheduler.py | near_sun 太阳角距排序后又被 RA 排序覆盖 | 使用明确角距和完整轮次可行性，检查东西侧、升落趋势 |
| survey/io_utils.py | TXT 第五列 W 硬编码；JSON 有 start_utc，但 TXT 没有等待/开始时刻 | 正确写滤片 token；将时间约束落实到真实执行接口或独立分组清单 |
| survey/run_daily.sh / followup.py | 默认开启后插 follow-up，按旧固定预算后移后续 start_utc | 新模式关闭这种无界后插；恢复任务在块内排，最终统一校验后发布 |
| history_update.py | 文件名历史只识别 MP 数字形式 | 保留兼容数字 field_id；program/campaign 单列；计划次数不当作已执行次数 |
| visualize_nightly.py | 依赖 footprint、cycle_id 和 start_utc | 新场写入足迹，展示场、目标、时间线、拒绝原因和未占用时间 |

已核查 TXT 的八列结构为：

```text
target_name  RA_hms  Dec_dms  0  FILTER_TOKEN  1  exptime_integer  0
```

只有 `FILTER_TOKEN` 要根据配置替换；其大小写/字符须取自现场已知合法命令或控制器说明。不能臆定 r 就是控制器的 r/R，也不能在缺失时静默回退到 W。配置中的 null 映射允许设计预览，但必须阻止导出可执行控制清单，直到实施者填入验证过的 token。

固定版本代码入口：

- https://github.com/xiaoyunao/SHARP/commit/2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757
- https://github.com/xiaoyunao/SHARP/blob/2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757/survey/config.py
- https://github.com/xiaoyunao/SHARP/blob/2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757/survey/scheduler.py
- https://github.com/xiaoyunao/SHARP/blob/2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757/survey/io_utils.py
- https://github.com/xiaoyunao/SHARP/blob/2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757/survey/run_daily.sh
- https://github.com/xiaoyunao/SHARP/blob/2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757/survey/followup.py

## 3. 天区设计的实际输出

### 3.1 晨昏场

为指定日期前后 3 个可用夜生成晚/晨各两个紧凑场组。以球面角距、真实升落和每次曝光高度为硬约束；优先 45–75° 太阳角距、50–65° 加权，黄纬候选 ±25°，低纬软加权。同组闭合循环每一跳都应符合邻场模型；初始邻场角距阈值设 3.5°，但这只是把约 3° 的场间距映射到预算的工程假设，需用现场日志更新。

不要默认晚、晨是同一批场，不要强行用全年固定黄道带。第 2、3 个可用夜每组优先保留 3 场，并允许几何不合格时更换。候选恢复通过更新星历/误差区生成独立覆盖需求。

每个建议场附：ICRS 场中心、焦面足迹、组内巡回次序、可开始区间、预计 5 次时刻、整段最低高度、真实太阳角距范围、月距、预计转场耗时、是否锚点、重叠后的独立面积。必须同时输出有解释的备选场，不能只给最终 8 个中心且无法复核。

### 3.2 Trojan 场

刷新 TJN 候选目录；对排名靠前的几十个目标交叉周期来源；按可测性、科学状态、两场共同 2 小时窗口、邻接关系以及未来 14 天轨迹联合选两个场。优先 2026 年秋冬的 L5，但最终服从当夜星历。

每场输出：全部候选和核心目标名单、r 或明确标记的 V 预测、周期来源/质量/备选解、科学角色、未来 14 天留场比例、坏区/星场污染提示、每晚可用块及最低高度。若以 V 作代理，不能偷偷在表头改称 r。

首个 14 天 campaign 场中心固定。重定心时，新 pointing_version 或新兼容数字 field_id 对应唯一坐标；保留 logical campaign 与前一版本的关联，不覆盖旧足迹。历史计数不得对 Trojan 使用旧的 3 夜冷却。

### 3.3 最小可追溯数据结构

```text
fields:
  field_id, pointing_version, program, campaign_id, group_id,
  ra_icrs_deg, dec_icrs_deg, footprint_version, position_angle_deg,
  valid_from_utc, valid_to_utc, priority, anchor, target_ids

targets:
  object_id, mpc_number, designation, population_source, population_snapshot_date,
  period_source, source_version, source_url, source_query_status,
  period_h, period_error_h, quality_raw, candidate_periods_h, science_role,
  brightness_value, brightness_band, brightness_source, brightness_is_proxy,
  ephemeris_source, ephemeris_epoch, retention_fraction, review_reason

exposures:
  exposure_id, block_id, group_id, program, campaign_id, field_id, pointing_version,
  cycle_id, visit_index, filter, controller_filter_token, exptime_s,
  planned_start_utc, planned_mid_utc, exposure_end_utc, readout_end_utc,
  not_before_utc, latest_start_utc, hard_stop_utc,
  preceding_slew_s, preceding_filter_change_s, center_alt_min_deg,
  footprint_alt_min_deg, sun_alt_max_deg, solar_elongation_deg, moon_sep_deg,
  planned_status, skip_reason
```

这是新模式推荐的逻辑字段，实施者可以适配现有 JSON，而不是创建不兼容的第二套含混字段。保留旧使用方需要的 `start_utc / mode / repeat` 等兼容字段并明确映射。`planned` 和 `actual` 必须分离。

## 4. 排程顺序与边界

1. 读取日期、站点、配置、版本化场表、目标星历、外部预约区间；生成太阳 −12°/−18° 穿越。
2. 生成晚间一小时、晨间最后一小时。已有外部预约是硬限制；若冲突则报告并缩短/跳过，不悄悄挤占。
3. 为晨昏分别挑选能完成五轮的组，晚间内侧先、晨间内侧后。不能只检查窗口中点或者“90 分钟有一半可见”。
4. 在去掉晨昏及外部预约后的 Sun<−18° 时段内寻找 Trojan 连续块；标准 120 min，90–119 min 两场，60–89 min 主场，少于 60 min 初版跳过。预算默认 7200 s/night，可为 0。
5. 每一曝光根据实际上一场/滤片计算开销，再检查曝光起止、必要读出和交接是否越过块截止；实际执行迟到时重新检查，不能机械追赶过期命令。
6. 同场不扣 slew，邻场按已验证模型，远场暂 60 s；滤片不变不扣换片时间。曝光中途不可因块截止主动截断，故开始前必须留足完成时间。
7. 月亮在地平线上时应用距月条件；明亮 Moon 门槛与高度余量均记录为初值，便于首轮调整。真实天气/设备约束由现场系统覆盖。
8. 出现高优先级恢复任务时，在预留余量或整组替代中处理；不使用旧后插模块把后面的整夜推迟。生成所有修改后的最终计划再统一校验、发布一次。
9. 输出未占用区间、已释放/跳过块、候选不足原因、预计科学目标数和实际曝光预算。即使没有合格 Trojan，也不能拿其他条带填满其预算。

## 5. 清单执行：先提供现场能够运行的版本

当前八列 TXT 不带时间。**不能把晚间、Trojan、晨间拼成一个 TXT，然后声称它会在中间等待其他项目，或者会等目标升起。**

初版最小执行交付可按执行组分别导出：`EV_A.txt`、`EV_B.txt`、`TJ.txt`、`AM_B.txt`、`AM_A.txt`，以及总 manifest。一个科学窗口内的两个组也要有独立的可开始时段；尤其清晨 A 组较低，不能因 B 组执行较快便提前开始。操作者按 UTC/北京时间 manifest 装载和交接即可开始试观测。

如现场已有经过验证的等待/暂停/停止/按需取下一条接口，Codex 接入该接口实现自动化，并在实际时钟上检查 `not_before/latest_start/hard_stop`。不发明新的望远镜指令。即使分组 TXT 手工执行，仍应记录实际 FITS 时间及错过的曝光，而不是把计划时刻当作完成时刻。

v0 脚本可一次生成全部分组文件和计划预览；人机交接方式不应阻碍首夜采集。但不得把“文件成功生成”报告成“已验证无人值守执行”。

## 6. 代表日期的窗口核算

下表北京时间，晨间均是次日；精确值由程序重新算。使用上述站点、Astropy 太阳位置、无折射，不含山体遮挡。

| 晚间日期 | 晚间窗口 | 次日清晨窗口 | 晚间 −12→−18 时长 |
|---|---|---|---:|
| 2026-11-03 | 18:05–19:05 | 04:42–05:42 | 31.8 min |
| 2026-12-21 | 17:51–18:51 | 05:24–06:24 | 33.2 min |
| 2027-03-20 | 19:20–20:20 | 04:14–05:14 | 32.4 min |
| 2027-06-21 | 20:59–21:59 | 02:24–03:24 | 49.7 min |
| 2027-09-23 | 19:05–20:05 | 04:00–05:00 | 32.3 min |

例如 11 月 3 日晚，扣除晨昏块后约 19:05–次日 04:42 可供共享，再取其中两小时给 Trojan，其他项目还可看到约 7 小时 37 分钟的未占用时段。这个数不扣天气、校准或已有预约，并不保证该时段全适合 Trojan。

## 7. 最小验收与三夜数据回收

这些验收针对可执行性，不要求先证明科学产率。

- 生成实际启动夜、春分和夏至代表夜；分别测试 Trojan 0/120 min、没有候选和外部时间冲突。所有曝光无重叠、太阳/高度合格、必要读出不越界。
- 手算核对同场、邻场、远场和换滤片四种开销；实际模型必须与 cadence/点数报告一致，循环回跳也计入。
- TXT 八列、坐标、曝光、已验证滤片 token 与 JSON 相符；缺失 token 阻止可执行导出，设计预览可继续。
- shell wrapper 后输出仍符合边界；新模式不会再走旧 follow-up 的无界后移。保留 legacy 基本生成能力。
- 场 ID/footprint/历史互相对应；真实 FITS 的 FILTER、EXPTIME 和时间记录与本次操作对应；使用实际帧更新历史。
- 第一夜正常执行并记录问题；第二夜保留主要场和参数；第三个有效夜后按执行完成率、恒星重复误差、背景、移动源恢复和目标留场情况调整一次。

每夜一页摘要：计划/实际窗口，计划/实际曝光，5 次 tracklet 完成数及短弧数，独立面积；Trojan 在场/可测/已查周期状态的对象数；按亮度的恒星重复误差；实际读出和邻/远场耗时；偏离计划及释放时间。可从既有日志和科学帧提取，不另起一个庞大的数据平台。

## 8. 实施完成的定义

Codex 应交付：版本化天区及核心对象表、可配置的新调度入口、JSON/分组 TXT/manifest/时间线、变更说明和有限验证结果。生成的未来日期清单明确标注有效夜和来源版本，失败时不把旧夜清单当作新结果发布。不得承诺本次规格已经得到某个新目标数、唯一周期数或发现率。

**按这套默认值开始 3 个有效夜，随后只根据已测到的瓶颈调整。** 具体参数改变要同时更新配置、预计 cadence 和场表，保持成果可比较。
