# 给 GPT Pro 的论文续写说明

请继续修改本压缩包中的完整 PASP LaTeX 工程 `manuscript/`。这是一份在
2026-08-03 数据快照和当前生产代码基础上完成的分析交接包。12 张正式图已经放入
`manuscript/figures/`，5 张正式表、统计摘要、代码冲突审计和作者待确认事项位于
`analysis_support/`。正式 PASP 编译入口是 `manuscript/manuscript.tex`；
`manuscript/manuscript_preview.tex` 只是预览入口。

## 你的任务

1. 通读 `manuscript/` 的全部正文、附录、数字宏和参考文献。
2. 先阅读下列证据文件，优先级从高到低：
   - `analysis_support/reports/FINAL_ANALYSIS_SUMMARY.md`
   - `analysis_support/reports/MANUSCRIPT_CODE_CONFLICTS.md`
   - `analysis_support/reports/AUTHOR_INPUTS_REQUIRED.md`
   - `analysis_support/reports/TASK_STATUS.md`
   - `analysis_support/tables/` 中的 Table 1--5
   - `analysis_support/qa/VALIDATION_REPORT.md`
3. 以这些冻结结果为准，系统修订正文、摘要、结论、附录、图注、表注和
   `paper_numbers.tex`。把已经有证据的旧 TODO 替换掉；真正需要作者、台站、上游或
   MPC 提供的信息仍保留为醒目的 TODO，不能自行猜测。
4. 把 12 张图按现有精确文件名接入正文。图的源 PDF 已在
   `manuscript/figures/`，不要重新绘图、改数或改变图中定义。
5. 根据 Table 1--5 决定正文表格的最终组织方式。可直接使用 `.tex`，也可为版式
   重排，但不能改变 CSV 中的统计粒度、状态和 guardrail。
6. 完成后请返回：
   - 可编译的完整修订 LaTeX 工程；
   - 编译后的 PDF；
   - 一份逐项修改记录；
   - 一份仍待作者回答的最短清单。

## 必须遵守的事实与表述边界

### 站点位置

- 当前代码并非“只有 960 m”。scheduler 和 known matcher 使用
  `117.575 deg E, 40.393 deg N, 960 m`；冻结的 short-arc orbit-confirmation
  使用 `117.575 deg E, 40.394239 deg N, 868.221 m`。
- L1/L2 header 中的 `GMG / 100.0313 / 26.6974 / 3227 m` 是上游模板污染，不能
  用作台站位置证据。
- 已完成 128/128 夜的替代 observer-location sensitivity run：全部 87,850 个
  linkage、正式 4,762 条和保留 58 条中，`fit_ok` 与 `is_good` 都是 0 次翻转。
- 这不是“只改海拔”的实验，因为纬度也改变了；也不能据此声称轨道元素稳定。
  部分三点极短弧轨道元素会大幅变化，因此禁止用这些 short-arc elements 做轨道
  种群推断。

### unknown pipeline 和人工复核

- 必须区分 `fit_ok` 与 `is_good`。正式 unknown catalog 的生产选择是
  `fit_ok && all_non_asteroid`；最终 4,762 条恰好同时全部为 `is_good=True`，但这是
  冻结样本上的等价，不能把两层定义写成同一个 gate。
- 人工链条是：`68 review-marked real -> 67 submission-selected -> 58 post-hoc
  retained`。唯一从 68 回撤到 67 的行是 `20260507 / 00001et`。
- `14,299` 是 linkage--detection memberships；全局唯一 detection keys 是
  `14,159`。不要再写成 “14,299 unique detections”。
- 58 行是 58 条 **single-night linkages**，不是 58 个确认的独立天体，更不是
  discovery count。
- 线性运动跨夜筛选得到 37 个 candidate components；它只是 screening result，
  不能等同于 37 个物理天体。
- 其中 6 条 linkage 在两遍 JPL 数值查询中都关联到 C/2025 Y1 (ATLAS)，但这不是
  MPC ingest、attribution、designation 或 discovery confirmation。权威 MPC 状态
  仍缺失。
- unknown ADES 冻结记录使用曝光开始时间；本地 179 行 midpoint correction 分析
  得到约 +15.001 s。是否向 MPC 沟通或更正必须由作者决定，不能在论文中暗示已经
  更正。

### known-object 结果

- known prediction 源行数为 13,311,871；按 exposure/object key 去重后名义分母为
  13,310,546，重复 1,325 行。
- 1 arcsec 匹配数为 534,780，名义匹配比例为 4.017716%。该分母是 `V<22.5` 加
  nominal WCS-projected detector rectangle，不是真实 photometrically detectable
  completeness denominator。
- 径向残差中位数为 0.263406 arcsec；这是受 1 arcsec 匹配阈值截断的 radial
  separation distribution，不能描述成无截断天测精度总体。
- blinded-known 的 75,899/102,347 = 0.741585 只是满足上游检测与几何条件后的
  conditional downstream-survival proxy，不是 image-level detection completeness。

### scheduler 与 operations

- 95.1743% 是 plan-active acquired cohort 内的 `7,455/7,833` acquired-frame plan
  compliance，不是 observatory efficiency，也不是完整 31,639 条计划的实现率。
- 历史 near-Sun 代码先按太阳距排序，随后在选择前把全候选重新按 RA 排序，因此
  实际是 RA-prefix 行为，不能写成“优先选择最小太阳距”。
- history feedback 只改变权重，不会禁止重复选择；而且生产 history 在
  2026-03-30 后基本停止更新。
- operations 的显式 normal-daily 完整事件链只有 2 夜。CPU/RAM、天气、设备故障、
  人工 override 和完整 observatory efficiency 没有充分证据，图表中的 N/A 不能
  改写成 0 或估计值。

## 写作要求

- 论文主线应定位为端到端巡天系统、实际运行、known-object 大样本天测/报告以及
  single-night unknown search 和人工复核，而不是“发现 58 个新天体”。
- 每个 headline number 必须保持其统计单位：exposure、detection、membership、
  tracklet、linkage、candidate component、object、night 不得混用。
- 对尚无外部证据的项目使用克制措辞，例如 `nominal`, `conditional proxy`,
  `candidate association`, `pending authoritative MPC reconciliation`, `N/A`。
- 不要直接照抄旧稿中的因果性表述。尤其不要把 twilight 分布归因于 near-Sun 模式，
  不要把计划缺口归因于天气或设备。
- 不要修改、重新计算或覆盖 `analysis_support/` 中的冻结证据。
- 若正文与证据冲突，以 `MANUSCRIPT_CODE_CONFLICTS.md` 和正式表格为准；若两个证据
  文件看似矛盾，请明确列为问题，不要自行选择更好看的数字。

## 尚需作者确认

至少保留以下未决项：正式 surveyed MPC 327 坐标、height type/datum；night-quality
mask 签字；正式台名/作者/单位/ORCID/CRediT/致谢；上游 reduction/depth/effective
mask/header 文档；58+9 的 MPC 逐行状态；6 条 C/2025 Y1 的论文处理政策；15 s 时间
更正政策；image-level injection；独立 CPU/RAM 和天气/设备日志。

请不要因为材料齐全就删除这些真正的外部 TODO。
