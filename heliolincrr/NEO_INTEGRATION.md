# digest2 / NEW NEO 接入设计

日期：2026-09-20。状态：已核对官方接口和仓库代码；尚未实现或部署。
本地 `main` 与 fetch 后 `origin/main` 同为 `2dffff5e`；服务器 SSH 连接超时，
因此尚不能确认服务器当前代码、环境和近期目标状态。

## 已确认的接口

- digest2 是可本地运行的短弧分类器。官方 Python 包提供 `Digest2`、
  `Observation` 和 `result.noid.NEO`；不需要每条轨迹请求在线评分服务。
- 按 MPC 来信及官方教程，使用 **NoID NEO score >= 65**，包括恰好 65；
  不用 Raw NEO 分数，也不把分数解释为真源概率或确认的 NEO 身份。
- 2026-09-20 GET 官方 PSV 提交表单确认以下实际 `obj_type` 值：
  - `neo candidate`：New NEO Candidate，对应来信要求的 NEW NEO。
  - `unclassified`：普通未分类对象。
  - `neocp`：已在 NEOCP 上的对象的后续观测。
  - `neo`：表单中单独的 NEO 选项，不等于 New NEO Candidate。
- 继续使用 `/submit_psv`，测试使用 `/submit_psv_test`。
  官方 curl 文档使用 `source=<file.psv`；现有实现使用 `source=@file.psv`，
  接入时应对齐并在测试接口验证，不把历史可用性当成当前格式保证。

来源：
- https://github.com/Smithsonian/mpc-public/tree/main/digest2
- https://docs.minorplanetcenter.net/tutorials/notebooks/mpc_tutorial_digest2/
- https://minorplanetcenter.net/submit_psv （直接 GET HTML 核对 radio 的 value）
- https://docs.minorplanetcenter.net/mpc-ops-docs/observations/command-line-submissions/
- https://docs.minorplanetcenter.net/mpc-ops-docs/astrometry/reporting-observations/

## 当前流程与最小改动

当前：known subtraction -> unknown catalog / trkSub -> GIF/review package ->
人工 submit CSV -> watcher -> export -> validate -> submit。
`export_unknown_ades.py` 默认 `--obj-type NEO`，整夜一个 PSV；
`submit_reviewed_unknown.py` 未提供分类分流；watcher 只记录整夜结果。

建议：

1. 新增 `score_unknown_digest2.py`，在 `run_single_night.sh` 分配 trkSub 后、
   打包 review 前为每个未知单夜目标评分。一个 trkSub 对应一组按时间排序、
   按原始 image/object ID 去重的观测，不对每一观测点单独评分。
   检查观测覆盖完整性、有限数值、至少两个不同历元和台址 327；
   缺行、未知台址、异常返回必须显式标为评分失败。
2. 与 `export_unknown_ades.py` 共用观测行构建逻辑，保证评分与最终提交的
   UTC 时间、RA/Dec、误差、星等和滤光片一致。可将这些行转换成 Python
   `Observation`，使用 ADES 误差模式；不要假设 `classify_file()` 支持 PSV，
   官方当前文档列出的是 MPC80 / ADES XML。
3. 在 `package_unknown_review.py` 的 manifest 添加 `digest2_neo_noid`、
   `digest2_neo_raw`、`neo_candidate`、`digest2_status`、评分版本及输入 hash；
   保持现有两列复核 CSV 兼容。网页消费者显示分数并优先展示高分候选；
   网页实现不在本次检查的入口代码中，部署时需确认实际消费者。
   高分不自动设置 `is_real=1`。
4. 提交前用最终观测快照再次核对评分输入 hash；发生变化必须重算。
   对人工接受的目标按下表拆批，每个目标的所有观测完整归入同一批：

   | 判定 | 动作 |
   | --- | --- |
   | is_real=0 / 尚未完成复核 | 不提交 |
   | is_real=1，未编号且 NoID NEO >=65 | `neo candidate` |
   | is_real=1，未编号且 NoID NEO <65 | `unclassified` |
   | is_real=1，评分失败或缺失 | 暂停该目标提交，显式报告错误，不能当作低分 |
   | 已获 MPC designation / NEOCP ID | 按现有身份走已知或 NEOCP 后续观测流程 |

5. `submit_reviewed_unknown.py` 协调拆批，`export_unknown_ades.py` 负责共用
   观测导出和请求；同时覆盖 direct exporter / single-night 的提交入口，
   防止绕过评分门控。普通批和 NEO 批分开 PSV、validate reply、submit reply。
   保留全部人工真源的统一 masked JSON/FITS，避免破坏 `survey` follow-up 摄入。
6. `watch_submit_reviews.py` 保留现有跨进程锁，新增不可变提交清单和逐批账本：
   night、trkSub 集合、review 内容 hash、PSV hash、obj_type、尝试状态、
   CurlID、时间与回执。逐批成功后立即落盘；另一批失败不重交成功批。
   网络超时且无法判断 MPC 是否已接收时标为待核对，不能盲目自动重发。
   review manifest 增加分数等元数据不应触发已提交观测重新上报。

## 环境和测量约束

- 在服务器 `heliolinc` 环境中检查 Python 兼容性后固定 digest2 版本；
  官方 Python 教程要求 Python >=3.9。若环境不兼容，评估调用 C CLI，
  不为此直接升级生产 heliolinc 环境。
- 缓存 ObsCodes（必须包含 327）、population model 和 MPC.config，记录
  版本/hash，固定 `repeatable=True`。不要每个目标实时下载配置。
- 保留官方误差配置/误差下限；本项目的质心误差不等于完整天测误差，
  不人为调小误差以得到更高分。
- 项目已有审计指出 unknown exporter 直接使用 L2 起始 MJD、比曝光中点
  早约 15 s。上线前确认当前服务器是否已修复，并统一新数据的时间约定；
  历史时间更正和已提交目标的处理应独立跟踪，不自动覆盖/重交。
- 先覆盖当前可提交 unknown 集合；随后审计上游 `fit_ok` 等筛选是否会漏掉
  真实高分短弧。不要把短弧最佳拟合轨道的 q 值当作 digest2 的替代。

## 验收与历史对象

- 先在临时目录回放，仅生成分数与拆批清单，不能覆盖生产复核记录。
- 检查 64.9 / 65 / 65.1、评分异常、空批、缺观测、真假混合、输入变化。
- 检查两批不重叠且并集恰为所有评分成功的人工接受目标。
- 模拟一批成功另一批失败、多 watcher 并发和连接结果不明确；确保成功批
  不重发，测试回执不能被标为正式 submission，HTTP 200 不等于业务接收成功。
- 两批分别通过 `/submit_psv_test`，保存并检查格式验证结果；正式端点不能
  用作测试。已有 submitted 状态应迁移保留，不能因新增评分重放历史夜次。
- 来信列出的七个 trkSub 作为追溯样本；需服务器数据才能计算本地分数。
  2026-09-04 WORKLOG 已记载 `00001nP`, `00001o2`, `00001oI`, `00001oj`
  为假源；`00001nF` 仅保留为可能真实。`00001vy` / `00001ua` 本次尚未核实。
  不因 MPC 高分提醒推翻人工假源结论，也不自动将历史目标重新提交。

实施顺序：服务器基线核对 -> 固定依赖/评分回放 -> review 展示 -> 拆批与账本
-> 测试接口验收 -> 新夜上线。首次上线不回扫已完成的历史提交。
