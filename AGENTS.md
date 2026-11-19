# Project AGENTS.md for `smt_asteroid`

## Scope

本文件只约束当前项目目录 `/Users/yunaoxiao/Desktop/smt_asteroid`。

## Source of truth

- 本仓库只保留服务器对齐版本，不再维护 `*_local` / `*_server` 双目录。
- 服务器基线路径固定为：
  - `/pipeline/xiaoyunao/survey`
  - `/pipeline/xiaoyunao/known_asteroid`
  - `/pipeline/xiaoyunao/heliolincrr`
- 仓库中的 `survey/`、`known_asteroid/`、`heliolincrr/` 默认分别对应上述服务器目录。
- 后续调试、重构和参数调整，都应优先保持与服务器行为一致。

## v2 服务器部署授权（2026-11-18 后续用户要求，优先于下方本地阶段限制）

- 用户明确要求在服务器完成v1归档/停止自动运行并部署首版v2，允许相应更新本地及服务器三个生产镜像目录；main分支仍不改。
- 新主入口为 `targeted_observing/server_service.py`，每日北京时间12点运行；每分钟tick处理持久事件与恢复。网站明确投递的单candidate JSON为最高优先级后随，名称 `FOLLOW_dddd`，当夜有效。北京时间17点或更早计划启动时刻截止受理，之后新请求保留当地时分转下一夜；开拍后不重排；服务器不自动从unknown产生后随、不做跨夜关联。
- v1完整归档及旧cron保存在 `/pipeline/xiaoyunao/v1/20261118_124436/`；v2配置/状态位于 `/pipeline/xiaoyunao/v2/`。公开列表、图件和小行星结果路径保持。部署/网站接口见 `targeted_observing/docs/SERVER_V2_20261118.md`。
- TJN按真实曝光/检出、跨夜基线、采样间隔和文献周期相位覆盖选场；计划曝光不得算入实拍历史。采样就绪不等于光变标定或唯一自转周期完成。
- MPC hold不解除，known恢复自动报需联系人确认及独立开关；unknown服务器自动报已退役，仅网站报。原始图测光、dashboard及无关EP任务保持。

## v2 运维更新（2026-11-19）

- scheduler维持北京时间12点自然运行；小行星处理独立于中午，在已发布列表end_local之后才监测L2，科学raw/L2对应及120秒稳定后入队。夜间不得把暂时无新文件当作观测结束。
- Windows同步script仅保留最近7个自然日的`YYYYMMDD_plan.txt`、`YYYYMMDD_plan_with_times.txt`及current_plan.txt；JSON/START_TIMES/分段等只放服务器存档，保持图件路径。
- 20261117已按用户要求标为异常夜，标记`/processed1/20261117/ASTEROID_DATA_INVALID.json`；L4删除、校验归档仅供审计，后续提取/历史/科学统计排除，raw/L1/L2保留。未经用户新指示不得自动解除。
- ASTORB/IERS周更使用reference_update：验证后原子切换、只留一份上一版，科学作业中不切换。服务器每周一10:20更新；项目独立定时任务每周一11:10新建对话核验与本地兜底，不绑定旧聊天。服务器网络失败可本地下载校验中转。JPL固定内核和冻结年度TJN轨迹的版本语义不混同ASTORB动态轨道表。
- digest2使用真实g/r/i/z，兼容历史G，保留官方解析器的波段换算；MPC hold不解除。仅收尾失败时允许原manifest验证后的resume-review，不重新提取整夜或重画已完成候选图。

## v2 方向与首阶段范围（2026-11-18 用户更新）

- 新策略为晨昏近太阳 MP + TJN + YSO。**YSO 例外，继续独立固定场**；此要求覆盖 11/16 的 YSO 共用网格约定。MP、TJN 及后续其他共用场任务采用下述 4397 场。
- 共用场名为 `MP_dddd`、`TJNL4_dddd`、`TJNL5_dddd`、`FOLLOW_dddd`；同一四位后缀必须对应相同 ICRS 中心。不得继续给 TJN 加 `MP_` 兼容前缀，不得直接沿用旧场后缀解释新位置。
- 冻结输入已纳入 `targeted_observing/data/shared_fields_20261116/`，与 11/16 用户选定 FITS 字节一致。YSO 保留 `YSO_ORI` / `YSO_CYG`。
- 未来服务器 v2 只做单夜已知提取和未知搜寻，输入覆盖全部合格科学数据（含 YSO/其他任务），不靠 `MP` 或四位数字通配符决定处理范围。服务器后随、跨夜关联和 unknown 自动上报迁到网站；单夜内轨迹关联保留。MPC 仍暂停；未来 known 服务器自动报，unknown 网站报，联系人待用户确认。
- 用户允许以后将服务器当前整套流程归档为 `/pipeline/xiaoyunao/` 下的 v1，新流程称 v2；这是对下方禁止平行版本目录的特定归档例外，不要求本地复制生产目录。
- 此前首阶段仅本地冻结/核查的范围已完成，后续部署以本文件上方“v2服务器部署授权”为准；首阶段历史见 `targeted_observing/docs/V2_FIELDS_AND_DAILY_LISTS_20261118.md`。

## 共用天区（2026-11-16 确认，11-18 更新 YSO 例外）

- 今后小行星巡天、TJN（木星特洛伊）监测和超新星监测共用兴隆 18° 筛选后的赤道平行网格；YSO 除外。基准表为 `targeted_observing/data/shared_fields_20261116/equatorial_footprints_xinglong18.fits`（原交付在 `outputs/equatorial_allsky_20261116/`），共 4397 场。
- 场身份由该表的 `field_id` 与 ICRS `center_ra`、`center_dec` 共同确定。`field_id` 是四字符文本 `0001`–`4397`（FITS `4A`，无前缀）。迁移旧任务时用目标坐标或所需覆盖范围重新匹配，历史编号需连同原天区表解释。
- 各任务可以选择各自的场子集和观测节奏，选中的场使用共享表的编号和中心。文件校验、复现命令和 TAN 审计见 [README 的共用天区说明](README.md#四类观测共用天区)。
- 这次确认适用于后续各任务更新。现有任务配置、已生成计划和服务器生产部署仍以各自更新记录为准。

## 本地新模式例外（2026-11-03 用户明确授权）

- `targeted_observing/` 是晨昏 + 木星特洛伊独立研发模式，在本地开发测试并同步 GitHub。
- 此目录不对应服务器路径，不连接/部署服务器，不改变现有三条生产流程。
- 允许此新增子目录；其余服务器镜像、MPC暂停、依赖及文档规则继续适用。
- 目标、验收和缺失现场输入以 `targeted_observing/docs/GOAL.md` 为准。

## 新模式生产隔离补充（2026-11-11 用户明确要求）

- 开发新模式时，不修改main分支，也不修改本地或服务器的survey/known_asteroid/heliolincrr生产代码。
- 速度/曝光中点等新适配只能写入targeted_observing/unknown_pipeline；新入口和新测试不得导入原生产实现。
- 380b777b曾在本地专用分支改动heliolincrr；已在后续隔离修正中恢复954e43d5原样。不得再次把新模式适配同步回生产目录。

## 临时调焦工具例外（2026-11-11 用户明确授权）

- 允许新增独立 `focus/`，保存程序和本地结果，不修改生产三目录或main。
- 可经默认SSH在服务器临时目录读取用户指定波段的调焦图并测量；原图/预处理图只读。
- 结果传回本地并校验后删除本次临时程序、日志、缓存和产物；不修改系统审计日志。
- 现场拍完确认前，排序须标明暂定，不能用无新文件自行推断结束。

## Environment rules

- `survey` 相关程序：
  - 服务器运行时使用默认 Python
  - 当前服务器默认 Python: `/home/smtpipeline/Softwares/miniconda3/bin/python`
- `known_asteroid` 相关程序：
  - 服务器运行时使用默认 Python
  - 当前服务器默认 Python: `/home/smtpipeline/Softwares/miniconda3/bin/python`
- `heliolincrr` 相关程序：
  - 必须使用 `heliolinc` 环境
  - 当前服务器解释器: `/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python`

## Layout rules

- 项目内只保留三套正式目录：
  - `survey`
  - `known_asteroid`
  - `heliolincrr`
- 不再新增 `*_local`、`*_server`、`v1`、`v2`、`server_v2` 之类平行目录。
- 如果服务器代码有变化，应直接同步到对应正式目录，而不是再建备份副本。

## Dependency rules

- `known_asteroid` 的大依赖文件按服务器布局放在：
  - `known_asteroid/astorb.dat`
  - `known_asteroid/de432s.bsp`
- `resources/known_asteroid/submit.xsd` 仅作为补充参考文件保留。
- 不要为了本地运行去改坏服务器绝对路径假设。
- 如果确实需要临时覆盖路径，优先通过命令行参数或环境变量覆盖，不要改默认服务器路径。

## Debugging rules

- 调试前先确认问题发生在：
  - `survey`
  - `known_asteroid`
  - `heliolincrr`
- 如果问题涉及真实行为、线上数据、定时任务、Slurm、运行环境或绝对路径：
  - 先检查服务器版本
  - 再把修复同步到本仓库
- 不要为了在本机上“更方便跑通”而引入偏离服务器行为的修改。

## Documentation rules

- 每次非平凡整理、调试、修复后，都更新：
  - `WORKLOG.md`
  - `PLAN.md`
- `WORKLOG.md` 和 `CHANGELOG.md` 都采用倒序维护：
  - 新日期块永远插在文件前面
  - 旧日期块必须保留，不要覆盖或替换
- `README.md` 保持项目总览，不写过程流水账。
- `WORKLOG.md` 写时间、目标、改动文件、关键命令、发现、验证、下一步。
- `PLAN.md` 维护当前目标、阻塞点、验证标准、近期步骤。

## Git rules

- 按正常 git 仓库管理，提交粒度保持在有意义的里程碑级别。
- 不提交大依赖文件：
  - `known_asteroid/astorb.dat`
  - `known_asteroid/de432s.bsp`
- 不提交运行产物、日志、缓存目录。
- 提交前确认没有把测试输出、旧目录或服务器临时文件误提交。

## Safety rules

- 不要把当前目录改成“便于本地运行但不再像服务器”的版本。
- 不要混用仓库相对路径和服务器绝对路径。
- 删除或重命名目录前，先确认引用已经同步更新。
