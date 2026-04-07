# Project AGENTS.md for `smt_asteroid`

## Scope

本文件只约束当前项目目录 `/Users/island/Desktop/smt_asteroid`。

## Source of truth

- 后续调试默认以服务器代码为准，不以本地副本为准。
- 服务器基线路径固定为：
  - `/pipeline/xiaoyunao/survey`
  - `/pipeline/xiaoyunao/known_asteroid`
  - `/pipeline/xiaoyunao/heliolincrr`
- 本地仓库内的 `*_server` 目录用于保存服务器版本备份。
- 本地仓库内的 `*_local` 目录用于保留本地可运行版本。

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

- 项目内始终保持两套目录：
  - `survey_local` / `survey_server`
  - `known_asteroid_local` / `known_asteroid_server`
  - `heliolincrr_local` / `heliolincrr_server`
- 不再新增 `v1`、`v2`、`server_v2` 之类平行目录。
- 新改动若来自服务器，应先更新 `*_server`，再决定是否同步到 `*_local`。

## Dependency rules

- `known_asteroid` 的本地依赖文件路径和服务器不同，必须明确区分。
- 本地依赖统一约定放在：
  - `resources/known_asteroid/astorb.dat`
  - `resources/known_asteroid/de432s.bsp`
  - `resources/known_asteroid/submit.xsd`
- 服务器版本脚本保留服务器原始路径假设，不要为了本地运行去改坏 `*_server`。
- 本地版本脚本要显式改成仓库内路径，不能继续默认写 `/processed1` 或 `/pipeline/xiaoyunao/...`，除非该参数就是给服务器用的覆盖项。

## Debugging rules

- 调试前先确认问题发生在：
  - `survey`
  - `known_asteroid`
  - `heliolincrr`
- 如果问题涉及真实行为、线上数据、定时任务、Slurm、运行环境或绝对路径：
  - 先检查服务器版本
  - 再决定本地副本是否需要同步
- 如果修改了服务器基线逻辑，本地 `*_server` 必须同步。
- 如果修改只为本地便捷运行，可只改 `*_local`，但要在 `WORKLOG.md` 记清楚。

## Documentation rules

- 每次非平凡整理、调试、修复后，都更新：
  - `WORKLOG.md`
  - `PLAN.md`
- `README.md` 保持项目总览，不写过程流水账。
- `WORKLOG.md` 写时间、目标、改动文件、关键命令、发现、验证、下一步。
- `PLAN.md` 维护当前目标、阻塞点、验证标准、近期步骤。

## Git rules

- 初始化后按正常 git 仓库管理。
- 不提交大依赖文件：
  - `resources/known_asteroid/astorb.dat`
  - `resources/known_asteroid/de432s.bsp`
- 不提交运行产物和本地测试数据目录。
- 提交前确认没有把旧版本目录、缓存目录、运行输出误提交。

## Safety rules

- 删除旧目录前，必须确认对应内容已经迁移到新的 `*_local` / `*_server`。
- 不要把 `*_server` 改成“看起来更好用但不再像服务器”的版本。
- 不要混用本地路径和服务器路径。
