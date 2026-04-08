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
