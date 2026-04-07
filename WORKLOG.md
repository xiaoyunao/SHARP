# WORKLOG

## 2026-04-07

- task: 初始化项目整理，建立本地/服务器双目录结构
- files_changed: `README.md`, `AGENTS.md`, `WORKLOG.md`, `PLAN.md`, `.gitignore`, `survey_local/*`, `known_asteroid_local/*`, `known_asteroid_server/*`, `heliolincrr_local/*`
- commands_run: `find`, `du -sh`, `sha256sum`, SSH 到服务器读取 `/pipeline/xiaoyunao/{survey,known_asteroid,heliolincrr}`，本地目录复制与整理
- key_findings:
  - `survey_v2_server` 的关键脚本与服务器哈希一致
  - `heliolincrr_v2` 的关键脚本与服务器哈希一致
  - `known_asteroid/pipeline_v2` 与服务器现版存在差异，已按服务器入口重建关键备份脚本
  - 服务器默认 Python 为 `Python 3.13.5`
  - `heliolinc` 环境 Python 为 `Python 3.10.20`
- validation:
  - 已完成服务器路径、脚本、解释器核对
  - `python3 -m py_compile survey_local/*.py heliolincrr_local/*.py known_asteroid_local/*.py known_asteroid_server/*.py survey_server/*.py heliolincrr_server/*.py` 通过
- remaining_issues:
  - 仍需初始化 git 并确认 GitHub 远端目标
  - `known_asteroid_server` 仍缺少部分服务器辅助脚本的完整回收
- next_step:
  - 本地提交第一次整理提交
  - 确认 GitHub 远端仓库
  - 后续调试时继续同步 `*_server` / `*_local`
