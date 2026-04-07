# smt_asteroid

本仓库整理了兴隆小行星相关三套流程，并明确拆成两份代码：

- `survey_local` / `known_asteroid_local` / `heliolincrr_local`
  - 用于本地开发和调试
  - 默认路径指向仓库内 `local_data/`
- `survey_server` / `known_asteroid_server` / `heliolincrr_server`
  - 作为服务器现版本备份
  - 保留服务器默认路径和运行假设

## 当前约定

- 后续调试默认以服务器版本为准：
  - `/pipeline/xiaoyunao/survey`
  - `/pipeline/xiaoyunao/known_asteroid`
  - `/pipeline/xiaoyunao/heliolincrr`
- Python 环境：
  - `survey` 和 `known_asteroid` 使用默认 Python
    - 服务器当前是 `/home/smtpipeline/Softwares/miniconda3/bin/python`
  - `heliolincrr` 使用 `heliolinc` 环境
    - 服务器当前是 `/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python`

## 目录结构

- `survey_local/`
- `survey_server/`
- `known_asteroid_local/`
- `known_asteroid_server/`
- `heliolincrr_local/`
- `heliolincrr_server/`
- `resources/known_asteroid/`
  - 放本地运行所需的大依赖文件
- `local_data/`
  - 本地测试数据根目录

## 本地运行

### survey_local

```bash
python3 -m survey_local.run_daily --date 2026-04-07
```

### known_asteroid_local

```bash
bash known_asteroid_local/run_pipeline.sh --batch false --submit-mpc false 20260406
```

### heliolincrr_local

```bash
bash heliolincrr_local/run_single_night.sh 20260406
```

如果本机没有 `heliolinc` 环境，需要显式指定：

```bash
PYTHON_BIN=/path/to/heliolinc/bin/python bash heliolincrr_local/run_single_night.sh 20260406
```

## 路径差异

- `*_server` 目录中的脚本保留服务器绝对路径，如 `/processed1`、`/pipeline/xiaoyunao/...`
- `*_local` 目录中的脚本改为仓库内路径
- `known_asteroid` 的大依赖文件本地统一放在 `resources/known_asteroid/`
  - 不提交到 git
  - 服务器版本仍按服务器路径假设运行

## 当前状态

- `survey_server` 的关键脚本已与服务器哈希核对一致
- `heliolincrr_server` 的关键脚本已与服务器哈希核对一致
- `known_asteroid_server` 已按 2026-04-07 服务器现版重建关键入口脚本

## GitHub

本地仓库会初始化为 git 仓库；推送到 GitHub 还需要一个明确的远端仓库目标。
