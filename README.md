# smt_asteroid

本仓库现在只保留一套代码，目标是直接镜像服务器上的实际运行版本，不再维护本地可运行分支。

## 当前约定

- 仓库内目录默认对应服务器路径：
  - `survey/` 对应 `/pipeline/xiaoyunao/survey`
  - `known_asteroid/` 对应 `/pipeline/xiaoyunao/known_asteroid`
  - `heliolincrr/` 对应 `/pipeline/xiaoyunao/heliolincrr`
- 后续修改以服务器行为为准，不再为了本地便捷运行改写路径假设。
- Python 环境：
  - `survey` 和 `known_asteroid` 使用服务器默认 Python
    - 当前为 `/home/smtpipeline/Softwares/miniconda3/bin/python`
  - `heliolincrr` 使用 `heliolinc` 环境
    - 当前为 `/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python`

## 目录结构

- `survey/`
- `known_asteroid/`
- `heliolincrr/`
- `resources/known_asteroid/`
  - 放补充参考文件，例如 `submit.xsd`

## 服务器运行

### survey

```bash
cd /pipeline/xiaoyunao/survey
/home/smtpipeline/Softwares/miniconda3/bin/python -m survey.run_daily --date 2026-04-07
```

或者：

```bash
cd /pipeline/xiaoyunao/survey
bash run_daily.sh 2026-04-07
```

### known_asteroid

```bash
cd /pipeline/xiaoyunao/known_asteroid
./submit_pipeline_slurm.sh --batch false --submit-mpc false 20260406
```

### heliolincrr

```bash
cd /pipeline/xiaoyunao/heliolincrr
PYTHON_BIN=/home/smtpipeline/Softwares/miniconda3/envs/heliolinc/bin/python bash run_single_night.sh 20260406
```

当前单夜流程已改为：

- `mask_gaia -> make_tracklet -> merge_tracklets_night`
- `run_linear_links_from_tracklets.py` 生成 `/pipeline/xiaoyunao/data/heliolincrr/<night>/rr_links`
- 默认单夜参数固定为 `speed=5 arcsec/h`、`direction=10 deg`、`require_shared_endpoint=True`

当前 15 夜流程仍使用 `run_rr_from_tracklets.py --profile w15`。

## 依赖与路径

- 代码中的服务器绝对路径保持不变，例如 `/processed1`、`/pipeline/xiaoyunao/...`
- `known_asteroid` 的大依赖文件按服务器布局放在 `known_asteroid/` 下：
  - `astorb.dat`
  - `de432s.bsp`
- 这两个大文件默认不提交到 git

## 当前状态

- `survey` 的关键脚本已与服务器哈希核对一致
- `heliolincrr` 的关键脚本已与服务器哈希核对一致
- `known_asteroid` 已按 2026-04-07 服务器现版重建关键入口脚本
- 2026-04-08 起，仓库已移除全部 `*_local` 目录

## 版本记录

- 仓库级更新说明见 `CHANGELOG.md`
