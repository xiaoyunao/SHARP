# known_asteroid_local

本目录是 `known_asteroid` 的本地可运行版本。

默认路径：

- 输入根目录：`local_data/processed1/`
- 本地依赖：`resources/known_asteroid/`
- runtime：`local_data/known_asteroid/runtime/`

说明：

- 本地优先使用顺序版 `run_pipeline.sh`
- `slurm_*` 和 `submit_pipeline_slurm.sh` 仍保留，主要用于和服务器流程对照

单夜运行：

```bash
bash known_asteroid_local/run_pipeline.sh --batch false --submit-mpc false 20260406
```
