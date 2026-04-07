# heliolincrr_local

本目录是 `heliolincrr` 的本地可运行版本。

默认路径：

- 原始输入根目录：`local_data/processed1/`
- 输出根目录：`local_data/heliolincrr/`
- Python：默认仍尝试 `heliolinc` 环境，可用 `PYTHON_BIN` 覆盖

单夜运行：

```bash
bash heliolincrr_local/run_single_night.sh 20260406
```

若本机环境不同：

```bash
PYTHON_BIN=/path/to/heliolinc/bin/python bash heliolincrr_local/run_single_night.sh 20260406
```
