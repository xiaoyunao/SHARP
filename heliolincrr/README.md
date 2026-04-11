# heliolincrr

本目录保存服务器当前 `heliolincrr` 代码，并作为仓库中的唯一正式版本。

说明：

- 保留服务器路径和环境假设
- 默认使用 `heliolinc` 环境
- 后续调试和调参都直接基于这一套代码进行
- 单夜正式入口为 `run_single_night.sh`
- 15 夜正式入口为 `run_pipeline_15.sh`
- 单夜统计脚本为 `summarize_single_night.py`
- 未知单夜 link 的可视化已并入 `run_single_night.sh`，底层脚本为 `plot_unknown_links.py`，输出到 `heliolincrr/plots/YYYYMMDD/`
