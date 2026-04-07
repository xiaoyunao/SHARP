# survey_local

本目录是 `survey` 的本地可运行版本。

默认路径：

- 包目录：仓库内 `survey_local/`
- footprints：`survey_local/footprints/equatorial_footprints.fits`
- processed 根目录：`local_data/processed1/`
- runtime：`local_data/survey/runtime/`
- publish 目录：`local_data/publish/survey/`

运行方式：

```bash
python3 -m survey_local.run_daily --date 2026-04-07
```

或：

```bash
bash survey_local/run_daily.sh 2026-04-07
```
