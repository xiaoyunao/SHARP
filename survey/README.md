# survey

服务器部署版每日调度程序。

## 行为

- 每天上午 `09:00` 运行一次。
- 当天运行日期记为当夜日期，例如 `2026-02-20`。
- footprint 默认路径：`/pipeline/xiaoyunao/survey/footprints/equatorial_footprints.fits`
- 历史曝光文件默认保存在 `workspace/history/exposure_history.fits`
- 如果历史文件不存在：
  - 直接新建空历史文件
  - 不回读 `/processed1`
  - 直接生成当夜观测列表
- 如果历史文件存在：
  - 尝试读取 `/processed1/YYYYMMDD/L2`
  - 只统计文件名里包含 `MP` 的文件
  - 从文件名中提取 `MP` 后面的编号，例如 `OBJ_MP_0023_0140_cat.fits.gz -> MP0023`
  - 每出现一次，对历史曝光次数加 `1`
  - 同一夜重复运行时，不会重复累计已经 ingest 过的文件
  - 如果目录不存在，或者没有 `L2`，或者没有任何 `MP` 文件，则跳过历史更新，直接生成当夜列表
- 如果设置了共享目录，还会额外复制：
  - `YYYYMMDD_plan.txt`
  - `current_plan.txt`
- `run_daily.sh` 默认会在基础巡天列表生成后调用 `survey.apply_followup`：
  - 只 ingest `FOLLOWUP_START_NIGHT` 及之后的 reviewed true unknown 源，默认 `20260624`
  - follow-up 目标写入 `MP_FU_<trkSub>_<field_id>`，仍会被后续 known/unknown 链条作为 MP 数据处理
  - 每个 active 源每晚插入 5 帧，满 2 个实际观测夜后完成
  - 若实际 L2 文件少于 5 帧，该夜记为 failed，后续继续排
  - 若从发现夜起满 10 天仍未完成，标记 abandoned
  - 插入失败时基础巡天列表保持可用，并在日志中写 WARN

## 输出

- `workspace/plans/YYYYMMDD_plan.json`
- `workspace/plans/YYYYMMDD_plan.txt`
- `workspace/history/exposure_history.fits`
- `workspace/history/l2_ingest_state.json`
- `workspace/followup/followup_state.json`
- `survey/logs/YYYYMMDD.log`
- `survey/plots/YYYYMMDD/`
  - `YYYYMMDD_cycle01.png`, `YYYYMMDD_cycle02.png`, ...
  - `YYYYMMDD_nightly_cycles.gif`
  - `YYYYMMDD_all_fields.png`
  - `YYYYMMDD_plots.json`

## 手动运行

```bash
python -m survey.run_daily \
  --date 2026-02-20 \
  --root-dir /pipeline/xiaoyunao/survey \
  --workspace /pipeline/xiaoyunao/survey/runtime \
  --footprints /pipeline/xiaoyunao/survey/footprints/equatorial_footprints.fits \
  --processed-root /processed1 \
  --publish-dir /path/to/shared/plan_dir
```

或者：

```bash
/bin/bash /pipeline/xiaoyunao/survey/run_daily.sh 2026-02-20
```

仅重跑 follow-up 插入：

```bash
PYTHONPATH=/pipeline/xiaoyunao \
/home/smtpipeline/Softwares/miniconda3/bin/python -m survey.apply_followup \
  --date 2026-06-25 \
  --workspace /pipeline/xiaoyunao/survey/runtime \
  --footprints /pipeline/xiaoyunao/survey/footprints/survey_fov_footprints_with_visibility.fits \
  --processed-root /processed1 \
  --review-root /pipeline/xiaoyunao/heliolincrr/review_packages \
  --publish-dir /pipeline/xiaoyunao/script \
  --state /pipeline/xiaoyunao/survey/runtime/followup/followup_state.json \
  --start-night 20260624
```

## 定时运行

参考 `survey/cron.example`：

```cron
0 9 * * * cd /pipeline/xiaoyunao && /bin/bash /pipeline/xiaoyunao/survey/run_daily.sh
```

## 部署说明

- 现在 `survey` 目录已经自包含，不再依赖其他 Python 包目录。
- 部署到服务器时，把 `equatorial_footprints.fits` 一起放到 `/pipeline/xiaoyunao/survey/footprints/`。
- 默认运行路径：
  - 包目录：`/pipeline/xiaoyunao/survey`
  - runtime：`/pipeline/xiaoyunao/survey/runtime`
  - processed：`/processed1`
  - 日志：`/pipeline/xiaoyunao/survey/logs/YYYYMMDD.log`
  - 图像：`/pipeline/xiaoyunao/survey/plots/YYYYMMDD/`

## 离线 slew 标定

`extract_slew_headers.py` 只读服务器原始 OBJ FITS 的主 header，不读取像素；
保留测试帧等完整序列，避免筛选科学帧后错误地跨帧配对。
`analyze_slew.py` 对连续科学曝光拟合两轴并行/串行线性、平方根、自由幂律和
加速—恒速—减速模型，输出分组交叉验证、块 bootstrap、最高速度剖面及诊断图。
这些工具不修改 scheduler 参数，也不由 daily 流程调用。

在仓库根目录执行（提取在服务器默认 Python 运行，拟合在本机运行）：

```bash
mkdir -p outputs/slew_calibration
ssh -p 20093 smtpipeline@www.xinglong-naoc.cn \
  '/home/smtpipeline/Softwares/miniconda3/bin/python - 20261021' \
  < survey/extract_slew_headers.py > outputs/slew_calibration/headers_20261021.json
python survey/analyze_slew.py outputs/slew_calibration/headers_20261021.json \
  --readout 14.47 --output outputs/slew_calibration/one_night
```

提取依赖 `astropy`；拟合另需本机已有的 `pandas`（表格和时间差）与 `scipy`
（稳健非线性优化），以及 `numpy`、`matplotlib`。这两个分析依赖不加入生产调度器的
必需依赖。可传入多份逐夜 JSON，输出中保留各夜与快照 SHA256。

时间差使用同一 `OBS_DATE` 字段，减前一帧 `EXPOSURE` 和平均读出时间；同一时间字段
的固定偏移会抵消，逐帧变化的延迟不会。`DATE` 整数秒和 `TEL_RA/DEC` 用于敏感性检查。
RA 小时乘 15 后处理跨零，不乘 `cos(Dec)`。得到的是相对跟踪的有效位移耗时关系，
不是编码器记录的机械运动曲线。最高速度剖面若在高速端平坦，不能把优化器返回值当作
测得的最高速度；常数项也不等同于独立测得的稳定时间。

小角度 RA 可能进入反复微调的慢分支。当前平方根模型强制单调，只描述样本中的正常转场；
其主导轴标签也是模型条件下的判断。小 RA 与大 Dec 同时移动时，微调延迟可能被
`max()` 遮蔽，良好的总耗时预测不能排除该分支。辨识它需要小 RA、小 Dec 的转场
覆盖，或控制/编码器时序；不能仅靠增加同类转场的重复数量。

多夜分析会额外做留整夜交叉验证，并按完整夜次 bootstrap，逐夜记录参数、主导轴
计数与预测误差；没有 OBJ 数据的日期可用空 `rows` 的提取 JSON 保留在来源清单中。
另做长科学曝光间隔的纳入敏感性检查，避免 180 秒常规转场上限决定拟合结论。
多夜最高速度敏感性会检查移除极少数大 RA 转场后参数是否仍受约束。

提取多夜 JSON 后可执行：

```bash
python survey/analyze_slew.py outputs/slew_calibration/headers_202610{15,16,17,18,19,20,21}.json \
  --output outputs/slew_calibration/week --micro-ra --bootstrap 500
```

`--micro-ra` 增加固定角尺度 0.1、0.3、1 度的候选附加延迟
`A*exp(-(dRA/scale)^2)`（仅正位移），允许 RA 耗时非单调；使用多个初值并做跨夜验证。
该函数只是探索性描述。`micro_delay_profile.csv` 记录固定尺度下的幅度剖面目标函数，
不是统计置信区间；被另一轴遮蔽时，优化器返回的幅度不能作为实测微调时间。
`small_ra_audit.json` 和 `small_ra_anomaly_repeats.csv` 同时保留小 RA 长耗时转场与
相同天区转场的重复对照，不能仅凭端点 header 将附加等待归因于某一轴。

`analyze_slew.py` 同时导出只含 RA、Dec 两个拟合子图的
`slew_axis_fits.png` / `slew_axis_fits.pdf`，保留两轴拟合图中的样本与曲线。

全年小 RA 检索使用冻结文件清单 `expanded/archive_inventory.json`
（`rows` 中保存 `night`、`directory`、`files`）。几何预筛选采用服务器的
`survey_fov_footprints_with_visibility.fits`，再以真实 header 核验，不能仅按文件名
推定实际指向。保留其他科学程序、同场重复、坏 header 和非标准指向各自的审计。

```bash
python survey/collect_slew_archive.py --source outputs/slew_calibration_20261022 --workers 8
python survey/investigate_small_ra.py --source outputs/slew_calibration_20261022 --use-archive --direct-ra-max 2 --output outputs/slew_calibration_20261022/expanded/ra_under2
```

收集器通过 SSH 只读主 header 的 2880 字节块直至 END，不读取图像像素；优先复用
已冻结快照，也可在前次收集停止后重跑以恢复，勿与仍在运行的同一收集任务重叠。
调查脚本的直接样本默认取非零 RA<2 度且 Dec<2 度；`--direct-ra-max 1` 可复现旧阈值，
`--output` 可另存图表而保留之前结果。扩大直接样本范围需要 `--use-archive`；不加时仅复用
最初14对RA<1度的定向探针。2026-10-22全年归档筛选得到70对、13夜、6种角度，
原RA<1度的14对完全保留；逐角度统计见 `isolated_small_ra_summary.csv`。
图中虚线为近期正常模型的参考外推，不是小角度数据强制单调拟合。
调查脚本使用同一夜、0.5 度 Dec 位移分层对照，分别比较 RA<2度与2–3度、1–2度与2–3度，
以及原来的 RA<1 度、<0.1 度和 <0.03 度对照，记录组内中位时间差、额外等待比例及整夜
bootstrap 区间。前两项本次分别为+0.029s（95%[-0.115,+0.182]）和-0.018s
（95%[-0.271,+0.187]）；没有检出总转场系统性额外变慢，不等于排除被Dec遮蔽的RA微调。
主比较仅使用
两端均约30秒曝光且符合正式天区指向的转场，长正间隔仍保留在等待分析中。
跨日期目录中相同时间、指向、曝光和通道的重复 header 会选择一个规范记录；
重复别名及其相邻转场不进入主比较，不在删除记录后重新连接前后帧。
这些是总转场时间的条件关联；不能把被 Dec 遮蔽的 RA 微调、人工/软件等待或
净位移为零的重复曝光直接解释为机械轴运动。零事件时不输出退化的零宽风险区间。
