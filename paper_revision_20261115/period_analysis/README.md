# SHARP 年度光变周期：离线移植与候选筛查

独立研究目录；不改变 survey / known_asteroid / heliolincrr、服务器任务、论文年度统计或已发表结论。这里只给出周期搜索假设，尚无确认的新自转周期。

## 交付来源与实际移植范围

找到了 `period_delivery_20261127` 的已解压副本：

`/Users/yunaoxiao/Desktop/gotta_asteroid_1/lightcurve_revision_20261127/handoff/input/period_delivery_20261127/`

Downloads 没找到同名 ZIP；macOS 拒绝终端读取 `.Trash`。没有恢复、删除或修改原交付。该副本含 97 个目标的点表、曲线、原程序和核验记录。11 个源文件的 SHA256 全部通过；指纹在 `validation.json`。

原程序链路是批处理 → 宽周期 LS/PDM → 加权傅里叶拟合 → P/2P 局部核验 → 人工选择 → 导出。原入口使用 TMTS/SiTian 绝对路径、已做几何和光行时修正的光变，以及每夜至少 60 点的筛选；不适合直接运行在 SHARP。

本目录重新实现独立的 NumPy/SciPy/Astropy 数值核心，沿用其加权正弦/余弦傅里叶模型。没有导入整个 16,127 行上游模块、运行发布脚本或接入其他项目数据库。`validate_port.py` 仅提取已审读的纯数值函数 `phase_fold_and_plot_P_only.fourier_fit_fixed`，与新实现对同样的 97 个目标、44,111 个最终拟合点逐一比较；最大系数差 `4.63e-14`。这验证固定周期拟合器的等价性，不表示搜索/人工判定流程完全相同。

为适应稀疏采样，这里有意作了以下调整：

- 每个试频都拟合各夜常数和二阶傅里叶系数，LS 仅提供额外候选。完整傅里叶网格参与搜索，避免先减夜均值后只依赖 LS 峰。
- 默认搜索 1–120 h，频率网格间隔不大于 `1/(8T)`；保留半频、倍频、±1 cycle/day 混叠和边界候选。不声称这一范围包含全部物理自转周期。
- 不使用目标已知周期给搜索提供种子；不随候选周期改变剔除点。误差按 `sqrt(mag_err²+0.03²)` 加固定探索性底噪，尚未由比较星重新估计。
- 对 2P 使用四阶模型与仅含第 2、4 谐波的嵌套模型比较，避免用阶数不足误排除倍周期。局部小误差和单个最优峰都不代表周期唯一。
- 没有直接移植 PDM 的 10-bin/3-cover 门槛：SHARP 每夜往往只有几次测量，默认采用能明确检查秩与夜零点退化的加权模型。PDM 可留作后续密采样样本的独立诊断。
- 留一夜检查只在训练夜内局部重拟合频率；测试夜额外拟合一个亮度常数，因此输出是**给定该夜均值的形状检验**，不是完全独立的预测。它也不是对全频率搜寻的假阳性校准。

## 年度输入与筛选

使用已经冻结的第一年标准 30 s 白光 MP 样本，UTC 曝光中点范围 `[2025-11-15 04:00, 2026-11-15 04:00)`。不是把全部一年数据强行折到同一光变形状上。`prepare_annual.py` 按目录夜/文件名连接年度曝光台账，采用台账 UTC 曝光中点和北京时间正午划分的物理夜；输出保留逐点源文件、obj_id、exposure_key 和输入行号。

预筛：年度 quality、编号身份、Flag=0、位置残差≤1″、正通量/误差且 SNR≥10、9.5≤mag≤22、0<magerr≤0.15、0.001≤帧零点散布≤0.1 mag；排除一检出多身份和同目标同曝光多匹配。每个目标只选时间上最密的≤30 d窗口，然后保留每夜≥3点的夜；至少40点、4夜、总点数≥夜数+12才进入搜索。窗口选择不看拟合结果。这里只检查每个目标一个窗口，其他季节和其他窗口仍可有科学价值。

825,431 条缓存关联中，346,645 条通过探索性预筛，涉及34,822个编号目标；520个达到窗口门槛，共30,160个窗口点。数量不代表成功测出周期。逐项首个剔除原因、输入指纹、采样目录和逐输入行剔除记录见 `results_20261127/`。未编号身份暂缓，不是判定它们不能测周期。

目前缓存没有 `Flag_ISO`、`Flag_ISO_Num`，且距离列全部为空，没有距离/相位角/光行时修正、比较星差分测光或图像检查。自由夜零点仅是探索性处理，会削弱长周期约束，也不能替代几何修正和逐帧标定。半周期/倍周期、时间窗口效应和系统误差仍需审查。所有输出都是观测者时刻的候选光变周期，不能当成精确恒星时自转周期。

## 文件与复现

依赖沿用年度 Python 环境：numpy、scipy、pandas、astropy、matplotlib、pyarrow。目录内运行：

```sh
python prepare_annual.py --output /tmp/sharp_period_inputs
python validate_port.py --output /tmp/sharp_period_validation.json
# 原交付仍存在时，可加 --delivery /absolute/path/period_delivery_20261127 做97目标等价核验
python run_pilot.py --inputs /tmp/sharp_period_inputs --output /tmp/sharp_period_pilot --limit 10
python run_pilot.py --inputs /tmp/sharp_period_inputs --output /tmp/sharp_period_screen --limit 520 --plots-limit 0 --cv-targets 0
# 所有盲搜结束后才可单独联网查询；不改变搜索结果
python check_references.py --results /tmp/sharp_period_pilot
```

`run_pilot.py --objects N:9723 N:7835 ...` 可选择采样目录内的指定对象。默认只为前10个生成 PDF/PNG、频谱及留一夜诊断；`--plots-limit` / `--cv-targets` 可调整。初筛全部目标不做逐目标图像复核。`period_hypotheses.csv` 保存全部搜索候选，`pilot_summary.csv` 的 `candidate_period_h` 仅是排名首位的光变假设；`competing_delta_bic_lt10=0` 也不表示排除了物理倍周期。`cv_nights=0` 表示没有做或无法做该检验，不能视为通过。

`check_references.py` 通过 [NASA/JPL SBDB Query API](https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html) 查询已有 `rot_per`，保留原始 JSON、URL、获取时刻和 API 版本。该库的空值只表示 SBDB 没有这项记录，不能证明首次测量；最终需逐一核对 [LCDB](https://www.minorplanet.info/PHP/lcdb.php) 及论文。公开周期仅用于盲搜后的对照，不用于挑选峰或事后改标可靠等级。

## 本轮完整结果

详见 [2026-11-28 结果报告](RESULTS_20261128.md)。`screening_20261127/` 保存520目标盲搜，`review_20261128/` 保存16目标的图件和条件留一夜复核；`results_20261127/` 是前10个试跑和全体输入采样/剔除台账。日期保留交付命名，复核在北京时间11/28完成。

完整筛查后执行 `python check_references.py --results /tmp/sharp_period_screen`，再执行 `python prioritize.py --results /tmp/sharp_period_screen`。后者的规则写在 `screening_summary.json`，只决定检查顺序；`missing_sbdb_review_queue.csv` 即16目标的选择依据。将其中 object_key 通过 `run_pilot.py --objects ... --limit 16 --plots-limit 16 --cv-targets 16` 传入即可重建重点复核。留一夜和图件不参与原先520目标的盲搜。

## 服务器恢复后的实际步骤

2026-11-27 使用项目 EasyTier SSH 入口试连一次，超时后停止。等待用户明确告知恢复；没有部署、提交 MPC 或更改服务器文件。

1. 只读检查现有 L2/L4：按这里的 exposure_key/source_file/obj_id 补全三个质量标记、全部测光孔径/误差、坏像元/饱和/混叠标记和真实曝光时刻。年度论文输入保持不动，另存新的带指纹周期输入。
2. 对候选与已知周期对照源提取周围比较星、图像切片，检查逐帧零点、位置/坏像元、背景与拥挤、孔径依赖和夜内系统误差。把可能污染点的排除理由固定后再比较所有周期。
3. 用冻结轨道/星历计算逐时刻距离、相位角和光行时；验证近似修正误差。分合理观测季节拟合，比较共享/夜零点和多个傅里叶阶数，检查谱窗、日混叠及P/2P；重做跨夜与时间子集检验。
4. 再核对 LCDB/原论文，分为已有周期恢复、改进/争议周期、首次候选与证据不足；对多解给候选集与原因。只有这些完成后，才能给出可用于论文的新周期与误差。
