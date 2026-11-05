# 2026-11-05 续接审查：结果与复现入口

已读取续接对话、执行清单和用户下载的 `SHARP_resume_20261105.zip`。本轮全部为本地开发和条件模拟；未连接服务器、未部署、未生成望远镜执行文件、未提交MPC。原基线固定为 `a5cb9454e56771ebe4211b7a544412056cce4c40`，原报告及缓存保留。

**推荐下一轮以 B 作为通用采样比较基准；如果明确优先 Thymbraeus／Paeon，可选 C 的10%科学份额。D 的20%目前代价偏高。晨间十访适合作为面积换短弧完整性的可选模式，没有证据表明它改善了低距日角覆盖。** 这些选择没有自动启用现场调度。

## 交付

- [365夜可按日期浏览的图索引](atlas/index.html)：本地打开；包含1张全年足迹、12张月图、365张高度与曝光时间图。GitHub HTML页面显示源码，下载仓库后打开即可。
- [13页全年／月度PDF](atlas/nominal_annual_atlas.pdf)：原基线的名义模拟排程，不是已完成观测。
- [原始365夜绘图数据](baseline_plot_data/manifest.json)：四个轻量文件，6137个夜次—块—场记录、37078条曝光；全部365夜输入哈希、配置及来源提交。无整夜零曝光，但没有Trojan的夜次也保留。
- [Trojan详细对比](TROJAN_CONTINUATION.md)、[晨昏详细对比](TWILIGHT_CONTINUATION.md)。
- [配对年度统计](paired_strategy_summary.csv)、[200次配对质量计数](paired_realizations.csv)、[特殊任务状态](special_task_summary.csv)、[16候选可见性及拒绝理由](special_availability_audit.csv)。
- `A/B/C/D_nightly_ledger.csv.gz`：各73000条时间与选择记录，含开始时刻、场ID、去重对象ID、特殊任务及实际开销；`*_science_tasks.csv.gz`记录全部16对象×200实现，未知物理完成状态留空。
- [全母表选择函数](full_parent_selection_function.csv.gz)：原16445对象、原筛选拒绝理由、各策略被采样次数与质量达标频率，零次也保留。未校正发现偏差，不是总体无偏推断。
- [统一时间账本](unified_time_ledger.csv.gz)：四策略×200×365，整窗晨昏、Trojan、外部预约、余量分列；单独校准额外时间为原基线的0，Trojan300秒交接已计入两小时。实际校准计划尚待现场输入。
- [51项回归及年度验收记录](validation/summary.json)、[菜单/星历/目录输入哈希](annual_input_manifest.json)。新B名义年365夜、37264曝光独立5秒验收全部通过。
- [配置快照](strategy_configs.json)、[对比来源哈希](comparison_provenance.json)、[晨昏来源哈希](twilight_provenance.json)、[原执行清单](PRO_REQUIREMENTS.md)、[原图册审查说明](ORIGINAL_PLOTTING_REVIEW.md)。原恢复图不冒充本轮计算图。

## 复现

使用本仓库已有 `requirements-annual.txt` 环境；绘图使用 matplotlib。PDF打包另外需要 Pillow、reportlab（本次使用应用附带的Python环境，不改变服务器依赖）。原始星历、IERS与年度菜单仍保存在忽略提交的 `cache/`、`products/`；精确数值复现需要保留同一输入快照。

```bash
python3 -m unittest discover -s targeted_observing/tests -v
# P0：相同旧菜单，仅修复去重排序
python3 -m targeted_observing.annual_simulation --fields 4 --dwells 2 --scenario baseline optimistic --policy-version 4 --save-curves --workers 2
# 文献权重与完成预测独立：--science-strength 0 不再隐含关闭预测
# 只有显式 --no-completion-forecast 才关闭预测
python3 -m targeted_observing.annual_visibility --free-centres --workers 3
python3 -m targeted_observing.annual_simulation --continuation A --fields 4 --dwells 2 --scenario baseline --workers 1 --save-curves
python3 -m targeted_observing.annual_simulation --continuation B --fields 4 --dwells 2 --scenario baseline --workers 1 --save-curves
python3 -m targeted_observing.annual_simulation --continuation C --fields 4 --dwells 2 --scenario baseline --workers 1 --save-curves
python3 -m targeted_observing.annual_simulation --continuation D --fields 4 --dwells 2 --scenario baseline --workers 1 --save-curves
python3 -m targeted_observing.continuation_report
python3 -m targeted_observing.continuation_figures
python3 -m targeted_observing.annual_morning_alignment --workers 3
python3 -m targeted_observing.annual_coverage --visits 10 --root-name twilight_morning10_aligned --skip-injections --nside 1024
# 将20分别替换为19、21，运行三种深度情景
python3 -m targeted_observing.annual_coverage --visits 10 --root-name twilight_morning10_aligned --causal --skip-area --limiting-mag 20
python3 -c "from targeted_observing.annual_effective_area import run; run(10,'twilight_morning10_aligned')"
python3 -m targeted_observing.continuation_twilight_report --geometry-only
python3 -m targeted_observing.continuation_twilight_report
# B的独立名义逐曝光合同；独立目录，不能覆盖原selected_nominal
python3 -m targeted_observing.annual_simulation --continuation B --fields 4 --dwells 2 --scenario baseline --workers 1 --nominal
python3 -m targeted_observing.annual_delivery --source-name nominal_continuation_B_baseline --output-name selected_nominal_B_20261105 --visibility-name visibility_free
```

P0-A使用恢复的四个 `review_tools`。原导出在HEAD仍为a5cb945时完成；以后不能用新HEAD伪装生成旧排程的代码。可在旧提交工作树中导出，或用独立产品副本及真实来源声明。导出器缺日期、缺坐标、缺来源时硬报错。图册打包命令为 `python3 -m targeted_observing.review_tools.publish_atlas PNG目录 导出目录 图册目录`；PNG需要先由 `plot_exported_plan --all-nights --workers 3 --iers ...` 生成。

## 未被这些模拟证明的内容

百颗目标保持不变，所有本轮标准情景仍为0/200次达到101；这不是总体概率严格为零的证明。质量光变、唯一周期、物理模型验证分别报告。16个科学候选不是16个已完成项目；`science_case_complete`／`period_solution_complete`均未知。尚缺真实r精度、跨夜校准、饱和/CCD污染与控制器输入；事件任务还缺本观测季有效事件预报、形状模型相位/极轴正演及真实测光。这些不阻止软件对比，但不能用采样代理代替科学结论。
