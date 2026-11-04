# 年度材料索引

- [Trojan 年度报告](TROJAN_ANNUAL_REPORT.md)
- [晨昏年度报告](TWILIGHT_ANNUAL_REPORT.md)
- [机器可读摘要](results.json)；其余同目录CSV和PNG由 `annual_report.py` 从完整运行结果生成。
- [年度目标](../GOAL.md)、[完整需求正文](../annual_requirements_pro_20261104.md)。

## 本地完整产品

这些是可重建运行产物，不提交Git。GitHub同步源码、配置、科学输入和本目录审查摘要，避免把大型缓存/日志加入仓库。

| 路径（相对 targeted_observing/products/annual） | 内容 |
|---|---|
| candidate_audit.csv / candidates.json | 全量16445行接受/拒绝理由、1063全年候选 |
| ephemeris_index.json / input_verification.json | 下载状态、原始星历hash、IERS与12历元抽验 |
| visibility/ | 365夜完整场菜单、版本、目标、共同可见窗口 |
| selected_nominal/ | 主策略365夜JSON、exposures.csv、time_ledger.csv、fields.json、retention.csv、独立verification.json |
| nominal_queue_*_science25_v3/ | 八种标准策略的全晴名义逐曝光表和365夜零计划账本 |
| queue_*_science25_v3/ | 各200天气实现的逐夜账本、逐对象产出；qualified_curves.npz保留误差/零点/光行时采样 |
| queue_4fields_dwell2_*_science25_v3/period_yield.csv | 全200天气的周期/振幅注入真值与结果，盲搜/先验分开 |
| twilight_5/ 与 twilight_10/ | 各365夜计划、历史备选菜单、causal_choices、146000行窗口账本、面积与移动源注入 |

`*_smoke`、旧的无science后缀队列、science v1/v2、fixed40诊断，以及无r深度后缀的旧moving_injections不属于本次最终结果。额外budget14400是申请敏感性，不是标准两小时执行计划。

## 复现顺序

在仓库根目录运行；需要Python、numpy/scipy/pandas/astropy/astroquery/requests/matplotlib/healpy。先运行基础fetch_sources/review_data取得母表与文献；年度计算另需annual_defaults.json所指定的IERS文件。当前本地已冻结iers_20261104.all；如在新机器获取，可从 https://maia.usno.navy.mil/ser7/finals2000A.all 下载并保存为新的带日期文件，确认MJD末端覆盖stop_date_exclusive再加4日，并更新配置，勿覆盖既有快照。后续用固定缓存和hash复跑；provenance.json记录版本与指纹。具体来源下载入口见 targeted_observing/README.md。默认日期、35日campaign、7日版本、200实现及随机种子全部在annual_defaults.json。

```bash
python3 -m targeted_observing.annual_sources --workers 2
python3 -m targeted_observing.annual_ephemeris_qa
python3 -m targeted_observing.annual_twilight --workers 4
python3 -m targeted_observing.annual_visibility --workers 4
python3 -m targeted_observing.annual_simulation --workers 4 --save-curves --scenario baseline optimistic conservative
python3 -m targeted_observing.annual_simulation --workers 4 --nominal --scenario baseline
python3 -m targeted_observing.annual_twilight_queue --visits 5 --workers 3
python3 -m targeted_observing.annual_twilight_queue --visits 10 --workers 3
python3 -m targeted_observing.annual_delivery
python3 -m targeted_observing.annual_effective_area
```

对两种visits分别运行annual_coverage的 `--causal --limiting-mag 19/20/21`；面积另用 `--nside 512 --skip-injections` 检查收敛。对八种baseline队列运行annual_periods，主四场驻场两张采用 `--objects 24 --phases 2`，其余 `--objects 8 --phases 1`。对主策略baseline与optimistic运行annual_period_yield。额外精度/时间实验分别用annual_simulation的 `--scenario precision_requirement --fields 4 --dwells 2` 和 `--duration-s 14400 --fields 8 --dwells 2 --scenario baseline optimistic precision_requirement`。全部完成后：

```bash
python3 -m targeted_observing.annual_report
python3 -m unittest discover -s targeted_observing/tests -v
```

报告生成器会检查24组365×200账本、时间闭合、年度去重质量数与累计完成单调性，缺最终结果即失败，不以早期诊断填补。两份报告包含模型假设、未建模项及未达百颗结论。没有任何远程发布命令。
