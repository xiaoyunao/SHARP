# 晨昏搜寻 + 木星特洛伊监测（本地研发模式）

本目录是独立 Python 包，不导入/改写旧 `survey` 调度，不包含 SSH、cron、生产发布或望远镜控制功能。所有计划标记 `LOCAL_SIMULATION_NOT_FOR_TELESCOPE`。

详细工作目标见 [docs/GOAL.md](docs/GOAL.md)，来源和样本解释见 [docs/SCIENCE.md](docs/SCIENCE.md)，2026-11-03实际本地结果见 [docs/RESULTS.md](docs/RESULTS.md)。原对话快照在 `docs/conversation_20261103.md`；四份原交接附件现已原样归档到 [docs/specification](docs/specification/manifest.json)。逐项差异与后续工作见 [附件对照审计](docs/SPEC_ALIGNMENT_20261103.md)。`defaults.json` 是当前程序配置，**不等于原 YAML 已完整接入**。

## 本地运行

从仓库根目录运行。推荐使用已有包含 Astropy 的 Python 或独立虚拟环境。

```bash
python3 -m pip install -r targeted_observing/requirements.txt
python3 -m targeted_observing.fetch_sources
python3 -m targeted_observing.review_data
python3 -m targeted_observing.plan --date 2026-11-03 --campaign targeted_observing/data/campaign_20261103.json --output targeted_observing/products/20261103
python3 -m targeted_observing.plan --date 2026-11-04 --campaign targeted_observing/data/campaign_20261103.json --previous targeted_observing/products/20261103/plan.json --shift-minutes 45 --output targeted_observing/products/20261104
python3 -m targeted_observing.plan --date 2026-11-06 --campaign targeted_observing/data/campaign_20261103.json --previous targeted_observing/products/20261104/plan.json --shift-minutes -45 --output targeted_observing/products/20261106
python3 -m unittest discover -s targeted_observing/tests -v
```

数据来源只在 `fetch_sources` 和重选场时联网；`plan` 使用给定冻结 campaign 与本地 IERS，不自动联网。缓存文件按URL和SHA256记录，重新获取更新资料请指定新的 `--cache` 目录。当前默认依赖是本机已有科学库；requests/astroquery用于公共数据获取，scipy用于太阳高度交点求根，matplotlib用于科学图。

按新日期重选14天 campaign：

```bash
python3 -m targeted_observing.design_campaign --start 2026-11-03 --output targeted_observing/products/campaign_20261103
python3 -m targeted_observing.campaign_visibility
python3 -m targeted_observing.seasonal
python3 -m targeted_observing.contamination
```

`contamination --campaign PATH --plans PLAN.json ...` 现在在计划访问时刻合并Hipparcos亮端补表和Gaia区域缓存，必须提供计划访问；不以白天小时星历极小值代替夜间污染。`seasonal` 的人口计数使用二体/H–G近似，只是季节方向检查；未来季节必须另取 Horizons 设计新场，不能直接复用11月冻结场。

## 默认行为和输出

- 晚晨分别60分钟，两个半小时槽；r30，4场×5访。太阳距角45–75°优先，75–90°扩展单独标记；无完整组释放时间。
- Trojan预算可设为0以停用；启用时仅一个连续最多120分钟块，r120、两场；扣除晨昏、−18°太阳界限和外部预约，按可见性降到90分钟两场/60分钟主场，少于60分钟跳过。
- 默认5分钟交接、每次15秒读出。转场取方案20/60秒与仓库实测正常分支模型加2秒余量的较大值；不是所有相邻场都强套20秒。
- 中心≥32.5°；采用完整足迹的球面包围角帽检验>30°，并保留时间采样间隔的保守余量。月距按30/40°工程初值；固定场逐曝光检查；目标逐访问记录在场状态，campaign几何留场门槛≥80%。
- 晨昏闭合巡回（包括末场回首场）每跳≤3.5°，阈值为附件暂定分类，不是实测角速度定律。优先保留4个旧场，失败尝试3个锚点加1新场；最终保留数写manifest。目标恢复需新星历，重拍原场不能当恢复成功。
- `plan.json`保存配置/EOP来源、时间和决策；`manifest.json`包含分块not-before/deadline、字段版本及UTC/北京时间未占用区间（扣整块及预约，含块内余量仍保留）；无外部预约时Trojan标记暂定分配；`EV_A.csv`等提供逐曝光时间；`fields.csv`汇总cadence；`validation.json`是独立5秒审计；`plan.png/pdf`展示时间轴与足迹。
- `data/campaign_20261103.json` 是版本化小型科学输入：两个固定场、5个核心的站心星历与文献证据；`target_audit_20261103.csv`保留60候选审查。缓存与运行输出不提交Git。

`defaults.json` 是可读JSON配置，不需要YAML依赖。v0只支持3°视场、晨昏4场5访/1小时；不支持的结构变化明确拒绝，不能只改配置就假定扩展模式已实现。曝光/开销/预约/月距等可调整；前三个有效夜建议保持r及30/120秒一致。

预约格式（UTC，追加到配置 `reservations`）：

```json
[{"start_utc":"2026-11-03T15:00:00Z","end_utc":"2026-11-03T16:30:00Z","program":"external"}]
```

## 执行契约

本地不生成可执行整夜TXT。旧控制器TXT没有日期/等待命令；仅有JSON时间戳不代表控制器会等待。`controller_filter_tokens` 默认空；`export_controller_drafts`缺合法映射与证据时拒绝，不能回退W。即使提供映射，也只输出带警示的分块草稿。

未来现场适配应按每次曝光实际到位时刻调用 `execution.gate`：早到等待；晚到对剩余完整组重算；超出截止/可见性则停止。该函数只返回决定，不发送命令。控制器需自行提供真实到位、读出完成、滤镜状态和序号；天气、地形、机械禁区及本地控制器异常不由天球几何证明安全。

新模式不调用旧无界follow-up后插，不读取旧三夜冷却历史，也不会把计划访问写成实际完成历史。恢复目标可在之后实现为整个块替换，必须走同一最终校验；本轮未实现真实NEOCP自动队列。

## 三有效夜复盘

填写 `data/feedback_template.csv` 的副本，记录实际曝光中点、有效标记、背景、FWHM、零点、深度、按亮度稳定恒星重复误差、跨夜零点、读出/转场及已知移动源恢复数。`valid`必须为0/1；计划行不是实际观测。

```bash
python3 -m targeted_observing.feedback --input my_actual_exposures.csv --output targeted_observing/products/feedback.json
```

首次复盘一次只调整主要瓶颈：场数、曝光、亮度样本或选场。r波段深度未知，不从白光20.8等推导实测性能；慢转测光不能任意减去每夜均值。生产短轨迹范围仍3–63角秒/小时，不能视为已验证快速NEO检索。

## 已知限制

- 实际r控制token、定时执行接口、CCD方向/坏区/亮星光晕待现场确认；没有服务器部署。
- 先筛60个高排名对象，场选择为启发式；不是全天全局最优或科学无偏样本。
- 场重定心必须新版本；14天留场不等于14晚均晴或一定可解周期。
- TESS 2025完整重分析表已接入，保留区间、备选周期与来源；18054未被该表收录不代表未知周期。
- IERS仅用于离线规划，预测EOP允许但记录来源/覆盖；超表范围拒绝。长期预排不等于临近执行时EOP仍新鲜。

实际反馈汇总增加 `complete_five_visits` 和 `usable_short_tracklet`：后者要求至少4次有效曝光且中点跨度≥720秒，包含满足条件的5次序列；`sequence_class`互斥分类。这里只统计采样条件，不证明某颗小行星已被检出或成功链接。Trojan行不套晨昏短弧分类。

## 最新审查落实

见 [GPT Pro审查修改与复现](docs/PRO_REVIEW_CHANGES_20261103.md)。晨昏采用A/B联合选择与三夜前瞻/复用，单场退化按预算取点，Trojan统一按计划访问≥80%留场。`design_campaign`现输出三组比较，完整结果位于输出目录各候选子目录；不自动覆盖旧campaign。当前角色分类、候选池和亮星缓存的限制均在报告中说明。

## 一年期模式（本地模拟已完成）

当前目标以 [年度GOAL](docs/GOAL.md) 和 [最新需求正文](docs/annual_requirements_pro_20261104.md) 为准。旧60对象/三夜方案只作为回归基线。`annual_defaults.json`固定日期、候选宽松门槛、2/3/4/6场、驻场张数、200天气实现与三种假定r性能；不是现场测量值。

另有 `precision_requirement` 灵敏度探索情景，不能称为已实现的仪器性能。年度主队列将文献权重以25%强度纳入完成概率评分，保留未核查状态；达到点数/夜数、但等待21天基线的对象先释放致密采样，届时回来补齐。评分不读取未来天气。`--science-strength 0`用于纯吞吐量对照，`--duration-s 10800/14400`仅作额外时间申请实验，不能由标准两小时模式执行。

年度去重面积使用HEALPix球面像素，需要 `pip install -r targeted_observing/requirements-annual.txt`；这一额外依赖仅用于面积评估，不改变原调度器依赖。

```bash
python3 -m targeted_observing.annual_sources --workers 2
python3 -m targeted_observing.annual_twilight --workers 4
python3 -m targeted_observing.annual_visibility --workers 4
python3 -m targeted_observing.annual_simulation --workers 4 --save-curves --scenario baseline optimistic conservative
python3 -m targeted_observing.annual_simulation --workers 4 --nominal --scenario baseline
python3 -m targeted_observing.annual_twilight_queue --visits 5 --workers 3
python3 -m targeted_observing.annual_twilight_queue --visits 10 --workers 3
python3 -m targeted_observing.annual_coverage --visits 5 --causal --limiting-mag 20
python3 -m targeted_observing.annual_coverage --visits 10 --causal --limiting-mag 20
python3 -m targeted_observing.annual_delivery
python3 -m targeted_observing.annual_effective_area
python3 -m unittest discover -s targeted_observing/tests -v
```

数据/计划保存在忽略提交的 `cache/annual` 和 `products/annual`；下载可续跑。`--days 3 --limit 60`仅用于Trojan/visibility冒烟，不能用于年度结论。晨昏全年输出包含两个策略各730窗口及零计划理由。星历必须完整；失败表不允许静默丢掉对象。新IERS `iers_20261104.all` 覆盖到2027-11-05，仍是预测规划数据。

名义计划使用计划成功反馈，不等于真实完成；`AnchorState.feedback`只应在实际有效短弧后更新三夜期限。年度天气模拟使用自己的逐对象累计状态；`annual_visibility`提供带保守时间余量的可选时窗，避免按未来天气设计队列。完成后的光变保存合格前缀，后续低质量伴随帧不会抹掉已有证书。

完整结果和其余周期/三深度/额外预算复现命令见 [年度材料索引](docs/annual_reports/README.md)、[Trojan报告](docs/annual_reports/TROJAN_ANNUAL_REPORT.md)、[晨昏报告](docs/annual_reports/TWILIGHT_ANNUAL_REPORT.md)。已生成365夜主方案、37078条曝光并独立验收；比较24组标准策略、各200天气实现。四场驻场两张的基准合格光变P10/P50/P90=36/42/47，当前两小时模式不能承诺百颗；周期与振幅另计。四小时/更高精度仅为额外申请实验，不改变正常预算。

`annual_trojan`是先前名义评分的参考实现；最终主方案由`annual_simulation --nominal`和`annual_delivery`生成。原三夜`defaults.json`及冻结campaign保留，年度推荐不静默替换它们。GitHub同步源码和审查摘要；大型年度产品、原始科学缓存留在本地。精确复现数字需要保留当前输入快照；日后重新下载的目录和IERS可能不同。

## 2026-11-05 续接：365夜图册与优化比较

见[本轮结果及复现入口](docs/continuation_20261105/README.md)、[逐日图索引](docs/continuation_20261105/atlas_B/index.html)。原a5cb945名义年四份轻量数据与378张图已交付；新增去重评分/独立完成预测、偏置中心、真实路线重评分、10%/20%预算内特殊科学状态，以及晨间十访后对齐。

当前 `annual_defaults.json` 的policy为4；`--continuation A`保留旧policy3，B/C/D显式启用新比较。基准测光200次年度中位分别42/47/44/39，质量达标不等于唯一周期或模型已验证；>100颗目标仍未达到。`--days`/`--realizations`的小试验使用独立目录后缀，`--nominal --continuation B`使用`nominal_`前缀，不覆盖天气实现。

科学任务使用文献候选周期和独立采样代理；普通质量达标后可继续长基线。完整物理解和实际r/饱和/污染仍待验证。图册PDF打包需要Pillow和reportlab，仅用于离线报告；现场与生产依赖未改变。原三夜campaign与生产路径继续保留。

2026-11-05 阶段收尾采用 B 名义年作为默认图册，包含 365 夜和 37,264 条曝光。图件与补充采样统计见[收尾说明](docs/continuation_20261105/FINAL_REVIEW.md)。
