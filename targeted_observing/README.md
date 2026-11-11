# 晨昏搜寻 + 木星特洛伊监测（本地研发模式）


## 当前晨昏入口：MP 离线固定列表（2026-11-11）

全年365夜已生成并校验，结果与局限见 [本轮结果](docs/offline_twilight_20261111/RESULTS.md)。

本入口适用于控制软件与服务器不联通的情况。每天提前生成两份八列观测列表；**没有实时反馈、没有列表内WAIT指令**。控制滤镜token由用户确认为 `R`，目标名沿用 `MP_0001` 格式。仅本地开发测试，不部署、不改cron。

```bash
python3 -m targeted_observing.offline_twilight --date 2026-11-11 \
  --output targeted_observing/products/mp_nights/2026-11-11
```

输出 `EV.txt`、`AM.txt`、`START_TIMES.json`、逐曝光时间CSV和完整 `plan.json`。**按START_TIMES中的北京时间启动各列表**，清晨属于次日。程序在单侧1小时窗口内比较开始时刻和晚15分钟两个预定启动时刻；启动时刻提前固定，不由现场反馈决定。空列表表示本次有限候选搜索未找到可行场组。

只使用生产原表1855个MP场，原编号、中心及角点保存在 `data/mp_footprints.json`，含原FITS来源hash。晨昏足迹高度25度（中心含边缘余量），月距沿用30/40度相位规则；从太阳−10度边界规划。先争取5访，不可行再降至4或3访。场数4–12，实际连续扫描的同场间隔要求5–12分钟，8分钟是软目标。近太阳优先，路线优化比较近赤纬和紧凑场组、RA顺序及闭合路线；不承诺全局数学最优。

曝光用Patat近天顶R背景先验、显式假设的低空/日距角/月光修正，按较亮的时序包络预计算5–30秒（5秒档位），−18度以下正常30秒。基准为太阳−10度、天区高度25度、日距角40度、月落时5秒；这是试观测初值，**未测得本机r的ADU/线性上限或极限星等**。使用 `offline_defaults.json` 调整初值，试观测数夜后重生成后续列表。

slew采用 `9.11+max(11.95*sqrt(dRA),12.30*sqrt(dDec))`，RA轴角不乘cos(dec)，跨零处理；读出预算15秒另加。首转场名义90秒、首次换片10秒。独立5秒网格检验还覆盖首转场更早及逐帧提前误差、合计300秒晚到包络；晚到和累计开销共同消耗这一预算，**不是允许任意延迟后继续执行的保证**。静态TXT不能自动截止，超出允许启动时间或严重延迟时须由观测端停止/换用重生成列表。天气回放会报告列表继续执行造成的不合格帧，不虚构实时中止或补拍。

全年固定列表及200条天气回放（汇总依赖沿用 `requirements-annual.txt` 中的 pandas/healpy）：

```bash
python3 -m targeted_observing.offline_annual --workers 8 --resume \
  --output targeted_observing/products/mp_offline_20261111
python3 -m targeted_observing.offline_delivery \
  --input targeted_observing/products/mp_offline_20261111 \
  --report targeted_observing/docs/offline_twilight_20261111 \
  --trojan-source targeted_observing/products/annual/selected_nominal_B_20261105
```

默认模拟2026-11-01至2027-10-31。`--resume`仅复用配置、输入grid/IERS哈希一致且验收通过的日计划；`--reuse`用于已有完整全年计划的回放。IERS沿用 `cache/iers_20261104.all`，缺缓存时按下方原年度说明获取覆盖全年的表。Trojan兼容审计可不传；该参数只核查旧B项目时间块仍在共享夜间内，不重跑或改变其科学策略。当前改的是晨昏，旧TW年度结果独立保留。

本地 masked catalogs 的新模式关联入口：

```bash
python3 -m targeted_observing.offline_processing --night 20261111 \
  --root /absolute/local/data --output /absolute/local/tracklets --workers 4
```

输入为 `<root>/<night>/mask_gaia/OBJ_MP_*_cat.fits.gz`。新入口从offline配置读取3–750″/小时，保留原下限与静态源/星等差筛选；原 `heliolincrr/run_single_night.sh` 默认仍为3–63，增加 `TRACKLET_VMIN/TRACKLET_VMAX` 环境覆盖接口。实际曝光中点优先EXPSTA/EXPEND，回退到开始＋EXPOSURE/EXPTIME的一半，支持primary HDU时间与extension WCS，不再固定加15秒。缺失/矛盾时间拒绝处理而非暗中假定30秒。

验证命令：`python3 -m unittest discover -s targeted_observing/tests -q`。近太阳速度注入仅验证配对范围及时间处理；全年条件检出模拟使用真值身份，**不是实际链接效率、定轨成功率或NEO发现数**。旧 `feedback.py` 与TW报告属于历史原型，不是新MP模式的第二套验收口径。

试观测后需保存每帧实际曝光起止、读出/转场、滤镜、源屏蔽后的分块背景及梯度、饱和比例、零点/极限星等和移动目标检出。离线修改理论曝光参数后重新生成下一晚列表；不要求现场实时传回。

## 历史 v0 入口与说明

以下直到后续年度审查章节为原 TW/v0 入口文档，仅供历史复现；其中固定30秒、−12度、实时执行闸门及4点/12分钟口径均不适用于上方当前MP入口。新MP计划标记 `OFFLINE_FIXED_LIST_LOCAL_TEST`，原v0标记 `LOCAL_SIMULATION_NOT_FOR_TELESCOPE`。两者均不连接服务器或直接控制望远镜。

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
