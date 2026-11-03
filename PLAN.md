# PLAN

## Current objective: 最新GPT Pro审查修改（2026-11-03）

- 继续仅本地研发、GitHub专用分支同步，不连接/部署服务器；原生产/MPC hold不动。
- 最新审查已通过read_thread取得并归档；官方镜像已取得TESS2025全表及Hipparcos亮端补表。9912行/52TJN/7候选命中复现，原5核心未被本表收录。
- 已修改：TESS完整区间与备选解、DAMIT一致比较、保守科学角色；参宿四/18054按计划访问风险；晨昏A/B联合候选和三夜前瞻/复用；原丢锚点案例恢复7/8；单场75/89分钟30/36点；统一计划访问80%留场及逐帧状态。
- 候选：原基线与两组替代方案，备选含22180/2223（另组加4792），远转场后17点/390秒；部分对象85.7%留场、覆盖12夜。完整报告见targeted_observing/docs/PRO_REVIEW_CHANGES_20261103.md；不自动覆盖旧冻结场。
- 验收：编译、21项测试、3组×3夜计划和3组×14夜Trojan独立审计；运行记录在products/pro_review，验收结果以verification.json为准。
- 待现场：r120精度—星等、亮星光晕/坏区/饱和、真实转场稳定时间、有效点率、合法滤镜token。尚无真实图像，不能宣称候选场测光合格。
- 仍待本地扩展：2–3场带等待回退、恢复队列、去重覆盖面积与数字ID/完整YAML兼容；当前保持禁用，不作为本轮已实现项。
- 下方历史生产记忆保持原日期口径，本轮不巡检服务器。

## Latest run check: natural cron completed (2026-10-28)

- 今日09:00:01自然cron启动、09:31:57完成，两段退出码均0；昨日5个部署文件hash一致，IERS新缓存有效，survey无IERS警告。
- 今晚504条计划已于09:00:30发布，runtime/日期TXT/current的SHA一致；观测电脑拉取和逐曝光高度策略未在本次检查。
- 昨晚10/27原始474/L2 474，全部474帧匹配完成，正式匹配3926条、mask15为4242条；41初选拒19（形态18+PSF1，46.3%），pass5/review17，22条人工候选/66次测量/22 GIF。评分、FITS/manifest、GIF解码、tar完整性通过。
- 无新增运行阻塞，19个旧MPC/hold保护文件不变、无新上报，hold保持。证据在`outputs/runtime_check_20261028/remote_report.json`；本次只读巡检完成，后续继续每日任务/人工复核。旧的5份上游缺失、逐曝光高度问题仍需独立处理。

## Latest readout: nightly screening counts (2026-10-27)

- 10/11串联上线后至10/26，8个有报告夜：1036初选、595拒绝（57.4%），其中形态556、PSF追加39；441可人工复核=pass91+review350，平均每个有报告夜过滤74.4条。
- 逐夜数、理由、raw SHA与人审队列已核对，来源/采集脚本在`outputs/screening_counts_20261027/`；无观测夜不计平均，不把历史回放当上线前实际表现，不用过滤数推断真实假源率。
- 本次查询完成，无生产变化；后续运行及待跟进问题沿用下节，MPC继续暂停。

## Current status: nightly recovery complete (2026-10-27)

- 已完成：本机下载USNO新IERS表（预测起点10/23）并备份后导入服务器；今晚501条计划已于16:14:39发布，三处TXT hash一致；survey失败与数据处理解耦已部署，本机/服务器各6项回归通过，MPC/hold保护文件未变。
- 数据恢复完成：16:44:29，396/396帧known无错误/警告，matched2748条、mask15为3036条。129条初选经串联筛选pass4/reject38/review87，91条人工候选/274次测量/91 GIF；评分、CSV/FITS、逐GIF解码、tar完整性验收通过。
- 执行与证据：Slurm array255223/finalize255224完成；串行绘图受控中止后在独立目录用原绘图器8进程恢复，129 GIF生成、10份串行已有文件hash一致。最终processing-only入口复用已有包成功，退出码0；19个既有MPC/hold保护文件不变、无新上报。凭据在`outputs/recovery_20261027/final_validation.json`，服务器对应`data/heliolincrr/recovery_20261027/`。
- 额外边界：独立QA发现既有调度50%窗口可见规则会留下14次高度<30deg曝光（最低24.234deg），月距全部通过；本轮保持原策略，逐行证据在`plan_validation.json`，后续可单独加强逐曝光可见性。尚未验证观测电脑是否已拉取新表。
- 额外输入边界：401原始帧中5帧缺上游L1/L2，本次覆盖全部396个可用L2；缺失文件清单在`input_completeness.json`，上游原因未在本次修改。
- 下一步：本次恢复已完成，继续保持MPC暂停；可按后续用户要求独立诊断5份上游产物缺失、加强逐曝光高度约束或核实观测电脑拉取。IERS手动离线更新及单独补跑命令已写README。

## Initial investigation (2026-10-27; superseded by recovery above)

- 巡检已完成：用户确认“计划未更新、昨晚无处理结果”，两者同源。10/25–27 IERS预测超过30天且主/备用下载失败，survey exit1导致总入口在known/unknown前退出；cron本身正常触发。
- 当前产物：10/25–27均无新计划；current_plan仍为10/24的495行。10/26有401个MP raw、396个MP L2，但L4空、known/unknown未启动。10/18–21已有known complete和95/73/24/53条复核包，缺GIF均0；22–25无MP输入。
- 当前范围：本次只读调查，无生产修复/缓存更新/补跑/发布。证据在 `outputs/runtime_check_20261027/remote_report.json` 及 `evidence/`。MPC hold持续有效。
- 后续里程碑：1)恢复可靠且预测起点足够新的IERS缓存；2)隔离生成、校验并发布当晚计划；3)保持hold补10/26 known/unknown；4)解耦survey失败与数据处理，并检查current_plan观测日期，防止静默复用旧表。
- 验收：IERS转换不再报错且使用新鲜预测；当晚JSON/TXT/runtime/发布目录一致；10/26逐帧known完成并生成通过筛选的复核包；任何计划故障仍明确报警且不阻断独立的数据处理；正式上报仍暂停。
- 注意：当前分支保留维护前astro_utils，生产基线应参照 `/Users/yunaoxiao/.codex/worktrees/nightly-maintenance/smt_asteroid` 的 `2bce0a18`，勿以旧分支代码覆盖服务器。尚未验证观测电脑实际拉取的本地副本。

## Latest completed task: Small-RA threshold expanded to 2 degrees (2026-10-22)

- 当前范围：按用户最新补充，小RA取非零<2deg；直接样本还要求Dec<2deg、两端30s及原有邻接/指向/唯一性条件，正时间无上限。
- 直接证据：70对、13夜、6种角度，较旧<1deg新增56对；实际Dec<0.81deg，各角度中位耗时随RA0.396..1.851deg从17.35增加到23.50s，整体16.11..26.03s，未出现明显小角度慢分支。
- 分层对照：RA<2 vs2..3deg：26夜483/300对，+0.029s、95%[-0.115,+0.182]；1..2 vs2..3deg：18夜241/211对，-0.018s、95%[-0.271,+0.187]。按整夜bootstrap；不能以总时间无增加排除被Dec遮蔽的微调。
- 交付：`outputs/slew_calibration_20261022/expanded/ra_under2/` 下的 `small_ra_evidence.png/.pdf`、`isolated_small_ra_pairs.csv`、`isolated_small_ra_summary.csv`、`small_ra_results.json`、`validation.json`；原结果完整保留，一周双子图路径不变。
- 验收：py_compile、端到端重跑、70对原始header算术及邻接检查、旧14对/原匹配队列/正常参数不变、新效应独立重算、图目视通过。
- 剩余边界：仅6种几何，<1deg仍14对，不能测出任意非单调小角度曲线或确认往复修正；后续成因核实依赖控制/编码器日志。
- 下一步：本次2deg范围更新已完成，按用户后续反馈继续；生产调度参数保持原状。

## Latest completed task: Two-axis figures and annual small-RA audit (2026-10-22)

- 用户要求：拟合图只保留RA/Dec两子图，并扩大范围直到有足够小RA证据；已完成双子图，并穷尽20251101..20261021可用巡天文件。微调机制的可辨识性仍受实际观测几何限制，不能宣称已测出完整小角度响应。
- 交付：`week/slew_axis_fits.png`、`.pdf` 为所要的一周双子图；`expanded/slew_axis_fits.*` 为近期30夜同款图；`expanded/small_ra_evidence.png` 和 `small_ra_results.json` 为全年调查结果，均在 `outputs/slew_calibration_20261022/`。
- 数据：184夜60503文件索引，MP为164夜51032份，51022可读/10不可读，50210对连续MP；273份跨目录重复header去重（50749唯一可读曝光键），不跨别名重接；5份非标准指向的相邻对从主比较隔离。现场1855天区表经全年header核对，不采用仓库旧天区表作此次几何结论。
- 小角度直接证据：仅14对小RA且小Dec正式30s曝光，跨5夜、RA0.396/0.763deg两组；耗时中位17.352/18.121s，无分钟级慢分支；不足以拟合RA<0.1deg的独立响应。
- 条件比较：同夜、0.5deg Dec分层，RA<1 vs1..3deg共92夜3882/2542对，时间差+0.0218s、整夜bootstrap95%[-0.0613,+0.0985]；RA<0.1对照49夜661/1559对，+0.0496s、区间[-0.0602,+0.1772]。额外>5s等待率未检出上升，不能据此识别或排除被Dec遮蔽的微调。
- 近期正常公式：30夜9835对T=9.24145+max(11.92277sqrtRA,12.22922sqrtDec)，留整夜中位误差0.497s/P90=1.259s，未改生产参数。
- 验收：全文header只读收集/完整性、20份独立读取对照、合成匹配效应与零事件不确定性边界、原始算术与全年统计、双子图与证据图目视通过。恢复入口和解释边界见survey/README。
- 剩余限制：header只有曝光端点，无法显示反复修正路径；其他程序2次小RA长等待与自动巡天分开。同场零位移也有长等待，不能直接归因于RA。现有归档缺少更多独立极小RA/小Dec几何，继续堆积同类样本不能解除遮蔽。
- 下一步：本次图与调查已交付；如继续追查物理原因，优先核对已有控制/到位/编码器日志，不以没有检出普遍效应声称所有小角度均无问题。

## Latest completed slew task: Weekly validation (2026-10-22)

- 目标已完成：按用户要求扩展20261015..20261021一周；10/16、17未找到OBJ，实际5夜2050份header、2037对有效科学转场。
- 当前正常转场公式：`T=9.1148+max(11.9542*sqrt(dRA),12.3044*sqrt(dDec))`，角度以deg、时间以s计；T已减曝光及14.47s读出。相对单夜系数稳定；不改生产调度参数。
- 模式核对：117对RA、1920对Dec主导（94.26%）；与按RA排序导致小RA、大Dec位移的巡天几何一致，主导标签仍有模型条件。
- 验收：10种模型比较；留整夜绝对误差中位0.530s、P90=1.293s，逐夜中位0.501..0.551s；500次整夜bootstrap及幅度/vmax剖面、坐标/时间/长间隔敏感性；完整端到端运行、合成参数恢复、独立重算预测和4张图目视通过。
- 小RA边界：1682对RA<1deg的Dec均>=2.496deg；最小RA0.00275deg，但两轴都小的对照仍缺失。非单调微调项没有实质跨夜预测提升，不能把拟合幅度当机械微调耗时。
- 偶发延迟：14个有效转场|residual|>5s均保留；10/20两组相同几何平时约38.4/43s、偶发135.1/161.1s，重复对照已保存。另两次>180s科学间隔单列，纳入后正常公式几乎不变。
- 机械峰速：RA>30deg仅5对，vmax对这几对极敏感；不报告优化器返回值为实测机械最高速度。正常模型不能解释全部长尾。
- 交付：`outputs/slew_calibration_20261022/week/{results.json,weekly_validation.png,slew_fit.png,slew_diagnostics.png,small_ra_diagnostics.png,small_ra_anomaly_repeats.csv,validation.json}`；7份header JSON保留在父目录。
- 下一步：本次周分析已完成；后续需要物理成因时，以已定位异常时间对照控制/编码器/到位日志；需要机械vmax或微调标定时优先补足几何与控制模式覆盖。

## Previous caveat: RA small-angle adjustment may be slow (2026-10-22)

- 用户补充：RA小角度微调可能反复修正/卡住；保留为控制模式假设，不能将旧测试图的长耗时分支直接当异常。
- 已核验：本夜393对RA<1deg转场全部同时移动Dec至少2.573deg，没有RA<1deg且Dec<2deg的对照；小RA时间可被max中的Dec分量遮蔽。旧模型的主导轴标签是单调假设下的结果。
- 模型边界：下面平方根模型仅为本夜正常转场的有效描述，不能据此排除非单调小RA分支，也不能用其小角度外推声称RA微调速度已测出。
- 后续方法：允许独立的小角度微调延迟/超时分支；若继续扩大到已允许的一周，应优先寻找两轴位移都小的样本、控制Dec位移和方向。只有端点header无法验证反复修正的具体机制。
- 产物：`one_night/small_ra_bins.csv`, `small_ra_pairs.csv`, `small_ra_identifiability.json`；额外延迟数字仅是遮蔽的构造示例，不是实测卡顿时间。

## Latest slew result: One-night pilot completed (2026-10-22)

- 用户范围：先做一夜，不够再做一周；本次20261021的有效耗时/平均速度模型已完成，未扩周。
- 数据：463份完整OBJ主header，462标准科学帧；461对连续科学曝光、110种有向天区转场；序号连续、无header读取失败。
- 公式：`T=9.102+max(12.086*sqrt(dRA),12.369*sqrt(dDec))`；角度以deg、时间以s计；T为同一OBS_DATE差减前帧EXPOSURE与14.47s平均读出。RA小时乘15、跨零、不乘cos(dec)。
- 验收：7种模型对照；天区转场分组/连续时间块各5折，平方根模型绝对误差中位0.523/0.532s、P90=1.400/1.373s；300次块bootstrap；时间字段/实际指向/读出均值敏感性；vmax剖面；合成边界检查和最终图目视通过。
- 产物：`outputs/slew_calibration_20261022/headers_20261021.json` 和 `one_night/{results.json,fitted_pairs.csv,slew_fit.png,slew_diagnostics.png}`；复现入口见 `survey/README.md`。
- 可用范围：RA主导27对、Dec主导434对；各自主导角范围RA2.27..62.05deg（长端很稀疏）、Dec2.55..20.00deg。有效运动平均速度为`sqrt(angle)/k`，含常数开销的转场速度应另算`angle/(t0+k*sqrt(angle))`。
- 未解决：vmax高速端无明确上限约束，不能把梯形模型最优值当实测峰速；常数项不是独立测得稳定时间。4次异常转场保留在所有拟合与验证中，长尾导致全样本分组CV RMSE5.284s；尚无跨夜验证或小于约2deg主导样本。
- 下一步：本次单夜交付完成；后续如要求跨夜泛化/机械最高速度，按已允许范围扩最近一周，重点增加大位移与恒速分支证据、核查异常时序。当前生产readout15s/slew35s未改。

## Latest completed task: Selected figures and GOTTA light-curve addition (2026-10-19)

- 当前交付：`outputs/exhibition_20261019/SHARP小行星巡天项目_中文介绍_20261019.docx`，12页/13图，继续保持纯中文正文加英文图。
- 已完成：删除上版图11/13/14并连续编号；图1换20260901天区；图10两个纵轴均0–1角秒；原图15（现图12）换为当前形态和PSF均pass、人工确认为真源的00001mX三帧。
- 光变内容：参考 `gotta_asteroid_1/revise1/main.tex` 补充GOTTA原型与兴隆施密特联合观测、周期分析方法、24目标先导样本与12/12分类；不自行修订该项目待合作者确认的具体周期和误差。
- 验收：只读核对现行筛选报告raw SHA及3帧双normal/人审is_real=1；13图内嵌hash与删除项检查通过；最终12页逐页目视通过。记录在 `expanded/qa/selected_final_validation.json`，16图旧稿已备份。
- 未解决事项：无本次文档阻塞；fetch因GitHub TLS失败未更新远端引用，本轮只提交本地记忆，不改生产或GOTTA项目。
- 下一步：以本次13图稿为当前版本，按用户后续意见继续修订。

## Previous draft: Expanded survey science with 16 manuscript-style figures (2026-10-19)

- 当时交付：14页普通中文段落与16幅英文图；按用户意见补回论文的科学统计，不做展板排版；已备份为 `outputs/exhibition_20261019/expanded/previous_sixteen_figure_document.docx`。
- 已完成：已知小行星天空/轨道类别/速度/太阳距角/距离/测光/名义匹配比例/天测残差/重复观测；未知候选数量/筛选/弧长与RMS/亮度速度/真实轨迹与三帧；新增月度观测、最新宽场和合成streak示意，保留多夜HelioLinC、机器学习、光变与相减研究。
- 当前数据：观测至20261018，49853帧/161夜/1444场/约10481平方度；known727133条/75967颗；初选9538条，当前保留84条/269测量/47夜；近期1873条筛选分为reject1021/pass303/review549。
- 风格要求：纯中文正文加图，图内英文，沿论文原字体/配色/图形参数；不添加来源、审计说明、空泛总结。新图使用当前数据，streak明确标为合成示意。
- 验收：论文12图与主要科学内容逐项对照完成；三人检查全部14页，无缺图/乱码/裁切；16张内嵌图、无表格/页眉页脚/手动分页；统计独立复核，修正太阳距角参考系和不同对象极值措辞。
- 未解决事项：无本次文档交付阻塞。提取/绘图脚本及验证在 `outputs/exhibition_20261019/expanded/`；运行图像产物不提交，生产问题和MPC暂停另由原主线处理。
- 下一步：以此扩充稿为准，按用户后续文字意见修订；数据变化时沿用原论文绘图方法重算。

## Previous draft: Plain Chinese project introduction with four figures (2026-10-19)

- 当时交付为4页普通段落加4图；取消初稿表格、来源说明、审计措辞和版式建议；现保留为 `outputs/exhibition_20261019/expanded/previous_four_figure_document.docx`。
- 用户风格要求：图内全部英文，字体/配色/图形参数沿用提供的20260830论文；后续修改继续遵守，不重新设计图表。
- 已完成：用截至20261018的全量known数据重画黄道天空分布、孔径测光和重访统计；最新Hapke三帧以原review风格重绘；补充15日HelioLinC多夜关联和streak机器学习研究。
- 验收：全量727133条/75967对象及HEALPix合计一致；字体/绘图代码与论文一致；DOCX四页逐页检查通过，无表格、页眉页脚或手动分页。
- 未解决事项：无本次文稿阻塞。运行快照、可复现脚本及检查文件均保留于`outputs/exhibition_20261019/revision/`，不提交运行产物。
- 下一步：按用户后续意见修改正文；数据变化时沿用论文原函数重画。生产问题和MPC暂停继续由原主线单独处理。

## Previous draft: Chinese exhibition document and fresh data snapshot (2026-10-19)

- 当前交付已完成：7页中文Word、6幅新图、可直接上展板的短文和来源说明；主文件在 `outputs/exhibition_20261019/SMT小行星巡天项目_中文展板内容与最新统计_20261019.docx`。
- 数据冻结：2026-10-19 15:34–15:38 CST只读提取；观测截至20261018，完整unknown筛选截至20261015。标准raw49853帧/161夜/1444场；known727133匹配行/75967目录身份/157夜；至少2夜43763个身份。
- 新版筛选范围：25个保存且完成报告夜，1873条=reject1021/pass303/review549；852为可人工队列，不是已核验/新发现数。10/18的255条数量门控拦截，单独列为未完成。
- 验收完成：known累计表与逐夜8列精确一致、25份筛选SHA一致、统计合计闭合；Word最终7页逐页看图通过，六图均新绘，Hapke运动图来自10/18 L1实际像素。
- 边界：原始标准曝光与含后随的匹配表口径不同；不报告确认发现、完备性或净覆盖面积。原有MPC暂停仍有效，未改生产代码/数据/人工决定。
- 尚未指定：最终展板尺寸、单位名称、团队署名和标识。当前Word可直接交给排版人员做内容依据；下一步按用户后续要求完成大幅版式。
- 延续事项：10/18数量门控、时间/星历参考维护及候选质量核验依照生产主线另行处理，本次展示材料不代表问题已修复或恢复上报。

## Current priority: human review routing and MPC correction request (2026-10-11)

- 用户最新规则：形态异常直接拒绝、不做 PSF；形态正常才拟合 PSF，PSF 异常拒绝。两关正常的 pass 和证据不足的 review 均进入人工核验；reject 不进。适用于全部 unknown。
- 已部署独立路由模块，不改变数值引擎/cache：747 条可人工核验，842 条拒绝；人工接受是最终提交的必要条件，旧标签不能绕过明确拒绝。
- 每天 09:00 和 reboot cron、单夜 Step6b、已有包续跑入口已核实；续跑会补筛选并备份刷新未上报旧包，保留人工决定；有提交记录的历史包保持。
- 验证：本机/服务器各 16 项门控，39 项隔离评分/提交回归，本机 3 项 daily/hold 已通过；真实复核包刷新/导出和 213 文件保护校验进行中。
- MPC 用户已授权直接联系、撤回或标错并请求台站 327 恢复处理条件；不代表恢复 astronomy submissions。已准备 33 条/105 次撤回/标错清单 + F5 的 3 次继续暂停，见 `outputs/mpc_request_20261011/`。
- 当前外部阻塞：Jira 浏览器操作持续超时；用户指定使用 wsgp2024@163.com，明确先不用 Gmail。163 尚未连接；已安排登录页，等待用户完成登录，不索取聊天明文密码/验证码。没有发出消息。
- 用户最新限制：先展示正文并获明确确认后才可发送邮件或 Jira；163 登录尚未完成，用户提供的登录凭据不得写入项目文件/日志/Git。
- 下一步：结束生产刷新/验收并提交里程碑；正文确认后再登录 163、核实 Reply-To、发送并记录回执。实际隔离 CLI pass+review 验收已 local_validated/network_requests=0。继续保留 MPC hold，不开展误拦率研究或声称天测全面修复。

## Previous milestone: serial morphology then PSF acceptance (2026-10-11)

- 最新明确要求：全部 unknown 先做形态；异常直接 reject 且不拟合，正常才进入 PSF；PSF 异常 reject，两关正常才 pass。证据不足 review 隔离，不进入常规人工候选；不限定 NEO。
- 数值阈值与至少 3 帧、>=2/3 异常规则保持；新 policy `serial-morphology-psf-20261011` 自动使旧缓存失效。上版“双异常才剔除”结果保留作对照。
- 已完成：本机/服务器各 12 项针对性回归，13:35 CST 备份后部署，约 14:25 CST 完成 20260830..20261010 的 22 夜重筛；8 个无保存 catalog 的夜目录明确保留为缺失，未重建 linking。
- 全体结果：1589 条中 reject842/pass274/review473；形态拒绝 761/不足 101，727 条正常进入 PSF 后拒绝 81/通过 274/不足 372。剔除率 53.0%，21 个非空夜平均 40.1 条、中位数 37；通过仍须人工确认。
- 已上报结果：47 条为 reject31/pass16/review0（规则命中 66.0%）；NEO 10/11、unclassified 4/4、历史未分流 unknown 17/32 被拒绝。旧假源 nP 新增被 PSF 拦住，oI/oj 仍双通过，不能宣称鬼像问题已完全解决。
- 路径：服务器部署/旧结果备份 `serial_morphology_deploy_20261011.s1t1upg5`、回放 `serial_morphology_replay_20261011`（均在 `/pipeline/xiaoyunao/data/heliolincrr/`）；本地执行凭据 `outputs/serial_morphology_20261011/`。
- 验收完成：1589 条串联路径重算，形态测量与上版一致；22 夜输入有效、202 个原始/人工/提交文件未变、3 个生产文件 hash 匹配。实际复核包 9 条/27 观测/9 GIF；旧人工选择导出 2 条/8 观测并完整本地校验。旧筛选备份、MPC hold 保留，无正式/测试上报或撤回消息。
- 交付：`outputs/serial_morphology_20261011/audit/` 的逐夜分关统计、47 条已上报明细、31 条规则命中待看图名单、逐帧指标及核验记录；完整 tar.gz SHA256 已核对，源码提交 `d0e70ce8`。
- 真源误拦率研究仍按用户要求暂不开展；MPC 撤回问题只核实官方流程，等待明确发出请求的授权。
- 下一步：本次工作完成；保留漏检证据，按后续用户要求决定 MPC 联系、规则改进及真源保留率研究。

## Superseded milestone: joint morphology/PSF rejection and history replay (2026-10-11)

- 最新授权：直接将轨迹级形态 AND PSF 规则接入主流程，从 20260830 起重筛已有 unknown；真源误拦概率研究按用户要求留到后续。
- 固定采用 FWHM ratio>1.5、Class_ANN<0.1 和 PSF width>1.3、delta_chi2>=9；两套判据分别要求至少 3 帧、>=2/3 异常，不要求异常帧完全相同。
- 代码已加入 raw 旁通过/隔离/审计输出，复核包、follow-up 关联、导出/正式 transport 共用门控；历史人工标签/PSV/账本保留，MPC hold 继续有效。
- 已完成：12:34 CST 原子部署并核对 hash；20260830..20261010 共 22 个有保存 catalog 的夜次回放成功。1589 条候选 = pass839/reject565/review185；8 个无保存 catalog 的夜目录明确记录，不重建 linking/orbit。生产源码与分支 `codex/unknown-morphology` 一致。
- 提交对账：正式收件证据与冻结 catalog/PSV 逐项对应，共 47 条 = pass25/reject22/review0；22/47=46.8% 是规则命中率。MPC 此次 11 个 NEO 候选中拦住 9 个，00001EM/00001F5 未触发，不据此确认真实。
- 已知局限：历史人工确认假源 00001nP/00001oI/00001oj 未满足形态多数，未进入 PSF 而通过；00001o2 拦住。明细保留人工标签/拟合执行情况，pass 不能当作真源。真源误拦概率与规则扩展按用户要求留待后续，不在本次继续试验。
- 验收完成：10 项新回归双平台、3 项 daily/hold、服务器独立 runtime 39 项评分/transport 回归；真实旧人工接受集导出由 6 条减至 4 条/12 观测且本地 ADES/digest2 校验通过。202 个原始/人工/提交文件未变、9 个部署文件 hash 一致、1589 个判定重新核算、22 份输入签名有效；本次无网络上报。
- 当前产物：`outputs/morphology_replay_20261011/` 含 CSV/逐帧报告/验证记录/完整审计归档；部署凭据 `outputs/morphology_deployment_20261011/`。服务器部署与备份 `/pipeline/xiaoyunao/data/heliolincrr/morphology_deploy_20261011.sm1jt1/`，回放 `/pipeline/xiaoyunao/data/heliolincrr/morphology_replay_20261011/`。
- 下一步：接入及回放任务已完成；保持 MPC hold，等待用户另行要求真源误拦研究或授权恢复上报。以下旧试验状态仅作历史依据，以本节部署结果为准。

## Current priority: known ADES rejection and daily recovery (2026-10-11)

- 本对话继续处理 known ADES 拒收与积压恢复；鬼像/PSF 筛选留给用户回到上一对话处理。用户已澄清封禁邮件系候选后随未确认，不再按已知身份误判方向追查。
- 已修复/部署 known 导出质量和本地 XSD 门控；CF58/CFL1 对应 20261005/06 各 4 条质量标志异常观测，正的小误差被格式化为 0.0000。14 项双平台回归通过。
- 两夜隔离纠正稿分别 6152/6659 条，全部保留行与旧报告一致，每对象 >=3；完整 XSD 通过。原 PSV/回执保留，MPC 状态仍为 rejected，未重报；纠正稿、审计、拒收记录在 incident_20261011_f1q_hiwq/known_ades_rejection_dsfe87d4。
- 已止血：服务器 `/pipeline/xiaoyunao/data/MPC_SUBMISSION_HOLD` 已启用；known/unknown 正式 transport 均阻断，提取和计划可继续。未发送 Jira/邮件，未恢复提交。
- 已修复并部署：digest2 `100.00000000000001` 端点舍入使整夜打包失败；daily 现在按夜隔离失败并保留非零退出/夜次错误。39 个评分/上报回归通过，89 条真实回放成功，88 个其他分数不变。
- 计划核查：20261010 471 条、20261011 474 条已生成并发布到 `/pipeline/xiaoyunao/script/`；待确认观测助手查看的是望远镜计划还是候选复核页面。
- 已恢复并验收：20261006 package，89 个候选/269 条 obs/89 GIF，全部评分成功且 tar 完整。正在恢复：20261007/08/10 先 known 后 unknown，保持提交暂停。20261009 无 MP L2。恢复记录 `/pipeline/xiaoyunao/data/heliolincrr/incident_20261011_f1q_hiwq/recovery_status.json`。
- 验收标准：每夜 known manifest 和完成日志一致、无 query/match 错误；复核包包含全部 GIF、manifest、评分；历史 MPC 账本和人工状态不变；故障夜不阻塞后续夜且失败明确可见。
- 尚未解决：MPC 候选测量质量与重新获准处理；旧 astorb 更新作为历史维护事项保留；已发现旧 known status 仅凭 partial merge 也可称 complete 的历史风险尚未改动，本次 recovery 另作逐帧日志审计。保留原始数据和上报证据，不以本次软件恢复宣称天测问题已解决。
- 下一步：完成补跑验收并记录结果；在明确整改证据和恢复上报授权前维持 hold。

## Historical experiment: PSF confirmation and joint morphology rule (2026-10-05)

- 当时完成隔离 PSF 对照；2026-10-11 用户已授权采用并部署，详见本文件首节。沿用普通恒星 PSF，按用户要求不另开发快速目标拖线 PSF；拖尾识别另有图像方案。
- 基线证据：`outputs/psf_feasibility_20261005/`。单帧 `FWHM/local_star_FWHM > 1.5 && Class_ANN < 0.1`，至少 3 帧且 >=2/3 帧异常；5 个原候选及 Cm 均触发。975 张 L2 / 16,982 条 known 检测中，314 条缺参考星；6353 个有效连续三帧窗口有 11 个误拦，涉及 8 个 object-night。
- 新证据：`outputs/psf_confirmation_20261005/`。全部 19 次候选观测 + 11 个误拦窗口的 27 次 known 观测；45/46 PSF 有效，225 次留出恒星拟合收敛。另一个 Honeycutt 曝光参考星不足，显式未评估；它不影响所属窗口多数判定。
- 有希望的组合是在轨迹级作 AND：先满足基线，再要求全部观测中 >=2/3 通过 PSF 宽化判据 `width_scale>1.3 && delta_chi2>=9`（探索值）。5 个原候选 + Cm 都仍触发，11 个已知误拦均消除；原先未触发的 6342 个窗口不可能被附加 AND 新增剔除，因此未对它们重做 PSF。
- 不能改成严格的“同一帧四条件全满足”多数票：同门限下 CB 会被放过。不能只用 reduced chi-square/qfit；与暗弱真源重叠，delta_chi2 也不是已校准显著性。
- 尚未解决：方案在同一小样本上探索且窗口重叠，不构成零误杀或前瞻完备率保证；Bx 一帧贴近宽度阈值，100 次模板重采样仅 82 次仍过阈值；宽度提高至 1.35 或 delta_chi2 至 16 分别放过 Bx/CB。5 个候选的 ghost 身份仍待独立确认。
- 后续研究（用户要求暂不开展）：在独立真实及伪源上检验；区分普通点源、拥挤混合、拖尾与 PSF 缺失，验证与已有拖尾识别的衔接及真实快速/暗弱目标保留率。当前固定参数已在 candidate 识别之后、人工复核之前接入。
- 保留原试验边界：需要独立留出样本和明确误拦/未评估分母才能估计保留率；跨帧映射、模板自包含排除、拟合失败处理已纳入本次软件验收，不把历史小样本表现等同生产精度保证。

## Paused by user (2026-10-11): reported-candidate image audit (2026-10-05)

- 已取回 `00001AE/BL/Bx/CB/Ck` 五张原始 GIF（16 帧），生成手机可看逐帧总图与 ZIP：`outputs/tracklet_audit_20261005/`。只读审计结束，生产状态未改。
- 五源的图像/编号/提交映射一致；局部恒星配准确认相对运动，排除“只因整场 WCS 漂移而把静止星连成链”的简单解释。
- 主要疑点是非恒星移动伪影：像斑多数偏宽，探测器位移与恒星漂移近乎反向，Bx 次夜又在 Cm 轨迹处近似重复；这些证据尚不能代替母星/反射几何或原始图减成因确认。
- 近期步骤：检查 Bx/Ck/Cm 同场全曝光及 L1/L0，追踪母星或随指向变化的 ghost；取得同行跟随目标/时间/覆盖范围后核对预测位置。
- 验收标准：用独立图像证据确认真伪与成因；任何形态/ghost 新门控先在隔离回放中验证，对已知小行星报告保留率，不能仅以 FWHM 阈值或高 digest2 分数定性。
- 尚未实施的后续操作：改变人工真假决定、撤回 MPC 候选、停用 follow-up、发送同行消息；先完成成因确认。

## Digest2 post-deployment observation (2026-09-23)

- 20260921：15 个 unknown 均完成本地评分，13 个 NoID NEO ≥65；人工最终全部判否，正式上报数 0，watcher 状态 `no_observations`。
- 20260922：unknown 数为 0，空评分 CSV/复核包正常生成，daily 正常结束。20260920 无观测目录，正常跳过。
- 两日相关日志无 ERROR/WARN/Traceback；9 个生产文件 hash 与部署记录一致。上一轮未完成的服务器说明文档同步已于 2026-09-23 完成。
- 当前待验：首个被人工接受的未知目标进入正式拆批、MPC 收件和后续邮件处理；不能把这两天空的上报结果当作正式接口已验收。
- 审计原始摘要：`outputs/digest2_health_20260923/remote_report.json`；下次有新提交时先查 `/processed1/<night>/L4/<night>_unknown_mpc_batches.json` 与 `review_submit_state.json`。

## Follow-up linkage investigation (2026-09-20)

- `00001nF` 的 20260909/10 后随 L2 数量达标，但没有完成跨晚目标确认：`associated_links=[]`，运动模型仍只有 20260902 的 3 个测量点。
- 主链只生成逐晚 tracklet/link；回挂步骤仅匹配 L4 unknown link，且只考虑 `active` 源。20260910 数据进入回挂前，本源已被 L2 数量判为 `complete`，产生时序漏检。
- 下一步先人工/离线确认后随图像中的同一移动源，并检验三夜轨道一致性；随后让回挂可处理刚 `complete` 的源，避免把“图像数量达标”误当成“链已关联”。

## Completed milestone: digest2 / NEW NEO production integration (2026-09-20)

- 最新授权已覆盖此前隔离测试限制；主流程接入和服务器部署完成，详见 `heliolincrr/DIGEST2_PIPELINE_UPDATE.md`。
- 流程：复核包全部预评分 -> 完整人工决定 -> 最终 ADES schema 校验与重算 -> NoID NEO >=65 使用 neo candidate，否则 unclassified -> 逐批正式回执账本。
- 保留真假门控、统一 masked JSON/FITS 和原后随入口；失败阻断，不以高分证明真源。--validate 为本地校验，不发请求。
- 36项双平台测试通过；真实20260910临时目录完整回放120目标372观测120 GIF，最终4目标18观测进入候选批；部署后真实入口本地验证通过。
- 历史7个submitted夜次保持终态；原共享state hash不变；此次正式POST=0。独立runtime与备份均位于 `/pipeline/xiaoyunao/data/heliolincrr/`。
- 交付：8页PPT在 `outputs/digest2_production_20260920/deliverables/`；MD记录本地原理、固定版本细节、运行路径与异常处理。
- 收尾：部署及烟测完成后 SSH 曾超时；额外服务器文档镜像已于 2026-09-23 同步。PPT/MD 本地交付已完成。
- 后续：首个新夜正式邮件/NEOCP状态需运行后核对；网页分数列、已有NEOCP身份自动路由、模型定期升级与曝光中点约定另行处理。
- 验收已完成：65边界、真假门控、完整观测与hash、schema错误、部分失败/并发去重、历史状态保护、部署路径烟测、PPT逐页检查。

## Previous objective: verify first follow-up night (2026-09-04)

- 20260904 最终计划包含 7 个 active source、每源 5 条，共 35 条 `MP_FU_*`；6 源来自 20260830，20260902 只保留可能真实的 `00001nF`。
- 下一验证：20260905 早晨逐源统计这 7 个源的昨晚 `MP_FU_*` L1/L2，确认每源是否达到 5 张并检查 reconcile 是否记录 success。
- 已修复：review watcher 调用 `survey.apply_followup` 时明确使用 survey 默认 Python；服务器 dry-run 和新 watcher 命令行已验证。
- 已重画：20260904 plan 图现在基于人工停用后的最终 433 行计划，单独标出 35 条 follow-up；绘图代码能把 `MP_FU_*` 映射回实际 selected footprint。
- 已修复：多个单夜 watcher 通过 state 专用跨进程锁串行执行；每轮在锁内 reload，并将完成检查、MPC 提交和状态写入作为一个事务。0830/0901/0902 首次有效 submission 记录已从日志恢复。
- 已验证：20260903 的 106 条候选全部判否，watcher 记录 `no_observations` 后正常退出；共享 state 同时保留此前五夜记录。
- 新发现：20260902 的最终人工 review 把 linkage 226/446 标为 0，但更早生成并实际提交的 submit snapshot 把两者标为 1；需给 review->submit 交接增加 accepted 清单确认、源文件 hash/mtime 和审计记录。
- 已确认并处理：`00001nP`, `00001o2`, `00001oI`, `00001oj` 均为假源，已从 follow-up state 与今晚计划停用；只有 `00001nF` 继续监测。
- 3000-row 门控在 20260903 成功拒绝 17 个同类污染帧，保留帧最大残留 2518；继续监测正常密场误拒率和 unknown 数量。
- 后续候选：给 follow-up 增加月距约束与线性外推不确定度门控；为陈旧的 `astorb.dat` 建立安全自动更新与版本记录。

## Completed milestone: Gaia-residual unknown protection (2026-09-02)

- 目标已完成：unknown 链生成前默认丢弃 Gaia 去除后残留数 `>2000` 的整帧星表。
- 已部署并与本地仓库对齐：
  - `heliolincrr/mask_gaia.py`
  - `heliolincrr/run_single_night.sh`
  - `heliolincrr/run_daily_unknown.sh`
- 历史夜恢复已完成：
  - `20260830`: 11 rejected frames，41 unknown links，41/41 GIF，125 review rows
  - `20260901`: 20 rejected frames，56 unknown links，56/56 GIF，169 review rows
- 验收已通过：保留帧残留数均不超过 2000；rejected 帧未进入最终 unknown；两份 review tar 可完整读取且 manifest 无 missing GIF。
- 当前剩余工作：
  - 人工填写两夜 review/submit CSV；本 milestone 未启动自动 MPC 提交
  - 上游 L1 astrometry 仍需增加 pixel scale、NMATCH、RMS 和跨曝光一致性质量门控
  - 20260829 需先修复/重跑因 IERS 过期而不完整的 known subtraction，再决定 unknown 重建
- 下一步：完成人工复核后使用现有 reviewed-unknown watcher 导出、验证并提交；继续追踪上游天测修复。

## Production incident follow-up (2026-09-01)

- 已确认 20260830/31 今早两次 known 上报分别对应两个观测夜，不是重复提交；20260830 是前一日 final export/submit 未完成后的 recovery 补交。
- 已排除：服务器断电、kernel OOM、Slurm preemption 和 13:40 才启动的图减 daily job 对 09:10 finalizer 的资源挤占；首次 finalizer 精确异常仍因日志丢弃而不可恢复。
- 已定位 unknown 爆量主因在上游 L1 天测质量：
  - `OBJ_MP_1431_0072` 和 `OBJ_MP_1430_0046` 是 pixel scale/NMATCH/photometric RMS 均异常的灾难性 WCS 错解
  - 其余异常帧有数角秒级空间畸变差，超过 Gaia mask 的 1.5 arcsec 半径
  - IERS 缓存过期是 20260829 known matching 的另一项独立故障；缓存已于 8 月 30 日自动刷新
- 当前待修复：
  - 为 known finalizer 保留独立 stdout/stderr、退出状态和 submission-attempt 记录
  - 为 MPC curl 增加有限超时、明确失败状态和安全重试，避免 daily 静默等待 8 小时
  - 在 known status 中校验 expected manifest parts 与 completed parts，不能把部分成功写成 `known_complete=true`
  - unknown 前增加 per-frame Gaia-mask 残留数异常门控，避免少数坏星表触发 tracklet/unknown 爆炸
  - 固定或自动刷新 IERS 数据，避免 20260829 的 192 帧 known query/match 静默失败
- 数据恢复：
  - 20260831 review package 已正常生成：57 links / 178 observations / 57 GIFs
  - 20260829 需先完成 known rematch；20260830 需先隔离异常 Gaia-mask 帧，再评估重建 unknown
- 验收标准：
  - 合成 finalizer 网络失败时能在持久日志/status 中看到原因，并在有限时间退出
  - 缺 parts 时 known status 明确失败，unknown 不启动
  - 异常单帧被隔离并报告，正常夜仍能自动发布 review package
  - recovery 不会对已有有效 reply 的同一夜重复提交
- 下一步：用户确认后实施上述最小修复，在临时目录回放 20260829/30，不直接覆盖生产结果。

## Corrected GOTTA-source figure bundle (2026-08-30)

- 已完成并验收独立目录 `SHARP_local_analysis_figures_20260830/`：
  - 12 张直接依据 GOTTA 源码重绘的 PDF/PNG 正式图；Figure 05c 已按作者要求删除
  - Figure 01 使用 GOTTA 原始 zscale + 61-pixel median background subtraction，无 N/E
  - Figure 02 恢复完整 SHARP production architecture 和 feedback loop
  - Figure 03 使用 GOTTA rainbow/LogNorm HEALPix density map
  - Figure 04 使用 GOTTA 原始八类 asteroid marker/color/order；04d/04e/04f 使用 GOTTA colorbar/marker grammar
  - Figure 06 仅为作者指定 GIF 的原始三帧 1x3，不增加内容
  - 58-linkage candidate table、58,482-row orbit join audit、orbit-class summary 和完整 figure data
  - GOTTA 源码快照、逐图映射、code/catalog/baseline/image-source audit、QA、运行脚本和全文件哈希
- 已确认的边界：
  - 不修改 manuscript 或生产代码
  - 当前 orbit confirmation 仍是 868.221 m，scheduler/known matcher 是 960 m
  - 无可靠 PHA 标志；NEO 仅按 `q < 1.3 AU` 派生
  - 缺时域 calibrated-source manifests，因此 Figures 07/08 合理缺席
- 验收：12/12 单页 PDF 通过 Poppler 渲染，12 张全尺寸渲染和 contact sheet 人工复核通过，candidate table 58 行，checksum manifest 全量复验通过。
- 下一步：上传新的桌面 `SHARP_local_analysis_figures_20260830.zip` 给 GPT Pro；旧 ZIP 不再使用。

## Completed side task (2026-08-27)

- 已生成 GPT Pro 论文续写交接包：
  - `/Users/yunaoxiao/Desktop/SHARP_PASP_GPT_handoff_20260827.zip`
  - 包含完整原始 LaTeX 工程、12 张正式图的 PDF/PNG、5 张正式表的 CSV/TeX、最终分析报告、紧凑证据摘要、QA 和 executed notebook
  - `START_HERE.md` 给出修订顺序、必须遵守的统计口径、代码冲突和外部 TODO
- 已重新解压验证：147 个受管文件哈希全通过，图表数量和关键文件齐全，ZIP 无压缩错误。
- 下一步：上传给 GPT Pro；论文修订以 `analysis_support/` 冻结结果为准，不让其自行补造 MPC、站点、天气或性能证据。

## Completed side task (2026-08-26)

- 已生成 PowerPoint 推荐 MP4：
  - `/Users/yunaoxiao/Desktop/known_asteroid_cutout_top40_mosaic_5x8_ppt_2560x1600_crf18_20260826.mp4`
  - H.264 High/yuv420p，`2560×1600`, `5 fps`, `60` 帧, `12 s`，约 `80 MB`
  - SSIM `0.995937`、PSNR `44.57 dB`；比 200 ms GIF 小 `58.3%`
  - PowerPoint 中设置 `Start: Automatically` 和 `Loop Until Stopped`
- 全分辨率 `3840×2400` CRF 16 试验版仍约 `196 MB`，仅作为备选，不推荐嵌入 PPT。
- 已在不重新编码画面的前提下生成 `200 ms/帧` 版本：
  - `/Users/yunaoxiao/Desktop/known_asteroid_cutout_top40_mosaic_5x8_60frame_200ms_20260826.gif`
  - `60` 帧，总循环 `12 s`；与 0 ms 版的 60 个解码画面逐帧一致
  - 清单：`/Users/yunaoxiao/Desktop/known_asteroid_cutout_top40_mosaic_5x8_60frame_200ms_20260826_manifest.csv`
- 已按修正后的逐帧取模规则生成 60 帧、0 ms 的完整合成 GIF：
  - `/Users/yunaoxiao/Desktop/known_asteroid_cutout_top40_mosaic_5x8_60frame_0ms_20260826.gif`
  - 每个全局帧均为完整 `3880×2425` 合成图；格子使用 `global_index % used_frames` 取帧
  - `60` 个 GIF image blocks 和 GCE delay=`0` 已逐块验证；2400 次格位映射检查通过
  - 取帧清单：`/Users/yunaoxiao/Desktop/known_asteroid_cutout_top40_mosaic_5x8_60frame_0ms_20260826_manifest.csv`
  - 静态预览：`/Users/yunaoxiao/Desktop/known_asteroid_cutout_top40_mosaic_5x8_60frame_0ms_20260826_preview.jpg`
- 已从 50 张 daily known asteroid cutout GIF 中逐帧精选 40 张，并生成无缝 `5×8` 循环大 GIF：
  - `/Users/yunaoxiao/Desktop/known_asteroid_cutout_top40_mosaic_5x8_20260826.gif`
  - 画布 `3880×2425`，`22` 帧，`500 ms/帧`，总循环 `11 s`，无限循环
  - 不同源帧数按完整播放相位统一映射；未简单截断
  - 筛选/格位清单：`/Users/yunaoxiao/Desktop/known_asteroid_cutout_top40_mosaic_5x8_20260826_manifest.csv`
  - 静态预览：`/Users/yunaoxiao/Desktop/known_asteroid_cutout_top40_mosaic_5x8_20260826_preview.jpg`
- 已将最近 10 个完整 daily known asteroid 夜次的 50 张 cutout GIF 打包到：
  - `/Users/yunaoxiao/Desktop/known_asteroid_cutout_gifs_latest50_20260826.tar.gz`
- 归档内含按夜次组织的 GIF、选择说明和逐文件 SHA-256 manifest；本地完整性验证通过。

## Current objective

SHARP PASP 的本地可完成分析已收口在 `paper_analysis_20260803/`。当前目标从“生成统计与图件”转为“等待作者/上游证据后冻结论文口径”。

- 已完成：
  - `20251115..20260715` raw/L2/known/unknown/plan/log/review 冻结清单与哈希
  - known recovery/astrometry/random-shift、unknown full funnel、cross-night/JPL、scheduler、operations、site sensitivity 和 ADES midpoint correction
  - 5 张 paper table、12 张 PDF/PNG 图、executed notebook、contact sheet 和机器验证报告
  - 555 个 artifact 的最终 snapshot manifest；63 PASS / 0 FAIL / 1 external SKIP
- 冻结解释：
  - 68 human-marked-real -> 67 submission-selected -> 58 post-audit-retained single-night linkages
  - 58 不能写成独立天体或 discovery count；37 只是线性运动候选 component
  - 6 条 C/2025 Y1 是 JPL 数值关联候选，不是 MPC attribution/designation
  - 960 m 替代站点重跑无 gate flip，但短弧轨道元素不可作 population inference
  - scheduler `95.1743%` 只适用于 plan-active acquired cohort，不是 observatory efficiency
- 当前外部阻塞：
  - 作者/台站确认 surveyed MPC 327 lon/lat/height/datum 和历史-as-run/未来代码修正政策
  - 作者签字确认 quality mask 与异常夜主文/附录策略
  - 上游 reduction、depth、effective detector mask、GMG header 修复和 image-injection 接口
  - MPC 58+9 逐行 ingest/identity/designation 记录，以及 C/2025 Y1 和 unknown 15 s midpoint correction 的处理政策
  - 独立 CPU/RAM、人工时间、天气/设备/override 日志
- 下一步：
  - 使用 `paper_analysis_20260803/reports/AUTHOR_INPUTS_REQUIRED.md` 收集最小回复
  - 收到外部证据后只更新相应冻结表、重跑受影响的表/图/validator/assembler
  - 论文文字修订以 `MANUSCRIPT_CODE_CONFLICTS.md` 为准；本分析任务不直接修改论文工程

## Deferred reference

2026-06-26 Agent Mail 通知候选方案，暂缓：

- 目标：当网页 check 产生 `<night>_submit.csv` 且某条 unknown link 被标记 `is_real=1` 后，向关注者邮件通知，邮件中包含 link 的轨道/观测摘要，方便团队及时关注和 follow-up。
- 推荐挂点：
  - `heliolincrr/watch_submit_reviews.py`
  - 在 `run_one()` 成功导出/validate/submit 后、`run_followup_update()` 附近触发通知
  - 使用 watcher state 对 submit CSV signature + review manifest signature 去重，避免重复发送
- 推荐邮件内容：
  - subject: `[SMT unknown] <night> <trkSub> passed review`
  - body: night, trkSub, linkage_id, n_obs, n_tracklets, rms/median/max residual, a/e/i, linear speed/direction, RA/Dec/MJD span, MPC validate/submit status, review package path, GIF path
  - 可选附件：对应 GIF；若附件策略复杂，先只写服务器绝对路径和 review package 路径
- 数据来源：
  - `/processed1/<night>/L4/<night>_unknown_links.json`
  - `/processed1/<night>/L4/<night>_unknown_links_submit_stats.json`
  - `/pipeline/xiaoyunao/heliolincrr/review_packages/<night>/<night>_unknown_review_manifest.json`
  - `/pipeline/xiaoyunao/heliolincrr/review_packages/<night>/gifs/<trkSub>_link....gif`
- 部署分歧：
  - 服务器当前 PATH 未发现 `agently-cli`, `npm`, `node`
  - 若要由 daily watcher 直接发信，需要在 `smtpipeline@www.xinglong-naoc.cn` 上安装并授权 Agent Mail CLI
  - 若不想把授权放服务器，可本机跑一个轮询 worker，SSH 读取 watcher state/manifest 后用本机 Agent Mail 发信
- 待决策：
  - 收件人列表
  - 每个 true link 一封，还是每夜 digest
  - 邮件是否附 GIF
  - 是否只通知 validate/submit 成功，还是网页 check 一通过就通知

2026-06-24 follow-up 调度入口已完成并部署：

- 新增 `survey.followup` 状态机和 `survey.apply_followup` 插入器
- follow-up 只处理 `20260624` 及之后新完成网页 check 的 true unknown 源，不回追历史 backlog
- 状态文件：
  - `/pipeline/xiaoyunao/survey/runtime/followup/followup_state.json`
- 调度规则：
  - 每个源需要两个实际满 `5` 帧的 follow-up 夜
  - 少于 `5` 帧的夜记为 failed，后续继续排
  - 从发现夜起满 `10` 天仍未完成则 abandoned
  - 目标名为 `MP_FU_<trkSub>_MP####`，既区别于普通 MP，又能进入 known/unknown 链条和历史曝光统计
  - 预测位置用当前 link detections 的线性 RA/Dec 模型；选择已有 footprint 中覆盖预测位置且中心距最近的 field
- 自动入口：
  - `survey/run_daily.sh` 生成基础计划后默认调用 `survey.apply_followup`
  - `heliolincrr/run_daily_unknown.sh` 在 unknown link 生成或已有 package skip 时调用 `associate_followup_links.py`
  - 单夜 `watch_submit_reviews.py` 由 daily unknown 启动时带 `--enable-followup`，网页 submit 处理完成后立即更新 follow-up 状态和当天 `current_plan.txt`
- 已部署到服务器：
  - `/pipeline/xiaoyunao/survey/followup.py`
  - `/pipeline/xiaoyunao/survey/apply_followup.py`
  - `/pipeline/xiaoyunao/survey/run_daily.sh`
  - `/pipeline/xiaoyunao/heliolincrr/associate_followup_links.py`
  - `/pipeline/xiaoyunao/heliolincrr/watch_submit_reviews.py`
  - `/pipeline/xiaoyunao/heliolincrr/run_daily_unknown.sh`
- 验证：
  - 本地和服务器 `py_compile`/`bash -n` 通过
  - 本地 smoke 合成 true link 后插入 `5` 条 `MP_FU_...` 行
  - 本地 smoke 合成 unknown link 后 `associate_followup_links.py` 返回 `matches=1`
  - 服务器真实计划 dry-run：`active_sources=0`, `inserted_rows=0`，未改写生产计划
- 待验证：
  - 等 `20260624` 之后第一条网页 submit true 源出现后，检查 `followup_state.json` 和 `current_plan.txt`
  - 第二天检查 L2 中 `MP_FU_...` 文件数是否触发 success/failed reconcile

2026-06-24 reviewed unknown 补报完成：

- `20251116..20260617` 历史 reviewed unknown 补报已经全部完成
- state `/pipeline/xiaoyunao/data/heliolincrr/review_submit_backlog_20251116_20260617.json`
- 最终 summary：
  - `review_packages=122`
  - `complete=122`
  - `pending=0`
  - `failed=0`
  - `submitted=35`
  - `no_observations=87`
- 顺序上从 `20251116` 连续完成到 `20260617`
- submit CSV 共 `116` 个、`4743` 行
  - `is_real=1`: `67`
  - `is_real=0`: `4676`
  - blank/invalid: `0`
- watcher 日志尾部已写出 `status=complete`，历史补报 watcher 进程已退出
- 后续只需要处理新 daily unknown check package，不再需要历史区间 backlog watcher 常驻

2026-06-23 断电恢复检查：

- 服务器 `2026-06-23 12:28 CST` 重启；上午断电属实
- `2026-06-23 09:00 CST` cron daily pipeline 在断电前正常跑完
- `2026-06-23 12:38 CST` `@reboot` daily pipeline 自动补跑完成
- `20260623` 观测脚本正常生成并发布，`n_exp=285`
- 昨晚目标夜 `20260622` 没有 processed data：
  - `/processed1/20260622` 不存在
  - known/unknown 均因 missing night dir 正常 skip
- 历史补报 watcher 已由 `@reboot` 自动恢复：
  - 当前 PID `12242`
  - 当前只有一个 watcher 进程，无重复运行
- 当前 `20251116..20260617` reviewed unknown 补报进度：
  - `review_packages=122`
  - `complete=76`
  - `submitted=19`
  - `no_observations=57`
  - `pending=46`
  - `failed=0`
  - submit CSV 共 `70` 个、`3138` 行
  - `is_real=1` 为 `35`，`is_real=0` 为 `3103`
  - 无空标签/非法值

2026-06-22 断电恢复检查：

- 服务器 `2026-06-22 10:20 CST` 重启，说明 `09:00` cron 当时错过
- `@reboot` daily pipeline 已在 `2026-06-22 10:31 CST` 自动补跑完成
- 今日观测脚本正常生成：
  - `/pipeline/xiaoyunao/survey/runtime/plans/20260622_plan.json`
  - `/pipeline/xiaoyunao/survey/runtime/plans/20260622_plan.txt`
  - `n_exp=282`
- 昨晚目标夜 `20260621` 没观测：
  - `/processed1/20260621` 不存在
  - known/unknown 均因 missing night dir 正常 skip
- recovery 回看 `20260615..20260621`：
  - `20260615/16/18/19/21` 缺 processed night dir，正常 skip
  - `20260617/20260620` 已有 known report 和 unknown review package，幂等 skip
- 历史补报 watcher 在断电后未自动恢复，已手动重启并加入 `@reboot`：
  - 当前 PID `201618`
  - crontab 新增 `@reboot sleep 900 && ... run_review_submit_backlog_watch.sh 20251116 20260617`
- `run_daily_unknown.sh` 已修复已有 review package 的恢复场景：
  - 已有正 unknown package 时不重跑 unknown
  - 会确保该夜 submit watcher 被启动
  - zero unknown package 不挂 watcher
- 当前 `20251116..20260617` reviewed unknown 补报进度：
  - `review_packages=122`
  - `complete=68`
  - `submitted=14`
  - `no_observations=54`
  - `pending=54`
  - `failed=0`
  - submit CSV 共 `62` 个，`2794` 行，`is_real=1` 为 `22`、`is_real=0` 为 `2772`，无空标签/非法值

2026-06-22 补充：

- 已生成项目状态英文单页 PPT：
  - `/Users/yunaoxiao/Desktop/smt_asteroid/outputs/smt_asteroid_project_status.pptx`
  - 桌面副本已修订：`/Users/yunaoxiao/Desktop/smt_asteroid_project_status.pptx`
  - 中文版本已按 DESI MWBP 参考 deck 风格重调：`/Users/yunaoxiao/Desktop/smt_asteroid_project_status_cn.pptx`
  - 版式参考 `/Users/yunaoxiao/Desktop/desi_mwbp_et.pptx` 的白底学术汇报风格
  - 内容概括 nightly scheduler、raw images、collaborator preprocessing/calibration/photometry、L1/L2 catalogs、known-object extraction/reporting、unknown search/review
  - 右侧放大嵌入服务器 known asteroid GIF 和本地 checked unknown candidate GIF
  - 底部状态采用高层统计：`354237` 次 known asteroid detections、`58224` 个 known objects、`22` 条 submitted reviewed unknown links
  - 桌面 known GIF 已替换为 `top03_Bakhchisaraj_20260528.gif`，并改为与 unknown GIF 一致的 `3` 帧、`3000 ms` 总时长
  - known/unknown 两个模块底部说明框已恢复上一版文字
- 若继续改 PPT，优先微调：
  - 标题和状态句措辞
  - 右侧示例 GIF 选择
  - 统计数字需按服务器最新产物刷新

2026-06-21 最新状态：

- Daily automation 已切到服务器 crontab：
  - `0 9 * * * cd /pipeline/xiaoyunao && /bin/bash /pipeline/xiaoyunao/heliolincrr/run_daily_pipeline.sh >> /pipeline/xiaoyunao/data/heliolincrr/daily_logs/cron_daily_pipeline.log 2>&1`
  - `@reboot sleep 600 && cd /pipeline/xiaoyunao && /bin/bash /pipeline/xiaoyunao/heliolincrr/run_daily_pipeline.sh >> /pipeline/xiaoyunao/data/heliolincrr/daily_logs/reboot_recovery_daily_pipeline.log 2>&1`
  - 旧的 `survey/run_daily.sh` 和 `known_asteroid/run_daily.sh` cron 已移除，避免 daily 入口不一致
- `2026-06-21 09:00` unknown 未触发的原因：
  - cron 还停留在旧的 survey/known 两条入口，没有安装 `run_daily_pipeline.sh`
  - known extraction 实际已触发并完成 `20260620`
  - known ADES/reply 初始未生成，是因为服务器 `known_asteroid/export_ades.py` 旧版本不支持 `slurm_merge_submit.sh` 传入的 `--dedup-history-root/--current-night` 参数，finalize 静默失败
- 已把正确的 `export_ades.py` 同步到服务器并通过 `py_compile`/`--help` 检查
- `heliolincrr/run_daily_pipeline.sh` 已更新：
  - survey 只跑 `RUN_DATE` 当天，用于生成今晚观测脚本
  - known/unknown 按 `RECOVERY_LOOKBACK_DAYS` 回看数据夜，适合断电后恢复
  - known daily 后等待 official known report 就绪，再运行 unknown
  - official known report ready 条件为 `<night>_matched_asteroids_ades.psv` + `<night>_mpc_reply.txt` 存在，或 status 明确 official matched `0`
- 自动 smoke 已完整跑通 `20260620`：
  - known official ADES/reply 已生成：`/processed1/20260620/L4/20260620_matched_asteroids_ades.psv`, `/processed1/20260620/L4/20260620_mpc_reply.txt`
  - unknown link `19`
  - GIF `19/19`
  - review package：`/pipeline/xiaoyunao/heliolincrr/review_packages/20260620/20260620_unknown_review_manifest.json`
  - review package manifest 中 `n_catalog_rows=19`, `n_gifs_copied=19`, `n_gifs_missing=0`
  - 网页已生成 `20260620_submit.csv`；19 条均为 `is_real=0`
  - 单夜 watcher 已自动写出 submit masked/ADES/stats，并以 `no_observations` 完成，未向 MPC 提交 unknown obsData
- 历史补报 watcher 仍在后台运行：
  - PID `2605614`
  - state `/pipeline/xiaoyunao/data/heliolincrr/review_submit_backlog_20251116_20260617.json`
  - 当前 summary：`complete=59`, `pending=63`, `failed=0`
  - `20260114` 的坏 submit CSV 已删除；网页随后重写了完整全假 submit，watcher 已处理为 `no_observations`
  - 不应与 daily 新夜 watcher 冲突，因为 state 文件和 night range 不同
- `watch_submit_reviews.py` 已避免对确定性坏输入反复 retry：
  - `blank/missing/unknown is_real`、重复 key 等错误会自动删除对应 `<night>_submit.csv`
  - 该夜回到 pending，等网页重新生成完整 submit CSV
  - submit CSV 或 review package mtime/size 变化后仍会重新处理

2026-06-20 状态：

- RA wrap 修复后的 full known/mask15 已完成并通过覆盖检查；unknown 后续应做 mask-only review rebuild，不应默认重画 GIF
- `unknown_remask_after_known_full_20260620_0346` 已停止；该流程错误地默认删除源 GIF 后重画，已完成到 `20251220`
- `remask_unknown_with_known.py` 已改为：
  - 默认保留 `/pipeline/xiaoyunao/heliolincrr/plots/<night>/unknown_link_*_<night>.gif`
  - 默认只清 `review_packages/<night>` 后重打 package
  - 当前 link 源 GIF 全存在时跳过 `plot_unknown_links.py`
  - `--skip-plots --require-existing-gifs` 可做严格 mask-only，缺图则不打不完整 package
- `20251220` smoke run 已验证：未重画 GIF，`unknown_count=20`, `n_gifs_missing=0`
- `20251221` 被旧残留 `plot_unknown_links.py` 清图后部分重画；已用 `plot_unknown_links.py --linkage-id` 只补缺失的 `8` 张 GIF，并重打 package，`n_gifs_missing=0`
- 严格 mask-only run 已完成：
  - command: `remask_unknown_with_known.py 20251222 20260617 --skip-plots --require-existing-gifs`
  - log: `/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_mask_only_strict_after_known_full_20260620_0945.log`
  - status: `/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_mask_only_strict_after_known_full_20260620_0945_status.tsv`
  - `done=91`, `skip=28`, `error=0`
  - done 夜 unknown link 合计 `3379`
  - 所有 done 夜均 `assigned_new=0`, `n_gifs_missing=0`
  - `20260103` 最终 mask 后 `unknown_count=1`
  - high-unknown skip: `20251226=366`, `20260111=378`, `20260528=926`, `20260611=653`
  - missing-input skip 共 `24` 夜，需后续分类确认
- ADES unknown 最终上报：
  - L2 catalog 已确认有 `Flux_Aper4`/`FluxErr_Aper4`
  - `logSNR` 采用 `log10(Flux_Aper4 / FluxErr_Aper4)`
  - `submit_reviewed_unknown.py` 默认透传 `--include-logsnr`
  - 如 MPC validation 临时不接受，可用 `submit_reviewed_unknown.py <night> --no-logsnr` 关闭
- Daily automation:
  - 服务器 known 匹配代码已确认包含 RA wrap 多中心查询；known Slurm 默认 `MASK_SEP_ARCSEC=1.5`
  - `known_asteroid/merge_night_parts.py` 会写 `/processed1/<night>/L4/<night>_known_asteroid_status.json`
  - `heliolincrr/summarize_single_night.py --require-mask15` 不再 fallback 到 1 角秒 matched；缺 mask15 时只允许 status 明确 empty mask
  - `heliolincrr/run_single_night.sh` 默认 `REQUIRE_KNOWN_MASK15=1`
  - `heliolincrr/run_daily_unknown.sh` 等 known mask15/status ready 后运行 unknown 并打 check 包；只有 check 包生成且 `n_catalog_rows > 0` 时才启动该夜 submit watcher
  - `heliolincrr/run_daily_pipeline.sh` 现在 survey 只跑 `RUN_DATE` 当天；known/unknown 按 `RECOVERY_LOOKBACK_DAYS` 回看实际数据夜补处理，适合 `0 9 * * *` 和 `@reboot` recovery
  - `heliolincrr/watch_submit_reviews.py` 扫描 `<night>_submit.csv`，要求已有 review manifest，用 state JSON 去重，并调用 `submit_reviewed_unknown.py` 导出/validate/submit
  - `heliolincrr/run_review_submit_backlog_watch.sh` 已用于历史补报区间 `20251116..20260617`
  - 已做 smoke：`20260103 --require-mask15` 得到 `unknown_count=1`；`20260605` known status 可表达 empty mask；watcher dry-run 可发现已有 submit CSV
  - `20260617` 真实 daily unknown 已完成：unknown link `0`，check package 已生成
  - 历史补报 watcher 正在后台运行：PID `2530505`，state `/pipeline/xiaoyunao/data/heliolincrr/review_submit_backlog_20251116_20260617.json`
  - 历史补报当前 summary：`review_packages=122`, `complete=42`, `pending=80`, `failed=0`
  - `20251116..20251231` 子区间 `37` 个 review package 已全部完成，真源 `20` 条 link / `60` 次探测
  - 该区间真源 GIF 和测光星表已导出到本地 `/Users/yunaoxiao/Desktop/true_unknown_20251116_20251231`
  - 已正式提交 reviewed unknown：`20251121` -> `2026-06-20T05:00:14.546_0000C0eY`; `20251124` -> `2026-06-20T05:02:22.435_0000C0ea`; `20251125` -> `2026-06-20T05:00:53.443_0000C0eZ`
  - cron 已在 `2026-06-21` 改为总入口 `run_daily_pipeline.sh`
- Known anomaly:
  - `20260111` known candidate 数正常（`all_asteroids=106042`），但 matched 极少（official `180`, mask15 `419`）
  - 相邻夜 matched 为几千到一万量级，且 `20260111` 候选星等分布并不更暗
  - 抽样显示 L2 catalog/WCS 相对 known ephemeris 有约 `8..11"` RA 方向系统偏移；Win/PSF 坐标一致
  - 初步判断不是 known 查询漏候选，而是该夜 L2 astrometry/WCS 有系统偏移；修复前不建议把该夜 known/unknown 作为可靠上报结果

回到 `heliolincrr` 主线，把单夜 unknown 搜索/提取自动化做稳，并在
`unknown_full_remask_20260514_171651` 全量重跑完成后做产物审计和 daily wrapper 收口。

最新全量重跑状态：

- 日志：`/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_full_remask_20260514_171651.log`
- 状态表：`/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_full_remask_20260514_171651_status.tsv`
- 范围：`20251119..20260514`
- 已完成时间：`2026-05-15T16:41:11+08:00`
- 状态表共 `177` 行：`done=115`, `skip=62`
- skip 原因：`no_l2=56`, `unknown_links_after_known_gt_200=3`, `no_mp_l2=2`, `missing_known_matched=1`
- 完成夜次累计 unknown link 数：`4764`
- 最后一行：`20260514 skip rc=0 no_l2`
- `2026-06-03` 复查：当前无 unknown 相关后台进程；combined link/detection 表仍与状态表一致

人工复核和上报状态：

- review package manifest/tar 共 `117` 个，包含早期 `20251116` 和 `20251118`
- review CSV 合计 `4848` 行，其中 `52` 行已填写 `is_real`，分布为 `1=4`, `0=48`
- 已填写 review 的夜次：`20260309`, `20260504`, `20260507`, `20260508`, `20260509`, `20260512`
- `export_unknown_ades.py` 已支持 `--review-csv`, `--require-review`, `--validate`, `--submit`
- `export_unknown_ades.py` 已支持网页最终提交文件 `--submit-csv`，默认读取 `/pipeline/xiaoyunao/heliolincrr/review_packages/<night>/<night>_submit.csv`
- `submit_reviewed_unknown.py` 已作为人工 check 完成后的单夜入口，默认生成 masked JSON/FITS、ADES PSV 和 stats，不自动 validate/submit
- `trkSub` 已从旧的 8 位 leading-zero managed base62 迁移为 MPC 要求的 7 位；服务器 history、unknown JSON/FITS/PSV、review/submit CSV、manifest/FITS、GIF 文件名、review tar.gz 和 combined 表均已迁移
- `run_single_night.sh` 默认关闭 unknown ADES export/validate/submit
- crontab 目前只有 known asteroid daily，没有 unknown daily wrapper
- 未发现 `*_unknown_mpc_reply.txt` 或 `*_unknown_validate_reply.txt`，人工筛选后的 unknown 上报尚无完成证据

`20260514` 之后新夜次处理：

- `unknown_after_20260514_20260603_162313` 已确认中断：状态表只有表头，日志停在 `20260528 mask_gaia`
- `unknown_after_20260514_20260617_161428` 已完成
- 日志：`/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_after_20260514_20260617_161428.log`
- 状态表：`/pipeline/xiaoyunao/data/heliolincrr/batch_logs/unknown_after_20260514_20260617_161428_status.tsv`
- 成功夜次：
  - `20260529`: `unknown=29`, `review_full_rows=87`, `ades_rows=87`, `n_gifs_missing=0`
  - `20260530`: `unknown=21`, `review_full_rows=63`, `ades_rows=63`, `n_gifs_missing=0`
  - `20260601`: `unknown=19`, `review_full_rows=57`, `ades_rows=57`, `n_gifs_missing=0`
- 跳过夜次：
  - `20260528`: `unknown_count=1168 > 200`
  - `20260611`: `unknown_count=653 > 200`
  - `20260605`: known-only 补跑无 matched detections，仍缺 matched FITS
- 本轮只生成结果/review package/ADES PSV；未 validate，未 submit

`2026-06-17` 复查补充：

- 当前有 unknown JSON/review manifest 的夜次共 `120` 个
- 其中 positive unknown 夜次 `114` 个，zero unknown 夜次 `6` 个，unknown link 合计 `4917`
- review CSV 共 `120` 个，`4917` 行；已填写 `110` 行，其中 `is_real=1` 为 `15` 条，`is_real=0` 为 `95` 条
- 2025 年子集共有 `37` 个夜次，其中 `35` 个 positive unknown、`2` 个 zero unknown，unknown link/review CSV 行合计 `1948`
- 2025 年 review CSV 尚未开始填写：`filled=0`, `is_real=1` 为 `0`, `is_real=0` 为 `0`
- 2025 年 review manifest、CSV、`*_unknown_review.tar.gz` 和 GIF 均完整，GIF 缺失合计 `0`
- `20260528` 高 unknown 主要来自异常密集 tracklet/link：`tracklets_total=92792`, `links_total=1170`, `unknown_count=1168`，集中在 group `47` 和 `OBJ_MP_1428_0236/0287/0309`
- `20260611` 高 unknown 同样来自密集场/文件组：`tracklets_total=267450`, `links_total=1319`, `unknown_count=653`；同时识别出 `all_same_asteroid=636` 个 known links
- `20260605` 未自动触发是因为缺少应在 `2026-06-06 09:00` 运行的 daily log；手动补跑后两个 MP 文件均无 matched detections，仍不能进入当前 unknown 流程
- 服务器 `2026-06-17 16:00:12` 重启，普通 cron 不会补跑当天 `09:00` 已错过的任务
- 已手动补跑 `2026-06-17` survey daily，生成并发布 `20260617_plan`
- 已手动补跑 `2026-06-17` known daily；target night `20260616` 因 `/processed1/20260616` 缺失而 skip
- 用户已确认 recovery/monitor 暂缓，待 unknown 人工 check 与上报链路稳定后再实现；当前重点是人工 check
- `20260512` 已存在网页生成的 `20260512_submit.csv`；dry run 成功生成 `1` 个 masked unknown link 和 `3` 行 ADES obsData，未 validate，未 submit
- `20260103` 的 6 条旧 unknown link 虽均被人工判为真源，但 JPL 临时检查确认其中至少 `5` 条为已知小行星；该夜应作为 known 提取漏检回归样本，不再作为 unknown submit 正例
- `2026-06-18` 迁移后 `trkSub` history 为 `4917` 条、全部 7 位，范围 `0000001..00001hj`；postcheck dry-run 无 8 位残留
- 早期 MPC unknown 正式 submit 为 `20260220`，submission ID `2026-05-14T03:00:19.564_0000Bx4c`；当时未经过人工 check，提交版本为 `34` 条 link、`102` 行 obsData。当前服务器 remask 后 `20260220` 版本为 `29` 条 link、`87` 行 obsData
- 8 位 trkSub 迁移备份 `/pipeline/xiaoyunao/data/heliolincrr/backups/trksub_8char_20260618_094807` 已按用户确认删除
- `2026-06-18` JPL 临时检查发现人工真源中至少 `6` 条实际为已知小行星；known 提取存在系统性漏检
- known 提取根因继续修复：L2 catalog BINTABLE 不再用 `NAXIS1/NAXIS2` 当图像尺寸；跨 RA=0 视场不能只改用 `RA=0` 查询中心，必须同时查原始 WCS 中心、`RA=0` 和 `RA=359` 后合并去重；official known matched 保持 `1.0"` 用于 ADES/MPC，上游另写 `*_matched_asteroids_mask15.fits` 用 `1.5"` 给 unknown 扣除
- 修复后临时验证：`3331 Kvistaberg`, `1522 Kokkola`, `2220 Hicks`, `168 Sibylla`, `1795 Woltjer`, `8219 (1996 JL)` 均可进入 known matched；临时验证目录已删除
- `2026-06-19` 再次用当前 `20260103` unknown 实时 RA/Dec 做 JPL/Horizons 复查，发现 remask 后新增的 `00001iF/00001iG/00001iH/00001iI` 仍分别是 `2125 Karl-Ontjes`, `2728 Yatskiv`, `3856 Lutskij`, `1989 Tatry`
- 已修复跨 RA=0 多中心查询并烟测 `OBJ_MP_1060_0055_cat.fits.gz`：上述四个对象均进入 all/matched/mask15，匹配到当前 unknown 的对应 `objID`
- `known_rematch_20260618_111935` 仍基于旧 RA=0 单中心逻辑生成，不能作为最终人工 check 输入；已废弃
- 已停止旧 GIF repair，并清空 review package/L4 中的 submit CSV 和 submit 派生产物
- 已启动新全量重跑：`RUN_ID=known_wrapfix_20260619_122524`
  - known/mask15 forced rematch：`20251116..20260617`
  - unknown remask：`20251116..20260617`
  - 参数：`FORCE_EXTRACT=1`, `MAX_PARALLEL=16`, `MASK_SEP_ARCSEC=1.5`, `MASK_MATCHED_SUFFIX=_mask15`
  - logs: `/pipeline/xiaoyunao/known_asteroid/runtime/logs/known_wrapfix_20260619_122524.log`, `/pipeline/xiaoyunao/known_asteroid/runtime/logs/known_wrapfix_20260619_122524_driver.log`
  - remask status: `/pipeline/xiaoyunao/data/heliolincrr/batch_logs/known_wrapfix_20260619_122524_unknown_remask_status.tsv`
- 本轮 remask 已改为每夜打包前默认清理旧错误 artifact：删除 `plots/<night>/unknown_link_*_<night>.gif` 和整个 `review_packages/<night>` 后再重画 GIF、重建 review package
- 批量 driver 已停用：旧 finalize 在 `FORCE_EXTRACT=1` 时会因为夜级 FITS 已存在而跳过 merge，导致 night-level known/mask15 仍是旧结果；已修复 `slurm_merge_submit.sh`，后续用逐夜手动流程推进
- 手动流程：每夜/每批提交 `ka_match` array，确认 array 离队后运行 `FORCE_MERGE=1 FORCE_EXTRACT=1 ./slurm_merge_submit.sh <night> false`，确认夜级 all/matched/mask15 时间戳为新 run，再进入下一夜
- 当前手动 merge 已完成到 `20251206`；`20251207..20260104` 是已提交但仍 active/pending 的 array，`20260105` 之后尚未继续提交
- 修复前生成的 known matched、unknown links 和 review packages 需要重跑后才能用于人工 check 后上报
- 已启动批量重跑：`RUN_ID=known_rematch_20260618_111935`
  - known forced rematch：中途因 Slurm `sbatch` 临时 `Resource temporarily unavailable` 停在 `20251228`；已给 submit/driver 脚本加 retry 并从 `20251228` 续跑
  - 最新进度口径改为三分类：有 `OBJ_MP` 输入的夜才作为应跑分母；无 `OBJ_MP` 或无 L2 的夜作为不适用/跳过
  - 当前应跑夜次 `128`，跳过/不适用 `22`；普通 known 和 mask15 本轮均已完成 `29` 夜，完成范围到 `20251217` 且包含 `20251221`
  - retry 后 `20251228` 已成功提交 array `190635` / finalize `190636`
  - `20260106` 曾因 `sbatch` stderr warning 污染 array job id，导致 finalize dependency 非法并停住；已修复 stdout/stderr 分离与 job id 解析，取消坏 array `191015`，并从 `20260106` 续跑成功，新 job 为 array `195437` / finalize `195438`
  - 后续 driver 因缺 L2 夜 `/processed1/20260203/L2` 退出；已改为缺 L2 时 skip 并继续，当前从 `20260205` 续跑
  - duplicate cleanup：并发残留导致的重复/坏 job `185697..185706` 已取消
  - `2026-06-19` known rematch 后的 unknown remask 已完成；旧人工 `<night>_submit.csv` 全部删除，review package 下当前 submit CSV 数量为 `0`，所有人工 check 需要基于新 review package 重做
  - unknown remask：已处理 `20251116..20260616`，不处理 `20260617` unknown
  - logs: `/pipeline/xiaoyunao/known_asteroid/runtime/logs/known_rematch_20260618_111935.log`, `/pipeline/xiaoyunao/known_asteroid/runtime/logs/known_rematch_20260618_111935_driver.log`
  - remask status: `/pipeline/xiaoyunao/data/heliolincrr/batch_logs/known_rematch_20260618_111935_unknown_remask_status.tsv`

短期目标：

- 暂停基于旧 unknown 产物的上报，等待 `known_rematch_20260618_111935` 的 known/remask 批处理完成
- 批处理完成后先抽查用户已人工检查的 `20260101..20260107` 和追加检查的 `20260110`，确认已知小行星不再进入 unknown
- 推进 unknown 人工 check 时，以修复后重建的 review package 为准
- 明确人工 check 完成判定：每个正 unknown 夜次存在 `<night>_submit.csv`，且每个 unknown `tracklet_id` 都有明确 `0/1`；zero unknown 夜次可直接视为完成
- 对已完成 check 的夜次用 `submit_reviewed_unknown.py <night>` 重新导出 filtered/masked 产物和 ADES，并先 validate，不自动 submit
- 抽查 review package tar 是否包含 GIF、review CSV、`*_unknown_review_full.fits`、`*_unknown_review_ades.fits`
- 对 skip 夜按原因分类记录，尤其是 `unknown_links_after_known_gt_200` 和 `missing_known_matched`
- 新增正式 daily unknown 入口，负责选择目标夜、防重复运行、日志记录和产物检查
- daily recovery/monitor 入口暂缓，等 unknown 人工 check 与上报链路稳定后再做
- 上报链路默认关闭；导出/提交必须显式打开
- GIF 可视化默认可跳过或限量，避免拖慢主计算

PPT 素材旁支已完成一批图件，统一在服务器
`/pipeline/xiaoyunao/ppt_assets_20260414` 管理，后续只按汇报版式做小幅微调。

## Milestones

1. `20260220` 已作为单夜 unknown 主测试夜完整跑通
2. `run_single_night.sh` 已串起 Gaia mask、tracklet、linear link、orbit confirm、summary、unknown GIF
3. `summarize_single_night.py` 已生成 unknown fit_ok catalog，并按 known asteroid match 结果扣除已知小行星
4. 已新增 `--trk-sub-map`、人工复核 CSV 过滤接口和 `export_unknown_ades.py` unknown ADES 导出器
5. 已实现全局 `trkSub` 历史分配：`/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl`
6. `20260220` unknown ADES 已通过 MPC test validation，并完成一次正式 submit
7. 0-link 夜次现在能生成 empty RR/orbit/summary/unknown JSON/FITS，并让 assign/plot 正常 0 行退出
8. 真实 0-link 夜 `20260128` 和 `20260224` 已重跑成功，均生成 unknown=0 的 summary 和 empty unknown FITS/JSON
9. `make_tracklet_linreproj.py` 已新增 `--max-tracklets-per-group`，默认 `100000`；`20260503` 重跑确认 dense group 被跳过
10. `20260503` 已重新生成 unknown GIF，并打包 review package，`n_gifs_missing=0`
11. `20251128`, `20260412`, `20260430` 已补齐空 unknown FITS/JSON
12. `mask_gaia.py` 已改为同时使用 `RA_Win/DEC_Win` 和 `RA_PSF/DEC_PSF` 匹配 Gaia，`run_single_night.sh` 使用 `--match-arcsec 1.5`
13. `package_unknown_review.py` 已新增 `*_unknown_review_full.fits` 和 `*_unknown_review_ades.fits`
14. `merge_tracklets_night.py --allow-empty` 和 `run_single_night.sh` 已支持 no-group/no-tracklet 夜生成 empty unknown 结果
15. `run_single_night.sh` 已支持 `MAX_UNKNOWN_LINKS_AFTER_KNOWN`，默认 `200`
16. `unknown_full_remask_20260514_171651` 已跑完整个区间，状态表统计为 `done=115`, `skip=62`
17. PPT 素材已新增已知小行星 detection histogram 和 orbit confirm diagnostic card 脚本
18. 已按人工复核包 `*_unknown_review_full.fits` 合并出全量 unknown link 表，桌面副本在 `/Users/island/Desktop/unknown_review_combined_20260516/`

## Outstanding issues

- 全量重跑完成后只做了状态表、combined 表和 manifest required outputs 的初步审计；还未逐夜深审 GIF/ADES/review full 内容
- `20260414` 目录名与 observing-night 归属不一致；本轮因 `missing_known_matched` 跳过
- `unknown_links_after_known_gt_200` 的 3 个 skip 夜需要记录具体 night 并决定是否人工排除或调参重跑
- 当前单夜自动化还没有“每日选择目标夜 + 防重复 + 日志 + 产物检查”的外层入口
- 当前没有“错过 09:00 cron 后自动补跑/报警”的 recovery 机制；服务器重启后需人工确认当天 survey/known/unknown 任务是否已执行
- GIF 可视化很慢，应在自动提取中默认可跳过或限量
- 15 夜流程尚未接入同一个 `trkSub` history
- unknown 真实提交策略仍需收紧；后续 daily wrapper 默认不应自动 submit，人工筛选后需先 validate 再显式 submit
- 人工复核目前只填写少量条目；修复 known 提取并重建对应 unknown/review 产物前，已筛选的 `is_real=1` 不应导出为正式 unknown ADES
- 修复前的 known 提取可能在 L2 BINTABLE 尺寸、RA=0 跨界、`1.0"` 半径三方面漏检已知小行星，需要确定重跑范围
- `20260528` 和 `20260611` unknown 数过高，需要单独复盘
- `20260605` known-only 补跑无 matched detections，是否允许 no matched 夜继续 unknown 需要决定
- PPT 侧 `known_object_detection_histogram` 已可出图，但 summary JSON 是否稳定记录该字段仍需补查

## Validation criteria

- 给定一个存在 `/processed1/<night>/L2` 和 known asteroid match 的夜次，自动入口能完整运行或明确跳过
- 成功夜必须存在，即使 unknown=0 也要写出 schema-only FITS/JSON：
  - `/pipeline/xiaoyunao/data/heliolincrr/<night>/analysis/<night>_single_night_summary.json`
  - `/pipeline/xiaoyunao/data/heliolincrr/<night>/analysis/<night>_single_night_summary.txt`
  - `/processed1/<night>/L4/<night>_unknown_links.json`
- summary 中 `counts.matched_detections_total > 0`，否则不能认为已完成已知小行星扣除
- 若启用 `trkSub`，每个 unknown link 的 `trk_sub` 必须符合 MPC 当前要求：不超过 7 个字符；自动分配器使用固定 7 位 `[0-9a-zA-Z]`
- `trkSub` 分配必须写入全局 history，并且重复运行同一 catalog 不新增重复记录
- 若启用人工复核或网页 submit CSV，`is_real=0` 的 `trk_sub` 不得进入 ADES PSV
- 网页 submit CSV 模式必须在导出前检查所有 unknown `tracklet_id` 都有明确 `0/1`
- 自动入口日志需记录 target night、skip/run reason、exit code、核心产物路径和 unknown count
- PPT 已知小行星图可稳定写出 `known_object_detection_histogram.png`
- PPT orbit confirm 图可稳定写出 `20260220_link234_orbit_fit_diagnostic*.png`
- 全量人工复核包合并表应满足 link-level 行数 `4764`、detection-level 行数 `14325`、逐夜 link 数与状态表无差异

## Next recommended steps

1. 明早检查 `cron_daily_pipeline.log`、`<run_date>_daily_pipeline.log` 和 `<run_date>_unknown_daily.log`，确认 09:00 cron 自然触发且不是人工 smoke
2. 继续让历史补报 watcher 等待 `20251116..20260617` 剩余 submit CSV，阶段性检查 `complete/pending/failed`
3. 遇到网页生成的不完整 submit CSV，确认 watcher 删除该文件并让该夜回到 pending
4. 抽查修复后 `20260102/20260103` 中 JPL 命中的已知小行星是否从 unknown 中消失
5. 以修复后产物继续网页筛选；旧 submit CSV 只能作为人工判断参考，不能直接上报
6. 对已完成 submit CSV 的夜次用 watcher/`submit_reviewed_unknown.py` 生成筛选后 PSV 并正式 submit；异常时先看 validate/submit reply
7. 对 `20260528`、`20260611` 单独复盘 high unknown 原因，重点查 known subtraction、mask 后源密度、dense group 和观测场区
8. 复盘 `20260111` L2 astrometry/WCS 约 `8..11"` RA 系统偏移；修复前不要把该夜 known/unknown 当可靠上报
9. 如 PPT 还要继续打磨，再补 detection histogram 的 cumulative 版或 orbit confirm 的 good/rejected 双栏对照图
