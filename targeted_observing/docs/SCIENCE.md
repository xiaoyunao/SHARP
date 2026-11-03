# 科学来源与选场解释（2026-11-03）

## 已核查的依据

- [Bolin et al. 2024，Palomar twilight survey](https://arxiv.org/abs/2409.15263)：支持r短曝光、同窗多次访问与恢复链；不把ZTF发现率缩放成本台预测。
- [ZTF官方Trojan实验](https://www.ztf.caltech.edu/ztf-experiments.html)：已有g/r、180秒、约小时采样实验；本项目应强调少数目标分钟级连续光变和多夜一致性。
- [Li et al. 2025](https://arxiv.org/abs/2504.15530)，[发表记录](https://authors.library.caltech.edu/records/8ecpv-f3d69)，[作者数据表](https://github.com/ZhuofuLi/ZTF-Jupiter-Trojan-Rotation-Periods)：读取216对象周期与可靠性列；低质量差异只能列为待检验候选。
- [Ďurech & Hanuš 2026](https://arxiv.org/abs/2610.06082)，[Zenodo结果表](https://doi.org/10.5281/zenodo.22812438)：下载完整/partial自转表，交叉247个Trojan背景；这改变了“无已知周期”的判断。
- [DAMIT数据导出](https://damit.cuni.cz/projects/damit/exports)：交叉asteroids与asteroid_models，以正式number连接；模型版本可能更新，与新论文不是独立样本。
- [LCDB官网](https://minorplanet.info/php/lcdb.php)：公开版2023-10-01。读取summary固定列、保留U；目录未命中不等于周期未知。
- [Pál et al. 2020 TSSYS-DR1](https://arxiv.org/abs/2001.05822)，[作者发布表](https://archive.konkoly.hu/pub/tssys/dr1/)：取得release.merge及列说明，检查本次全部60个候选。
- [TESS 2025重分析](https://arxiv.org/abs/2501.07189)，[CDS目录](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A%2BA/693/A66?format=html&tex=true)：文章可定位，FTP主站和HTTPS镜像下载分别SSL EOF/502失败；**未完成逐对象重分析表匹配**，DR1不能替代这份表。
- [TESS EDR2](https://arxiv.org/abs/2503.13332)：已获取论文，44目标的特定样本不代表全部TESS覆盖。
- [Mottola et al. 2011](https://doi.org/10.1088/0004-6256/141/5/170)：直接检查5476的5.78小时和5511的长周期讨论，后者依赖跨夜比较星一致测光。
- [SBDB Query文档](https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html)及[Horizons API文档](https://ssd-api.jpl.nasa.gov/doc/horizons.html)：TJN目录用于初筛，实际目标轨迹使用站心观测星历，RA/Dec为ICRF天测坐标、时间UTC。

## 当前两场的科学用途

| 核心对象 | 公共周期证据 | 本轮用途 |
|---|---|---|
| 5476 Mulius | LCDB U=3，5.78 h；2026模型5.778835 h | 可靠周期/时间系统控制 |
| 6002 Eetion | LCDB U=3−，12.918 h；新模型12.800626 h | 第二个交叉核对对象；约1%差别不能称巨大冲突 |
| 4754 Panthoos | LCDB U=2+，27.68 h；新模型66.859319 h | 检查旧解与新解预测；新工作可能已解决旧歧义，不能预先宣称未解问题 |
| 5511 Cloanthus | LCDB U=2−，336 h；新模型322.773350 h | 慢转/跨夜零点试验，14天不足保证唯一周期 |
| 18054 (1999 SW7) | 本次五类表无可用周期记录；定向网页搜索未找到一手周期结果 | **未核查充分**的候选，不标“首次测周期” |

上述5个对象在TSSYS-DR1均未命中；不推断其他TESS发布没有结果。主场4个核心、邻场1个核心，V约17.45–18.32，V仅为亮度代理，不是r极限或已测r星等。样本有两个控制对象而不是强行满足未知占比；有新数据后可按科学价值调整。

两场位于低银纬（约−8°），因此明亮目标数不是唯一评价。Gaia亮星检查与轨迹最近距离写入独立污染报告；完整深星密度、CCD坏区/鬼像/光晕和目标掠星测光需要首轮图像。不能声称已得到无污染样本。

## 对生产处理的约束

`heliolincrr/run_single_night.sh` 当前短轨迹参数 `vmin=3`、`vmax=63`，单位角秒/小时。本轮只读本地源码，不改变搜索范围或形态筛选；拖线可能被点源形态门控拒绝。需要未来用注入源/已知高速NEO验证速度、拖线长度、测光及恢复延迟的选择函数，才能量化晨昏发现能力。

数据中的月距阈值、V代理门槛和目标组成是初始工程/科学选择；r30/r120的测光精度和深度只能由实际观测给出。
