# 三幅论文式配图及中文图注

## 图 1　SHARP 实拍天区与观测量

2025 年 11 月 16 日夜至 2026 年 11 月 23 日夜，SHARP 在 176 个有效观测夜取得 55,605 张合格曝光，快门累计开放 490.5 小时，实拍几何覆盖约 11,199 平方度。图 (a) 使用 ICRS Mollweide 投影，颜色为覆盖次数，橙线为黄道。图 (b) 为各月合格曝光数，末月截至 11 月 23 日夜。

图源 **SHARP 巡天项目**。实际曝光来自兴隆原始头与 L2 WCS，年度面积使用 NSIDE 1024，地图显示使用同一曝光样本的 NSIDE 512 网格。 [PNG](fig01_annual_coverage.png) · [PDF](fig01_annual_coverage.pdf) · [年度统计](annual_summary.csv)

## 图 2　小行星 (282917) 2007 LA33 的三次实拍

背景恒星在三幅图中保持固定，黄色箭头标出小行星 (282917) 2007 LA33 的位置。2026 年 9 月 1 日三幅图的首末曝光中点相隔约 61 分钟，目标相对背景恒星移动约 29 角秒。

图源 **中国科学院国家天文台兴隆观测基地 / SHARP 巡天项目**。三帧源为 `/processed1/20260901/L1/OBJ_MP_1849_0201.fits.gz`、`OBJ_MP_1849_0223.fits.gz` 和 `OBJ_MP_1849_0245.fits.gz`，均为 W 波段，约 30 秒。UTC 中点分别为 **16:14:02.102、16:44:42.572、17:15:18.395**。按实测 WCS 重投影为北上东左、1.15 角秒每像素的同一网格。用各帧背景中值和噪声尺度归一，再使用同一灰度拉伸。静态图和 GIF 采用同三帧。MPC 用内部轨迹 00001mX 逐次匹配到该已知身份。 [PNG](fig02_asteroid_motion.png) · [PDF](fig02_asteroid_motion.pdf) · [GIF](fig02_asteroid_motion.gif) · [MPC 核对表](candidate_review_crosswalk.csv)

## 图 3　巡天检出的已知小行星轨道分布

图 (a) 为半长轴与偏心率，图 (b) 为半长轴与轨道倾角。颜色和点型区分七个轨道类别。两图横轴和右图倾角轴采用对数刻度。左图虚线标出近日点距 1.3 天文单位的边界。

图源 **SHARP 巡天项目**。每个点对应一颗已知小行星，两面板各绘出同一批 81,909 颗。轨道要素来自本次冻结的 Lowell ASTORB 快照，按轨道要素分组。图例的特洛伊区域和远距类别采用此分组定义。 [PNG](fig03_orbital_distribution.png) · [PDF](fig03_orbital_distribution.pdf) · [逐对象轨道表](known_object_summary.csv.gz)
