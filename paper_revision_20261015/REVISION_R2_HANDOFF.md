# SHARP R2 handoff — 2026-11-25

起点：科学输入/上一轮核验提交 `2f5f71d5a074c4cfd712b0e3fbef240c40745e2a`，GPT Pro 2026-11-24 R2 全文工程。按作者最新指示只核读当前上游实现，不再追查更早备份。工作和交付均在论文分支；服务器全程只读。

## 已完成的稿件与新图

- `manuscript/manuscript.tex`、各章节、文献、10 幅图和本地重新编译的 `manuscript.pdf` 直接纳入仓库。17 页，8 正文图+2 附图。Survey/Performance/Initial Results 为 2/3/2 个小节，Methods/Discussion 连续段落，Transient Sources 文件只有标题和 label。
- Fig.2、5、6 已由随稿脚本从原冻结 parquet/ledger 重新生成紧凑数据再绘图。最终 PDF/320 dpi PNG 同步到 `figures/`；论文只引用同一份 PDF 内容。其余七幅图保持旧 SHA。
- 实际选样：Fig.5b 500,092；Fig.5c 531,416/531,612，星等范围外 196；Fig.5d 全部 531,612/128 夜。Fig.6 全局及速度 75,899/102,347，星等相关 75,897/102,315；排除的 32 中 2 个保留；比例网格 N<30 遮罩。
- `scripts/make_figures.py` 接入 R2 绘图，避免再次生成旧版三图；`scripts/finalize_revision_r2.py` 更新图指纹及本轮数字来源，校验筛选摘要、输入哈希、结构和图件一致性。

## 已修改、供 GPT Pro 直接采用的方法句子

Section 2.1 已补写：

> Reference ICRS positions at catalog epoch J2016.0 are propagated to the Julian year of the UTC exposure start using the catalog proper motions and Astropy's apply_space_motion, with the right-ascension proper motion interpreted as mu_alpha cos(delta). Parallax and radial velocity are not supplied; sources with nonfinite propagated coordinates retain their catalog positions. The propagated coordinates form the reference catalog consumed by SCAMP, whose fitted TPV solution supplies the final image WCS.

并补入以下实际调用差别：

> The photometric reference reader also computes propagated coordinates, but the white-light calibration matches stars using their original catalog coordinates.

这是本轮最需要保留的方法细节：`cal_pos.get_gaiacat` -> X_WORLD/Y_WORLD -> SCAMP 的链路成立；`calib_flux.gaia_mags` 实际以原始 ra/dec 建立测光匹配坐标。不要把“WCS 传播”概括为“所有参考星匹配都传播”。不修改程序、不重算零点。

三张冻结曝光对应 2025-11-21、2026-02-18、2026-05-08，均与冻结台账尺寸及 mtime 一致；每张 4 个高自行参考星+2 个缺失自行回退星，共 18 行明细。使用传播坐标时的亮参考星位置残差中位数为 0.1151/0.0755/0.0741 arcsec，原坐标为 0.1903/0.1748/0.1823 arcsec；此为局部链路印证，不增加新的巡天性能样本，也未重做 WCS。

调用链、源代码/瓦片 SHA、配置、曝光时间、坐标明细、两段可直接采用的英文和证据边界见 `inputs/calibration_method_note_20261124_r2.md`、`inputs/calibration_chain_evidence_r2.json` 与 `tables/calibration_*_r2.csv`。Gaia EDR3 依据作者确认，瓦片本身的精简头未带发行版声明；不再作为作者待答问题。

## 图注与数字同步

R2 图注定义均保留；本地字体可用性使 PDF 字节和字体选择与 GPT Pro 包不同，新 SHA 已据实更新。Fig.5c 是径向 O−C 密度；signed 夜偏移只在 d。Fig.6 的计数与比例同网格且未平滑，1D 使用原分箱及 95% Wilson 区间。`figure_revision/pdf_validation.json` 更新为实际生成文件的尺寸和嵌入字体。

摘要继续引用合格曝光、物理夜重建且沿用历史状态政策的 75,467/101,784=74.14%。正文另列历史 74.16% 和进一步限制为无隔离且阶段归档可读子样本的 73,787/99,426=74.21%。58/179 全部通过；4,762 自动 links 保持历史处理产物定义。远距 79/9/31 是关联级分类；身份 7/7/3 对应组总关联 107/9/3；31 未决中 28 属于已有其他运动证据的身份。深度 20.25 mag、18 px 参考零点、FWHM/2 噪声面积、无额外包围能量补偿均不变。

## 验收与后续材料

- 三图原冻结输入 SHA、绘图计数和一维分母核对通过；10 幅论文/支撑图 PDF 一致，嵌入字体通过；7 个小节和空暂现源章节自动核对通过。
- Tectonic 0.16.9 本地重新编译成功。无未定义引用、无 overfull/underfull box；仅上游 AASTeX 类名提示和 Tectonic 读取旧 lineno 包注释时的 UTF-8 警告。17 页渲染后检查，新图及空章节正常。
- 不再重派上一轮五项核验，不扩展 52 条候选身份研究。
- 后续按作者/同事材料合并系统响应曲线、进一步仪器资料、完整署名致谢/in-prep 书目信息与暂现源正文。本轮无需作者先答新问题。
