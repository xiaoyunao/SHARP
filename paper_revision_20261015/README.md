# SHARP 启动期论文分析交付

兴隆 60/90 cm Schmidt 的宽黄道带重复成像，用实拍几何覆盖、场夜重访、测量质量与单夜移动源关联来量化巡天的启动期表现。

正文继续采用 **2025-11-15 至 2026-07-15** 的原八月冻结期间，并在图表内注明各自的质量和处理条件。历史代码阅读入口为 `2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757`，旧产物的 August provenance 将生产算法基线记为 `d2f0057bd4cdd4c7a3c6c2431f6449a72c327284`。后续区间仅按可信原图曝光中点盘点，边界为北京时间 2026-10-16 00:00，保留了边界当夜的截止前后计数。

重建了 **39,590 张合格曝光和 129 夜**。实拍 WCS 覆盖约 **10,928 deg²**，曾达到同夜三访的面积约 **10,471 deg²**。实拍相邻重访中位 37.65 分钟，目录模型深度中位 20.25 mag，原冻结条件关联保留率 74.16%。

## 2026-11-24 R2 稿件整合（2026-11-25 完成）

完整 [LaTeX 源码与编译 PDF](manuscript/)、[三幅新版图及紧凑复现数据](figure_revision/)、[定标方法核查](inputs/calibration_method_note_20261124_r2.md) 和 [给 GPT Pro 的交接](REVISION_R2_HANDOFF.md) 已直接纳入本分支。保持 8 主图+2 附图。Gaia WCS 使用传播坐标；W 测光匹配实际读取原始坐标，这一区别已写入正文。

```bash
python figure_revision/redraw_figures.py --prepare --source .
python scripts/finalize_revision_r2.py --source .
cd manuscript
tectonic --keep-logs manuscript.tex
# 或 latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript.tex
```

`--prepare` 只重建三图的数据/绘图，不运行生产程序；省略它可直接用紧凑图数据重绘。当前图指纹与数字来源使用 `scripts/finalize_revision_r2.py` 更新。旧 `build_manuscript_source_map.py` 及 `verify_revision_geometry.py` 对应前一稿的历史 SHA，仅用于旧版本核对。

## 2026-11-24 二次核验

[核验摘要](REVISION_CHECK_SUMMARY_20261124.md) 包含合格曝光条件连接率、58 条候选的逐成员对照、模型深度及远距目标复核。[稿件数字来源](tables/manuscript_number_source_map.csv) 记录正文与图注的样本定义。[L1 局部切片](audit_figures/distant_l1_local_evidence.pdf) 为远距目标核验材料。

可重跑新增统计及切片图。L2 的局部搜寻结果、L1 像素切片和冻结产品状态均附在仓库中。

```bash
python scripts/reconcile_historical_cohorts.py
python scripts/review_distant_associations.py
python scripts/render_distant_stamps.py
```

仅对上一轮历史稿，可在论文目录中运行 `python scripts/build_manuscript_source_map.py <源工程根路径>` 和 `python scripts/verify_revision_geometry.py <源工程根路径>`。两项校验使用源工程中的 LaTeX、PDF 和图指纹。

## 给正文的文件

- [figures/](figures/) 提供可直接使用的正文图和附图，每幅含 PDF 和 300 dpi PNG（R2 新绘三图 320 dpi）。[figure_data/](figure_data/) 提供绘图数据。
- [tables/main_results.csv](tables/main_results.csv) 为正文主结果表。
- [MANUSCRIPT_FACTS.md](MANUSCRIPT_FACTS.md) 提供可引用的数字、方法事实、图注要点和来源。
- [FIGURE_CHANGELOG.md](FIGURE_CHANGELOG.md) 对应旧稿十二图与新版八幅主图、两幅附图。
- [AUTHOR_TODO.md](AUTHOR_TODO.md) 列出影响论文结论的待补证事项。
- [tables/baseline_reproduction.csv](tables/baseline_reproduction.csv) 核对旧快照计数。[tables/observing_summary.json](tables/observing_summary.json) 和 [tables/performance_summary.json](tables/performance_summary.json) 给出本次计算的统计。
- [tables/exposure_ledger.csv.gz](tables/exposure_ledger.csv.gz) 包含原曝光、夜次、质量、L1/L2 可用状态和冻结 L4 product 状态。L4 的夜级存在字段与逐曝光预测行数单独记录。
- [tables/extension_by_night.csv](tables/extension_by_night.csv) 汇总七月后档案。后续存在不同过滤器和处理配置，边界夜的夜级 L4 包含截止后的曝光，故独立交付。
- [freeze_manifest.json](freeze_manifest.json) 列出输入、分析环境和 SHA-256。

## 本地复现

从 GitHub 检出论文分支后，在 `paper_revision_20261015` 目录中用 Python 3.12 运行以下命令。绘图数据、观测归档、已筛选的星表列和地球定向表在仓库中。DE432s 星历由 JPL 发布，下载后核对 SHA-256。[JPL NAIF 星历目录](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/)图 1 和图 8 为附带的已验收原版成品。图 1 的原全幅 L1 约 278 MB，位置与校验值在源清单中。

```bash
curl --fail --location --output inputs/de432s.bsp https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp
shasum -a 256 -c <<'EOF'
363f32e14f5255359ac32c4d38080cf28ab55564a5e16696a75f63394b666e9b  inputs/de432s.bsp
EOF
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/exposures.py
python scripts/known_performance.py
python scripts/stellar_repeatability.py
python scripts/make_figures.py
python scripts/finalize_revision_r2.py
python scripts/verify_bundle.py
python scripts/freeze_delivery.py
```

可从冻结匹配及链接记录再运行一次条件关联分析。输出中有对象夜逐行数据、身份映射与各箱的 Wilson 区间。输出目录须为新路径。

```bash
python scripts/replay_frozen_recovery.py \
  --products-dir inputs \
  --l2-manifest inputs/l2_manifest.csv.gz \
  --output-dir replayed_recovery
```

每夜 ≥3 次的几何面积先对当夜各曝光的足迹叠加，再取所有夜的达标天区并集。没有坏区、遮挡与检出概率的全量空间掩膜。点源深度为该星表记录的 5σ 光度模型估计。关联曲线复现原回顾性已检出样本，按天体夜计数，和 1″ 已知产物的样本定义分别记录。
