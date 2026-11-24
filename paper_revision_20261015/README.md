# SHARP 启动期论文分析交付

兴隆 60/90 cm Schmidt 的宽黄道带重复成像，用实拍几何覆盖、场夜重访、测量质量与单夜移动源关联来量化巡天的启动期表现。

正文继续采用 **2025-11-15 至 2026-07-15** 的原八月冻结期间，并在图表内注明各自的质量和处理条件。历史代码阅读入口为 `2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757`，旧产物的 August provenance 将生产算法基线记为 `d2f0057bd4cdd4c7a3c6c2431f6449a72c327284`。后续区间仅按可信原图曝光中点盘点，边界为北京时间 2026-10-16 00:00，保留了边界当夜的截止前后计数。

重建了 **39,590 张合格曝光和 129 夜**。实拍 WCS 覆盖约 **10,928 deg²**，曾达到同夜三访的面积约 **10,471 deg²**。实拍相邻重访中位 37.65 分钟，目录模型深度中位 20.25 mag，原冻结条件关联保留率 74.16%。

## 给正文的文件

- [MANUSCRIPT_FACTS.md](MANUSCRIPT_FACTS.md) 提供可引用的数字、方法事实、图注要点和来源。
- [FIGURE_CHANGELOG.md](FIGURE_CHANGELOG.md) 对应旧稿十二图与新版八幅主图、两幅附图。
- [AUTHOR_TODO.md](AUTHOR_TODO.md) 列出影响论文结论的待补证事项。
- [tables/baseline_reproduction.csv](tables/baseline_reproduction.csv) 核对旧快照计数。[tables/observing_summary.json](tables/observing_summary.json) 和 [tables/performance_summary.json](tables/performance_summary.json) 给出本次计算的统计。
- [tables/exposure_ledger.csv.gz](tables/exposure_ledger.csv.gz) 包含原曝光、夜次、质量、L1/L2 可用状态和冻结 L4 product 状态。L4 的夜级存在字段与逐曝光预测行数单独记录。
- [tables/extension_by_night.csv](tables/extension_by_night.csv) 汇总七月后档案。后续存在不同过滤器和处理配置，边界夜的夜级 L4 包含截止后的曝光，故独立交付。
- [freeze_manifest.json](freeze_manifest.json) 列出输入、分析环境和 SHA-256。

## 本地复现

在 ZIP 解压出的 `paper_revision_20261015` 目录中，用 Python 3.12 运行以下命令。绘图数据、观测归档、已筛选的星表列、冻结星历和地球定向表均在包内。图 1 和图 8 为附带的已验收原版成品。图 1 的原全幅 L1 约 278 MB，位置与校验值在源清单中。

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/exposures.py
python scripts/known_performance.py
python scripts/stellar_repeatability.py
python scripts/make_figures.py
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
