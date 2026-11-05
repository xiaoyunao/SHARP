# SHARP 实际计划绘图工具

这些独立脚本不修改 SHARP、不部署、不联网，不需要重跑 200 条天气实现。

## 1. 从原始年度计划导出小型 bundle

在保存了 `targeted_observing/products/annual/selected_nominal` 的机器上运行：

```bash
python3 export_plot_data.py /path/to/SHARP \
  --output /path/to/SHARP/targeted_observing/docs/annual_reports/plot_data \
  --start-date 2026-11-01 \
  --stop-date-exclusive 2027-11-01 \
  --source-commit a5cb9454e56771ebe4211b7a544412056cce4c40
```

输入也可以直接指向 `selected_nominal`。默认要求 365 个连续夜次，包括零计划夜。其他完整区间可显式设置 `--expected-nights`；单日烟测用 `--expected-nights 1`。已有输出默认不覆盖，必要时传 `--overwrite`。

导出 `fields.csv`、`exposures.csv.gz`、`night_index.json`、`manifest.json`。这些小文件可与报告一同提交 Git，原始星历缓存与天气实现不必提交。CSV 的 `*_json` 列保存原始列表/字典；缺必需字段会报出文件与 JSON 路径。

## 2. 用 bundle 画实际指向图

依赖现有 Python 环境中的 `numpy`、`matplotlib`、`astropy`。脚本本身不会安装或下载依赖。

```bash
python3 plot_exported_plan.py /path/to/plot_data \
  --output /path/to/annual_plan_figures
```

输出：

- `pointings_aitoff.png`：真实计划中的 EV / AM / Trojan 指向，以不同颜色显示球面视场边框；同一 field/version 的重复访问仅画一遍。
- `monthly_pointings_YYYY-MM.png`：每月三个模式分别裁剪的局部 RA/Dec 图。每个面板的坐标范围独立，RA 向左增加。
- `plot_provenance.json`：源 bundle 哈希、源码版本、计划状态、字段/目标计数及绘图范围。

视场宽度来自 manifest 配置，位置角来自字段本身，边框在球面 sky-offset 坐标系计算；不会用 `RA±1.5°` 代替真实球面边框。SHARP 当前产品 PA 为 0。小图中的短标签是 field ID 的尾部，不是目标位置。

**计数必须分开：**field/version 数、独立指向几何数、独立计划目标 ID 数和曝光数分别保存。字段可能互相重叠，不能把字段数乘 9 当作独立面积；目标 ID 来自主/伴随成员列表，不能当作成功光变或新周期数。没有目标星历时不把目标画在视场中心，不生成虚构源位置。

## 3. 画指定夜次的实际曝光高度与时间表

```bash
python3 plot_exported_plan.py /path/to/plot_data \
  --output /path/to/annual_plan_figures \
  --night 2026-12-15 \
  --night 2027-04-15 \
  --iers /path/to/SHARP/targeted_observing/cache/iers_20261104.all
```

`--night` 可重复；日期必须在真实 bundle 中。一次批量生成全部夜图可用：

```bash
python3 plot_exported_plan.py /path/to/plot_data \
  --output /path/to/annual_plan_figures \
  --all-nights \
  --iers /path/to/SHARP/targeted_observing/cache/iers_20261104.all
```

`--all-nights` 按真实 `night_index.json` 循环，不补造日期；零计划夜仍生成太阳/月亮与明确的空计划/跳过原因。不能与 `--night` 同时使用。脚本会在绘图前核对本地 IERS 覆盖所有请求夜次。加 `--no-monthly` 可跳过月图。

`night_YYYY-MM-DD.png` 上半图显示太阳/月亮、太阳 −12°/−18°、源配置的最低/中心高度参考线（当前为 30°/32.5°，未另加采样余量），以及每帧计划曝光中点的指向高度；字段淡线只画在原计划 held block 内。下半图按块/场显示逐曝光彩条、读出及转场/换滤镜开销，横轴为北京时间，午夜后属于次日。

台址必须来自 `manifest.configurations[config_id].site`，或同时显式指定：

```bash
--site-lat 40.393 --site-lon 117.575 --site-height 960
```

不会静默假设台址。夜图必须提供覆盖该夜的本地 IERS-A 文件，自动下载始终关闭；高度为无大气折射的几何高度，不是望远镜执行或测光成功认证。源计划的空闲时间原样显示，不由脚本重新分配。

如运行环境将依赖安装在额外目录，可给绘图脚本传 `--deps /path/to/packages`。

## 4. 特殊目标代表夜数据（与实际年度计划图分开）

`compute_special_visibility.py` 可对已有 Horizons 索引计算每月代表夜及两晚的高度轨迹。它不生成原年度计划，只用于候选对象的几何比较：

```bash
python3 compute_special_visibility.py \
  --snapshot /path/to/SHARP \
  --index /path/to/special_ephemeris_index.json \
  --iers /path/to/SHARP/targeted_observing/cache/iers.all \
  --output /path/to/special_plot_data
```

输出的 provenance 明确区分代表夜与全年排程，V 与 r、高度可见与高精度可测、可行窗口与获分配时间。

## 本次烟测范围

只用已有重建代表夜做一日烟测：16 场记录、80 帧；没有声称重建了原始 365 夜产品。样例图显著标注 `RECONSTRUCTED TEST INPUT - not the original annual schedule`。正式图需在持有真实年度产品的机器上先导出 bundle。缺 bundle、缺字段、哈希不一致、缺指定日期或缺台址/IERS 时明确报错，不补造坐标。
