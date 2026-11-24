# SHARP 一周年通讯稿交付

- [事实表](NEWS_FACTS.md)、[Word 替换段落](WORD_REPLACEMENTS.md)、[中文图注和署名](FIGURE_CAPTIONS.md)
- [年度汇总](annual_summary.csv)、[逐月汇总](monthly_summary.csv)、[候选阶段与身份汇总](candidate_status_summary.csv)
- [观测天区和月度进度](fig01_annual_coverage.png)、[小行星移动三帧](fig02_asteroid_motion.png)、[轨道分布](fig03_orbital_distribution.png)、[循环动图](fig02_asteroid_motion.gif)
- [逐曝台账](exposure_ledger.csv.gz)、[已知对象及轨道要素](known_object_summary.csv.gz)、[逐候选身份核验](candidate_review_crosswalk.csv)、[覆盖分辨率对照](coverage_summary.csv)

图 1 和图 2 的字体、白底及配色参考 Pal 5 New 项目。图 3 按作者提供的 SHARP 论文图样绘制，使用衬线字体、用不同的色彩及点型区分七类轨道、用对数轴展示半长轴和倾角范围。三幅图按论文版式绘制，使用白底、统一刻度和紧凑轴标，均有不少于 3000 像素长边的 PNG 及矢量文字 PDF。图内仅留必要英文、坐标和数字。三帧动图与静态版采用相同图像和时刻。

## 复现

在 SHARP 仓库根目录，用具有 NumPy、Pandas、Astropy、Healpy、SciPy、Matplotlib、Pillow 的 Python 运行。交付目录附带覆盖网格、选取的切图和轨道结果表，可据此重绘三图。

```bash
news=outputs/news_release_20261124
python3 "$news/scripts/plot_news.py"
python3 "$news/scripts/plot_news.py" --orbit-only
```

使用完整的兴隆与论文快照可重算逐曝、已知和候选统计。汇总另外需要 PyArrow。完整输入快照保存在原工作区 `evidence/`，波段切图和 HEALPix 网格随交付分支提供。原始大 FITS 继续保存在兴隆。

```bash
python3 "$news/scripts/summarize.py"
python3 "$news/scripts/summarize_candidates.py"
```

`paper_revision_20261015/` 取自提交 `e71f9f6a2ea95205e0a30230039ea65b50949914`，并用 `2f5f71d5` 复核冻结候选归组。`nov09/` 保留早期后审记录。`server_*.jsonl.gz` 保存 L2 足迹、原始头、单夜处理清单、known 选列、候选目录及复核状态。名义排程只作为逐帧任务标签的证据。`wamo.json`、`jpl_*.json` 和 `astorb_elements.csv.gz` 保存身份查询与轨道分组依据。

如需复查兴隆现场输入，使用项目默认 EasyTier SSH，把脚本通过标准输入交给服务器 Python。重查会得到查询时刻的存档状态，保留原版证据目录即可同时对照。MPC 和 JPL 脚本只读查询。

```bash
news=outputs/news_release_20261124
ssh_sharp() {
  ssh -o BatchMode=yes -o StrictHostKeyChecking=yes \
    -o 'HostKeyAlias=[www.xinglong-naoc.cn]:20093' \
    -p 22 smtpipeline@100.100.1.93 "$@"
}
science_python=/home/smtpipeline/Softwares/miniconda3/bin/python
mkdir -p "$news/recheck"
for kind in exposures known candidates; do
  ssh_sharp "set -o pipefail; $science_python - $kind --end 20261123 | gzip -1" \
    < "$news/scripts/collect_server.py" > "$news/recheck/server_${kind}.jsonl.gz"
done
ssh_sharp "$science_python -" < "$news/scripts/collect_images.py" \
  > "$news/recheck/image_pixels.npz"
ssh_sharp "set -o pipefail; $science_python - | gzip -1" < "$news/scripts/collect_orbits.py" \
  > "$news/recheck/astorb_elements.csv.gz" 2> "$news/recheck/astorb_metadata.json"
```

`query_candidates.py` 以相同内部编号向 MPC WAMO 查询，`query_jpl.py` 用 JPL SBDB 固定分组条件查询。两者将更新工作区的相应 JSON，复现已交付数字请保留原快照。

`quality` 表示曝光时间、星表和几何足迹通过相应时期的质量门槛。`counted_raw` 标明档案去重选中的代表记录。`invalid_directory` 保留用于审计，年度统计另行排除。`input_known` 和 `input_unknown` 标记由现存处理产品可确认的科学输入。缺少正式未知候选目录的夜次可能触发历史候选数量门控。原始夜标签、经头时刻与目录一致性核实的夜数、合格夜数分别列出。年度主样本对旧版时间可疑的匹配做保守处理。

精简交付压缩包包含全部图件、文案、统计和脚本。复算的大体积输入快照随工作区目录留存。
