# SHARP TJN固定天区交接包

日期：2026-11-11。先读 `CODEX_TJN_FIXED_FIELDS_20261111.md`，长期任务见 `SHARP_PROJECT_MEMORY_20261111.md`。

本包交付明确的网格方案、旧MP/年度B输入、1063候选ID、分云结果、几何草案与复现脚本。只涉及TJN天区，不修改其他观测策略。

## 主要文件

| 文件 | 内容 |
|---|---|
| `tjn_fixed_grid_seed_v0.csv` | 522场草案：L4 241，L5 281；未冻结 |
| `TJL4_seed_v0.csv`、`TJL5_seed_v0.csv` | 相同草案的两云分表 |
| `tjn_grid_proposal.json` | 母网格、编号与模板契约 |
| `01_tjn_fixed_seed_atlas.png/pdf` | 两云旧中心与固定场草案 |
| `02_tj_test_grid_comparison.png/pdf` | 今晚TJ_test与恒Dec网格对照 |
| `seed_summary.json` | 样本、草案范围、名义覆盖核算与限制 |
| `legacy_to_seed_cells.csv` | 旧场到草案单元的多对多关系 |
| `inputs/` | 原始输入、候选/分云中间数据 |
| `audit/` | 独立名义TAN检查、来源与交付核验 |
| `scripts/` | 三个独立复现脚本，不接望远镜 |

522场来源于旧B足迹包络，完整1,063对象全年轨迹和实际WCS须由Codex从已有本地缓存补充，之后才能冻结正式场库。草案中的旧对象成员只是来源记录，不是新场成员预测；本包没有可直接上望远镜的执行列表。

所有人类可读编号在本包中都是草案，正式冻结后才实行不可重排/不可移动规则。数据坐标单位为度，中心为ICRS；名义足迹为保守1.15″/pix的TAN模型，实际几何尚待核实。

依赖：Python、NumPy、pandas、matplotlib。复现命令见主交接文档第11节。
