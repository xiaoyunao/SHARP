# PLAN

## Current objective

保持当前单夜正式链路稳定：

- `mask_gaia -> make_tracklet -> merge_tracklets_night -> run_linear_links_from_tracklets.py -> orbit_confirm_links.py -> summarize_single_night.py -> plot_unknown_links.py`
- 单夜默认 linking 参数固定为 `speed=5 arcsec/h`、`direction=10 deg`、`require_shared_endpoint=True`
- 单夜默认 orbit confirm 参数固定为 `max_v_kms=35`
- 单夜正式输出固定写入 `/pipeline/xiaoyunao/data/heliolincrr/<night>/rr_links`

## Milestones

1. 单夜正式流程已切换到 `run_linear_links_from_tracklets.py`
2. 15 夜流程已拆分为 `merge_tracklets_15.py -> run_rr_from_tracklets_15.py -> orbit_confirm_links_15.py`
3. `20260220` 单夜历史试验目录与分析结果已清理，只保留正式 `rr_links`
4. `20260220` 已从空目录按新单夜正式链路完整重跑到 summary、unknown catalog 和 GIF

## Outstanding issues

- 单夜正式结果当前为：`n_links=582`、`n_member_rows=1166`、`fit_ok=569/582`、`is_good=569/582`
- 当前已知检测覆盖为：`in_rr_link=1590/2304`、`rr_given_tracklet=69.01%`
- 当前失败 link 只有 `13` 条，主失败原因仍是 `max_v`
- 将 `max_v=30` 提到 `35` 后，只新增了 `4` 条边缘 unknown links，已知结果不变
- 当前已统一写出 `20260220_single_night_summary.{json,txt}`、`/processed1/20260220/L4/20260220_unknown_links.{fits,json}`、`heliolincrr/plots/20260220/*.gif`
- `summarize_single_night.py` 还需要把用户明确关心的 6 组统计口径固化为固定字段，避免后续每次手算
- 当前单夜 purity 已足够高，后续如果继续优化，核心问题只剩 completeness
- 单夜 formal path 已清理完成，后续不要再在正式目录旁边长期保留 `linear_links_*`、`rr_links_*` 试验目录
- 15 夜 `w15` 流程仍需继续与单夜完全分开维护

## Validation criteria

- `20260220` 能稳定从空目录重跑 `mask_gaia` 和纯 `make tracklet`
- `run_linear_links_from_tracklets.py` 能稳定写出 `links_tracklets.fits`, `linkage_members.fits`, `link_edges.fits`, `rr_summary.json`
- `orbit_confirm_links.py` 能稳定在正式 `rr_links` 下写出 `orbit_confirm/orbit_links.fits` 与 `orbit_obs_residuals.fits`
- `summarize_single_night.py` 能稳定写出 summary、unknown catalog、以及固定的 6 组 requested metrics
- `plot_unknown_links.py` 能稳定按 unknown catalog 产出 GIF 和 summary
- 正式单夜目录中不保留旧试验目录与旧分析副本

## Next recommended steps

1. 先补全 `summarize_single_night.py` 的 requested metrics 字段并验证到现有夜次
2. 保持当前单夜正式参数 `speed=5`, `direction=10`, `require_shared_endpoint=True`, `orbit max_v=35` 不变
3. 若后续继续提升单夜召回，只围绕线性 linking 的边规则做小范围扫描
4. 15 夜流程继续只通过 `run_pipeline_15.sh` 及其 `_15` 脚本单独维护和记录
