# PLAN

## Current objective

为项目汇报准备一批可直接上 PPT 的图和动图，并把素材统一放到服务器
`/pipeline/xiaoyunao/ppt_assets_20260414`：

- `survey` 全年 `365` 晚模拟覆盖动画
- `known_asteroid` 的 `asteroid_orbits` 与 `nside64_counts`，且要与 `/Volumes/Foundation/Asteroid/sitian_stats.ipynb` 保持一致
- `heliolincrr` 的 Gaia masking + 静止源剔除流程图

## Milestones

1. 已在服务器建立独立素材目录 `/pipeline/xiaoyunao/ppt_assets_20260414`
2. `known_asteroid` 最新总历史表已补齐轨道元素，并已生成 notebook 直读版 `/Users/yunaoxiao/Desktop/sitian.fits`
3. 已按 `sitian_stats.ipynb` 原代码在本机写出 `asteroid_orbits.png` 与 `sitian_ra_dec_healpix_nside64_counts_linear_clim0_200.png`
4. `heliolincrr` 已产出 Gaia masking + 静止源剔除流程图
5. `survey` 12 晚测试动画已跑通并成功写出 GIF
6. `survey` 365 晚正式动画已在后台运行

## Outstanding issues

- `survey_year_full` 仍在后台运行，尚未完成全部 `365` 帧
- 服务器无 `Times New Roman`，因此服务器端 `gaia_static_flow` 只能近似 notebook 风格，不能保证字形完全一致
- 若后续还要微调 `known_asteroid` 图，原则上应直接改桌面 notebook 或其输入 `sitian.fits`，不要再随意改绘图参数
- `Gaia masking` 当前已升级为包含静止源剔除的流程图；若希望更“动图化”，还可再做局部放大或 before/after 闪烁版
- 这批汇报图暂时放在独立目录，还未正式并入 README 或长期统计流程

## Validation criteria

- 桌面 `/Users/yunaoxiao/Desktop/sitian.fits` 存在且包含：
  - `orbit_elements_a`
  - `orbit_elements_e`
  - `orbit_elements_i`
  - `object_orbit_class_name`
- 桌面存在按 notebook 原代码写出的：
  - `/Users/yunaoxiao/Desktop/asteroid_orbits.png`
  - `/Users/yunaoxiao/Desktop/sitian_ra_dec_healpix_nside64_counts_linear_clim0_200.png`
- `/pipeline/xiaoyunao/ppt_assets_20260414/outputs/gaia_masking_test_fast/` 能稳定写出：
  - `20260220_gaia_masking_story.png`
  - `20260220_gaia_masking_story.gif`
  - `20260220_gaia_masking_story.json`
- `/pipeline/xiaoyunao/ppt_assets_20260414/outputs/gaia_static_flow/` 能稳定写出：
  - `20260220_gaia_static_flow.png`
  - `20260220_gaia_static_flow.json`
- `/pipeline/xiaoyunao/ppt_assets_20260414/outputs/survey_test/` 能稳定写出测试 GIF
- `/pipeline/xiaoyunao/ppt_assets_20260414/logs/survey_year_full.log` 持续推进直到写出正式 `365` 晚 GIF

## Next recommended steps

1. 盯 `survey_year_full.log`，待 `365` 晚动画写完后抽查头尾帧和 GIF 节奏
2. 若用户还想继续调 `known_asteroid`，优先直接在桌面 notebook 上改，因为当前图已经是按原代码生成
3. 若需要更强视觉冲击，再补 `Gaia masking` 的 zoom-in 或闪烁对比动图
4. 用户确认最终想要的导出格式后，再决定是否额外输出 mp4 或更大尺寸 PNG
