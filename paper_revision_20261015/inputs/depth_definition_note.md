# 历史 MLIM、参考零点和孔径星等

## 校准链路

历史 `calib_flux.py` 按 CCD9k 的固定孔径直径列表找到 18 px 项。在参考星的 Gaia 派生 W 星等与仪器 18 px 孔径星等之差上做 2.5σ 裁剪，取中位数作为逐帧 `CCDZP`。源代码的分支是优先参考星等小于 15 的匹配星，数量少于 50 时放宽到 18。各固定孔径同步保存 `APERdd`，其符号定义是参考 18 px 仪器星等减 dd px 仪器星等的中位数。

设 $f_d$ 是直径 $d$ 的背景扣除孔径通量，以 ADU 计，$Z$ 是 `CCDZP`，$C_d$ 是 `APERdd`。历史 `photometry/phot.py` 先按零点换算到目录通量单位，再应用孔径修正。以 `Mag_Aper4` 为例，$d=8$，

```text
F_d = f_d × 10^[0.4 × (22.5 − Z − C_d)]
m_d = 22.5 − 2.5 log10(F_d) = Z − 2.5 log10(f_d) + C_d
```

同一因子作用于通量和其误差，因此不会改变相应固定孔径的通量信噪比。目录单位使用 22.5 的 nanomaggy 数值约定，光度本身采用代码定义的 Gaia 参考 W 标尺。

## 5σ 模型量

```text
A_noise = π × (FWHM_pix / 2)^2
F5_ADU = 5 × sqrt[A_noise × (sky_ADU_per_pix × gain_e_per_ADU + RON_e²)] / gain_e_per_ADU
MLIM = CCDZP − 2.5 log10(F5_ADU)
```

`getmlim` 接受的是 18 px 参考零点、背景 ADU、FWHM 像素、增益 e⁻/ADU 和读出噪声电子数。FWHM/2 只用于设置噪声面积。代码既没有使用 `APER08`，也没有对 FWHM/2 半径内的包围能量另做因子或星等修正。

18 份历史 L2 的内容 SHA-256 与上一轮一致。用头值复算 `MLIM`，最大差 0.00483 mag。图像 FWHM 在 2.200–3.085 px，参考直径均为 18 px，`APER18=0`，`APER08` 在 −0.2121 至 −0.0990 mag。对每帧 582–1,192 个无标记、高星状可靠度的亮源，把目录通量逆去光度孔径修正，恢复固定孔径的相对增长曲线。包夹各帧 FWHM 直径的 2 px 与 4 px 孔径，跨帧中位增长比 $f_d/f_{18}$ 分别约为 0.197 和 0.542。FWHM 直径处的具体比值未被直接采样，保留这两端而不采用假设 PSF。18 帧的 `Mag_Aper4`、归档通量、逐帧零点和 `APER08` 符合上述代数关系。

可用于正文的英文表述：

> The median catalog background-noise model depth is 20.25 mag on the Gaia-referenced white-light scale. The calculation uses the 18-pixel reference-aperture zero point and a circular noise area of radius FWHM/2, with no separate encircled-energy factor.

[18 帧头值复算](../tables/depth_aperture_reconciliation_18.csv) · [固定孔径增长曲线](../tables/depth_aperture_growth_curves_18.csv) · [上游文件指纹](calibration_provenance.md)
