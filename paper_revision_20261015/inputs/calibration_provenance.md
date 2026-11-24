# 上游定标与星历来源

SHARP CCD9k 的原始读出为 10,000 × 9,400 像素。上游 overscan 裁剪后为 9,216 × 9,232 像素。历史 L2 头的 `FILTER=W`、`BAND=w`、`CALI_REF=Gaia` 与 CCD9k 配置一致。SExtractor 的孔径直径列表以 2、4、6、8 像素开始，因此 `Mag_Aper4` 的测量直径是 8 像素。参考定标使用直径 18 像素的孔径。W 颜色变换和深度算法见[论文事实表](../MANUSCRIPT_FACTS.md)，逐曝光核对见[光度头复算](../tables/sample_photometric_header_and_depth_check.csv)。

Aleph 1.0.0 中的 `query_mixed_cat` 先从小行星延迟历元的日心赤道位置扣除将兴隆站址由地固坐标换算得到的日心位置，再输出 RA、Dec。旧运行读取冻结 DE432s。论文分析在同一站点原点重算日距角，并用 DE432s 提供的独立行星矢量交叉核对。

| 核对文件 | SHA-256 |
|---|---|
| `calib_flux.py` | `d1c728ccbc962aa0e16cbf25998cad168098d2893c831af12131c11ea6afb014` |
| `ccd9k.py` | `23ef8a2279bc5cea67a7949ab1f0b4d394b698eba92a7c1d50261ff5eac9bb0e` |
| `photometry/phot.py` | `955d3060de929a438e7705c779b80d1d916334379392ff1ef7ac1d88bf0b2824` |
| `photometry/config1.sex` | `4a63ac6544582fe0f0b299d48c27e0a5a691e097d548594286dc6a85b433c095` |
| `aleph/Query.py` | `137ae6e8572f9f33c07d777b17e21d372a6e9e091db764ffc1a8dda47fb331d7` |
| `aleph/aleph_utils.py` | `80523e7464ee6478f1676728d05c580df3eaf896656dc185ba21da51f9d2eb6e` |
