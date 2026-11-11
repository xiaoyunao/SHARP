# SMT r-band focus measurement

独立工具；用户于 2026-11-11 授权新增 `focus/` 并在服务器临时运行。生产目录、原始 FITS、定时任务均不修改。程序和结果保留本地；服务器临时目录在结果传回并逐文件 SHA256 验证后删除。SSH/系统审计日志不属于本工具清理范围。

## 方法

- 读取全幅 16 放大器，不抽样选星；SEP 3σ、至少12像素检测，逐源局部背景扣除。
- 实测0083原图亮星约53800 ADU平顶（全图最大54422），头SATURATE=65535不可信；去overscan后平台约51000，保守屏蔽>45000 ADU像素及包含它们的源，避免饱和尺寸随亮度增大污染排序。`--saturation`可覆盖，不能未经检查直接信头值。
- 使用所有符合可靠测量要求的源：S/N≥30，剔除饱和/坏像素、边界、混叠标志、明显宇宙线尺度、邻源污染和增长曲线不收敛。逐放大器记录每种剔除原因。
- 主量为 D50（包围50%光的直径）及 D80（80%），不把 D50 冒称 FWHM。光度孔径≥12像素或6倍长轴二阶矩；扩大20%后两个半径变化均须≤10%。环状源仍能测量；环状源不报告单峰径向 FWHM。
- 比较时平移配准、双向唯一匹配同一批星；各4×4等面积区域内求对数尺寸比中位数，再跨区域取中位数，避免星密度不均主导排序。至少8区、每区10星，否则拒绝给出匹配排序。
- D50 排序，D80、亮星 S/N>80 子样本、椭率和空间分布交叉检查。直接比较领先图与每幅图，不仅比较各自对第一幅的误差条。
- bootstrap 区间只表示源采样误差，不能消除不同曝光的 seeing 变化。尺寸差很小或 D50/D80 不一致时须结合复拍判断。曝光内星像质量可比较，不能从单波段焦扫独立认定 band 间焦差的原因。
- 源坐标为去overscan后图像的0基像素坐标；尺寸CSV单位pixel，摘要同时给出头中 `PIXSCAL1` 换算的角秒。

源为未经星表分类的恒星候选；中位数及同源比值抵抗少量星系/未分辨双星污染。低S/N源无法可靠测尺寸，保留拒绝计数。未给出精确相机gain时，S/N使用背景噪声、局部背景估计误差和保守1 e/ADU泊松项，属于筛选估计。

## 服务器环境及原图处理

解释器：`/home/smtpipeline/Softwares/miniconda3/bin/python`，已确认 NumPy 2.3.3、SciPy 1.16.2、Astropy 7.1.0、SEP 1.4.1、Matplotlib 3.10.6；不安装依赖。

`.fits.gz` 的 subscan 图为 9232×9216。新 `focus-...fits` 不含 `OBJ`，生产预处理会跳过，故用独立的16通道行overscan算法直接读 9400×10000 原图。布局/异常overscan处理根据只读核查 `/pipeline/smtpipeline/ccd9k.py` 实现，不导入生产程序。

首幅 raw 与 subscan 两条路径独立运行：修正饱和阈值后1,754个源的所有13列逐项完全一致，证据在结果目录 `raw_processed_validation_corrected.json`（初版1,832源记录仍留档）。

## 本地发起远程运行

在本项目目录执行（只给明确的 r 波段调焦曝光；程序也检查FILTER）：

```bash
python3 focus/run_remote.py \
  /processed1/20261111/subscan/OBJ/OBJ_MP_1831_0080.fits.gz \
  /raw1/20261111/OBJ/focus-5.84_0081.fits \
  /raw1/20261111/OBJ/focus-5.86_0082.fits \
  /raw1/20261111/OBJ/focus-5.88_0083.fits \
  /raw1/20261111/OBJ/focus-5.90_0084.fits \
  /raw1/20261111/OBJ/focus-5.92_0085.fits \
  --out focus/results/20261111
```

后续仅传新帧即可：已有本地源表/摘要被复制到本次临时目录，测量后在服务器重新比较整个本地批次；旧结果不会因只新增一帧而消失。FITS读前后大小/mtime须一致；尚在写入则拒绝本次测量。单帧失败时不宣称已有完整排序。

`last_remote_run.json` 保存本次远程临时路径、SHA256和清理状态。中断/断网时检查 `cleaned`；若为false，只处理该记录中的确切临时目录，不清理其他任务。

结果包括：逐源CSV、逐帧JSON、全场尺寸图、固定31×31像素源切图、`comparison.json`、`comparison.png`、`RESULTS.md`。未指定范围时结果保持 provisional，最终比较范围必须由现场确认，不能靠一段时间无新文件推断结束。用户已确认本次用到0088，最终选择0083。

仅按已测数据定稿（不重测FITS）：

```bash
python3 focus/run_remote.py --out focus/results/20261111 --final-through 88
```

定稿要求源表序号连续、末帧等于指定值、所有帧都能排序；在服务器临时目录重算比较、校验传回并清理。

## 验证

`python3 -m py_compile focus/*.py`；在上述服务器环境中运行 `python -B validate_synthetic.py`。解析高斯D50误差分别0.88%/0.18%，Moffat为2.65%；环状源识别、宽度排序、饱和/宇宙线剔除通过。相近帧排序还以真实同源、亮星和空间分区复核。

SEP包围能量实现参考：[flux_radius API](https://sep.readthedocs.io/en/stable/api/sep.flux_radius.html)。
