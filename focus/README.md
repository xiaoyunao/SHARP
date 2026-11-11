# SMT 调焦主程序

独立读取服务器图像并测量星像。生产目录、原始FITS和定时任务保持原状；结果传回并校验SHA256后删除本次服务器临时目录。今晚结果索引见 [results/README.md](results/README.md)。

## 保留的主程序

- `run_remote.py`：本地发起SSH、复制必要模块、运行、传回校验和清理。
- `measure_focus.py`：全幅16通道overscan、源筛选、D50/D80与形态测量。
- `compare_focus.py`：同波段焦扫的同星/分区比较，供主入口调用。
- `compare_passbands.py`：不同波段完全共同源的图像质量对照，不推导滤镜焦差。

测试和一次性复核脚本已清理，需要恢复时查看git历史。结果中保留逐源CSV、逐帧摘要、最终图表与必要波段来源信息；校验结论汇总于`results/provenance.json`。

## 运行

在本项目目录本地执行，实际计算在默认服务器：

```bash
python3 focus/run_remote.py /raw1/20261111/OBJ/focus-5.88-G_0092.fits \
  --band g --allow-header-mismatch --out focus/results/new_g_run
```

`--band`支持g/r/i/z，实际滤镜已确认但图头陈旧时才使用`--allow-header-mismatch`；未确认用`--band unconfirmed`。摘要分别记录实际band与header_band，不改原图头。不同波段应使用独立输出目录。已有目录会带入其源表/摘要继续比较；后续新焦扫使用新目录。

仅对已有同波段结果定稿：`python3 focus/run_remote.py --out <结果目录> --final-through <最后序号>`。不指定截止时报告保留provisional。`last_remote_run.json`在运行时记录传输SHA和清理状态；断网/中断后以该记录核查本次临时目录。

跨波段对照：`python3 focus/compare_passbands.py <参考摘要.json> <其他摘要.json> ... --out <输出目录>`，第一个输入为参考，必须每band一帧；可将程序及源表传到服务器临时目录执行。

## 测量口径与环境

- SEP全幅检测，使用所有通过S/N≥30、饱和/坏像素/边界/混叠/增长曲线收敛筛选的源。D50/D80为包围50%/80%光的直径，D50不是FWHM。
- 同星双向匹配，4×4分区等权；D80、亮星子样本、椭率和空间分布交叉检查。bootstrap只含源采样误差，不含seeing系统变化。
- 2026-11-11实际饱和平台低于65535头值；去overscan后采用45000 ADU保守阈值。源图采用完整峰值线性显示，避免99%截白造成假胖核。
- 服务器解释器`/home/smtpipeline/Softwares/miniconda3/bin/python`；NumPy、SciPy、Astropy、SEP、Matplotlib。不导入生产代码，不安装新依赖。
- 支持9400×10000 raw及9232×9216 subscan，原图在内存独立去overscan。CSV坐标为去overscan后0基像素，宽度为pixel；摘要按PIXSCAL1另给角秒。
- 单一焦位不同波段的图像大小不等于各自最优焦点；曝光时长、seeing和成像位置可能影响结果。
