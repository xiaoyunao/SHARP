import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

def str_to_float(arr, unit=None):
    """
    将 SBDB 保存的字符串数组转换为 float。
    如果字符串中有单位（如 'AU', 'deg'），去掉单位再转 float。
    无效数值（空字符串、'--'）自动转 np.nan
    """
    out = []
    for x in arr:
        if x is None:
            out.append(np.nan)
            continue
        s = str(x).strip()
        if s in ('', '--', 'nan', 'NaN'):
            out.append(np.nan)
            continue
        if unit and s.endswith(unit):
            s = s[:-len(unit)].strip()
        try:
            out.append(float(s))
        except ValueError:
            out.append(np.nan)
    return np.array(out, dtype=float)

# -----------------------------
# 1. 读取 FITS
# -----------------------------
t = Table.read("/Users/yunaoxiao/Desktop/sitian.fits")

# -----------------------------
# 2. 提取并转换数据
# -----------------------------
a = str_to_float(t['orbit_elements_a'], unit='AU')
e = str_to_float(t['orbit_elements_e'])        # 离心率本身无单位
i = str_to_float(t['orbit_elements_i'], unit='deg')
classes = np.array(t['object_orbit_class_name'], dtype=str)

# -----------------------------
# 3. 分类准备
# -----------------------------
unique_classes = np.unique(classes)
markers = ['o', 's', 'D', '^', 'v', '*', 'p', 'h', '+', 'x']
colors = plt.cm.tab10.colors  # 默认10色循环

class_marker_color = {}
for idx, cls in enumerate(unique_classes):
    class_marker_color[cls] = (markers[idx % len(markers)],
                               colors[idx % len(colors)])

# -----------------------------
# 4. 创建画布和子图
# -----------------------------
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 30

fig, axes = plt.subplots(1, 2, figsize=(18, 9))
ax1, ax2 = axes

# -----------------------------
# 5. 绘制 a vs e（左图）
# -----------------------------
for cls in unique_classes:
    mask = classes == cls
    if cls == 'TransNeptunian Object':
        alpha_val = 0.8
        color_val = 'r'
        size_val = 100
    elif cls == 'Main-belt Asteroid':
        alpha_val = 0.1
        color_val = class_marker_color[cls][1]
        size_val = 10
    else:
        alpha_val = 0.8
        color_val = class_marker_color[cls][1]
        size_val = 40

    ax1.scatter(a[mask], e[mask],
                marker=class_marker_color[cls][0],
                color=color_val,
                s=size_val,
                alpha=alpha_val)

# q = 1.3 AU 虚线
a_line = np.logspace(np.log10(1.3), 2, 100)
e_line = 1 - 1.3 / a_line
ax1.plot(a_line, e_line, 'k--', label='q=1.3 AU')
ax1.legend()
ax1.set_xscale('log')
ax1.set_xlabel('Semimajor axis [AU]')
ax1.set_ylabel('Eccentricity')

# -----------------------------
# 6. 绘制 a vs i（右图）
# -----------------------------
for cls in unique_classes:
    mask = classes == cls
    if cls == 'TransNeptunian Object':
        alpha_val = 0.8
        color_val = 'r'
        size_val = 100
    elif cls == 'Main-belt Asteroid':
        alpha_val = 0.1
        color_val = class_marker_color[cls][1]
        size_val = 20
    else:
        alpha_val = 0.8
        color_val = class_marker_color[cls][1]
        size_val = 40

    ax2.scatter(a[mask], i[mask],
                marker=class_marker_color[cls][0],
                color=color_val,
                s=size_val,
                alpha=alpha_val,
                label=cls)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Semimajor axis [AU]')
ax2.set_ylabel('Inclination [°]')
ax2.legend(fontsize=18, loc='best')

plt.tight_layout()
plt.savefig('/Users/yunaoxiao/Desktop/asteroid_orbits.png', dpi=300)
plt.show()
