"""Paper-style fixed-field atlas and a self-contained technical handoff."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from astropy.time import Time
from .fixed_fields import Grid,boundary
from .geometry import offset,coord
from .fixed_field_build import OUTPUT,HANDOFF,write_json
from .plot_style import apply_style,save_pair,BLUE,ORANGE,GREEN


def draw(output=OUTPUT):
    apply_style();figdir=output/'figures';figdir.mkdir(exist_ok=True)
    doc=json.loads((output/'registry.json').read_text());g=Grid(**doc['grid']);fields=pd.DataFrame(doc['fields'])
    archive=np.load(output/'annual_samples.npz');samples={k:archive[k] for k in ('ids','strict','ra','dec')};candidates=pd.read_csv(HANDOFF/'inputs/tjn_annual_candidates_1063.csv',dtype={'object_id':str}).set_index('object_id')
    old=pd.read_csv(HANDOFF/'inputs/tjn_B_unique_fields_with_cloud.csv');colors={'L4':BLUE,'L5':ORANGE}
    def outlines(ax,rows,color,margin=0,lw=.45,alpha=.5,ls='solid'):
        segments=[boundary(r.ra_deg,r.dec_deg,g,margin=margin,samples=17) for r in rows.itertuples()]
        ax.add_collection(LineCollection(segments,colors=color,linewidths=lw,alpha=alpha,linestyles=ls))
    def tracks(ax,cloud):
        for i,k in enumerate(samples['ids']):
            if candidates.loc[k,'cloud_proxy']!=cloud:continue
            select=np.flatnonzero(samples['strict'][i])[::12]
            ax.scatter(samples['ra'][i,select],samples['dec'][i,select],s=.7,color=colors[cloud],alpha=.35,rasterized=True)
    fig,ax=plt.subplots(figsize=(14,6),layout='constrained')
    for cloud in colors:
        fs=fields[fields.planning_cloud==cloud];outlines(ax,fs,colors[cloud],alpha=.4);tracks(ax,cloud)
        ax.scatter([],[],s=30,color=colors[cloud],label=f'{cloud}: {len(fs)} fixed fields')
    ax.scatter(old.ra_deg,old.dec_deg,marker='x',s=8,c='.2',lw=.6,label='383 legacy B centres')
    ax.set(xlim=(270,20),ylim=(-25,61),xlabel='ICRS RA [deg]',ylabel='ICRS Dec [deg]',title='SHARP fixed TJN fields | 2026 Nov — 2027 Oct')
    ax.legend(loc='upper center',ncol=3,fontsize=11);ax.grid(alpha=.15)
    fig.text(.5,-.02,'2.50° Dec rows; RA spacing = 360° / ceil(360° cos Dec / 2.50°). Instrument acceptance pending.',ha='center',fontsize=11)
    save_pair(fig,figdir/'01_annual_atlas');plt.close(fig)
    cases={'L4':'L4_spring_2027-04-15.json','L5':'L5_winter_2026-12-15.json'}
    for n,cloud in enumerate(colors,2):
        fs=fields[fields.planning_cloud==cloud];b=old[old.cloud_proxy==cloud]
        plan=json.loads((output/'representative_nights'/cases[cloud]).read_text());active={f['field_id'] for block in plan['blocks'] for f in block['fields']}
        fig,ax=plt.subplots(figsize=(12,8),layout='constrained');outlines(ax,fs,'.55',alpha=.45);tracks(ax,cloud)
        outlines(ax,fs[fs.field_id.isin(active)],colors[cloud],lw=2,alpha=1)
        ax.scatter(b.ra_deg,b.dec_deg,s=8,marker='x',c='.25',lw=.5)
        ax.set(xlim=(fs.ra_deg.max()+4,fs.ra_deg.min()-4),ylim=(fs.dec_deg.min()-4,fs.dec_deg.max()+4),xlabel='ICRS RA [deg]',ylabel='ICRS Dec [deg]',title=f'{cloud}: {len(fs)} fixed fields | active example {plan["date"]}')
        ax.legend(handles=[Line2D([],[],color='.55',label='Static nominal footprints'),Line2D([],[],color=colors[cloud],lw=2,label='Four scheduled fields'),Line2D([],[],color=colors[cloud],marker='.',ls='',label='Annual visible trajectories')],loc='upper center',ncol=3,fontsize=10)
        ax.grid(alpha=.15);save_pair(fig,figdir/f'0{n}_{cloud}_detail');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(14,6),layout='constrained');ax=axes[0]
    top=fields.loc[fields.dec_deg.idxmax()];fs=fields[(abs(fields.ra_deg-top.ra_deg)<6)&(fields.dec_deg>=top.dec_deg-5)]
    outlines(ax,fs,ORANGE,lw=1,alpha=.9);outlines(ax,fs,BLUE,margin=g.margin+g.guard,lw=.8,alpha=.9,ls='dashed')
    ax.scatter(fs.ra_deg,fs.dec_deg,s=14,c='k')
    ax.set(xlim=(top.ra_deg+6,top.ra_deg-6),ylim=(top.dec_deg-6.6,top.dec_deg+2),xlabel='ICRS RA [deg]',ylabel='ICRS Dec [deg]',title='High-Dec rows: nominal / quality + guard')
    ax.grid(alpha=.15)
    row=pd.read_csv(output/'row_geometry_audit.csv');ax=axes[1]
    for key,color,label in [('EQ25-v1',BLUE,'2.50° adopted'),('MP-pitch-proposal',ORANGE,'2.67° proposed')]:
        d=row[row.grid==key];ax.plot(d.dec_deg,d.minimum_edge_deg*60,'o-',ms=3,color=color,label=label)
    ax.axhline(9.6,color='k',ls='--',label='0.15° quality + 0.01° guard')
    ax.set(xlabel='Row centre Dec [deg]',ylabel='Minimum cell-to-edge margin [arcmin]',title='Dense TAN ownership-cell boundary audit')
    ax.legend(fontsize=10);ax.grid(alpha=.15);save_pair(fig,figdir/'04_edge_geometry');plt.close(fig)
    test=pd.read_csv(output/'TJ_test_tracks.csv.gz',dtype={'object_id':str});prov=json.loads((output/'TJ_test_provenance.json').read_text())
    fig,ax=plt.subplots(figsize=(12,8),layout='constrained')
    fs=fields[(fields.ra_deg>71)&(fields.ra_deg<85)&(fields.dec_deg>19)&(fields.dec_deg<29)]
    outlines(ax,fs,'.6',alpha=.55);outlines(ax,fs,'.6',margin=g.margin+g.guard,alpha=.6,ls='dotted')
    q=np.linspace(-1.5,1.5,65);xx=np.r_[q,np.full(65,1.5),-q,np.full(65,-1.5)];yy=np.r_[np.full(65,-1.5),q,np.full(65,1.5),-q]
    bp=offset(coord(79.002262,23.719287),xx,yy);oldb=np.c_[bp.ra.deg,bp.dec.deg];ax.plot(oldb[:,0],oldb[:,1],color=GREEN,ls='--',lw=2,label='Unchanged TJ_test nominal footprint')
    for k,color in [('17365',BLUE),('4722',ORANGE)]:
        t=test[(test.object_id==k)&(test.period=='35days')];n=test[(test.object_id==k)&(test.period=='tonight')]
        ax.plot(t.ra_deg,t.dec_deg,color=color,lw=2,label=k+' | next 35 days')
        ax.scatter(n.ra_deg,n.dec_deg,s=5,color=color,zorder=5)
    involved=set(';'.join(test.quality_field_ids.fillna('')).split(';'))
    for f in fs[fs.field_id.isin(involved)].itertuples():ax.text(f.ra_deg,f.dec_deg+.2,f.field_id,fontsize=9,ha='center')
    ax.scatter(79.002262,23.719287,marker='x',color=GREEN,s=80)
    ax.set(xlim=(83,71),ylim=(20,28),xlabel='ICRS RA [deg]',ylabel='ICRS Dec [deg]',title='TJ_test transition | the fixed grid splits tonight’s pair')
    ax.legend(loc='upper left',fontsize=10);ax.grid(alpha=.15);save_pair(fig,figdir/'05_TJ_test_transition');plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(12,6),sharex=True,layout='constrained')
    for ax,k,color in zip(axes,['17365','4722'],[BLUE,ORANGE]):
        t=test[(test.object_id==k)&(test.period=='35days')].copy();day=(Time(t.time_utc.tolist()).unix-Time('2026-11-11T04:00:00').unix)/86400
        names=sorted(set(';'.join(t.quality_field_ids.fillna('')).split(';'))- {''})
        for j,fid in enumerate(names):
            inside=t.quality_field_ids.fillna('').map(lambda s:fid in s.split(';')).to_numpy()
            ax.scatter(day[inside],np.full(inside.sum(),j),color=color,s=3)
        ax.set(yticks=range(len(names)),yticklabels=names,title=f'{k}: geometric quality-footprint membership (daylight included)')
        ax.grid(alpha=.15)
    axes[-1].set(xlabel='Days since 2026-11-11 12:00 Beijing time',xlim=(0,35));save_pair(fig,figdir/'06_target_field_crossings');plt.close(fig)


def report(output=OUTPUT):
    s=json.loads((output/'summary.json').read_text());a=json.loads((output/'independent_audit.json').read_text());w=json.loads((output/'wcs_audit.json').read_text())
    fields=pd.read_csv(output/'registry.csv');objects=pd.read_csv(output/'objects.csv');special=pd.read_csv(output/'special_targets.csv');nights=pd.read_csv(output/'representative_nights.csv')
    lines=['# TJN 固定监测天区：2026-11-11 候选发行版', '',
        '**已完成全年轨迹支持的名义固定场划分；状态为 RC1，尚未正式冻结为仪器认证场库。** 实测 TPV 模型检查通过，但真实坏像素/渐晕遮罩、模板公共裁边和高低赤纬重复指向/PA 分布尚缺。今晚控制脚本及既有生产代码保持原状。','',
        '## 场库和几何选择','',
        '| 分区 | 场数 | 编号 | 中心 RA 范围 | 中心 Dec 范围 |','|---|---:|---|---|---|']
    for cloud,g in fields.groupby('planning_cloud'):
        lines.append(f'| {cloud} | {len(g)} | TJ{cloud}_001—TJ{cloud}_{len(g):03d} | {g.ra_deg.min():.2f}°—{g.ra_deg.max():.2f}° | {g.dec_deg.min():+.1f}°—{g.dec_deg.max():+.1f}° |')
    lines += ['', '这些范围仅为分离区域的外包络，不能当成满铺矩形。每晚只取其中少量场。静态注册表不写目标名单；逐时归属另外保存。', '',
        '采用 RA/Dec=0°锚点、Dec 行距 **2.50°**，同行 `N=ceil(360*cos(Dec)/2.50)`、`RA=k*360/N`。相机设计轴约定为 +x 北、+y 东，PA=0°，宇称 −1；这是一项长期几何约定，硬件实际 PA 仍需标定。', '',
        f'原建议 2.67264°/2.67496°节距在保守 TAN 模型下，保留单侧 0.15°测光留边和 0.01°轨迹余量后，{s["proposed_pitch_quality_loss_samples"]:,}/{s["visible_samples"]:,} 个可见采样点（{100*s["proposed_pitch_quality_loss_samples"]/s["visible_samples"]:.2f}%）即使允许邻场也无法进入安全区。故在首次定版前整体缩小节距；没有移动任何已正式冻结的旧场。未优化全局相位，保留建议的零锚点。', '',
        f'采用 9216×1.15″=2.944°的保守内接正方形；不把约 3°、1.16″、9232/9224 像素口径混用。逐行每边 2,049 点，最小名义留边 {s["nominal_min_cell_edge_deg"]*60:.3f}′，超过所需 9.6′。旧 `in_field(... margin=.15)` 是移动目标质量条件，并非实测模板裁边。', '',
        '## 全年覆盖、选择函数和模板代价','',
        '年度为北京时间观测夜 2026-11-01—2027-10-31。原 16,445 母表及 1,063 候选保持不变；L4/L5 来自原日心黄经代理标签，不按视赤经猜云。逐文件核验全部 6 小时站心 ICRF/UTC Horizons 缓存与既有哈希，未重新下载母表或星历。', '',
        '30 分钟采样采用 RA unwrap 插值；用太阳/高度 3.88°时间包络和星历运动余量保留可能的观测片段，月亮不永久裁掉场。严格统计口径为目标高度≥32.5°、太阳≤−18°、V≤19.2；核心代理 V≤18.5。这里的“可见天数”表示该夜至少一个采样点，不能等同完整两小时或有效测光夜。', '',
        '| 样本 | L4 | L5 | 合计 |','|---|---:|---:|---:|',
        '| 原年度候选 | 610 | 453 | 1063 |']
    for label,col in [('至少一次严格可见','visible_samples'),('至少一次核心亮度代理','core_samples')]:
        count=objects[objects[col]>0].groupby('cloud').size();lines.append(f'| {label} | {count.get("L4",0)} | {count.get("L5",0)} | {sum(count)} |')
    visible=objects[objects.visible_samples>0]
    lines += ['',f'- 严格可见 {s["visible_samples"]:,} 点全部在新库质量模型内；保守选择包络 {s["conservative_envelope_samples"]:,} 点。未可见对象保留在逐对象表，不计为“覆盖失败”。',
        f'- 旧 522 场草案只覆盖严格样本的 {a["seed_nominal_covered"]:,} 点（{100*a["seed_nominal_covered"]/a["strict_samples"]:.2f}%）；{a["seed_objects_with_any_nominal_miss"]} 个对象有部分轨迹不在草案名义足迹内。',
        f'- 独立检查 {a["swept_interior_points"]:,} 个间插点，最小最佳留边 {a["swept_minimum_best_edge_deg"]*60:.3f}′；旧 B 的 383 场另以 61×61 共 {a["legacy_footprint_points"]:,} 点检查，无质量区漏点。',
        f'- 可见对象按每天中间采样的主场统计，年内切场中位 {visible.primary_day_field_changes.median():.0f} 次、最大 {visible.primary_day_field_changes.max()} 次；连续可见/核心天数及逐日多场关系见 objects/object_days。该主场计数不是强制的实际切场策略。',
        f'- 单场活跃采样时刻同时可用目标均值 {a["overlap_mean_simultaneous_targets"]:.2f}、最多 {a["overlap_max_simultaneous_targets"]} 个，计入重叠伴随对象；全库覆盖按 object_id 去重。',
        f'- r 模板几何槽位至少 {s["r_template_count"]} 个，四滤镜 {s["griz_template_count"]} 个；实际还需乘以模板版本/几何变体，不能理解为所需曝光数或已获得模板。',
        '- `field_occupancy.csv` 的 core_active 指至少20个严格可见的主场采样；其余为候补/边界标签。标签可更新，不改变中心。原 B 的质量光变、唯一周期产额没有继承，也没有开展新天气模拟。','',
        '## 实测 WCS 的证据范围','',
        f'恢复本地 4 个历史场的 12 份完整 TPV WCS，裁剪 CRPIX 加回裁剪原点。目录 X_Win/Y_Win 从1起算，转换为0起算后重投影误差≤{w["origin_zero_reconstruction_max_pix"]:.2g}像素。每份边缘1,028点正反变换通过。像元尺度约1.155″，+x PA约−0.336°至−0.271°，+y向东，宇称−1。', '',
        f'将这些完整畸变模型平移到所有候选行的情景检查，最小边缘余量 {w["minimum_transported_edge_pix"]:.1f}像素，超过约500.9像素的0.16°要求。这些原始样本实际赤纬仅约−3°至+8°，不能冒称高赤纬硬件PA和模板遮罩已经实测通过。`in_usable`/`usable_mask_ref` 保留空值，正式冻结状态明确为 false。', '',
        '使用完整 [Astropy WCS 变换](https://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html)，而非只凭四角或RA–Dec外包矩形。','',
        '## 今晚 TJ_test 与35天留场','',
        '**两颗目标在新固定网格中会分开，不能把今晚中心直接替换成附近新场。** 17365今晚102帧落入 `TJL5_436`，4722落入 `TJL5_468`，没有共同质量场。附近新中心 `TJL5_468` 为 (79.694656488550°, +25.000000000000°)，17365今晚在该场内的比例为0。', '',
        '未来35天，两者在固定库任意质量场中的几何覆盖都是100%，但17365依次涉及434/435/436，4722涉及466/467/468。4722在468场内约26.90%的全时采样、23.76%的可见采样；17365始终不在468内。逐帧表和轨迹图单独交付。单场组合效率发生损失，需要未来队列优化，不能用旧两颗同场产额代替。', '',
        '今晚 `20261111.txt` SHA256仍为 `e43db50712f924f382b9f776f62d825abca5d8dc54d46e598427b29daf78147b`；旧图像保留 `TJ_test / AF_148ecdaae5` 几何身份。','',
        '## 特殊目标','', '| 对象 | 原年度样本 | 可见夜数 | 核心代理夜数 |','|---|---|---:|---:|']
    for r in special.itertuples():
        lines.append(f'| {r.object_id} {r.name} | {"是" if r.in_annual_1063 else "否"} | {int(r.visible_days) if pd.notna(r.visible_days) else "未有全年轨迹"} | {int(r.core_days) if pd.notna(r.core_days) else "—"} |')
    lines += ['', '38615/228155原预筛因亮度或赤纬拒绝，本轮没有静默加入1063比较样本；其余未可见对象仍保留。Hektor等受本站高度条件限制；这不代表科学价值不足。', '',
        '## 有限接入验证','', '| 日期/案例 | 场数 | r120帧 | 独立目标 | 转场秒 | 验证 |','|---|---:|---:|---:|---:|---|']
    for r in nights.itertuples():lines.append(f'| {r.date} / {r.case} | {r.fields} | {r.exposures} | {r.unique_targets} | {r.slew_s:.1f} | 通过 |')
    lines += ['', '仅新增 opt-in `fixed_trojan` 查询入口，原调度入口保持。采用既有2小时预算、r120、驻场两张、该次配置的读出时间、5分钟交接和原转场模型，逐曝光5秒几何审查；RA机械轴位移仍使用wrap后的RA差，未乘cos(Dec)。每夜4场40帧是这些代表夜的结果，不是全年承诺。','',
        '## 文件和复现','',
        '- `registry.csv` / `registry.fits` / `registry.json` 与 `TJL4.csv`、`TJL5.csv`：静态场库；编号在RC后续增补时也只追加，正式发行仍需现场几何验收。',
        '- `object_samples.csv.gz`、`object_field_intervals.csv.gz`、`object_days.csv.gz`、`objects.csv`：动态主场、留边和可见性；重叠全部成员可用查询API获得。',
        '- `legacy_to_fixed.csv`、`seed_vs_final_objects.csv`、`neighbors.csv`、`special_targets.csv`：对照与完整性。',
        '- `row_geometry_audit.csv`、`wcs_*`、`measured_wcs_boundaries.json`、`independent_audit.json`：名义/实测情景分层证据。',
        '- `figures/01`—`06`：全年、两云局部、边缘、今晚对照、跨场图，各有PNG/PDF。',
        '- `provenance.json` 与 `MANIFEST_SHA256.json`：输入与交付哈希。原522场草案保存在原交接包，不混用RC场号。','',
        '```bash','python3 -m targeted_observing.fixed_field_build','python3 -m targeted_observing.fixed_field_wcs','python3 -m targeted_observing.fixed_trojan','python3 -m targeted_observing.fixed_field_audit','python3 -m targeted_observing.fixed_field_report','python3 -m unittest targeted_observing.tests.test_fixed_fields -v','```','',
        '完整重算需本地已冻结的年度缓存；GitHub保留代码、候选发行数据、来源/验收和紧凑图表，原始逐时产物留本地交付包。正式冻结前的下一步是提供科学遮罩/模板公共区域并实测跨夜指向与PA；六模块重构、overhead重拟合、晨昏实测曝光与后随升级仅记入PROJECT_MEMORY，不在本次实施。','']
    (output/'REPORT.md').write_text('\n'.join(lines))

if __name__=='__main__':draw();report()
