"""Build curated annual science tables and figures from completed local runs.

Raw caches/plans remain under products; only these compact report artifacts are
versioned. Re-running this module requires the final simulations to be complete.
"""
import json
import hashlib
import platform
from importlib.metadata import version
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import healpy as hp
from .catalog import BASE
from .annual_sources import dates
from .annual_coverage import pixels
from .geometry import Sky, setup_iers, coord, unix

ROOT=BASE/'products/annual'
OUT=BASE/'docs/annual_reports'
PRIMARY='queue_4fields_dwell2_baseline_science25_v3'


def read(path):return json.loads(path.read_text())


def save_table(name,frame):
    frame.to_csv(OUT/f'{name}.csv',index=False)
    return frame


def table(frame):
    def cell(x):
        if isinstance(x,float):return f'{x:.4f}' if 0<abs(x)<.01 else f'{x:.2f}'
        return str(x)
    rows=['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |']
    rows+=['| '+' | '.join(cell(x) for x in row)+' |' for row in frame.itertuples(index=False,name=None)]
    return '\n'.join(rows)


def plot_save(fig,name):
    fig.savefig(OUT/f'{name}.png',dpi=160,bbox_inches='tight');plt.close(fig)


def provenance():
    spec=read(BASE/'annual_defaults.json')
    files=[BASE/'annual_defaults.json',BASE/'defaults.json',BASE/'data/weather_prior.json',
           BASE/'cache/trojans.json',BASE/'cache/tess2025.tsv',BASE/'cache'/spec['iers_file'],
           ROOT/'input_verification.json',ROOT/'screening_manifest.json',ROOT/'selected_nominal/verification.json',
           ROOT/'selected_nominal/exposures.csv']
    files+=list(BASE.glob('annual_*.py'))
    files+=list(ROOT.glob('queue_*_science25_v3/summary.json'))
    files+=list(ROOT.glob('queue_4fields_dwell2_*_science25_v3/period_yield_summary.json'))
    result={'scope':'frozen input/code/product fingerprints for this local annual assessment',
        'python':platform.python_version(),'packages':{k:version(k) for k in ['numpy','scipy','pandas','astropy','astroquery','healpy','matplotlib']},
        'configuration':spec,'sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    (OUT/'provenance.json').write_text(json.dumps(result,indent=2)+'\n')


def strategy_tables():
    summaries=[];cadences=[];period_comparison=[]
    for n in [2,3,4,6]:
        for d in [1,2]:
            for s in ['conservative','baseline','optimistic']:
                root=ROOT/f'queue_{n}fields_dwell{d}_{s}_science25_v3';x=read(root/'summary.json')
                assert x['nights']==365 and x['realizations']==200 and x['queue_policy_version']==3
                ledger=pd.read_csv(root/'nightly_execution_ledger.csv')
                assert len(ledger)==73000 and ledger.groupby('realization').size().eq(365).all()
                assert (ledger.groupby('realization').quality_objects_cumulative.diff().dropna()>=0).all()
                used=ledger.allocated_s>0
                parts=['exposure_s','readout_s','slew_s','filter_s','delay_s','handover_s','slack_s']
                assert np.max(abs(ledger.loc[used,parts].sum(axis=1)-ledger.loc[used,'allocated_s']))<1e-5
                assert ledger.slack_s.min()>=-1e-5 and ledger.allocated_s.max()<=7200
                objs=pd.read_csv(root/'object_yield.csv')
                counts=objs.groupby('realization').quality_lightcurve.sum().reindex(range(200),fill_value=0)
                np.testing.assert_array_equal(counts.to_numpy(),x['counts'])
                hours=ledger.groupby('realization').allocated_s.sum()/3600
                p10,p50,p90=x['quality_count_P10_P50_P90']
                summaries.append({'fields':n,'dwell':d,'scenario':s,'P10':p10,'P50':p50,'P90':p90,
                                  'P_ge101':x['probability_ge101'],'allocated_hours_P50':hours.median()})
            rows=pd.read_csv(ROOT/f'queue_{n}fields_dwell{d}_baseline_science25_v3/execution_realization0.csv')
            gaps=rows.groupby(['date','field_id']).mid_unix.diff().dropna()
            points=rows.groupby(['date','field_id']).size()
            cadences.append({'fields':n,'dwell':d,'points_per_field_P10':points.quantile(.1),
                'points_per_field_P50':points.median(),'points_per_field_P90':points.quantile(.9),
                'gap_s_min':gaps.min(),'gap_s_P10':gaps.quantile(.1),'gap_s_P50':gaps.median(),
                'gap_s_P90':gaps.quantile(.9),'gap_s_max':gaps.max()})
            grid=pd.read_csv(ROOT/f'queue_{n}fields_dwell{d}_baseline_science25_v3/period_injection.csv')
            for mode,g in grid.groupby('mode'):
                period_comparison.append({'fields':n,'dwell':d,'mode':mode,'trials':len(g),
                    'unique_correct_fraction':g.period_unique_recovered.mean(),'false_unique_fraction':g.false_unique.mean(),
                    'amplitude_fraction':g.amplitude_recovered.mean()})
    summary=save_table('strategy_comparison',pd.DataFrame(summaries));save_table('cadence',pd.DataFrame(cadences))
    save_table('period_cadence_comparison',pd.DataFrame(period_comparison))
    fig,ax=plt.subplots(figsize=(10,4.6));colors=['#aa493b','#2676b8','#288e68']
    for j,s in enumerate(['conservative','baseline','optimistic']):
        f=summary[summary.scenario==s];x=np.arange(8)+(j-1)*.2
        ax.errorbar(x,f.P50,yerr=[f.P50-f.P10,f.P90-f.P50],fmt='o',capsize=3,label=s,color=colors[j])
    ax.axhline(101,color='black',ls='--',label='science target: 101')
    ax.set_xticks(np.arange(8),[f'{n} fields\ndwell {d}' for n in [2,3,4,6] for d in [1,2]])
    ax.set(ylabel='Unique quality lightcurves / year',title='200 shared weather realizations · P10 / P50 / P90 · 2 h maximum')
    ax.legend(ncol=4,loc='upper center',bbox_to_anchor=(.5,1.23));ax.grid(axis='y',alpha=.2);fig.tight_layout()
    plot_save(fig,'trojan_yield')
    sens=[]
    for folder in [PRIMARY,'queue_4fields_dwell2_optimistic_science25_v3','queue_4fields_dwell2_precision_requirement_science25_v3']+[f'queue_8fields_dwell2_{s}_budget14400_science25_v3' for s in ['baseline','optimistic','precision_requirement']]:
        x=read(ROOT/folder/'summary.json');sens.append({'fields':4 if '4fields' in folder else 8,'hours_per_block':x['allocation_s']/3600,
            'sigma_r18':x['scenario']['sigma_r18_120s'],'P10':x['quality_count_P10_P50_P90'][0],
            'P50':x['quality_count_P10_P50_P90'][1],'P90':x['quality_count_P10_P50_P90'][2],'P_ge101':x['probability_ge101']})
    save_table('budget_precision',pd.DataFrame(sens))
    return summary,pd.DataFrame(cadences),pd.DataFrame(sens)


def trojan_figures():
    ledger=pd.read_csv(ROOT/PRIMARY/'nightly_execution_ledger.csv')
    array=ledger.pivot(index='realization',columns='date',values='quality_objects_cumulative').to_numpy()
    q=np.percentile(array,[10,50,90],axis=0);frame=pd.DataFrame({'date':sorted(ledger.date.unique()),'P10':q[0],'P50':q[1],'P90':q[2]})
    save_table('cumulative_quality',frame)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    axes[0].fill_between(np.arange(365),q[0],q[2],alpha=.2,color='#2676b8');axes[0].plot(q[1],color='#2676b8')
    axes[0].set(xlabel='Night since 2026-11-01',ylabel='Completed unique objects',title='4 fields × dwell 2 · baseline')
    objects=pd.read_csv(ROOT/PRIMARY/'object_yield.csv');first=objects[objects.realization==0]
    bins=np.arange(14,20.6,.4)
    axes[1].hist(first.median_V_proxy,bins=bins,label='sampled',alpha=.4)
    axes[1].hist(first[first.quality_lightcurve].median_V_proxy,bins=bins,label='quality LC',alpha=.7)
    axes[1].set(xlabel='Median predicted V (r=V scenario)',ylabel='Objects',title='Weather realization 0');axes[1].legend()
    fig.tight_layout();plot_save(fig,'trojan_progress')
    roles=read(ROOT/PRIMARY/'science_roles.json');objects['role']=objects.object_id.astype(str).map(lambda k:roles[k]['role'])
    role_rows=[]
    for role in sorted(objects.role.unique()):
        counts=objects[(objects.role==role)&objects.quality_lightcurve].groupby('realization').size().reindex(range(200),fill_value=0)
        role_rows.append({'role':role,'P10':counts.quantile(.1),'P50':counts.median(),'P90':counts.quantile(.9)})
    save_table('science_role_yield',pd.DataFrame(role_rows))
    # Every sampled object remains visible in the compact annual index, including zero-quality objects.
    compact=objects.groupby('object_id').agg(simulations_sampled=('realization','size'),quality_realizations=('quality_lightcurve','sum'),
        median_V=('median_V_proxy','median'),median_good_points=('good_points','median'),median_good_nights=('good_nights','median'),
        median_baseline_days=('baseline_days','median'),median_error_mag=('median_random_error_mag','median')).reset_index()
    compact['quality_probability']=compact.quality_realizations/200
    compact['science_role']=compact.object_id.astype(str).map(lambda k:roles[k]['role'])
    save_table('object_summary',compact)
    grid=pd.read_csv(ROOT/PRIMARY/'period_injection.csv')
    grouped=grid.groupby(['mode','shape','period_injected_h','amplitude_injected_mag']).agg(
        trials=('object_id','size'),unique_correct=('period_unique_recovered','mean'),false_unique=('false_unique','mean'),amplitude_correct=('amplitude_recovered','mean')).reset_index()
    save_table('period_grid',grouped)
    fig,axes=plt.subplots(1,3,figsize=(12,3.7),sharey=True)
    for ax,shape in zip(axes,['asymmetric','symmetric','complex']):
        g=grid[(grid['shape']==shape)&(grid['mode']=='blind')].pivot_table(index='amplitude_injected_mag',columns='period_injected_h',values='period_unique_recovered')
        im=ax.imshow(g,origin='lower',aspect='auto',vmin=0,vmax=1,cmap='viridis')
        ax.set_xticks(range(len(g.columns)),[f'{p:g}' for p in g.columns],rotation=60)
        ax.set_yticks(range(3),[.1,.3,.6]);ax.set(title=shape,xlabel='Injected physical period (h)')
    axes[0].set_ylabel('Peak-to-peak amplitude (mag)');fig.colorbar(im,ax=list(axes),label='Blind correct + unique fraction',shrink=.8)
    plot_save(fig,'period_recovery')
    period=[]
    for scenario in ['baseline','optimistic']:
        folder=ROOT/f'queue_4fields_dwell2_{scenario}_science25_v3';x=read(folder/'period_yield_summary.json')
        assert 'mixture' in x['population']
        p=pd.read_csv(folder/'period_yield.csv')
        for mode,g in p.groupby('mode'):
            counts=g.groupby('realization').period_unique_recovered.sum().reindex(range(200),fill_value=0)
            period.append({'scenario':scenario,'mode':mode,'P10':counts.quantile(.1),'P50':counts.median(),'P90':counts.quantile(.9),
                'false_unique_fraction':g.false_unique.mean(),'amplitude_correct_fraction':g.amplitude_recovered.mean()})
    save_table('period_yield',pd.DataFrame(period))
    return pd.DataFrame(role_rows),pd.DataFrame(period)


def monthly_fields():
    spec=read(BASE/'annual_defaults.json');days=dates(spec['start_date'],spec['stop_date_exclusive'])
    audit=pd.read_csv(ROOT/'candidate_audit.csv',dtype={'object_id':str}).set_index('object_id')
    monthly=defaultdict(lambda:{'objects':set(),'fields':set(),'pointings':{},'eligible_nights':0})
    for day in days:
        menu=read(ROOT/'visibility'/f'{day}.json');data=np.load(ROOT/'visibility'/f'{day}.npz')
        active=np.flatnonzero(np.any(data['valid'],axis=1));m=monthly[day[:7]]
        m['eligible_nights']+=bool(len(active))
        for i in active:
            f=menu['fields'][i];m['objects'].update(f['core_targets']);m['fields'].add(f['field_id'])
        p=read(ROOT/'selected_nominal'/f'{day}.json')
        for b in p['blocks']:
            if b['block_id']=='TJ':
                for f in b['fields']:m['pointings'][f['field_id']]=f
    rows=[];fig,axes=plt.subplots(3,4,figsize=(13,8),sharex=True,sharey=True)
    for ax,(month,m) in zip(axes.flat,monthly.items()):
        ids=list(m['objects']);counts=audit.loc[ids].cloud_proxy.value_counts() if ids else {}
        fs=list(m['pointings'].values());ax.scatter([f['ra_deg'] for f in fs],[f['dec_deg'] for f in fs],s=12,color='#2676b8')
        ax.set(title=month,xlim=(0,360),ylim=(-35,55));ax.grid(alpha=.15)
        rows.append({'month':month,'eligible_nights':m['eligible_nights'],'eligible_objects':len(ids),
            'L4_proxy':counts.get('L4',0),'L5_proxy':counts.get('L5',0),'selected_field_versions':len(fs)})
    fig.supxlabel('RA (deg)');fig.supylabel('Dec (deg)');fig.suptitle('Selected nominal pointings · 35-day campaigns / 7-day versions');fig.tight_layout()
    plot_save(fig,'monthly_pointings');monthly_table=save_table('monthly_visibility',pd.DataFrame(rows))
    nominal=pd.read_csv(ROOT/'selected_nominal/time_ledger.csv');nominal['month']=nominal.date.str[:7]
    hours=nominal.groupby('month')[['EV_allocated_s','TJ_allocated_s','AM_allocated_s','free_s']].sum()/3600
    save_table('monthly_nominal_hours',hours.rename(columns=lambda x:x.removesuffix('_s')+'_h').reset_index())
    fig,ax=plt.subplots(figsize=(10,4));hours.plot.bar(stacked=True,ax=ax,color=['#e29e4d','#2676b8','#ad679d','#d6dce3'])
    ax.set(ylabel='Hours',xlabel='',title='Nominal allocation within Sun < -12° night · weather losses not deducted')
    ax.legend(['Evening','Trojan','Morning','Unallocated'],ncol=4);fig.tight_layout();plot_save(fig,'monthly_time')
    return monthly_table,hours


def twilight_tables():
    spec=read(BASE/'annual_defaults.json');setup_iers(BASE/'cache'/spec['iers_file'])
    nights=dates(spec['start_date'],spec['stop_date_exclusive']);area_rows=[];bins=[];moving=[];convergence=[]
    for visits in [5,10]:
        root=ROOT/f'twilight_{visits}';summary=read(root/'coverage_nside1024.json');low=read(root/'coverage_nside512.json')
        assert summary['night_count']==low['night_count']==365
        area_rows.append({'visits':visits,**{k:v for k,v in summary.items() if isinstance(v,(int,float)) and k not in ['nside','night_count']}})
        convergence.append({'visits':visits,'unique_relative_difference':abs(low['unique_year_deg2']/summary['unique_year_deg2']-1),
                            'qualified_relative_difference':abs(low['qualified_short_arc_unique_year_deg2']/summary['qualified_short_arc_unique_year_deg2']-1)})
        groups=defaultdict(set)
        for day in nights:
            p=read(root/f'{day}.json');sky=Sky(p['config'])
            rows=[r for b in p['blocks'] for r in b['rows']]
            if not rows:continue
            ts=np.array([unix(r['mid_utc']) for r in rows]);frame,sun,_,_=sky.environment(ts)
            centers=coord([r['ra_deg'] for r in rows],[r['dec_deg'] for r in rows]).transform_to(frame)
            elong=centers.separation(sun).deg
            for r,e in zip(rows,elong):
                label='45–55' if e<55 else '55–75' if e<75 else '75–90'
                groups[(day[:7],r['block_id'][:2],label)].update(pixels(r['ra_deg'],r['dec_deg'],1024).tolist())
        for (month,side,label),pix in groups.items():
            bins.append({'visits':visits,'month':month,'side':side,'elongation_deg':label,'unique_area_deg2':len(pix)*hp.nside2pixarea(1024,degrees=True)})
        for depth in [19,20,21]:
            df=pd.read_csv(root/f'moving_injections_causal_r{depth}.csv',
                dtype={'tracklet':'boolean','recovered_plus1':'boolean','recovered_plus3':'boolean'})
            for (speed,mag),g in df.groupby(['speed_deg_day','magnitude']):
                first=g[g.tracklet];plus1=g[g.recovered_plus1.notna()];plus3=g[g.recovered_plus3.notna()]
                moving.append({'visits':visits,'limiting_mag':depth,'speed_deg_day':speed,'magnitude':mag,'injected':len(g),
                    'tracklet_fraction':g.tracklet.mean(),'plus1_joint_fraction':plus1.recovered_plus1.mean(),
                    'plus3_joint_fraction':plus3.recovered_plus3.mean(),
                    'plus1_given_tracklet':first.recovered_plus1.mean(),'plus3_given_tracklet':first.recovered_plus3.mean(),
                    'within_existing_rate_range':g.within_existing_rate_range.iloc[0]})
    a=save_table('twilight_area',pd.DataFrame(area_rows));b=save_table('twilight_monthly_elongation',pd.DataFrame(bins))
    m=pd.DataFrame(moving)
    for key in ['tracklet_fraction','plus1_joint_fraction','plus3_joint_fraction','plus1_given_tracklet','plus3_given_tracklet']:
        m[key]=pd.to_numeric(m[key],errors='coerce').astype(float)
    save_table('moving_recovery',m);save_table('healpix_convergence',pd.DataFrame(convergence))
    fig,axes=plt.subplots(2,2,figsize=(12,6),sharey=True)
    for i,v in enumerate([5,10]):
        for j,side in enumerate(['EV','AM']):
            g=b[(b.visits==v)&(b.side==side)].pivot(index='month',columns='elongation_deg',values='unique_area_deg2').fillna(0)
            g.reindex(columns=['45–55','55–75','75–90'],fill_value=0).plot.bar(stacked=False,ax=axes[i,j],color=['#2a9d8f','#e9c46a','#e76f51'])
            axes[i,j].set(title=f'{v} visits / field · {side}',xlabel='',ylabel='Monthly unique footprint (deg²)')
            axes[i,j].tick_params(axis='x',rotation=60);axes[i,j].legend(fontsize=8)
    fig.tight_layout();plot_save(fig,'twilight_elongation')
    plot_moving(m)
    return a,b,m


def plot_moving(m):
    fig,axes=plt.subplots(2,3,figsize=(11,6),sharex=True,sharey=True,layout='constrained')
    for i,v in enumerate([5,10]):
        for j,key in enumerate(['tracklet_fraction','plus1_given_tracklet','plus3_given_tracklet']):
            g=m[(m.visits==v)&(m.limiting_mag==20)].pivot(index='magnitude',columns='speed_deg_day',values=key)
            im=axes[i,j].imshow(g,origin='lower',aspect='auto',vmin=0,vmax=1,cmap='viridis')
            axes[i,j].set_xticks(range(5),g.columns);axes[i,j].set_yticks(range(5),g.index)
            axes[i,j].set(title=f'{v} visits · '+['short arc','+1 | short arc','+3 | short arc'][j],
                          xlabel='Angular speed (deg/day)' if i==1 else '')
        axes[i,0].set_ylabel('Injected r magnitude')
    fig.colorbar(im,ax=list(axes.flat),label='Conditional fraction',shrink=.75)
    fig.suptitle('Dark r30 depth = 20 mag scenario · white: no initial short arc')
    plot_save(fig,'moving_recovery')


def write_reports(payload,summary,cadence,sensitivity,roles,period,monthly,hours,area,moving):
    grid=pd.read_csv(OUT/'period_grid.csv')
    effective=[{'visits':v,**read(ROOT/f'twilight_{v}/effective_coverage_realization0.json')} for v in [5,10]]
    effective_table=pd.DataFrame([{k:r[k] for k in ['visits','unique_year_deg2','qualified_short_arc_unique_year_deg2','area_on_ge2_nights_deg2','area_on_ge3_nights_deg2']} for r in effective])
    save_table('twilight_effective_area',effective_table)
    period_table=period.rename(columns={'false_unique_fraction':'错误但宣称唯一的比例','amplitude_correct_fraction':'振幅恢复比例'})
    source_text='''- 完整需求正文：[年度要求归档](../annual_requirements_pro_20261104.md)。其额外附件未下载；实现依据是完整正文及已归档的原四份附件。
- [JPL SBDB Query API](https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html) 与 [Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html)：全量 TJN 分类与兴隆站心星历。
- [TESS 2025 重分析，A&A 693 A66](https://arxiv.org/abs/2501.07189)：保留周期区间、备选解和质量信息。完整表已取得，9912行，不再列为缺失输入。
- [Li 等 2025，Icarus 438, 116609](https://doi.org/10.1016/j.icarus.2025.116609)：216颗 Trojan 周期，其中80颗此前无已知周期；不能把目录未命中自动当作新周期。
- [Huang 等 2026，KMTNet](https://arxiv.org/abs/2606.26184)：长期高重复采样的科学参考，三站和目标样本不同，不能按面积换算 SHARP 产额。
- [兴隆历史天气研究，MNRAS 511, 5363](https://academic.oup.com/mnras/article/511/4/5363/6543517)：使用已发表月度可观测率；可观测夜不等于高精度测光夜。逐年数字见 `data/weather_prior.json`。
- [Durech 与 Hanus 2026 模型](https://arxiv.org/abs/2610.06082)、DAMIT、冻结 LCDB2023、TESS DR1 与人工文献审查共同保留；并非最新文献的完备普查。'''
    trojan=f'''# SHARP 木星特洛伊年度监测方案

版本：2026-11-04；规划区间：2026-11-01—2027-10-31（北京时间夜次，内部UTC）。状态：本地模拟与项目申请材料，尚未现场执行。

## 结论与拟申请模式

建议将 **每夜最多两小时、r120、3—4场、每场连续两张** 作为年度候选主模式。当前交付的365夜完整曝光表采用四场驻场两张；三场同样值得在前三个有效夜实测后比较。两者基准年度合格光变的中位数均约42颗，四场在乐观精度情景下更高。没有为了达到百颗而降低质量门槛。

**两小时方案尚不能支撑“超过100颗高质量光变”的申请承诺。** 四场基准 P10/P50/P90 为36/42/47，乐观为55/62.5/70。全部24组标准对照均为0/200次达到101颗；这是本模型下的经验频率，不能解释成真实概率严格为零。二项抽样的一侧95%上限约1.49%，也不包含测光、天气和调度模型误差。

四小时、八场且假设r18单张随机误差约0.01mag时，P10/P50/P90为102.9/117/128.1，达到101的频率92%；这同时要求更高性能与更多时间，是追加时间的敏感性实验。尚未证明本仪器达到该精度，也未将四小时写入正常执行模式。

## 科学目的和独立计数

目标是可公开复用、支持自转研究的多夜时序，重点包括弱周期/别名检验、不同观测几何下的振幅变化、慢转候选和可靠控制源。主带对象只能是额外成果，未计入这里的 Trojan 数量。Lucy 的科学背景不等于这些场直接承担任务保障。

合格光变必须同时满足：≥100合格r点、≥10有效夜、每有效夜≥5点且跨度≥1800秒、基线≥21日；清洁点误差中位数≤0.03mag、至少80%≤0.05mag；单点夜零点偏移≤0.02mag，并检查按夜加权的零点RMS。严重极亮星场先排除，其余污染/坏点采用显式损失情景。所有数值是模拟判据，真实图像上的污染、饱和和系统误差仍须验证。

光变合格、周期正确且唯一、振幅恢复三者分别统计。对象第一次通过时保存合格前缀和完成证书，后续偶然落入场内的低质量帧不会删除已有合格数据；所有原始模拟点仍保留。如真实校准被推翻，必须另行重评，不能沿用模拟证书。

## 输入、季节与场版本

全量16445行JPL TJN母表逐行保留筛选结果；每周二体粗筛、H<14、V代理≤19.7、Dec≥−29°后得到1063颗。所有1063颗均有6小时采样的站心Horizons星历，覆盖全年及末端+3日恢复。12个独立逐小时查询的插值抽验最大位置差0.511角秒；不是对所有目标的逐点保证。缺H和粗筛拒绝项可审查，不能宣称人口样本完备。两云标签是瞬时黄经代理，不替代长期共振积分。

场组按35日科学campaign管理，场中心每7日产生新版本；跨版本按对象ID积累点数、有效夜、基线和完成状态。全足迹高度>30°，中心≥32.5°；月在地平线上时普通30°、照明≥75%时40°避让。可见菜单采用300秒网格并加时间角速度余量，实际名义曝光独立5秒复验。无共同窗口时正式输出零计划并释放时间。

{table(monthly)}

表中的“eligible_objects”是该月至少一个候选场有完整两小时合法窗口的几何对象数，不是光变产额；不同月不可相加。L4/L5是上述代理标签。可见季由星历计算，不强制全年每月占用时间。

![月度选中天区](monthly_pointings.png)

每场保留主目标与伴随目标。在每个场版本的实际计划访问上检查80%留场，未达到者标为secondary而非删除曝光记录。年内同一对象跨场去重。V<2恒星距场中心4.122°内的场采用保守排除；Hipparcos V<6补足Gaia亮端，参宿四不会再因Gaia缺失漏检。其余亮星、拥挤、CCD方向和坏区尚没有完整实测模型。

## 调度、cadence与时间账本

转场使用已有header正常分支模型并加2秒余量，与20/60秒预算取大；同场零转场，读出15秒，真正换滤镜10秒，进入块按未知前一指向60秒，末尾另留300秒交接。只有完整轮巡进入计划。额外延迟、读出/转场抖动及预约取消都进入年度执行账本。

{table(cadence)}

上表来自基准天气实现0的实际曝光，按场逐夜统计。同场连拍会产生约135秒的短间隔和较长回访间隔，不能用一个平均cadence代替。不同场组可能相距很远，因此这里的长间隔常比最初165秒/帧的手算大。周期对照见 `period_cadence_comparison.csv`；除主方案24对象×2相位外，其他策略使用8对象×1相位、每个都有三种形状，样本构成不同，仅作诊断，不能当作严格同目标的因果排序。

全年名义预留：Trojan {hours.TJ_allocated_s.sum():.1f}小时，黄昏 {hours.EV_allocated_s.sum():.1f}小时，清晨 {hours.AM_allocated_s.sum():.1f}小时；−12°夜界内尚未分配 {hours.free_s.sum():.1f}小时。预留时间包含交接和空隙，不等于曝光积分时间，也不保证当夜可执行。

![月度时间](monthly_time.png)

## 共同天气与测光情景

200条全年天气序列供各策略共同使用，随机种子20261104。月可观测率来自历史资料；潜变量AR(1)相关系数0.55、部分云概率0.25、外部取消概率0.15、最大延时180秒、额外读出均值0.5秒、额外转场均值2秒、seeing中位2角秒均是规划情景。当前天气关闭可知，未来随机天气不可用于选场；仅使用已发生反馈、可预知几何与历史月天气先验。不同策略共用潜在扰动，不代表按真实绝对时刻重放同一片云。

三种精度以120秒、r18为参考：保守随机误差0.040/夜零点0.020/有效比例0.70；基准0.025/0.012/0.85；乐观0.015/0.008/0.95。逐点误差再依亮度、气团、月光和seeing变化，含0.005mag底噪。V暂作为r；颜色偏移±0.3mag可使源/背景主导项约改变15%—32%，不能用白光20.8极限替代r精度。精度情景覆盖较宽但不等价于经过校准的颜色误差传播。

{table(summary)}

![年度光变产出](trojan_yield.png)
![累计完成和亮度](trojan_progress.png)

科学评分采用测光可用性×(0.75+0.25×文献权重)×剩余完成机会及已观测进度。达到点数/夜数但尚未满21日时降低占用、到期补齐；完成后不再获取主目标优先权。原始文献质量保持不变：2223是控制源，22180保留别名复核，未查到证据保持literature_unchecked。

{table(roles)}

各角色分位数不能相加得到总分位数。未核查对象不计为“首次周期”。年度逐对象表 `object_summary.csv` 给出采样和质量概率，完整200实现逐对象明细仍在本地产品目录。

## 周期与振幅恢复

实际采样注入2—1000小时、0.1/0.3/0.6mag、随机相位，以及非对称双峰、完全对称双峰、含四次谐波三种形状。主网格8640次拟合含盲搜和±5%外部先验。另对所有200天气实现中合格对象进行总体情景注入：可靠目录/模型周期优先，其余对象固定抽取对数均匀周期，三形状权重0.6/0.2/0.2。这是透明的示例人口先验，并非真实自转分布测量。

![周期恢复](period_recovery.png)

{table(period_table)}

可靠恢复同时要求与真值相差<1%、竞争频率Δχ²≥9、周期移动±1%重新拟合后Δχ²≥9。后两项只是候选集合内的启发式，不是校准后的后验概率。对称双峰可能只能识别半周期，因此显式报告“错误但宣称唯一”的比例。算法先筛12个一阶峰及半/倍频，再局部二阶拟合，不能替代最终多阶Fourier/别名人工审查。振幅容差为max(0.05mag,20%)。

逐对象周期表另保存20相位箱覆盖率与最大相位空隙，均针对声明的注入周期，不把未知真周期的相位覆盖视为实测。总体注入采用Horizons光行时修正；代表网格以实际接收采样时刻同时生成与拟合信号。假设距离/相角已理想校正、形状和自转恒定，没有逐夜去均值。误差由平均亮度预测，未把自转亮度变化反馈到探测截断。先验和盲搜结果相互重叠，不能相加；这些也不是已经获得的周期或新发现数量。

## 百颗缺口与可执行决策

{table(sensitivity)}

标准两小时、四场在假定0.01mag性能下中位80颗，达101频率0.5%；只增加精度仍不足以稳健达到目标。四小时八场基准中位62，乐观96.5；只有额外精度实验中位117、概率92%，仍略低于设计中位120的愿望。正常调度器继续拒绝超过两小时，这些额外预算实验不能直接导出执行。

固定基准精度、全年无天气/预约/坏帧损失的同一策略名义实现也只有79颗，说明天气并非唯一瓶颈；这是该策略的对照，不是望远镜理论上限。主要限制包括可测亮对象密度、南天L4完整足迹窗口、月距、21日基线、≥10个有效夜、远转场和测光质量。尚未穷举所有场中心和全年全局最优解，不能据本结果证明其他调度不可能超过百颗。

建议申请文本先采用诚实的两小时情景区间，同时将“百颗”保留为需要额外资源与实测性能共同支持的升级目标。前三个有效夜固定r，测量r120误差—星等、夜零点、光晕、坏区、真实转场和有效点率；之后更新参数，重跑剩余年度，并按一项主要参数一次调整。无需等试观测才拥有完整年度设计。

## 交付、复现和边界

详见 [材料索引与复现](README.md)。365夜最终组合计划共{payload['verification']['exposures']}条曝光，逐条独立5秒几何检查通过；预算与未占用时间闭合。全套35项回归测试通过。当前所有产出都是本地模拟，不含控制器合法token、真实图像的测光认证或服务器部署。原生产目录与MPC暂停状态未改。

{source_text}
'''
    (OUT/'TROJAN_ANNUAL_REPORT.md').write_text(trojan)
    compact_area=area[['visits','exposure_area_deg2','unique_year_deg2','qualified_short_arc_unique_year_deg2','area_on_ge2_nights_deg2','area_on_ge3_nights_deg2']]
    compact_moving=moving[(moving.limiting_mag==20)&(moving.magnitude==18)][['visits','speed_deg_day','tracklet_fraction','plus1_given_tracklet','plus3_given_tracklet','within_existing_rate_range']]
    twilight=f'''# SHARP 晨昏年度观测方案

版本：2026-11-04；规划区间2026-11-01—2027-10-31，共365夜、730个晨昏窗口。状态：本地设计与条件效率模拟。

## 方案与科学范围

保留 **r30、每侧最多一小时、两组各四场五访** 为覆盖基线，同时交付一组四场十访的全年对照。晚窗从太阳降到−12°开始，晨窗在太阳升到−12°结束。这一小时包含部分暗夜，不等同于−12°到−18°的天文暮光。高度门槛使极小距角不可达，因此报告把45—55°、55—75°与75—90°分别列出，后者称扩展搜索。

科学目标包括低距角小天体搜索、短弧确认与恢复机会的量化。没有使用归一化NEO轨道人口模型，**不报告年度绝对发现数**。当前生产处理速度约3—63角秒/小时，不因为天区已排好就自动具备高速NEO发现能力；快于这个范围的模拟只说明几何/测光潜力。

## 时间、几何与组间切换

读出15秒、换片10秒、转场采用实测正常分支下限加余量；晨昏相邻阈值3.5°并检查闭合回跳，未知进入指向预算60秒。全足迹高度>30°、中心≥32.5°；黄纬±25°、月距30/40°、太阳高度≤−12°。每组完整曝光才保留，剩余不足不塞入不完整组。晚上先小距角后大距角，早晨相反。

五访模式每侧最多40张，四场轮巡产生约20分钟弧段；十访也是最多40张，单组持续更久。最终间隔以逐曝光时间为准。十访对漏帧更稳健，但低高度场更难维持一个完整窗口。45°附近必要高度关系只是下限，五/十访全过程还须同时合法。

## 全年场与三有效夜锚点

名义年度每夜都有主组或明确零计划。天气队列使用当晚和过去14日的合法场组菜单，原锚点只有形成≥4次、跨度≥720秒的有效采样后才累积一个有效夜；达到三有效夜后优先新场。场组被局部修复时按field ID保留已完成有效夜，不能无限重置年龄。几何失效或主备不兼容时记录替换原因。

已有原11月6日丢锚点回归保留，年度版本进一步保存每个窗口的选中组、保留数和替换理由。队列没有提前使用未来随机天气；可预知几何是允许的。菜单只含历史/当前生成组，未穷举所有连续天区；真实待恢复候选及不确定区尚未接入，因此下述恢复只是扫描策略的机会，不是候选驱动后随的完整系统。

## 覆盖面积：必须区分四种口径

名义全年几何覆盖，HEALPix nside1024，约0.00328平方度/像素：

{table(compact_area)}

exposure_area是9平方度×曝光张数，仅是曝光面积之和。unique_year按球面足迹去重；qualified_short_arc采用缩进0.15°边界后的足迹，要求同窗口累计≥4访/≥12分钟；两夜/三夜是不同北京时间夜次重复。单月面积不能直接相加得到全年独立面积。nside512与1024的相对差异见 `healpix_convergence.csv`，有限像素化误差保留。

天气实现0的实际有效帧重放（与因果队列逐窗口good_exposures和qualified_fields逐项核对）：

{table(effective_table)}

这张表已扣除天气关闭、预约取消、云、迟到、几何失效和坏帧，尚未施加星等相关检出概率。它是一个共同天气实现的面积，不能当作200次面积分布的中位数。200次实现分别有146000条窗口记录；五访/十访的合格场窗口数P10/P50/P90为1922.5/2125.5/2336与1144.8/1256/1370.6，同一场不同夜重复计数，不能直接乘9当作全年面积。

![月度太阳距角覆盖](twilight_elongation.png)

图按每次曝光中点的场中心太阳距角分箱，每月/早晚/箱内去重足迹。一个场可能跨箱，所以不同箱也不能相加后声称完全独立面积。本次五访方案45—55°仅有少量覆盖，十访方案没有该箱覆盖；大部分属于55—75°及75—90°，不能将全数包装成极低距角巡天。完整表 `twilight_monthly_elongation.csv` 保留月份、晨/昏及策略；没有低距角合格场时允许只排扩展区或零计划。

## 移动源注入与恢复

速度取0.1/0.3/0.7/1.5/5度/日，星等17—21，每个已排窗口/速度/星等单元12个随机位置方向。注入位置在已选足迹内均匀分布，沿大圆恒角速移动三日；不是轨道人口模型，没有加速度或真实NEO轨道偏好。每个初始窗口轮换使用200天气实现之一，后续+1/+3日沿用同一实现及相应因果选场，不偷看未来随机天气。

检出概率采用logistic宽度0.25mag；r30暗夜深度分别19/20/21mag三情景，加入最长0.8mag暮光罚项、气团消光、月光、seeing及拖线损失。移动源检出未另加真实CCD坏区掩膜，也未生成图像或假源；这些额外损失需实测，因此效率仍偏理想化。窗口初始已关闭的情况单列在执行账本，未混入“给定已排窗口”的效率分母。末年+1/+3日超出申请区间的恢复右删失，不当失败计数。输入星历/IERS仍留有+3日余量。

以下为r18、暗夜深度20mag情景；所有数值均为条件比例：

{table(compact_moving)}

![移动目标恢复](moving_recovery.png)

tracklet_fraction要求同夜≥4次检出且跨度≥12分钟。plus1/3_given_tracklet的分母为初夜已成功短弧且后续未删失，后夜允许晨或昏任一窗口形成短弧。完整三深度×五星等×五速度结果见 `moving_recovery.csv`，另有以全注入为分母的联合成功率。

**同场复拍从未直接算作目标恢复。** 每次都推进目标位置并检查实际新足迹；5度/日移动天体通常迅速离开3度视场。这里的后续成功是利用注入真值关联的上限，没有假源、轨道歧义或实际链接算法，不能称为真实连接成功率。当前处理范围只包含表中的0.1、0.3度/日两档；0.7、1.5、5必须单列，未来改造检测程序还需要真实图像注入与误报评估。

## 选择与三夜校准

两组五访在同样40张上限下通常覆盖更广；十访可为较高漏检率提供更多冗余，但不一定提高跨夜恢复。申请先保留五访基线及十访对照，不能只看曝光数或某一星等成功率作决定。因未归一化人口，当前不能据效率表断言哪一种绝对NEO产量更高。

首三个有效夜保持r30和主要读出设置，实测天空背景随太阳高度/气团/月距变化、检出完整率、拖线损失、转场、坏区和tracklet误报；选择已知移动目标做真实再获取。然后以实际深度替换19/20/21情景，校准三有效夜更新规则和快速确认时效，重新生成余下年度。明亮NEO确认、彗星活动等机会任务可作为后续独立优先队列，当前未暗中占用YSO/SN/public项目时间。

## 交付、限制与来源

两种模式各365份名义夜计划/730窗口；200条因果执行序列；主方案与Trojan组合后365夜共{payload['verification']['exposures']}条曝光，独立5秒复验通过。各夜显式零计划和未分配区间可查。控制器token/定时接口、实际快速目标检测和现场发布尚未启用，未连接服务器。

所有表和复现命令见 [材料索引](README.md)。[Bolin等2024 ZTF晨昏研究](https://arxiv.org/abs/2409.15263)提供分组多访和后随的科学参考；不同地平限制不可直接移植其覆盖或发现率。共同天气来源及输入与Trojan报告相同：

{source_text}
'''
    (OUT/'TWILIGHT_ANNUAL_REPORT.md').write_text(twilight)
    index='''# 年度材料索引

- [Trojan 年度报告](TROJAN_ANNUAL_REPORT.md)
- [晨昏年度报告](TWILIGHT_ANNUAL_REPORT.md)
- [机器可读摘要](results.json)；其余同目录CSV和PNG由 `annual_report.py` 从完整运行结果生成。
- [年度目标](../GOAL.md)、[完整需求正文](../annual_requirements_pro_20261104.md)。

## 本地完整产品

这些是可重建运行产物，不提交Git。GitHub同步源码、配置、科学输入和本目录审查摘要，避免把大型缓存/日志加入仓库。

| 路径（相对 targeted_observing/products/annual） | 内容 |
|---|---|
| candidate_audit.csv / candidates.json | 全量16445行接受/拒绝理由、1063全年候选 |
| ephemeris_index.json / input_verification.json | 下载状态、原始星历hash、IERS与12历元抽验 |
| visibility/ | 365夜完整场菜单、版本、目标、共同可见窗口 |
| selected_nominal/ | 主策略365夜JSON、exposures.csv、time_ledger.csv、fields.json、retention.csv、独立verification.json |
| nominal_queue_*_science25_v3/ | 八种标准策略的全晴名义逐曝光表和365夜零计划账本 |
| queue_*_science25_v3/ | 各200天气实现的逐夜账本、逐对象产出；qualified_curves.npz保留误差/零点/光行时采样 |
| queue_4fields_dwell2_*_science25_v3/period_yield.csv | 全200天气的周期/振幅注入真值与结果，盲搜/先验分开 |
| twilight_5/ 与 twilight_10/ | 各365夜计划、历史备选菜单、causal_choices、146000行窗口账本、面积与移动源注入 |

`*_smoke`、旧的无science后缀队列、science v1/v2、fixed40诊断，以及无r深度后缀的旧moving_injections不属于本次最终结果。额外budget14400是申请敏感性，不是标准两小时执行计划。

## 复现顺序

在仓库根目录运行；需要Python、numpy/scipy/pandas/astropy/astroquery/requests/matplotlib/healpy。先运行基础fetch_sources/review_data取得母表与文献；年度计算另需annual_defaults.json所指定的IERS文件。当前本地已冻结iers_20261104.all；如在新机器获取，可从 https://maia.usno.navy.mil/ser7/finals2000A.all 下载并保存为新的带日期文件，确认MJD末端覆盖stop_date_exclusive再加4日，并更新配置，勿覆盖既有快照。后续用固定缓存和hash复跑；provenance.json记录版本与指纹。具体来源下载入口见 targeted_observing/README.md。默认日期、35日campaign、7日版本、200实现及随机种子全部在annual_defaults.json。

```bash
python3 -m targeted_observing.annual_sources --workers 2
python3 -m targeted_observing.annual_ephemeris_qa
python3 -m targeted_observing.annual_twilight --workers 4
python3 -m targeted_observing.annual_visibility --workers 4
python3 -m targeted_observing.annual_simulation --workers 4 --save-curves --scenario baseline optimistic conservative
python3 -m targeted_observing.annual_simulation --workers 4 --nominal --scenario baseline
python3 -m targeted_observing.annual_twilight_queue --visits 5 --workers 3
python3 -m targeted_observing.annual_twilight_queue --visits 10 --workers 3
python3 -m targeted_observing.annual_delivery
python3 -m targeted_observing.annual_effective_area
```

对两种visits分别运行annual_coverage的 `--causal --limiting-mag 19/20/21`；面积另用 `--nside 512 --skip-injections` 检查收敛。对八种baseline队列运行annual_periods，主四场驻场两张采用 `--objects 24 --phases 2`，其余 `--objects 8 --phases 1`。对主策略baseline与optimistic运行annual_period_yield。额外精度/时间实验分别用annual_simulation的 `--scenario precision_requirement --fields 4 --dwells 2` 和 `--duration-s 14400 --fields 8 --dwells 2 --scenario baseline optimistic precision_requirement`。全部完成后：

```bash
python3 -m targeted_observing.annual_report
python3 -m unittest discover -s targeted_observing/tests -v
```

报告生成器会检查24组365×200账本、时间闭合、年度去重质量数与累计完成单调性，缺最终结果即失败，不以早期诊断填补。两份报告包含模型假设、未建模项及未达百颗结论。没有任何远程发布命令。
'''
    (OUT/'README.md').write_text(index)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    verification=read(ROOT/'selected_nominal/verification.json')
    assert verification['passed'] and verification['nights']==365
    summary,cadence,sensitivity=strategy_tables();roles,period=trojan_figures()
    monthly,hours=monthly_fields();area,elong,moving=twilight_tables()
    payload={'strategy':summary.to_dict('records'),'cadence':cadence.to_dict('records'),'sensitivity':sensitivity.to_dict('records'),
             'roles':roles.to_dict('records'),'period':period.to_dict('records'),'monthly':monthly.to_dict('records'),
             'nominal_hours':{k.removesuffix('_s')+'_h':v for k,v in hours.sum().items()},'coverage':area.to_dict('records'),
             'verification':{k:v for k,v in verification.items() if k!='records'}}
    (OUT/'results.json').write_text(json.dumps(payload,indent=2)+'\n')
    write_reports(payload,summary,cadence,sensitivity,roles,period,monthly,hours,area,moving)
    provenance()
    print(json.dumps({k:v for k,v in payload.items() if k in ['nominal_hours','verification','period']},indent=2),flush=True)


if __name__=='__main__':main()
