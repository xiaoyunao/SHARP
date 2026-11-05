"""Actual scheduled-time sky coverage and same-weather twilight comparisons."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import AltAz,GeocentricMeanEcliptic,get_sun
from astropy.time import Time
import astropy.units as u
from .catalog import BASE,write_csv
from .geometry import Sky,setup_iers,coord,unix
from .annual_sources import dates
from .continuation_report import archive_csv,sha,OUT,ROOT


def read_injections(path):
    return pd.read_csv(path,dtype={k:'boolean' for k in ['tracklet','recovered_plus1','recovered_plus3']})


def geometry():
    spec=json.loads((BASE/'annual_defaults.json').read_text());setup_iers(BASE/'cache'/spec['iers_file'])
    nights=dates(spec['start_date'],spec['stop_date_exclusive']);allrows=[]
    for name in ['twilight_5','twilight_morning10_aligned']:
        rows=[]
        for day in nights:
            p=json.loads((ROOT/name/f'{day}.json').read_text());sky=Sky(p['config'])
            for b in p['blocks']:
                if not b['block_id'].startswith('AM'):continue
                for r in b['rows']:rows.append({'strategy':name,'date':day,'field_id':r['field_id'],'mid_unix':unix(r['mid_utc']),'ra_deg':r['ra_deg'],'dec_deg':r['dec_deg']})
        d=pd.DataFrame(rows);t=Time(d.mid_unix.to_numpy(),format='unix');c=coord(d.ra_deg.to_numpy(),d.dec_deg.to_numpy())
        aa=AltAz(obstime=t,location=sky.location,pressure=0*u.hPa)
        d['elongation_deg']=c.transform_to(aa).separation(get_sun(t).transform_to(aa)).deg
        d['beta_deg']=c.transform_to(GeocentricMeanEcliptic(equinox=Time('J2000'),obstime=t)).lat.deg
        allrows.append(d)
    data=pd.concat(allrows,ignore_index=True);archive_csv(data,OUT/'morning_scheduled_geometry.csv.gz')
    plot_twilight()


def plot_twilight():
    from .plot_style import BLUE, ORANGE, apply_style, save_pair
    apply_style()
    data = pd.read_csv(OUT/'morning_scheduled_geometry.csv.gz')
    xbins, ybins = np.arange(44, 92, 2), np.arange(-26, 28, 2)
    density = [np.histogram2d(g.elongation_deg, g.beta_deg, bins=[xbins,ybins])[0].T
               for _, g in data.groupby('strategy', sort=False)]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), sharey=True, layout='constrained')
    for ax, panel, label in zip(axes, density, ['AM5','AM10 · end aligned']):
        image = ax.pcolormesh(xbins, ybins, np.ma.masked_equal(panel,0), cmap='viridis',
                              vmin=0, vmax=max(p.max() for p in density))
        ax.set_title(label)
        ax.set_xlabel('Solar elongation [deg]')
        ax.grid()
    axes[0].set_ylabel('Ecliptic latitude [deg]')
    fig.colorbar(image, ax=axes, label='Exposure centres [count]')
    fig.suptitle('365-night nominal morning schedule')
    save_pair(fig, OUT/'morning_geometry')
    plt.close(fig)

    recovery = pd.read_csv(OUT/'twilight_conditional_recovery.csv')
    recovery = recovery.query("side == 'AM' and dark_depth == 20 and magnitude == 19")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.9), sharey=True, layout='constrained')
    measures = [('tracklet_fraction','Initial short arc'), ('plus1_conditional_fraction','+1 night'),
                ('plus3_conditional_fraction','+3 nights')]
    for ax,(key,title) in zip(axes, measures):
        for strategy,label,color,style in [('twilight_5','AM5',BLUE,'o-'),
                ('twilight_morning10_aligned','AM10 end aligned',ORANGE,'s--')]:
            group = recovery[recovery.strategy==strategy].sort_values('speed_deg_day')
            ax.plot(group.speed_deg_day, group[key]*100, style, color=color, label=label)
        ax.set_xscale('log')
        ax.set_xlabel('Motion [deg/day]')
        ax.set_title(title)
        ax.grid()
    axes[0].set_ylabel('Conditional efficiency [%]')
    axes[0].legend()
    fig.suptitle('Injected r=19 · dark r30 depth 20 · common 200-weather ensemble')
    save_pair(fig, OUT/'morning_recovery')
    plt.close(fig)


def summarize():
    rows=[];injections=[];sources={}
    for name in ['twilight_5','twilight_morning10_aligned']:
        root=ROOT/name;d=pd.read_csv(root/'causal_execution_ledger.csv')
        if len(d)!=146000:raise ValueError('Incomplete twilight causal ledger')
        nominal=json.loads((root/'coverage_nside1024.json').read_text())
        effective=json.loads((root/'effective_coverage_realization0.json').read_text())
        for side,g in d.groupby('side'):
            counts=g.groupby('realization').qualified_fields.sum()
            rows.append({'strategy':name,'side':side,'qualified_field_windows_p10':np.percentile(counts,10),
                'qualified_field_windows_p50':np.median(counts),'qualified_field_windows_p90':np.percentile(counts,90),
                'mean_good_exposures':g.groupby('realization').good_exposures.sum().mean(),
                'nominal_year_union_both_sides_deg2':nominal['unique_year_deg2'],
                'nominal_repeated_ge2_nights_deg2':nominal['area_on_ge2_nights_deg2'],
                'effective_year_union_realization0_deg2':effective['unique_year_deg2'],
                'qualified_plus1_area_sum_deg2_night_pairs':effective['qualified_area_day_plus_1_sum_deg2'],
                'qualified_plus3_area_sum_deg2_night_pairs':effective['qualified_area_day_plus_3_sum_deg2']})
        # Persist each realization's actual totals; identical weather IDs permit pairing.
        compact=d.groupby(['realization','side'])[['scheduled_exposures','good_exposures','qualified_fields','retained_anchors']].sum().reset_index()
        compact.to_csv(OUT/f'{name}_paired_window_totals.csv',index=False)
        for depth in [19,20,21]:
            p=root/f'moving_injections_causal_r{depth}.csv';x=read_injections(p)
            for (side,mag,speed),g in x.groupby(['side','magnitude','speed_deg_day']):
                # Nullable last-year followups are excluded from the denominator.
                follow={}
                for gap in [1,3]:
                    valid=g[f'recovered_plus{gap}'].notna();qualified=g.tracklet&valid
                    success=g.loc[qualified,f'recovered_plus{gap}'].astype(bool)
                    follow.update({f'plus{gap}_eligible_tracklets':int(qualified.sum()),f'plus{gap}_conditional_fraction':float(success.mean()) if len(success) else None})
                injections.append({'strategy':name,'side':side,'dark_depth':depth,'magnitude':mag,'speed_deg_day':speed,
                                   'injected':len(g),'tracklet_fraction':g.tracklet.mean(),**follow})
            sources[str(p.relative_to(BASE))]=sha(p)
        sources[str((root/'causal_execution_ledger.csv').relative_to(BASE))]=sha(root/'causal_execution_ledger.csv')
    write_csv(OUT/'twilight_comparison.csv',rows);write_csv(OUT/'twilight_conditional_recovery.csv',injections)
    (OUT/'twilight_provenance.json').write_text(json.dumps({'input_hashes':sources,
        'scope':'same weather ensemble; conditional truth-associated moving-source recovery, not NEO discovery rate',
        'geometry':'nominal actual exposure centres; separate HEALPix independent area estimates; no PNG-derived coordinates',
        'end_alignment':'30s nominal guard; simulated additive delays are not absorbed by waiting and may lose late visits'},indent=2)+'\n')
    # Common night time partition includes whole held twilight windows and TJ overhead.
    from .annual_weather import ensemble
    spec=json.loads((BASE/'annual_defaults.json').read_text());nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    weather=ensemble(nights,200,spec['seed']);total={}
    for day in nights:
        p=json.loads((ROOT/'twilight_5'/f'{day}.json').read_text());total[day]=p['night']['morning'][1]-p['night']['evening'][0]
    tw=pd.read_csv(ROOT/'twilight_5/causal_execution_ledger.csv')
    tw['held_s']=np.where(tw.scheduled_exposures>0,3600.,0.)
    held=tw.groupby(['realization','date']).held_s.sum();partitions=[]
    for label in 'ABCD':
        tj=pd.read_csv(ROOT/f'continuation_{label}_baseline/nightly_execution_ledger.csv')
        for r in tj.to_dict('records'):
            ri=r['realization'];ni=nights.index(r['date']);flags=weather['reserved'][ri,ni]
            external=float(np.dot(flags,[3600,7200,3600]));tw_s=held.loc[(ri,r['date'])]
            # External reservations are counted even on weather-closed nights.
            free=total[r['date']]-tw_s-r['allocated_s']-external
            if free<-.01:raise ValueError('Unified budget overrun')
            partitions.append({'strategy':label,'realization':ri,'date':r['date'],'night_s':total[r['date']],
                               'twilight_held_s':tw_s,'trojan_held_s':r['allocated_s'],'external_reserved_s':external,
                               'unallocated_s':free,'separate_calibration_s':0,
                               'trojan_handover_already_in_held_s':r['handover_s']})
    archive_csv(pd.DataFrame(partitions),OUT/'unified_time_ledger.csv.gz')


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--geometry-only',action='store_true');p.add_argument('--plots-only',action='store_true');a=p.parse_args()
    if a.plots_only:plot_twilight()
    elif a.geometry_only:geometry()
    else:summarize()
