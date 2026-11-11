"""Review a saved fixed night: original MP footprints, timing and sampling."""
import argparse
import hashlib
import json
import shutil
import zipfile
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import astropy.units as u
from astropy.coordinates import BarycentricMeanEcliptic, SkyCoord
from astropy.time import Time
from .catalog import BASE, write_csv
from .geometry import coord, setup_iers, local_time
from .offline_twilight import load_grid, make_sky, verify_window
from .plot_style import apply_style, save_pair


TZ = ZoneInfo('Asia/Shanghai')


def stamp(t):
    return datetime.fromtimestamp(t, TZ)


def report(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib.patches import Polygon

    plan = json.loads((source/'plan.json').read_text())
    output.mkdir(parents=True, exist_ok=True)
    setup_iers(BASE/'cache'/plan['config']['iers_file'])
    sky = make_sky()
    for block in plan['blocks']:
        result = verify_window(sky, block, plan['config'])
        if not result['passed']:
            raise ValueError(result)
    for name in ('EV.txt', 'AM.txt', 'START_TIMES.json', 'EV_timing.csv', 'AM_timing.csv', 'plan.json'):
        if source.resolve() != output.resolve():
            shutil.copyfile(source/name, output/name)
    all_fields, _ = load_grid()
    apply_style()
    field_fig, maps = plt.subplots(2, 2, figsize=(15, 11), layout='constrained')
    stats_fig, axes = plt.subplots(4, 2, figsize=(15, 15), layout='constrained')
    field_stats = []
    summaries = []
    ecliptic = SkyCoord(lon=np.linspace(0, 360, 1441)*u.deg, lat=np.zeros(1441)*u.deg,
                        frame=BarycentricMeanEcliptic(equinox=Time('J2000'))).icrs
    for column, block in enumerate(plan['blocks']):
        rows, fields, side = block['rows'], block['fields'], block['side']
        if not rows:
            for ax in [maps[0, column], maps[1, column], *axes[:, column]]:
                ax.text(.5, .5, side+': no feasible list', transform=ax.transAxes, ha='center')
            summaries.append({'side': side, 'frames': 0})
            continue
        colors = plt.get_cmap('tab10')(np.arange(len(fields)))
        reference = float(np.mean([f['ra_deg'] for f in fields]))
        wrap = lambda ra: (np.asarray(ra)-reference+180) % 360-180+reference
        ax = maps[0, column]
        xs = np.concatenate([wrap(f['corners_ra_deg']) for f in fields])
        ys = np.concatenate([f['corners_dec_deg'] for f in fields])
        ax.scatter(wrap([f['ra_deg'] for f in all_fields]), [f['dec_deg'] for f in all_fields],
                   s=7, color='.82', label='Original MP grid')
        ex = wrap(ecliptic.ra.deg)
        order = np.argsort(ex)
        ax.plot(ex[order], ecliptic.dec.deg[order], color='.4', ls='--', lw=1, label='Ecliptic J2000' if np.any((ecliptic.dec.deg>ys.min()-1)&(ecliptic.dec.deg<ys.max()+1.5)&(ex>xs.min()-1.1)&(ex<xs.max()+1.1)) else '_nolegend_')
        centers = np.array([[wrap(f['ra_deg']), f['dec_deg']] for f in fields])
        for i, (f, color) in enumerate(zip(fields, colors)):
            ax.add_patch(Polygon(np.c_[wrap(f['corners_ra_deg']), f['corners_dec_deg']],
                                 facecolor=color, edgecolor=color, alpha=.18))
            ax.scatter(*centers[i], color=color, s=35)
            ax.annotate(f"{i+1}: {f['field_id']}", centers[i], xytext=(3, 6),
                        textcoords='offset points', fontsize=10, color=color)
            following = (i+1) % len(fields)
            ax.annotate('', centers[following], centers[i],
                        arrowprops={'arrowstyle': '->', 'color': '.25', 'lw': 1.1,
                                    'linestyle': '--' if following == 0 else '-'})
        ax.set(xlim=(xs.max()+1.1, xs.min()-1.1), ylim=(ys.min()-1, ys.max()+1.5),
               xlabel='ICRS RA (deg; increasing left)', ylabel='ICRS Dec (deg)',
               title=f'{side}: original footprints / repeated route')
        ax.set_aspect(1/np.cos(np.deg2rad(np.mean(ys))))
        ax.legend(loc='lower left', fontsize=9)

        # All objects in this panel use one common epoch, not each field's visit time.
        snapshot = rows[0]['mid_unix']
        frame, sun, moon, fraction = sky.environment([snapshot])
        aa = coord([f['ra_deg'] for f in fields], [f['dec_deg'] for f in fields]).transform_to(frame[0])
        ax = maps[1, column]
        ax.axhspan(-25, 25, color='.93')
        ax.axhline(25, color='red', ls='--', lw=1, label='25 deg altitude')
        ax.scatter(sun.az.deg, sun.alt.deg, s=170, marker='*', color='goldenrod', label='Sun')
        ax.scatter(moon.az.deg, moon.alt.deg, s=70, marker='o', color='.45', label=f'Moon ({fraction[0]:.0%} lit)')
        for i, (f, color) in enumerate(zip(fields, colors)):
            fp = coord(f['corners_ra_deg'], f['corners_dec_deg']).transform_to(frame[0])
            ax.add_patch(Polygon(np.c_[fp.az.deg, fp.alt.deg], facecolor=color, edgecolor=color, alpha=.25))
            ax.scatter(aa.az.deg[i], aa.alt.deg[i], color=color, s=30)
            ax.annotate(str(i+1), (aa.az.deg[i], aa.alt.deg[i]), xytext=(3, 6), textcoords='offset points', fontsize=10)
        az_values = [*aa.az.deg, float(sun.az.deg[0])]
        if moon.alt.deg[0] >= -25:
            az_values.append(float(moon.az.deg[0]))
        else:
            ax.text(.98, .04, f'Moon altitude {moon.alt.deg[0]:.1f} deg (below panel)',
                    transform=ax.transAxes, ha='right', fontsize=10, color='.4')
        ax.set(xlim=(min(az_values)-10, max(az_values)+10), ylim=(-25, 50),
               xlabel='Azimuth (deg; E=90, S=180, W=270)', ylabel='Altitude (deg)',
               title=f'{side} snapshot: {stamp(snapshot):%m-%d %H:%M:%S} UTC+8')
        ax.legend(loc='upper right', fontsize=10)

        start, stop = block['reserved_window_unix']
        times = np.linspace(start, stop, 181)
        metrics = sky.metrics([f['ra_deg'] for f in fields], [f['dec_deg'] for f in fields], times)
        gaps_all = []
        for i, (f, color) in enumerate(zip(fields, colors)):
            rr = [r for r in rows if r['field_id'] == f['field_id']]
            mids = np.array([r['mid_unix'] for r in rr])
            m = sky.metrics([f['ra_deg']], [f['dec_deg']], mids)
            gaps = np.diff(mids)/60
            gaps_all.extend(gaps)
            field_stats.append({'side': side, 'field_id': f['field_id'], 'order': i+1,
                'ra_deg': f['ra_deg'], 'dec_deg': f['dec_deg'], 'visits': len(rr),
                'first_mid_local': local_time(mids[0]), 'last_mid_local': local_time(mids[-1]),
                'exposure_s': '/'.join(str(r['exposure_s']) for r in rr),
                'cadence_min_minutes': float(min(gaps)), 'cadence_median_minutes': float(np.median(gaps)),
                'cadence_max_minutes': float(max(gaps)), 'arc_span_minutes': float(np.ptp(mids)/60),
                'center_alt_at_mid_min_deg': float(m['alt'].min()),
                'footprint_alt_bound_at_mid_min_deg': float(m['alt'].min()-f['radius_deg']),
                'solar_elongation_mid_min_deg': float(m['elong'].min()),
                'solar_elongation_mid_max_deg': float(m['elong'].max()),
                'moon_sep_mid_min_deg': float(m['moon_sep'].min())})
            axes[0, column].plot([stamp(t) for t in times], metrics['alt'][i]-f['radius_deg'],
                                  color=color, lw=1.4, label=f['field_id'])
            axes[0, column].scatter([stamp(t) for t in mids], m['alt'][0]-f['radius_deg'], color=color, s=25)
            axes[1, column].barh(np.full(len(rr), i+1), [r['exposure_s']/86400 for r in rr],
                left=mdates.date2num([stamp(r['start_unix']) for r in rr]), height=.6, color=color)
            axes[3, column].scatter(i+1+np.linspace(-.15,.15,len(gaps)), gaps, color=color, s=30)
        axes[0, column].axhline(25, color='red', ls='--', lw=1)
        axes[0, column].set(ylabel='Footprint altitude\nlower bound (deg)', title=f'{side}: altitude; dots = exposure midpoints')
        axes[0, column].legend(fontsize=9, ncol=3, loc='upper right')
        axes[1, column].set(yticks=np.arange(1, len(fields)+1),
            yticklabels=[f['field_id'] for f in fields], title=f'{side}: continuous exposure timeline', ylabel='Field')
        for ax in axes[:2, column]:
            ax.set_xlim(stamp(start), stamp(stop))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=TZ))
            ax.set_xlabel(f'{stamp(start):%Y-%m-%d} Beijing time (UTC+8)')
            ax.axvline(stamp(block['window_unix'][0]), color='.3', ls=':', lw=1)
        axes[3, column].axhline(8, color='.5', ls='--', label='8 min preference')
        axes[3, column].axhline(5, color='red', ls=':', label='5 min floor')
        axes[3, column].set(xticks=np.arange(1, len(fields)+1), xticklabels=[f['field_id'] for f in fields],
            ylim=(4.5, max(9., max(gaps_all)+.5)), ylabel='Successive midpoints (min)',
            title=f'{side}: successive visits to the same field')
        axes[3, column].tick_params(axis='x', rotation=30)
        axes[3, column].legend(fontsize=10)
        ax = axes[2, column]
        ax.step([stamp(r['mid_unix']) for r in rows], [r['exposure_s'] for r in rows],
                where='mid', color='navy', lw=1.3, label='Fixed exposure')
        ax.scatter([stamp(r['mid_unix']) for r in rows], [r['exposure_s'] for r in rows], s=10, color='navy')
        ax.set(ylabel='Exposure (s)', yticks=[5,10,15,20,25,30], ylim=(0,35),
               title=f'{side}: exposure follows predicted twilight')
        ax.set_xlim(stamp(start), stamp(stop))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=TZ))
        ax.set_xlabel(f'{stamp(start):%Y-%m-%d} Beijing time (UTC+8)')
        sun_axis = ax.twinx()
        sun_axis.plot([stamp(t) for t in times], metrics['sun_alt'], color='goldenrod', lw=1, label='Sun altitude')
        sun_axis.axhline(-18, color='goldenrod', ls=':', lw=1)
        sun_axis.set_ylabel('Sun altitude (deg)', color='goldenrod')
        ax.legend(loc='upper left', fontsize=9)
        sun_axis.legend(loc='upper right', fontsize=9)
        spent = sum(r['exposure_s']+r['readout_s']+r['slew_s']+r['filter_change_s'] for r in rows)
        summaries.append({'side': side, 'fields': len(fields), 'frames': len(rows), 'filter': plan['config']['controller_filter_token'],
            'dispatch_local': local_time(block['window_unix'][0]),
            'first_shutter_local': local_time(rows[0]['start_unix']),
            'last_shutter_end_local': local_time(rows[-1]['end_unix']),
            'nominal_readout_complete_local': local_time(rows[-1]['release_unix']),
            'hard_end_local': local_time(stop),
            'window_seconds': stop-start, 'tail_unused_seconds': stop-rows[-1]['release_unix'],
            'visits_per_field': {f['field_id']: sum(r['field_id']==f['field_id'] for r in rows) for f in fields},
            'exposure_counts': {str(t): sum(r['exposure_s']==t for r in rows) for t in sorted({r['exposure_s'] for r in rows})},
            'exposures_seconds': sorted({r['exposure_s'] for r in rows}),
            'cadence_minutes_min_median_max': [float(min(gaps_all)), float(np.median(gaps_all)), float(max(gaps_all))],
            'integration_seconds': sum(r['exposure_s'] for r in rows), 'execution_seconds': spent,
            'slew_seconds': sum(r['slew_s'] for r in rows), 'readout_seconds': sum(r['readout_s'] for r in rows),
            'verified': True})
    for ax in [*maps.flat, *axes.flat]:
        ax.grid(alpha=.15)
    field_fig.suptitle(f"SHARP MP night {plan['date']} | fixed lists, {plan['config']['controller_filter_token']} | original grid", fontsize=20)
    stats_fig.suptitle(f"SHARP MP night {plan['date']} | planned sampling, not observed data", fontsize=20)
    save_pair(field_fig, output/'fields', 'Original MP footprint corners and same-epoch sky snapshots')
    save_pair(stats_fig, output/'statistics', 'Fixed-list altitude bounds, exposures and revisit cadence')
    plt.close('all')
    write_csv(output/'field_statistics.csv', field_stats)
    summary = {'date': plan['date'], 'source_plan_sha256': hashlib.sha256((source/'plan.json').read_bytes()).hexdigest(),
               'blocks': summaries, 'note': 'Nominal predictions only; midnight rolls into the next local date. No server operations.'}
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (output/'READ_ME.txt').write_text(
        'SHARP '+plan['date']+' night review\n'
        'Read START_TIMES.json. Start EV.txt and AM.txt separately at their Beijing times; AM is the next day.\n'
        'These are twilight-only eight-column lists, filter r, no WAIT or automatic cutoff.\n'
        'Exposure is a theoretical prior. Start at the fixed boundary; full-hour lists have no reserved 300-second delay budget.\n'
        'Fields/statistics figures are predictions, not observed data.\n')
    files=[p for p in sorted(output.iterdir()) if p.is_file() and p.suffix!='.zip' and p.name!='manifest.json']
    manifest={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    suffix = '_full_hour_review.zip' if plan['config'].get('fill_window') else '_tonight.zip'
    with zipfile.ZipFile(output/('SHARP_'+plan['date']+suffix), 'w', zipfile.ZIP_DEFLATED) as z:
        for file in files+[output/'manifest.json']:
            z.write(file, file.name)
    print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    report(Path(a.source), Path(a.output))


if __name__ == '__main__':
    main()
