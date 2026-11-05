"""Summarize B/C cadence, effective nights, science gates and role yield."""
import json
import gzip
from datetime import date

import numpy as np
import pandas as pd

from .annual_weather import lightcurve_quality
from .annual_tasks import load_tasks
from .continuation_report import OUT, ROOT, sha

FOCUS = [17365, 22180, 4902, 5511, 5012, 7543, 5284, 617]


def cadence():
    with gzip.open(OUT/'B_plot_data/exposures.csv.gz', 'rt') as file:
        nominal = pd.read_csv(file)
    nominal = nominal[nominal['mode'] == 'trojan'].copy()
    nominal['mid_unix'] = pd.to_datetime(nominal.mid_utc, format='ISO8601').astype('int64') / 1e9
    scenarios = [('B nominal', nominal)]
    for strategy in 'BC':
        path = ROOT/f'continuation_{strategy}_baseline/execution_realization0.csv'
        scenarios.append((f'{strategy} baseline realization 0', pd.read_csv(path)))
    visits, gaps = [], []
    for scenario, data in scenarios:
        for (night, field), group in data.groupby(['date', 'field_id']):
            mid = np.sort(group.mid_unix.to_numpy())
            if len(mid) % 2:
                raise ValueError(f'Incomplete two-frame dwell: {scenario} {night} {field}')
            visits.append(dict(scope=scenario, date=night, field_id=field, exposures=len(mid),
                               midpoint_span_s=np.ptp(mid)))
            for kind, values in [('within_dwell', mid[1::2]-mid[::2]),
                                 ('same_field_first_to_first', np.diff(mid[::2])),
                                 ('same_field_last_to_next_first', mid[2::2]-mid[1:-1:2])]:
                gaps.extend(dict(scope=scenario, date=night, field_id=field, kind=kind, seconds=value)
                            for value in values)
    pd.DataFrame(visits).to_csv(OUT/'cadence_field_nights.csv', index=False)
    series = pd.DataFrame(gaps)
    series.to_csv(OUT/'cadence_intervals.csv.gz', index=False,
                  compression={'method':'gzip', 'mtime':0})
    summary = series.groupby(['scope', 'kind']).seconds
    table = summary.quantile([.1, .5, .9]).unstack().rename(columns={.1:'p10_s', .5:'p50_s', .9:'p90_s'})
    table['intervals'] = summary.size()
    table.to_csv(OUT/'cadence_quantiles.csv')


def effective_nights(state, quality, mask):
    times, nights, errors, offsets, valid = state
    eligible = mask & valid & np.isfinite(times) & (errors <= .05) & (np.abs(offsets) <= quality['max_zeropoint_mag'])
    return np.array([night for night in np.unique(nights[eligible])
        if (eligible & (nights == night)).sum() >= quality['points_per_night']
        and np.ptp(times[eligible & (nights == night)]) >= quality['night_span_s']], dtype=int)


def object_statistics():
    cross_night, campaign_rows, focus, precision, ensemble_rows = [], [], [], [], []
    manifests = {}
    definitions = load_tasks()
    for strategy in 'BC':
        root = ROOT/f'continuation_{strategy}_baseline'
        for name in ['run_config.json', 'samples_realization0.json', 'execution_realization0.csv',
                     'science_task_status.csv', 'summary.json', 'object_yield.csv']:
            manifests[f'{root.name}/{name}'] = sha(root/name)
        spec = json.loads((root/'run_config.json').read_text())['spec']
        quality = spec['quality']
        length = spec['campaign_days']
        start = date.fromisoformat(spec['start_date'])
        certificates = pd.read_csv(root/'object_yield.csv')
        counts = certificates[certificates.quality_lightcurve]
        for column in ['good_points', 'good_nights', 'baseline_days']:
            p10, p50, p90 = counts[column].quantile([.1, .5, .9])
            ensemble_rows.append(dict(strategy=strategy, realizations=spec['weather_realizations'],
                qualified_object_realizations=len(counts), metric=column, p10=p10, p50=p50, p90=p90))
        history = json.loads((root/'samples_realization0.json').read_text())
        for object_id, record in history.items():
            state = [np.asarray(record[key]) for key in ['times','night_indices','errors','zeropoints','valid']]
            nights = state[1]
            prefix = np.arange(len(nights)) < record.get('completion_index', len(nights))
            good = effective_nights(state, quality, np.ones(len(nights), dtype=bool))
            certificate = effective_nights(state, quality, prefix)
            passed = record['quality']['quality_lightcurve']
            if len(certificate) != record['quality']['good_nights']:
                raise ValueError(f'Effective-night count disagrees with saved certificate: {object_id}')
            # The frozen year starts on a campaign boundary. Campaign indices
            # reproduce visibility_free's consecutive 35-night menus.
            single_campaign_pass = False
            for campaign in np.unique(nights // length):
                mask = nights // length == campaign
                subset = good[good//length == campaign]
                passes = lightcurve_quality(*(a[mask & prefix] for a in state), quality)['quality_lightcurve']
                single_campaign_pass |= passes
                campaign_rows.append(dict(strategy=strategy, realization=0, object_id=object_id,
                    campaign_start=date.fromordinal(start.toordinal()+int(campaign)*length).isoformat(),
                    effective_nights=len(subset), baseline_days=np.ptp(subset) if len(subset) else 0,
                    quality_pass_inside_certificate_prefix=passes))
            gaps = np.diff(certificate)
            breaks = np.flatnonzero(gaps >= length) + 1
            linked_long_gap = False
            if passed and len(breaks):
                starts = certificate[breaks]
                segments = np.searchsorted(starts, nights)
                linked_long_gap = not any(lightcurve_quality(*(a[(segments==i)&prefix] for a in state),
                    quality)['quality_lightcurve'] for i in range(len(starts)+1))
            cross_night.append(dict(strategy=strategy, realization=0, object_id=object_id,
                quality_lightcurve=passed, effective_nights=len(good), annual_baseline_days=np.ptp(good) if len(good) else 0,
                longest_gap_days=np.diff(good).max(initial=0), certificate_effective_nights=len(certificate),
                certificate_baseline_days=np.ptp(certificate) if len(certificate) else 0,
                certificate_longest_gap_days=gaps.max(initial=0),
                requires_cross_campaign_certificate=passed and not single_campaign_pass,
                requires_linking_across_gap_ge_35d=linked_long_gap))
        tasks = pd.read_csv(root/'science_task_status.csv')
        for object_id, group in tasks[tasks.object_id.isin(FOCUS)].groupby('object_id'):
            fail = dict(points=group.good_task_points < 100, nights=group.task_nights < 14,
                        baseline=group.task_baseline_days < group.required_baseline_days,
                        phase=group.min_phase_coverage < .8)
            if len(group) != spec['weather_realizations'] or not (
                    group.sampling_ready == ~np.logical_or.reduce(list(fail.values()))).all():
                raise ValueError(f'Saved task gates disagree: {strategy} {object_id}')
            record = history.get(str(object_id), {})
            times = np.asarray(record.get('times', []))
            calibrated = np.asarray(record.get('valid', []), dtype=bool) & (
                np.abs(record.get('zeropoints', [])) <= quality['max_zeropoint_mag'])
            limit = group.required_single_point_error_mag.iloc[0]
            usable = calibrated & (np.asarray(record.get('errors', [])) <= limit)
            for period in definitions[str(object_id)]['periods_h']:
                phase_boxes = np.floor((times[usable] / (period * 3600) % 1) * 20)
                precision.append(dict(strategy=strategy, realization=0, object_id=object_id,
                    period_hours=period, phase_coverage=len(np.unique(phase_boxes))/20,
                    valid_within_zeropoint_limit=calibrated.sum(),
                    good_after_task_error_limit=usable.sum(), error_limit_mag=limit))
            focus.append(dict(strategy=strategy, object_id=object_id, name=group.name.iloc[0],
                realizations=len(group), error_limit_mag=group.required_single_point_error_mag.iloc[0],
                min_required_days=group.required_baseline_days.iloc[0],
                median_good_points=group.good_task_points.median(), median_effective_nights=group.task_nights.median(),
                median_baseline_days=group.task_baseline_days.median(), median_min_phase=group.min_phase_coverage.median(),
                photometric_complete=group.photometric_complete.mean(), sampling_ready=group.sampling_ready.mean(),
                period_solution_complete=None, science_case_complete=None,
                **{f'fails_{key}_fraction':value.mean() for key,value in fail.items()}))
    pd.DataFrame(cross_night).to_csv(OUT/'BC_object_timing_realization0.csv', index=False)
    pd.DataFrame(campaign_rows).to_csv(OUT/'BC_campaign_timing_realization0.csv', index=False)
    pd.DataFrame(ensemble_rows).to_csv(OUT/'BC_ensemble_certificates.csv', index=False)
    pd.DataFrame(precision).to_csv(OUT/'BC_focus_precision_realization0.csv', index=False)
    pd.DataFrame(focus).to_csv(OUT/'BC_focus_science_gates.csv', index=False)
    path = OUT/'object_selection_function.csv'
    selection = pd.read_csv(path)
    role_yield = selection[selection.strategy.isin(['B','C'])].groupby(['strategy','role']).quality_probability.sum()
    role_yield.rename('expected_quality_objects').to_csv(OUT/'BC_role_yield.csv')
    manifests[path.name] = sha(path)
    manifests['B_plot_data/manifest.json'] = sha(OUT/'B_plot_data/manifest.json')
    (OUT/'BC_statistics_inputs.json').write_text(json.dumps(manifests, indent=2)+'\n')


if __name__ == '__main__':
    cadence()
    object_statistics()
