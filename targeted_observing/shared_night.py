"""Local v2 grid integration probe: verified MP lists and TJN candidate membership."""
import argparse
import json
from pathlib import Path

import numpy as np

from .annual_trojan import load_tracks
from .catalog import write_csv
from .offline_twilight import load_config, make_plan, write_plan
from .shared_fields import SharedFields, cloud_labels


def build(day, output):
    output = Path(output)
    catalog = SharedFields()
    config = load_config()
    config['version'] = 'shared-grid-integration-20261118'
    config['assumptions'] = [s for s in config['assumptions'] if not s.startswith('Original MP coordinates')]
    config['assumptions'].append('Shared EQ-Xinglong18-20261116 centres; MP centre ecliptic latitude remains within 20 degrees. Representative TAN enclosing cap.')
    plan = make_plan(day, config, grid=catalog.twilight_grid())
    plan['status'] = 'LOCAL_V2_GRID_PROBE_NOT_FULL_NIGHT_SCHEDULE'
    write_plan(plan, output / 'twilight')
    # Broad per-night screening; final TJN allocation is deliberately not the old
    # fixed 21:52/01:40 November test timetable.
    start = max(plan['dark_crossings'][0], plan['unallocated_interval_unix'][0])
    end = min(plan['dark_crossings'][1], plan['unallocated_interval_unix'][1])
    if start >= end:
        raise ValueError('No dark interval remains for candidate membership')
    times = np.linspace(start, end, int(np.ceil((end-start)/600))+1)
    rows = catalog.memberships(load_tracks(), times, cloud_labels())
    # Explicit header even for an empty candidate night.
    (output/'tjn_memberships.csv').write_text('object_id,cloud,field_id,execution_name,ra_deg,dec_deg,retained_fraction,V_min,V_max\n')
    write_csv(output/'tjn_memberships.csv', rows)
    result = dict(date=day, grid=catalog.manifest, mp_blocks=[dict(
        side=b['side'], fields=len(b['fields']), frames=len(b['rows']),
        passed=b['verification']['passed']) for b in plan['blocks']],
        tjn_memberships=len(rows), tjn_unique_objects=len({r['object_id'] for r in rows}),
        tjn_unique_fields=len({r['field_id'] for r in rows}),
        membership_sample_start_unix=start, membership_sample_end_unix=end,
        status='GRID_AND_MEMBERSHIP_PROBE_ONLY',
        outstanding='TJN visibility/priority/allocation, YSO reservation and full-night replay must be integrated before a v2 full-night TXT is delivered.')
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--date', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.date, args.output), indent=2))
