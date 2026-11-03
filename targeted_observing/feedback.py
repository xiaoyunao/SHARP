"""Summarize actual exposure logs without confusing planned visits with completed visits."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import numpy as np
from .geometry import unix

COLUMNS=['date','block_id','field_id','visit','actual_start_utc','actual_end_utc','filter','valid','background','fwhm_arcsec','zeropoint','limiting_mag','stable_star_rms_mag','zero_point_offset_mag','readout_s','slew_s','known_moving_expected','known_moving_recovered','measured_r_mag','notes']


def summarize(rows):
    groups=defaultdict(list);seen=set()
    for r in rows:
        key=tuple(r[k] for k in ['date','block_id','field_id','visit'])
        if key in seen:raise ValueError('Duplicate actual exposure key '+str(key))
        seen.add(key)
        if r['valid'] not in {'0','1'}:raise ValueError('valid must be explicit 0/1')
        if r['valid']=='1':
            if unix(r['actual_end_utc'])<=unix(r['actual_start_utc']):raise ValueError('Invalid actual exposure times')
            groups[key[:3]].append(r)
        else:groups[key[:3]]
    output=[]
    for key,rr in groups.items():
        mids=sorted((unix(r['actual_start_utc'])+unix(r['actual_end_utc']))/2 for r in rr)
        output.append({'date':key[0],'block_id':key[1],'field_id':key[2],'valid_visits':len(rr),
            'cadence_s':float(np.median(np.diff(mids))) if len(mids)>1 else None,
            'arc_span_s':mids[-1]-mids[0] if mids else 0,
            'stable_star_rms_median_mag':float(np.median([float(r['stable_star_rms_mag']) for r in rr if r['stable_star_rms_mag']])) if any(r['stable_star_rms_mag'] for r in rr) else None})
    return {'fields':output,'note':'No asteroid RMS substituted for measurement error; empty rows are not completed observations'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--template');p.add_argument('--input');p.add_argument('--output');a=p.parse_args()
    if a.template:
        with Path(a.template).open('w',newline='') as f:csv.writer(f,lineterminator='\n').writerow(COLUMNS)
    elif a.input and a.output:
        rows=list(csv.DictReader(Path(a.input).open()));Path(a.output).write_text(json.dumps(summarize(rows),indent=2)+'\n')
    else:p.error('Use --template PATH or --input CSV --output JSON')

if __name__=='__main__':main()
