"""Assemble and independently audit the selected causal policy's nominal year.

This converts the numeric experiment back into the normal exposure contract.
It does not enable telescope export, install cron jobs, or contact a server.
"""
import copy
import csv
import json
from collections import defaultdict
from datetime import date, timedelta
import numpy as np
from .catalog import BASE, interpolate, write_csv
from .annual_sources import dates
from .annual_trojan import load_tracks
from .geometry import Sky, setup_iers, unix, isot, in_field
from .scheduler import visits, block_record, free_intervals
from .retention import report
from .validate import validate_plan


def assemble():
    spec=json.loads((BASE/'annual_defaults.json').read_text())
    root=BASE/'products/annual';source=root/'nominal_queue_4fields_dwell2_baseline_science25_v3'
    out=root/'selected_nominal';out.mkdir(exist_ok=True)
    tracks=load_tracks();setup_iers(BASE/'cache'/spec['iers_file'])
    c=json.loads((BASE/'defaults.json').read_text());c['trojan']['exposures_per_dwell']=2
    sky=Sky(c);nights=dates(spec['start_date'],spec['stop_date_exclusive'])
    by_day=defaultdict(list)
    with (source/'execution_realization0.csv').open() as handle:
        for row in csv.DictReader(handle):by_day[row['date']].append(row)
    with (source/'nightly_execution_ledger.csv').open() as handle:
        ledger={r['date']:r for r in csv.DictReader(handle)}
    plans=[];versions=defaultdict(list)
    for day in nights:
        menu=json.loads((root/'visibility'/f'{day}.json').read_text())
        p=json.loads((root/'twilight_5'/f'{day}.json').read_text())
        p['config']=copy.deepcopy(c);p['status']='LOCAL_ANNUAL_NOMINAL_NOT_FOR_TELESCOPE'
        p['source_policy']=source.name;p['observations_completed']=False
        p['nominal_feedback']='All-clear realization of causal completion policy; not an operational weather forecast.'
        numeric=by_day[day]
        if numeric:
            ids=list(dict.fromkeys(r['field_id'] for r in numeric))
            lookup={f['field_id']:f for f in menu['fields']};fs=[copy.deepcopy(lookup[k]) for k in ids]
            a=float(numeric[0]['block_start_unix']);b=float(numeric[0]['block_deadline_unix'])
            rows=visits(fs,a,b-c['trojan']['handover_s'],120,100,c,'TJ','trojan')
            assert len(rows)==len(numeric),(day,len(rows),len(numeric))
            np.testing.assert_allclose([unix(r['mid_utc']) for r in rows],
                                       [float(r['mid_unix']) for r in numeric],atol=.001,rtol=0)
            for f in fs:
                rr=[r for r in rows if r['field_id']==f['field_id']];ts=[unix(r['mid_utc']) for r in rr]
                for k in f['core_targets']:
                    ra,dec=interpolate(tracks[k],ts)
                    inside=in_field(ra,dec,f['ra_deg'],f['dec_deg'],3,.15)
                    for r,good in zip(rr,inside):r['target_states'][k]={'in_field':bool(good),'photometric_usable':None}
            p['blocks'].append(block_record('TJ',fs,rows,a,b,{}))
            version=fs[0]['version'];p['field_version_start']=version
            p['field_version_stop']=(date.fromisoformat(version)+timedelta(days=7)).isoformat()
            versions[version].append(p)
        else:
            p.setdefault('skips',[]).append('TJ: '+ledger[day]['reason'])
        p['blocks'].sort(key=lambda b:unix(b['not_before_utc']))
        occupied=[{'start_utc':b['not_before_utc'],'end_utc':b['hard_deadline_utc']} for b in p['blocks']]
        p['free_intervals']=[{'start_utc':isot(a),'end_utc':isot(b),'duration_s':b-a}
            for a,b in free_intervals((p['night']['evening'][0],p['night']['morning'][1]),occupied)]
        plans.append(p)
    audits=[];field_table=[];retention=[]
    campaigns={}
    for version,ps in versions.items():
        fields={f['field_id']:f for p in ps for b in p['blocks'] if b['block_id']=='TJ' for f in b['fields']}
        scheduled={p['date']:[r for b in p['blocks'] if b['block_id']=='TJ' for r in b['rows']] for p in ps}
        sampling=report(list(fields.values()),tracks,scheduled,c)
        for f in fields.values():
            secondary=[k for k in f['core_targets'] if sampling[f['field_id']+':'+k]['geometric_retention']<.8]
            f['secondary_targets']=secondary;f['core_targets']=[k for k in f['core_targets'] if k not in secondary]
            field_table.append(f)
        retention.extend({'version':version,'key':k,**{a:b for a,b in v.items() if a!='nights'}} for k,v in sampling.items())
        campaigns[version]={'start_date':version,'stop_date':ps[0]['field_version_stop'],
                            'fields':list(fields.values()),'tracks':tracks,'planned_visits':scheduled}
        for p in ps:
            for b in p['blocks']:
                if b['block_id']=='TJ':b['fields']=[fields[f['field_id']] for f in b['fields']]
    flat=[];time_rows=[]
    for i,p in enumerate(plans):
        audit=validate_plan(p,campaigns.get(p.get('field_version_start')))
        if not audit['passed']:raise RuntimeError((p['date'],audit))
        p['validation']=audit;audits.append({'date':p['date'],**audit})
        for b in p['blocks']:
            flat.extend({'date':p['date'],**r,'target_states':json.dumps(r['target_states'])} for r in b['rows'])
        allocated={kind:sum(unix(b['hard_deadline_utc'])-unix(b['not_before_utc']) for b in p['blocks'] if b['block_id'].startswith(kind)) for kind in ('EV','TJ','AM')}
        free=sum(x['duration_s'] for x in p['free_intervals']);total=p['night']['morning'][1]-p['night']['evening'][0]
        assert abs(sum(allocated.values())+free-total)<.01
        time_rows.append({'date':p['date'],**{k+'_allocated_s':v for k,v in allocated.items()},'free_s':free,'night_s':total,'trojan_zero_reason':ledger[p['date']]['reason']})
        (out/f"{p['date']}.json").write_text(json.dumps(p)+'\n')
        if (i+1)%30==0:print(f'Assembled and independently validated {i+1}/365 nights',flush=True)
    write_csv(out/'exposures.csv',flat);write_csv(out/'time_ledger.csv',time_rows);write_csv(out/'retention.csv',retention)
    (out/'fields.json').write_text(json.dumps(field_table)+'\n')
    result={'nights':len(plans),'twilight_windows':2*len(plans),'exposures':len(flat),'passed':all(a['passed'] for a in audits),
            'source':source.name,'independent_geometry_step_s':5,'records':audits}
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k!='records'},flush=True)


if __name__=='__main__':assemble()
