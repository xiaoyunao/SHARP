"""Compact report and portable ZIP for the offline MP experiment."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import healpy as hp
import astropy.units as u
from .catalog import BASE, write_csv
from .geometry import coord, unix
from .offline_twilight import load_grid


def audit_exports(root, plans):
    """Read actual controller files back, independently of the export writer."""
    from datetime import datetime, timedelta
    from astropy.coordinates import SkyCoord
    fields, provenance = load_grid()
    lookup = {f['field_id']: f for f in fields}
    frames = 0
    for plan in plans:
        directory = root/'nights'/plan['date']
        starts = json.loads((directory/'START_TIMES.json').read_text())
        assert starts['date'] == plan['date']
        assert plan['grid'] == provenance
        for block, dispatch in zip(plan['blocks'], starts['lists'], strict=True):
            assert block['verification']['passed']
            assert dispatch['side'] == block['side']
            start = datetime.fromisoformat(dispatch['start_local'])
            end = datetime.fromisoformat(dispatch['hard_end_local'])
            assert abs(start.timestamp()-block['window_unix'][0]) < 1
            assert abs(end.timestamp()-block['window_unix'][1]) < 1
            expected_day = datetime.fromisoformat(plan['date']).date()
            if block['side'] == 'AM':
                expected_day += timedelta(days=1)
            assert start.date() == expected_day
            path = directory/dispatch['list']
            assert hashlib.sha256(path.read_bytes()).hexdigest() == dispatch['sha256']
            lines = path.read_text().splitlines()
            assert len(lines) == len(block['rows']) == dispatch['exposures']
            for line, row in zip(lines, block['rows'], strict=True):
                parts = line.split()
                assert len(parts) == 8
                assert parts[0] == 'MP_'+row['field_id'][2:]
                assert parts[3:6] == ['0', plan['config']['controller_filter_token'], '1']
                assert int(parts[6]) == row['exposure_s'] and parts[7] == '0'
                field = lookup[row['field_id']]
                assert (row['ra_deg'], row['dec_deg']) == (field['ra_deg'], field['dec_deg'])
                exported = SkyCoord(parts[1], parts[2], unit=(u.hourangle, u.deg))
                assert exported.separation(coord(field['ra_deg'], field['dec_deg'])).arcsec < .02
                frames += 1
    return {'nights': len(plans), 'lists': 2*len(plans), 'frames': frames,
            'eight_columns_original_MP_R_exposure_and_hashes_passed': True}


def coverage(root, plans):
    fields,_=load_grid();lookup={f['field_id']:f for f in fields};pixels={}
    per_night=[]
    for p in plans:
        used={r['field_id'] for b in p['blocks'] for r in b['rows']}
        for fid in used-pixels.keys():
            f=lookup[fid]
            vertices=hp.ang2vec(f['corners_ra_deg'],f['corners_dec_deg'],lonlat=True)
            pixels[fid]=hp.query_polygon(512,vertices)
        per_night.append(np.unique(np.concatenate([pixels[f] for f in used])) if used else np.array([],int))
    all_pixels=np.concatenate(per_night);unique,n=np.unique(all_pixels,return_counts=True)
    unit=hp.nside2pixarea(512,degrees=True)
    result={'nside':512,'unique_nominal_area_deg2':len(unique)*unit,
            'nominal_area_ge3_nights_deg2':int(np.sum(n>=3))*unit,
            'method':'HEALPix pixel centres within original four-corner spherical polygons; nominal, not weather-adjusted'}
    (root/'coverage.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def audit_trojan_compatibility(root, plans, source):
    """Preserve the B block, only check that expanded shared night still contains it."""
    rows=[]
    for p in plans:
        path=source/(p['date']+'.json')
        old=json.loads(path.read_text())
        for b in old['blocks']:
            if b['mode']!='trojan':continue
            start=unix(b['not_before_utc']);end=unix(b['hard_deadline_utc'])
            free=p['unallocated_interval_unix']
            rows.append({'date':p['date'],'block_id':b['block_id'],'exposures':len(b['rows']),
                         'fits_new_shared_night':free[0]<=start<end<=free[1],
                         'source_file_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    if not all(r['fits_new_shared_night'] for r in rows):raise ValueError('Trojan overlap')
    write_csv(root/'existing_trojan_compatibility.csv',rows)
    return {'existing_B_blocks':len(rows),'all_fit_new_shared_night':True,
            'scope':'Existing B Trojan scheduling and science yields are preserved, not rerun by this twilight experiment'}


def figures(root, report, ledger, exposures):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .plot_style import apply_style,save_pair,MODES
    apply_style()
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    for side in ('EV','AM'):
        d=ledger[ledger.side==side];color=MODES[side]
        axes[0,0].plot(pd.to_datetime(d.date),d.fields,color=color,label=side,lw=1)
        valid=d.dropna(subset=['elongation_median_deg'])
        axes[0,1].scatter(pd.to_datetime(valid.date),valid.elongation_median_deg,s=12,color=color,label=side,alpha=.65)
        axes[1,0].hist(d.cadence_median_s.dropna()/60,bins=np.arange(4.5,12.6,.5),
                       histtype='step',lw=2,color=color,label=side)
    axes[0,0].set(ylabel='MP fields per window',title='Fixed MP lists')
    axes[0,1].set(ylabel='Median solar elongation (deg)',title='Near-Sun geometry')
    axes[1,0].set(xlabel='Median revisit interval (min)',ylabel='Windows',title='No in-list waiting')
    bins=np.arange(5,36,5);counts=[int((exposures.exposure_s==t).sum()) for t in bins[:-1]]
    axes[1,1].bar(bins[:-1],counts,width=3.5,color=MODES['EV'])
    axes[1,1].set(xlabel='Fixed exposure (s)',ylabel='Exposures',title='Theoretical exposure prior')
    for ax in axes.flat:ax.grid(alpha=.15)
    for ax in (axes[0,0],axes[0,1],axes[1,0]):ax.legend(frameon=False)
    for ax in axes[0]:ax.tick_params(axis='x',rotation=25)
    save_pair(fig,report/'annual_overview','MP offline nominal year, theoretical background; no discovery forecast')
    plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--report',required=True)
    p.add_argument('--trojan-source');a=p.parse_args();root=Path(a.input);report=Path(a.report);report.mkdir(parents=True,exist_ok=True)
    summary=json.loads((root/'summary.json').read_text())
    plans=[json.loads(p.read_text()) for p in sorted((root/'nights').glob('*/plan.json'))]
    if len(plans)!=summary['nights']:raise ValueError('Night count mismatch')
    summary['export_audit']=audit_exports(root,plans)
    summary['delivery_source_hashes']={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()
        for name in ('offline_delivery.py','annual_weather.py','data/weather_prior.json')}
    ledger=pd.read_csv(root/'time_ledger.csv');exposures=pd.read_csv(root/'exposures.csv')
    summary['coverage']=coverage(root,plans)
    if a.trojan_source:summary['trojan_compatibility']=audit_trojan_compatibility(root,plans,Path(a.trojan_source))
    summary['cadence_seconds_min_median_max']=[float(ledger.cadence_min_s.min()),float(ledger.cadence_median_s.median()),float(ledger.cadence_max_s.max())]
    summary['mean_optimized_cycle_slew_s']=float(ledger.optimized_cycle_slew_s.mean())
    summary['mean_RA_order_cycle_slew_s']=float(ledger.ra_cycle_slew_s.mean())
    summary['visit_count_windows']={str(k):int(v) for k,v in ledger.visits.value_counts().sort_index().items()}
    summary['solar_elongation_P10_P50_P90']=np.percentile(exposures.elongation_deg,[10,50,90]).tolist()
    ledger['month']=ledger.date.str[:7]
    monthly=ledger.groupby(['month','side']).agg(windows=('date','count'),
        active_windows=('fields',lambda x:int((x>0).sum())),field_windows=('fields','sum'),
        frames=('exposures','sum'),median_cadence_s=('cadence_median_s','median'),
        median_elongation_deg=('elongation_median_deg','median')).reset_index()
    monthly.to_csv(report/'monthly.csv',index=False)
    for name in ('time_ledger.csv','conditional_detections.csv'):
        (report/name).write_bytes((root/name).read_bytes())
    (report/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (root/'delivery_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (root/'READ_BEFORE_USE.txt').write_text(
        'SHARP MP offline fixed lists — local trial products\n'
        'Each date labels the evening; AM belongs to the next local calendar day.\n'
        'Read nights/DATE/START_TIMES.json, then start EV.txt and AM.txt separately at the stated Beijing times.\n'
        'The eight-column lists use filter R and contain NO WAIT, date, realtime feedback, or automatic cutoff commands.\n'
        'Late start and accumulated timing delays share a 300-second allowance. Stop/replan beyond it; TXT cannot enforce this.\n'
        'An empty list means no feasible group was found by this bounded search; do not substitute another date.\n'
        'Exposure times are an uncalibrated 5–30 second theoretical prior. Inspect trial images before adopting them routinely.\n'
        'These are twilight lists only. Existing Trojan blocks are checked separately, not appended here.\n'
        'summary.json is the simulation summary; delivery_summary.json additionally includes coverage, cadence and compatibility.\n', encoding='utf-8')
    figures(root,report,ledger,exposures)
    manifest=[]
    for file in sorted((root/'nights').glob('*/*')):
        if file.is_file():manifest.append({'file':str(file.relative_to(root)),
                            'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
    (root/'delivery_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    bundle=root/'SHARP_MP_offline_night_lists_365.zip'
    with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
        for entry in manifest:z.write(root/entry['file'],entry['file'])
        for name in ('summary.json','delivery_summary.json','READ_BEFORE_USE.txt','time_ledger.csv',
                     'coverage.json','delivery_manifest.json','conditional_detections.csv'):
            z.write(root/name,name)
        if a.trojan_source:
            z.write(root/'existing_trojan_compatibility.csv','existing_trojan_compatibility.csv')
        z.writestr('offline_defaults.json',json.dumps(plans[0]['config'],indent=2)+'\n')
    (report/'delivery.json').write_text(json.dumps({'bundle_name':bundle.name,
        'bundle_sha256':hashlib.sha256(bundle.read_bytes()).hexdigest(),'files':len(manifest),
        'raw_products':'Ignored products directory; reproduce from committed code/config/grid and pinned IERS.'},indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
