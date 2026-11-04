"""Independent hourly Horizons checks of the annual six-hour interpolation."""
import hashlib
import json
from datetime import date,timedelta
import numpy as np
from astropy.time import Time
from .catalog import BASE,interpolate,write_csv
from .annual_sources import fetch_one
from .annual_trojan import load_tracks
from .geometry import coord,setup_iers


def main():
    spec=json.loads((BASE/'annual_defaults.json').read_text());config=json.loads((BASE/'defaults.json').read_text())
    tracks=load_tracks();keys=sorted(tracks);rows=[]
    for month in range(12):
        k=keys[int((month+.5)*len(keys)/12)]
        start=date(2026+(month+10)//12,(month+10)%12+1,15)
        check=dict(spec,start_date=start.isoformat(),stop_date_exclusive=(start+timedelta(days=2)).isoformat(),lookahead_days=0,horizons_step='1h')
        path=fetch_one(k,check,config);reference=json.loads(__import__('pathlib').Path(path).read_text())
        seconds=Time(reference['datetime_jd'],format='jd').unix
        ra,dec=interpolate(tracks[k],seconds)
        separation=coord(ra,dec).separation(coord(reference['RA'],reference['DEC'])).arcsec
        dv=np.interp(reference['datetime_jd'],tracks[k]['datetime_jd'],tracks[k]['V'])-np.array(reference['V'])
        rows.append({'object_id':k,'reference_start':start.isoformat(),'hourly_samples':len(seconds),
                     'max_position_residual_arcsec':float(separation.max()),'max_V_residual_mag':float(abs(dv).max()),
                     'reference_sha256':hashlib.sha256(__import__('pathlib').Path(path).read_bytes()).hexdigest()})
        print(k,start,rows[-1]['max_position_residual_arcsec'],flush=True)
    root=BASE/'products/annual';write_csv(root/'ephemeris_interpolation_qa.csv',rows)
    iers=setup_iers(BASE/'cache'/spec['iers_file'])
    hashes={k:hashlib.sha256(__import__('pathlib').Path(path).read_bytes()).hexdigest()
            for k,path in json.loads((root/'ephemeris_index.json').read_text())['paths'].items()}
    result={'source':'JPL Horizons independent one-hour queries; same topocentric site',
            'sampled_objects':len(rows),'coverage':'12 epochs; not an exhaustive guarantee for all trajectories',
            'max_position_residual_arcsec':max(r['max_position_residual_arcsec'] for r in rows),
            'passed':all(r['max_position_residual_arcsec']<3 and r['max_V_residual_mag']<.02 for r in rows),
            'iers':iers,'ephemeris_sha256':hashes}
    (root/'input_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise SystemExit('Annual interpolation needs a finer step')


if __name__=='__main__':main()
