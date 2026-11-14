"""Read-only server-side catalog extraction; CONFIG is supplied by the caller.

Execute over SSH stdin using the server default Python. No remote files written.
Returns target detections and identical sky reference stars across all nights.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree


def xyz(ra, dec):
    a, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack((np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)))


def main(config):
    columns = ['objID', 'RA_Win', 'DEC_Win', 'X_Win', 'Y_Win', 'Flag',
               'Flag_ISO_Num', 'FWHM', 'Class_ANN', 'Mag_PSF', 'MagErr_PSF',
               'Mag_Kron', 'MagErr_Kron', 'Mag_Aper3', 'MagErr_Aper3',
               'Mag_Aper4', 'MagErr_Aper4', 'Mag_Aper5', 'MagErr_Aper5']
    references, measurements, inventory = {}, [], []
    targets = config['tracks']
    for night in config['nights']:
        root = Path('/processed1')/night
        matched = fits.getdata(root/'L4'/f'{night}_matched_asteroids.fits', 1)
        files = sorted((root/'L2').glob('*TJ*cat.fits*'))
        for f in files:
            d, h = fits.getdata(f, 1, header=True)
            start = datetime.fromisoformat(h['EXPSTA']).replace(tzinfo=timezone.utc).timestamp()
            end = datetime.fromisoformat(h['EXPEND']).replace(tzinfo=timezone.utc).timestamp()
            jd = (start+end)/2/86400+2440587.5
            tree = cKDTree(xyz(d['RA_Win'], d['DEC_Win']))
            inv = dict(night=night, file=f.name, start=start, end=end, jd=jd,
                       filter=h.get('FILTER'), n_sources=len(d), targets=[])
            for number in config['numbers']:
                m = matched[(matched['number']==int(number)) & (matched['source_file']==f.name)]
                if len(m)>1:
                    raise ValueError('Duplicate target identity in '+f.name)
                if len(m):
                    indices = np.flatnonzero(d['objID']==m['objID'][0])
                    assert len(indices)==1
                    idx = int(indices[0])
                    sep = 3600*np.hypot((d['RA_Win'][idx]-m['ra'][0])*np.cos(np.deg2rad(m['dec'][0])), d['DEC_Win'][idx]-m['dec'][0])
                    identity = 'production known match'
                elif night=='20261111' and number in ('17365','4722'):
                    track = targets[number]
                    assert track['datetime_jd'][0] <= jd <= track['datetime_jd'][-1]
                    ra, dec = [np.interp(jd,track['datetime_jd'],track[c]) for c in ['RA','DEC']]
                    dist, ii = tree.query(xyz([ra],[dec]), k=2)
                    sep = float(dist[0,0]*206264.806)
                    if sep>2 or dist[0,1]*206264.806<3:
                        continue
                    idx = int(ii[0,0]); identity = 'cached topocentric Horizons, 2 arcsec unique match'
                else:
                    continue
                inv['targets'].append(number)
                row = {c:float(d[c][idx]) for c in columns}
                if number not in references:
                    ra, dec = row['RA_Win'],row['DEC_Win']
                    distance = np.hypot((d['RA_Win']-ra)*np.cos(np.deg2rad(dec)),d['DEC_Win']-dec)
                    valid = ((distance<.45)&(d['Flag']==0)&(d['Flag_ISO_Num']==0)&
                             (d['Class_ANN']>.8)&(d['Mag_PSF']>14)&(d['Mag_PSF']<18)&
                             (d['MagErr_PSF']<.025)&(d['MagErr_PSF']>0))
                    inds = np.flatnonzero(valid)
                    nearest = tree.query(xyz(d['RA_Win'][inds],d['DEC_Win'][inds]),k=2)[0][:,1]*206264.806
                    inds = inds[nearest>10]
                    # Remove all stars close to the target's entire three-night path.
                    if number in targets:
                        tr=targets[number]; tt=np.asarray(tr['datetime_jd'])
                        take=(tt>=2461355.5)&(tt<=2461358.5)
                        tracktree=cKDTree(xyz(np.asarray(tr['RA'])[take],np.asarray(tr['DEC'])[take]))
                        keep=tracktree.query(xyz(d['RA_Win'][inds],d['DEC_Win'][inds]))[0]*206264.806>20
                        inds=inds[keep]
                    else:
                        inds=inds[distance[inds]>.08]
                    inds=inds[np.argsort(d['MagErr_PSF'][inds])[:150]]
                    references[number]=[dict(ra=float(d['RA_Win'][i]),dec=float(d['DEC_Win'][i])) for i in inds]
                refs=references[number]
                dd, ii=tree.query(xyz([r['ra'] for r in refs],[r['dec'] for r in refs]),k=1)
                rows=[]
                for dist,i in zip(dd,ii):
                    rows.append({c:float(d[c][i]) for c in columns} if dist*206264.806<1 else None)
                measurements.append(dict(night=night,file=f.name,jd=jd,exposure_s=end-start,
                    number=number,identity=identity,sep_arcsec=float(sep),target=row,refs=rows))
            inventory.append(inv)
        print(f'{night}: {len(files)} catalogs read',file=sys.stderr,flush=True)
    print(json.dumps(dict(config=config,references=references,measurements=measurements,inventory=inventory),allow_nan=True))


if __name__ == '__main__':
    main(CONFIG)
