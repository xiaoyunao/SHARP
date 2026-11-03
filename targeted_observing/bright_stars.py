"""Hipparcos bright-end supplement plus Gaia; geometric risks, not measured halo radii."""
import json
import numpy as np
from astropy.time import Time
from .catalog import BASE,interpolate
from .geometry import coord,in_field,unix
from .review_data import table

def load(cache=BASE/'cache',epoch=2026.84):
    rows=[]
    for r in table(cache/'hip_bright_pm.tsv'):
        if not r['RAICRS'] or not r['DEICRS']:continue
        ra,dec=float(r['RAICRS']),float(r['DEICRS']);dt=epoch-1991.25
        ra+=float(r['pmRA'] or 0)*dt/3600000/np.cos(np.deg2rad(dec));dec+=float(r['pmDE'] or 0)*dt/3600000
        rows.append({'id':'HIP '+r['HIP'],'ra_deg':ra%360,'dec_deg':dec,'mag':float(r['Vmag']),'band':'V','source':'I/239/hip_main'})
    gp=cache/'gaia_bright_complete.json'
    if gp.exists():
        for k,ra,dec,g in json.loads(gp.read_text())['data']:
            rows.append({'id':'Gaia '+str(k),'ra_deg':ra,'dec_deg':dec,'mag':g,'band':'G','source':'Gaia DR3 cached region only'})
    return rows

def assess(campaign,visits,cache=BASE/'cache'):
    stars=load(cache);results=[]
    for f in campaign['fields']:
        # Include a two-degree exterior buffer for off-detector scattered light.
        nearby=[s for s in stars if coord(s['ra_deg'],s['dec_deg']).separation(coord(f['ra_deg'],f['dec_deg'])).deg<4.122]
        for k in f['core_targets']:
            for date,rows in visits.items():
                rr=[r for r in rows if r['field_id']==f['field_id']]
                if not rr:continue
                ra,dec=interpolate(campaign['tracks'][k],[unix(r['mid_utc']) for r in rr])
                for s in nearby:
                    sep=coord(ra,dec).separation(coord(s['ra_deg'],s['dec_deg'])).arcmin
                    if sep.min()>60 and s['mag']>=2:continue
                    results.append({'night':date,'field_id':f['field_id'],'object_id':k,'star':s,
                        'star_in_footprint':bool(in_field(s['ra_deg'],s['dec_deg'],f['ra_deg'],f['dec_deg'],3,0)),
                        'separation_arcmin':sep.tolist(),'min_arcmin':float(sep.min()),'max_arcmin':float(sep.max()),
                        'status':'geometric_proximity_requires_images; not a calibrated exclusion radius',
                        'photometric_usable_fraction':None})
    return {'risks':results,'hip_rows':sum(s['source']=='I/239/hip_main' for s in stars),
        'limits':'Hipparcos V<6 global bright supplement; Gaia cached region only, V/G not interchangeable; no clean-field certification'}
