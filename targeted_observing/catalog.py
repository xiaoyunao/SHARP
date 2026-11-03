"""Cached public catalogue acquisition and planning-only two-body preselection."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import requests
import astropy.units as u
from astropy.coordinates import get_body_barycentric, CartesianRepresentation, SkyCoord, BarycentricMeanEcliptic
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from .geometry import coord, in_field, isot

BASE=Path(__file__).parent


def download(url, path, params=None):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        r=requests.get(url,params=params,timeout=90); r.raise_for_status()
        path.write_bytes(r.content)
        path.with_suffix(path.suffix+'.source.json').write_text(json.dumps({
            'url':r.url,'retrieved_utc':datetime.now(timezone.utc).isoformat(),
            'sha256':hashlib.sha256(r.content).hexdigest()},indent=2)+'\n')
    return path


def approximate(catalog, date):
    """Two-body ranking only. Never used for final field retention or executable pointing."""
    t=Time(date)
    records=[dict(zip(catalog['fields'],r)) for r in catalog['data']]
    records=[r for r in records if r['H'] is not None and float(r['H'])<13.0]
    vals=lambda k:np.array([float(r[k]) for r in records])
    a,e=vals('a'),vals('e')
    ma=np.deg2rad(vals('ma'))+.01720209895/a**1.5*(t.tdb.jd-vals('epoch'))
    E=ma.copy()
    for _ in range(15): E-=(E-e*np.sin(E)-ma)/(1-e*np.cos(E))
    x=a*(np.cos(E)-e);y=a*np.sqrt(1-e*e)*np.sin(E)
    om,w,inc=[np.deg2rad(vals(k)) for k in ['om','w','i']]
    X=(np.cos(om)*np.cos(w)-np.sin(om)*np.sin(w)*np.cos(inc))*x+(-np.cos(om)*np.sin(w)-np.sin(om)*np.cos(w)*np.cos(inc))*y
    Y=(np.sin(om)*np.cos(w)+np.cos(om)*np.sin(w)*np.cos(inc))*x+(-np.sin(om)*np.sin(w)+np.cos(om)*np.cos(w)*np.cos(inc))*y
    Z=np.sin(w)*np.sin(inc)*x+np.cos(w)*np.sin(inc)*y
    ob=np.deg2rad(23.439291111)
    xyz=np.array([X,Y*np.cos(ob)-Z*np.sin(ob),Y*np.sin(ob)+Z*np.cos(ob)])
    earth=(get_body_barycentric('earth',t)-get_body_barycentric('sun',t)).xyz.to_value(u.au)
    geo=xyz-earth[:,None];r=np.linalg.norm(xyz,axis=0);d=np.linalg.norm(geo,axis=0)
    phase=np.arccos(np.clip(np.sum(xyz*geo,axis=0)/r/d,-1,1))
    phi1=np.exp(-3.33*np.tan(phase/2)**.63);phi2=np.exp(-1.87*np.tan(phase/2)**1.22)
    V=vals('H')+5*np.log10(r*d)-2.5*np.log10(.85*phi1+.15*phi2)
    c=SkyCoord(CartesianRepresentation(geo*u.au),frame='icrs')
    jup=(get_body_barycentric('jupiter',t)-get_body_barycentric('sun',t)).xyz.to_value(u.au)
    jlon=np.arctan2(jup[1]*np.cos(ob)+jup[2]*np.sin(ob),jup[0])
    diff=(np.arctan2(Y,X)-jlon+np.pi)%(2*np.pi)-np.pi
    for n,row in enumerate(records):
        row.update(ra_deg=float(c.spherical.lon.deg[n]),dec_deg=float(c.spherical.lat.deg[n]),
                   V_proxy=float(V[n]),cloud_proxy='L4' if diff[n]>0 else 'L5')
    return records


def horizons_track(pdes, start, stop, config, cache):
    key=f'{pdes}_{start}_{stop}_1h'
    path=Path(cache)/f'horizons_{key}.json'
    if path.exists(): return json.loads(path.read_text())
    s=config['site']
    h=Horizons(id=str(pdes),id_type='smallbody',location={'lon':s['lon_deg'],'lat':s['lat_deg'],'elevation':s['height_m']/1000},
               epochs={'start':start,'stop':stop,'step':'1h'})
    table=h.ephemerides(quantities='1,3,9,19,20,23,24',cache=False)
    keys=['datetime_jd','RA','DEC','RA_rate','DEC_rate','V','r','delta','elong','alpha']
    data={k:np.ma.filled(table[k],np.nan).tolist() for k in keys}
    data.update(pdes=str(pdes),targetname=str(table['targetname'][0]),
                source='JPL Horizons observer ephemeris; topocentric astrometric ICRF; UTC',
                query=h.uri,retrieved_utc=datetime.now(timezone.utc).isoformat(),site=config['site'])
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    return data


def interpolate(track, seconds):
    jd=np.asarray(track['datetime_jd']);q=Time(np.atleast_1d(seconds),format='unix').jd
    if np.min(q)<jd[0] or np.max(q)>jd[-1]: raise ValueError('Ephemeris outside coverage')
    ra=np.rad2deg(np.interp(q,jd,np.unwrap(np.deg2rad(track['RA']))))%360
    dec=np.interp(q,jd,track['DEC'])
    return ra,dec


def field_members(field, tracks, width=3, margin=.15):
    return [k for k,t in tracks.items() if np.all(in_field(t['RA'],t['DEC'],field['ra_deg'],field['dec_deg'],width,margin))]


def write_csv(path, rows):
    if not rows: return
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
