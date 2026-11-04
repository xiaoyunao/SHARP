"""Conservative exposure-wide geometry after simulated timing shifts."""
import numpy as np
from astropy.coordinates import GeocentricMeanEcliptic
from astropy.time import Time
from .geometry import unix,coord


def prepare(block,sky):
    if all(k in block for k in ['timing_grid','timing_alt','timing_moon_loss']):return
    c=sky.config;a=unix(block['not_before_utc']);b=unix(block['hard_deadline_utc'])
    grid=np.linspace(a,b,int(np.ceil((b-a)/30))+1)
    fields=block['fields'];m=sky.metrics([f['ra_deg'] for f in fields],[f['dec_deg'] for f in fields],grid)
    guard=.0043*15+.01;radius=2.122
    valid=(m['alt']>=max(c['center_alt_deg'],c['min_alt_deg']+radius)+guard)
    moonlim=np.where(m['moon_fraction']>=c['bright_moon_fraction']-.001,c['bright_moon_sep_deg'],c['moon_sep_deg'])
    valid&=(m['moon_sep']>=moonlim[None,:]+radius+guard)|(m['moon_alt'][None,:]<-guard)
    valid&=m['sun_alt'][None,:]<=c['twilight']['sun_alt_deg']+.002
    hi=c['twilight']['extended_elongation_deg'][1] if block['extended_search'] else c['twilight']['elongation_deg'][1]
    valid&=(m['elong']>=c['twilight']['elongation_deg'][0]+.001)&(m['elong']<=hi-.001)
    ecliptic=GeocentricMeanEcliptic(equinox=Time('J2000'),obstime=Time(grid,format='unix'))
    beta=coord([f['ra_deg'] for f in fields],[f['dec_deg'] for f in fields])[:,None].transform_to(ecliptic).lat.deg
    valid&=abs(beta)<=c['twilight']['beta_max_deg']-.00001
    block['timing_grid']=grid.tolist()
    block['timing_valid']={f['field_id']:valid[i].tolist() for i,f in enumerate(fields)}
    block['timing_alt']={f['field_id']:m['alt'][i].tolist() for i,f in enumerate(fields)}
    illuminated=np.where(m['moon_alt']>0,m['moon_fraction'],0)
    loss=1.25*np.log10(1+8*illuminated[None,:]*np.exp(-m['moon_sep']/40))
    block['timing_moon_loss']={f['field_id']:loss[i].tolist() for i,f in enumerate(fields)}
    block.pop('_timing_matrix',None)


def mask(block,field_ids,midpoints,exposure_s=30):
    """Accept only when all neighboring grid intervals covering exposure are valid."""
    if '_timing_matrix' not in block:
        block['_timing_index']={fid:i for i,fid in enumerate(block['timing_valid'])}
        block['_timing_matrix']=np.array(list(block['timing_valid'].values()),dtype=bool)
    grid=np.asarray(block['timing_grid']);times=np.asarray(midpoints);out=np.ones(times.shape,bool)
    lo=np.searchsorted(grid,times-exposure_s/2,side='right')-1
    hi=np.searchsorted(grid,times+exposure_s/2,side='left')
    out&=(lo>=0)&(hi<len(grid))
    # With a 30s exposure and <=30s grid spacing, at most three grid nodes occur.
    # Require every enclosing node rather than interpolating boolean validity.
    row_codes=np.array([block['_timing_index'][fid] for fid in field_ids])
    for offset in range(4):
        ix=np.minimum(lo+offset,hi)
        out&=block['_timing_matrix'][row_codes,np.clip(ix,0,len(grid)-1)]
    return out


def row_times(block):
    if '_row_times' not in block:block['_row_times']=np.array([unix(r['mid_utc']) for r in block['rows']])
    return block['_row_times']
