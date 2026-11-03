"""Gaia DR3 bright-star screening; not a validated optical halo or completeness model."""
import json
from pathlib import Path
import numpy as np
from .catalog import BASE,download,write_csv
from .geometry import coord,in_field


def main():
    camp=json.loads((BASE/'data/campaign_20261103.json').read_text())
    query="SELECT source_id,ra,dec,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',89,10,4)) AND phot_g_mean_mag<8"
    params={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'json','QUERY':query,'MAXREC':'10000'}
    p=download('https://gea.esac.esa.int/tap-server/tap/sync',BASE/'cache/gaia_bright_complete.json',params)
    data=json.loads(p.read_text())['data']
    if len(data)>=10000:raise ValueError('Possible truncated catalogue; do not report complete bright-star check')
    stars=coord([r[1] for r in data],[r[2] for r in data]);rows=[]
    for f in camp['fields']:
        inside=[r for r in data if in_field(r[1],r[2],f['ra_deg'],f['dec_deg'],3,0)]
        for k in f['core_targets']:
            t=camp['tracks'][k];seps=coord(t['RA'],t['DEC'])[:,None].separation(stars[None,:]).deg
            idx=np.unravel_index(seps.argmin(),seps.shape)
            rows.append({'pdes':k,'field_id':f['field_id'],'galactic_lat_deg':float(coord(f['ra_deg'],f['dec_deg']).galactic.b.deg),
                'field_stars_G_lt8':len(inside),'field_brightest_G':min(r[3] for r in inside) if inside else None,
                'nearest_track_star_arcmin':float(seps.min()*60),'nearest_star_G':data[idx[1]][3],
                'status':'requires_on_sky_halo_and_blending_QA; not a clean-field certification'})
    write_csv(BASE/'products/campaign_20261103/contamination.csv',rows)
    print(json.dumps({'query':query,'gaia_rows':len(data),'targets':rows},indent=2))

if __name__=='__main__':main()
