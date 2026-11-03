"""Geometric retention at planned exposure midpoints, separate from photometric usability."""
import numpy as np
from .catalog import interpolate
from .geometry import in_field,unix

def report(fields,tracks,visits_by_night,config):
    output={}
    for f in fields:
        for k in f['core_targets']:
            nights={};total=inside=0
            for date,rows in visits_by_night.items():
                ts=[unix(r['mid_utc']) for r in rows if r['field_id']==f['field_id']]
                if not ts:continue
                ra,dec=interpolate(tracks[k],ts)
                good=in_field(ra,dec,f['ra_deg'],f['dec_deg'],config['fov_deg'],config['edge_margin_deg'])
                nights[date]={'planned_points':len(ts),'geometric_points':int(good.sum()),'in_field':good.tolist()}
                total+=len(ts);inside+=int(good.sum())
            output[f['field_id']+':'+k]={'object_id':k,'field_id':f['field_id'],'planned_points':total,
                'geometric_points':inside,'geometric_retention':inside/total if total else 0.,
                'nights_with_geometric_points':sum(r['geometric_points']>0 for r in nights.values()),
                'nights':nights,'photometric_usable_fraction':None,'photometric_status':'requires_actual_images'}
    return output
