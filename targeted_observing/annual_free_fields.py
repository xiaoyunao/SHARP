"""Finite offset-centre search on the unchanged annual parent sample.

Pair midpoints and local footprint centroids supplement the frozen anchors.
This is a reproducible candidate search, not continuous global optimization.
"""
import numpy as np
from scipy.spatial import cKDTree
from .annual_trojan import snapshot
from .geometry import coord, in_field, unix
from .scheduler import field


def supplement(fields, tracks, start, stop, stars):
    ids,ra,dec,mag=snapshot(tracks,(unix(start)+unix(stop))/2)
    eligible=np.flatnonzero((mag<=19.2)&(dec>=-26)&np.isfinite(mag))
    xyz=coord(ra[eligible],dec[eligible]).cartesian.xyz.value.T
    centres=[]
    for a,b in sorted(cKDTree(xyz).query_pairs(2*np.sin(np.deg2rad(3.8)/2))):
        v=xyz[a]+xyz[b];v/=np.linalg.norm(v)
        centres.append((np.rad2deg(np.arctan2(v[1],v[0]))%360,np.rad2deg(np.arcsin(v[2]))))
    for f in fields:
        members=eligible[in_field(ra[eligible],dec[eligible],f['ra_deg'],f['dec_deg'],3,.15)]
        if len(members)>1:
            v=coord(ra[members],dec[members]).cartesian.xyz.value.mean(axis=1);v/=np.linalg.norm(v)
            centres.append((np.rad2deg(np.arctan2(v[1],v[0]))%360,np.rad2deg(np.arcsin(v[2]))))
    seen={(round(f['ra_deg'],4),round(f['dec_deg'],4)) for f in fields}
    star_coord=coord([s['ra_deg'] for s in stars],[s['dec_deg'] for s in stars]);mags=np.array([s['mag'] for s in stars])
    result=list(fields)
    for r,d in centres:
        key=round(r,4),round(d,4)
        if key in seen:continue
        seen.add(key);members=eligible[in_field(ra[eligible],dec[eligible],r,d,3,.15)]
        if not len(members):continue
        f=field(r,d,'AF',start);f['core_targets']=[ids[j] for j in members]
        f['campaign_members']=f['core_targets'];f['anchor_object']=None
        f['centre_method']='spherical pair midpoint or local centroid'
        separation=coord(r,d).separation(star_coord).deg
        f['very_bright_star_risk']=bool(np.any((mags<2)&(separation<4.122)))
        f['bright_stars_within_cap']=int(np.sum(separation<2.122));result.append(f)
    return result
