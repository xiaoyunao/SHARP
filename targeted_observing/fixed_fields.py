"""Long-lived equatorial TJN cells, independent of nightly target membership.

No production imports. Registry IDs are append-only, including retired fields.
The nominal design is distinct from measured detector/template certification.
"""
from dataclasses import dataclass, asdict
import copy
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


REGISTRY_PATH = Path(__file__).resolve().parent/'data/tjn_fixed_20261111/registry.json'


@dataclass(frozen=True)
class Grid:
    sx: float = 2.5
    sy: float = 2.5
    version: str = 'EQ25-v1'
    # Inscribed square of the smaller detector dimension at 1.15 arcsec/pixel.
    width: float = 9216 * 1.15 / 3600
    margin: float = .15
    guard: float = .01

    def row(self, j):
        dec = np.asarray(j) * self.sy
        if np.any(abs(dec) >= 90):
            raise ValueError('Polar rows are outside this regional grid')
        n = np.ceil(360 * np.cos(np.deg2rad(dec)) / self.sx).astype(int)
        return dec, n, 360 / n

    def assign(self, ra, dec):
        ra, dec = np.broadcast_arrays(ra, dec)
        j = np.floor(dec / self.sy + .5).astype(int)
        _, n, step = self.row(j)
        k = np.floor((ra % 360) / step + .5).astype(int) % n
        return j, k

    def key(self, j, k):
        return f'{self.version}_R{int(j):+04d}_C{int(k):03d}'

    def centre(self, j, k):
        dec, _, step = self.row(j)
        return np.asarray(k) * step, dec


def tan_xy(ra, dec, ra0, dec0):
    da = np.deg2rad((np.asarray(ra) - ra0 + 180) % 360 - 180)
    d, d0 = np.deg2rad(dec), np.deg2rad(dec0)
    den = np.sin(d0)*np.sin(d) + np.cos(d0)*np.cos(d)*np.cos(da)
    x = np.cos(d)*np.sin(da)/den
    y = (np.cos(d0)*np.sin(d)-np.sin(d0)*np.cos(d)*np.cos(da))/den
    return np.where(den > 0, np.rad2deg(x), np.inf), np.where(den > 0, np.rad2deg(y), np.inf)


def from_tan(ra0, dec0, x, y):
    a, d = np.deg2rad(ra0), np.deg2rad(dec0)
    x, y = np.deg2rad(x), np.deg2rad(y)
    vx = np.cos(d)*np.cos(a)-x*np.sin(a)-y*np.sin(d)*np.cos(a)
    vy = np.cos(d)*np.sin(a)+x*np.cos(a)-y*np.sin(d)*np.sin(a)
    vz = np.sin(d)+y*np.cos(d)
    return np.rad2deg(np.arctan2(vy,vx)) % 360, np.rad2deg(np.arctan2(vz,np.hypot(vx,vy)))


def edge_distance(ra, dec, ra0, dec0, grid=Grid()):
    x, y = tan_xy(ra,dec,ra0,dec0)
    return grid.width/2 - np.maximum(abs(x),abs(y))


def boundary(ra, dec, grid=Grid(), margin=0, samples=65):
    q = np.linspace(-1,1,samples); h = grid.width/2-margin
    x = h*np.r_[q,np.ones(samples),-q,-np.ones(samples)]
    y = h*np.r_[-np.ones(samples),q,np.ones(samples),-q]
    return np.column_stack(from_tan(ra,dec,x,y))


def cell_boundary(j, grid=Grid(), samples=2049):
    d, _, step = grid.row(j);q=np.linspace(-.5,.5,samples)
    return np.r_[q*step,np.full(samples,step/2),-q*step,np.full(samples,-step/2)], np.r_[np.full(samples,d-grid.sy/2),d+q*grid.sy,np.full(samples,d+grid.sy/2),d-q*grid.sy]


def unit_vectors(ra,dec):
    a,d=np.deg2rad(ra),np.deg2rad(dec)
    return np.column_stack([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])


def registry(cells, grid=Grid(), previous=None, release='20261111-rc1'):
    """cells maps (row,column) to original planning cloud; never duplicate geometry."""
    old = copy.deepcopy(previous) if previous else {'grid':asdict(grid),'fields':[]}
    if old['grid'] != asdict(grid):
        raise ValueError('Cannot change geometry in an existing registry')
    result=old['fields']; by_cell={r['grid_cell_id']:r for r in result}
    if len(by_cell)!=len(result) or len({r['field_id'] for r in result})!=len(result):
        raise ValueError('Duplicate field or cell identity')
    counters={c:max([int(r['field_id'].split('_')[1]) for r in result if r['planning_cloud']==c]+[0]) for c in ('L4','L5')}
    for (j,k),cloud in sorted(cells.items(),key=lambda item:(item[1],item[0])):
        if cloud not in counters:raise ValueError('Missing/invalid planning cloud')
        dec,n,step=grid.row(j)
        if not 0 <= k < n:raise ValueError('Invalid column')
        key=grid.key(j,k)
        if key in by_cell:continue
        if previous and release==previous['release']:
            raise ValueError('Appending cells requires a new release identifier')
        counters[cloud]+=1;ra=float(k*step)
        row={'field_id':f'TJ{cloud}_{counters[cloud]:03d}','grid_cell_id':key,'planning_cloud':cloud,
             'row_index':int(j),'column_index':int(k),'columns_in_full_row':int(n),
             'ra_deg':ra,'dec_deg':float(dec),'frame':'ICRS','pa_deg':0.,
             'orientation':'design +x north, +y east; PA of +x east of north',
             'parity':-1,'geometry_version':grid.version+'-TAN-inscribed',
             'nominal_boundary':{'model':'TAN','width_deg':grid.width,'height_deg':grid.width,
                                 'pa_deg':0.,'boundary_generator':'fixed_fields.boundary(samples=65)'},
             'usable_mask_ref':None,'template_common_region_ref':None,
             'created_release':release,'status':'CANDIDATE_PENDING_INSTRUMENT_ACCEPTANCE'}
        result.append(row);by_cell[key]=row
    return {'release':release,'grid':asdict(grid),'fields':result,
            'freeze_status':'nominal RC identities stable; formal instrument release pending masks/PA/repeatability'}


class FixedFields:
    def __init__(self, document):
        self.document=document;self.grid=Grid(**document['grid']);self.fields=document['fields']
        self.tree=cKDTree(unit_vectors([f['ra_deg'] for f in self.fields],[f['dec_deg'] for f in self.fields]))
        self.by_id={f['field_id']:f for f in self.fields}

    @classmethod
    def read(cls,path):
        with open(path) as f:return cls(json.load(f))

    def locate(self,ra,dec,quality=False):
        """All overlapping fields; includes both planning clouds, never infers object cloud."""
        points=unit_vectors(np.atleast_1d(ra),np.atleast_1d(dec))
        neighbors=self.tree.query_ball_point(points,2*np.sin(np.deg2rad(self.grid.width)/2))
        threshold=self.grid.margin+self.grid.guard if quality else 0
        return [[self.fields[i]['field_id'] for i in ns if self.fields[i]['status']!='RETIRED' and
                 edge_distance(a,d,self.fields[i]['ra_deg'],self.fields[i]['dec_deg'],self.grid)>=threshold]
                for a,d,ns in zip(np.atleast_1d(ra),np.atleast_1d(dec),neighbors)]

    def at(self,tracks,seconds,mag_limit=19.2):
        """Opt-in TJN query adapter; static geometry plus ephemeral object memberships."""
        from .catalog import interpolate
        from astropy.time import Time
        out={};jd=Time(seconds,format='unix').jd
        for k,t in tracks.items():
            v=float(np.interp(jd,t['datetime_jd'],t['V']))
            if v>mag_limit:continue
            ra,dec=interpolate(t,seconds)
            for fid in self.locate(ra,dec,quality=True)[0]:
                f=out.setdefault(fid,{**self.by_id[fid],'core_targets':[],'secondary_targets':[],
                                     'version':self.document['release'],
                                     'very_bright_star_risk':False,'bright_stars_within_cap':0})
                f['core_targets'].append(k)
        return sorted(out.values(),key=lambda f:f['field_id'])

    def template_key(self,fid,physical_filter,template_version,instrument='SHARP'):
        if physical_filter not in ('g','r','i','z'):raise ValueError('Explicit physical filter required')
        f=self.by_id[fid]
        return (instrument,f['geometry_version'],fid,physical_filter,template_version)
