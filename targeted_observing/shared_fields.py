"""Frozen common pointing identities for MP and TJN; no production imports."""
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time
from scipy.spatial import cKDTree

from .catalog import interpolate
from .fixed_fields import tan_xy, from_tan, unit_vectors

ROOT = Path(__file__).resolve().parent / 'data/shared_fields_20261116'
NAME = re.compile(r'(MP|TJNL4|TJNL5)_([0-9]{4})')


class SharedFields:
    def __init__(self):
        self.manifest = json.loads((ROOT / 'manifest.json').read_text())
        path = ROOT / self.manifest['catalog']
        if hashlib.sha256(path.read_bytes()).hexdigest() != self.manifest['catalog_sha256']:
            raise ValueError('Shared catalog checksum changed')
        self.table = Table.read(path)
        ids = [str(x) for x in self.table['field_id']]
        if ids != [f'{i:04d}' for i in range(1, 4398)]:
            raise ValueError('Invalid shared field identities')
        self.by_id = {str(r['field_id']): {
            'field_id': str(r['field_id']), 'ra_deg': float(r['center_ra']),
            'dec_deg': float(r['center_dec']), 'frame': 'ICRS',
            'grid_release': self.manifest['release'],
        } for r in self.table}
        self.ids = ids
        self.tree = cKDTree(unit_vectors(self.table['center_ra'], self.table['center_dec']))

    def name(self, prefix, field_id):
        if prefix not in ('MP', 'TJNL4', 'TJNL5') or field_id not in self.by_id:
            raise ValueError('Expected a program prefix and a four-character shared ID')
        return prefix + '_' + field_id

    def resolve(self, name):
        match = NAME.fullmatch(name)
        if match is None or match[2] not in self.by_id:
            raise ValueError('Not a shared-grid observation name: ' + name)
        return {**self.by_id[match[2]], 'program': match[1], 'execution_name': name}

    def edge_distance(self, ra, dec, field_id):
        f = self.by_id[field_id]
        x, y = tan_xy(ra, dec, f['ra_deg'], f['dec_deg'])
        angle = np.deg2rad(self.manifest['clockwise_roll_deg'])
        # Match audit_equatorial_tan.TanCoverage's rotated east/north basis.
        xx, yy = np.cos(angle)*x-np.sin(angle)*y, np.sin(angle)*x+np.cos(angle)*y
        hx, hy = np.tan(np.deg2rad([
            self.manifest['tan_width_east_deg']/2,
            self.manifest['tan_height_north_deg']/2]))
        return np.minimum(np.rad2deg(hx)-abs(xx), np.rad2deg(hy)-abs(yy))

    def locate(self, ra, dec, quality=True):
        ra, dec = np.broadcast_arrays(np.atleast_1d(ra), np.atleast_1d(dec))
        if not np.all(np.isfinite(ra) & np.isfinite(dec) & (abs(dec) <= 90)):
            raise ValueError('Invalid sky coordinates')
        neighbors = self.tree.query_ball_point(unit_vectors(ra, dec), 2*np.sin(np.deg2rad(2.2)/2))
        margin = (self.manifest['quality_margin_deg'] + self.manifest['quality_guard_deg']) if quality else 0
        return [sorted(self.ids[i] for i in ns
                       if self.edge_distance(a, d, self.ids[i]) >= margin)
                for a, d, ns in zip(ra, dec, neighbors)]

    def twilight_grid(self, ecliptic_limit=20):
        """Keep the existing MP science latitude cut; change pointing grid only."""
        rows = self.table[abs(self.table['ecliptic_lat_deg']) <= ecliptic_limit]
        fields = []
        angle = np.deg2rad(self.manifest['clockwise_roll_deg'])
        hx, hy = np.tan(np.deg2rad([
            self.manifest['tan_width_east_deg']/2,
            self.manifest['tan_height_north_deg']/2]))
        xx, yy = np.array([-hx, hx, hx, -hx]), np.array([-hy, -hy, hy, hy])
        x, y = np.cos(angle)*xx+np.sin(angle)*yy, -np.sin(angle)*xx+np.cos(angle)*yy
        radius = float(np.rad2deg(np.arctan(np.hypot(hx, hy))))
        for r in rows:
            f = self.by_id[str(r['field_id'])]
            ra, dec = from_tan(f['ra_deg'], f['dec_deg'], np.rad2deg(x), np.rad2deg(y))
            # The existing twilight planner's internal MPdddd adapter exports MP_dddd.
            fields.append({**f, 'field_id': 'MP'+f['field_id'], 'radius_deg': radius,
                           'corners_ra_deg': ra.tolist(), 'corners_dec_deg': dec.tolist()})
        return fields, {**self.manifest, 'mp_ecliptic_limit_deg': ecliptic_limit,
                        'mp_candidate_fields': len(fields)}

    def memberships(self, tracks, times, clouds, mag_limit=19.2, retained_fraction=.8):
        """All overlapping fields retaining a target over the supplied visit times.

        This is candidate membership, not a telescope schedule or a detection.
        The caller must query again at the final scheduled exposure midpoints.
        """
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or not len(times) or not np.all(np.isfinite(times)):
            raise ValueError('Explicit finite visit times required')
        jd = Time(times, format='unix').jd
        output = []
        for key, track in tracks.items():
            cloud = clouds[key]
            if cloud not in ('L4', 'L5'):
                raise ValueError('Missing frozen cloud classification: '+key)
            ra, dec = interpolate(track, times)  # Reject out-of-range cached ephemerides.
            mag = np.interp(jd, track['datetime_jd'], track['V'])
            if mag.min() > mag_limit:
                continue
            hits = self.locate(ra, dec)
            for fid in sorted({fid for group in hits for fid in group}):
                good = np.array([fid in group for group in hits]) & (mag <= mag_limit)
                if good.mean() >= retained_fraction:
                    output.append(dict(object_id=key, cloud=cloud, field_id=fid,
                        execution_name=self.name('TJN'+cloud, fid), **{
                            k: self.by_id[fid][k] for k in ('ra_deg', 'dec_deg')},
                        retained_fraction=float(good.mean()), V_min=float(mag.min()),
                        V_max=float(mag.max())))
        return output


def cloud_labels():
    path = ROOT.parent / 'tjn_handoff_20261111/inputs/tjn_annual_candidates_1063.csv'
    with path.open() as stream:
        return {r['object_id']: r['cloud_proxy'] for r in csv.DictReader(stream)}
