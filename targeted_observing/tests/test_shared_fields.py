import unittest

import numpy as np
from astropy.time import Time

from targeted_observing.audit_equatorial_tan import TanCoverage
from targeted_observing.full_night_test import command, read_lines
from targeted_observing.shared_fields import SharedFields


class SharedGridContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = SharedFields()

    def test_programs_share_geometry_and_text_suffix(self):
        fields = [self.catalog.resolve(prefix+'_0010') for prefix in ('MP', 'TJNL4', 'TJNL5', 'FOLLOW')]
        self.assertEqual({f['field_id'] for f in fields}, {'0010'})
        self.assertEqual(len({(f['ra_deg'], f['dec_deg']) for f in fields}), 1)
        for bad in ('MP_10', 'MP_0000', 'MP_4398', 'MP_TJL5_400', 'YSO_ORI', 'TJNL5_0010_extra'):
            with self.assertRaises(ValueError):
                self.catalog.resolve(bad)

    def test_all_centers_are_inside_their_own_quality_footprints(self):
        c = self.catalog
        hits = c.locate(c.table['center_ra'], c.table['center_dec'])
        self.assertTrue(all(fid in group for fid, group in zip(c.ids, hits)))
        self.assertEqual(c.locate([0, 360], [0, 0]), c.locate([360, 0], [0, 0]))

    def test_nominal_membership_agrees_with_existing_tan_audit(self):
        c = self.catalog
        model = TanCoverage(c.table, c.manifest['tan_width_east_deg'],
                            c.manifest['tan_height_north_deg'], c.manifest['clockwise_roll_deg'])
        rng = np.random.default_rng(18)
        ra = np.r_[rng.uniform(0, 360, 600), 0, 360, 0, 90]
        dec = np.r_[rng.uniform(-35, 90, 600), 85, 85, 90, -90]
        self.assertEqual([bool(group) for group in c.locate(ra, dec, quality=False)],
                         (model.margin(ra, dec) >= 0).tolist())

    def test_controller_roundtrip_preserves_shared_center(self):
        c = self.catalog
        for prefix in ('MP', 'TJNL5'):
            f = c.resolve(prefix+'_0010')
            text = command({**f, 'field_id': f['execution_name']}, 'r', 5, 120)
            row = read_lines(text, 'TEST')[0]
            self.assertEqual(row['field_id'], prefix+'_0010')
            self.assertAlmostEqual(row['ra_deg'], f['ra_deg'], places=5)
            self.assertAlmostEqual(row['dec_deg'], f['dec_deg'], places=5)

    def test_track_crossing_is_not_assigned_only_by_first_position(self):
        c = self.catalog
        f, other = c.by_id['2000'], c.by_id['2010']
        times = np.array([1800000000., 1800003600., 1800007200.])
        track = dict(datetime_jd=Time(times, format='unix').jd.tolist(),
                     RA=[f['ra_deg'], other['ra_deg'], other['ra_deg']],
                     DEC=[f['dec_deg'], other['dec_deg'], other['dec_deg']], V=[18., 18., 18.])
        self.assertEqual(c.memberships({'test': track}, times, {'test': 'L5'}), [])
        with self.assertRaisesRegex(ValueError, 'coverage'):
            c.memberships({'test': track}, times+86400, {'test': 'L5'})

    def test_mp_subset_preserves_science_latitude_cut(self):
        c = self.catalog
        fields, provenance = c.twilight_grid()
        self.assertEqual(len(fields), int(np.count_nonzero(abs(c.table['ecliptic_lat_deg']) <= 20)))
        self.assertEqual(provenance['catalog_sha256'], c.manifest['catalog_sha256'])
        for f in fields:
            identity = c.resolve('MP_'+f['field_id'][2:])
            self.assertEqual((identity['ra_deg'], identity['dec_deg']), (f['ra_deg'], f['dec_deg']))


if __name__ == '__main__':
    unittest.main()
