import unittest
from targeted_observing.catalog import BASE
from targeted_observing.full_night_test import read_lines, replay, movement


class FullNightContracts(unittest.TestCase):
    def test_supplied_m31_counts_and_exposures(self):
        text=(BASE/'data/night_inputs/20261111_M31_GRIZ.txt').read_text()
        commands=read_lines(text,'M31')
        self.assertEqual(sum(r['count'] for r in commands),104)
        self.assertEqual(sum(r['count']*r['exposure_s'] for r in commands),10245)
        self.assertEqual([r['filter'] for r in commands],['G']*3+['R']*3+['I']*3+['Z']*2)

    def test_repeat_counts_and_filter_only_change_do_not_add_slew(self):
        text=('A_G 00:42:44.330 +41:16:07.50 0 G 3 5 0\n'
              'A_R 00:42:44.330 +41:16:07.50 0 R 1 60 0\n')
        rows=replay(read_lines(text,'TEST'),0)
        self.assertEqual(len(rows),4)
        self.assertEqual(rows[-1]['release_unix'],245)
        self.assertEqual(rows[-1]['slew_s'],0)
        self.assertEqual(rows[-1]['filter_change_s'],10)

    def test_same_name_different_coordinates_requires_slew(self):
        a={'field_id':'A','ra_deg':359,'dec_deg':20}
        b={'field_id':'A','ra_deg':1,'dec_deg':20}
        self.assertGreater(movement(a,b),20)
        self.assertEqual(movement(a,a),0)


if __name__=='__main__':unittest.main()
