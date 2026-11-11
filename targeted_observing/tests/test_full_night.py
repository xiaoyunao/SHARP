import unittest
from targeted_observing.catalog import BASE
from targeted_observing.full_night_test import read_lines, replay, movement, lowercase_filters, TJ, YSO, YSO_FIELDS


class FullNightContracts(unittest.TestCase):
    def test_supplied_m31_counts_and_exposures(self):
        text=(BASE/'data/night_inputs/20261111_M31_GRIZ.txt').read_text()
        commands=read_lines(text,'M31')
        self.assertEqual(sum(r['count'] for r in commands),104)
        self.assertEqual(sum(r['count']*r['exposure_s'] for r in commands),10245)
        self.assertEqual([r['filter'] for r in commands],['g']*3+['r']*3+['i']*3+['z']*2)

    def test_repeat_counts_and_filter_only_change_do_not_add_slew(self):
        text=('A_G 00:42:44.330 +41:16:07.50 0 G 3 5 0\n'
              'A_R 00:42:44.330 +41:16:07.50 0 R 1 60 0\n')
        rows=replay(read_lines(text,'TEST'),0)
        self.assertEqual(len(rows),4)
        self.assertEqual(rows[-1]['release_unix'],245)
        self.assertEqual(rows[-1]['slew_s'],0)
        self.assertEqual(rows[-1]['filter_change_s'],10)

    def test_field_names_and_filter_normalization_preserve_other_columns(self):
        self.assertEqual(TJ['field_id'],'TJN_17365')
        self.assertEqual(YSO['field_id'],'YSO_ORI')
        self.assertEqual(YSO_FIELDS['YSO_CYG'],{'ra_deg':314.,'dec_deg':44.1})
        text=(BASE/'data/night_inputs/20261111_M31_GRIZ.txt').read_text()
        converted=lowercase_filters(text)
        for old,new in zip(text.splitlines(),converted.splitlines()):
            a,b=old.split(),new.split()
            self.assertEqual(b[4],a[4].lower())
            self.assertEqual(a[:4]+a[5:],b[:4]+b[5:])

    def test_same_name_different_coordinates_requires_slew(self):
        a={'field_id':'A','ra_deg':359,'dec_deg':20}
        b={'field_id':'A','ra_deg':1,'dec_deg':20}
        self.assertGreater(movement(a,b),20)
        self.assertEqual(movement(a,a),0)


if __name__=='__main__':unittest.main()
