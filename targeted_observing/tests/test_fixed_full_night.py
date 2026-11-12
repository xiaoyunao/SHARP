import fnmatch
import re
import unittest

from targeted_observing.controller_fields import controller_name, canonical_name, is_fixed_name
from targeted_observing.fixed_full_night import make_command, fill_to_handover, rows_for, compress_dwells


class FixedControllerTests(unittest.TestCase):
    def test_alias_preserves_cloud_number_and_matches_pipeline_globs(self):
        for field in ('TJL4_001','TJL5_436','TJL5_1234'):
            name=controller_name(field)
            self.assertEqual(canonical_name(name),field)
            self.assertTrue(is_fixed_name(name))
            filename='OBJ_'+name+'_0031_cat.fits.gz'
            self.assertTrue(fnmatch.fnmatch(filename,'OBJ_MP_*_cat.fits.gz'))
            self.assertIsNotNone(re.fullmatch(r'.*MP.*\.(fits|fits\.gz)$',filename))
        self.assertFalse(is_fixed_name('MP_1570'))
        for bad in ('MP_436','TJL6_001','TJL5_436_extra','MP_TJL5_436'):
            with self.assertRaises(ValueError):controller_name(bad)

    def test_five_frame_dwells_switch_only_after_complete_visits(self):
        from itertools import groupby
        a={'field_id':'MP_TJL5_400','ra_deg':79.4117647,'dec_deg':20.}
        b={'field_id':'MP_TJL5_436','ra_deg':77.9104478,'dec_deg':22.5}
        nxt=make_command({'field_id':'YSO_ORI','ra_deg':84.,'dec_deg':-4.8},'g',60,'YSO',72)
        commands=[make_command(a,'r',5,'EV')]
        fill_to_handover(commands,0,[a,b],nxt,4000,'TJN_MAIN',dwell=5)
        rr=rows_for(commands,0)
        groups=[list(g) for _,g in groupby([r for r in rr if r['block']=='TJN_MAIN'],key=lambda r:r['field_id'])]
        self.assertGreater(len(groups),3)
        for group in groups[:-1]:
            self.assertEqual(len(group),5)
            self.assertTrue(all(r['exposure_s']==120 for r in group))
        for group in groups:
            self.assertTrue(all(r['slew_s']==0 for r in group[1:]))
        self.assertTrue(all(5<=r['exposure_s']<=120 for group in groups for r in group))
        exported=compress_dwells(commands)
        self.assertEqual(exported[1]['count'],5)
        expanded=rows_for(exported,0)
        self.assertEqual([(r['field_id'],r['start_unix'],r['exposure_s']) for r in expanded],
                         [(r['field_id'],r['start_unix'],r['exposure_s']) for r in rr])
        self.assertTrue(all(c['count']==1 for c in commands))
        for invalid in (0,-1,1.5,True):
            with self.assertRaises(ValueError):fill_to_handover([],0,[a,b],nxt,4000,'TJN_MAIN',dwell=invalid)

    def test_small_tail_and_wrap_handover_have_no_wait_or_short_invalid_exposure(self):
        field={'field_id':'MP_TJL5_436','ra_deg':359.,'dec_deg':20.}
        nxt={'field_id':'YSO_ORI','ra_deg':1.,'dec_deg':20.}
        next_command=make_command(nxt,'g',60,'YSO',72)
        for extra in (0,1,5,19,20,134,150,.6):
            commands=[make_command(field,'r',5,'EV')]
            base=rows_for(commands+[make_command(field,'r',120,'TJN'),next_command],0)
            target=next(r['start_unix'] for r in base if r['block']=='YSO')+extra
            fill_to_handover(commands,0,[field],next_command,target,'TJN')
            rows=rows_for(commands+[next_command],0)
            got=next(r['start_unix'] for r in rows if r['block']=='YSO')
            self.assertGreaterEqual(target-got,0)
            self.assertLess(target-got,1.01)
            self.assertTrue(all(5<=r['exposure_s']<=120 for r in rows if r['block']=='TJN'))
            self.assertTrue(all(abs(a['release_unix']-b['overhead_start_unix'])<1e-8 for a,b in zip(rows,rows[1:])))


if __name__=='__main__':unittest.main()
