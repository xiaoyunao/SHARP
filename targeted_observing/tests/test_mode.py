import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from targeted_observing.catalog import BASE,interpolate
from targeted_observing.geometry import coord,offset,in_field,setup_iers,unix
from targeted_observing.scheduler import field,slew,visits,free_intervals,trojan
from targeted_observing.products import export_controller_drafts


class CoreTests(unittest.TestCase):
    def setUp(self):self.c=json.loads((BASE/'defaults.json').read_text())

    def test_ra_wrap_and_high_declination(self):
        f=coord(359.9,75);p=offset(f,[0,1.3,1.6],[0,1.3,0])
        np.testing.assert_array_equal(in_field(p.ra.deg,p.dec.deg,359.9,75),[True,True,False])
        self.assertLess(slew(field(359.9,0),field(.1,0),self.c),30)

    def test_budget_and_no_partial_cycle(self):
        fs=[field(10,20),field(12.9,20)]
        rows=visits(fs,1800000000,1800006900,120,22,self.c,'TJ','trojan')
        self.assertEqual(len(rows)%2,0)
        self.assertLessEqual(unix(rows[-1]['release_utc']),1800006900)
        self.assertEqual(rows[0]['filter_change_s'],10)
        self.assertTrue(all(r['filter_change_s']==0 for r in rows[1:]))
        self.assertGreater(rows[1]['slew_s'],20)  # Calibrated floor replaces optimistic nominal timing.

    def test_long_slew_and_repeated_field(self):
        a=field(10,10);self.assertEqual(slew(a,a,self.c),0)
        self.assertGreaterEqual(slew(a,field(120,20),self.c),60)

    def test_reservations(self):
        a=unix('2026-11-03T12:00:00Z');b=a+7200
        rs=[{'start_utc':'2026-11-03T12:30:00Z','end_utc':'2026-11-03T13:30:00Z'}]
        self.assertEqual(free_intervals((a,b),rs),[(a,a+1800),(a+5400,b)])
        self.assertEqual(free_intervals((a,b),[{'start_utc':'2026-11-03T11:00:00Z','end_utc':'2026-11-03T15:00:00Z'}]),[])

    def test_missing_filter_mapping_refuses_without_files(self):
        plan={'config':self.c,'blocks':[{'rows':[{'filter':'r'}]}]}
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'controller'
            with self.assertRaisesRegex(ValueError,'Verified filter'):export_controller_drafts(plan,p)
            self.assertFalse(p.exists())

    def test_ephemeris_wrap_and_no_extrapolation(self):
        from astropy.time import Time
        t={'datetime_jd':[2461348.5,2461349.5],'RA':[359.,1.],'DEC':[10.,10.]}
        ra,dec=interpolate(t,Time(2461349.,format='jd').unix)
        self.assertAlmostEqual(ra[0]%360,0)
        with self.assertRaisesRegex(ValueError,'coverage'):interpolate(t,0)


class IntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (BASE/'cache/iers.all').exists():raise unittest.SkipTest('Run fetch_sources first for offline geometry tests')
        setup_iers(BASE/'cache/iers.all')

    def setUp(self):self.c=json.loads((BASE/'defaults.json').read_text())

    def test_solar_windows_and_season(self):
        from targeted_observing.geometry import Sky,local_time
        s=Sky(self.c);n=s.night('2026-11-03')
        self.assertEqual(local_time(n['evening'][0])[:16],'2026-11-03T18:05')
        self.assertEqual(local_time(n['morning'][1])[:16],'2026-11-04T05:41')
        for d in ['2026-12-21','2027-03-20','2027-06-21','2027-09-23']:
            n=s.night(d);self.assertAlmostEqual(n['evening'][1]-n['evening'][0],3600)
            self.assertLess(n['dark'][0],n['dark'][1])

    def test_footprint_rejection_and_sun(self):
        from targeted_observing.geometry import Sky
        s=Sky(self.c);a=unix('2026-11-03T14:00:00Z')
        self.assertFalse(s.check_field(field(0,-80),a,a+120,'trojan')['ok'])
        self.assertFalse(s.check_field(field(0,40),unix('2026-11-03T04:00:00Z'),unix('2026-11-03T04:02:00Z'),'trojan')['ok'])

    def test_trojan_fallback_and_reservations(self):
        from targeted_observing.geometry import Sky
        s=Sky(self.c);a=unix('2026-11-03T15:30:00Z');fs=[field(55,20,'TJ'),field(57.9,20,'TJ')]
        for minutes,count in [(120,2),(100,2),(75,1),(59,0)]:
            b=trojan(s,(a,a+minutes*60),fs,{})
            self.assertEqual(len(b['fields']) if b else 0,count)
        c=copy.deepcopy(self.c);c['reservations']=[{'start_utc':'2026-11-03T15:00:00Z','end_utc':'2026-11-03T19:00:00Z'}]
        self.assertIsNone(trojan(Sky(c),(a,a+7200),fs,{}))

    def test_independent_audit_detects_corrupt_exposure(self):
        from targeted_observing.geometry import Sky
        from targeted_observing.validate import validate_plan
        s=Sky(self.c);a=unix('2026-11-03T15:30:00Z');fs=[field(55,20,'TJ'),field(57.9,20,'TJ')]
        b=trojan(s,(a,a+7200),fs,{})
        plan={'config':self.c,'date':'2026-11-03','blocks':[b]}
        campaign={'start_date':'2026-11-03','stop_date':'2026-11-17','tracks':{}}
        self.assertTrue(validate_plan(plan,campaign)['passed'])
        plan['blocks'][0]['rows'][-1]['end_utc']=plan['blocks'][0]['rows'][-1]['start_utc']
        self.assertFalse(validate_plan(plan,campaign)['passed'])

    def test_execution_gate_early_and_deadline(self):
        from targeted_observing.geometry import Sky,isot
        from targeted_observing.execution import gate
        s=Sky(self.c);a=unix('2026-11-03T15:30:00Z');fs=[field(55,20,'TJ'),field(57.9,20,'TJ')]
        b=trojan(s,(a,a+7200),fs,{})
        plan={'config':self.c,'date':'2026-11-03','blocks':[b]};campaign={'tracks':{}}
        self.assertEqual(gate(plan,'TJ',0,isot(a),campaign)['action'],'wait')
        self.assertEqual(gate(plan,'TJ',0,b['rows'][0]['start_utc'],campaign)['action'],'dry_run_exposure_eligible')
        self.assertEqual(gate(plan,'TJ',0,isot(a+7100),campaign)['action'],'abort_remaining_block')

if __name__=='__main__':unittest.main()
