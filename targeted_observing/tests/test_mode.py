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
        self.assertEqual(rows[0]['filter_change_s'],24)
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

    def test_closed_tour_includes_return_jump(self):
        from targeted_observing.scheduler import closed_neighbor_tour
        # Consecutive forward hops are legal; the last-to-first hop is not.
        self.assertFalse(closed_neighbor_tour([field(x,0) for x in [0,2.8,5.6,8.4]],self.c))
        self.assertTrue(closed_neighbor_tour([field(0,0),field(2.8,0),field(2.8,2.8),field(0,2.8)],self.c))
        self.assertTrue(closed_neighbor_tour([field(0,0),field(3.4,0)],self.c))

    def test_actual_short_tracklet_boundary(self):
        from targeted_observing.feedback import summarize
        from targeted_observing.geometry import isot
        a=unix('2026-11-03T10:30:00Z')
        for count,span,usable,complete in [(3,900,False,False),(4,719,False,False),(4,720,True,False),(5,1040,True,True)]:
            rows=[{'date':'2026-11-03','block_id':'EV_A','field_id':'test','visit':str(i+1),
                   'valid':'1','actual_start_utc':isot(a+t),'actual_end_utc':isot(a+t+30),'stable_star_rms_mag':''}
                  for i,t in enumerate(np.linspace(0,span,count))]
            result=summarize(rows)['fields'][0]
            self.assertEqual(result['usable_short_tracklet'],usable)
            self.assertEqual(result['complete_five_visits'],complete)
            for row in rows:row['valid']='0'
            self.assertFalse(summarize(rows)['fields'][0]['usable_short_tracklet'])

    def test_missing_filter_mapping_refuses_without_files(self):
        plan={'config':self.c,'blocks':[{'rows':[{'filter':'r'}]}]}
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'controller'
            with self.assertRaisesRegex(ValueError,'Verified filter'):export_controller_drafts(plan,p)
            self.assertFalse(p.exists())

    def test_empty_plan_clears_stale_generated_lists_only(self):
        from targeted_observing.products import write_products
        plan={'date':'2026-11-03','status':'LOCAL_SIMULATION_NOT_FOR_TELESCOPE','blocks':[],
              'skips':['TJ: disabled by zero nightly budget'],'free_intervals':[],
              'free_intervals_scope':'test'}
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for name in ['TJ.csv','EV_A.csv','exposures.csv','fields.csv','observer_notes.csv']:
                (p/name).write_text('previous run')
            write_products(plan,{'passed':True},p)
            self.assertFalse((p/'TJ.csv').exists())
            self.assertFalse((p/'EV_A.csv').exists())
            self.assertEqual((p/'exposures.csv').read_text(),'')
            self.assertEqual((p/'fields.csv').read_text(),'')
            self.assertEqual((p/'observer_notes.csv').read_text(),'previous run')

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

    def test_zero_budget_and_free_interval_partition(self):
        from unittest.mock import patch
        from targeted_observing.config import validate_config
        from targeted_observing.scheduler import make_plan
        from targeted_observing.geometry import Sky,isot
        self.c['trojan']['block_s']=0
        validate_config(self.c)
        sky=Sky(self.c);night=sky.night('2026-11-03')
        self.assertIsNone(trojan(sky,night['dark'],[],{}))
        a,b=night['evening']
        self.c['reservations']=[{'start_utc':isot(a+1800),'end_utc':isot(a+5400)}]
        block={'block_id':'EV_A','not_before_utc':isot(a),'hard_deadline_utc':isot(b)}
        with patch('targeted_observing.scheduler.twilight',side_effect=[[block],[]]):
            plan=make_plan(self.c,'2026-11-03')
        self.assertTrue(any('zero nightly budget' in r for r in plan['skips']))
        self.assertEqual(len(plan['free_intervals']),1)
        self.assertAlmostEqual(sum(r['duration_s'] for r in plan['free_intervals']),night['morning'][1]-a-5400,places=2)
        # Zero must be explicit; malformed or unsupported budgets still fail.
        for budget in [-1,1800,float('nan')]:
            self.c['trojan']['block_s']=budget
            with self.assertRaises(ValueError):validate_config(self.c)

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
