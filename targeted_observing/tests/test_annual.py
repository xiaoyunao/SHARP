import json
import unittest
import numpy as np
from targeted_observing.catalog import BASE
from targeted_observing.annual_sources import dates
from targeted_observing.annual_weather import ensemble,lightcurve_quality
from targeted_observing.annual_twilight import AnchorState


class AnnualContracts(unittest.TestCase):
    def test_full_year_and_leap_calendar(self):
        self.assertEqual(len(dates('2026-11-01','2027-11-01')),365)
        self.assertEqual(len(dates('2024-01-01','2025-01-01')),366)

    def test_shared_weather_is_reproducible_and_correlated(self):
        nights=dates('2026-11-01','2027-11-01')
        a=ensemble(nights,200,42);b=ensemble(nights,200,42)
        for key in a:np.testing.assert_array_equal(a[key],b[key])
        np.testing.assert_array_equal(a['observable'][:,:,0],a['observable'][:,:,2])
        self.assertGreater(np.corrcoef(a['observable'][:,:-1,0].ravel(),a['observable'][:,1:,0].ravel())[0,1],.25)

    def test_quality_is_not_point_count_or_period(self):
        q=json.loads((BASE/'annual_defaults.json').read_text())['quality']
        ns=np.repeat(np.arange(10)*3,10);ts=ns*86400+np.tile(np.linspace(0,3600,10),10)
        err=np.full(100,.02);zp=np.full(100,.01);valid=np.ones(100,dtype=bool)
        result=lightcurve_quality(ts,ns,err,zp,valid,q)
        self.assertTrue(result['quality_lightcurve']);self.assertIsNone(result['period_reliable'])
        self.assertFalse(lightcurve_quality(ts,ns,np.full(100,.04),zp,valid,q)['quality_lightcurve'])
        # Identical point count concentrated into fewer nights fails.
        self.assertFalse(lightcurve_quality(ts,ns//6,err,zp,valid,q)['quality_lightcurve'])

    def test_anchor_ages_only_on_successful_feedback(self):
        state=AnchorState();b={'block_id':'EV_A','fields':[{'field_id':str(i)} for i in range(4)]}
        ids={str(i) for i in range(4)}
        state.feedback([b],'night0',set())
        self.assertEqual(len(state.previous('EV')),1)
        for night in ['n1','n2','n3']:state.feedback([b],night,ids)
        self.assertEqual(state.previous('EV'),[])
        self.assertEqual(state.previous('AM'),[])

    def test_numeric_execution_matches_scheduler_and_accounts_time(self):
        from targeted_observing.annual_simulation import numeric_visits
        from targeted_observing.scheduler import field,visits
        from targeted_observing.geometry import unix
        c=json.loads((BASE/'defaults.json').read_text());c['trojan']['exposures_per_dwell']=2
        fields=[field(359,20,'TJ'),field(15,28,'TJ'),field(19,22,'TJ')]
        a=1800000000.
        expected=visits(fields,a,a+6900,120,100,c,'TJ','trojan')
        actual=numeric_visits(fields,2,c)
        self.assertEqual(len(actual),len(expected))
        np.testing.assert_allclose(actual[:,1]+a,[unix(r['mid_utc']) for r in expected],atol=.001,rtol=0)
        delayed=numeric_visits(fields,2,c,180,1,3)
        self.assertLessEqual(delayed[-1,2],6900)
        self.assertEqual(len(delayed)%6,0)

    def test_period_recovers_asymmetric_double_peaked_signal(self):
        from targeted_observing.annual_periods import injected,search
        rng=np.random.default_rng(1);t=np.sort(rng.uniform(0,800,200));error=np.full(200,.02)
        result=search(t,injected(t,7.31,.3,.8)+rng.normal(0,error),error)
        self.assertLess(abs(result['period_h']-7.31),.01)
        self.assertTrue(result['unique_in_searched_candidates'])
        self.assertLess(abs(result['amplitude_mag']-.3),.05)

    def test_spherical_area_does_not_count_repeat_exposures_twice(self):
        from targeted_observing.annual_coverage import pixels,window_pixels
        from targeted_observing.scheduler import field
        from targeted_observing.geometry import isot
        import healpy as hp
        f=field(359.9,60);pix=pixels(f['ra_deg'],f['dec_deg'],1024)
        self.assertAlmostEqual(len(pix)*hp.nside2pixarea(1024,degrees=True),9,delta=.1)
        b={'fields':[f],'rows':[{'field_id':f['field_id'],'mid_utc':isot(1800000000+i*300)} for i in range(5)]}
        area,qualified=window_pixels([b],512)
        self.assertLess(len(qualified),len(area))
        self.assertTrue(set(qualified).issubset(set(area)))
        repeated,_=window_pixels([b,b],512)
        np.testing.assert_array_equal(area,repeated)

    def test_shifted_exposure_checks_all_enclosing_grid_nodes(self):
        from targeted_observing.annual_timing import mask
        block={'timing_grid':[0,30,60,90,120],'timing_valid':{'A':[True,True,False,True,True]}}
        self.assertFalse(mask(block,['A'],np.array([60.]))[0])
        self.assertFalse(mask(block,['A'],np.array([-10.]))[0])
        self.assertTrue(mask(block,['A'],np.array([15.]))[0])

    def test_queue_releases_baseline_wait_then_returns_to_finish(self):
        from targeted_observing.annual_simulation import pick
        fields=[{'core_targets':['old']},{'core_targets':['new']}]
        menu={'fields':fields};data={'valid':np.ones((2,1),dtype=bool)}
        history={'old':{'times':[0],'night_indices':[0],'quality':{'good_points':120,'good_nights':10,'first_good_night':0}}}
        completion=({'old':20,'new':20},{'old':60,'new':60},10)
        self.assertEqual(pick(menu,data,{'old':1,'new':.9},1,history,18,completion)[1],[1])
        self.assertEqual(pick(menu,data,{'old':1,'new':.9},1,history,21,completion)[1],[0])

    def test_accepted_curve_survives_later_noisy_companion_points(self):
        from targeted_observing.annual_simulation import update_quality
        q=json.loads((BASE/'annual_defaults.json').read_text())['quality']
        nights=np.repeat(np.arange(10)*3,10);times=nights*86400+np.tile(np.linspace(0,3600,10),10)
        s={'times':times.tolist(),'night_indices':nights.tolist(),'errors':[.02]*100,'zeropoints':[.01]*100,'valid':[True]*100}
        update_quality(s,q);self.assertTrue(s['quality']['quality_lightcurve'])
        for key,value in [('times',40*86400),('night_indices',40),('errors',.2),('zeropoints',.01),('valid',True)]:s[key].extend([value]*100)
        update_quality(s,q)
        self.assertTrue(s['quality']['quality_lightcurve']);self.assertEqual(s['completion_index'],100)
        self.assertEqual(len(s['times']),200)
