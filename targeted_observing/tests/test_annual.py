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
