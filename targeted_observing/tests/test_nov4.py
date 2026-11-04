"""Regressions identified by the second independent review."""
import copy
import json
import unittest
import numpy as np
from astropy.time import Time
from targeted_observing.catalog import BASE
from targeted_observing.geometry import Sky, setup_iers, unix, isot
from targeted_observing.scheduler import field, visits, trojan, twilight
from targeted_observing.execution import gate
from targeted_observing.validate import validate_plan


class ReviewRegressions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (BASE/'cache/iers.all').exists():raise unittest.SkipTest('IERS required')
        setup_iers(BASE/'cache/iers.all')

    def setUp(self):self.c=json.loads((BASE/'defaults.json').read_text())

    def test_aabb_actual_budget_and_nonuniform_cadence(self):
        self.c['trojan']['exposures_per_dwell']=2
        self.c['use_header_slew_floor']=False
        a=unix('2026-11-03T15:30:00Z')
        fs=[field(55,20,'TJ'),field(72,20,'TJ')]
        b=trojan(Sky(self.c),(a,a+7200),fs,{})
        self.assertEqual(len(b['fields']),2)
        for f in fs:
            rows=[r for r in b['rows'] if r['field_id']==f['field_id']]
            self.assertEqual(len(rows),20)
            self.assertEqual([r['visit'] for r in rows],list(range(1,21)))
            self.assertEqual(set(np.diff([unix(r['mid_utc']) for r in rows])),{135.,525.})
        campaign={'start_date':'2026-11-03','stop_date':'2026-11-17','fields':fs,'tracks':{}}
        self.assertTrue(validate_plan({'config':self.c,'date':'2026-11-03','blocks':[b]},campaign)['passed'])

    def test_delayed_first_exposure_updates_target_state(self):
        a=unix('2026-11-03T15:30:00Z');f=field(55,20,'TJ');f['core_targets']=['edge']
        rows=visits([f],a,a+4200,120,30,self.c,'TJ','trojan')
        ts=np.array([unix(r['mid_utc']) for r in rows])
        # Moving outward at first frame, then inward: delay changes only first flag.
        track={'datetime_jd':Time([ts[0]-120,ts[0],ts[0]+60,ts[1],ts[-1]+120],format='unix').jd.tolist(),
               'RA':[55,55,57,55,55],'DEC':[20]*5}
        for r in rows:r['target_states']={'edge':{'in_field':True,'photometric_usable':None}}
        b={'block_id':'TJ','mode':'trojan','fields':[f],'rows':rows,'extended_search':False,'hard_deadline_utc':isot(a+4800)}
        out=gate({'config':self.c,'blocks':[b]},'TJ',0,isot(unix(rows[0]['start_utc'])+60),{'tracks':{'edge':track}})
        self.assertEqual(out['action'],'dry_run_exposure_eligible')
        self.assertFalse(out['exposure']['target_states']['edge']['in_field'])
        self.assertGreater(out['target_sampling_at_actual_times'][f['field_id']][f['field_id']+':edge']['geometric_retention'],.9)

    def test_old_field_revalidated_against_stricter_beta(self):
        old=json.loads((BASE/'data/anchor_regression_20261104.json').read_text())['blocks']
        sky=Sky(self.c);blocks=twilight(sky,sky.night('2026-11-06')['evening'],'EV',old)
        # Make every nonzero-latitude existing field illegal without redesigning it.
        strict=copy.deepcopy(self.c);strict['twilight']['beta_max_deg']=0
        b=blocks[0];r=b['rows'][0]
        q=Sky(strict).check_field(b['fields'][0],unix(r['start_utc']),unix(r['end_utc']),'twilight',b['extended_search'])
        self.assertIn('ecliptic_latitude',q['reasons'])
        result=validate_plan({'config':strict,'date':'2026-11-06','blocks':blocks})
        self.assertTrue(any('ecliptic latitude' in x for x in result['errors']))

    def test_manual_control_preserves_conflicting_quality(self):
        from targeted_observing.science import assess
        review=json.loads((BASE/'data/manual_science_review.json').read_text())['2223']
        evidence={'LCDB':{'period_h':22.741,'U':'2+'},'manual_review':review}
        result=assess(evidence)
        self.assertEqual(result['role'],'control')
        self.assertEqual(result['periods'][0]['U'],'2+')
        self.assertEqual(result['manual_review']['measurements'][0]['paper_quality'],3)
