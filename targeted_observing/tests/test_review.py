import json
import unittest
import numpy as np
from astropy.time import Time
from targeted_observing.catalog import BASE
from targeted_observing.geometry import setup_iers,unix,isot,Sky
from targeted_observing.scheduler import field,trojan,twilight
from targeted_observing.science import assess
from targeted_observing.retention import report

class EvidenceTests(unittest.TestCase):
    def test_damit_discrepancy_is_evidence_not_priority_boost(self):
        e={'LCDB':{'period_h':8.47,'U':'2'},'DAMIT':[{'period_h':10.3792}]}
        a=assess(e)
        self.assertEqual({p['source'] for p in a['periods']},{'LCDB','DAMIT'})
        self.assertTrue(a['review_reasons']);self.assertEqual(a['weight'],.1)
        self.assertEqual(assess(None)['role'],'literature_unchecked')

    def test_tess_full_table_and_ambiguous_period(self):
        from targeted_observing.review_data import tess_records
        path=BASE/'cache/tess2025.tsv'
        if not path.exists():self.skipTest('Run review_data for source integration test')
        rows,meta=tess_records(path)
        self.assertEqual((meta['rows'],meta['numbered'],meta['zero_width']),(9912,9868,34))
        self.assertGreater(rows['22180']['upper_h'],500)
        self.assertTrue(rows['24453']['other_periods_h'])
        self.assertNotIn('18054',rows)

    def test_retention_uses_visits_not_daytime_grid(self):
        c=json.loads((BASE/'defaults.json').read_text());f=field(0,0,'TJ');f['core_targets']=['1']
        ts=unix('2026-11-03T15:00:00Z')+np.arange(6)*120
        t={'datetime_jd':Time(ts,format='unix').jd.tolist(),'RA':[0,0,0,0,2,2],'DEC':[0]*6}
        rows=[{'field_id':f['field_id'],'mid_utc':isot(x)} for x in ts[:5]]
        r=report([f],{'1':t},{'night':rows},c)[f['field_id']+':1']
        self.assertEqual(r['geometric_retention'],.8)
        self.assertIsNone(r['photometric_usable_fraction'])
        self.assertEqual(r['nights_with_geometric_points'],1)

class GeometryReviewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (BASE/'cache/iers.all').exists():raise unittest.SkipTest('IERS required')
        setup_iers(BASE/'cache/iers.all')
    def setUp(self):self.c=json.loads((BASE/'defaults.json').read_text())

    def test_single_field_fills_budget(self):
        camp=json.loads((BASE/'data/campaign_20261103.json').read_text());a=unix('2026-11-03T17:00:00Z')
        for minutes,n in [(75,30),(89,36)]:
            b=trojan(Sky(self.c),(a,a+minutes*60),camp['fields'],camp['tracks'])
            self.assertIsNotNone(b);self.assertEqual(len(b['rows']),n)
            self.assertLessEqual(unix(b['rows'][-1]['release_utc']),a+minutes*60-300)

    def test_old_anchor_regression_retains_seven(self):
        old=json.loads((BASE/'data/anchor_regression_20261104.json').read_text())['blocks']
        sky=Sky(self.c);blocks=twilight(sky,sky.night('2026-11-06')['evening'],'EV',old)
        self.assertEqual(len(blocks),2)
        self.assertGreaterEqual(sum(b['retained_anchor_count'] for b in blocks),7)

    def test_hipparcos_catches_betelgeuse(self):
        from targeted_observing.bright_stars import load
        from targeted_observing.geometry import in_field
        if not (BASE/'cache/hip_bright_pm.tsv').exists():self.skipTest('Run review_data')
        star=next(s for s in load() if s['id']=='HIP 27989')
        self.assertTrue(in_field(star['ra_deg'],star['dec_deg'],89.349997,8.554764))
        self.assertLess(star['mag'],1)
