import copy
import unittest
import numpy as np
from targeted_observing.annual_selection import pick,DEFAULT_QUALITY
from targeted_observing.annual_simulation import pick_legacy


class SelectionReview(unittest.TestCase):
    def fixture(self):
        return ({'fields':[{'core_targets':[k]} for k in ['A','A','B','C']]},
                {'valid':np.array([[1,0],[1,0],[0,1],[0,1]],bool)}, {'A':10,'B':6,'C':6})

    def test_review_counterexample_with_and_without_forecast(self):
        for forecast in [None,({'A':30,'B':30,'C':30},{'A':60,'B':60,'C':60},10)]:
            menu,data,values=self.fixture()
            self.assertEqual(pick_legacy(menu,data,values,2,{},0,forecast),(0,[0]))
            self.assertEqual(pick(menu,data,values,2,{},0,forecast),(1,[2,3]))

    def test_empty_and_no_available_fields(self):
        self.assertIsNone(pick({'fields':[]},{'valid':np.zeros((0,2),bool)},{},2,{},0))
        menu,data,values=self.fixture();data['valid'][:]=False
        self.assertIsNone(pick(menu,data,values,2,{},0))

    def test_single_field_and_completed_object(self):
        menu={'fields':[{'core_targets':['A']}]};data={'valid':np.ones((1,2),bool)}
        self.assertEqual(pick(menu,data,{'A':1},4,{},0),(0,[0]))
        self.assertIsNone(pick(menu,data,{'A':1},4,{'A':{'quality':{'quality_lightcurve':True}}},0))

    def test_equal_union_scores_rotate_only_among_ties(self):
        menu={'fields':[{'core_targets':[str(i)]} for i in range(5)]};data={'valid':np.eye(5,dtype=bool)}
        values={str(i):1 for i in range(5)}
        self.assertEqual(pick(menu,data,values,1,{},0)[0],0)
        self.assertNotEqual(pick(menu,data,values,1,{},1)[0],0)

    def test_nondefault_quality_forecast_and_baseline_wait(self):
        q=dict(DEFAULT_QUALITY,points=12,nights=2,baseline_days=7)
        menu={'fields':[{'core_targets':['A']},{'core_targets':['B']}]};data={'valid':np.ones((2,1),bool)}
        history={'A':{'times':[0],'night_indices':[0],'quality':{'good_points':12,'good_nights':2,'first_good_night':0}}}
        forecast=({'A':10,'B':10},{'A':14,'B':14},6)
        self.assertEqual(pick(menu,data,{'A':1,'B':.9},1,history,6,forecast,q)[1],[1])
        self.assertEqual(pick(menu,data,{'A':1,'B':.9},1,history,7,forecast,q)[1],[0])

    def test_overlapping_fields_do_not_hide_lower_ranked_distinct_target(self):
        menu={'fields':[{'core_targets':['A']} for _ in range(40)]+[{'core_targets':['B']}]}
        data={'valid':np.ones((41,1),bool)}
        self.assertEqual(pick(menu,data,{'A':10,'B':6},2,{},0),(0,[0,40]))

    def test_completed_prefix_does_not_affect_other_targets(self):
        menu,data,values=self.fixture()
        history={'A':{'quality':{'quality_lightcurve':True}}}
        self.assertEqual(pick(menu,data,values,2,history,0),(1,[2,3]))
