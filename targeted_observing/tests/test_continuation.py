import copy
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from targeted_observing.catalog import BASE
from targeted_observing.annual_selection import pick_route,DEFAULT_QUALITY
from targeted_observing.annual_tasks import status,weights
from targeted_observing.review_tools import export_plot_data as export
from targeted_observing.annual_twilight_queue import group_names


class ContinuationContracts(unittest.TestCase):
    def test_route_cost_changes_start_choice(self):
        c=json.loads((BASE/'defaults.json').read_text());c.update(far_slew_s=600,use_header_slew_floor=False)
        fields=[{'field_id':str(i),'ra_deg':r,'dec_deg':0,'core_targets':[str(i)]} for i,r in enumerate([0,90,180,270,10,11,12,13])]
        menu={'fields':fields};data={'valid':np.array([[1,0]]*4+[[0,1]]*4,bool),'starts':np.array([0.,1000.])}
        result=pick_route(menu,data,{str(i):1 for i in range(8)},4,{},0,None,DEFAULT_QUALITY,c,2)
        self.assertEqual(result[0],1)
        self.assertEqual(set(result[1]),{4,5,6,7})

    def test_science_survives_general_photometric_completion(self):
        task={'tier':'A: candidate','periods_h':[738.],'amplitude_mag':None}
        s={'times':[0.,30*86400.],'night_indices':[0,30],'errors':[.02,.02],'valid':[True,True],
           'zeropoints':[0,0],'quality':{'quality_lightcurve':True}}
        result=status(task,s,DEFAULT_QUALITY)
        self.assertFalse(result['sampling_ready']);self.assertTrue(result['photometric_complete'])
        self.assertIsNone(result['science_case_complete']);self.assertIsNone(result['period_solution_complete'])
        self.assertGreater(weights({'4902':task},{'4902':s},{'4902':1},DEFAULT_QUALITY,34)['4902'],0)
        self.assertGreater(result['required_baseline_days'],90)

    def test_low_amplitude_science_has_separate_precision_requirement(self):
        task={'periods_h':[19.39,26.1284],'amplitude_mag':.0925}
        result=status(task,None,DEFAULT_QUALITY)
        self.assertAlmostEqual(result['required_single_point_error_mag'],.0185)
        self.assertFalse(result['sampling_ready'])

    def test_mixed_twilight_group_order(self):
        m={'options':{'EV_A':[],'EV_B':[],'AM_A':[]}}
        self.assertEqual(group_names(m,'EV'),['EV_A','EV_B']);self.assertEqual(group_names(m,'AM'),['AM_A'])

    def test_export_missing_coordinates_and_nonfinite_values_are_errors(self):
        for obj in [{'dec_deg':0},{'ra_deg':float('nan'),'dec_deg':0},{'ra_deg':360,'dec_deg':0}]:
            with self.assertRaises(export.InputError):export.coordinates(obj,'fixture')

    def test_export_missing_calendar_day_is_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            for day in ['2026-11-01','2026-11-03']:(p/(day+'.json')).write_text('{}')
            args=SimpleNamespace(input=p,expected_nights=2)
            with self.assertRaisesRegex(export.InputError,'missing nights'):export.run(args)

    def test_export_sha_mismatch_is_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'.git').mkdir()
            with patch.object(export.subprocess,'run',return_value=SimpleNamespace(stdout='abc\n')):
                with self.assertRaisesRegex(export.InputError,'differs'):export.commit_metadata(p,'def')

    def test_export_requires_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(export.InputError,'source commit'):export.commit_metadata(Path(tmp),None)

    def test_right_censored_followups_preserve_false(self):
        from io import StringIO
        from targeted_observing.continuation_twilight_report import read_injections
        data=read_injections(StringIO('tracklet,recovered_plus1,recovered_plus3\nTrue,False,\nTrue,True,False\n'))
        self.assertEqual(float(data.recovered_plus1.mean()),.5)
        self.assertEqual(float(data.recovered_plus3.mean()),0.)
        self.assertTrue(data.recovered_plus3.isna().iloc[0])
