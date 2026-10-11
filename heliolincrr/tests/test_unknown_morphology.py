"""Regression checks for routing, provenance, and the production gate; no network."""
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import unknown_morphology as morphology
import export_unknown_ades as exporter
import digest2_pipeline as pipeline
import package_unknown_review as package


class MorphologyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.night = '20261002'
        self.paths = morphology.paths_for(self.night, self.root)
        self.rows = [{'linkage_id': i, 'trk_sub': f'TEST{i:03d}'} for i in range(3)]
        morphology.atomic_json(self.paths['raw'], self.rows)
        morphology.atomic_json(self.paths['passed'], self.rows[:1])
        morphology.atomic_json(self.paths['quarantine'], self.rows[1:])
        self.input = self.root / 'input'
        self.input.write_text('pixels')
        self.report = {
            'night': self.night, 'policy': morphology.POLICY,
            'engine_sha256': morphology.digest(Path(morphology.__file__).read_bytes()),
            'raw_sha256': morphology.digest(self.paths['raw'].read_bytes()),
            'inputs': [morphology.signature(self.input)],
            'counts': {'pass': 1, 'reject': 1, 'review': 1},
            'outputs': {k: {'sha256': morphology.digest(self.paths[k].read_bytes())}
                        for k in ('passed', 'quarantine')},
            'candidates': [{'row_sha256': morphology.row_hash(r), 'status': s}
                           for r, s in zip(self.rows, ['pass', 'reject', 'review'])],
        }
        morphology.atomic_json(self.paths['report'], self.report)

    def test_independent_majorities_keep_cb_evidence(self):
        self.assertEqual(morphology.decide([True, True, False], [True, False, True]), 'reject')

    def test_three_of_four_required(self):
        self.assertEqual(morphology.decide([True]*4, [True, True, False, False]), 'pass')
        self.assertEqual(morphology.decide([True]*4, [True, True, True, False]), 'reject')

    def test_missing_evidence_is_not_a_normal_measurement(self):
        self.assertEqual(morphology.decide([True, None, False], [True]*3), 'review')
        self.assertEqual(morphology.decide([True]*3, [True, None, False]), 'review')
        self.assertEqual(morphology.decide([True]*3, [False, None, False]), 'pass')
        self.assertEqual(morphology.decide([True, True], [True, True]), 'review')

    def test_fixed_and_broadened_psf_fit(self):
        rng = np.random.default_rng(42)
        template = np.exp(-morphology.RADIUS**2/8)
        template /= template.sum()
        for width in (1.0, 1.7):
            shape = morphology.map_coordinates(template,
                [morphology.Y/width+25, morphology.X/width+25], order=3,
                mode='constant', cval=0)/width**2
            stamp = 1000 + 10000*shape + rng.normal(0, 1, shape.shape)
            result = morphology.fit_stamp(stamp, 0, 0, 2, template)
            self.assertTrue(result['fit_success'])
            self.assertFalse(result['at_bound'])
            self.assertAlmostEqual(result['width_scale'], width, delta=.03)

    def test_passed_subset_and_strict_submission(self):
        rows, _ = morphology.approved_rows(self.night, self.root, self.rows)
        self.assertEqual(rows, self.rows[:1])
        with self.assertRaisesRegex(ValueError, 'rejected or unresolved'):
            morphology.approved_rows(self.night, self.root, self.rows, strict=True)
        with self.assertRaisesRegex(ValueError, 'does not match'):
            morphology.approved_rows(self.night, self.root, [{'linkage_id': 99}])

    def test_changed_inputs_block_gate(self):
        self.input.write_text('changed pixels')
        with self.assertRaisesRegex(ValueError, 'inputs changed'):
            morphology.read_current_report(self.night, self.root)

    def test_changed_catalog_and_output_block_gate(self):
        self.paths['passed'].write_text('[]')
        with self.assertRaisesRegex(ValueError, 'catalog changed'):
            morphology.read_current_report(self.night, self.root)
        morphology.atomic_json(self.paths['passed'], self.rows[:1])
        self.paths['raw'].write_text('[]')
        with self.assertRaisesRegex(ValueError, 'current raw catalog'):
            morphology.read_current_report(self.night, self.root)

    def test_incomplete_review_is_checked_after_automatic_filter(self):
        review = self.root / 'review.csv'
        review.write_text('tracklet_id,is_real\nTEST000,0\n')
        args = SimpleNamespace(catalog=str(self.paths['raw']), processed_root=str(self.root),
                               night=self.night, require_morphology=True, review_csv=str(review),
                               require_review=True, require_complete_review=True)
        observations, stats, rows = exporter.build_obs_rows(args, return_catalog=True)
        self.assertEqual((observations, rows), ([], []))
        self.assertEqual(stats['catalog_links'], 1)
        self.assertEqual(stats['review_missing_links'], 0)

    def test_production_process_cannot_score_or_send_quarantine(self):
        args = SimpleNamespace(night=self.night, processed_root=str(self.root))
        with patch.object(pipeline, 'score_snapshot') as score:
            with self.assertRaises(ValueError):
                pipeline.process_export(args, self.rows[1:], self.paths['raw'].parent)
            score.assert_not_called()

    def test_package_cannot_bypass_gate_with_raw_catalog(self):
        args = package.build_argparser().parse_args([
            self.night, '--processed-root', str(self.root), '--catalog', str(self.paths['raw'])])
        with self.assertRaises(ValueError):
            package.build_package(args)


if __name__ == '__main__':
    unittest.main()
