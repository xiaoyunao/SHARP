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
import unknown_review_gate as review_gate


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

    def test_abnormal_catalog_rejects_without_psf(self):
        self.assertEqual(morphology.decide([True, True, False], []), 'reject')
        self.assertEqual(morphology.decide([True]*3, [False]*3), 'reject')

    def test_normal_catalog_requires_psf_pass(self):
        self.assertEqual(morphology.decide([False]*3, [False]*3), 'pass')
        self.assertEqual(morphology.decide([False]*3, [True]*3), 'reject')
        self.assertEqual(morphology.decide([False]*3, []), 'review')
        self.assertEqual(morphology.decide([False]*3, [None]*3), 'review')

    def test_three_of_four_required(self):
        self.assertEqual(morphology.decide([False]*4, [True, True, False, False]), 'pass')
        self.assertEqual(morphology.decide([False]*4, [True, True, True, False]), 'reject')
        self.assertEqual(morphology.decide([True, True, False, False], [False]*4), 'pass')
        self.assertEqual(morphology.decide([True, True, True, False], []), 'reject')

    def test_missing_evidence_is_not_a_normal_measurement(self):
        self.assertEqual(morphology.decide([True, None, False], [True]*3), 'review')
        self.assertEqual(morphology.decide([False]*3, [True, None, False]), 'review')
        self.assertEqual(morphology.decide([False]*3, [False, None, False]), 'pass')
        self.assertEqual(morphology.decide([True, True], [True, True]), 'review')

    def test_screen_routes_only_catalog_normal_candidates_to_psf(self):
        images = [f'frame{i}.fits.gz' for i in range(3)]
        rows = [{'linkage_id': oid, 'trk_sub': f'TEST{oid:03d}', 'n_obs': 3,
                 'image_names': ';'.join(images), 'objids': ';'.join([str(oid)]*3),
                 'mjds': '61000;61000.01;61000.02'} for oid in range(1, 6)]
        morphology.atomic_json(self.paths['raw'], rows)
        calls = []

        def measured_frames(function, tasks, workers):
            calls.append(function.__name__)
            frames = []
            for _, _, name, ids in tasks:
                if function is morphology.catalog_frame:
                    measurements = {str(oid): {'catalog_bad': {1: True, 2: False, 3: False,
                                                               4: None, 5: False}[oid]} for oid in ids}
                else:
                    self.assertEqual(ids, [2, 3, 5])
                    measurements = {str(oid): {'psf_bad': {2: True, 3: False, 5: None}[oid]}
                                    for oid in ids}
                frames.append({'image_name': name, 'observations': measurements, 'inputs': []})
            return frames

        with patch.object(morphology, 'run_frames', side_effect=measured_frames):
            report = morphology.screen_night(self.night, self.root, force=True)
        self.assertEqual(calls, ['catalog_frame', 'psf_frame'])
        candidates = report['candidates']
        self.assertEqual([r['status'] for r in candidates], ['reject', 'reject', 'pass', 'review', 'review'])
        self.assertEqual([r['reason'] for r in candidates], [
            'catalog_majority_abnormal', 'psf_majority_abnormal', 'both_stages_normal',
            'catalog_evidence_insufficient', 'psf_evidence_insufficient'])
        self.assertTrue(all('psf_bad' not in o for r in (candidates[0], candidates[3])
                            for o in r['observations']))
        self.assertEqual(json.loads(self.paths['passed'].read_text()), [rows[2]])
        self.assertEqual(morphology.read_current_report(self.night, self.root)['counts'],
                         {'pass': 1, 'reject': 2, 'review': 2})

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
        review.write_text('tracklet_id,is_real\nTEST000,0\nTEST002,0\n')
        args = SimpleNamespace(catalog=str(self.paths['raw']), processed_root=str(self.root),
                               night=self.night, require_morphology=True, review_csv=str(review),
                               require_review=True, require_complete_review=True)
        observations, stats, rows = exporter.build_obs_rows(args, return_catalog=True)
        self.assertEqual((observations, rows), ([], []))
        self.assertEqual(stats['catalog_links'], 2)
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

    def test_unresolved_goes_to_review_but_reject_cannot_override(self):
        rows, _ = review_gate.reviewable_rows(self.night, self.root, self.rows)
        self.assertEqual(rows, [self.rows[0], self.rows[2]])
        with self.assertRaisesRegex(ValueError, 'rejected'):
            review_gate.reviewable_rows(self.night, self.root, self.rows, strict=True)
        review = self.root / 'review.csv'
        review.write_text('tracklet_id,is_real\nTEST000,1\nTEST001,1\n')
        with self.assertRaisesRegex(ValueError, 'human acceptance'):
            review_gate.require_human_acceptance(rows, review)
        review.write_text('tracklet_id,is_real\nTEST000,1\nTEST002,1\n')
        review_gate.require_human_acceptance(rows, review)

    def test_missing_unresolved_decision_blocks_complete_export(self):
        review = self.root / 'review.csv'
        review.write_text('tracklet_id,is_real\nTEST000,0\n')
        args = SimpleNamespace(catalog=str(self.paths['raw']), processed_root=str(self.root),
                               night=self.night, require_morphology=True, review_csv=str(review),
                               require_review=True, require_complete_review=True)
        with self.assertRaises(ValueError):
            exporter.build_obs_rows(args, return_catalog=True)

    def test_reviewed_unresolved_can_reach_local_validation(self):
        review = self.root / 'review.csv'
        review.write_text('tracklet_id,is_real\nTEST002,1\n')
        args = SimpleNamespace(night=self.night, processed_root=str(self.root),
                               review_csv=str(review), out=str(self.root / 'out.psv'),
                               submit=False, response_out='')
        with patch.object(pipeline, 'score_snapshot', return_value=(self.root / 'plan', {})) as score:
            pipeline.process_export(args, self.rows[2:], self.paths['raw'].parent)
            score.assert_called_once()
        review.write_text('tracklet_id,is_real\nTEST002,0\n')
        with patch.object(pipeline, 'score_snapshot') as score:
            with self.assertRaisesRegex(ValueError, 'human acceptance'):
                pipeline.process_export(args, self.rows[2:], self.paths['raw'].parent)
            score.assert_not_called()

    def test_package_includes_unresolved_and_preserves_decisions(self):
        out = self.root / 'package'
        out.mkdir()
        (out / 'gifs').mkdir()
        (out / 'gifs/TEST001_link0001_20261002.gif').write_bytes(b'old gif')
        review = out / f'{self.night}_unknown_review.csv'
        review.write_text('tracklet_id,is_real\nTEST000,0\nTEST001,1\nTEST002,1\n')
        args = package.build_argparser().parse_args([
            self.night, '--processed-root', str(self.root), '--out-dir', str(out),
            '--plot-root', str(self.root / 'plots')])
        with patch.object(package, 'write_review_ades_fits', return_value=('', 0, {'scores': {}})), \
             patch.object(package, 'write_review_detail_fits', return_value=('', 0)):
            result = package.build_package(args)
        manifest = json.loads(Path(result['manifest_json']).read_text())
        self.assertEqual([r['morphology_status'] for r in manifest['rows']], ['pass', 'review'])
        self.assertEqual(review.read_text(), 'tracklet_id,is_real\nTEST000,0\nTEST002,1\n')
        self.assertEqual(manifest['n_catalog_rows'], 2)
        self.assertFalse(list((out / 'gifs').glob('*.gif')))
        self.assertEqual(len(list(self.root.glob(f'{self.night}_excluded_gifs_*/*.gif'))), 1)


if __name__ == '__main__':
    unittest.main()
