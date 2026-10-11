"""Incident regressions without production data, Slurm, or network calls."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]


class DailyIncidentTests(unittest.TestCase):
    def test_hold_runs_known_extraction_without_submit_and_does_not_rerun(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wrapper = root / 'run_daily.sh'
            wrapper.write_bytes((ROOT / 'known_asteroid/run_daily.sh').read_bytes())
            runner = root / 'submit_pipeline_slurm.sh'
            runner.write_text('#!/bin/bash\nprintf "%s\\n" "$@" >> "$(dirname "$0")/calls"\n')
            runner.chmod(0o755)
            l2 = root / 'processed/20261007/L2'
            l2.mkdir(parents=True)
            (l2 / 'OBJ_MP_001_cat.fits').touch()
            hold = root / 'HOLD'
            hold.touch()
            env = dict(os.environ, PROCESSED_ROOT=str(root / 'processed'), TARGET_NIGHT='20261007',
                       WORKSPACE=str(root / 'runtime'), MPC_SUBMISSION_HOLD_FILE=str(hold), PYTHON=sys.executable)
            cmd = ['bash', str(wrapper), '2026-10-08']
            first = subprocess.run(cmd, env=env, text=True, capture_output=True)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            calls = (root / 'calls').read_text()
            self.assertIn('--submit-mpc\nfalse\n', calls)
            l4 = l2.parent / 'L4'
            l4.mkdir()
            (l4 / '20261007_known_asteroid_status.json').write_text(json.dumps({
                'known_complete': True, 'matched_asteroids': {'rows': 123}}))
            second = subprocess.run(cmd, env=env, text=True, capture_output=True)
            self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
            self.assertEqual((root / 'calls').read_text(), calls)

    def test_failed_night_stops_its_steps_but_not_remaining_nights(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = (ROOT / 'heliolincrr/run_daily_pipeline.sh').read_text()
            source = source.removesuffix('main 2>&1 | tee -a "${LOG_PATH}"\n')
            self.assertNotIn('main 2>&1', source)
            source += '''
night_for_offset() { printf '202610%02d\\n' "$((11 - $1))"; }
run_processing_for_night() {
  echo "started=$1"
  if [[ "$1" == "20261006" ]]; then
    false
    echo "UNSAFE_CONTINUATION"
  fi
  echo "finished=$1"
}
main 2>&1 | tee -a "${LOG_PATH}"
'''
            env = dict(os.environ, LOG_DIR=tmp, SURVEY_DIR=tmp, RECOVERY_LOOKBACK_DAYS='7')
            result = subprocess.run(['bash'], input=source, text=True, capture_output=True, env=env)
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertNotIn('UNSAFE_CONTINUATION', result.stdout)
            self.assertNotIn('finished=20261006', result.stdout)
            self.assertIn('finished=20261007', result.stdout)
            self.assertIn('finished=20261010', result.stdout)
            self.assertIn('incomplete: 1 failed nights', result.stdout)

    def test_known_submission_hold_blocks_transport(self):
        spec = importlib.util.spec_from_file_location('known_export_incident', ROOT / 'known_asteroid/export_ades.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            hold = Path(tmp) / 'HOLD'
            hold.write_text('327 suspended')
            with patch.dict(os.environ, {'MPC_SUBMISSION_HOLD_FILE':str(hold)}), patch.object(module.subprocess, 'run') as run:
                with self.assertRaisesRegex(RuntimeError, 'submission paused'):
                    module.run_curl(module.SUBMIT_URL, SimpleNamespace())
                run.assert_not_called()


if __name__ == '__main__':
    unittest.main()
