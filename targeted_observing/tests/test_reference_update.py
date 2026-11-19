from datetime import datetime
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from targeted_observing import reference_update as reference
from targeted_observing.server_common import load_config, runtime, write_json, TZ


class ReferenceUpdateTests(unittest.TestCase):
    def config(self, folder):
        c=load_config();c['base']=str(Path(folder)/'base')
        write_json(Path(c['base'])/'v2/config.json',c)
        p=Path(c['base'])/'known_asteroid/astorb.dat';p.parent.mkdir();p.write_bytes(b'old orbit data')
        package=Path(folder)/'package';(package/'cache').mkdir(parents=True)
        return c,package,p

    def download(self, urls, path):
        path.write_bytes(gzip.compress(b'new orbit data') if path.name.endswith('.gz') else b'new IERS')
        return dict(url=urls[0],sha256=reference.sha(path))

    def test_validated_update_is_monthly_and_keeps_one_previous_version(self):
        with tempfile.TemporaryDirectory() as folder:
            c,package,orb=self.config(folder)
            with patch.object(reference,'__file__',str(package/'reference_update.py')),\
                 patch.object(reference,'download',side_effect=self.download),\
                 patch.object(reference,'validate_astorb',return_value=dict(objects=1_600_000)),\
                 patch.object(reference,'validate_iers',return_value=dict(last_mjd=62000)),\
                 patch.object(reference.subprocess,'check_output',return_value=''),\
                 patch('astropy.utils.data.import_file_to_cache'):
                self.assertEqual(reference.update(c)['status'],'complete')
                self.assertEqual(reference.update(c)['status'],'already_updated_this_month')
            self.assertEqual(orb.read_bytes(),b'new orbit data')
            self.assertEqual((runtime(c)/'reference_updates/previous/astorb.dat').read_bytes(),b'old orbit data')
            self.assertEqual(json.loads((Path(c['base'])/'v2/config.json').read_text())['iers_file'],'iers_current.all')
            self.assertFalse(list((runtime(c)/'reference_updates').glob('.staging_*')))

    def test_validation_failure_preserves_production_and_does_not_mark_month_done(self):
        with tempfile.TemporaryDirectory() as folder:
            c,package,orb=self.config(folder)
            with patch.object(reference,'download',side_effect=self.download),\
                 patch.object(reference,'validate_astorb',side_effect=ValueError('invalid orbital file')):
                with self.assertRaisesRegex(ValueError,'invalid orbital file'):reference.update(c)
            self.assertEqual(orb.read_bytes(),b'old orbit data')
            self.assertTrue(reference.due(c))
            self.assertEqual(json.loads((runtime(c)/'reference_updates/status.json').read_text())['status'],'retry_required')

    def test_queue_depends_on_active_processing_and_deduplicates_itself(self):
        with tempfile.TemporaryDirectory() as folder:
            c,package,orb=self.config(folder)
            with patch.object(reference.subprocess,'check_output',side_effect=['123|sharp_v2_20261118\n','124;cluster']) as run:
                self.assertEqual(reference.enqueue(c)['after_processing_jobs'],['123'])
                self.assertIn('afterany:123',run.call_args.args[0])
            with patch.object(reference.subprocess,'check_output',return_value='124|sharp_ref_update\n'):
                self.assertEqual(reference.enqueue(c)['status'],'already_queued')


if __name__=='__main__':unittest.main()
