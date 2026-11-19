from datetime import datetime, timedelta
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from targeted_observing.server_common import TZ, load_config, runtime, write_json
from targeted_observing.server_l2_watch import readiness, recent_nights, science_stem
from targeted_observing import server_service as service


class L2WatchTests(unittest.TestCase):
    def inputs(self, folder):
        config = load_config()
        config.update(base=str(Path(folder)/'base'), processed_root=str(Path(folder)/'processed'),
                      raw_root=str(Path(folder)/'raw'))
        revision = Path(folder)/'revision'
        write_json(revision/'plan.json', dict(end_local='2026-11-19T06:08:05+08:00'))
        write_json(runtime(config)/'published/20261118.json', dict(revision=str(revision)))
        raw = Path(config['raw_root'])/'20261118/OBJ/OBJ_YSO_ORI_r_0001.fits'
        l2 = Path(config['processed_root'])/'20261118/L2/OBJ_YSO_ORI_r_0001_cat.fits.gz'
        raw.parent.mkdir(parents=True);l2.parent.mkdir(parents=True)
        raw.write_bytes(b'raw');l2.write_bytes(b'catalog')
        now = datetime.now(TZ)+timedelta(seconds=121)
        return config, raw, l2, now

    def test_quiet_midnight_gap_cannot_start_processing(self):
        with tempfile.TemporaryDirectory() as folder:
            c,raw,l2,now = self.inputs(folder)
            value = readiness('20261118',c,datetime(2026,11,19,1,tzinfo=TZ))
            self.assertEqual(value['status'],'waiting_for_observing_end')

    def test_science_names_and_missing_catalog_block_the_whole_night(self):
        with tempfile.TemporaryDirectory() as folder:
            c,raw,l2,now = self.inputs(folder)
            self.assertTrue(readiness('20261118',c,now)['ready'])
            missing=raw.with_name('OBJ_OtherProgram_i_0002.fits');missing.write_bytes(b'raw')
            (raw.parent/'OBJ_start.fits').write_bytes(b'non-science')
            value=readiness('20261118',c,now)
            self.assertEqual(value['missing_l2'],[missing.name])
            self.assertFalse(value['ready'])
            self.assertEqual(science_stem(l2),'OBJ_YSO_ORI_r_0001')
            self.assertIsNone(science_stem(raw.with_name('OBJ_focus_0001.fits')))

    def test_replaced_raw_requires_new_l2_and_settling_after_last_write(self):
        with tempfile.TemporaryDirectory() as folder:
            c,raw,l2,now = self.inputs(folder)
            raw.write_bytes(b'replaced raw')
            self.assertEqual(readiness('20261118',c,now)['stale_l2'],[raw.name])
            l2.write_bytes(b'new catalog')
            self.assertEqual(readiness('20261118',c)['status'],'waiting_for_stable_input')
            self.assertTrue(readiness('20261118',c,now)['ready'])

    def test_watch_only_recent_previous_nights_from_deployment(self):
        c=load_config()
        self.assertEqual(recent_nights(c,datetime(2026,11,19,1,tzinfo=TZ)),['20261118'])
        self.assertEqual(recent_nights(c,datetime(2026,11,22,1,tzinfo=TZ)),['20261121','20261120','20261119'])

    def test_morning_tick_enqueues_without_running_scheduler(self):
        class Morning(datetime):
            @classmethod
            def now(cls,tz=None):return datetime(2026,11,19,7,tzinfo=TZ)
        with tempfile.TemporaryDirectory() as folder:
            c,raw,l2,now=self.inputs(folder)
            with patch.object(service,'datetime',Morning),patch.object(service,'enqueue_processing') as enqueue,\
                 patch.object(service,'daily') as daily,patch.object(service,'plan') as plan:
                service.tick(c)
            enqueue.assert_called_once_with('20261118',c)
            daily.assert_not_called();plan.assert_not_called()

    def test_noon_entry_only_schedules_and_does_not_launch_processing(self):
        with tempfile.TemporaryDirectory() as folder:
            c,raw,l2,now=self.inputs(folder)
            write_json(runtime(c)/'published/20261119.json',dict(revision='scheduled'))
            with patch.object(service,'enqueue_processing') as enqueue:
                self.assertEqual(service.daily(c,'20261119'),dict(plan=dict(revision='scheduled')))
            enqueue.assert_not_called()

    def test_queue_is_idempotent_and_requeues_late_catalogs(self):
        with tempfile.TemporaryDirectory() as folder:
            c,raw,l2,now=self.inputs(folder)
            with patch('targeted_observing.server_l2_watch.readiness',return_value=dict(ready=True)),\
                 patch.object(service.subprocess,'check_output',side_effect=['','123;cluster','123']):
                value=service.enqueue_processing('20261118',c)
                self.assertEqual(value['job_id'],'123')
                self.assertEqual(service.enqueue_processing('20261118',c),value)
                write_json(runtime(c)/'processing/20261118.json',dict(status='complete'))
                self.assertEqual(service.enqueue_processing('20261118',c),value)
                l2.with_name('OBJ_EP_r_0002_cat.fits').write_bytes(b'late catalog')
                self.assertEqual(service.enqueue_processing('20261118',c)['status'],'already_running')


if __name__=='__main__':
    unittest.main()
