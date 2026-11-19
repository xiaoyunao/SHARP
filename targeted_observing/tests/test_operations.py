from contextlib import closing
from datetime import datetime
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from astropy.table import Table

from targeted_observing.server_common import load_config, runtime, write_json, TZ
from targeted_observing.server_quality import MARKER
from targeted_observing.server_history import connect, states, ingest
from targeted_observing.server_pipeline import process
from targeted_observing.science_manifest import build
from targeted_observing.server_publish import observing_times_text
from targeted_observing import server_service as service
from known_asteroid.update_all_matched_history import purge_excluded, iter_matched_files


class OperationsTests(unittest.TestCase):
    def config(self, folder):
        c=load_config();c.update(base=str(Path(folder)/'base'),processed_root=str(Path(folder)/'processed'))
        return c

    def test_excluded_night_cannot_recreate_l4_or_reenter_sampling(self):
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);root=Path(c['processed_root'])
            with closing(connect(c)) as db,db:
                db.execute('INSERT INTO exposures VALUES (?,?,?,?,?,?,?,?,?)',('bad','20261117','x',1,'r',30,80,20,'MP'))
                db.execute('INSERT INTO detections VALUES (?,?,?,?)',('bad','4722',.02,1))
                db.execute('INSERT INTO ingested VALUES (?,?)',('20261117','old'))
            write_json(root/'20261117'/MARKER,dict(reason='wrong-night duplicate data'))
            self.assertEqual(states(c),{})
            self.assertEqual(ingest(c,['20261117'])[0]['status'],'excluded_invalid_night')
            self.assertEqual(process('20261117',c)['status'],'excluded_invalid_night')
            self.assertFalse((root/'20261117/L4').exists())
            self.assertEqual(build('20261117',root,Path(folder)/'manifest.json')['accepted'],[])
            with closing(connect(c)) as db:
                self.assertEqual(db.execute('SELECT count(*) FROM ingested').fetchone()[0],0)

    def test_global_history_removes_invalid_night_even_without_new_input(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)/'processed';history=Path(folder)/'history';history.mkdir()
            path=history/'all.fits';state=history/'state.json'
            Table(dict(source_night=['20261117','20261118'],number=[1,2])).write(path)
            write_json(state,dict({'20261117':1,'20261118':1}))
            write_json(root/'20261117'/MARKER,dict(reason='invalid'))
            l4=root/'20261117/L4';l4.mkdir()
            Table(dict(number=[1])).write(l4/'20261117_matched_asteroids.fits')
            purge_excluded(root,path,state)
            self.assertEqual(list(Table.read(path)['source_night']),['20261118'])
            self.assertNotIn('20261117',json.loads(state.read_text()))
            self.assertEqual(iter_matched_files(root),[])

    def test_night_tick_does_not_even_enter_l2_poll_before_plan_end(self):
        class NightClock(datetime):
            @classmethod
            def now(cls,tz=None):return datetime(2026,11,19,1,tzinfo=TZ)
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);revision=Path(folder)/'revision'
            write_json(revision/'plan.json',dict(end_local='2026-11-19T06:08:05+08:00'))
            write_json(runtime(c)/'published/20261118.json',dict(revision=str(revision)))
            with patch.object(service,'datetime',NightClock),patch.object(service,'enqueue_processing') as enqueue:
                service.tick(c)
            enqueue.assert_not_called()

    def test_human_time_file_preserves_a_required_gap(self):
        def row(t):
            return dict(field_id='MP_0010',filter='r',exposure_s=30,mode='EV',candidate_id=None,
                        ra_deg=80.,dec_deg=20.,overhead_start_unix=t,start_unix=t+10,release_unix=t+58)
        plan=dict(night='20261118',start_local='2026-11-18T17:42:00+08:00',
                  end_local='2026-11-19T06:08:05+08:00',rows=[row(1794994920),row(1794995020)])
        text=observing_times_text(plan)
        self.assertEqual(text.count('开始本段'),2)
        self.assertEqual(text.count('MP_0010'),2)
        self.assertIn('没有 WAIT',text)

    def test_review_resume_rejects_inputs_changed_since_search(self):
        from targeted_observing import server_pipeline as pipeline
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);p=Path(folder)/'catalog.fits';p.write_bytes(b'searched catalog')
            stat=p.stat()
            write_json(runtime(c)/'manifests/20261118_science.json',
                       dict(input_signature='same',accepted=[dict(path=str(p),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns)]))
            write_json(runtime(c)/'processing/20261118.json',dict(known_complete=True,search_complete=True,input_signature='same'))
            p.write_bytes(b'new catalog with changed measurements')
            with patch.object(pipeline,'finish_outputs') as finish:
                with self.assertRaisesRegex(RuntimeError,'Input changed'):pipeline.resume_review('20261118',c)
            finish.assert_not_called()

    def test_reactivation_preserves_noon_and_has_one_reference_check(self):
        from targeted_observing import activate_server_v2 as activation
        class FixedClock(datetime):
            @classmethod
            def now(cls,tz=None):return datetime(2026,11,19,11,tzinfo=TZ)
        with tempfile.TemporaryDirectory() as folder:
            base=Path(folder);write_json(base/'v2/runtime/published/20261119.json',dict(verification=dict(passed=True)))
            write_json(base/'v2/deployment.json',{})
            (base/'data').mkdir();(base/'data/MPC_SUBMISSION_HOLD').touch()
            old='0 7 * * * independent-ep\n# BEGIN SHARP_REFERENCE_UPDATE\nold-reference-entry\n# END SHARP_REFERENCE_UPDATE\n'
            with patch.object(activation,'datetime',FixedClock),\
                 patch.object(activation.subprocess,'check_output',return_value=old),\
                 patch.object(activation.subprocess,'run') as install,patch('builtins.print'):
                activation.activate(base)
            cron=install.call_args.kwargs['input']
            self.assertIn('0 12 * * *',cron)
            self.assertIn('independent-ep',cron)
            self.assertEqual(cron.count('reference_update --enqueue'),1)
            self.assertNotIn('old-reference-entry',cron)


if __name__=='__main__':unittest.main()
