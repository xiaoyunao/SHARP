from datetime import datetime, timedelta
from contextlib import closing
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

from targeted_observing.server_common import load_config, TZ, runtime, write_json, followup_cutoff
from targeted_observing.server_followup import validate, submit, queue_snapshot
from targeted_observing.server_history import sampling, ingest, states
from targeted_observing.science_manifest import build, load


def catalog(path,start='2026-11-18T11:00:00',seconds=120):
    primary=fits.PrimaryHDU();primary.header['EXPSTA']=start;primary.header['EXPTIME']=seconds
    primary.header['FILTER']='r'
    table=Table(dict(objID=[1],RA_Win=[80.],DEC_Win=[20.],MagErr_PSF=[.02],Flag=[0]))
    ext=fits.BinTableHDU(table)
    for key,value in dict(CTYPE1='RA---TAN',CTYPE2='DEC--TAN',CRVAL1=80.,CRVAL2=20.,
                          CRPIX1=4608.,CRPIX2=4616.,CD1_1=-.00032,CD2_2=.00032).items():ext.header[key]=value
    path.parent.mkdir(parents=True,exist_ok=True)
    fits.HDUList([primary,ext]).writeto(path,overwrite=True)


class V2Contracts(unittest.TestCase):
    def config(self,folder):
        value=load_config();value['base']=str(Path(folder)/'base');value['processed_root']=str(Path(folder)/'processed')
        return value

    def request(self,key='candidateA',night='20261118'):
        return dict(candidate_id=key,night=night,field_id='2630',start_time='2026-11-18T23:30:00+08:00',frames=5,exposure_s=30,filter='r')

    def test_request_contract_requires_timezone_and_four_text_digits(self):
        self.assertEqual(validate(self.request())['execution_name'],'FOLLOW_2630')
        for updates in ({'field_id':2630},{'start_time':'2026-11-18T23:00:00'},
                        {'frames':True},{'filter':'W'},{'night':'20261117'},{'candidate_id':'../escape'}):
            with self.assertRaises((ValueError,TypeError)):validate({**self.request(),**updates})

    def test_every_candidate_creates_event_and_whole_night_queue_is_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder)
            for hour,key in ((11,'A'),(13,'B'),(14,'C')):
                result=submit(self.request(key),c,datetime(2026,11,18,hour,tzinfo=TZ))
                self.assertEqual(result['scheduling'],'noon run' if hour==11 else 'replan queued')
            queue,errors=queue_snapshot(c,'20261118')
            self.assertEqual({r['candidate_id'] for r in queue},{'A','B','C'});self.assertFalse(errors)
            self.assertEqual(len(list((runtime(c)/'events/20261118').glob('*.json'))),3)
            self.assertEqual(queue_snapshot(c,'20261119'),([],[]))
            with self.assertRaises(ValueError):submit(self.request('expired'),c,datetime(2026,11,19,13,tzinfo=TZ))

    def test_cutoff_routes_late_commands_to_next_night_preserving_clock_time(self):
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder)
            ontime=submit(self.request('early'),c,datetime(2026,11,18,16,59,59,tzinfo=TZ))
            self.assertEqual(ontime['night'],'20261118');self.assertFalse(ontime['deferred'])
            for key,hour in (('boundary',17),('late',23)):
                result=submit(self.request(key),c,datetime(2026,11,18,hour,tzinfo=TZ))
                self.assertEqual(result['night'],'20261119');self.assertTrue(result['deferred'])
                self.assertEqual(result['start_time'],'2026-11-19T23:30:00+08:00')
                self.assertEqual(result['scheduling'],'noon run')
            overnight=submit(self.request('overnight'),c,datetime(2026,11,19,1,tzinfo=TZ))
            self.assertEqual(overnight['night'],'20261119')
            self.assertEqual({r['candidate_id'] for r in queue_snapshot(c,'20261118')[0]},{'early'})

    def test_earlier_execution_start_closes_acceptance_early(self):
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);revision=runtime(c)/'plans/probe'
            write_json(revision/'plan.json',dict(start_local='2026-11-18T16:30:00+08:00',window_start_local='2026-11-18T16:30:00+08:00'))
            write_json(runtime(c)/'published/20261118.json',dict(revision=str(revision)))
            self.assertEqual(followup_cutoff(c,'20261118'),datetime(2026,11,18,16,30,tzinfo=TZ))
            value=submit(self.request(),c,datetime(2026,11,18,16,30,tzinfo=TZ))
            self.assertEqual(value['night'],'20261119')

    def test_shorter_ev_reservation_does_not_close_before_actual_execution(self):
        from targeted_observing import server_service as service
        class PrestartClock(datetime):
            @classmethod
            def now(cls,tz=None):return datetime(2026,11,18,16,45,tzinfo=TZ)
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);revision=runtime(c)/'plans/probe'
            write_json(revision/'plan.json',dict(start_local='2026-11-18T17:30:00+08:00',window_start_local='2026-11-18T16:30:00+08:00'))
            write_json(runtime(c)/'published/20261118.json',dict(revision=str(revision)))
            self.assertEqual(followup_cutoff(c,'20261118'),datetime(2026,11,18,17,tzinfo=TZ))
            self.assertFalse(submit(self.request(),c,PrestartClock.now())['deferred'])
            with patch.object(service,'datetime',PrestartClock),\
                 patch('targeted_observing.server_scheduler.build',return_value=dict(start_local='2026-11-18T16:55:00+08:00')) as build_mock,\
                 patch('targeted_observing.server_publish.publish',return_value={'revision':'test'}):
                service.plan('20261118',c,'prestart_test')
            self.assertGreater(build_mock.call_args.kwargs['now'],PrestartClock.now().timestamp())
            self.assertNotIn('previous_plan',build_mock.call_args.kwargs)

    def test_frozen_plan_never_invokes_scheduler_and_tick_defers_pending(self):
        from targeted_observing import server_service as service
        class FrozenClock(datetime):
            @classmethod
            def now(cls,tz=None):return datetime(2026,11,18,19,tzinfo=TZ)
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);revision=runtime(c)/'plans/probe'
            write_json(revision/'plan.json',dict(start_local='2026-11-18T17:42:00+08:00'))
            write_json(runtime(c)/'published/20261118.json',dict(revision=str(revision),followups=[]))
            submit(self.request(),c,datetime(2026,11,18,16,59,59,tzinfo=TZ))
            with patch.object(service,'datetime',FrozenClock),patch('targeted_observing.server_scheduler.build') as build_mock:
                for finalize in (False,True):
                    with self.assertRaisesRegex(ValueError,'frozen'):service.plan('20261118',c,'test',finalize=finalize)
                with patch.object(service,'enqueue_processing'):
                    service.tick(c)
                build_mock.assert_not_called()
            self.assertEqual(queue_snapshot(c,'20261119')[0][0]['candidate_id'],'candidateA')
            results=list((runtime(c)/'events/20261118').glob('*.result.json'))
            self.assertEqual(len(results),1)
            self.assertEqual(json.loads(results[0].read_text())['status'],'deferred_before_execution_freeze')

    def test_direct_sftp_late_request_is_moved_to_tomorrow_queue(self):
        from targeted_observing import server_service as service
        class ArrivalClock(datetime):
            @classmethod
            def now(cls,tz=None):return datetime(2026,11,18,17,1,tzinfo=TZ)
            @classmethod
            def fromtimestamp(cls,value,tz=None):return cls.now(tz)
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder)
            source=runtime(c)/'followups/20261118/candidateA.json'
            write_json(source,self.request())
            with patch.object(service,'datetime',ArrivalClock):
                self.assertEqual(service.discover_events('20261118',c),[])
            self.assertFalse(source.exists())
            self.assertEqual(queue_snapshot(c,'20261119')[0][0]['start_time'],'2026-11-19T23:30:00+08:00')
            self.assertTrue((runtime(c)/'deferred_followups/20261118/candidateA.json').exists())

    def test_slow_render_cannot_replace_public_list_after_start(self):
        from targeted_observing import server_publish as publisher
        class FinishedTooLate(datetime):
            @classmethod
            def now(cls,tz=None):return datetime(2026,11,18,17,42,tzinfo=TZ)
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);public=Path(c['base'])/'script/current_plan.txt'
            public.parent.mkdir(parents=True);public.write_text('original final list')
            value=dict(night='20261118',start_local='2026-11-18T17:42:00+08:00',rows=[])
            with patch.object(publisher,'plots',return_value={}),patch.object(publisher,'datetime',FinishedTooLate):
                with self.assertRaisesRegex(ValueError,'rendering'):publisher.publish(value,c,'late_render_test')
            self.assertEqual(public.read_text(),'original final list')
            self.assertFalse((runtime(c)/'published/20261118.json').exists())

    def test_plot_history_matches_real_centers_across_ra_zero_not_old_names(self):
        from targeted_observing.server_history import connect
        from targeted_observing.server_publish import plot_history
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder)
            with closing(connect(c)) as db, db:
                for uid,ra,name in [('one',359.999,'MP_legacy_0010'),('two',0.04,'OtherProgram'),('far',90.,'MP_2630')]:
                    db.execute('INSERT INTO exposures VALUES (?,?,?,?,?,?,?,?,?)',(uid,'20261117',uid,1.,'r',30.,ra,20.,name))
            fields=Table(dict(field_id=['2630','YSO_ORI'],center_ra=[.001,80.],center_dec=[20.,20.]))
            history=plot_history(c,fields)
            self.assertEqual(list(history['exposure_count']),[2,0])
            self.assertEqual(history.meta['unmatched_exposures'],1)

    def test_sync_script_retention_keeps_ten_evenings_current_and_unrelated_files(self):
        from targeted_observing.script_retention import cleanup
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);root=Path(c['base'])/'script';root.mkdir(parents=True)
            for day in range(18,30):
                for suffix in ('plan.txt','plan.json','START_TIMES.json','segment01.txt'):
                    (root/f'202611{day}_{suffix}').write_text('list')
            for name in ('current_plan.txt','operator_notes.txt','20261001_unknown.txt'):
                (root/name).write_text('keep')
            cleanup(c,datetime(2026,11,28,12,tzinfo=TZ))
            self.assertFalse((root/'20261118_plan.txt').exists())
            self.assertTrue((root/'20261119_plan.txt').exists())
            self.assertTrue((root/'20261128_plan.txt').exists())
            self.assertTrue((root/'current_plan.txt').exists())
            self.assertTrue((root/'operator_notes.txt').exists())
            self.assertTrue((root/'20261001_unknown.txt').exists())

    def test_initial_script_cleanup_starts_today_not_older_retention_window(self):
        from targeted_observing.script_retention import cleanup
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);root=Path(c['base'])/'script';root.mkdir(parents=True)
            for night in ('20261117','20261118'):(root/(night+'_plan.txt')).write_text('list')
            cleanup(c,datetime(2026,11,18,14,tzinfo=TZ))
            self.assertFalse((root/'20261117_plan.txt').exists())
            self.assertTrue((root/'20261118_plan.txt').exists())

    def test_unavailable_twilight_segment_is_explicitly_skipped_without_masking_errors(self):
        from targeted_observing import offline_twilight as twilight
        config=twilight.load_config()
        empty=lambda sky,fields,bounds,side,c:dict(side=side,rows=[],fields=[],reason='No full-hour route found')
        failed=lambda *a:dict(passed=False,errors=['no full-hour schedule'],samples=0)
        with patch.object(twilight,'setup_iers',return_value={}),patch.object(twilight,'make_sky'),\
             patch.object(twilight,'night_windows',return_value={'EV':[0,3600],'AM':[10000,13600],'dark_crossings':[3600,10000]}),\
             patch.object(twilight,'plan_window',side_effect=empty),patch.object(twilight,'verify_window',side_effect=failed) as check:
            value=twilight.make_plan('2026-11-02',config,grid=([],{}),allow_empty=True)
            self.assertTrue(all(b['verification']['skipped'] for b in value['blocks']))
            with self.assertRaises(ValueError):twilight.make_plan('2026-11-02',config,grid=([],{}))
            check.side_effect=lambda *a:dict(passed=False,errors=['footprint altitude'],samples=1)
            with self.assertRaises(ValueError):twilight.make_plan('2026-11-02',config,grid=([],{}),allow_empty=True)

    def test_shorter_twilight_fallback_keeps_cadence_and_minimum_visits(self):
        from targeted_observing import offline_twilight as twilight
        config=twilight.load_config();checked=[]
        def plan_window(sky,fields,bounds,side,c):
            return dict(side=side,rows=[] if c['fill_window'] else [{'visit':i} for i in range(3)],fields=[],
                        reason='test shorter route',window_unix=bounds,reserved_window_unix=bounds)
        def check(sky,block,c):
            checked.append(c)
            return dict(passed=True,errors=[],samples=3)
        with patch.object(twilight,'setup_iers',return_value={}),patch.object(twilight,'make_sky'),\
             patch.object(twilight,'night_windows',return_value={'EV':[0,3600],'AM':[10000,13600],'dark_crossings':[3600,10000]}),\
             patch.object(twilight,'plan_window',side_effect=plan_window),patch.object(twilight,'verify_window',side_effect=check):
            value=twilight.make_plan('2026-11-02',config,grid=([],{}),allow_empty=True,shorter_fallback=True)
        self.assertTrue(all(b['shorter_fallback'] for b in value['blocks']))
        for c in checked:
            self.assertFalse(c['fill_window']);self.assertEqual(c['minimum_visits'],3)
            self.assertEqual(c['cadence_min_s'],config['cadence_min_s'])
            self.assertEqual(c['min_footprint_alt_deg'],config['min_footprint_alt_deg'])

    def test_twilight_selection_rejects_hour_angle_before_final_validation(self):
        from targeted_observing.offline_twilight import WindowGeometry,load_config
        g=WindowGeometry.__new__(WindowGeometry);g.bounds=[0,3600];g.times=np.array([0,1800,3600])
        g.fields=[dict(radius_deg=2.1)];g.c=load_config();shape=(1,3)
        g.m=dict(alt=np.full(shape,50.),moon_sep=np.full(shape,120.),elong=np.full(shape,50.),
                 moon_alt=np.full(3,-30.),moon_fraction=np.zeros(3),sun_alt=np.full(3,-20.))
        g.ha=np.full(shape,4.8);self.assertTrue(g.valid(0,0,3600))
        g.ha=np.full(shape,5.1);self.assertFalse(g.valid(0,0,3600))

    def test_extra_short_cadence_frames_preserve_an_existing_valid_sequence(self):
        from targeted_observing.server_scheduler import cadence_visits
        self.assertEqual(cadence_visits([0,600,1200,1248,1288.5]),3)
        self.assertEqual(cadence_visits([0,48,96,144]),1)
        self.assertEqual(cadence_visits([0,720.0000001,1440]),3)

    def test_science_inputs_include_yso_and_follow_but_reject_stale_and_focus(self):
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);root=Path(c['processed_root']);l2=root/'20261118/L2'
            for i,name in enumerate(('OBJ_YSO_ORI_0001_cat.fits','OBJ_FOLLOW_2630_0002_cat.fits','OBJ_EP_0003_cat.fits','OtherProgram_0004_cat.fits')):
                catalog(l2/name,start=f'2026-11-18T11:{i:02d}:00')
            catalog(l2/'OBJ_MP_2630_0004_cat.fits',start='2026-11-17T11:00:00')
            catalog(l2/'OBJ_focus-5.88_0005_cat.fits')
            out=Path(folder)/'manifest.json';value=build('20261118',root,out)
            self.assertEqual(len(value['accepted']),4);self.assertEqual(len(value['rejected']),2)
            self.assertAlmostEqual(value['accepted'][0]['mid_unix'],Time('2026-11-18T11:03:00').unix,places=4)
            path=Path(value['accepted'][0]['path']);path.touch()
            with self.assertRaises(RuntimeError):load(out)

    def test_burst_counts_do_not_replace_night_span_and_baseline(self):
        c=load_config();rows=[];epoch=datetime(2026,11,1,tzinfo=TZ).timestamp()
        for i in range(200):rows.append(dict(mid=epoch+i,night='20261101',band='r',error=.02,good=True))
        self.assertEqual(sampling(rows,c)['good_nights'],0)
        self.assertFalse(sampling(rows,c)['sampling_ready'])
        rows=[]
        for day in range(0,30,3):
            for i in range(12):rows.append(dict(mid=epoch+day*86400+i*300,night=(datetime(2026,11,1)+timedelta(days=day)).strftime('%Y%m%d'),band='r',error=.02,good=True))
        self.assertTrue(sampling(rows,c)['sampling_ready'])
        for r in rows:r['band']='g'
        self.assertEqual(sampling(rows,c)['good_points'],0)

    def test_history_is_actual_deduplicated_and_replaces_reprocessed_measurements(self):
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);root=Path(c['processed_root']);file='OBJ_TJNL5_2630_0001_cat.fits'
            catalog(root/'20261118/L2'/file)
            l4=root/'20261118/L4';l4.mkdir()
            path=l4/'20261118_matched_asteroids.fits'
            Table(dict(source_file=[file],number=[4722],MagErr_PSF=[.02],Flag=[0])).write(path)
            ingest(c,['20261118']);ingest(c,['20261118'])
            self.assertEqual(states(c)['4722']['good_points'],1)
            Table(dict(source_file=[file],number=[4722],MagErr_PSF=[.2],Flag=[0])).write(path,overwrite=True)
            ingest(c,['20261118'])
            self.assertEqual(states(c)['4722']['good_points'],0)

    def test_unknown_export_uses_actual_midpoint_and_physical_filter(self):
        from heliolincrr.export_unknown_ades import load_detection_row
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'20261118/L2/OBJ_YSO_ORI_0001_cat.fits'
            catalog(path)
            with fits.open(path,mode='update') as hdus:
                table=Table(hdus[1].data)
                for key,value in dict(MJD=1.,RAErr_Win=.00001,DECErr_Win=.00001,Mag_Aper4=18.,MagErr_Aper4=.02).items():
                    table[key]=[value]
                h=hdus[1].header.copy();hdus[1]=fits.BinTableHDU(table,header=h)
                hdus[0].header['FILTER']='g'
            row=load_detection_row(root,'20261118','OBJ_YSO_ORI_0001.fits',1)
            self.assertAlmostEqual(row['MJD'],Time('2026-11-18T11:01:00').mjd,places=9)
            self.assertEqual(row['physical_filter'],'g')

    def test_known_history_reprocessing_replaces_same_count_rows(self):
        from known_asteroid.update_all_matched_history import append_one
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);combined=root/'all.fits';state=root/'state.json';source=root/'source.fits'
            Table(dict(epoch=[1.],number=[4722])).write(source)
            value={};append_one(combined,state,'20261118',source,value)
            Table(dict(epoch=[2.],number=[4722])).write(source,overwrite=True)
            append_one(combined,state,'20261118',source,value)
            table=Table.read(combined);self.assertEqual(len(table),1);self.assertEqual(table['epoch'][0],2.)

    def test_invalid_night_is_not_empty_known_complete_and_keeps_old_evidence(self):
        from targeted_observing.server_pipeline import process
        with tempfile.TemporaryDirectory() as folder:
            c=self.config(folder);l4=Path(c['processed_root'])/'20261118/L4';l4.mkdir(parents=True)
            (l4/'20261118_known_asteroid_status.json').write_text('{"known_complete":true}')
            raw=l4/'20261118_unknown_links.json';raw.write_text('[{"old_evidence":true}]')
            value=process('20261118',c)
            self.assertFalse(value['known_complete']);self.assertFalse(value['unknown_complete'])
            self.assertFalse(json.loads((l4/'20261118_known_asteroid_status.json').read_text())['known_complete'])
            self.assertEqual(json.loads(raw.read_text()),[{'old_evidence':True}])


if __name__=='__main__':unittest.main()
