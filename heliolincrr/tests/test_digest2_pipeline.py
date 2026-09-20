"""Production integration gates; requests never leave this test process."""
import json
import multiprocessing
import os
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import test_digest2_rehearsal as rehearsal
import digest2_pipeline as pipeline
import score_unknown_digest2 as scoring
import watch_submit_reviews as watcher

RECEIPT = 'Submission ID is 2026-09-20T00:00:00.001_00000TEST'


def concurrent_submit(plan, ledger, calls):
    def transport(batch):
        with open(calls, 'a') as f:
            f.write(batch['id'] + '\n')
        return RECEIPT
    pipeline.submit_batches(Path(plan), Path(ledger), transport)


class ProductionTests(unittest.TestCase):
    def setUp(self):
        rehearsal.RehearsalTests.setUp(self)
        self.ledger = self.root / 'formal.json'

    prepare = rehearsal.RehearsalTests.prepare

    def plan(self):
        report = self.prepare()
        scoring.write_json(self.args.out_dir / 'schema_validation.json', {
            'status': 'local_validated', 'schema_sha256': scoring.sha256(pipeline.SCHEMA.read_bytes()),
            'batches': {b['id']: b['sha256'] for b in report['batches']}})
        return self.args.out_dir / 'routing.json'

    def valid_psv(self):
        header = ('# version=2017\n# observatory\n! mpcCode 327\n'
                  '# submitter\n! name SMT\n# observers\n! name SMT\n# measurers\n! name SMT\n'
                  '# telescope\n! design Schmidt\n! aperture 0.6\n! detector CCD\n')
        columns = self.columns.replace('|mag|band', '|astCat|mag|band|photCat')
        lines = [line.replace('|20|G', '|Gaia3E|20|G|Gaia3E') for line in self.lines]
        return (header + columns + '\n' + '\n'.join(lines) + '\n').encode()

    def test_schema_valid_and_corrupt_ra_header(self):
        data = self.valid_psv()
        pipeline.validate_schema(data)
        for bad in (data.replace(b'130.0', b'NOT_A_NUMBER'),
                    data.replace(b'# submitter\n! name SMT\n', b'')):
            with self.assertRaises(Exception):
                pipeline.validate_schema(bad)

    def test_success_rerun_no_duplicate(self):
        calls = []
        plan = self.plan()
        transport = lambda b: calls.append(b['id']) or RECEIPT
        pipeline.submit_batches(plan, self.ledger, transport)
        state = pipeline.submit_batches(plan, self.ledger, transport)
        self.assertEqual(len(calls), 2)
        self.assertTrue(all(b['status'] == 'received' for b in state['batches'].values()))

    def test_partial_failure_and_retry_does_not_resend(self):
        calls = []
        plan = self.plan()
        def transport(b):
            calls.append(b['id'])
            if len(calls) == 2:
                raise TimeoutError('response lost')
            return RECEIPT
        for _ in range(2):
            with self.assertRaisesRegex(RuntimeError, 'reconciliation'):
                pipeline.submit_batches(plan, self.ledger, transport)
        self.assertEqual(len(calls), 2)
        self.assertFalse(watcher.failure_is_retryable({"stderr_tail": "MPC attempt unresolved; manual reconciliation required"}))
        state = json.loads(self.ledger.read_text())
        self.assertEqual([b['status'] for b in state['batches'].values()].count('received'), 1)

    def test_interrupted_attempt_no_network(self):
        plan = self.plan()
        b = json.loads(plan.read_text())['batches'][0]
        scoring.write_json(self.ledger, {'mode': 'formal_mpc', 'batches': {b['id']: {**b, 'status': 'sending'}}})
        with self.assertRaisesRegex(RuntimeError, 'reconciliation'):
            pipeline.submit_batches(plan, self.ledger, lambda b: self.fail('must not send'))

    def test_changed_batch_preflight_blocks_before_any_send(self):
        plan = self.plan()
        b = json.loads(plan.read_text())['batches'][-1]
        scoring.write_json(self.ledger, {'mode': 'formal_mpc', 'batches': {'old-hash': {**b, 'status': 'received'}}})
        with self.assertRaisesRegex(RuntimeError, 'changed batch'):
            pipeline.submit_batches(plan, self.ledger, lambda b: self.fail('must not send'))

    def test_missing_or_stale_schema_proof_blocks(self):
        plan = self.plan()
        proof = plan.parent / 'schema_validation.json'
        proof.write_text('{}')
        with self.assertRaises(ValueError):
            pipeline.submit_batches(plan, self.ledger, lambda b: self.fail('must not send'))

    def test_nonreceipt_does_not_mark_received(self):
        plan = self.plan()
        with self.assertRaises(RuntimeError):
            pipeline.submit_batches(plan, self.ledger, lambda b: 'Submission format valid.')
        self.assertEqual(next(iter(json.loads(self.ledger.read_text())['batches'].values()))['status'], 'needs_reconciliation')

    def test_concurrent_formal_submission_exactly_once(self):
        plan = self.plan()
        calls = self.root / 'calls'
        ctx = multiprocessing.get_context('spawn')
        processes = [ctx.Process(target=concurrent_submit, args=(str(plan),str(self.ledger),str(calls))) for _ in range(2)]
        for p in processes:p.start()
        for p in processes:
            p.join(20)
            self.assertEqual(p.exitcode, 0)
        self.assertEqual(len(calls.read_text().splitlines()), 2)

    def test_historical_submitted_terminal_after_review_change(self):
        self.assertTrue(watcher.already_done({'status':'submitted','submit_csv':{'old':1}},
                                            {'new':2},{'new':2},True))

    def test_legacy_state_or_reply_blocks_direct_export(self):
        state = self.root / 'state.json'
        scoring.write_json(state, {'nights':{'20260910':{'status':'submitted'}}})
        with patch.dict(os.environ, {'SMT_REVIEW_SUBMIT_STATE':str(state)}):
            with self.assertRaisesRegex(RuntimeError,'Already submitted'):
                pipeline.guard_legacy('20260910',self.root)
            state.write_text('{}')
            (self.root/'20260910_unknown_mpc_reply.txt').write_text(RECEIPT)
            with self.assertRaisesRegex(RuntimeError,'legacy MPC reply'):
                pipeline.guard_legacy('20260910',self.root)

    def test_curl_exact_payload_no_redirect_or_retry(self):
        from types import SimpleNamespace
        batch = {'payload':b'exact bytes', 'obj_type':'neo candidate'}
        reply = SimpleNamespace(returncode=0, stdout=(RECEIPT+'\nSMT_HTTP_STATUS:200\n').encode(),stderr=b'')
        with patch.object(pipeline.subprocess,'run',return_value=reply) as run:
            self.assertEqual(pipeline.post_batch(batch,'a@b.c','test'),RECEIPT)
        args=run.call_args
        self.assertEqual(args.kwargs['input'],batch['payload'])
        self.assertIn('source=<-',args.args[0])
        self.assertNotIn('-L',args.args[0])
        self.assertEqual(args.kwargs['timeout'],100)


if __name__ == '__main__':
    unittest.main()
