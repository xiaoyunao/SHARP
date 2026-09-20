"""Run with the isolated digest2 environment; all network transports are mocked."""
from __future__ import annotations

import argparse
import json
import multiprocessing
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import score_unknown_digest2 as scoring
import validate_neo_batches as validation


class FakeEngine:
    def __init__(self, **kwargs):
        assert kwargs["repeatable"] is True
        assert kwargs["no_threshold"] is False

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def classify_tracklet(self, observations, *, classes, is_ades):
        assert is_ades and classes == ["NEO"]
        # A deliberately conflicting Raw score detects use of the wrong column.
        score = {130: 64.9, 131: 65.0, 132: 65.1}[int(observations[0].ra)]
        return SimpleNamespace(noid=SimpleNamespace(NEO=score), raw=SimpleNamespace(NEO=100-score))


def concurrent_validate(plan, calls):
    def transport(batch):
        with open(calls, "a") as stream:
            stream.write(batch["id"] + "\n")
        return "Submission format valid."
    validation.validate(Path(plan), transport, simulation=True)


class RehearsalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.columns = "trkSub|mode|stn|obsTime|ra|dec|rmsRA|rmsDec|mag|band"
        self.header = "# version=2017\n# observatory\n! mpcCode 327\n" + self.columns + "\n"
        self.lines = [f"T{i}|CCD|327|2026-09-10T12:{minute:02d}:00.000Z|{130+i+j*.001}|10|0.2|0.2|20|G"
                      for i in range(3) for j, minute in enumerate((0, 10, 20))]
        self.psv = self.root / "source.psv"
        self.psv.write_text(self.header + "\n".join(self.lines) + "\n")
        self.catalog = self.root / "catalog.json"
        self.catalog.write_text(json.dumps([{"trk_sub": f"T{i}", "n_obs": 3} for i in range(3)]))
        self.review = self.root / "review.csv"
        self.review.write_text("tracklet_id,is_real\nT0,1\nT1,1\nT2,1\n")
        self.obscodes = self.root / "sites.txt"
        self.obscodes.write_text("327 117.5750 0.76278 +0.64470 Peking Observatory, Xinglong Station\n")
        self.model = self.root / "model.csv"
        self.model.write_text("fake model used only with mocked C engine\n")
        self.config = self.root / "MPC.config"
        self.config.write_text("obserr=1.0\n")
        self.args = argparse.Namespace(psv=self.psv, catalog=self.catalog, review_csv=self.review,
                                       obscodes=self.obscodes, model=self.model, config=self.config,
                                       out_dir=self.root / "output")

    def prepare(self):
        with patch("digest2.Digest2", FakeEngine):
            return scoring.prepare(self.args)

    def plan(self):
        self.prepare()
        return self.args.out_dir / "routing.json"

    def test_threshold_raw_score_and_complete_disjoint_batches(self):
        report = self.prepare()
        self.assertEqual([report["scores"][f"T{i}"]["route"] for i in range(3)],
                         ["unclassified", "neo candidate", "neo candidate"])
        batches = validation.checked_batches(self.args.out_dir / "routing.json")
        self.assertEqual(sum(b["n_obs"] for b in batches), 9)
        self.assertEqual({k for b in batches for k in b["trkSub"]}, {"T0", "T1", "T2"})

    def test_rejected_high_score_and_pending_are_not_exported(self):
        self.review.write_text("tracklet_id,is_real\nT0,\nT1,0\nT2,1\n")
        report = self.prepare()
        self.assertEqual(report["scores"]["T0"]["route"], "pending_review")
        self.assertEqual(report["scores"]["T1"]["route"], "rejected")
        self.assertEqual(report["batches"][0]["trkSub"], ["T2"])

    def test_no_review_scores_but_never_exports(self):
        self.args.review_csv = None
        report = self.prepare()
        self.assertEqual(len(report["scores"]), 3)
        self.assertEqual(report["batches"], [])

    def test_missing_entire_tracklet_and_missing_row_are_blocked(self):
        self.psv.write_text(self.header + "\n".join(self.lines[:2] + self.lines[6:]) + "\n")
        report = self.prepare()
        self.assertEqual(report["blocked"], ["T0", "T1"])
        self.assertEqual(report["batches"][0]["trkSub"], ["T2"])

    def test_engine_failure_never_becomes_low_score(self):
        with patch.object(FakeEngine, "classify_tracklet", side_effect=RuntimeError("engine failed")):
            report = self.prepare()
        self.assertEqual(report["blocked"], ["T0", "T1", "T2"])
        self.assertEqual(report["batches"], [])

    def test_nonfinite_engine_result_blocked(self):
        result = SimpleNamespace(noid=SimpleNamespace(NEO=float("nan")), raw=SimpleNamespace(NEO=60))
        with patch.object(FakeEngine, "classify_tracklet", return_value=result):
            self.assertEqual(len(self.prepare()["blocked"]), 3)

    def test_official_parser_dropping_rows_is_detected(self):
        with patch("digest2.observation.parse_ades_psv", return_value={}):
            report = self.prepare()
        self.assertEqual(report["blocked"], ["T0", "T1", "T2"])
        self.assertTrue(all("parser dropped" in s["error"] for s in report["scores"].values()))

    def test_observation_failures(self):
        _, _, groups = scoring.read_psv(self.psv.read_text())
        for field, value in (("ra", "nan"), ("ra", "360"), ("dec", "91"), ("rmsRA", "0"),
                             ("stn", "999"), ("obsTime", "2026-09-10T12:00:00"),
                             ("band", "X"), ("mode", "CMO"), ("provID", "2026 AA")):
            with self.subTest(field=field, value=value):
                rows = [dict(r) for r in groups["T0"]]
                rows[0][field] = value
                with self.assertRaises(ValueError):
                    scoring.check_rows(rows, 3, {"327"})

    def test_duplicate_epoch_and_long_arc_rejected(self):
        _, _, groups = scoring.read_psv(self.psv.read_text())
        rows = groups["T0"]
        rows[1]["obsTime"] = rows[0]["obsTime"]
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            scoring.check_rows(rows, 3, {"327"})
        rows[1]["obsTime"] = "2026-09-12T00:00:00Z"
        with self.assertRaisesRegex(ValueError, "single-night"):
            scoring.check_rows(rows, 3, {"327"})

    def test_parser_never_silently_drops_bad_columns(self):
        with self.assertRaisesRegex(ValueError, "column count"):
            scoring.read_psv(self.psv.read_text() + "Tbad|CCD\n")

    def test_duplicate_review_and_duplicate_catalog_rejected(self):
        with self.assertRaises(ValueError):
            scoring.read_decisions("tracklet_id,is_real\nT1,0\nT1,1\n")
        with self.assertRaises(ValueError):
            scoring.read_catalog('[{"trk_sub":"T1","n_obs":3},{"trk_sub":"T1","n_obs":3}]')
        with self.assertRaises(ValueError):
            scoring.read_catalog('[{"trk_sub":"T1","n_obs":2,"image_names":"a;a","objids":"1;1"}]')

    def test_empty_catalog_no_requests(self):
        self.psv.write_text(self.header)
        self.catalog.write_text("[]")
        plan = self.plan()
        with patch("subprocess.run") as run:
            state = validation.validate(plan, lambda b: self.fail("empty request"), simulation=True)
        self.assertEqual(state["batches"], {})
        run.assert_not_called()

    def test_existing_output_never_overwritten(self):
        self.prepare()
        with self.assertRaises(FileExistsError):
            self.prepare()

    def test_snapshot_mutation_blocks_all_requests(self):
        plan = self.plan()
        (plan.parent / "inputs/review.csv").write_text("changed")
        with self.assertRaisesRegex(ValueError, "snapshot"):
            validation.validate(plan, lambda b: self.fail("must not request"), simulation=True)

    def test_payload_mutation_blocks_all_requests(self):
        plan = self.plan()
        (plan.parent / "unclassified.psv").write_text("changed")
        with self.assertRaisesRegex(ValueError, "content"):
            validation.validate(plan, lambda b: self.fail("must not request"), simulation=True)

    def test_partial_failure_and_restart_never_resend_success_or_uncertain(self):
        plan = self.plan()
        calls = []
        def transport(batch):
            calls.append(batch["obj_type"])
            if batch["obj_type"] == "unclassified":
                raise TimeoutError("timeout after possible receive")
            return "Submission format valid."
        state = validation.validate(plan, transport, simulation=True)
        self.assertEqual({b["status"] for b in state["batches"].values()}, {"validated", "needs_reconciliation"})
        validation.validate(plan, transport, simulation=True)
        self.assertEqual(len(calls), 2)

    def test_interrupted_attempt_does_not_retry(self):
        plan = self.plan()
        batches = validation.checked_batches(plan)
        state = {"mode": "simulation", "batches": {b["id"]: {"status": "validating"} for b in batches}}
        validation.save_state(plan.parent / "simulation_state.json", state)
        result = validation.validate(plan, lambda b: self.fail("resend"), simulation=True)
        self.assertEqual({b["status"] for b in result["batches"].values()}, {"needs_reconciliation"})

    def test_http_success_or_curlid_does_not_imply_validation(self):
        plan = self.plan()
        result = validation.validate(plan, lambda b: "Submission ID is test. Error: bad format", simulation=True)
        self.assertEqual({b["status"] for b in result["batches"].values()}, {"rejected"})

    def test_concurrent_workers_send_each_batch_once(self):
        plan = self.plan()
        calls = self.root / "calls.txt"
        ctx = multiprocessing.get_context("spawn")
        workers = [ctx.Process(target=concurrent_validate, args=(str(plan), str(calls))) for _ in range(2)]
        for worker in workers:
            worker.start()
        for worker in workers:
            worker.join(15)
            self.assertEqual(worker.exitcode, 0)
        entries = calls.read_text().splitlines()
        self.assertEqual(len(entries), 2)
        self.assertEqual(len(set(entries)), 2)

    def test_transport_only_uses_test_url_and_exact_hashed_bytes(self):
        batch = {"obj_type": "neo candidate", "payload": b"exact payload\n"}
        with patch("subprocess.run", return_value=subprocess.CompletedProcess([], 0, b"Submission format valid.\nSMT_HTTP_STATUS:200\n", b"")) as run:
            validation.test_transport(batch, "observer@example.invalid")
        args, kwargs = run.call_args
        self.assertIn("https://minorplanetcenter.net/submit_psv_test", args[0])
        self.assertNotIn("https://minorplanetcenter.net/submit_psv", args[0])
        self.assertIn("source=<-", args[0])
        self.assertEqual(kwargs["input"], batch["payload"])

    def test_actual_mpc_receipt_is_neither_format_pass_nor_rejection(self):
        reply = "[SMT digest2 offline integration format test].  Submission ID is 2026-09-20T03:52:14.232_00000kOY\n"
        result = validation.classify_test_reply(reply)
        self.assertEqual(result["status"], "test_received")
        self.assertFalse(result["format_confirmed"])
        self.assertEqual(result["test_request_id"], "2026-09-20T03:52:14.232_00000kOY")

    def test_receipt_only_batch_is_not_resent(self):
        plan = self.plan()
        reply = "Submission ID is 2026-09-20T03:52:14.232_00000kOY\n"
        state = validation.validate(plan, lambda b: reply, simulation=True)
        self.assertEqual({b["status"] for b in state["batches"].values()}, {"test_received"})
        validation.validate(plan, lambda b: self.fail("must not resend receipt-only batch"), simulation=True)

    def test_actual_short_payload_http_400_is_definite_rejection(self):
        body = b"error: data size (117 bytes) not in allowed range (200 to 1000000000) bytes\n\nSMT_HTTP_STATUS:400\n"
        with patch("subprocess.run", return_value=subprocess.CompletedProcess([], 22, body, b"HTTP 400")):
            with self.assertRaisesRegex(validation.TestRequestRejected, "HTTP 400.*data size"):
                validation.test_transport({"obj_type": "unclassified", "payload": b"bad"}, "observer@example.invalid")

    def test_explicit_rejection_and_unknown_reply_have_distinct_states(self):
        plan = self.plan()
        def transport(batch):
            if batch["obj_type"] == "neo candidate":
                raise validation.TestRequestRejected("HTTP 400: invalid")
            return "unexpected service response"
        state = validation.validate(plan, transport, simulation=True)
        self.assertEqual({b["status"] for b in state["batches"].values()}, {"rejected", "needs_reconciliation"})

    def test_default_validator_cli_has_no_network(self):
        plan = self.plan()
        result = subprocess.run([sys.executable, str(Path(validation.__file__)), str(plan)],
                                text=True, capture_output=True, check=True)
        self.assertEqual(json.loads(result.stdout)["network_requests"], 0)
        self.assertFalse((plan.parent / "validation_state.json").exists())


if __name__ == "__main__":
    unittest.main()
