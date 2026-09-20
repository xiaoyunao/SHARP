#!/usr/bin/env python3
"""Validate rehearsal batches ONLY via MPC's test endpoint, with a local ledger.

There is deliberately no production URL, submit option, or 'submitted' state.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from score_unknown_digest2 import json_hash, read_psv, sha256

TEST_URL = "https://minorplanetcenter.net/submit_psv_test"


class TestRequestRejected(RuntimeError):
    """A definite HTTP client-error response, distinct from a transport timeout."""


def classify_test_reply(reply: str) -> dict:
    """A test CurlID acknowledges receipt; even malformed RA can receive one."""
    valid = bool(re.search(r"\bSubmission format valid\.", reply, re.I))
    failure = bool(re.search(r"\b(?:invalid|error|rejected|not valid)\b", reply, re.I))
    receipt = re.search(
        r"\bSubmission ID is\s+(\d{4}-\d{2}-\d{2}T[\d:.]+_[A-Za-z0-9]+)\b", reply,
    )
    if failure:
        status = "rejected"
    elif valid:
        status = "validated"
    elif receipt:
        status = "test_received"
    else:
        status = "needs_reconciliation"
    result = {"status": status, "format_confirmed": status == "validated"}
    if receipt:
        result["test_request_id"] = receipt.group(1)
    return result


def save_state(path: Path, data: dict) -> None:
    tmp = path.with_suffix(".tmp")
    with tmp.open("w") as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    tmp.replace(path)
    fd = os.open(str(path.parent), os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def checked_batches(plan_path: Path) -> list[dict]:
    plan = json.loads(plan_path.read_text())
    if plan.get("mode") != "offline_rehearsal" or plan.get("production_enabled") is not False:
        raise ValueError("Expected an offline rehearsal plan")
    base = plan_path.parent.resolve()
    # Validate the frozen inputs as well as the final payloads before ANY request.
    provenance = plan["provenance"]
    snapshots = [("inputs/observations.psv", provenance["input_psv_sha256"]),
                 ("inputs/catalog.json", provenance["catalog_sha256"])]
    if provenance["review_sha256"]:
        snapshots.append(("inputs/review.csv", provenance["review_sha256"]))
    for name in ("model", "config", "obscodes"):
        snapshots.append((provenance[name]["snapshot"], provenance[name]["sha256"]))
    for relative, digest in snapshots:
        path = (base / relative).resolve()
        if not path.is_relative_to(base) or sha256(path.read_bytes()) != digest:
            raise ValueError(f"Changed or unsafe input snapshot: {relative}")
    seen, batches = set(), []
    for batch in plan["batches"]:
        obj_type = batch["obj_type"]
        if obj_type not in ("neo candidate", "unclassified"):
            raise ValueError("Unsupported object type")
        path = (base / batch["file"]).resolve()
        if not path.is_relative_to(base):
            raise ValueError("Batch path escapes rehearsal directory")
        payload = path.read_bytes()
        digest = sha256(payload)
        if digest != batch["sha256"] or batch["id"] != json_hash({"obj_type": obj_type, "psv_sha256": digest}):
            raise ValueError("Batch content or identity changed after scoring")
        _, _, groups = read_psv(payload.decode())
        if not groups or set(groups) != set(batch["trkSub"]) or seen.intersection(groups):
            raise ValueError("Empty, duplicate or inconsistent batch objects")
        if sum(map(len, groups.values())) != batch["n_obs"]:
            raise ValueError("Batch observation count mismatch")
        for key, rows in groups.items():
            score = plan["scores"][key]
            if score["status"] != "ok" or score["decision"] != 1 or score["route"] != obj_type:
                raise ValueError("Batch bypasses review/scoring gate")
            route = "neo candidate" if score["neo_noid"] >= 65 else "unclassified"
            if route != obj_type or json_hash(sorted(rows, key=lambda r: r["obsTime"])) != score["observations_sha256"]:
                raise ValueError("Payload differs from scored observations or threshold route")
        seen.update(groups)
        batches.append({**batch, "payload": payload})
    expected = {k for k, s in plan["scores"].items() if s["route"] in ("neo candidate", "unclassified")}
    if seen != expected:
        raise ValueError("Batch union differs from accepted, scored objects")
    return batches


def test_transport(batch: dict, email: str) -> str:
    # source=<- streams the already-hashed bytes, avoiding a read/send file race.
    proc = subprocess.run(
        ["curl", "-q", "--silent", "--show-error", "--fail-with-body", "--connect-timeout", "15",
         "--max-time", "90", "--write-out", "\nSMT_HTTP_STATUS:%{http_code}\n", TEST_URL,
         "--form-string", "ack=SMT digest2 offline integration format test",
         "--form-string", f"ac2={email}", "--form-string", f"obj_type={batch['obj_type']}",
         "-F", "source=<-"], input=batch["payload"], capture_output=True,
        timeout=100, check=False,
    )
    output = proc.stdout.decode(errors="replace")
    marker = re.search(r"\nSMT_HTTP_STATUS:(\d{3})\n$", output)
    http_status = int(marker.group(1)) if marker else 0
    reply = output[:marker.start()] if marker else output
    if 400 <= http_status < 500:
        raise TestRequestRejected(f"MPC test HTTP {http_status}: {reply}")
    if proc.returncode or not 200 <= http_status < 300:
        raise RuntimeError(f"Test request uncertain: curl={proc.returncode}; "
                           + proc.stderr.decode(errors="replace")[-2000:]
                           + reply[-2000:])
    return reply


def validate(plan_path: Path, transport, *, simulation: bool = False) -> dict:
    """Persist each result before starting another batch; serialize all workers."""
    base = plan_path.parent
    name = "simulation_state.json" if simulation else "validation_state.json"
    state_path = base / name
    with (base / (name + ".lock")).open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        batches = checked_batches(plan_path)
        state = json.loads(state_path.read_text()) if state_path.exists() else {
            "mode": "simulation" if simulation else "mpc_test_only", "batches": {}}
        expected_mode = "simulation" if simulation else "mpc_test_only"
        if state.get("mode") != expected_mode:
            raise ValueError("Ledger mode mismatch")
        for batch in batches:
            old = state["batches"].get(batch["id"], {})
            if old.get("status") in ("validated", "test_received", "needs_reconciliation", "rejected"):
                continue
            if old.get("status") == "validating":
                old["status"] = "needs_reconciliation"
                old["error"] = "Interrupted attempt; do not automatically resend"
                save_state(state_path, state)
                continue
            entry = {k: v for k, v in batch.items() if k != "payload"}
            entry.update(status="validating", started_utc=datetime.now(timezone.utc).isoformat())
            state["batches"][batch["id"]] = entry
            save_state(state_path, state)
            try:
                reply = transport(batch)
                entry["reply"] = reply
                entry.update(classify_test_reply(reply))
            except TestRequestRejected as exc:
                entry.update(status="rejected", format_confirmed=False, error=str(exc))
            except Exception as exc:
                entry.update(status="needs_reconciliation", error=str(exc))
            entry["finished_utc"] = datetime.now(timezone.utc).isoformat()
            save_state(state_path, state)
        return state


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path, help="routing.json from score_unknown_digest2.py")
    parser.add_argument("--test-endpoint", action="store_true", help="Actually call MPC test validator; default is local checks only")
    parser.add_argument("--email", default="", help="Acknowledgment address required only for test endpoint")
    args = parser.parse_args()
    if not args.test_endpoint:
        batches = checked_batches(args.plan)
        print(json.dumps({"mode": "local_checks_only", "batches": len(batches), "network_requests": 0}))
        return
    if not args.email or args.email.startswith("-") or "@" not in args.email:
        parser.error("--email is required for MPC test validation")
    state = validate(args.plan, lambda batch: test_transport(batch, args.email))
    print(json.dumps({"mode": state["mode"], "statuses": {k: v["status"] for k, v in state["batches"].items()}}, indent=2))
    if any(v["status"] != "validated" for v in state["batches"].values()):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
