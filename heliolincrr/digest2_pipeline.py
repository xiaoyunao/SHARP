#!/usr/bin/env python3
"""Production unknown scoring and durable per-batch MPC submission.

The parent pipeline keeps its original Python environment. Only the worker uses
SMT_DIGEST2_PYTHON. Scoring and --validate are entirely local (no MPC requests).
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import re
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from score_unknown_digest2 import prepare, read_psv, sha256
from validate_neo_batches import checked_batches, classify_test_reply, save_state

RUNTIME = Path('/pipeline/xiaoyunao/data/heliolincrr/digest2_runtime')
SUBMIT_URL = 'https://minorplanetcenter.net/submit_psv'
SCHEMA = Path(__file__).with_name('ades_submit.xsd')


def require_submission_enabled() -> None:
    hold = Path(os.environ.get('MPC_SUBMISSION_HOLD_FILE',
                               '/pipeline/xiaoyunao/data/MPC_SUBMISSION_HOLD'))
    if hold.exists():
        raise RuntimeError(f'MPC submission paused: {hold}; no request sent')


def validate_schema(payload: bytes) -> None:
    """Lossless single-block conversion for our exporter, then ADES 2017 XSD."""
    from lxml import etree
    metadata, columns, groups = read_psv(payload.decode())
    root = etree.Element('ades', version='2017')
    block = etree.SubElement(root, 'obsBlock')
    context = etree.SubElement(block, 'obsContext')
    parent = None
    for line in metadata[1:]:
        parts = line[1:].strip().split(None, 1)
        name, value = parts[0], parts[1] if len(parts) == 2 else ''
        if line.startswith('#'):
            parent = etree.SubElement(context, name)
            if value:
                parent.text = value
        else:
            if parent is None:
                raise ValueError('ADES metadata child without parent')
            etree.SubElement(parent, name).text = value
    data = etree.SubElement(block, 'obsData')
    for rows in groups.values():
        for row in rows:
            optical = etree.SubElement(data, 'optical')
            for name in columns:
                if row[name]:
                    etree.SubElement(optical, name).text = row[name]
    schema = etree.XMLSchema(etree.parse(str(SCHEMA)))
    schema.assertValid(root)


def score_snapshot(psv: Path, catalog: list[dict], root: Path,
                   review_source: Path | None = None) -> tuple[Path, dict]:
    """Freeze caller inputs; verified accepted catalog is normalized to trkSub."""
    root.mkdir(parents=True, exist_ok=True)
    job = Path(tempfile.mkdtemp(prefix='run_', dir=root))
    (job / 'source.psv').write_bytes(psv.read_bytes())
    (job / 'catalog.json').write_text(json.dumps(catalog, allow_nan=False))
    cmd = [os.environ.get('SMT_DIGEST2_PYTHON', str(RUNTIME / 'env/bin/python')),
           str(Path(__file__).resolve()), '--worker', '--psv', str(job / 'source.psv'),
           '--catalog', str(job / 'catalog.json'), '--out-dir', str(job / 'score'),
           '--obscodes', os.environ.get('SMT_DIGEST2_OBSCODES', str(RUNTIME / 'ObsCodes.html'))]
    if review_source is not None:
        original = review_source.read_bytes()
        (job / 'source_review.csv').write_bytes(original)
        with (job / 'accepted.csv').open('w') as stream:
            writer = csv.writer(stream)
            writer.writerow(['tracklet_id', 'is_real'])
            writer.writerows([(r.get('trk_sub') or r['trkSub'], 1) for r in catalog])
        cmd.extend(['--review-csv', str(job / 'accepted.csv')])
    proc = subprocess.run(cmd, text=True, capture_output=True, timeout=900)
    if proc.returncode:
        raise RuntimeError(f'digest2/schema gate failed; audit={job}\n{proc.stdout}\n{proc.stderr}')
    plan_path = job / 'score/routing.json'
    report = json.loads(plan_path.read_text())
    report['source_review_sha256'] = sha256(original) if review_source is not None else None
    save_state(job / 'audit.json', report)
    return plan_path, report


def guard_legacy(night: str, l4: Path) -> None:
    """A historical receipt is terminal even if review/manifest files changed."""
    ledger = l4 / f'{night}_unknown_mpc_batches.json'
    old_state = Path(os.environ.get('SMT_REVIEW_SUBMIT_STATE',
                         '/pipeline/xiaoyunao/data/heliolincrr/review_submit_state.json'))
    if old_state.exists():
        status = json.loads(old_state.read_text()).get('nights', {}).get(night, {}).get('status')
        if status == 'submitted':
            raise RuntimeError(f'Already submitted night {night}; historical resubmission blocked')
    reply = l4 / f'{night}_unknown_mpc_reply.txt'
    if not ledger.exists() and reply.exists() and reply.read_text().strip():
        raise RuntimeError(f'Existing legacy MPC reply: {reply}; reconcile before any submission')


def post_batch(batch: dict, email: str, ack: str) -> str:
    require_submission_enabled()
    proc = subprocess.run(
        ['curl', '-q', '--silent', '--show-error', '--fail-with-body',
         '--connect-timeout', '15', '--max-time', '90',
         '--write-out', '\nSMT_HTTP_STATUS:%{http_code}\n', SUBMIT_URL,
         '--form-string', f'ack={ack}', '--form-string', f'ac2={email}',
         '--form-string', f"obj_type={batch['obj_type']}", '-F', 'source=<-'],
        input=batch['payload'], capture_output=True, timeout=100)
    output = proc.stdout.decode(errors='replace')
    marker = re.search(r'\nSMT_HTTP_STATUS:(\d{3})\n$', output)
    code = int(marker.group(1)) if marker else 0
    reply = output[:marker.start()] if marker else output
    if proc.returncode or not 200 <= code < 300:
        raise RuntimeError(f'MPC request requires reconciliation: HTTP={code} curl={proc.returncode}\n'
                           + reply + '\n' + proc.stderr.decode(errors='replace'))
    return reply


def submit_batches(plan_path: Path, ledger_path: Path, transport) -> dict:
    """Persist before I/O. An ambiguous attempt is never automatically retried."""
    require_submission_enabled()
    with Path(str(ledger_path) + '.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        batches = checked_batches(plan_path)
        proof = json.loads((plan_path.parent / 'schema_validation.json').read_text())
        if proof.get('status') != 'local_validated' or proof.get('schema_sha256') != sha256(SCHEMA.read_bytes()):
            raise ValueError('Missing/current ADES schema validation required')
        if proof.get('batches') != {b['id']: b['sha256'] for b in batches}:
            raise ValueError('Schema validation does not match submission batches')
        if not batches or json.loads(plan_path.read_text())['blocked']:
            raise ValueError('Empty or blocked submission plan')
        state = json.loads(ledger_path.read_text()) if ledger_path.exists() else {
            'mode': 'formal_mpc', 'batches': {}}
        if state.get('mode') != 'formal_mpc':
            raise ValueError('Submission ledger mode mismatch')
        # Preflight ALL overlaps before posting either batch. Review changes cannot
        # repackage an already sent object under a new payload/batch hash.
        for batch in batches:
            for key, old in state['batches'].items():
                if set(batch['trkSub']).intersection(old['trkSub']) and key != batch['id']:
                    raise RuntimeError('Previously attempted objects changed batch; manual reconciliation required')
        for batch in batches:
            old = state['batches'].get(batch['id'])
            if old:
                if old['status'] == 'received':
                    continue
                raise RuntimeError('MPC attempt unresolved; manual reconciliation required; do not resend')
            require_submission_enabled()
            entry = {k: v for k, v in batch.items() if k != 'payload'}
            entry.update(status='sending', plan=str(plan_path),
                         started_utc=datetime.now(timezone.utc).isoformat())
            state['batches'][batch['id']] = entry
            save_state(ledger_path, state)
            try:
                reply = transport(batch)
                entry['reply'] = reply
                parsed = classify_test_reply(reply)
                if parsed['status'] != 'test_received':
                    raise RuntimeError('Formal endpoint did not return an unambiguous Submission ID')
                entry.update(status='received', submission_id=parsed['test_request_id'])
            except Exception as exc:
                entry.update(status='needs_reconciliation', error=str(exc))
                save_state(ledger_path, state)
                raise RuntimeError('MPC attempt unresolved; manual reconciliation required') from exc
            entry['finished_utc'] = datetime.now(timezone.utc).isoformat()
            save_state(ledger_path, state)
        return state


def process_export(args, kept_rows: list[dict], l4: Path) -> None:
    plan, report = score_snapshot(Path(args.out), kept_rows, l4 / 'digest2', Path(args.review_csv))
    print('[DIGEST2]', plan)
    print('[LOCAL-VALIDATED] ADES schema, complete observations, human acceptance and digest2 routes')
    if args.submit:
        ledger = l4 / f'{args.night}_unknown_mpc_batches.json'
        result = submit_batches(plan, ledger, lambda batch: post_batch(batch, args.ac2_email, args.ack))
        print('[MPC-RECEIVED]', json.dumps(result, ensure_ascii=False))
        if args.response_out:
            save_state(Path(args.response_out), result)
    elif args.response_out:
        save_state(Path(args.response_out), {'status': 'local_validated', 'plan': str(plan),
                                            'network_requests': 0})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true', required=True)
    parser.add_argument('--psv', type=Path, required=True)
    parser.add_argument('--catalog', type=Path, required=True)
    parser.add_argument('--review-csv', type=Path)
    parser.add_argument('--obscodes', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    args.model = args.config = None
    validate_schema(args.psv.read_bytes())
    report = prepare(args)
    if any(s['status'] != 'ok' for s in report['scores'].values()):
        raise RuntimeError('Unknown scoring failed; inspect routing.json')
    for batch in checked_batches(args.out_dir / 'routing.json'):
        validate_schema(batch['payload'])
    save_state(args.out_dir / 'schema_validation.json', {
        'status': 'local_validated', 'schema_sha256': sha256(SCHEMA.read_bytes()),
        'payload_sha256': report['provenance']['input_psv_sha256'],
        'batches': {b['id']: b['sha256'] for b in report['batches']}})


if __name__ == '__main__':
    main()
