#!/usr/bin/env python3
"""Replay saved unknown candidates without touching reviews or MPC ledgers.

Received submissions are identified by formal batch receipts, or legacy receipt
+ matching exported PSV/catalog snapshots. No MPC endpoint is called.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from score_unknown_digest2 import read_psv
from unknown_morphology import atomic_json, digest, paths_for, row_hash, screen_night


def submission_evidence(night, root):
    l4 = root / night / 'L4'
    ledger = l4 / f'{night}_unknown_mpc_batches.json'
    evidence = {}
    if ledger.exists():
        data = json.loads(ledger.read_text())
        if data.get('mode') != 'formal_mpc':
            raise ValueError(f'{night}: unexpected ledger mode')
        for batch in data['batches'].values():
            if batch['status'] != 'received':
                if batch['status'] in ('sending', 'needs_reconciliation'):
                    raise ValueError(f'{night}: unresolved submission attempt')
                continue
            plan = Path(batch['plan'])
            payload = (plan.parent / batch['file']).read_bytes()
            if digest(payload) != batch['sha256']:
                raise ValueError(f'{night}: submitted payload hash mismatch')
            _, _, groups = read_psv(payload.decode())
            if set(groups) != set(batch['trkSub']) or sum(map(len, groups.values())) != batch['n_obs']:
                raise ValueError(f'{night}: submitted identities/counts mismatch')
            frozen = json.loads((plan.parent.parent / 'catalog.json').read_text())
            snapshots = {r.get('trk_sub') or r.get('trkSub'): r for r in frozen}
            for trk in groups:
                if trk in evidence:
                    raise ValueError(f'{night}: duplicate received trkSub')
                evidence[trk] = {'submission_id': batch['submission_id'], 'obj_type': batch['obj_type'],
                                 'source': 'formal_batch_and_frozen_catalog', 'row': snapshots[trk]}
        return evidence
    reply = l4 / f'{night}_unknown_mpc_reply.txt'
    if not reply.exists():
        return evidence
    receipts = re.findall(r'Submission ID is\s+([0-9T:.\-_A-Za-z]+)', reply.read_text())
    if not receipts:
        raise ValueError(f'{night}: legacy reply has no recognizable receipt')
    payload = l4 / f'{night}_unknown_links_submit_ades.psv'
    _, _, groups = read_psv(payload.read_text())
    snapshots = json.loads((l4 / f'{night}_unknown_links_submit_masked.json').read_text())
    by_id = {r.get('trk_sub') or r.get('trkSub'): r for r in snapshots}
    if set(groups) != set(by_id):
        raise ValueError(f'{night}: legacy PSV/catalog identities differ')
    for trk, observations in groups.items():
        if len(observations) != int(by_id[trk]['n_obs']):
            raise ValueError(f'{night}: legacy observation count differs')
        evidence[trk] = {'submission_id': ';'.join(receipts), 'obj_type': 'legacy_unknown',
                         'source': 'legacy_receipt_and_export_snapshot', 'row': by_id[trk]}
    return evidence


def protected_files(nights, root, review_root, state_path, hold_path):
    files = [p for p in [state_path, hold_path] if p.exists()]
    for night in nights:
        l4 = root / night / 'L4'
        for path in l4.glob(f'{night}_unknown*'):
            if (path.name.endswith(('_links.json', '_links.fits')) or
                '_submit_' in path.name or '_mpc_' in path.name) and path.suffix != '.lock':
                files.append(path)
        review = review_root / night
        files.extend(review.glob('*.csv'))
        files.extend(review.glob('*manifest.json'))
    return {str(p): digest(p.read_bytes()) for p in sorted(set(files)) if p.is_file()}


def csv_output(path, rows, fields):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--start', default='20260830')
    ap.add_argument('--end', required=True)
    ap.add_argument('--processed-root', type=Path, default=Path('/processed1'))
    ap.add_argument('--review-root', type=Path, default=Path('/pipeline/xiaoyunao/heliolincrr/review_packages'))
    ap.add_argument('--state', type=Path, default=Path('/pipeline/xiaoyunao/data/heliolincrr/review_submit_state.json'))
    ap.add_argument('--hold', type=Path, default=Path('/pipeline/xiaoyunao/data/MPC_SUBMISSION_HOLD'))
    ap.add_argument('--out-dir', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()
    if not (len(args.start) == len(args.end) == 8 and args.start.isdigit() and args.end.isdigit()
            and args.start <= args.end and 1 <= args.workers <= 16):
        ap.error('Invalid date range or worker count')
    args.out_dir.mkdir(parents=True, exist_ok=True)
    existing = sorted(p.name for p in args.processed_root.iterdir()
                      if p.is_dir() and p.name.isdigit() and args.start <= p.name <= args.end)
    nights = [n for n in existing if paths_for(n, args.processed_root)['raw'].exists()]
    before_path = args.out_dir / 'protected_files_before.json'
    before = protected_files(nights, args.processed_root, args.review_root, args.state, args.hold)
    if before_path.exists():
        original = json.loads(before_path.read_text())
        if original != before:
            raise RuntimeError('Protected inputs changed since the original replay; reconcile before resuming')
    else:
        atomic_json(before_path, before)
    submitted = {n: submission_evidence(n, args.processed_root) for n in nights}
    status = {'start': args.start, 'end': args.end, 'started_utc': datetime.now(timezone.utc).isoformat(),
              'status': 'running', 'nights': {}, 'no_saved_catalog': [n for n in existing if n not in nights]}
    atomic_json(args.out_dir / 'status.json', status)
    all_candidates, nightly, errors = [], [], []
    for night in nights:
        print(f'[NIGHT] {night}', flush=True)
        try:
            report = screen_night(night, args.processed_root, args.workers)
            raw = json.loads(paths_for(night, args.processed_root)['raw'].read_text())
            raw_by_id = {r.get('trk_sub') or r.get('trkSub'): r for r in raw}
            for trk, evidence in submitted[night].items():
                if trk not in raw_by_id or row_hash(evidence['row']) != row_hash(raw_by_id[trk]):
                    raise ValueError(f'Submitted snapshot differs from raw candidate {trk}')
            received_counts = Counter()
            for candidate in report['candidates']:
                trk = candidate['trk_sub']
                evidence = submitted[night].get(trk)
                item = {k: candidate[k] for k in ['linkage_id', 'trk_sub', 'status', 'reason',
                                                 'catalog_bad_count', 'psf_bad_count']}
                item.update(night=night, n_obs=len(candidate['observations']), submitted=bool(evidence),
                            submission_id=evidence['submission_id'] if evidence else '',
                            obj_type=evidence['obj_type'] if evidence else '',
                            submission_evidence=evidence['source'] if evidence else '')
                all_candidates.append(item)
                if evidence:
                    received_counts[candidate['status']] += 1
            record = {'night': night, 'unknown': report['n_candidates'], **report['counts'],
                      'submitted': len(submitted[night]),
                      **{'submitted_' + k: received_counts[k] for k in ('pass', 'reject', 'review')}}
            nightly.append(record)
            status['nights'][night] = {'status': 'complete', **record}
            print(json.dumps(record), flush=True)
        except Exception as exc:
            error = {'night': night, 'error': f'{type(exc).__name__}: {exc}'}
            errors.append(error)
            status['nights'][night] = {'status': 'failed', **error}
            print(json.dumps(error), flush=True)
        atomic_json(args.out_dir / 'status.json', status)
    after = protected_files(nights, args.processed_root, args.review_root, args.state, args.hold)
    changed = [p for p in set(before) | set(after) if before.get(p) != after.get(p)]
    columns = ['night', 'trk_sub', 'linkage_id', 'status', 'reason', 'n_obs',
               'catalog_bad_count', 'psf_bad_count', 'submitted', 'submission_id', 'obj_type',
               'submission_evidence']
    csv_output(args.out_dir / 'all_candidates.csv', all_candidates, columns)
    csv_output(args.out_dir / 'submitted_candidates.csv', [r for r in all_candidates if r['submitted']], columns)
    csv_output(args.out_dir / 'nightly_summary.csv', nightly,
               ['night', 'unknown', 'pass', 'reject', 'review', 'submitted',
                'submitted_pass', 'submitted_reject', 'submitted_review'])
    summary = {'start': args.start, 'end': args.end, 'nights': len(nights),
               'completed_nights': len(nightly), 'errors': errors,
               'unknown': len(all_candidates), 'counts': dict(Counter(r['status'] for r in all_candidates)),
               'submitted': sum(r['submitted'] for r in all_candidates),
               'submitted_counts': dict(Counter(r['status'] for r in all_candidates if r['submitted'])),
               'no_saved_catalog': status['no_saved_catalog'],
               'protected_file_count': len(before), 'protected_files_changed': changed,
               'interpretation': 'reject means flagged by the adopted morphology rule, not independently confirmed false; pass does not establish reality.',
               'nightly': nightly}
    atomic_json(args.out_dir / 'summary.json', summary)
    status['status'] = 'failed' if errors or changed else 'complete'
    status['finished_utc'] = datetime.now(timezone.utc).isoformat()
    atomic_json(args.out_dir / 'status.json', status)
    print(json.dumps(summary, indent=2), flush=True)
    if errors or changed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
