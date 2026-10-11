"""Route measured candidates to human review without changing PSF cache provenance."""
from pathlib import Path

from unknown_morphology import approved_rows, row_hash

ROUTING_VERSION = 'pass-and-unresolved-to-human-20261011'


def reviewable_rows(night: str, processed_root: Path, rows: list[dict], strict=False):
    # Reuse the numerical engine's input/output and raw-row integrity checks.
    _, report = approved_rows(night, processed_root, rows)
    allowed = {r['row_sha256'] for r in report['candidates']
               if r['status'] in ('pass', 'review')}
    selected = [row for row in rows if row_hash(row) in allowed]
    if strict and len(selected) != len(rows):
        raise ValueError('Candidate contains morphology-rejected rows')
    return selected, report


def require_human_acceptance(rows: list[dict], review_csv: Path):
    from export_unknown_ades import read_review_csv, review_decision_for_row
    review, _ = read_review_csv(str(review_csv), require_values=False)
    if any(review_decision_for_row(row, review) != 1 for row in rows):
        raise ValueError('Every submitted candidate requires explicit human acceptance')
