"""Persistent observing-night exclusions, kept outside disposable L4 results."""
import json
from pathlib import Path

MARKER = 'ASTEROID_DATA_INVALID.json'


def exclusion(processed_root, night):
    path = Path(processed_root)/night/MARKER
    return json.loads(path.read_text()) if path.exists() else None


def purge_ledger(db, processed_root):
    for (night,) in db.execute('SELECT night FROM exposures UNION SELECT night FROM ingested').fetchall():
        if not exclusion(processed_root, night):
            continue
        db.execute('DELETE FROM detections WHERE uid IN (SELECT uid FROM exposures WHERE night=?)', (night,))
        db.execute('DELETE FROM exposures WHERE night=?', (night,))
        db.execute('DELETE FROM ingested WHERE night=?', (night,))
    db.commit()
