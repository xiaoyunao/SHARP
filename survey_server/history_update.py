from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from astropy.table import Table


MP_PATTERN = re.compile(r"MP[_-]?(\d+)", re.IGNORECASE)


def extract_mp_field_id(path: Path) -> str | None:
    match = MP_PATTERN.search(path.name)
    if match is None:
        return None
    number = int(match.group(1))
    return f"MP{number:04d}"


def count_l2_mp_exposures(l2_dir: Path) -> dict[str, int]:
    if not l2_dir.exists() or not l2_dir.is_dir():
        return {}
    counts: dict[str, int] = {}
    for path in sorted(l2_dir.iterdir()):
        if not path.is_file():
            continue
        field_id = extract_mp_field_id(path)
        if field_id is None:
            continue
        counts[field_id] = counts.get(field_id, 0) + 1
    return counts


def apply_counts_to_history(history: Table, counts: dict[str, int], obs_night: str) -> Table:
    if not counts:
        return history
    out = history.copy()
    field_ids = np.asarray(out["field_id"]).astype("U16")
    idx_by_field = {field_id: idx for idx, field_id in enumerate(field_ids)}
    for field_id, count in counts.items():
        idx = idx_by_field.get(field_id)
        if idx is None:
            continue
        out["exposure_count"][idx] += count
        out["last_observed_night"][idx] = obs_night
    return out
