from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.table import Table


def create_history_template(footprints: Table) -> Table:
    history = Table()
    history["field_id"] = np.asarray(footprints["field_id"]).astype("U16")
    history["exposure_count"] = np.zeros(len(footprints), dtype=np.int32)
    history["last_observed_night"] = np.array([""] * len(footprints), dtype="U8")
    return history


def load_history(history_path: Path, footprints: Table) -> Table:
    if history_path.exists():
        return Table.read(history_path)
    return create_history_template(footprints)


def save_history(history: Table, history_path: Path) -> None:
    history_path.parent.mkdir(parents=True, exist_ok=True)
    history.write(history_path, overwrite=True)
