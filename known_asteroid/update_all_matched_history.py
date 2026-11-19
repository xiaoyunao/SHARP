#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Maintain an incrementally updated all_matched_asteroids history table.")
    parser.add_argument("--processed-root", default="/processed1")
    parser.add_argument("--history-dir", default="/pipeline/xiaoyunao/known_asteroid/runtime/history")
    parser.add_argument("--night", default=None, help="Optional YYYYMMDD to prioritize updating one night.")
    return parser.parse_args()


def load_state(path: Path) -> dict[str, int]:
    if not path.exists():
        return {}
    return {str(k): int(v) for k, v in json.loads(path.read_text(encoding="utf-8")).items()}


def save_state(path: Path, state: dict[str, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary=path.with_suffix('.tmp.json')
    temporary.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def write_table(table,path):
    temporary=path.with_suffix('.tmp.fits')
    table.write(temporary,overwrite=True)
    temporary.replace(path)


def iter_matched_files(processed_root: Path) -> list[tuple[str, Path]]:
    rows: list[tuple[str, Path]] = []
    for path in sorted(processed_root.glob("20*/L4/*_matched_asteroids.fits")):
        night = path.name.split("_", 1)[0]
        if len(night) == 8 and night.isdigit() and not (processed_root/night/'ASTEROID_DATA_INVALID.json').exists():
            rows.append((night, path))
    return rows


def purge_excluded(processed_root: Path, all_path: Path, state_path: Path) -> None:
    excluded = {p.parent.name for p in processed_root.glob('20*/ASTEROID_DATA_INVALID.json')}
    if not excluded:
        return
    if all_path.exists():
        old = Table.read(all_path, memmap=True)
        if 'source_night' in old.colnames:
            keep = ~np.isin(np.asarray(old['source_night']).astype(str), list(excluded))
            if not keep.all():
                write_table(old[keep], all_path)
    if state_path.exists():
        state = load_state(state_path)
        save_state(state_path, {k:v for k,v in state.items() if k not in excluded})


def read_with_night(path: Path, night: str) -> Table:
    table = Table.read(path, memmap=True)
    if "source_night" not in table.colnames:
        table["source_night"] = np.array([night] * len(table), dtype="U8")
    else:
        table["source_night"] = np.asarray([str(x).strip() for x in table["source_night"]]).astype("U8")
    return table


def rebuild_all(all_path: Path, state_path: Path, matched_files: list[tuple[str, Path]]) -> None:
    tables: list[Table] = []
    state: dict[str, int] = {}
    for night, path in matched_files:
        try:
            table = read_with_night(path, night)
        except Exception:
            continue
        tables.append(table)
        state[night] = len(table)
    if tables:
        combined = vstack(tables, metadata_conflicts="silent")
        all_path.parent.mkdir(parents=True, exist_ok=True)
        write_table(combined,all_path)
    save_state(state_path, state)


def append_one(all_path: Path, state_path: Path, night: str, matched_path: Path, state: dict[str, int]) -> None:
    new_table = read_with_night(matched_path, night)
    if all_path.exists():
        old = Table.read(all_path, memmap=True)
        old=old[np.asarray(old['source_night']).astype(str)!=night]
        combined = vstack([old, new_table], metadata_conflicts="silent")
    else:
        combined = new_table
    all_path.parent.mkdir(parents=True, exist_ok=True)
    write_table(combined,all_path)
    state[night] = len(new_table)
    save_state(state_path, state)


def main() -> None:
    args = parse_args()
    processed_root = Path(args.processed_root)
    history_dir = Path(args.history_dir)
    history_dir.mkdir(parents=True, exist_ok=True)
    all_path = history_dir / "all_matched_asteroids.fits"
    state_path = history_dir / "all_matched_ingest_state.json"
    purge_excluded(processed_root, all_path, state_path)

    matched_files = iter_matched_files(processed_root)
    if not matched_files:
        return

    state = load_state(state_path)
    if args.night:
        matched_files = sorted(matched_files, key=lambda item: item[0] != args.night)

    need_rebuild = not all_path.exists()
    if not need_rebuild and state:
        for night, path in matched_files:
            if night in state:
                try:
                    row_count = len(Table.read(path, memmap=True))
                except Exception:
                    continue
                if int(state[night]) != int(row_count):
                    need_rebuild = True
                    break

    if need_rebuild:
        rebuild_all(all_path, state_path, matched_files)
        return

    # Reprocessing can change timing/matches without changing the number of rows.
    if args.night and args.night in state:
        for night,path in matched_files:
            if night==args.night:
                append_one(all_path,state_path,night,path,state)
                break

    changed = False
    for night, path in matched_files:
        if night in state:
            continue
        append_one(all_path, state_path, night, path, state)
        changed = True

    if not changed and not state_path.exists():
        save_state(state_path, state)


if __name__ == "__main__":
    main()
