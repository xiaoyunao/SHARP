#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
from astropy.table import Table


OLD_MANAGED_RE = re.compile(r"^0[0-9A-Za-z]{7}$")
TEXT_SUFFIXES = {".json", ".jsonl", ".csv", ".psv", ".txt", ".md", ".tsv"}
FITS_SUFFIXES = {".fits", ".fit"}
DEFAULT_ROOTS = [
    "/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl",
    "/pipeline/xiaoyunao/data/heliolincrr/20*/analysis/*single_night_summary*.json",
    "/pipeline/xiaoyunao/data/heliolincrr/20*/analysis/*single_night_summary*.txt",
    "/pipeline/xiaoyunao/heliolincrr/review_packages/**/*.csv",
    "/pipeline/xiaoyunao/heliolincrr/review_packages/**/*.json",
    "/pipeline/xiaoyunao/heliolincrr/review_packages/**/*.fits",
    "/pipeline/xiaoyunao/heliolincrr/review_packages/**/*.psv",
    "/pipeline/xiaoyunao/heliolincrr/review_packages/**/*.txt",
    "/processed1/20*/L4/*unknown*.json",
    "/processed1/20*/L4/*unknown*.fits",
    "/processed1/20*/L4/*unknown*.psv",
    "/processed1/20*/L4/*unknown*reply*.txt",
]
DEFAULT_RENAME_ROOTS = ["/pipeline/xiaoyunao/heliolincrr/review_packages"]


def convert_value(value: Any) -> Any:
    if isinstance(value, str):
        return value[1:] if OLD_MANAGED_RE.fullmatch(value) else value
    if isinstance(value, list):
        return [convert_value(item) for item in value]
    if isinstance(value, dict):
        return {key: convert_value(item) for key, item in value.items()}
    return value


def converted_text(text: str) -> tuple[str, int]:
    count = 0

    def repl(match: re.Match[str]) -> str:
        nonlocal count
        count += 1
        return match.group(0)[1:]

    return re.sub(r"(?<![0-9A-Za-z])0[0-9A-Za-z]{7}(?![0-9A-Za-z])", repl, text), count


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def migrate_json(path: Path, dry_run: bool) -> int:
    data = json.loads(path.read_text(encoding="utf-8"))
    new_data = convert_value(data)
    if new_data == data:
        return 0
    if not dry_run:
        write_json(path, new_data)
    return 1


def migrate_jsonl(path: Path, dry_run: bool) -> int:
    changed = 0
    out_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            out_lines.append(line)
            continue
        data = json.loads(line)
        new_data = convert_value(data)
        changed += int(new_data != data)
        out_lines.append(json.dumps(new_data, ensure_ascii=False, sort_keys=True))
    if changed and not dry_run:
        path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return changed


def migrate_text(path: Path, dry_run: bool) -> int:
    text = path.read_text(encoding="utf-8", errors="ignore")
    new_text, count = converted_text(text)
    if count and not dry_run:
        path.write_text(new_text, encoding="utf-8")
    return count


def migrate_fits(path: Path, dry_run: bool) -> int:
    table = Table.read(path)
    changed = 0
    for col in table.colnames:
        dtype = table[col].dtype
        if dtype.kind not in {"U", "S", "O"}:
            continue
        values = np.asarray(table[col]).astype(str)
        new_values = np.asarray([convert_value(value) for value in values])
        mask = new_values != values
        if np.any(mask):
            table[col] = new_values
            changed += int(np.count_nonzero(mask))
    if changed and not dry_run:
        table.write(path, overwrite=True)
    return changed


def migrate_csv_with_reader(path: Path, dry_run: bool) -> int:
    text = path.read_text(encoding="utf-8", errors="ignore")
    sample = text[:2048]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t|") if sample.strip() else csv.excel
    except csv.Error:
        return migrate_text(path, dry_run)
    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    if not rows or not rows[0]:
        return migrate_text(path, dry_run)
    changed = 0
    new_rows: list[dict[str, Any]] = []
    for row in rows:
        new_row = {key: convert_value(value) for key, value in row.items()}
        changed += int(new_row != row)
        new_rows.append(new_row)
    if changed and not dry_run:
        with path.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), dialect=dialect)
            writer.writeheader()
            writer.writerows(new_rows)
    return changed


def candidate_files(roots: list[Path]) -> list[Path]:
    out: list[Path] = []
    for root in roots:
        root_text = str(root)
        if any(char in root_text for char in "*?["):
            out.extend(Path(item) for item in glob.glob(root_text, recursive=True) if Path(item).is_file())
            continue
        if root.is_file():
            out.append(root)
            continue
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            suffix = path.suffix.lower()
            if suffix in TEXT_SUFFIXES or suffix in FITS_SUFFIXES:
                out.append(path)
    return sorted(out)


def migrate_filenames(roots: list[Path], dry_run: bool) -> tuple[int, int]:
    changed_files = 0
    changed_items = 0
    for root in roots:
        if not root.exists():
            continue
        paths = [root] if root.is_file() else sorted(root.rglob("*"), reverse=True)
        for path in paths:
            new_name, count = converted_text(path.name)
            if not count or new_name == path.name:
                continue
            target = path.with_name(new_name)
            if target.exists():
                raise FileExistsError(f"Cannot rename {path} -> {target}: target exists")
            changed_files += 1
            changed_items += count
            print(f"[RENAME] {path} -> {target}")
            if not dry_run:
                path.rename(target)
    return changed_files, changed_items


def migrate_file(path: Path, dry_run: bool) -> int:
    suffix = path.suffix.lower()
    if suffix == ".json":
        return migrate_json(path, dry_run)
    if suffix == ".jsonl":
        return migrate_jsonl(path, dry_run)
    if suffix == ".csv":
        return migrate_csv_with_reader(path, dry_run)
    if suffix in TEXT_SUFFIXES:
        return migrate_text(path, dry_run)
    if suffix in FITS_SUFFIXES:
        return migrate_fits(path, dry_run)
    return 0


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Migrate managed unknown trkSub values from 8-char leading-zero form to 7-char MPC form."
    )
    ap.add_argument("roots", nargs="*", default=DEFAULT_ROOTS, help="Files or directories to scan")
    ap.add_argument(
        "--rename-root",
        action="append",
        default=DEFAULT_RENAME_ROOTS,
        help="Directory whose filenames should also be migrated; repeatable",
    )
    ap.add_argument("--no-rename-files", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="Optional file count limit for testing")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    files = candidate_files([Path(root) for root in args.roots])
    if args.limit:
        files = files[: args.limit]
    changed_files = 0
    changed_items = 0
    for path in files:
        try:
            changed = migrate_file(path, dry_run=bool(args.dry_run))
        except Exception as exc:
            print(f"[WARN] skip {path}: {exc}")
            continue
        if changed:
            changed_files += 1
            changed_items += changed
            print(f"[CHANGE] {path} items={changed}")
    renamed_files = 0
    renamed_items = 0
    if not args.no_rename_files:
        renamed_files, renamed_items = migrate_filenames([Path(root) for root in args.rename_root], bool(args.dry_run))
    print(
        json.dumps(
            {
                "dry_run": bool(args.dry_run),
                "files_scanned": len(files),
                "changed_files": changed_files,
                "changed_items": changed_items,
                "renamed_files": renamed_files,
                "renamed_items": renamed_items,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
