#!/usr/bin/env python3
from __future__ import annotations

import argparse
import fcntl
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from astropy.table import Table


ALPHABET = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
TRK_SUB_WIDTH = 7
TRK_SUB_RE = re.compile(rf"^[0-9A-Za-z]{{{TRK_SUB_WIDTH}}}$")
DEFAULT_HISTORY = "/pipeline/xiaoyunao/data/heliolincrr/trksub_history.jsonl"


def int_to_base62(value: int) -> str:
    if value <= 0:
        raise ValueError("trkSub sequence value must be positive")
    digits: list[str] = []
    base = len(ALPHABET)
    while value:
        value, rem = divmod(value, base)
        digits.append(ALPHABET[rem])
    out = "".join(reversed(digits)).rjust(TRK_SUB_WIDTH, "0")
    if len(out) > TRK_SUB_WIDTH:
        raise OverflowError(f"trkSub sequence exhausted {TRK_SUB_WIDTH} base62 digits")
    return out


def base62_to_int(value: str) -> int:
    value = str(value).strip()
    if not TRK_SUB_RE.fullmatch(value):
        raise ValueError(f"Invalid managed trkSub {value!r}; expected exactly {TRK_SUB_WIDTH} base62 characters")
    out = 0
    base = len(ALPHABET)
    for char in value:
        out = out * base + ALPHABET.index(char)
    return out


def split_semicolon(value: Any) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def load_unknown_rows(path: Path) -> list[dict[str, Any]]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError(f"Unknown catalog must contain a JSON list: {path}")
    return rows


def write_unknown_rows(rows: list[dict[str, Any]], json_path: Path, fits_path: Path | None) -> None:
    json_path.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if fits_path:
        if rows:
            Table(rows=rows).write(fits_path, overwrite=True)
        elif fits_path.exists():
            Table.read(fits_path)[:0].write(fits_path, overwrite=True)


def detection_signature(row: dict[str, Any]) -> str:
    image_names = split_semicolon(row.get("image_names"))
    objids = split_semicolon(row.get("objids"))
    mjds = split_semicolon(row.get("mjds"))
    ras = split_semicolon(row.get("ras_deg"))
    decs = split_semicolon(row.get("decs_deg"))
    parts: list[str] = []
    for idx, image_name in enumerate(image_names):
        objid = objids[idx] if idx < len(objids) else ""
        mjd = mjds[idx] if idx < len(mjds) else ""
        ra = ras[idx] if idx < len(ras) else ""
        dec = decs[idx] if idx < len(decs) else ""
        parts.append(f"{image_name}:{objid}:{mjd}:{ra}:{dec}")
    return "|".join(sorted(parts))


def identity_key(row: dict[str, Any], mode: str, night: str) -> str:
    return f"{mode}:{night}:{detection_signature(row)}"


def load_history(path: Path) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], int]:
    records: list[dict[str, Any]] = []
    by_identity: dict[str, dict[str, Any]] = {}
    max_value = 0
    if not path.exists():
        return records, by_identity, max_value
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        record = json.loads(line)
        trk_sub = record.get("trk_sub", "")
        value = base62_to_int(trk_sub)
        max_value = max(max_value, value)
        key = str(record.get("identity_key", ""))
        if key:
            by_identity[key] = record
        records.append(record)
    return records, by_identity, max_value


def record_for_row(row: dict[str, Any], args: argparse.Namespace, trk_sub: str, key: str) -> dict[str, Any]:
    return {
        "history_version": 1,
        "trk_sub": trk_sub,
        "sequence_value": base62_to_int(trk_sub),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "mode": args.mode,
        "night": args.night,
        "identity_key": key,
        "source_catalog": str(Path(args.catalog).resolve()),
        "source_summary": str(Path(args.summary_json).resolve()) if args.summary_json else "",
        "linkage_id": int(row["linkage_id"]) if str(row.get("linkage_id", "")).strip() else None,
        "internal_tracklet_ids": split_semicolon(row.get("tracklet_ids")),
        "image_names": split_semicolon(row.get("image_names")),
        "objids": split_semicolon(row.get("objids")),
        "mjds": split_semicolon(row.get("mjds")),
        "ras_deg": split_semicolon(row.get("ras_deg")),
        "decs_deg": split_semicolon(row.get("decs_deg")),
        "n_tracklets": row.get("n_tracklets"),
        "n_obs": row.get("n_obs"),
        "fit_ok": row.get("fit_ok"),
        "is_good": row.get("is_good"),
        "rms_arcsec": row.get("rms_arcsec"),
        "a_au": row.get("a_au"),
        "ecc": row.get("ecc"),
        "inc_deg": row.get("inc_deg"),
    }


def assign(args: argparse.Namespace) -> dict[str, int]:
    catalog = Path(args.catalog)
    fits_path = Path(args.fits) if args.fits else None
    history = Path(args.history)
    history.parent.mkdir(parents=True, exist_ok=True)
    lock_path = Path(str(history) + ".lock")

    rows = load_unknown_rows(catalog)
    with lock_path.open("a+", encoding="utf-8") as lock_fh:
        fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX)
        records, by_identity, max_value = load_history(history)
        next_value = max(max_value + 1, 1)
        new_records: list[dict[str, Any]] = []
        reused = 0
        assigned = 0

        for row in rows:
            key = identity_key(row, args.mode, args.night)
            existing = by_identity.get(key)
            if existing is not None:
                row["trk_sub"] = existing["trk_sub"]
                reused += 1
                continue

            current = str(row.get("trk_sub", "")).strip()
            if current:
                if not TRK_SUB_RE.fullmatch(current):
                    raise ValueError(f"Existing trk_sub {current!r} is not managed {TRK_SUB_WIDTH}-character base62")
                trk_sub = current
                next_value = max(next_value, base62_to_int(trk_sub) + 1)
            else:
                trk_sub = int_to_base62(next_value)
                next_value += 1
                row["trk_sub"] = trk_sub
            rec = record_for_row(row, args, trk_sub, key)
            new_records.append(rec)
            by_identity[key] = rec
            assigned += 1

        if not args.dry_run and new_records:
            with history.open("a", encoding="utf-8") as fh:
                for rec in new_records:
                    fh.write(json.dumps(rec, ensure_ascii=False, sort_keys=True) + "\n")
        if not args.dry_run:
            write_unknown_rows(rows, catalog, fits_path)

    return {
        "catalog_rows": len(rows),
        "assigned_new": assigned,
        "reused_existing": reused,
        "history_records_before": len(records),
        "history_records_added": len(new_records),
    }


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Assign globally unique MPC trkSub values to unknown catalog rows.")
    ap.add_argument("night", help="Target night in YYYYMMDD format")
    ap.add_argument("--mode", default="single-night", help="Pipeline mode stored in history, e.g. single-night or window15")
    ap.add_argument("--catalog", required=True, help="Unknown catalog JSON to update in place")
    ap.add_argument("--fits", default="", help="Optional unknown catalog FITS to update in place")
    ap.add_argument("--summary-json", default="", help="Optional summary JSON path recorded in history")
    ap.add_argument("--history", default=DEFAULT_HISTORY, help=f"Global trkSub history JSONL (default: {DEFAULT_HISTORY})")
    ap.add_argument("--dry-run", action="store_true")
    return ap


def main() -> None:
    args = build_argparser().parse_args()
    stats = assign(args)
    stats["history"] = str(Path(args.history))
    stats["catalog"] = str(Path(args.catalog))
    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
