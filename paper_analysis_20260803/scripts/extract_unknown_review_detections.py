#!/usr/bin/env python3
"""Freeze compact unknown review detections with L2 photometry.

Designed for read-only execution on the production server.  The output root
must not overlap the production inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from astropy.io import fits


COLUMNS = {
    "trk_sub": "trk_sub",
    "linkage_id": "linkage_id",
    "detection_index": "detection_index",
    "image_name": "image_name",
    "catalog_name": "catalog_name",
    "obj_id": "objID",
    "source_tracklet_ids": "source_tracklet_ids",
    "link_tracklet_ids": "link_tracklet_ids",
    "groups": "groups",
    "exp_pairs": "exp_pairs",
    "mjd": "mjd",
    "ra_win_deg": "RA_Win",
    "dec_win_deg": "DEC_Win",
    "ra_psf_deg": "RA_PSF",
    "dec_psf_deg": "DEC_PSF",
    "ra_err_deg": "RAErr_Win",
    "dec_err_deg": "DECErr_Win",
    "x_win_px": "X_Win",
    "y_win_px": "Y_Win",
    "mag_aper4": "Mag_Aper4",
    "magerr_aper4": "MagErr_Aper4",
    "mag_psf": "Mag_PSF",
    "magerr_psf": "MagErr_PSF",
    "flag": "Flag",
    "n_obs": "n_obs",
    "n_tracklets": "n_tracklets",
    "rms_arcsec": "rms_arcsec",
    "med_arcsec": "med_arcsec",
    "max_arcsec": "max_arcsec",
    "a_au": "a_au",
    "ecc": "ecc",
    "inc_deg": "inc_deg",
    "best_v1_kms": "best_v1_kms",
    "lin_rms_arcsec": "lin_rms_arcsec",
    "lin_speed_arcsec_per_day": "lin_speed_arcsec_per_day",
    "lin_dir_deg": "lin_dir_deg",
}

UNITS = {
    "mjd": "d",
    "ra_win_deg": "deg",
    "dec_win_deg": "deg",
    "ra_psf_deg": "deg",
    "dec_psf_deg": "deg",
    "ra_err_deg": "deg",
    "dec_err_deg": "deg",
    "x_win_px": "pix",
    "y_win_px": "pix",
    "mag_aper4": "mag",
    "magerr_aper4": "mag",
    "mag_psf": "mag",
    "magerr_psf": "mag",
    "rms_arcsec": "arcsec",
    "med_arcsec": "arcsec",
    "max_arcsec": "arcsec",
    "a_au": "AU",
    "inc_deg": "deg",
    "best_v1_kms": "km/s",
    "lin_rms_arcsec": "arcsec",
    "lin_speed_arcsec_per_day": "arcsec/day",
    "lin_dir_deg": "deg",
}


def iter_nights(start: str, end: str):
    value = datetime.strptime(start, "%Y-%m-%d").date()
    finish = datetime.strptime(end, "%Y-%m-%d").date()
    while value <= finish:
        yield value.strftime("%Y%m%d")
        value += timedelta(days=1)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def native(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values)
    if array.dtype.byteorder not in ("=", "|"):
        array = array.byteswap().view(array.dtype.newbyteorder("="))
    if array.dtype.kind in {"S", "U", "O"}:
        return np.asarray(
            [
                item.decode("utf-8", "replace").strip()
                if isinstance(item, (bytes, bytearray))
                else str(item).strip()
                for item in array
            ],
            dtype=object,
        )
    return array


def with_unit_metadata(table: pa.Table) -> pa.Table:
    """Attach explicit physical units without changing column values or types."""

    fields = []
    for field in table.schema:
        metadata = dict(field.metadata or {})
        unit = UNITS.get(field.name)
        if unit is not None:
            metadata[b"unit"] = unit.encode("ascii")
        fields.append(field.with_metadata(metadata or None))
    return table.cast(pa.schema(fields, metadata=table.schema.metadata))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--review-root",
        type=Path,
        default=Path("/pipeline/xiaoyunao/heliolincrr/review_packages"),
    )
    parser.add_argument("--start", default="2025-11-15")
    parser.add_argument("--end", default="2026-07-15")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    review_root = args.review_root.resolve()
    if output == review_root or review_root in output.parents or output in review_root.parents:
        raise ValueError("output directory must not overlap review root")
    output.mkdir(parents=True, exist_ok=False)
    temp = output / "unknown_review_detections.parquet.inprogress"
    final = output / "unknown_review_detections.parquet"
    writer = None
    statuses = []
    total = 0
    try:
        for night in iter_nights(args.start, args.end):
            path = review_root / night / f"{night}_unknown_review_full.fits"
            if not path.exists():
                statuses.append({"night": night, "status": "missing", "path": str(path), "rows": 0})
                continue
            before = path.stat()
            try:
                with fits.open(path, memmap=True) as hdul:
                    data = hdul[1].data
                    names = set(hdul[1].columns.names)
                    missing = [source for source in COLUMNS.values() if source not in names]
                    if missing:
                        raise KeyError(f"missing columns: {missing}")
                    if len(data) == 0:
                        after = path.stat()
                        statuses.append(
                            {
                                "night": night,
                                "status": "ok_zero",
                                "path": str(path),
                                "rows": 0,
                                "size_bytes": before.st_size,
                                "mtime_utc": datetime.fromtimestamp(before.st_mtime, timezone.utc).isoformat(),
                                "sha256": sha256(path),
                                "changed_during_read": (before.st_size, before.st_mtime_ns)
                                != (after.st_size, after.st_mtime_ns),
                            }
                        )
                        continue
                    payload = {
                        "night": np.repeat(night, len(data)),
                        "provenance_path": np.repeat(str(path), len(data)),
                        "provenance_row": np.arange(len(data), dtype=np.int64),
                    }
                    payload.update({target: native(data[source]) for target, source in COLUMNS.items()})
                    table = with_unit_metadata(pa.Table.from_pydict(payload))
                if writer is None:
                    writer = pq.ParquetWriter(temp, table.schema, compression="zstd", compression_level=6)
                writer.write_table(table)
                total += len(table)
                after = path.stat()
                changed = (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns)
                statuses.append(
                    {
                        "night": night,
                        "status": "changed_during_read" if changed else "ok",
                        "path": str(path),
                        "rows": len(table),
                        "size_bytes": before.st_size,
                        "mtime_utc": datetime.fromtimestamp(before.st_mtime, timezone.utc).isoformat(),
                        "sha256": sha256(path),
                    }
                )
            except Exception as exc:
                statuses.append(
                    {"night": night, "status": "error", "path": str(path), "rows": 0, "error": f"{type(exc).__name__}: {exc}"}
                )
        if writer is None:
            raise RuntimeError("no review detections were found")
        writer.close()
        writer = None
        os.replace(temp, final)
    finally:
        if writer is not None:
            writer.close()
        if temp.exists():
            temp.unlink()

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "rows": total,
        "nights_present": sum(row["status"] in {"ok", "ok_zero"} for row in statuses),
        "nonempty_nights": sum(row["status"] == "ok" for row in statuses),
        "zero_nights": sum(row["status"] == "ok_zero" for row in statuses),
        "errors": sum(row["status"] in {"error", "changed_during_read"} for row in statuses),
        "parquet_sha256": sha256(final),
        "parquet_size_bytes": final.stat().st_size,
        "files": statuses,
    }
    (output / "unknown_review_detections_status.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: summary[key] for key in ("rows", "nights_present", "errors", "parquet_size_bytes")}, indent=2))


if __name__ == "__main__":
    main()
