#!/usr/bin/env python3
"""Stream selected frozen production products into compact Parquet tables.

This program is intended to run on the SHARP production server.  All source
FITS files are opened read-only.  Output is permitted only beneath the explicit
``--output-dir`` and is rejected if that directory is inside either production
input root.  Existing outputs are never overwritten.

Six normalized Parquet tables are produced:

* ``known_all.parquet``
* ``known_matched.parquet``
* ``known_mask15.parquet``
* ``unknown_catalog.parquet``
* ``orbit_links.parquet``
* ``orbit_obs_residuals.parquet``

Every row carries the observing night, product name, input-product provenance,
source-row number, and a ``source_file`` value.  For known-object products,
``source_file`` is the originating L2 catalog.  For the unknown catalog it is
the semicolon-delimited exposure list.  Orbit products do not carry exposure
names, so their ``source_file`` is the input FITS basename; ``provenance_path``
always identifies the exact input product without this semantic fallback.

The script scans every expected product for every inclusive night in the date
range.  Missing files are normal audit records.  Read or schema failures are
recorded and scanning continues, but the final process exit status is nonzero.
No production file or directory is created, renamed, or modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import socket
import sys
from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from astropy import __version__ as astropy_version
from astropy.io import fits


SCHEMA_VERSION = "1.0"
DEFAULT_START = "2025-11-15"
DEFAULT_END = "2026-07-15"
DEFAULT_PROCESSED_ROOT = "/processed1"
DEFAULT_HELIOLINC_ROOT = "/pipeline/xiaoyunao/data/heliolincrr"
DEFAULT_PROTECTED_ROOTS = (
    Path("/raw1"),
    Path("/processed1"),
    Path("/pipeline/xiaoyunao"),
)


@dataclass(frozen=True)
class FieldSpec:
    """One normalized output field and its accepted FITS-column aliases."""

    name: str
    aliases: tuple[str, ...]
    kind: str
    required: bool = False
    unit: str | None = None
    description: str | None = None


@dataclass(frozen=True)
class ProductSpec:
    """Location and output schema for one nightly production product."""

    name: str
    output_name: str
    root_kind: str
    relative_template: str
    fields: tuple[FieldSpec, ...]
    source_aliases: tuple[str, ...] = ()


def field(
    name: str,
    *aliases: str,
    kind: str = "float64",
    required: bool = False,
    unit: str | None = None,
    description: str | None = None,
) -> FieldSpec:
    return FieldSpec(
        name=name,
        aliases=tuple(aliases),
        kind=kind,
        required=required,
        unit=unit,
        description=description,
    )


PREDICTED_ASTEROID_FIELDS = (
    field("asteroid_name", "name", "astorb_name", kind="string"),
    field("asteroid_number", "number", "num", "astorb_number", kind="identifier"),
    field("epoch_mjd", "epoch", "mjd", "MJD", required=True, unit="d"),
    field("pred_ra_deg", "ra", "RA", required=True, unit="deg"),
    field("pred_dec_deg", "dec", "DEC", required=True, unit="deg"),
    field("pred_v_mag", "mag", "V", unit="mag"),
    field(
        "pred_heliocentric_distance_au",
        "r",
        "r_au",
        "heliocentric_distance",
        unit="AU",
    ),
    field(
        "pred_geocentric_distance_au",
        "delta",
        "delta_au",
        "geocentric_distance",
        unit="AU",
    ),
    field(
        "pred_solar_elongation_deg",
        "elongation",
        "elong",
        "solar_elongation",
        unit="deg",
    ),
    field("pred_phase_angle_deg", "phase", "alpha", "phase_angle", unit="deg"),
)


MATCHED_DETECTION_FIELDS = (
    field("obs_ra_win_deg", "RA_Win", required=True, unit="deg"),
    field("obs_dec_win_deg", "DEC_Win", required=True, unit="deg"),
    field("obs_ra_err_deg", "RAErr_Win", unit="deg"),
    field("obs_dec_err_deg", "DECErr_Win", unit="deg"),
    field("obs_ra_psf_deg", "RA_PSF", unit="deg"),
    field("obs_dec_psf_deg", "DEC_PSF", unit="deg"),
    field("x_win_px", "X_Win", "XWIN_IMAGE", unit="pix"),
    field("y_win_px", "Y_Win", "YWIN_IMAGE", unit="pix"),
    field("obj_id", "objID", "NUMBER", "ID", kind="identifier"),
    field("mag_kron", "Mag_Kron", "MAG_KRON", unit="mag"),
    field("magerr_kron", "MagErr_Kron", "MAGERR_KRON", unit="mag"),
    field("mag_aper4", "Mag_Aper4", "MAG_APER4", unit="mag"),
    field("magerr_aper4", "MagErr_Aper4", "MAGERR_APER4", unit="mag"),
    field("mag_aper8", "Mag_Aper8", "MAG_APER8", unit="mag"),
    field("magerr_aper8", "MagErr_Aper8", "MAGERR_APER8", unit="mag"),
    field("mag_psf", "Mag_PSF", "MAG_PSF", unit="mag"),
    field("magerr_psf", "MagErr_PSF", "MAGERR_PSF", unit="mag"),
    field("flux_aper4", "Flux_Aper4", "FLUX_APER4"),
    field("fluxerr_aper4", "FluxErr_Aper4", "FLUXERR_APER4"),
    field("fwhm_px", "FWHM", "FWHM_IMAGE", unit="pix"),
    field("flag", "Flag", "FLAGS", kind="int64"),
    field("flag_iso_num", "Flag_ISO_Num", kind="int64"),
)


UNKNOWN_CATALOG_FIELDS = (
    field("trk_sub", "trk_sub", "trkSub", kind="string"),
    field("linkage_id", "linkage_id", kind="int64", required=True),
    field("n_tracklets", "n_tracklets", kind="int64"),
    field("n_obs", "n_obs", kind="int64"),
    field("fit_ok", "fit_ok", kind="bool", required=True),
    field("is_good", "is_good", kind="bool"),
    field("rms_arcsec", "rms_arcsec", unit="arcsec"),
    field("med_arcsec", "med_arcsec", unit="arcsec"),
    field("max_arcsec", "max_arcsec", unit="arcsec"),
    field("a_au", "a_au", unit="AU"),
    field("ecc", "ecc"),
    field("inc_deg", "inc_deg", unit="deg"),
    field("raan_deg", "raan_deg", unit="deg"),
    field("argp_deg", "argp_deg", unit="deg"),
    field("nu_deg", "nu_deg", unit="deg"),
    field("best_v1_kms", "best_v1_kms", unit="km/s"),
    field("lin_rms_arcsec", "lin_rms_arcsec", unit="arcsec"),
    field(
        "lin_speed_arcsec_per_day",
        "lin_speed_arcsec_per_day",
        unit="arcsec/day",
    ),
    field("lin_dir_deg", "lin_dir_deg", unit="deg"),
    field("tracklet_ids", "tracklet_ids", kind="string"),
    field("image_names", "image_names", kind="string"),
    field("objids", "objids", kind="string"),
    field("mjds", "mjds", kind="string", unit="d"),
    field("ras_deg", "ras_deg", kind="string", unit="deg"),
    field("decs_deg", "decs_deg", kind="string", unit="deg"),
    field("groups", "groups", kind="string"),
    field("exp_pairs", "exp_pairs", kind="string"),
)


ORBIT_LINK_FIELDS = (
    field("linkage_id", "linkage_id", kind="int64", required=True),
    field("n_tracklets", "n_tracklets", kind="int64"),
    field("n_nights", "n_nights", kind="int64"),
    field("n_obs", "n_obs", kind="int64"),
    field("fit_ok", "fit_ok", kind="bool", required=True),
    field("hypo_r_au", "hypo_r_au", unit="AU"),
    field("hypo_rdot_au_day", "hypo_rdot_au_day", unit="AU/day"),
    field("hypo_rdd_au_day2", "hypo_rdd_au_day2", unit="AU/day2"),
    field("rms_arcsec", "rms_arcsec", unit="arcsec"),
    field("med_arcsec", "med_arcsec", unit="arcsec"),
    field("max_arcsec", "max_arcsec", unit="arcsec"),
    field("is_good", "is_good", kind="bool", required=True),
    field("a_au", "a_au", unit="AU"),
    field("ecc", "ecc"),
    field("inc_deg", "inc_deg", unit="deg"),
    field("raan_deg", "raan_deg", unit="deg"),
    field("argp_deg", "argp_deg", unit="deg"),
    field("nu_deg", "nu_deg", unit="deg"),
    field("pred_ra_deg", "pred_ra_deg", unit="deg"),
    field("pred_dec_deg", "pred_dec_deg", unit="deg"),
    field("lin_rms_arcsec", "lin_rms_arcsec", unit="arcsec"),
    field(
        "lin_speed_arcsec_per_day",
        "lin_speed_arcsec_per_day",
        unit="arcsec/day",
    ),
    field("lin_dir_deg", "lin_dir_deg", unit="deg"),
    field("best_v1_kms", "best_v1_kms", unit="km/s"),
    field("min_rejected_max_v_kms", "min_rejected_max_v_kms", unit="km/s"),
    field("best_rejected_max_v_kms", "best_rejected_max_v_kms", unit="km/s"),
    field("fail_reason", "fail_reason", kind="string"),
    field("fail_counts", "fail_counts", kind="string"),
    field("seed_tracklets", "seed_tracklets", kind="string"),
    field("inlier_tracklets", "inlier_tracklets", kind="string"),
    field("rejected_tracklets", "rejected_tracklets", kind="string"),
    field("seed_max_v_kms", "seed_max_v_kms", unit="km/s"),
    field("final_max_v_kms", "final_max_v_kms", unit="km/s"),
)


ORBIT_RESIDUAL_FIELDS = (
    field("linkage_id", "linkage_id", kind="int64", required=True),
    field("obs_key", "obs_key", kind="string"),
    field("tracklet_id", "tracklet_id", kind="string"),
    field("mjd", "mjd", required=True, unit="d"),
    field("ra_deg", "ra_deg", required=True, unit="deg"),
    field("dec_deg", "dec_deg", required=True, unit="deg"),
    field("resid_arcsec", "resid_arcsec", required=True, unit="arcsec"),
    field("used", "used", kind="bool", required=True),
)


PRODUCTS = (
    ProductSpec(
        name="known_all",
        output_name="known_all.parquet",
        root_kind="processed",
        relative_template="{night}/L4/{night}_all_asteroids.fits",
        fields=PREDICTED_ASTEROID_FIELDS,
        source_aliases=("source_file", "srcfile"),
    ),
    ProductSpec(
        name="known_matched",
        output_name="known_matched.parquet",
        root_kind="processed",
        relative_template="{night}/L4/{night}_matched_asteroids.fits",
        fields=PREDICTED_ASTEROID_FIELDS + MATCHED_DETECTION_FIELDS,
        source_aliases=("source_file", "srcfile"),
    ),
    ProductSpec(
        name="known_mask15",
        output_name="known_mask15.parquet",
        root_kind="processed",
        relative_template="{night}/L4/{night}_matched_asteroids_mask15.fits",
        fields=PREDICTED_ASTEROID_FIELDS + MATCHED_DETECTION_FIELDS,
        source_aliases=("source_file", "srcfile"),
    ),
    ProductSpec(
        name="unknown_catalog",
        output_name="unknown_catalog.parquet",
        root_kind="processed",
        relative_template="{night}/L4/{night}_unknown_links.fits",
        fields=UNKNOWN_CATALOG_FIELDS,
        source_aliases=("image_names",),
    ),
    ProductSpec(
        name="orbit_links",
        output_name="orbit_links.parquet",
        root_kind="heliolinc",
        relative_template="{night}/{rr_subdir}/orbit_confirm/orbit_links.fits",
        fields=ORBIT_LINK_FIELDS,
    ),
    ProductSpec(
        name="orbit_obs_residuals",
        output_name="orbit_obs_residuals.parquet",
        root_kind="heliolinc",
        relative_template="{night}/{rr_subdir}/orbit_confirm/orbit_obs_residuals.fits",
        fields=ORBIT_RESIDUAL_FIELDS,
    ),
)


PROVENANCE_FIELDS = (
    pa.field("night", pa.string(), nullable=False),
    pa.field("product", pa.string(), nullable=False),
    pa.field("provenance_path", pa.string(), nullable=False),
    pa.field("provenance_basename", pa.string(), nullable=False),
    pa.field("provenance_size_bytes", pa.int64(), nullable=False),
    pa.field("provenance_mtime_utc", pa.string(), nullable=False),
    pa.field("source_row", pa.int64(), nullable=False),
    pa.field("source_file", pa.string(), nullable=True),
)


def arrow_type(kind: str) -> pa.DataType:
    if kind in {"string", "identifier"}:
        return pa.string()
    if kind == "float64":
        return pa.float64()
    if kind == "int64":
        return pa.int64()
    if kind == "bool":
        return pa.bool_()
    raise ValueError(f"Unsupported field kind: {kind}")


def product_schema(product: ProductSpec) -> pa.Schema:
    fields = list(PROVENANCE_FIELDS)
    for spec in product.fields:
        metadata: dict[bytes, bytes] = {}
        if spec.unit:
            metadata[b"unit"] = spec.unit.encode("utf-8")
        if spec.description:
            metadata[b"description"] = spec.description.encode("utf-8")
        fields.append(
            pa.field(
                spec.name,
                arrow_type(spec.kind),
                nullable=not spec.required,
                metadata=metadata or None,
            )
        )
    metadata = {
        b"schema_version": SCHEMA_VERSION.encode("ascii"),
        b"source_product": product.name.encode("utf-8"),
        b"extractor": Path(__file__).name.encode("utf-8"),
    }
    return pa.schema(fields, metadata=metadata)


class StreamingParquetWriter:
    """Write row groups incrementally, then atomically publish one Parquet file."""

    def __init__(
        self,
        output_path: Path,
        schema: pa.Schema,
        compression: str | None,
        compression_level: int | None,
    ) -> None:
        self.output_path = output_path
        self.inprogress_path = output_path.with_name(output_path.name + ".inprogress")
        self.schema = schema
        self.compression = compression
        self.compression_level = compression_level
        self.writer: pq.ParquetWriter | None = None
        self.rows_written = 0

    def _open(self) -> None:
        if self.writer is not None:
            return
        kwargs: dict[str, Any] = {
            "where": self.inprogress_path,
            "schema": self.schema,
            "compression": self.compression,
            "use_dictionary": True,
            "write_statistics": True,
        }
        if self.compression_level is not None:
            kwargs["compression_level"] = self.compression_level
        self.writer = pq.ParquetWriter(**kwargs)

    def write(self, table: pa.Table) -> None:
        if table.num_rows == 0:
            return
        self._open()
        assert self.writer is not None
        self.writer.write_table(table, row_group_size=table.num_rows)
        self.rows_written += table.num_rows

    def finish(self) -> None:
        self._open()
        assert self.writer is not None
        self.writer.close()
        self.writer = None
        os.replace(self.inprogress_path, self.output_path)

    def abort(self) -> None:
        if self.writer is not None:
            self.writer.close()
            self.writer = None


def parse_iso_date(value: str) -> date:
    try:
        return datetime.strptime(value, "%Y-%m-%d").date()
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected YYYY-MM-DD, got {value!r}") from exc


def iter_nights(start: date, end: date) -> Iterable[str]:
    current = start
    while current <= end:
        yield current.strftime("%Y%m%d")
        current += timedelta(days=1)


def is_same_or_descendant(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def iso_mtime(stat_result: os.stat_result) -> str:
    return datetime.fromtimestamp(stat_result.st_mtime, tz=timezone.utc).isoformat()


def file_sha256(path: Path, block_bytes: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_bytes)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def resolve_column_map(names: Sequence[str], fields: Sequence[FieldSpec]) -> dict[str, str | None]:
    exact = {name: name for name in names}
    folded: dict[str, str] = {}
    for name in names:
        folded.setdefault(name.casefold(), name)

    resolved: dict[str, str | None] = {}
    for spec in fields:
        actual = None
        for alias in spec.aliases:
            if alias in exact:
                actual = exact[alias]
                break
            if alias.casefold() in folded:
                actual = folded[alias.casefold()]
                break
        resolved[spec.name] = actual
    return resolved


def resolve_first_column(names: Sequence[str], aliases: Sequence[str]) -> str | None:
    mapping = resolve_column_map(
        names,
        (FieldSpec("value", tuple(aliases), "string"),),
    )
    return mapping["value"]


def scalar_to_text(value: Any, identifier: bool = False) -> str | None:
    if value is None or value is np.ma.masked:
        return None
    if isinstance(value, (bytes, np.bytes_)):
        text = bytes(value).decode("utf-8", errors="replace").strip()
        return text or None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float):
        if not np.isfinite(value):
            return None
        if identifier and value.is_integer():
            return str(int(value))
    text = str(value).strip()
    return text or None


def mask_for(values: Any) -> np.ndarray:
    if np.ma.isMaskedArray(values):
        return np.asarray(np.ma.getmaskarray(values), dtype=bool)
    return np.zeros(len(values), dtype=bool)


def arrow_array(values: Any, spec: FieldSpec) -> pa.Array:
    """Convert one FITS column slice to the fixed Arrow type."""

    masked = np.ma.asarray(values)
    data = np.ma.getdata(masked)
    mask = mask_for(masked)

    if spec.kind in {"string", "identifier"}:
        identifier = spec.kind == "identifier"
        converted = [
            None if mask[index] else scalar_to_text(value, identifier=identifier)
            for index, value in enumerate(data)
        ]
        return pa.array(converted, type=pa.string())

    if spec.kind == "float64":
        converted = np.asarray(data, dtype=np.float64)
        return pa.array(converted, mask=mask, type=pa.float64())

    if spec.kind == "int64":
        converted: list[int | None] = []
        for index, value in enumerate(data):
            if mask[index] or value is None or value is np.ma.masked:
                converted.append(None)
                continue
            try:
                numeric = float(value)
            except (TypeError, ValueError):
                converted.append(None)
                continue
            converted.append(int(numeric) if np.isfinite(numeric) else None)
        return pa.array(converted, type=pa.int64())

    if spec.kind == "bool":
        converted_bool: list[bool | None] = []
        for index, value in enumerate(data):
            if mask[index] or value is None or value is np.ma.masked:
                converted_bool.append(None)
                continue
            if isinstance(value, (bytes, np.bytes_)):
                value = bytes(value).decode("ascii", errors="ignore")
            if isinstance(value, str):
                text = value.strip().lower()
                if text in {"1", "true", "t", "yes", "y"}:
                    converted_bool.append(True)
                elif text in {"0", "false", "f", "no", "n"}:
                    converted_bool.append(False)
                else:
                    converted_bool.append(None)
            else:
                converted_bool.append(bool(value))
        return pa.array(converted_bool, type=pa.bool_())

    raise ValueError(f"Unsupported field kind: {spec.kind}")


def constant_array(value: Any, length: int, data_type: pa.DataType) -> pa.Array:
    return pa.array([value] * length, type=data_type)


def build_chunk_table(
    *,
    product: ProductSpec,
    schema: pa.Schema,
    data: Any,
    start: int,
    stop: int,
    night: str,
    path: Path,
    stat_result: os.stat_result,
    column_map: dict[str, str | None],
    source_column: str | None,
) -> pa.Table:
    length = stop - start
    arrays: list[pa.Array] = [
        constant_array(night, length, pa.string()),
        constant_array(product.name, length, pa.string()),
        constant_array(str(path), length, pa.string()),
        constant_array(path.name, length, pa.string()),
        constant_array(int(stat_result.st_size), length, pa.int64()),
        constant_array(iso_mtime(stat_result), length, pa.string()),
        pa.array(np.arange(start, stop, dtype=np.int64), type=pa.int64()),
    ]

    if source_column is None:
        arrays.append(constant_array(path.name, length, pa.string()))
    else:
        source_spec = FieldSpec("source_file", (source_column,), "string")
        arrays.append(arrow_array(data[source_column][start:stop], source_spec))

    for spec in product.fields:
        actual = column_map[spec.name]
        if actual is None:
            arrays.append(pa.nulls(length, type=arrow_type(spec.kind)))
        else:
            arrays.append(arrow_array(data[actual][start:stop], spec))

    return pa.Table.from_arrays(arrays, schema=schema)


def product_path(
    product: ProductSpec,
    *,
    night: str,
    processed_root: Path,
    heliolinc_root: Path,
    rr_subdir: str,
) -> Path:
    root = processed_root if product.root_kind == "processed" else heliolinc_root
    relative = product.relative_template.format(night=night, rr_subdir=rr_subdir)
    return root / relative


def empty_status(product: ProductSpec, night: str, path: Path) -> dict[str, Any]:
    return {
        "night": night,
        "product": product.name,
        "path": str(path),
        "exists": False,
        "status": "missing",
        "error_type": "",
        "error": "",
        "size_bytes": 0,
        "mtime_utc": "",
        "mtime_ns": 0,
        "sha256": "",
        "fits_checksum": "",
        "fits_datasum": "",
        "hdu": 1,
        "input_rows": 0,
        "rows_written": 0,
        "source_column": "",
        "column_map": {},
        "missing_required_columns": [],
        "missing_optional_columns": [],
    }


def extract_one_file(
    *,
    product: ProductSpec,
    night: str,
    path: Path,
    writer: StreamingParquetWriter,
    chunk_rows: int,
    calculate_sha256: bool,
) -> dict[str, Any]:
    status = empty_status(product, night, path)
    if not path.is_file():
        return status

    status["exists"] = True
    before = path.stat()
    status["size_bytes"] = int(before.st_size)
    status["mtime_utc"] = iso_mtime(before)
    status["mtime_ns"] = int(before.st_mtime_ns)

    try:
        with fits.open(
            path,
            mode="readonly",
            memmap=True,
            lazy_load_hdus=False,
        ) as hdul:
            if len(hdul) <= 1:
                raise ValueError(f"missing binary-table HDU 1; HDUs={len(hdul)}")
            hdu = hdul[1]
            status["fits_checksum"] = str(hdu.header.get("CHECKSUM", ""))
            status["fits_datasum"] = str(hdu.header.get("DATASUM", ""))
            data = hdu.data
            if data is None:
                names: list[str] = []
                input_rows = 0
            else:
                names = list(data.names or [])
                input_rows = len(data)
            status["input_rows"] = int(input_rows)

            column_map = resolve_column_map(names, product.fields)
            status["column_map"] = column_map
            missing_required = [
                spec.name
                for spec in product.fields
                if spec.required and column_map[spec.name] is None
            ]
            missing_optional = [
                spec.name
                for spec in product.fields
                if not spec.required and column_map[spec.name] is None
            ]
            status["missing_required_columns"] = missing_required
            status["missing_optional_columns"] = missing_optional
            source_column = resolve_first_column(names, product.source_aliases)
            status["source_column"] = source_column or ""

            if missing_required:
                status["status"] = "schema_error"
                status["error_type"] = "MissingRequiredColumns"
                status["error"] = ", ".join(missing_required)
                return status

            for start in range(0, input_rows, chunk_rows):
                stop = min(input_rows, start + chunk_rows)
                table = build_chunk_table(
                    product=product,
                    schema=writer.schema,
                    data=data,
                    start=start,
                    stop=stop,
                    night=night,
                    path=path,
                    stat_result=before,
                    column_map=column_map,
                    source_column=source_column,
                )
                writer.write(table)
                status["rows_written"] += table.num_rows

        if calculate_sha256:
            status["sha256"] = file_sha256(path)

        after = path.stat()
        if before.st_size != after.st_size or before.st_mtime_ns != after.st_mtime_ns:
            status["status"] = "changed_during_read"
            status["error_type"] = "InputChangedDuringRead"
            status["error"] = (
                f"before(size={before.st_size},mtime_ns={before.st_mtime_ns}); "
                f"after(size={after.st_size},mtime_ns={after.st_mtime_ns})"
            )
        elif status["missing_optional_columns"]:
            status["status"] = "ok_with_missing_optional"
        else:
            status["status"] = "ok"
    except Exception as exc:
        status["status"] = "error"
        status["error_type"] = type(exc).__name__
        status["error"] = str(exc)

    return status


STATUS_CSV_FIELDS = (
    "night",
    "product",
    "path",
    "exists",
    "status",
    "error_type",
    "error",
    "size_bytes",
    "mtime_utc",
    "mtime_ns",
    "sha256",
    "fits_checksum",
    "fits_datasum",
    "hdu",
    "input_rows",
    "rows_written",
    "source_column",
    "column_map",
    "missing_required_columns",
    "missing_optional_columns",
)


def atomic_write_text(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".inprogress")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def write_status_csv(path: Path, statuses: Sequence[dict[str, Any]]) -> None:
    temporary = path.with_name(path.name + ".inprogress")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for status in statuses:
            row = dict(status)
            for key in ("column_map", "missing_required_columns", "missing_optional_columns"):
                row[key] = json.dumps(row[key], ensure_ascii=False, sort_keys=True)
            writer.writerow(row)
    os.replace(temporary, path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Read-only extraction of frozen known/unknown/orbit FITS products "
            "into six compact Parquet tables."
        )
    )
    parser.add_argument("--start", type=parse_iso_date, default=parse_iso_date(DEFAULT_START))
    parser.add_argument("--end", type=parse_iso_date, default=parse_iso_date(DEFAULT_END))
    parser.add_argument("--processed-root", type=Path, default=Path(DEFAULT_PROCESSED_ROOT))
    parser.add_argument("--heliolinc-root", type=Path, default=Path(DEFAULT_HELIOLINC_ROOT))
    parser.add_argument(
        "--rr-subdir",
        default="rr_links",
        help="Nightly heliolinc linkage directory (default: rr_links).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Destination directory; must not be inside either production input root.",
    )
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=100_000,
        help="Maximum FITS rows per in-memory Arrow chunk/Parquet row group.",
    )
    parser.add_argument(
        "--compression",
        choices=("zstd", "snappy", "gzip", "none"),
        default="zstd",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=6,
        help="Used for zstd/gzip; ignored for snappy/none.",
    )
    parser.add_argument(
        "--skip-sha256",
        action="store_true",
        help="Skip whole-file SHA-256 hashing to avoid a second sequential read.",
    )
    parser.add_argument(
        "--fail-on-missing",
        action="store_true",
        help="Return nonzero if any expected nightly product is absent.",
    )
    return parser


def validate_args(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    if args.start > args.end:
        raise ValueError(f"start {args.start} is after end {args.end}")
    if args.chunk_rows <= 0:
        raise ValueError("--chunk-rows must be positive")
    if not args.rr_subdir or "/" in args.rr_subdir or args.rr_subdir in {".", ".."}:
        raise ValueError("--rr-subdir must be one plain directory name")

    processed_root = args.processed_root.expanduser().resolve(strict=False)
    heliolinc_root = args.heliolinc_root.expanduser().resolve(strict=False)
    output_dir = args.output_dir.expanduser().resolve(strict=False)

    protected_roots = {
        processed_root,
        heliolinc_root,
        *(root.resolve(strict=False) for root in DEFAULT_PROTECTED_ROOTS),
    }
    for input_root in sorted(protected_roots, key=str):
        if is_same_or_descendant(output_dir, input_root) or is_same_or_descendant(
            input_root, output_dir
        ):
            raise ValueError(
                f"refusing output that overlaps a protected production root: "
                f"output={output_dir}, protected_root={input_root}"
            )

    expected = [product.output_name for product in PRODUCTS]
    expected.extend(("file_status.json", "file_status.csv", "row_counts.json"))
    collisions = [
        output_dir / name
        for name in expected
        if (output_dir / name).exists()
        or (output_dir / f"{name}.inprogress").exists()
    ]
    if collisions:
        joined = "\n  ".join(str(path) for path in collisions)
        raise FileExistsError(
            "refusing to overwrite existing extraction outputs; use a new output directory:\n  "
            + joined
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    return processed_root, heliolinc_root, output_dir


def run(args: argparse.Namespace) -> int:
    processed_root, heliolinc_root, output_dir = validate_args(args)
    nights = list(iter_nights(args.start, args.end))
    compression = None if args.compression == "none" else args.compression
    compression_level = (
        args.compression_level if args.compression in {"zstd", "gzip"} else None
    )

    writers = {
        product.name: StreamingParquetWriter(
            output_path=output_dir / product.output_name,
            schema=product_schema(product),
            compression=compression,
            compression_level=compression_level,
        )
        for product in PRODUCTS
    }
    statuses: list[dict[str, Any]] = []
    started = datetime.now(timezone.utc)

    try:
        for night_index, night in enumerate(nights, start=1):
            for product in PRODUCTS:
                path = product_path(
                    product,
                    night=night,
                    processed_root=processed_root,
                    heliolinc_root=heliolinc_root,
                    rr_subdir=args.rr_subdir,
                )
                status = extract_one_file(
                    product=product,
                    night=night,
                    path=path,
                    writer=writers[product.name],
                    chunk_rows=args.chunk_rows,
                    calculate_sha256=not args.skip_sha256,
                )
                statuses.append(status)
            if night_index % 10 == 0 or night_index == len(nights):
                print(
                    f"[scan] nights={night_index}/{len(nights)} "
                    f"files={len(statuses)} rows="
                    + ",".join(
                        f"{name}:{writer.rows_written}"
                        for name, writer in writers.items()
                    ),
                    flush=True,
                )
    except BaseException:
        for writer in writers.values():
            writer.abort()
        raise

    for writer in writers.values():
        writer.finish()

    finished = datetime.now(timezone.utc)
    counts_by_product = {
        product.name: {
            "rows_written": writers[product.name].rows_written,
            "parquet_path": str(output_dir / product.output_name),
            "parquet_size_bytes": int((output_dir / product.output_name).stat().st_size),
            "parquet_sha256": file_sha256(output_dir / product.output_name),
            "files_ok": sum(
                status["status"] in {"ok", "ok_with_missing_optional"}
                for status in statuses
                if status["product"] == product.name
            ),
            "files_missing": sum(
                status["status"] == "missing"
                for status in statuses
                if status["product"] == product.name
            ),
            "files_error": sum(
                status["status"] in {"error", "schema_error", "changed_during_read"}
                for status in statuses
                if status["product"] == product.name
            ),
        }
        for product in PRODUCTS
    }
    row_counts = {
        "schema_version": SCHEMA_VERSION,
        "date_range_inclusive": {"start": args.start.isoformat(), "end": args.end.isoformat()},
        "night_count": len(nights),
        "expected_file_count": len(statuses),
        "present_file_count": sum(status["exists"] for status in statuses),
        "missing_file_count": sum(status["status"] == "missing" for status in statuses),
        "error_file_count": sum(
            status["status"] in {"error", "schema_error", "changed_during_read"}
            for status in statuses
        ),
        "products": counts_by_product,
    }
    status_payload = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": finished.isoformat(),
        "started_utc": started.isoformat(),
        "elapsed_seconds": (finished - started).total_seconds(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version.replace("\n", " "),
        "dependencies": {
            "numpy": np.__version__,
            "astropy": astropy_version,
            "pyarrow": pa.__version__,
        },
        "arguments": {
            "start": args.start.isoformat(),
            "end": args.end.isoformat(),
            "processed_root": str(processed_root),
            "heliolinc_root": str(heliolinc_root),
            "rr_subdir": args.rr_subdir,
            "output_dir": str(output_dir),
            "chunk_rows": args.chunk_rows,
            "compression": args.compression,
            "compression_level": compression_level,
            "sha256": not args.skip_sha256,
            "fail_on_missing": args.fail_on_missing,
        },
        "row_counts": row_counts,
        "files": statuses,
    }

    atomic_write_text(
        output_dir / "file_status.json",
        json.dumps(status_payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
    )
    write_status_csv(output_dir / "file_status.csv", statuses)
    atomic_write_text(
        output_dir / "row_counts.json",
        json.dumps(row_counts, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
    )

    errors = row_counts["error_file_count"]
    missing = row_counts["missing_file_count"]
    print(
        f"[done] output={output_dir} present={row_counts['present_file_count']} "
        f"missing={missing} errors={errors}",
        flush=True,
    )
    if errors:
        return 2
    if args.fail_on_missing and missing:
        return 3
    return 0


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        exit_code = run(args)
    except (FileExistsError, ValueError) as exc:
        parser.error(str(exc))
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
