#!/usr/bin/env python3
"""Compare frozen orbit-confirmation products with an alternate-site rerun.

The comparison is keyed by ``(night, linkage_id)`` and never reads or writes
the live production tree.  It reports selection-status flips separately for
all orbit links, the formal unknown catalog, and the 58-link reviewed sample.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table


KEYS = ["night", "linkage_id"]
NUMERIC_COLUMNS = [
    "rms_arcsec",
    "med_arcsec",
    "max_arcsec",
    "a_au",
    "ecc",
    "inc_deg",
    "raan_deg",
    "argp_deg",
    "nu_deg",
    "pred_ra_deg",
    "pred_dec_deg",
    "lin_rms_arcsec",
    "best_v1_kms",
    "final_max_v_kms",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_unique(frame: pd.DataFrame, label: str) -> None:
    duplicate = frame.duplicated(KEYS, keep=False)
    if duplicate.any():
        sample = frame.loc[duplicate, KEYS].head(5).to_dict("records")
        raise ValueError(f"{label} has duplicate keys; examples={sample}")


def load_alternate(root: Path) -> tuple[pd.DataFrame, list[dict], list[Path]]:
    parts: list[pd.DataFrame] = []
    nights: list[dict] = []
    paths: list[Path] = []
    for path in sorted(root.glob("[0-9]" * 8 + "/orbit_links.fits")):
        night = path.parent.name
        table = Table.read(path, hdu=1)
        keep = [
            column
            for column in ["linkage_id", "fit_ok", "is_good", *NUMERIC_COLUMNS]
            if column in table.colnames
        ]
        part = table[keep].to_pandas()
        part.insert(0, "night", night)
        parts.append(part)
        paths.append(path)
        nights.append(
            {
                "night": night,
                "alternate_path": str(path),
                "alternate_sha256": sha256(path),
                "alternate_rows": int(len(part)),
            }
        )
    if not parts:
        raise FileNotFoundError(f"no alternate orbit_links.fits below {root}")
    frame = pd.concat(parts, ignore_index=True)
    frame["night"] = frame["night"].astype("string").str.zfill(8)
    frame["linkage_id"] = pd.to_numeric(frame["linkage_id"], errors="raise").astype("int64")
    ensure_unique(frame, "alternate rerun")
    return frame, nights, paths


def key_index(path: Path | None, label: str) -> pd.MultiIndex:
    if path is None:
        return pd.MultiIndex.from_arrays([[], []], names=KEYS)
    if path.suffix.lower() == ".parquet":
        frame = pd.read_parquet(path, columns=KEYS)
    else:
        frame = pd.read_csv(path, usecols=KEYS, dtype={"night": "string"})
    frame["night"] = frame["night"].astype("string").str.zfill(8)
    frame["linkage_id"] = pd.to_numeric(frame["linkage_id"], errors="raise").astype("int64")
    ensure_unique(frame, label)
    return pd.MultiIndex.from_frame(frame[KEYS])


def finite_summary(values: pd.Series) -> dict[str, float | int | None]:
    array = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    array = array[np.isfinite(array)]
    if not len(array):
        return {"n": 0, "median": None, "p16": None, "p84": None, "p95_abs": None, "max_abs": None}
    return {
        "n": int(len(array)),
        "median": float(np.median(array)),
        "p16": float(np.percentile(array, 16)),
        "p84": float(np.percentile(array, 84)),
        "p95_abs": float(np.percentile(np.abs(array), 95)),
        "max_abs": float(np.max(np.abs(array))),
    }


def status_counts(frame: pd.DataFrame, mask: pd.Series) -> dict[str, int]:
    selected = frame.loc[mask]
    return {
        "n": int(len(selected)),
        "matched_both_products_n": int(selected["present_both"].sum()),
        "fit_ok_baseline_n": int(selected["fit_ok_baseline"].fillna(False).sum()),
        "fit_ok_alternate_n": int(selected["fit_ok_alternate"].fillna(False).sum()),
        "fit_ok_flip_n": int(selected["fit_ok_flip"].fillna(False).sum()),
        "is_good_baseline_n": int(selected["is_good_baseline"].fillna(False).sum()),
        "is_good_alternate_n": int(selected["is_good_alternate"].fillna(False).sum()),
        "is_good_flip_n": int(selected["is_good_flip"].fillna(False).sum()),
        "missing_in_alternate_n": int(selected["missing_in_alternate"].sum()),
        "new_in_alternate_n": int(selected["missing_in_baseline"].sum()),
    }


def json_value(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-parquet", type=Path, required=True)
    parser.add_argument("--alternate-root", type=Path, required=True)
    parser.add_argument("--formal-unknown", type=Path)
    parser.add_argument("--high-confidence", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    baseline_path = args.baseline_parquet.resolve()
    alternate_root = args.alternate_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)

    baseline_columns = ["night", "linkage_id", "fit_ok", "is_good", *NUMERIC_COLUMNS]
    available = pd.read_parquet(baseline_path).columns
    baseline_columns = [column for column in baseline_columns if column in available]
    baseline = pd.read_parquet(baseline_path, columns=baseline_columns)
    baseline["night"] = baseline["night"].astype("string").str.zfill(8)
    baseline["linkage_id"] = pd.to_numeric(baseline["linkage_id"], errors="raise").astype("int64")
    ensure_unique(baseline, "baseline")

    alternate, alternate_nights, alternate_paths = load_alternate(alternate_root)
    comparison = baseline.merge(alternate, on=KEYS, how="outer", suffixes=("_baseline", "_alternate"), indicator=True)
    comparison["present_both"] = comparison["_merge"].eq("both")
    comparison["missing_in_baseline"] = comparison["_merge"].eq("right_only")
    comparison["missing_in_alternate"] = comparison["_merge"].eq("left_only")
    for field in ("fit_ok", "is_good"):
        left = comparison.get(f"{field}_baseline", pd.Series(pd.NA, index=comparison.index)).astype("boolean")
        right = comparison.get(f"{field}_alternate", pd.Series(pd.NA, index=comparison.index)).astype("boolean")
        comparison[f"{field}_flip"] = comparison["present_both"] & left.ne(right)
    for field in NUMERIC_COLUMNS:
        left_name = f"{field}_baseline"
        right_name = f"{field}_alternate"
        if left_name in comparison and right_name in comparison:
            comparison[f"delta_{field}"] = (
                pd.to_numeric(comparison[right_name], errors="coerce")
                - pd.to_numeric(comparison[left_name], errors="coerce")
            )

    index = pd.MultiIndex.from_frame(comparison[KEYS])
    formal = key_index(args.formal_unknown.resolve() if args.formal_unknown else None, "formal unknown")
    reviewed = key_index(args.high_confidence.resolve() if args.high_confidence else None, "high confidence")
    comparison["in_formal_unknown_catalog"] = index.isin(formal)
    comparison["in_high_confidence_58"] = index.isin(reviewed)
    comparison = comparison.sort_values(KEYS).reset_index(drop=True)

    night_rows: list[dict] = []
    for night, group in comparison.groupby("night", sort=True):
        row = {"night": night, **status_counts(group, pd.Series(True, index=group.index))}
        row["formal_unknown_n"] = int(group["in_formal_unknown_catalog"].sum())
        row["formal_fit_ok_flip_n"] = int((group["in_formal_unknown_catalog"] & group["fit_ok_flip"]).sum())
        row["formal_is_good_flip_n"] = int((group["in_formal_unknown_catalog"] & group["is_good_flip"]).sum())
        row["high_confidence_n"] = int(group["in_high_confidence_58"].sum())
        row["high_confidence_fit_ok_flip_n"] = int((group["in_high_confidence_58"] & group["fit_ok_flip"]).sum())
        row["high_confidence_is_good_flip_n"] = int((group["in_high_confidence_58"] & group["is_good_flip"]).sum())
        night_rows.append(row)
    night_frame = pd.DataFrame(night_rows)

    all_mask = pd.Series(True, index=comparison.index)
    formal_mask = comparison["in_formal_unknown_catalog"]
    reviewed_mask = comparison["in_high_confidence_58"]
    delta_summaries = {
        column.removeprefix("delta_"): finite_summary(comparison[column])
        for column in comparison.columns
        if column.startswith("delta_")
    }
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "comparison_key": KEYS,
        "baseline": {
            "path": str(baseline_path),
            "sha256": sha256(baseline_path),
            "rows": int(len(baseline)),
            "observer": {"longitude_deg_east": 117.575, "latitude_deg_north": 40.394239, "height_m": 868.221},
        },
        "alternate": {
            "root": str(alternate_root),
            "nights": int(alternate["night"].nunique()),
            "rows": int(len(alternate)),
        },
        "all_orbit_links": status_counts(comparison, all_mask),
        "formal_unknown_catalog": status_counts(comparison, formal_mask),
        "high_confidence_58": status_counts(comparison, reviewed_mask),
        "numeric_deltas_alternate_minus_baseline": delta_summaries,
        "interpretation_guardrail": (
            "This is a deterministic observer-location sensitivity comparison, not a new MPC submission. "
            "Canonical surveyed MPC 327 coordinates still require author/observatory confirmation."
        ),
    }
    run_summary = alternate_root / "run_summary.json"
    if run_summary.exists():
        summary["alternate_run_summary"] = json.loads(run_summary.read_text(encoding="utf-8"))

    comparison.to_csv(output / "orbit_site_sensitivity_by_link.csv", index=False)
    night_frame.to_csv(output / "orbit_site_sensitivity_by_night.csv", index=False)
    (output / "orbit_site_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=json_value) + "\n", encoding="utf-8"
    )
    inputs = [baseline_path, *alternate_paths]
    for optional in (args.formal_unknown, args.high_confidence):
        if optional is not None:
            inputs.append(optional.resolve())
    with (output / "hashes.sha256").open("w", encoding="utf-8") as handle:
        for path in inputs:
            handle.write(f"{sha256(path)}  {path}\n")
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_value))


if __name__ == "__main__":
    main()
