#!/usr/bin/env python3
"""Assemble a machine-readable, hash-locked paper data snapshot.

The assembler reads only frozen products below a supplied project root.  It
does not consult production/live directories, contact the network, or invent
values for absent inputs.  Every component is labelled ``complete``,
``incomplete``, or ``blocked`` and every available input is SHA256 hashed.

Two generated files are written by default at the project root:

* ``paper_data_snapshot.json`` -- component inventory, recomputed metrics and
  cross-source QA checks;
* ``hashes.sha256`` -- a standard, sorted sha256sum-compatible manifest.

Existing outputs are refused.  ``--overwrite-generated`` may replace only
those two explicitly selected generated files; it never removes input data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_JSON = PROJECT_ROOT / "paper_data_snapshot.json"
DEFAULT_HASH_MANIFEST = PROJECT_ROOT / "hashes.sha256"
SCHEMA_VERSION = "1.0"


FIGURE_STEMS = (
    "fig01_system_architecture",
    "fig02_footprint_coverage",
    "fig03_data_accounting",
    "fig04_scheduler_example",
    "fig05_known_method_and_residuals",
    "fig06_unknown_pipeline_examples",
    "fig07_nightly_exposure_history",
    "fig08_known_results",
    "fig09_unknown_funnel",
    "fig10_high_confidence_distributions",
    "fig11_twilight_context",
    "fig12_operations_timeline",
)


IGNORED_NAMES = {".DS_Store"}
IGNORED_SUFFIXES = (".inprogress", ".tmp", ".pyc")


@dataclass(frozen=True)
class ComponentSpec:
    """A frozen analysis component and its minimum completeness contract."""

    name: str
    directory: Path
    required_files: tuple[str, ...] = ()
    require_any_file: bool = False
    description: str = ""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if hasattr(value, "item"):
        return json_safe(value.item())
    if isinstance(value, Mapping):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return str(value)


def path_label(path: Path, project_root: Path) -> str:
    resolved = path.resolve(strict=False)
    try:
        return resolved.relative_to(project_root.resolve()).as_posix()
    except ValueError:
        return resolved.as_posix()


def usable_file(path: Path, generated_outputs: set[Path]) -> bool:
    resolved = path.resolve(strict=False)
    if resolved in generated_outputs or not path.is_file():
        return False
    if path.name in IGNORED_NAMES or path.name.startswith("._"):
        return False
    if any(path.name.endswith(suffix) for suffix in IGNORED_SUFFIXES):
        return False
    if "__pycache__" in path.parts:
        return False
    return True


def collect_files(directory: Path, generated_outputs: set[Path]) -> list[Path]:
    if not directory.is_dir():
        return []
    return sorted(
        (path for path in directory.rglob("*") if usable_file(path, generated_outputs)),
        key=lambda path: path.as_posix(),
    )


def parquet_profile(path: Path) -> dict[str, Any]:
    try:
        import pyarrow.parquet as pq

        metadata = pq.ParquetFile(path).metadata
        schema = pq.ParquetFile(path).schema_arrow
        return {
            "row_count": int(metadata.num_rows),
            "column_count": int(len(schema.names)),
            "columns": list(schema.names),
            "profile_status": "ok",
        }
    except Exception as exc:  # pragma: no cover - environment-dependent backend errors
        return {
            "row_count": None,
            "column_count": None,
            "columns": None,
            "profile_status": "error",
            "profile_error": f"{type(exc).__name__}: {exc}",
        }


def csv_profile(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, [])
            row_count = sum(1 for _ in reader)
        return {
            "row_count": int(row_count),
            "column_count": int(len(header)),
            "columns": header,
            "profile_status": "ok",
        }
    except Exception as exc:
        return {
            "row_count": None,
            "column_count": None,
            "columns": None,
            "profile_status": "error",
            "profile_error": f"{type(exc).__name__}: {exc}",
        }


def json_profile(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            return {
                "json_type": "object",
                "top_level_keys": sorted(str(key) for key in payload),
                "profile_status": "ok",
            }
        if isinstance(payload, list):
            return {
                "json_type": "array",
                "item_count": len(payload),
                "profile_status": "ok",
            }
        return {"json_type": type(payload).__name__, "profile_status": "ok"}
    except Exception as exc:
        return {
            "json_type": None,
            "profile_status": "error",
            "profile_error": f"{type(exc).__name__}: {exc}",
        }


def file_profile(path: Path, project_root: Path, digest: str) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": path_label(path, project_root),
        "size_bytes": int(path.stat().st_size),
        "sha256": digest,
    }
    suffix = path.suffix.lower()
    if suffix in {".csv", ".tsv"}:
        result.update(csv_profile(path))
    elif suffix == ".parquet":
        result.update(parquet_profile(path))
    elif suffix == ".json":
        result.update(json_profile(path))
    return result


def audit_component(
    spec: ComponentSpec,
    *,
    project_root: Path,
    generated_outputs: set[Path],
    digest_cache: dict[Path, str],
) -> dict[str, Any]:
    directory = spec.directory.resolve(strict=False)
    files = collect_files(directory, generated_outputs)
    missing = [name for name in spec.required_files if not (directory / name).is_file()]
    issues: list[str] = []
    if not directory.is_dir():
        status = "blocked"
        issues.append("component directory does not exist")
    elif spec.require_any_file and not files:
        status = "blocked"
        issues.append("component directory contains no usable files")
    elif missing:
        status = "incomplete" if files else "blocked"
        issues.append("required files are missing")
    else:
        status = "complete"

    profiles: list[dict[str, Any]] = []
    for path in files:
        resolved = path.resolve()
        if resolved not in digest_cache:
            digest_cache[resolved] = sha256_file(resolved)
        profile = file_profile(resolved, project_root, digest_cache[resolved])
        profiles.append(profile)
        if profile.get("profile_status") == "error":
            status = "blocked"
            issues.append(f"cannot profile {profile['path']}")

    return {
        "status": status,
        "description": spec.description,
        "directory": path_label(directory, project_root),
        "required_files": list(spec.required_files),
        "missing_required_files": missing,
        "artifact_count": len(files),
        "artifacts": profiles,
        "issues": sorted(set(issues)),
    }


def read_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def read_csv(path: Path, *, usecols: Sequence[str] | None = None) -> pd.DataFrame | None:
    if not path.is_file():
        return None
    try:
        return pd.read_csv(path, usecols=usecols, low_memory=False)
    except (OSError, ValueError, pd.errors.ParserError):
        return None


def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.fillna(False).astype(bool)
    return (
        series.astype("string")
        .str.strip()
        .str.lower()
        .isin({"1", "true", "yes", "y", "t"})
    )


def metric(
    metrics: dict[str, dict[str, Any]],
    name: str,
    value: Any,
    *,
    source: Path,
    project_root: Path,
    grain: str,
    unit: str,
    method: str,
) -> None:
    metrics[name] = {
        "value": json_safe(value),
        "unit": unit,
        "grain": grain,
        "source": path_label(source, project_root),
        "method": method,
    }


def qa_check(
    checks: list[dict[str, Any]],
    name: str,
    left: Any,
    right: Any,
    *,
    left_source: Path,
    right_source: Path,
    project_root: Path,
    tolerance: float = 0.0,
) -> None:
    if left is None or right is None:
        status = "unavailable"
        difference = None
    else:
        try:
            difference = float(left) - float(right)
            status = "pass" if math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=tolerance) else "fail"
        except (TypeError, ValueError):
            difference = None
            status = "pass" if left == right else "fail"
    checks.append(
        {
            "name": name,
            "status": status,
            "left": json_safe(left),
            "right": json_safe(right),
            "difference_left_minus_right": json_safe(difference),
            "tolerance": tolerance,
            "left_source": path_label(left_source, project_root),
            "right_source": path_label(right_source, project_root),
        }
    )


def scalar_summary_metrics(
    payload: Mapping[str, Any] | None,
    *,
    prefixes: Sequence[str] = (),
) -> dict[str, Any]:
    """Return numeric/boolean result leaves without copying raw row arrays."""

    if payload is None:
        return {}
    selected: Mapping[str, Any] = payload
    for prefix in prefixes:
        value = selected.get(prefix)
        if not isinstance(value, Mapping):
            return {}
        selected = value

    output: dict[str, Any] = {}

    def visit(value: Any, key: str, depth: int) -> None:
        if isinstance(value, (bool, int, float)) and not isinstance(value, str):
            output[key] = json_safe(value)
        elif isinstance(value, Mapping) and depth < 3:
            for child_key, child in value.items():
                child_name = f"{key}.{child_key}" if key else str(child_key)
                visit(child, child_name, depth + 1)

    for key, value in selected.items():
        visit(value, str(key), 0)
    return output


def add_top_level_numeric_metrics(
    metrics: dict[str, dict[str, Any]],
    payload: Mapping[str, Any] | None,
    *,
    prefix: str,
    source: Path,
    project_root: Path,
) -> None:
    if payload is None:
        return
    for key, value in payload.items():
        if isinstance(value, bool) or (isinstance(value, (int, float)) and not isinstance(value, bool)):
            metric(
                metrics,
                f"{prefix}.{key}",
                value,
                source=source,
                project_root=project_root,
                grain="summary-defined",
                unit="source-defined",
                method="read from frozen summary JSON",
            )


def recompute_metrics(
    project_root: Path,
    directories: Mapping[str, Path],
    components: Mapping[str, dict[str, Any]],
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    metrics: dict[str, dict[str, Any]] = {}
    checks: list[dict[str, Any]] = []
    summary_scalars: dict[str, Any] = {}

    config_path = project_root / "config" / "snapshot.json"
    config = read_json(config_path)
    if config:
        for key in ("snapshot_label", "observation_start", "observation_end", "production_algorithm_commit", "repository_commit_reviewed"):
            if key in config:
                metric(
                    metrics,
                    f"config.{key}",
                    config[key],
                    source=config_path,
                    project_root=project_root,
                    grain="snapshot",
                    unit="identifier",
                    method="read from frozen configuration",
                )

    coverage_dir = directories["coverage"]
    coverage_csv = coverage_dir / "coverage_by_night.csv"
    fields_csv = coverage_dir / "field_visit_counts.csv"
    coverage_json = coverage_dir / "coverage_summary.json"
    coverage = read_csv(coverage_csv)
    fields = read_csv(fields_csv)
    coverage_summary = read_json(coverage_json)
    recomputed_coverage: dict[str, Any] = {}
    if coverage is not None and not coverage.empty:
        if "raw_exposure_n" in coverage:
            recomputed_coverage["raw_exposure_n"] = int(pd.to_numeric(coverage["raw_exposure_n"], errors="coerce").sum())
        if "night" in coverage:
            recomputed_coverage["raw_night_n"] = int(coverage["night"].astype("string").nunique())
        for column, key in (
            ("cumulative_open_shutter_hours", "open_shutter_hours"),
            ("cumulative_unique_field_n", "observed_field_n"),
            ("cumulative_healpix_pixel_n", "observed_pixel_n"),
            ("cumulative_area_deg2", "unique_area_deg2"),
        ):
            if column in coverage:
                values = pd.to_numeric(coverage[column], errors="coerce").dropna()
                if len(values):
                    recomputed_coverage[key] = float(values.iloc[-1])
    if fields is not None:
        if "field_id" in fields:
            recomputed_coverage["observed_field_n_from_fields"] = int(fields["field_id"].nunique())
        if "exposure_n" in fields:
            recomputed_coverage["raw_exposure_n_from_fields"] = int(pd.to_numeric(fields["exposure_n"], errors="coerce").sum())
        if "open_shutter_s" in fields:
            recomputed_coverage["open_shutter_hours_from_fields"] = float(pd.to_numeric(fields["open_shutter_s"], errors="coerce").sum() / 3600.0)
    for key, value in recomputed_coverage.items():
        source = fields_csv if key.endswith("_from_fields") else coverage_csv
        metric(
            metrics,
            f"coverage.{key}",
            value,
            source=source,
            project_root=project_root,
            grain="snapshot",
            unit="count_or_hours_or_deg2_as_named",
            method="recomputed from frozen CSV rows",
        )
    if coverage_summary:
        for key in ("raw_exposure_n", "raw_night_n", "observed_field_n", "observed_pixel_n", "open_shutter_hours", "unique_area_deg2"):
            qa_check(
                checks,
                f"coverage.{key}.summary_vs_recomputed",
                coverage_summary.get(key),
                recomputed_coverage.get(key),
                left_source=coverage_json,
                right_source=coverage_csv,
                project_root=project_root,
                tolerance=1e-9 if key in {"open_shutter_hours", "unique_area_deg2"} else 0.0,
            )

    inventory_dir = directories["inventory"]
    for filename, prefix in (("raw_manifest.csv", "raw"), ("l2_manifest.csv", "l2")):
        path = inventory_dir / filename
        frame = read_csv(path)
        if frame is None:
            continue
        metric(metrics, f"inventory.{prefix}_row_n", len(frame), source=path, project_root=project_root, grain="file row", unit="rows", method="count frozen CSV rows")
        if "night" in frame:
            metric(metrics, f"inventory.{prefix}_night_n", frame["night"].astype("string").nunique(), source=path, project_root=project_root, grain="night", unit="nights", method="count distinct night")
        if prefix == "raw" and "field_id" in frame:
            metric(metrics, "inventory.raw_field_n", frame["field_id"].nunique(), source=path, project_root=project_root, grain="field", unit="fields", method="count distinct field_id")
        if "exposure_s" in frame:
            metric(metrics, f"inventory.{prefix}_open_shutter_hours", pd.to_numeric(frame["exposure_s"], errors="coerce").sum() / 3600.0, source=path, project_root=project_root, grain="exposure", unit="hours", method="sum exposure_s / 3600")

    frozen_dir = directories["frozen_products"]
    frozen_counts_path = frozen_dir / "row_counts.json"
    frozen_counts = read_json(frozen_counts_path)
    if frozen_counts and isinstance(frozen_counts.get("products"), dict):
        for product, values in frozen_counts["products"].items():
            if not isinstance(values, dict) or "rows_written" not in values:
                continue
            metric(metrics, f"frozen_products.{product}.rows", values["rows_written"], source=frozen_counts_path, project_root=project_root, grain=f"{product} row", unit="rows", method="read frozen extraction row count")
            parquet = frozen_dir / f"{product}.parquet"
            if parquet.is_file():
                profile = parquet_profile(parquet)
                qa_check(checks, f"frozen_products.{product}.row_count_json_vs_parquet", values["rows_written"], profile.get("row_count"), left_source=frozen_counts_path, right_source=parquet, project_root=project_root)

    known_path = directories["derived_known"] / "known_population_summary.json"
    unknown_path = directories["derived_unknown"] / "unknown_population_summary.json"
    blinded_path = directories["blinded_proxy"] / "blinded_known_summary.json"
    known_summary = read_json(known_path)
    unknown_summary = read_json(unknown_path)
    blinded_summary = read_json(blinded_path)
    add_top_level_numeric_metrics(metrics, known_summary, prefix="known", source=known_path, project_root=project_root)
    add_top_level_numeric_metrics(metrics, unknown_summary, prefix="unknown", source=unknown_path, project_root=project_root)
    add_top_level_numeric_metrics(metrics, blinded_summary, prefix="blinded_known_proxy", source=blinded_path, project_root=project_root)

    unknown_dir = directories["derived_unknown"]
    link_path = unknown_dir / "unknown_all_links.csv"
    member_path = unknown_dir / "unknown_all_link_memberships.csv"
    unique_path = unknown_dir / "unknown_unique_detections.csv"
    links = read_csv(link_path)
    members = read_csv(member_path)
    unique = read_csv(unique_path)
    if links is not None and unknown_summary:
        qa_check(checks, "unknown.catalog_linkages.summary_vs_rows", unknown_summary.get("catalog_linkages"), len(links), left_source=unknown_path, right_source=link_path, project_root=project_root)
        for column, summary_key in (("initial_human_selected", "initial_human_selected_linkages"), ("posthoc_retained", "posthoc_retained_linkages")):
            if column in links:
                qa_check(checks, f"unknown.{summary_key}.summary_vs_rows", unknown_summary.get(summary_key), int(bool_series(links[column]).sum()), left_source=unknown_path, right_source=link_path, project_root=project_root)
    if members is not None and unknown_summary:
        qa_check(checks, "unknown.memberships.summary_vs_rows", unknown_summary.get("linkage_detection_memberships"), len(members), left_source=unknown_path, right_source=member_path, project_root=project_root)
    if unique is not None and unknown_summary:
        qa_check(checks, "unknown.unique_detections.summary_vs_rows", unknown_summary.get("globally_unique_detection_keys"), len(unique), left_source=unknown_path, right_source=unique_path, project_root=project_root)

    summary_sources: dict[str, tuple[Path, tuple[str, ...]]] = {
        "scheduler": (directories["scheduler"] / "scheduler_mode_summary.json", ()),
        "operations": (directories["operations"] / "operations_latency_summary.json", ()),
        "orbit_site": (directories["orbit_site"] / "orbit_site_sensitivity_summary.json", ()),
        "cross_night_primary": (unknown_dir / "unknown_linear_motion_candidate_group_summary.json", ("primary_results",)),
        "cross_night_pairs": (unknown_dir / "unknown_linear_motion_candidate_group_summary.json", ("pair_accounting",)),
        "jpl_first_pass": (directories["jpl"] / "jpl_identification_summary.json", ()),
        "jpl_second_pass": (directories["jpl"] / "second_pass" / "jpl_second_pass_summary.json", ()),
        "paper_tables": (directories["paper_tables"] / "paper_tables_summary.json", ()),
    }
    for name, (source, nested) in summary_sources.items():
        summary_scalars[name] = scalar_summary_metrics(read_json(source), prefixes=nested)

    jpl_dir = directories["jpl"]
    first_summary_path = jpl_dir / "jpl_identification_summary.json"
    first_summary = read_json(first_summary_path)
    first_status_path = jpl_dir / "jpl_query_status.csv"
    first_status = read_csv(first_status_path)
    if first_summary and first_status is not None:
        qa_check(checks, "jpl.first_pass.query_links_summary_vs_rows", first_summary.get("query_linkages"), len(first_status), left_source=first_summary_path, right_source=first_status_path, project_root=project_root)
        if "status" in first_status:
            successful = int(first_status["status"].astype("string").str.lower().eq("ok").sum())
            qa_check(checks, "jpl.first_pass.success_summary_vs_rows", first_summary.get("successful_queries"), successful, left_source=first_summary_path, right_source=first_status_path, project_root=project_root)
    second_summary_path = jpl_dir / "second_pass" / "jpl_second_pass_summary.json"
    second_summary = read_json(second_summary_path)
    confirmations_path = jpl_dir / "second_pass" / "jpl_second_pass_confirmations.csv"
    confirmations = read_csv(confirmations_path)
    if second_summary and confirmations is not None:
        qa_check(checks, "jpl.second_pass.rows_summary_vs_csv", second_summary.get("second_pass_result_rows"), len(confirmations), left_source=second_summary_path, right_source=confirmations_path, project_root=project_root)
        if "numerically_confirmed_candidate" in confirmations:
            qa_check(checks, "jpl.second_pass.confirmed_summary_vs_csv", second_summary.get("numerically_confirmed_linkages"), int(bool_series(confirmations["numerically_confirmed_candidate"]).sum()), left_source=second_summary_path, right_source=confirmations_path, project_root=project_root)

    table_summary_path = directories["paper_tables"] / "paper_tables_summary.json"
    table_summary = read_json(table_summary_path)
    if table_summary and isinstance(table_summary.get("outputs"), dict):
        for filename, values in table_summary["outputs"].items():
            if not isinstance(values, dict) or "row_count" not in values:
                continue
            table_path = directories["paper_tables"] / filename
            profile = csv_profile(table_path) if table_path.is_file() else {}
            qa_check(checks, f"paper_tables.{filename}.summary_vs_csv", values.get("row_count"), profile.get("row_count"), left_source=table_summary_path, right_source=table_path, project_root=project_root)

    figures = components["figures"]
    present_pairs = sum(
        1
        for stem in FIGURE_STEMS
        if (directories["figures"] / f"{stem}.pdf").is_file()
        and (directories["figures"] / f"{stem}.png").is_file()
    )
    metric(metrics, "figures.complete_pdf_png_pairs", present_pairs, source=directories["figures"], project_root=project_root, grain="planned figure", unit="figure pairs", method="count planned stems with both PDF and PNG")
    qa_check(checks, "figures.required_files_vs_present", len(figures["required_files"]), len(figures["required_files"]) - len(figures["missing_required_files"]), left_source=directories["figures"], right_source=directories["figures"], project_root=project_root)

    return metrics, checks, summary_scalars


def choose_directory(project_root: Path, explicit: Path | None, candidates: Sequence[str]) -> Path:
    if explicit is not None:
        return explicit.expanduser().resolve(strict=False)
    paths = [(project_root / candidate).resolve(strict=False) for candidate in candidates]
    existing = [path for path in paths if path.is_dir()]
    return existing[0] if existing else paths[0]


def build_directories(args: argparse.Namespace, project_root: Path) -> dict[str, Path]:
    return {
        "inventory": choose_directory(project_root, args.inventory_dir, ("snapshot/inventory",)),
        "inventory_partial": choose_directory(project_root, None, ("snapshot/inventory_partial",)),
        "tables": choose_directory(project_root, args.tables_dir, ("snapshot/tables",)),
        "coverage": choose_directory(project_root, args.coverage_dir, ("snapshot/coverage",)),
        "frozen_products": choose_directory(project_root, args.frozen_products_dir, ("snapshot/frozen_products",)),
        "derived_known": choose_directory(project_root, args.derived_known_dir, ("snapshot/derived_known", "snapshot/known_population")),
        "derived_unknown": choose_directory(project_root, args.derived_unknown_dir, ("snapshot/derived_unknown",)),
        "blinded_proxy": choose_directory(
            project_root,
            args.blinded_proxy_dir,
            ("snapshot/blinded_known", "snapshot/blinded_known_proxy", "snapshot/blinded_proxy"),
        ),
        "scheduler": choose_directory(project_root, args.scheduler_dir, ("snapshot/scheduler", "snapshot/scheduler_analysis")),
        "scheduler_plans": choose_directory(project_root, None, ("snapshot/plans",)),
        "operations": choose_directory(project_root, args.operations_dir, ("snapshot/operations", "snapshot/operations_analysis")),
        "operations_logs": choose_directory(project_root, None, ("snapshot/operations_logs",)),
        "orbit_site": choose_directory(
            project_root,
            args.orbit_site_dir,
            ("snapshot/orbit_site_comparison", "snapshot/orbit_site_sensitivity"),
        ),
        "jpl": choose_directory(project_root, args.jpl_dir, ("snapshot/jpl_identification",)),
        "paper_tables": choose_directory(project_root, args.paper_tables_dir, ("snapshot/paper_tables", "paper_tables", "tables")),
        "review_sample": choose_directory(project_root, None, ("snapshot/review_sample",)),
        "tracklet_stage": choose_directory(project_root, None, ("snapshot/tracklet_stage_counts",)),
        "random_shift": choose_directory(project_root, None, ("snapshot/random_shift",)),
        "provenance": choose_directory(project_root, None, ("snapshot/provenance",)),
        "figures": choose_directory(project_root, args.figures_dir, ("figures",)),
        "figure_data": choose_directory(project_root, args.figure_data_dir, ("figure_data",)),
    }


def component_specs(project_root: Path, directories: Mapping[str, Path]) -> list[ComponentSpec]:
    figure_required = tuple(
        filename
        for stem in FIGURE_STEMS
        for filename in (f"{stem}.pdf", f"{stem}.png")
    )
    return [
        ComponentSpec("config", project_root / "config", ("snapshot.json",), description="Frozen snapshot definition and executed configuration."),
        ComponentSpec("inventory", directories["inventory"], ("raw_manifest.csv", "raw_engineering_manifest.csv", "l2_manifest.csv", "nightly_products.csv", "mask_gaia_stage_counts.csv", "collector_summary.json"), description="Final server inventory; partial inventory is never substituted."),
        ComponentSpec("inventory_partial", directories["inventory_partial"], require_any_file=True, description="Explicitly partial collection retained for provenance only."),
        ComponentSpec("night_tables", directories["tables"], ("night_status.csv", "quality_mask.csv", "unknown_stage_counts_by_night.csv", "unknown_stage_definitions.csv", "unknown_stage_counts_total.json", "snapshot_summary.json"), description="Reconciled night and funnel accounting tables."),
        ComponentSpec("coverage", directories["coverage"], ("coverage_by_night.csv", "coverage_summary.json", "field_visit_counts.csv", "healpix_coverage.fits"), description="Frozen field/HEALPix coverage products."),
        ComponentSpec("frozen_products", directories["frozen_products"], ("row_counts.json", "known_all.parquet", "known_matched.parquet", "known_mask15.parquet", "orbit_links.parquet", "orbit_obs_residuals.parquet", "unknown_catalog.parquet", "unknown_review_detections.parquet"), description="Frozen production tables used by downstream analyses."),
        ComponentSpec("derived_known", directories["derived_known"], ("known_population_summary.json", "known_recovery_by_detection.parquet", "known_astrometric_residuals.parquet", "known_recovery_binned.csv", "known_recovery_by_night.csv", "known_astrometric_residuals_binned.csv"), description="Known-object recovery and astrometric analysis."),
        ComponentSpec("derived_unknown", directories["derived_unknown"], ("unknown_population_summary.json", "unknown_all_links.csv", "unknown_all_link_memberships.csv", "unknown_unique_detections.csv", "unknown_detection_multiplicity.csv", "unknown_high_confidence_links_recomputed.csv", "unknown_high_confidence_detections_recomputed.csv"), description="Unknown-object link/detection accounting."),
        ComponentSpec("blinded_known_proxy", directories["blinded_proxy"], ("blinded_known_summary.json", "blinded_known_link_classification.csv", "blinded_known_recovery_by_object.csv", "blinded_known_recovery_binned.csv"), description="Identity-blind known-object survival proxy."),
        ComponentSpec("scheduler", directories["scheduler"], ("scheduler_plan_realization_by_night.csv", "scheduler_exposure_matches.csv", "scheduler_mode_summary.csv", "scheduler_mode_summary.json"), description="Scheduler realization analysis."),
        ComponentSpec("scheduler_plans", directories["scheduler_plans"], require_any_file=True, description="Frozen scheduler plan inputs."),
        ComponentSpec("operations", directories["operations"], ("pipeline_latency_by_night.csv", "operations_event_evidence.csv", "operations_log_hashes.csv", "operations_logs.sha256", "operations_latency_summary.json"), description="Reconstructed operations latency analysis."),
        ComponentSpec("operations_source_logs", directories["operations_logs"], require_any_file=True, description="Frozen operations logs used by reconstruction."),
        ComponentSpec("orbit_site_sensitivity", directories["orbit_site"], ("orbit_site_sensitivity_by_link.csv", "orbit_site_sensitivity_by_night.csv", "orbit_site_sensitivity_summary.json", "hashes.sha256"), description="Orbit-confirmation observer-site sensitivity analysis."),
        ComponentSpec("cross_night", directories["derived_unknown"], ("review_and_mpc_status_with_linear_motion_candidate_groups.csv", "unknown_cross_night_pair_metrics.csv", "unknown_linear_motion_candidate_group_assignments.csv", "unknown_linear_motion_candidate_group_summary.json"), description="Linear-motion repeat-candidate screening; not orbit identification."),
        ComponentSpec("jpl_two_pass", directories["jpl"], ("jpl_identification_summary.json", "jpl_query_status.csv", "jpl_candidate_matches.csv", "jpl_best_match_by_link.csv", "mpc_reconciliation.csv", "second_pass/jpl_second_pass_summary.json", "second_pass/jpl_second_pass_confirmations.csv", "second_pass/jpl_second_pass_query_status.csv"), description="JPL first- and second-pass candidate association evidence."),
        ComponentSpec("paper_tables", directories["paper_tables"], ("table1_configuration_environment.csv", "table2_data_accounting.csv", "table3_known_recovery_astrometry.csv", "table4_unknown_funnel_retention.csv", "table5_retained_links.csv", "paper_tables_summary.json"), description="Publication-ready tables generated from frozen analysis inputs."),
        ComponentSpec("review_sample", directories["review_sample"], require_any_file=True, description="Frozen reviewed sample, status rows, and cutouts."),
        ComponentSpec("tracklet_stage_counts", directories["tracklet_stage"], require_any_file=True, description="Frozen tracklet-stage log accounting."),
        ComponentSpec("known_random_shift", directories["random_shift"], require_any_file=True, description="Frozen known-matcher random-shift control."),
        ComponentSpec("provenance", directories["provenance"], require_any_file=True, description="Production code/environment provenance."),
        ComponentSpec("figures", directories["figures"], figure_required, description="Twelve planned figures, each as vector PDF and 300-dpi PNG."),
        ComponentSpec("figure_data", directories["figure_data"], require_any_file=True, description="Figure-specific frozen plotting data; coverage is recorded by filename prefix."),
    ]


def write_atomic_text(path: Path, text: str, *, allow_replace: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not allow_replace:
        raise FileExistsError(f"refusing to overwrite generated output: {path}")
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, prefix=path.name + ".", suffix=".inprogress", delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def assemble(args: argparse.Namespace) -> dict[str, Any]:
    project_root = args.project_root.expanduser().resolve(strict=True)
    output_json = args.output_json.expanduser().resolve(strict=False) if args.output_json else project_root / "paper_data_snapshot.json"
    hash_manifest = args.hash_manifest.expanduser().resolve(strict=False) if args.hash_manifest else project_root / "hashes.sha256"
    if output_json == hash_manifest:
        raise ValueError("JSON output and hash manifest must be different files")
    existing = [path for path in (output_json, hash_manifest) if path.exists()]
    if existing and not args.overwrite_generated:
        raise FileExistsError("refusing to overwrite generated output(s): " + ", ".join(map(str, existing)))

    directories = build_directories(args, project_root)
    generated_outputs = {output_json.resolve(strict=False), hash_manifest.resolve(strict=False)}
    digest_cache: dict[Path, str] = {}
    components: dict[str, dict[str, Any]] = {}
    for spec in component_specs(project_root, directories):
        components[spec.name] = audit_component(spec, project_root=project_root, generated_outputs=generated_outputs, digest_cache=digest_cache)

    metrics, checks, summary_scalars = recompute_metrics(project_root, directories, components)
    failed_checks = [check["name"] for check in checks if check["status"] == "fail"]
    unavailable_checks = [check["name"] for check in checks if check["status"] == "unavailable"]
    blocked = [name for name, value in components.items() if value["status"] == "blocked"]
    incomplete = [name for name, value in components.items() if value["status"] == "incomplete"]
    if failed_checks or blocked:
        overall_status = "blocked"
    elif incomplete or unavailable_checks:
        overall_status = "incomplete"
    else:
        overall_status = "complete"

    figure_data_files = collect_files(directories["figure_data"], generated_outputs)
    figure_data_coverage = {
        f"fig{index:02d}": [path_label(path, project_root) for path in figure_data_files if path.name.startswith(f"fig{index:02d}_")]
        for index in range(1, 13)
    }
    config = read_json(project_root / "config" / "snapshot.json") or {}
    payload = {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": overall_status,
        "project_root": str(project_root),
        "snapshot_identity": {
            key: config.get(key)
            for key in ("snapshot_label", "observation_start", "observation_end", "production_algorithm_commit", "repository_commit_reviewed")
        },
        "assembly_policy": {
            "live_data_read": False,
            "network_access": False,
            "missing_value_policy": "never fill from expected values; mark the component incomplete or blocked",
            "partial_inventory_used_for_headlines": False,
            "hash_manifest_excludes_generated_outputs": [path_label(output_json, project_root), path_label(hash_manifest, project_root)],
            "overwrite_policy": "only the two explicitly selected generated outputs may be atomically replaced with --overwrite-generated",
        },
        "status_summary": {
            "blocked_components": blocked,
            "incomplete_components": incomplete,
            "failed_qa_checks": failed_checks,
            "unavailable_qa_checks": unavailable_checks,
        },
        "components": components,
        "headline_metrics": metrics,
        "component_summary_scalars": summary_scalars,
        "data_quality_checks": checks,
        "figure_data_coverage": figure_data_coverage,
        "hash_manifest": {
            "path": path_label(hash_manifest, project_root),
            "algorithm": "SHA256",
            "format": "sha256sum: <hex><two spaces><path>",
            "artifact_count": len(digest_cache),
        },
    }

    manifest_lines = []
    for path, digest in sorted(digest_cache.items(), key=lambda item: path_label(item[0], project_root)):
        label = path_label(path, project_root)
        if "\n" in label or "\r" in label:
            raise ValueError(f"cannot write newline-containing path to hash manifest: {label!r}")
        manifest_lines.append(f"{digest}  {label}\n")

    write_atomic_text(hash_manifest, "".join(manifest_lines), allow_replace=args.overwrite_generated)
    try:
        write_atomic_text(output_json, json.dumps(json_safe(payload), indent=2, ensure_ascii=False, sort_keys=True) + "\n", allow_replace=args.overwrite_generated)
    except Exception:
        # A first-run partial write must not leave only one of the two outputs.
        if not existing and hash_manifest.exists():
            hash_manifest.unlink()
        raise
    return payload


def add_directory_args(parser: argparse.ArgumentParser) -> None:
    for flag in (
        "inventory", "tables", "coverage", "frozen-products", "derived-known",
        "derived-unknown", "blinded-proxy", "scheduler", "operations",
        "orbit-site", "jpl", "paper-tables", "figures", "figure-data",
    ):
        parser.add_argument(f"--{flag}-dir", type=Path)


def write_csv_fixture(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def run_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="assemble-paper-snapshot-") as temporary:
        root = Path(temporary) / "paper"
        root.mkdir()
        config = {
            "snapshot_label": "fixture-snapshot",
            "observation_start": "2099-01-01",
            "observation_end": "2099-01-02",
            "production_algorithm_commit": "fixture-production-commit",
            "repository_commit_reviewed": "fixture-reviewed-commit",
        }
        (root / "config").mkdir()
        (root / "config" / "snapshot.json").write_text(json.dumps(config), encoding="utf-8")

        parser = build_parser()
        seed_args = parser.parse_args(["--project-root", str(root)])
        directories = build_directories(seed_args, root)
        specs = component_specs(root, directories)
        for spec in specs:
            spec.directory.mkdir(parents=True, exist_ok=True)
            for name in spec.required_files:
                path = spec.directory / name
                path.parent.mkdir(parents=True, exist_ok=True)
                if path.suffix == ".json":
                    path.write_text("{}\n", encoding="utf-8")
                elif path.suffix == ".csv":
                    write_csv_fixture(path, [{"fixture": 1}])
                elif path.suffix == ".parquet":
                    pd.DataFrame({"fixture": [1]}).to_parquet(path, index=False)
                else:
                    path.write_bytes(b"fixture\n")
            if spec.require_any_file and not any(spec.directory.rglob("*")):
                (spec.directory / "fixture.txt").write_text("fixture\n", encoding="utf-8")

        coverage = directories["coverage"]
        write_csv_fixture(coverage / "coverage_by_night.csv", [
            {"night": "20990101", "raw_exposure_n": 2, "cumulative_open_shutter_hours": 1.0, "cumulative_unique_field_n": 2, "cumulative_healpix_pixel_n": 10, "cumulative_area_deg2": 5.5},
            {"night": "20990102", "raw_exposure_n": 3, "cumulative_open_shutter_hours": 2.5, "cumulative_unique_field_n": 3, "cumulative_healpix_pixel_n": 20, "cumulative_area_deg2": 9.5},
        ])
        write_csv_fixture(coverage / "field_visit_counts.csv", [
            {"field_id": 1, "exposure_n": 2, "open_shutter_s": 3600},
            {"field_id": 2, "exposure_n": 2, "open_shutter_s": 3600},
            {"field_id": 3, "exposure_n": 1, "open_shutter_s": 1800},
        ])
        (coverage / "coverage_summary.json").write_text(json.dumps({"raw_exposure_n": 5, "raw_night_n": 2, "observed_field_n": 3, "observed_pixel_n": 20, "open_shutter_hours": 2.5, "unique_area_deg2": 9.5}), encoding="utf-8")

        inventory = directories["inventory"]
        write_csv_fixture(inventory / "raw_manifest.csv", [{"night": "20990101", "field_id": 1, "exposure_s": 30}, {"night": "20990102", "field_id": 2, "exposure_s": 40}])
        write_csv_fixture(inventory / "l2_manifest.csv", [{"night": "20990101", "exposure_s": 30}])

        frozen = directories["frozen_products"]
        products = ("known_all", "known_matched", "known_mask15", "orbit_links", "orbit_obs_residuals", "unknown_catalog", "unknown_review_detections")
        for product in products:
            pd.DataFrame({"fixture": [1]}).to_parquet(frozen / f"{product}.parquet", index=False)
        (frozen / "row_counts.json").write_text(json.dumps({"products": {product: {"rows_written": 1} for product in products}}), encoding="utf-8")

        unknown = directories["derived_unknown"]
        write_csv_fixture(unknown / "unknown_all_links.csv", [{"initial_human_selected": True, "posthoc_retained": True}, {"initial_human_selected": False, "posthoc_retained": False}])
        write_csv_fixture(unknown / "unknown_all_link_memberships.csv", [{"x": 1}, {"x": 2}, {"x": 3}])
        write_csv_fixture(unknown / "unknown_unique_detections.csv", [{"x": 1}, {"x": 2}])
        (unknown / "unknown_population_summary.json").write_text(json.dumps({"catalog_linkages": 2, "initial_human_selected_linkages": 1, "posthoc_retained_linkages": 1, "linkage_detection_memberships": 3, "globally_unique_detection_keys": 2}), encoding="utf-8")

        jpl = directories["jpl"]
        write_csv_fixture(jpl / "jpl_query_status.csv", [{"status": "ok"}, {"status": "ok"}])
        (jpl / "jpl_identification_summary.json").write_text(json.dumps({"query_linkages": 2, "successful_queries": 2}), encoding="utf-8")
        write_csv_fixture(jpl / "second_pass/jpl_second_pass_confirmations.csv", [{"numerically_confirmed_candidate": True}])
        (jpl / "second_pass/jpl_second_pass_summary.json").write_text(json.dumps({"second_pass_result_rows": 1, "numerically_confirmed_linkages": 1}), encoding="utf-8")

        paper_tables = directories["paper_tables"]
        outputs: dict[str, Any] = {}
        for filename in ("table1_configuration_environment.csv", "table2_data_accounting.csv", "table3_known_recovery_astrometry.csv", "table4_unknown_funnel_retention.csv", "table5_retained_links.csv"):
            write_csv_fixture(paper_tables / filename, [{"fixture": 1}])
            outputs[filename] = {"row_count": 1}
        (paper_tables / "paper_tables_summary.json").write_text(json.dumps({"outputs": outputs}), encoding="utf-8")

        payload = assemble(seed_args)
        assert payload["status"] == "complete", payload["status_summary"]
        assert payload["headline_metrics"]["coverage.raw_exposure_n"]["value"] == 5
        manifest = (root / "hashes.sha256").read_text(encoding="utf-8").splitlines()
        assert manifest and all(len(line.split("  ", 1)[0]) == 64 for line in manifest)
        manifest_labels = [line.split("  ", 1)[1] for line in manifest]
        assert "paper_data_snapshot.json" not in manifest_labels
        assert "hashes.sha256" not in manifest_labels

        try:
            assemble(seed_args)
        except FileExistsError:
            pass
        else:  # pragma: no cover - assertion path
            raise AssertionError("self-test: overwrite refusal did not trigger")

        overwrite_args = parser.parse_args(["--project-root", str(root), "--overwrite-generated"])
        overwritten = assemble(overwrite_args)
        assert overwritten["status"] == "complete"

        (directories["scheduler"] / "scheduler_mode_summary.json").unlink()
        incomplete_root = Path(temporary) / "incomplete.json"
        incomplete_hashes = Path(temporary) / "incomplete.sha256"
        incomplete_args = parser.parse_args(["--project-root", str(root), "--output-json", str(incomplete_root), "--hash-manifest", str(incomplete_hashes)])
        incomplete = assemble(incomplete_args)
        assert incomplete["status"] in {"incomplete", "blocked"}
        assert incomplete["components"]["scheduler"]["status"] == "incomplete"
    print("[self-test] assemble_paper_snapshot: PASS")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--hash-manifest", type=Path)
    parser.add_argument("--overwrite-generated", action="store_true", help="Atomically replace only the selected JSON and SHA256 outputs.")
    parser.add_argument("--self-test", action="store_true")
    add_directory_args(parser)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.self_test:
            run_self_test()
            return
        payload = assemble(args)
    except (FileExistsError, FileNotFoundError, ValueError, AssertionError) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": payload["status"], "status_summary": payload["status_summary"], "hash_artifact_count": payload["hash_manifest"]["artifact_count"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
