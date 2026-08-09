#!/usr/bin/env python3
"""Build publication tables from explicitly supplied, frozen analysis products.

The script never discovers or modifies production data and never contacts
external services.  Every value is copied from, or arithmetically checked
against, an explicitly provided JSON/CSV input.  Missing optional evidence
becomes a visible ``pending``/``not_available`` state rather than an inferred
value.

JPL candidate associations and authoritative MPC submission/designation states
are kept in separate column families.  A JPL first- or second-pass association
can never populate an MPC ingest state, identification, or designation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd


SCHEMA_VERSION = "1.0"
DEFAULT_EXPECTED_RETAINED = 58

TABLE_FILENAMES = {
    "table1": "table1_configuration_environment.csv",
    "table2": "table2_data_accounting.csv",
    "table3": "table3_known_recovery_astrometry.csv",
    "table4": "table4_unknown_funnel_retention.csv",
    "table5": "table5_retained_links.csv",
}
SUMMARY_FILENAME = "paper_tables_summary.json"

TABLE1_COLUMNS = (
    "row_order",
    "section",
    "parameter",
    "value",
    "unit",
    "evidence_status",
    "definition",
    "source_input",
    "input_manifest_hash",
)
TABLE2_COLUMNS = (
    "row_order",
    "metric",
    "value",
    "unit",
    "grain",
    "evidence_status",
    "definition",
    "source_input",
    "input_manifest_hash",
)
TABLE3_COLUMNS = (
    "row_order",
    "section",
    "metric",
    "value",
    "unit",
    "denominator",
    "evidence_status",
    "definition",
    "caveat",
    "source_input",
    "input_manifest_hash",
)
TABLE4_COLUMNS = (
    "stage_order",
    "stage",
    "value",
    "unit",
    "scope",
    "evidence_status",
    "denominator_stage",
    "retention_fraction",
    "definition",
    "interpretation_guardrail",
    "source_input",
    "input_manifest_hash",
)
TABLE5_COLUMNS = (
    "row_order",
    "link_key",
    "origin_night",
    "trk_sub",
    "linkage_id",
    "n_obs",
    "n_tracklets",
    "unique_detection_n",
    "fit_ok",
    "is_good",
    "rms_arcsec",
    "median_residual_arcsec",
    "max_arcsec",
    "a_au",
    "ecc",
    "inc_deg",
    "raan_deg",
    "argp_deg",
    "nu_deg",
    "best_v1_kms",
    "linear_rms_arcsec",
    "speed_arcsec_per_hour",
    "direction_deg",
    "first_mjd",
    "median_mjd",
    "last_mjd",
    "median_ra_deg",
    "median_dec_deg",
    "median_mag_aper4",
    "time_span_hours",
    "ecliptic_lon_j2000_deg",
    "ecliptic_lat_j2000_deg",
    "solar_elongation_deg",
    "nearest_astronomical_twilight",
    "nearest_twilight_abs_hours",
    "nearest_twilight_signed_hours",
    "final_paper_status",
    "review_notes",
    "source_gif",
    "linear_motion_candidate_group_id",
    "candidate_group_link_count",
    "candidate_group_distinct_night_count",
    "candidate_group_cross_night_primary_edge_count",
    "candidate_group_has_cross_night_primary_edge",
    "candidate_group_is_singleton",
    "candidate_group_interpretation",
    "jpl_first_pass_status",
    "jpl_first_pass_object_name",
    "jpl_first_pass_separation_arcsec",
    "jpl_first_pass_strict_candidate",
    "jpl_first_pass_plausible_candidate",
    "jpl_second_pass_status",
    "jpl_second_pass_object_name",
    "jpl_second_pass_separation_arcsec",
    "jpl_second_pass_numerically_confirmed_candidate",
    "jpl_evidence_scope",
    "mpc_evidence_status",
    "mpc_evidence_source",
    "mpc_submission_id",
    "mpc_ingest_state",
    "mpc_identification",
    "provisional_designation",
    "known_object_id",
    "identity_guardrail",
    "retained_link_source_path",
    "input_manifest_hash",
)

JPL_GUARDRAIL = (
    "JPL positional/rate association only; not evidence of MPC ingest, MPC "
    "identification, provisional designation, or discovery confirmation."
)
MPC_PENDING = "pending_no_authoritative_mpc_evidence"
RETAINED_GUARDRAIL = (
    "One retained single-night linkage row; not necessarily one independent object."
)
UNKNOWN_RETENTION_GUARDRAIL = "Stage retention/post-audit precision proxy only; not a completeness or discovery rate."


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def combined_hash(hashes: dict[str, str]) -> str:
    text = "".join(f"{name}:{digest}\n" for name, digest in sorted(hashes.items()))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def normalize_text(value: Any) -> str:
    if value is None or value is pd.NA:
        return ""
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace").strip()
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "null", "<na>"} else text


def normalize_night(value: Any) -> str:
    digits = re.sub(r"\D", "", normalize_text(value))
    return digits[:8] if len(digits) >= 8 else ""


def normalize_linkage_id(value: Any) -> str:
    text = normalize_text(value)
    if re.fullmatch(r"[-+]?\d+(?:\.0+)?", text):
        return str(int(float(text)))
    return text


def truthy(value: Any) -> bool:
    return normalize_text(value).lower() in {"1", "true", "t", "yes", "y"}


def numeric(value: Any) -> float | None:
    try:
        result = float(normalize_text(value))
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def display_value(value: Any) -> Any:
    if value is None or value is pd.NA:
        return ""
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.floating, float)):
        if not math.isfinite(float(value)):
            return ""
        if math.isclose(float(value), round(float(value)), abs_tol=1.0e-12):
            return int(round(float(value)))
        return float(value)
    if isinstance(value, bool):
        return value
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, ensure_ascii=False)
    return normalize_text(value)


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: JSON root must be an object")
    return payload


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False)


def require_columns(frame: pd.DataFrame, required: Iterable[str], path: Path) -> None:
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{path}: missing required columns {missing}")


def deep_get(payload: dict[str, Any] | None, dotted: str) -> Any:
    if payload is None:
        return None
    current: Any = payload
    for part in dotted.split("."):
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    return current


def record_check(
    checks: list[dict[str, Any]],
    *,
    check_id: str,
    status: str,
    detail: str,
    severity: str = "critical",
) -> None:
    checks.append(
        {
            "check_id": check_id,
            "status": status,
            "severity": severity,
            "detail": detail,
        }
    )


def assert_close(
    checks: list[dict[str, Any]],
    check_id: str,
    left: Any,
    right: Any,
    *,
    left_label: str,
    right_label: str,
    tolerance: float = 1.0e-9,
) -> None:
    left_value = numeric(left)
    right_value = numeric(right)
    if left_value is None or right_value is None:
        record_check(
            checks,
            check_id=check_id,
            status="not_run_missing_value",
            severity="medium",
            detail=f"{left_label}={left!r}; {right_label}={right!r}",
        )
        return
    if not math.isclose(left_value, right_value, rel_tol=tolerance, abs_tol=tolerance):
        raise ValueError(
            f"data-quality check {check_id} failed: "
            f"{left_label}={left_value}, {right_label}={right_value}"
        )
    record_check(
        checks,
        check_id=check_id,
        status="pass",
        detail=f"{left_label}={left_value:g}; {right_label}={right_value:g}",
    )


def table1_configuration_environment(
    config: dict[str, Any],
    environment: dict[str, Any] | None,
    production_hashes: pd.DataFrame | None,
    orbit_site_summary: dict[str, Any] | None,
    *,
    input_hash: str,
    checks: list[dict[str, Any]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    def add(
        section: str,
        parameter: str,
        value: Any,
        unit: str,
        definition: str,
        source: str,
        status: str | None = None,
    ) -> None:
        text = display_value(value)
        available = text != ""
        rows.append(
            {
                "row_order": len(rows) + 1,
                "section": section,
                "parameter": parameter,
                "value": text,
                "unit": unit,
                "evidence_status": status
                or ("available" if available else "not_available"),
                "definition": definition,
                "source_input": source,
                "input_manifest_hash": input_hash,
            }
        )

    config_specs = (
        ("facility", "facility_name", "facility.name", "", "Facility name"),
        (
            "facility",
            "observatory_name",
            "facility.observatory_name",
            "",
            "Observatory name",
        ),
        (
            "facility",
            "telescope_aperture",
            "facility.telescope_aperture_m",
            "m",
            "Telescope aperture",
        ),
        ("facility", "detector_name", "facility.detector_name", "", "Science detector"),
        (
            "facility",
            "pixel_scale",
            "facility.pixel_scale_arcsec",
            "arcsec pixel^-1",
            "Detector pixel scale",
        ),
        (
            "facility",
            "field_of_view",
            "facility.field_of_view_deg2",
            "deg2",
            "Single-pointing field of view",
        ),
        (
            "facility",
            "unfiltered_response",
            "facility.unfiltered_response",
            "",
            "Unfiltered effective response or calibration convention",
        ),
        (
            "facility",
            "upstream_pipeline_version",
            "upstream.pipeline_version",
            "",
            "Upstream image-processing version",
        ),
        ("snapshot", "snapshot_label", "snapshot_label", "", "Frozen analysis label"),
        (
            "snapshot",
            "observation_start",
            "observation_start",
            "UTC date",
            "Inclusive start",
        ),
        ("snapshot", "observation_end", "observation_end", "UTC date", "Inclusive end"),
        ("snapshot", "pipeline_timezone", "timezone", "", "Operational timezone"),
        (
            "version",
            "production_algorithm_commit",
            "production_algorithm_commit",
            "git SHA",
            "Frozen algorithm revision",
        ),
        (
            "version",
            "repository_commit_reviewed",
            "repository_commit_reviewed",
            "git SHA",
            "Repository revision reviewed",
        ),
        (
            "quality",
            "quality_rules_status",
            "quality_rules.status",
            "",
            "Author sign-off state",
        ),
        (
            "site",
            "scheduler_longitude",
            "site_configuration_as_executed.scheduler.longitude_deg_east",
            "deg E",
            "Scheduler site longitude",
        ),
        (
            "site",
            "scheduler_latitude",
            "site_configuration_as_executed.scheduler.latitude_deg_north",
            "deg N",
            "Scheduler site latitude",
        ),
        (
            "site",
            "scheduler_height",
            "site_configuration_as_executed.scheduler.height_m",
            "m",
            "Scheduler site height",
        ),
        (
            "site",
            "known_matcher_longitude",
            "site_configuration_as_executed.known_matcher.longitude_deg_east",
            "deg E",
            "Known matcher site longitude",
        ),
        (
            "site",
            "known_matcher_latitude",
            "site_configuration_as_executed.known_matcher.latitude_deg_north",
            "deg N",
            "Known matcher site latitude",
        ),
        (
            "site",
            "known_matcher_height",
            "site_configuration_as_executed.known_matcher.height_m",
            "m",
            "Known matcher site height",
        ),
        (
            "site",
            "single_night_orbit_longitude",
            "site_configuration_as_executed.single_night_orbit_confirm.longitude_deg_east",
            "deg E",
            "Single-night orbit-confirmation longitude",
        ),
        (
            "site",
            "single_night_orbit_latitude",
            "site_configuration_as_executed.single_night_orbit_confirm.latitude_deg_north",
            "deg N",
            "Single-night orbit-confirmation latitude",
        ),
        (
            "site",
            "single_night_orbit_height",
            "site_configuration_as_executed.single_night_orbit_confirm.height_m",
            "m",
            "Single-night orbit-confirmation height",
        ),
        (
            "site",
            "multi_night_orbit_longitude",
            "site_configuration_as_executed.multi_night_orbit_confirm.longitude_deg_east",
            "deg E",
            "Multi-night orbit-confirmation longitude",
        ),
        (
            "site",
            "multi_night_orbit_latitude",
            "site_configuration_as_executed.multi_night_orbit_confirm.latitude_deg_north",
            "deg N",
            "Multi-night orbit-confirmation latitude",
        ),
        (
            "site",
            "multi_night_orbit_height",
            "site_configuration_as_executed.multi_night_orbit_confirm.height_m",
            "m",
            "Multi-night orbit-confirmation height",
        ),
        (
            "site",
            "legacy_twilight_longitude",
            "site_configuration_as_executed.legacy_twilight_analysis_not_final.longitude_deg_east",
            "deg E",
            "Legacy twilight-analysis longitude, not final",
        ),
        (
            "site",
            "legacy_twilight_latitude",
            "site_configuration_as_executed.legacy_twilight_analysis_not_final.latitude_deg_north",
            "deg N",
            "Legacy twilight-analysis latitude, not final",
        ),
        (
            "site",
            "legacy_twilight_height",
            "site_configuration_as_executed.legacy_twilight_analysis_not_final.height_m",
            "m",
            "Legacy twilight-analysis height, not final",
        ),
        (
            "site",
            "canonical_site_status",
            "site_configuration_as_executed.final_canonical_site_status",
            "",
            "Final surveyed-site confirmation state",
        ),
        (
            "site_metadata_warning",
            "upstream_header_observatory",
            "upstream_header_site_metadata.observatory",
            "",
            "Upstream header observatory label",
        ),
        (
            "site_metadata_warning",
            "upstream_header_longitude",
            "upstream_header_site_metadata.longitude_deg_east",
            "deg E",
            "Upstream header longitude; use only according to its status",
        ),
        (
            "site_metadata_warning",
            "upstream_header_latitude",
            "upstream_header_site_metadata.latitude_deg_north",
            "deg N",
            "Upstream header latitude; use only according to its status",
        ),
        (
            "site_metadata_warning",
            "upstream_header_height",
            "upstream_header_site_metadata.height_m",
            "m",
            "Upstream header height; use only according to its status",
        ),
        (
            "site_metadata_warning",
            "upstream_header_site_status",
            "upstream_header_site_metadata.status",
            "",
            "Validity of upstream header site metadata",
        ),
        (
            "known",
            "predicted_v_limit",
            "known.predicted_v_limit_mag",
            "mag",
            "Nominal known-object prediction limit",
        ),
        (
            "known",
            "official_match_radius",
            "known.official_match_radius_arcsec",
            "arcsec",
            "Official known-object match radius",
        ),
        (
            "known",
            "unknown_mask_radius",
            "known.unknown_mask_radius_arcsec",
            "arcsec",
            "Known subtraction/mask radius",
        ),
        (
            "unknown",
            "gaia_match_radius",
            "unknown.gaia_match_radius_arcsec",
            "arcsec",
            "Gaia static-source match radius",
        ),
        (
            "unknown",
            "minimum_speed",
            "unknown.minimum_speed_arcsec_per_hour",
            "arcsec h^-1",
            "Tracklet speed lower bound",
        ),
        (
            "unknown",
            "maximum_speed",
            "unknown.maximum_speed_arcsec_per_hour",
            "arcsec h^-1",
            "Tracklet speed upper bound",
        ),
        (
            "unknown",
            "maximum_magnitude_difference",
            "unknown.maximum_magnitude_difference",
            "mag",
            "Tracklet magnitude consistency",
        ),
        (
            "unknown",
            "maximum_speed_difference",
            "unknown.maximum_speed_difference_arcsec_per_hour",
            "arcsec h^-1",
            "Link speed consistency",
        ),
        (
            "unknown",
            "maximum_direction_difference",
            "unknown.maximum_direction_difference_deg",
            "deg",
            "Link direction consistency",
        ),
        (
            "unknown",
            "skip_common_area_as_executed",
            "unknown.production_skip_common_area",
            "boolean",
            "Production common-area option",
        ),
        (
            "unknown",
            "catalog_selector_as_executed",
            "unknown.single_night_catalog_selector_as_executed",
            "",
            "Single-night catalog selection",
        ),
        (
            "unknown",
            "ADES_time_as_executed",
            "unknown.unknown_ades_time_as_executed",
            "",
            "Unknown ADES epoch convention",
        ),
        (
            "unknown",
            "orbit_tracklet_time_reference",
            "unknown.orbit_and_tracklet_time_reference",
            "",
            "Orbit/tracklet epoch convention",
        ),
        (
            "unknown",
            "standard_exposure_midpoint_offset",
            "unknown.standard_exposure_midpoint_offset_s",
            "s",
            "Standard start-to-midpoint offset",
        ),
        (
            "unknown",
            "quarantine_threshold",
            "unknown.quarantine_threshold_links",
            "linkages night^-1",
            "Unknown-night quarantine threshold",
        ),
        (
            "unknown_orbit_quality",
            "minimum_tracklets",
            "unknown.single_night_is_good_requires.minimum_tracklets",
            "tracklets",
            "Single-night is_good minimum tracklets",
        ),
        (
            "unknown_orbit_quality",
            "minimum_used_observations",
            "unknown.single_night_is_good_requires.minimum_used_observations",
            "observations",
            "Single-night is_good minimum used observations",
        ),
        (
            "unknown_orbit_quality",
            "minimum_used_fraction",
            "unknown.single_night_is_good_requires.minimum_used_fraction",
            "fraction",
            "Single-night is_good minimum used fraction",
        ),
        (
            "unknown_orbit_quality",
            "rms_limit",
            "unknown.single_night_is_good_requires.rms_arcsec_lte",
            "arcsec",
            "Single-night is_good RMS limit",
        ),
        (
            "unknown_orbit_quality",
            "maximum_residual_limit",
            "unknown.single_night_is_good_requires.max_residual_arcsec_lte",
            "arcsec",
            "Single-night is_good maximum residual limit",
        ),
        (
            "unknown_orbit_quality",
            "linear_rms_limit",
            "unknown.single_night_is_good_requires.linear_rms_arcsec_lte",
            "arcsec",
            "Single-night is_good linear-motion RMS limit",
        ),
    )
    for section, parameter, dotted, unit, definition in config_specs:
        add(section, parameter, deep_get(config, dotted), unit, definition, "config")

    if environment is None:
        for parameter, unit, definition in (
            ("host", "", "Production host"),
            ("platform", "", "Production platform"),
            ("python_version", "", "Default production Python"),
            ("short_arc_python_version", "", "Short-arc environment Python"),
            ("slurm_accounting_status", "", "Resource-accounting availability"),
        ):
            add("environment", parameter, None, unit, definition, "environment")
    else:
        add(
            "environment",
            "host",
            environment.get("host"),
            "",
            "Production host",
            "environment",
        )
        add(
            "environment",
            "platform",
            environment.get("platform"),
            "",
            "Production platform",
            "environment",
        )
        add(
            "environment",
            "python_executable",
            environment.get("python_executable"),
            "",
            "Default interpreter",
            "environment",
        )
        add(
            "environment",
            "python_version",
            environment.get("python_version"),
            "",
            "Default Python",
            "environment",
        )
        add(
            "model",
            "known_propagation",
            deep_get(environment, "model_notes.known_propagation"),
            "",
            "Known-object propagation model",
            "environment",
        )
        add(
            "model",
            "short_arc_orbit",
            deep_get(environment, "model_notes.short_arc_orbit"),
            "",
            "Short-arc orbit model and ephemeris",
            "environment",
        )
        for package, version in sorted((environment.get("packages") or {}).items()):
            add(
                "environment",
                f"default_package_{package}",
                version,
                "version",
                "Frozen package version",
                "environment",
            )
        short_arc = environment.get("short_arc_environment") or {}
        short_payload: dict[str, Any] = {}
        try:
            short_payload = json.loads(normalize_text(short_arc.get("stdout")))
            if not isinstance(short_payload, dict):
                short_payload = {}
        except (AttributeError, json.JSONDecodeError):
            short_payload = {}
        add(
            "environment",
            "short_arc_python_executable",
            short_payload.get("executable"),
            "",
            "Short-arc interpreter",
            "environment",
        )
        add(
            "environment",
            "short_arc_python_version",
            short_payload.get("python"),
            "",
            "Short-arc Python",
            "environment",
        )
        for package, version in sorted((short_payload.get("packages") or {}).items()):
            add(
                "environment",
                f"short_arc_package_{package}",
                version,
                "version",
                "Short-arc package version",
                "environment",
            )
        accounting = environment.get("slurm_accounting_probe") or {}
        accounting_status = (
            "unavailable_slurm_accounting_disabled"
            if numeric(accounting.get("returncode")) not in {None, 0.0}
            else "available"
        )
        add(
            "environment",
            "slurm_accounting_status",
            accounting_status,
            "",
            "Resource-accounting availability",
            "environment",
        )
        add(
            "reference",
            "footprint_rows",
            deep_get(environment, "footprint.n_rows"),
            "fields",
            "Footprint library rows",
            "environment",
        )
        add(
            "reference",
            "gaia_tile_count",
            deep_get(environment, "gaia_tiles.count"),
            "tiles",
            "Frozen Gaia HEALPix tile count",
            "environment",
        )
        add(
            "reference",
            "gaia_release_provenance",
            deep_get(environment, "gaia_tiles.release_provenance_in_files"),
            "",
            "Gaia release provenance status",
            "environment",
        )

    if production_hashes is None:
        add(
            "provenance",
            "production_files_hashed",
            None,
            "files",
            "Production file hash count",
            "production_hashes",
        )
    else:
        require_columns(
            production_hashes,
            {"category", "hash_status", "sha256", "path"},
            Path("production_hashes"),
        )
        ok = production_hashes["hash_status"].map(normalize_text).eq("ok")
        add(
            "provenance",
            "production_files_listed",
            len(production_hashes),
            "files",
            "Frozen production provenance rows",
            "production_hashes",
        )
        add(
            "provenance",
            "production_files_hashed_ok",
            int(ok.sum()),
            "files",
            "Rows with a recorded SHA-256",
            "production_hashes",
        )
        add(
            "provenance",
            "production_code_files_hashed_ok",
            int((ok & production_hashes["category"].eq("production_code")).sum()),
            "files",
            "Production-code rows with SHA-256",
            "production_hashes",
        )

    orbit_specs = (
        (
            "baseline_longitude",
            "baseline.observer.longitude_deg_east",
            "deg E",
            "As-executed orbit baseline longitude",
        ),
        (
            "baseline_latitude",
            "baseline.observer.latitude_deg_north",
            "deg N",
            "As-executed orbit baseline latitude",
        ),
        (
            "baseline_height",
            "baseline.observer.height_m",
            "m",
            "As-executed orbit baseline height",
        ),
        (
            "alternate_960m_longitude",
            "alternate_run_summary.observer.longitude_deg_east",
            "deg E",
            "960 m sensitivity-run longitude",
        ),
        (
            "alternate_960m_latitude",
            "alternate_run_summary.observer.latitude_deg_north",
            "deg N",
            "960 m sensitivity-run latitude",
        ),
        (
            "alternate_960m_height",
            "alternate_run_summary.observer.height_m",
            "m",
            "960 m sensitivity-run height",
        ),
        (
            "all_link_fit_ok_flips",
            "all_orbit_links.fit_ok_flip_n",
            "linkages",
            "fit_ok classification changes in the 960 m sensitivity run",
        ),
        (
            "all_link_is_good_flips",
            "all_orbit_links.is_good_flip_n",
            "linkages",
            "is_good classification changes in the 960 m sensitivity run",
        ),
        (
            "formal_unknown_fit_ok_flips",
            "formal_unknown_catalog.fit_ok_flip_n",
            "linkages",
            "Formal-unknown fit_ok changes in the 960 m sensitivity run",
        ),
        (
            "formal_unknown_is_good_flips",
            "formal_unknown_catalog.is_good_flip_n",
            "linkages",
            "Formal-unknown is_good changes in the 960 m sensitivity run",
        ),
        (
            "retained58_fit_ok_flips",
            "high_confidence_58.fit_ok_flip_n",
            "linkages",
            "Retained-sample fit_ok changes in the 960 m sensitivity run",
        ),
        (
            "retained58_is_good_flips",
            "high_confidence_58.is_good_flip_n",
            "linkages",
            "Retained-sample is_good changes in the 960 m sensitivity run",
        ),
        (
            "site_sensitivity_guardrail",
            "interpretation_guardrail",
            "",
            "Interpretation boundary for the site sensitivity run",
        ),
    )
    for parameter, dotted, unit, definition in orbit_specs:
        add(
            "orbit_site_sensitivity",
            parameter,
            deep_get(orbit_site_summary, dotted),
            unit,
            definition,
            "orbit_site_summary",
        )
    if orbit_site_summary is not None:
        for check_id, config_key, summary_key in (
            (
                "table1.orbit_baseline_height_vs_config",
                "site_configuration_as_executed.single_night_orbit_confirm.height_m",
                "baseline.observer.height_m",
            ),
            (
                "table1.orbit_alternate_height_vs_scheduler_config",
                "site_configuration_as_executed.scheduler.height_m",
                "alternate_run_summary.observer.height_m",
            ),
        ):
            assert_close(
                checks,
                check_id,
                deep_get(config, config_key),
                deep_get(orbit_site_summary, summary_key),
                left_label=f"config.{config_key}",
                right_label=f"orbit_site_summary.{summary_key}",
            )
    return pd.DataFrame(rows, columns=TABLE1_COLUMNS)


def product_row_count(row_counts: dict[str, Any] | None, product: str) -> Any:
    return deep_get(row_counts, f"products.{product}.rows_written")


def table2_data_accounting(
    snapshot: dict[str, Any] | None,
    coverage: dict[str, Any] | None,
    row_counts: dict[str, Any] | None,
    night_status: pd.DataFrame | None,
    *,
    input_hash: str,
    checks: list[dict[str, Any]],
) -> pd.DataFrame:
    night_frame: pd.DataFrame | None = None
    status_sums: dict[str, int] = {}
    if night_status is not None:
        require_columns(night_status, {"night"}, Path("night_status"))
        night_frame = night_status.copy()
        night_frame["night"] = night_frame["night"].map(normalize_night)
        if night_frame["night"].eq("").any() or night_frame["night"].duplicated().any():
            raise ValueError("night_status night keys must be valid and unique")
        for column in (
            "raw_science_n",
            "raw_all_mp_n",
            "l2_mp_n",
            "known_predicted_n",
            "known_match1_n",
            "known_mask15_n",
            "known_ades_n",
            "unknown_n",
            "review_real_n",
            "submit_real_n",
            "audit_initial_n",
            "audit_real_n",
        ):
            if column not in night_frame:
                continue
            values = pd.to_numeric(night_frame[column], errors="coerce")
            invalid = (
                values.isna() | values.lt(0) | ~np.isclose(values, np.round(values))
            )
            if invalid.any():
                examples = (
                    night_frame.loc[invalid, ["night", column]]
                    .head()
                    .to_dict("records")
                )
                raise ValueError(
                    f"night_status {column} must contain non-negative integer counts: {examples}"
                )
            status_sums[column] = int(values.sum())

    candidates: dict[str, list[tuple[str, Any]]] = {
        "strict_raw_frames": [
            ("snapshot_summary", deep_get(snapshot, "strict_raw_frames")),
            ("coverage_summary", deep_get(coverage, "raw_exposure_n")),
            ("night_status", status_sums.get("raw_science_n")),
        ],
        "strict_raw_nights": [
            ("snapshot_summary", deep_get(snapshot, "strict_raw_nights")),
            ("coverage_summary", deep_get(coverage, "raw_night_n")),
            (
                "night_status",
                (
                    int((pd.to_numeric(night_frame["raw_science_n"]) > 0).sum())
                    if night_frame is not None and "raw_science_n" in night_frame
                    else None
                ),
            ),
        ],
        "observed_unique_fields": [
            ("snapshot_summary", deep_get(snapshot, "raw_fields")),
            ("coverage_summary", deep_get(coverage, "observed_field_n")),
        ],
        "open_shutter_hours": [
            ("snapshot_summary", deep_get(snapshot, "open_shutter_hours")),
            ("coverage_summary", deep_get(coverage, "open_shutter_hours")),
        ],
        "unique_healpix_area": [
            ("coverage_summary", deep_get(coverage, "unique_area_deg2"))
        ],
        "all_mp_fits": [
            ("snapshot_summary", deep_get(snapshot, "all_mp_fits")),
            ("night_status", status_sums.get("raw_all_mp_n")),
        ],
        "strict_l2_catalogs": [
            ("snapshot_summary", deep_get(snapshot, "l2_catalogs")),
            ("night_status", status_sums.get("l2_mp_n")),
        ],
        "strict_l2_nights": [
            (
                "night_status",
                (
                    int((pd.to_numeric(night_frame["l2_mp_n"]) > 0).sum())
                    if night_frame is not None and "l2_mp_n" in night_frame
                    else None
                ),
            )
        ],
        "known_prediction_rows": [
            ("snapshot_summary", deep_get(snapshot, "known_predictions")),
            ("frozen_row_counts", product_row_count(row_counts, "known_all")),
            ("night_status", status_sums.get("known_predicted_n")),
        ],
        "known_matches_1arcsec": [
            ("snapshot_summary", deep_get(snapshot, "known_matches_1arcsec")),
            ("frozen_row_counts", product_row_count(row_counts, "known_matched")),
            ("night_status", status_sums.get("known_match1_n")),
        ],
        "known_matches_1p5arcsec": [
            ("snapshot_summary", deep_get(snapshot, "known_matches_1p5arcsec")),
            ("frozen_row_counts", product_row_count(row_counts, "known_mask15")),
            ("night_status", status_sums.get("known_mask15_n")),
        ],
        "known_ades_observation_rows": [
            ("night_status", status_sums.get("known_ades_n"))
        ],
        "unknown_catalog_linkages": [
            ("snapshot_summary", deep_get(snapshot, "unknown_catalog_linkages")),
            ("frozen_row_counts", product_row_count(row_counts, "unknown_catalog")),
            ("night_status", status_sums.get("unknown_n")),
        ],
        "human_review_marked_real_linkages": [
            (
                "snapshot_summary",
                deep_get(snapshot, "human_review_marked_real_linkages"),
            ),
            ("night_status", status_sums.get("review_real_n")),
        ],
        "submission_selected_linkages": [
            ("snapshot_summary", deep_get(snapshot, "submission_selected_linkages")),
            ("night_status", status_sums.get("submit_real_n")),
        ],
        "initial_human_selected_linkages": [
            (
                "snapshot_summary",
                deep_get(snapshot, "initial_review_selected_linkages"),
            ),
            ("night_status", status_sums.get("audit_initial_n")),
        ],
        "posthoc_retained_linkages": [
            ("snapshot_summary", deep_get(snapshot, "posthoc_retained_linkages")),
            ("night_status", status_sums.get("audit_real_n")),
        ],
        "posthoc_retained_detection_memberships": [
            (
                "snapshot_summary",
                deep_get(snapshot, "posthoc_retained_detection_memberships"),
            )
        ],
    }
    specs = (
        (
            "strict_raw_frames",
            "frames",
            "exposure",
            "Strict standard raw science frames",
        ),
        ("strict_raw_nights", "nights", "night", "Nights containing strict raw frames"),
        (
            "observed_unique_fields",
            "fields",
            "field",
            "Distinct scheduled field identifiers acquired",
        ),
        (
            "open_shutter_hours",
            "h",
            "exposure time",
            "Sum of archived strict-raw exposure times",
        ),
        (
            "unique_healpix_area",
            "deg2",
            "sky pixel",
            "Union area of observed footprint polygons",
        ),
        (
            "all_mp_fits",
            "files",
            "file",
            "All MP-named raw FITS including nonstandard/engineering names",
        ),
        (
            "strict_l2_catalogs",
            "catalogs",
            "exposure catalog",
            "Strict standard MP L2 catalogs",
        ),
        (
            "strict_l2_nights",
            "nights",
            "night",
            "Nights containing strict standard MP L2 catalogs",
        ),
        (
            "known_prediction_rows",
            "predictions",
            "prediction-exposure row",
            "Nominal known-object prediction rows",
        ),
        (
            "known_matches_1arcsec",
            "matched detections",
            "prediction-exposure row",
            "Known matches inside the official 1 arcsec radius",
        ),
        (
            "known_matches_1p5arcsec",
            "matched detections",
            "prediction-exposure row",
            "Known detections in the 1.5 arcsec subtraction mask",
        ),
        (
            "known_ades_observation_rows",
            "observation rows",
            "ADES observation row",
            "Known-object ADES rows generated in the frozen status ledger",
        ),
        (
            "unknown_catalog_linkages",
            "linkages",
            "single-night linkage",
            "Rows in frozen post-known unknown catalogs",
        ),
        (
            "human_review_marked_real_linkages",
            "linkages",
            "single-night linkage",
            "Rows marked real during the operational human review",
        ),
        (
            "submission_selected_linkages",
            "linkages",
            "single-night linkage",
            "Rows selected for unknown-object submission products",
        ),
        (
            "initial_human_selected_linkages",
            "linkages",
            "single-night linkage",
            "Initial human-selected linkage rows",
        ),
        (
            "posthoc_retained_linkages",
            "linkages",
            "single-night linkage",
            "Post-hoc retained high-confidence linkage rows",
        ),
        (
            "posthoc_retained_detection_memberships",
            "memberships",
            "linkage-detection membership",
            "Detection memberships in retained linkages",
        ),
    )
    rows: list[dict[str, Any]] = []
    for metric, unit, grain, definition in specs:
        available = [
            (source, value)
            for source, value in candidates[metric]
            if numeric(value) is not None
        ]
        if len(available) > 1:
            for other_source, other_value in available[1:]:
                assert_close(
                    checks,
                    f"table2.{metric}.{available[0][0]}_vs_{other_source}",
                    available[0][1],
                    other_value,
                    left_label=available[0][0],
                    right_label=other_source,
                )
        if available:
            source, value = available[0]
            evidence_status = "available"
        else:
            source, value, evidence_status = "", None, "pending_missing_frozen_metric"
        rows.append(
            {
                "row_order": len(rows) + 1,
                "metric": metric,
                "value": display_value(value),
                "unit": unit,
                "grain": grain,
                "evidence_status": evidence_status,
                "definition": definition,
                "source_input": source,
                "input_manifest_hash": input_hash,
            }
        )

    if night_frame is not None:
        for metric, column, definition in (
            (
                "calendar_nights_in_table",
                "night",
                "Calendar-night rows in frozen status table",
            ),
            (
                "primary_science_included_nights",
                "primary_science_included",
                "Nights included in primary survey statistics",
            ),
            (
                "unknown_science_included_nights",
                "unknown_science_included",
                "Nights included in unknown-pipeline statistics",
            ),
        ):
            if column == "night":
                value = len(night_frame)
            elif column in night_frame:
                value = int(night_frame[column].map(truthy).sum())
            else:
                value = None
            rows.append(
                {
                    "row_order": len(rows) + 1,
                    "metric": metric,
                    "value": display_value(value),
                    "unit": "nights",
                    "grain": "night",
                    "evidence_status": (
                        "available"
                        if value is not None
                        else "pending_missing_frozen_metric"
                    ),
                    "definition": definition,
                    "source_input": "night_status" if value is not None else "",
                    "input_manifest_hash": input_hash,
                }
            )
    return pd.DataFrame(rows, columns=TABLE2_COLUMNS)


def table3_known_recovery_astrometry(
    known: dict[str, Any] | None,
    random_shift: dict[str, Any] | None,
    known_residuals: pd.DataFrame | None,
    *,
    input_hash: str,
    checks: list[dict[str, Any]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    residual_quantiles: dict[str, float] = {}
    if known_residuals is not None:
        require_columns(
            known_residuals,
            {"radial_residual_arcsec"},
            Path("known_residuals"),
        )
        radial = pd.to_numeric(
            known_residuals["radial_residual_arcsec"], errors="coerce"
        ).to_numpy(dtype=float)
        if not len(radial) or not np.isfinite(radial).all() or (radial < 0).any():
            raise ValueError(
                "known_residuals radial_residual_arcsec must be non-empty, finite, and non-negative"
            )
        p50, p84, p90, p95 = np.percentile(radial, [50, 84, 90, 95])
        residual_quantiles = {
            "median_radial_residual_arcsec": float(p50),
            "p84_radial_residual_arcsec": float(p84),
            "p90_radial_residual_arcsec": float(p90),
            "p95_radial_residual_arcsec": float(p95),
        }
        expected_rows = numeric(deep_get(known, "matched_1arcsec_n"))
        if expected_rows is not None:
            assert_close(
                checks,
                "table3.residual_row_count",
                len(radial),
                expected_rows,
                left_label="known residual rows",
                right_label="matched_1arcsec_n",
            )
        truncation = numeric(deep_get(known, "radial_distribution_truncated_at_arcsec"))
        if truncation is not None and float(np.max(radial)) > truncation + 1.0e-9:
            raise ValueError(
                "known residuals exceed the reported radial truncation: "
                f"max={float(np.max(radial))}, threshold={truncation}"
            )
    denominator_text = normalize_text(
        deep_get(known, "prediction_denominator_definition")
    )
    nominal_caveat = (
        denominator_text or "pending: nominal denominator definition not supplied"
    )
    specs = (
        (
            "recovery",
            "source_prediction_rows_including_duplicates",
            "source_prediction_rows_including_duplicates",
            "rows",
            "",
            "Known prediction rows before duplicate-key collapse",
            nominal_caveat,
        ),
        (
            "recovery",
            "predicted_nominal_n",
            "predicted_nominal_n",
            "predictions",
            "predicted_nominal_n",
            "Unique nominal prediction keys",
            nominal_caveat,
        ),
        (
            "recovery",
            "matched_1arcsec_n",
            "matched_1arcsec_n",
            "matched detections",
            "predicted_nominal_n",
            "Predictions with a production match inside 1 arcsec",
            "Nearest-neighbor thresholded association",
        ),
        (
            "recovery",
            "match_fraction_nominal",
            "match_fraction_nominal",
            "fraction",
            "predicted_nominal_n",
            "matched_1arcsec_n / predicted_nominal_n",
            "Nominal recovery proxy; not detectable-object completeness",
        ),
        (
            "recovery",
            "matched_unique_identity_keys",
            "matched_unique_identity_keys",
            "identity keys",
            "",
            "Distinct frozen known identity keys among matches",
            "Identity-key count is not a unique-object survey denominator",
        ),
        (
            "data_quality",
            "prediction_key_duplicates",
            "prediction_key_duplicates",
            "rows",
            "",
            "Duplicate prediction keys removed from denominator",
            "Duplicates are retained as a reported quality diagnostic",
        ),
        (
            "data_quality",
            "matched_prediction_key_duplicates",
            "matched_prediction_key_duplicates",
            "rows",
            "",
            "Duplicate matched prediction keys",
            "Expected to be zero",
        ),
        (
            "astrometry",
            "median_radial_residual_arcsec",
            "median_radial_residual_arcsec",
            "arcsec",
            "matched_1arcsec_n",
            "Median radial predicted-measured separation",
            "Distribution is truncated by the 1 arcsec production match threshold",
        ),
        (
            "astrometry",
            "p84_radial_residual_arcsec",
            "p84_radial_residual_arcsec",
            "arcsec",
            "matched_1arcsec_n",
            "84th percentile radial separation",
            "Distribution is threshold-truncated",
        ),
        (
            "astrometry",
            "p90_radial_residual_arcsec",
            "p90_radial_residual_arcsec",
            "arcsec",
            "matched_1arcsec_n",
            "90th percentile radial separation",
            "Pending unless frozen by the known-population analysis; distribution is threshold-truncated",
        ),
        (
            "astrometry",
            "p95_radial_residual_arcsec",
            "p95_radial_residual_arcsec",
            "arcsec",
            "matched_1arcsec_n",
            "95th percentile radial separation",
            "Distribution is threshold-truncated",
        ),
        (
            "astrometry",
            "median_dra_cosdec_arcsec",
            "median_dra_cosdec_arcsec",
            "arcsec",
            "matched_1arcsec_n",
            "Median signed RA*cos(Dec) residual",
            "Signed component",
        ),
        (
            "astrometry",
            "median_ddec_arcsec",
            "median_ddec_arcsec",
            "arcsec",
            "matched_1arcsec_n",
            "Median signed Dec residual",
            "Signed component",
        ),
        (
            "astrometry",
            "robust_sigma_dra_arcsec",
            "robust_sigma_dra_arcsec",
            "arcsec",
            "matched_1arcsec_n",
            "0.7413 times RA residual IQR",
            "Robust component scale",
        ),
        (
            "astrometry",
            "robust_sigma_ddec_arcsec",
            "robust_sigma_ddec_arcsec",
            "arcsec",
            "matched_1arcsec_n",
            "0.7413 times Dec residual IQR",
            "Robust component scale",
        ),
        (
            "astrometry",
            "radial_distribution_truncated_at_arcsec",
            "radial_distribution_truncated_at_arcsec",
            "arcsec",
            "",
            "Production radial truncation",
            "Selection threshold, not an uncertainty estimate",
        ),
    )
    for section, metric, key, unit, denominator, definition, caveat in specs:
        summary_value = deep_get(known, key)
        residual_value = residual_quantiles.get(key)
        if numeric(summary_value) is not None and residual_value is not None:
            assert_close(
                checks,
                f"table3.{key}.summary_vs_residuals",
                summary_value,
                residual_value,
                left_label="known summary",
                right_label="recomputed residual percentile",
                tolerance=1.0e-12,
            )
        if numeric(summary_value) is not None:
            value = summary_value
            source = "known_summary"
        elif residual_value is not None:
            value = residual_value
            source = "known_residuals"
        else:
            value = None
            source = ""
        rows.append(
            {
                "row_order": len(rows) + 1,
                "section": section,
                "metric": metric,
                "value": display_value(value),
                "unit": unit,
                "denominator": denominator,
                "evidence_status": (
                    "available"
                    if numeric(value) is not None
                    else "pending_missing_known_summary"
                ),
                "definition": definition,
                "caveat": caveat,
                "source_input": source,
                "input_manifest_hash": input_hash,
            }
        )
    random_specs = (
        ("random_shift_prediction_trials", "trials", "Shifted prediction trials"),
        (
            "random_shift_matches_lt_1arcsec",
            "matches",
            "Shifted trials matched inside 1 arcsec",
        ),
        (
            "random_shift_chance_fraction",
            "fraction",
            "Position-shift chance-match fraction",
        ),
    )
    for key, unit, definition in random_specs:
        value = deep_get(random_shift, key)
        rows.append(
            {
                "row_order": len(rows) + 1,
                "section": "chance_match_control",
                "metric": key,
                "value": display_value(value),
                "unit": unit,
                "denominator": (
                    "random_shift_prediction_trials"
                    if key != "random_shift_prediction_trials"
                    else ""
                ),
                "evidence_status": (
                    "available"
                    if numeric(value) is not None
                    else "pending_missing_random_shift_summary"
                ),
                "definition": definition,
                "caveat": normalize_text(deep_get(random_shift, "method")),
                "source_input": (
                    "random_shift_summary" if numeric(value) is not None else ""
                ),
                "input_manifest_hash": input_hash,
            }
        )

    if known is not None:
        predicted = numeric(known.get("predicted_nominal_n"))
        matched = numeric(known.get("matched_1arcsec_n"))
        fraction = numeric(known.get("match_fraction_nominal"))
        if predicted is not None and matched is not None and fraction is not None:
            assert_close(
                checks,
                "table3.nominal_match_fraction_closure",
                fraction,
                matched / predicted if predicted else float("nan"),
                left_label="reported fraction",
                right_label="matched/predicted",
            )
    return pd.DataFrame(rows, columns=TABLE3_COLUMNS)


def nonnegative_count(value: Any, label: str) -> int | None:
    parsed = numeric(value)
    if parsed is None:
        return None
    if parsed < 0 or not math.isclose(parsed, round(parsed), abs_tol=1.0e-9):
        raise ValueError(f"{label} must be a non-negative integer count, got {value!r}")
    return int(round(parsed))


def table4_unknown_funnel_retention(
    stage_totals: dict[str, Any] | None,
    unknown: dict[str, Any],
    blinded_known: dict[str, Any] | None,
    *,
    stage_scope: str,
    input_hash: str,
    checks: list[dict[str, Any]],
) -> pd.DataFrame:
    scoped = deep_get(stage_totals, stage_scope) or {}
    if not isinstance(scoped, dict):
        raise ValueError(f"unknown stage totals {stage_scope} must be an object")
    specs = (
        (
            1,
            "l2_source_detection",
            "l2_detection_n",
            "detection",
            None,
            "L2 catalog source rows",
        ),
        (
            2,
            "magnitude_flag_survivor_detection",
            "mag_flag_prefilter_detection_n",
            "detection",
            None,
            "Magnitude/flag prefilter survivors",
        ),
        (
            3,
            "gaia_survivor_detection",
            "gaia_survivor_n",
            "detection",
            None,
            "Detections surviving Gaia static-source removal",
        ),
        (
            4,
            "grouped_gaia_input_detection",
            "grouped_gaia_input_detection_n",
            "detection",
            None,
            "Detections entering grouped Gaia-mask processing",
        ),
        (
            5,
            "common_area_survivor_detection",
            "common_area_survivor_detection_n",
            "detection",
            None,
            "Common-area stage survivors as executed",
        ),
        (
            6,
            "edge_shell_survivor_detection",
            "edge_shell_survivor_detection_n",
            "detection",
            None,
            "Detections surviving the edge-shell filter",
        ),
        (
            7,
            "reference_anchored_static_survivor_detection",
            "static_survivor_detection_n",
            "detection",
            None,
            "Detections surviving the reference-anchored static filter",
        ),
        (
            8,
            "two_point_tracklet",
            "tracklet_n",
            "tracklet",
            None,
            "Two-point tracklets; a new counting grain",
        ),
        (
            9,
            "shared_endpoint_linkage",
            "link_n",
            "linkage",
            None,
            "Shared-endpoint linkages; a new counting grain",
        ),
        (
            10,
            "orbit_fit_ok_linkage",
            "orbit_fit_n",
            "linkage",
            None,
            "Linkages with a successful numerical orbit fit",
        ),
        (
            11,
            "orbit_is_good_linkage",
            "orbit_is_good_n",
            "linkage",
            None,
            "Diagnostic parallel orbit-quality threshold",
        ),
        (
            12,
            "fit_ok_all_non_known_linkage",
            "orbit_fit_non_known_n",
            "linkage",
            None,
            "fit_ok linkages whose member detections are all non-known",
        ),
        (
            13,
            "is_good_all_non_known_linkage",
            "orbit_is_good_non_known_n",
            "linkage",
            None,
            "Diagnostic is_good and all-non-known intersection",
        ),
        (
            14,
            "post_known_catalog_linkage",
            "unknown_n",
            "linkage",
            "catalog_linkages",
            "Rows in the frozen post-known unknown catalogs",
        ),
        (
            15,
            "initial_human_selected_linkage",
            "audit_initial_n",
            "linkage",
            "initial_human_selected_linkages",
            "Initial human-selected audit-ledger linkages",
        ),
        (
            16,
            "posthoc_retained_linkage",
            "audit_real_n",
            "linkage",
            "posthoc_retained_linkages",
            "Post-hoc retained high-confidence linkages",
        ),
    )
    denominator_map = {
        "magnitude_flag_survivor_detection": "l2_source_detection",
        "gaia_survivor_detection": "magnitude_flag_survivor_detection",
        "grouped_gaia_input_detection": "gaia_survivor_detection",
        "common_area_survivor_detection": "grouped_gaia_input_detection",
        "edge_shell_survivor_detection": "common_area_survivor_detection",
        "reference_anchored_static_survivor_detection": "edge_shell_survivor_detection",
        "orbit_fit_ok_linkage": "shared_endpoint_linkage",
        "orbit_is_good_linkage": "orbit_fit_ok_linkage",
        "fit_ok_all_non_known_linkage": "orbit_fit_ok_linkage",
        "is_good_all_non_known_linkage": "fit_ok_all_non_known_linkage",
        "post_known_catalog_linkage": "fit_ok_all_non_known_linkage",
        "initial_human_selected_linkage": "post_known_catalog_linkage",
        "posthoc_retained_linkage": "initial_human_selected_linkage",
    }
    values: dict[str, int | None] = {}
    sources: dict[str, str] = {}
    for _, stage, stage_key, _, summary_key, _ in specs:
        stage_value = scoped.get(stage_key) if stage_key else None
        summary_value = unknown.get(summary_key) if summary_key else None
        # The population summary covers the complete frozen catalog, whereas
        # ``included_unknown_nights`` deliberately removes quarantined nights.
        # Compare the sources only when their scopes are the same.
        if (
            stage_scope == "all_nights"
            and stage_value is not None
            and summary_value is not None
        ):
            assert_close(
                checks,
                f"table4.{stage}.stage_vs_population_summary",
                stage_value,
                summary_value,
                left_label=f"unknown_stage_totals.{stage_scope}.{stage_key}",
                right_label=f"unknown_population_summary.{summary_key}",
            )
        if stage_value is not None:
            values[stage] = nonnegative_count(stage_value, stage)
            sources[stage] = "unknown_stage_totals"
        elif summary_value is not None:
            values[stage] = nonnegative_count(summary_value, stage)
            sources[stage] = "unknown_population_summary"
        else:
            values[stage] = None
            sources[stage] = ""

    rows: list[dict[str, Any]] = []
    for order, stage, _, unit, _, definition in specs:
        value = values[stage]
        denominator_stage = denominator_map.get(stage, "")
        denominator = values.get(denominator_stage) if denominator_stage else None
        value_scope = (
            stage_scope
            if sources[stage] == "unknown_stage_totals"
            else "frozen_catalog_or_review_snapshot"
        )
        denominator_scope = (
            stage_scope
            if sources.get(denominator_stage) == "unknown_stage_totals"
            else "frozen_catalog_or_review_snapshot"
        )
        scopes_are_comparable = (
            not denominator_stage
            or value_scope == denominator_scope
            or (
                stage_scope == "all_nights"
                and {sources[stage], sources.get(denominator_stage)}
                == {"unknown_stage_totals", "unknown_population_summary"}
            )
        )
        fraction = (
            value / denominator
            if value is not None and denominator and scopes_are_comparable
            else None
        )
        if value is not None and denominator is not None and scopes_are_comparable:
            if value > denominator:
                raise ValueError(
                    f"unknown funnel is non-monotonic at {stage}: "
                    f"{value} > {denominator_stage}={denominator}"
                )
            record_check(
                checks,
                check_id=f"table4.{stage}.monotonic_count",
                status="pass",
                detail=f"{stage}={value}; {denominator_stage}={denominator}",
            )
        row_guardrail = UNKNOWN_RETENTION_GUARDRAIL
        if denominator_stage and not scopes_are_comparable:
            row_guardrail += (
                " Retention fraction intentionally blank because numerator and "
                "denominator scopes differ."
            )
        elif denominator_stage and denominator == 0:
            row_guardrail += (
                " Retention fraction undefined because the denominator is zero."
            )
        rows.append(
            {
                "stage_order": order,
                "stage": stage,
                "value": display_value(value),
                "unit": unit,
                "scope": value_scope,
                "evidence_status": (
                    "available"
                    if value is not None
                    else "pending_missing_frozen_stage_count"
                ),
                "denominator_stage": denominator_stage,
                "retention_fraction": display_value(fraction),
                "definition": definition,
                "interpretation_guardrail": row_guardrail,
                "source_input": sources[stage],
                "input_manifest_hash": input_hash,
            }
        )

    accounting_specs = (
        (
            20,
            "linkage_detection_membership",
            "linkage_detection_memberships",
            "memberships",
            "Detection memberships across catalog linkages",
        ),
        (
            21,
            "globally_unique_detection_key",
            "globally_unique_detection_keys",
            "detections",
            "Unique (night,image,objID) detection keys",
        ),
        (
            22,
            "duplicate_linkage_membership",
            "duplicate_memberships",
            "memberships",
            "Memberships beyond globally unique detection keys",
        ),
        (
            23,
            "posthoc_retained_detection_membership",
            "posthoc_retained_memberships",
            "memberships",
            "Detection memberships in retained linkages",
        ),
        (
            24,
            "posthoc_retained_unique_detection_key",
            "posthoc_retained_unique_detection_keys",
            "detections",
            "Unique detection keys in retained linkages",
        ),
    )
    for order, stage, key, unit, definition in accounting_specs:
        value = nonnegative_count(unknown.get(key), stage)
        rows.append(
            {
                "stage_order": order,
                "stage": stage,
                "value": display_value(value),
                "unit": unit,
                "scope": "frozen_unknown_catalog_memberships",
                "evidence_status": (
                    "available"
                    if value is not None
                    else "pending_missing_unknown_population_metric"
                ),
                "denominator_stage": "",
                "retention_fraction": "",
                "definition": definition,
                "interpretation_guardrail": (
                    "Membership, detection, linkage, and independent-object grains are distinct."
                ),
                "source_input": (
                    "unknown_population_summary" if value is not None else ""
                ),
                "input_manifest_hash": input_hash,
            }
        )
    memberships = numeric(unknown.get("linkage_detection_memberships"))
    unique = numeric(unknown.get("globally_unique_detection_keys"))
    duplicates = numeric(unknown.get("duplicate_memberships"))
    if memberships is not None and unique is not None and duplicates is not None:
        assert_close(
            checks,
            "table4.membership_duplicate_closure",
            memberships - unique,
            duplicates,
            left_label="memberships-unique detections",
            right_label="reported duplicate memberships",
        )

    blinded_scope = "retrospective_blinded_known_proxy_ge3_and_rate"
    blinded_guardrail = (
        normalize_text(deep_get(blinded_known, "guardrail"))
        or "Conditional survival proxy only; not L2 completeness or image injection."
    )
    eligible = nonnegative_count(
        deep_get(blinded_known, "survival_for_ge3_and_rate.eligible_object_n"),
        "blinded-known eligible object-night rows",
    )
    top_level_eligible = nonnegative_count(
        deep_get(blinded_known, "eligible_ge3_and_rate_n"),
        "blinded-known top-level eligible object-night rows",
    )
    if eligible is not None and top_level_eligible is not None:
        assert_close(
            checks,
            "table4.blinded_known_eligible_closure",
            eligible,
            top_level_eligible,
            left_label="nested eligible object-night rows",
            right_label="top-level eligible_ge3_and_rate_n",
        )
    blinded_specs = (
        (30, "blinded_known_eligible_ge3_rate_object_n", None, None),
        (
            31,
            "blinded_known_strict_linked_object_n",
            "strict_linked_n",
            "strict_linked_fraction",
        ),
        (
            32,
            "blinded_known_fit_ok_survived_object_n",
            "fit_ok_survived_n",
            "fit_ok_survived_fraction",
        ),
        (
            33,
            "blinded_known_is_good_survived_object_n",
            "is_good_survived_n",
            "is_good_survived_fraction",
        ),
    )
    for order, stage, count_key, fraction_key in blinded_specs:
        count = (
            eligible
            if count_key is None
            else nonnegative_count(
                deep_get(blinded_known, f"survival_for_ge3_and_rate.{count_key}"),
                stage,
            )
        )
        reported_fraction = (
            None
            if fraction_key is None
            else numeric(
                deep_get(
                    blinded_known,
                    f"survival_for_ge3_and_rate.{fraction_key}",
                )
            )
        )
        calculated_fraction = (
            count / eligible
            if count_key is not None and count is not None and eligible
            else None
        )
        if reported_fraction is not None and calculated_fraction is not None:
            assert_close(
                checks,
                f"table4.{stage}.fraction_closure",
                reported_fraction,
                calculated_fraction,
                left_label="reported blinded-known fraction",
                right_label="count/eligible",
            )
        if (
            count_key is not None
            and count is not None
            and eligible is not None
            and count > eligible
        ):
            raise ValueError(f"{stage}={count} exceeds eligible denominator={eligible}")
        rows.append(
            {
                "stage_order": order,
                "stage": stage,
                "value": display_value(count),
                "unit": "object-night rows",
                "scope": blinded_scope,
                "evidence_status": (
                    "available"
                    if count is not None
                    else "pending_missing_blinded_known_summary"
                ),
                "denominator_stage": (
                    "blinded_known_eligible_ge3_rate_object_n" if count_key else ""
                ),
                "retention_fraction": display_value(
                    reported_fraction
                    if reported_fraction is not None
                    else calculated_fraction
                ),
                "definition": (
                    "Eligible known object-night rows with >=3 detections and a 3--63 arcsec h^-1 rate"
                    if count_key is None
                    else "Identity-blind downstream linkage/orbit survival proxy"
                ),
                "interpretation_guardrail": blinded_guardrail,
                "source_input": "blinded_known_summary" if count is not None else "",
                "input_manifest_hash": input_hash,
            }
        )
    return pd.DataFrame(rows, columns=TABLE4_COLUMNS)


def normalize_link_keys(
    frame: pd.DataFrame,
    path: Path,
    *,
    night_column: str,
) -> pd.DataFrame:
    require_columns(frame, {night_column, "trk_sub", "linkage_id"}, path)
    result = frame.copy()
    result["_night"] = result[night_column].map(normalize_night)
    result["_trk_sub"] = result["trk_sub"].map(normalize_text)
    result["_linkage_id"] = result["linkage_id"].map(normalize_linkage_id)
    invalid = result[["_night", "_trk_sub", "_linkage_id"]].eq("").any(axis=1)
    if invalid.any():
        examples = (
            result.loc[invalid, [night_column, "trk_sub", "linkage_id"]]
            .head()
            .to_dict("records")
        )
        raise ValueError(f"{path}: invalid link keys: {examples}")
    result["_link_key"] = (
        result["_night"] + "|" + result["_trk_sub"] + "|" + result["_linkage_id"]
    )
    if result["_link_key"].duplicated().any():
        duplicates = (
            result.loc[result["_link_key"].duplicated(False), "_link_key"]
            .head()
            .tolist()
        )
        raise ValueError(f"{path}: duplicate link keys: {duplicates}")
    return result


def frame_lookup(frame: pd.DataFrame | None) -> dict[str, pd.Series]:
    if frame is None:
        return {}
    return {row["_link_key"]: row for _, row in frame.iterrows()}


def get_first(row: pd.Series | None, names: Sequence[str]) -> Any:
    if row is None:
        return ""
    for name in names:
        if name in row.index and normalize_text(row[name]) != "":
            return row[name]
    return ""


def verified_mpc_rows(
    frame: pd.DataFrame | None,
    *,
    checks: list[dict[str, Any]],
) -> pd.DataFrame | None:
    if frame is None:
        record_check(
            checks,
            check_id="table5.authoritative_mpc_evidence",
            status="not_run_input_not_supplied",
            severity="high",
            detail="All retained-link MPC states remain pending.",
        )
        return None
    require_columns(frame, {"evidence_source"}, Path("mpc_evidence"))
    if "evidence_verified" in frame:
        verified = frame["evidence_verified"].map(truthy)
    elif "evidence_status" in frame:
        verified = (
            frame["evidence_status"]
            .map(normalize_text)
            .str.lower()
            .isin({"verified", "authoritative", "confirmed", "evidence_verified"})
        )
    else:
        verified = pd.Series(False, index=frame.index)
    normalized_source = frame["evidence_source"].map(normalize_text)
    jpl_only_source = normalized_source.str.contains(
        r"jpl|horizons", case=False, regex=True
    )
    verified &= normalized_source.ne("") & ~jpl_only_source
    if "evidence_authority" in frame:
        authority = frame["evidence_authority"].map(normalize_text).str.lower()
        verified &= authority.isin(
            {"mpc", "minor planet center", "minor_planet_center"}
        )
    selected = frame.loc[verified].copy()
    record_check(
        checks,
        check_id="table5.authoritative_mpc_evidence",
        status="pass" if len(selected) else "no_verified_rows",
        severity="high",
        detail=(
            f"verified_rows={len(selected)}; supplied_rows={len(frame)}; "
            f"jpl_only_rows_rejected={int(jpl_only_source.sum())}"
        ),
    )
    return selected


def table5_retained_links(
    review: pd.DataFrame,
    retained_links: pd.DataFrame,
    cross_groups: pd.DataFrame | None,
    jpl_first: pd.DataFrame | None,
    jpl_second: pd.DataFrame | None,
    mpc_evidence: pd.DataFrame | None,
    *,
    expected_retained: int,
    input_hash: str,
    checks: list[dict[str, Any]],
) -> pd.DataFrame:
    require_columns(
        review,
        {"origin_night", "trk_sub", "linkage_id", "final_paper_status", "n_obs"},
        Path("review_status"),
    )
    review = normalize_link_keys(
        review, Path("review_status"), night_column="origin_night"
    )
    review = review.loc[
        review["final_paper_status"]
        .map(normalize_text)
        .eq("retained_after_posthoc_audit")
    ].copy()
    if len(review) != expected_retained:
        raise ValueError(
            f"review_status retained rows={len(review)}, expected {expected_retained}"
        )
    record_check(
        checks,
        check_id="table5.retained_review_row_count",
        status="pass",
        detail=f"retained_rows={len(review)}",
    )

    retained_links = normalize_link_keys(
        retained_links, Path("retained_links"), night_column="night"
    )
    retained_lookup = frame_lookup(retained_links)
    missing_retained = sorted(set(review["_link_key"]) - set(retained_lookup))
    extra_retained = sorted(set(retained_lookup) - set(review["_link_key"]))
    if missing_retained or extra_retained:
        raise ValueError(
            "retained-link/review key mismatch: "
            f"missing_metrics={missing_retained[:5]}, extra_metrics={extra_retained[:5]}"
        )
    record_check(
        checks,
        check_id="table5.retained_link_key_coverage",
        status="pass",
        detail=f"one-to-one keys={len(review)}",
    )

    if cross_groups is not None:
        cross_groups = normalize_link_keys(
            cross_groups, Path("cross_night_groups"), night_column="night"
        )
        missing = sorted(set(review["_link_key"]) - set(cross_groups["_link_key"]))
        if missing:
            raise ValueError(
                f"cross-night assignments miss retained keys: {missing[:5]}"
            )
        require_columns(
            cross_groups,
            {
                "linear_motion_candidate_group_id",
                "candidate_group_link_count",
                "candidate_group_distinct_night_count",
            },
            Path("cross_night_groups"),
        )
        retained_groups = cross_groups.loc[
            cross_groups["_link_key"].isin(set(review["_link_key"]))
        ].copy()
        if (
            retained_groups["linear_motion_candidate_group_id"]
            .map(normalize_text)
            .eq("")
            .any()
        ):
            raise ValueError(
                "cross-night assignments contain a blank candidate group ID"
            )
        for group_id, group in retained_groups.groupby(
            "linear_motion_candidate_group_id"
        ):
            claimed_links = {
                nonnegative_count(value, f"{group_id} candidate_group_link_count")
                for value in group["candidate_group_link_count"]
            }
            claimed_nights = {
                nonnegative_count(
                    value, f"{group_id} candidate_group_distinct_night_count"
                )
                for value in group["candidate_group_distinct_night_count"]
            }
            if claimed_links != {len(group)}:
                raise ValueError(
                    f"{group_id}: claimed link counts {claimed_links} != {len(group)}"
                )
            actual_nights = group["_night"].nunique()
            if claimed_nights != {actual_nights}:
                raise ValueError(
                    f"{group_id}: claimed night counts {claimed_nights} != {actual_nights}"
                )
        record_check(
            checks,
            check_id="table5.linear_motion_candidate_group_closure",
            status="pass",
            detail=(
                f"assigned_links={len(retained_groups)}; "
                f"candidate_components={retained_groups['linear_motion_candidate_group_id'].nunique()}"
            ),
        )
    if jpl_first is not None:
        first_night_column = "night" if "night" in jpl_first else "origin_night"
        jpl_first = normalize_link_keys(
            jpl_first, Path("jpl_first_pass"), night_column=first_night_column
        )
        extra = sorted(set(jpl_first["_link_key"]) - set(review["_link_key"]))
        if extra:
            raise ValueError(
                f"JPL first-pass rows fall outside retained sample: {extra[:5]}"
            )
    if jpl_second is not None:
        second_night_column = "night" if "night" in jpl_second else "origin_night"
        jpl_second = normalize_link_keys(
            jpl_second, Path("jpl_second_pass"), night_column=second_night_column
        )
        extra = sorted(set(jpl_second["_link_key"]) - set(review["_link_key"]))
        if extra:
            raise ValueError(
                f"JPL second-pass rows fall outside retained sample: {extra[:5]}"
            )
    if mpc_evidence is not None:
        mpc_night_column = "night" if "night" in mpc_evidence else "origin_night"
        mpc_evidence = normalize_link_keys(
            mpc_evidence, Path("mpc_evidence"), night_column=mpc_night_column
        )
    mpc_evidence = verified_mpc_rows(mpc_evidence, checks=checks)

    cross_lookup = frame_lookup(cross_groups)
    first_lookup = frame_lookup(jpl_first)
    second_lookup = frame_lookup(jpl_second)
    mpc_lookup = frame_lookup(mpc_evidence)
    rows: list[dict[str, Any]] = []
    review = review.sort_values(["_night", "_trk_sub", "_linkage_id"])
    for _, review_row in review.iterrows():
        key = review_row["_link_key"]
        link = retained_lookup[key]
        group = cross_lookup.get(key)
        first = first_lookup.get(key)
        second = second_lookup.get(key)
        mpc = mpc_lookup.get(key)

        review_n = numeric(review_row.get("n_obs"))
        link_n = numeric(link.get("n_obs"))
        if (
            review_n is not None
            and link_n is not None
            and not math.isclose(review_n, link_n)
        ):
            raise ValueError(
                f"{key}: review n_obs={review_n:g}, retained-link n_obs={link_n:g}"
            )

        first_name = normalize_text(
            get_first(first, ("jpl_object_name", "first_pass_object_name"))
        )
        first_strict = truthy(get_first(first, ("strict_candidate_association",)))
        first_plausible = truthy(get_first(first, ("plausible_candidate_association",)))
        if first is None:
            first_status = (
                "not_supplied" if jpl_first is None else "not_in_first_pass_table"
            )
        elif first_strict:
            first_status = "strict_candidate_association"
        elif first_plausible:
            first_status = "plausible_candidate_association"
        elif first_name:
            first_status = "candidate_within_search_radius_not_associated"
        else:
            first_status = (
                normalize_text(get_first(first, ("jpl_query_state",)))
                or "no_candidate_row"
            )

        second_confirmed = truthy(
            get_first(second, ("numerically_confirmed_candidate",))
        )
        if second is None:
            second_status = (
                "not_supplied"
                if jpl_second is None
                else "not_selected_or_no_second_pass_row"
            )
        elif second_confirmed:
            second_status = "numerically_confirmed_jpl_candidate_association"
        else:
            second_status = "second_pass_row_not_confirmed"
        second_name = normalize_text(
            get_first(second, ("second_pass_object_name", "jpl_object_name"))
        )

        if mpc is None:
            mpc_status = MPC_PENDING
            mpc_source = ""
            submission_id = ""
            ingest_state = "pending"
            mpc_identification = ""
            designation = ""
            known_object_id = ""
        else:
            mpc_status = "verified_authoritative_mpc_evidence"
            mpc_source = normalize_text(mpc.get("evidence_source"))
            submission_id = normalize_text(mpc.get("submission_id"))
            ingest_state = normalize_text(mpc.get("mpc_ingest_state"))
            mpc_identification = normalize_text(mpc.get("mpc_identification"))
            designation = normalize_text(mpc.get("provisional_designation"))
            known_object_id = normalize_text(mpc.get("known_object_id"))

        rows.append(
            {
                "row_order": len(rows) + 1,
                "link_key": key,
                "origin_night": review_row["_night"],
                "trk_sub": review_row["_trk_sub"],
                "linkage_id": review_row["_linkage_id"],
                "n_obs": display_value(link_n if link_n is not None else review_n),
                "n_tracklets": display_value(get_first(link, ("n_tracklets",))),
                "unique_detection_n": display_value(
                    get_first(link, ("unique_detection_n",))
                ),
                "fit_ok": display_value(get_first(link, ("fit_ok",))),
                "is_good": display_value(get_first(link, ("is_good",))),
                "rms_arcsec": display_value(get_first(link, ("rms_arcsec",))),
                "median_residual_arcsec": display_value(
                    get_first(link, ("med_arcsec", "median_residual_arcsec"))
                ),
                "max_arcsec": display_value(get_first(link, ("max_arcsec",))),
                "a_au": display_value(get_first(link, ("a_au",))),
                "ecc": display_value(get_first(link, ("ecc",))),
                "inc_deg": display_value(get_first(link, ("inc_deg",))),
                "raan_deg": display_value(get_first(link, ("raan_deg",))),
                "argp_deg": display_value(get_first(link, ("argp_deg",))),
                "nu_deg": display_value(get_first(link, ("nu_deg",))),
                "best_v1_kms": display_value(get_first(link, ("best_v1_kms",))),
                "linear_rms_arcsec": display_value(
                    get_first(link, ("lin_rms_arcsec", "linear_rms_arcsec"))
                ),
                "speed_arcsec_per_hour": display_value(
                    get_first(link, ("speed_arcsec_per_hour",))
                ),
                "direction_deg": display_value(
                    get_first(link, ("lin_dir_deg", "direction_deg"))
                ),
                "first_mjd": display_value(get_first(link, ("first_mjd",))),
                "median_mjd": display_value(get_first(link, ("median_mjd",))),
                "last_mjd": display_value(get_first(link, ("last_mjd",))),
                "median_ra_deg": display_value(get_first(link, ("median_ra_deg",))),
                "median_dec_deg": display_value(get_first(link, ("median_dec_deg",))),
                "median_mag_aper4": display_value(
                    get_first(link, ("median_mag_aper4",))
                ),
                "time_span_hours": display_value(get_first(link, ("time_span_hours",))),
                "ecliptic_lon_j2000_deg": display_value(
                    get_first(link, ("ecliptic_lon_j2000_deg",))
                ),
                "ecliptic_lat_j2000_deg": display_value(
                    get_first(link, ("ecliptic_lat_j2000_deg",))
                ),
                "solar_elongation_deg": display_value(
                    get_first(link, ("solar_elongation_deg",))
                ),
                "nearest_astronomical_twilight": display_value(
                    get_first(link, ("nearest_astronomical_twilight",))
                ),
                "nearest_twilight_abs_hours": display_value(
                    get_first(link, ("nearest_twilight_abs_hours",))
                ),
                "nearest_twilight_signed_hours": display_value(
                    get_first(link, ("nearest_twilight_signed_hours",))
                ),
                "final_paper_status": "retained_after_posthoc_audit",
                "review_notes": normalize_text(review_row.get("notes")),
                "source_gif": normalize_text(review_row.get("source_gif")),
                "linear_motion_candidate_group_id": display_value(
                    get_first(
                        group,
                        ("linear_motion_candidate_group_id", "cross_night_group_id"),
                    )
                ),
                "candidate_group_link_count": display_value(
                    get_first(group, ("candidate_group_link_count",))
                ),
                "candidate_group_distinct_night_count": display_value(
                    get_first(group, ("candidate_group_distinct_night_count",))
                ),
                "candidate_group_cross_night_primary_edge_count": display_value(
                    get_first(
                        group, ("candidate_group_cross_night_primary_edge_count",)
                    )
                ),
                "candidate_group_has_cross_night_primary_edge": display_value(
                    get_first(group, ("candidate_group_has_cross_night_primary_edge",))
                ),
                "candidate_group_is_singleton": display_value(
                    get_first(group, ("candidate_group_is_singleton",))
                ),
                "candidate_group_interpretation": normalize_text(
                    get_first(group, ("candidate_interpretation",))
                )
                or "pending_cross_night_candidate_screen",
                "jpl_first_pass_status": first_status,
                "jpl_first_pass_object_name": first_name,
                "jpl_first_pass_separation_arcsec": display_value(
                    get_first(first, ("separation_arcsec", "jpl_separation_arcsec"))
                ),
                "jpl_first_pass_strict_candidate": (
                    first_strict if first is not None else ""
                ),
                "jpl_first_pass_plausible_candidate": (
                    first_plausible if first is not None else ""
                ),
                "jpl_second_pass_status": second_status,
                "jpl_second_pass_object_name": second_name,
                "jpl_second_pass_separation_arcsec": display_value(
                    get_first(second, ("second_pass_separation_arcsec",))
                ),
                "jpl_second_pass_numerically_confirmed_candidate": (
                    second_confirmed if second is not None else ""
                ),
                "jpl_evidence_scope": JPL_GUARDRAIL,
                "mpc_evidence_status": mpc_status,
                "mpc_evidence_source": mpc_source,
                "mpc_submission_id": submission_id,
                "mpc_ingest_state": ingest_state,
                "mpc_identification": mpc_identification,
                "provisional_designation": designation,
                "known_object_id": known_object_id,
                "identity_guardrail": RETAINED_GUARDRAIL,
                "retained_link_source_path": normalize_text(
                    get_first(link, ("provenance_path", "source_file"))
                ),
                "input_manifest_hash": input_hash,
            }
        )
    result = pd.DataFrame(rows, columns=TABLE5_COLUMNS)
    jpl_confirmed = result["jpl_second_pass_numerically_confirmed_candidate"].map(
        truthy
    )
    pending_mpc = result["mpc_evidence_status"].eq(MPC_PENDING)
    pending_identity = result.loc[
        pending_mpc,
        [
            "mpc_submission_id",
            "mpc_identification",
            "provisional_designation",
            "known_object_id",
        ],
    ].apply(lambda column: column.map(normalize_text))
    if pending_identity.ne("").any().any():
        raise ValueError(
            "pending MPC rows contain populated authoritative identity fields"
        )
    record_check(
        checks,
        check_id="table5.jpl_does_not_promote_mpc_state",
        status="pass",
        detail=(
            f"JPL_second_pass_confirmed={int(jpl_confirmed.sum())}; "
            f"MPC_pending={int(pending_mpc.sum())}; no JPL-derived MPC fields"
        ),
    )
    return result


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (pd.Timestamp, datetime)):
        return value.isoformat()
    if value is pd.NA:
        return None
    return value


def latex_view(name: str, frame: pd.DataFrame) -> pd.DataFrame:
    columns = {
        "table1": ["section", "parameter", "value", "unit", "evidence_status"],
        "table2": ["metric", "value", "unit", "grain", "evidence_status"],
        "table3": ["section", "metric", "value", "unit", "caveat"],
        "table4": ["stage", "value", "unit", "retention_fraction", "evidence_status"],
        "table5": [
            "origin_night",
            "trk_sub",
            "linkage_id",
            "n_obs",
            "rms_arcsec",
            "linear_motion_candidate_group_id",
            "jpl_second_pass_object_name",
            "mpc_evidence_status",
            "provisional_designation",
        ],
    }[name]
    return frame[columns]


def write_outputs(
    output_dir: Path,
    tables: dict[str, pd.DataFrame],
    summary_base: dict[str, Any],
    *,
    emit_latex: bool,
) -> None:
    if output_dir.exists():
        raise FileExistsError(
            f"refusing to overwrite existing output directory: {output_dir}"
        )
    temporary = output_dir.with_name(output_dir.name + ".inprogress")
    if temporary.exists():
        raise FileExistsError(f"staged output directory already exists: {temporary}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary.mkdir()
    output_hashes: dict[str, str] = {}
    for name, frame in tables.items():
        csv_path = temporary / TABLE_FILENAMES[name]
        frame.to_csv(csv_path, index=False)
        output_hashes[csv_path.name] = sha256_file(csv_path)
        if emit_latex:
            tex_path = temporary / TABLE_FILENAMES[name].replace(".csv", ".tex")
            latex_view(name, frame).to_latex(
                tex_path,
                index=False,
                escape=True,
                longtable=name == "table5",
                na_rep="",
            )
            output_hashes[tex_path.name] = sha256_file(tex_path)
    summary = dict(summary_base)
    summary["outputs"] = {
        name: {
            "row_count": int(len(tables[key])),
            "sha256": output_hashes[name],
        }
        for key, name in TABLE_FILENAMES.items()
    }
    summary["output_hashes"] = output_hashes
    summary_path = temporary / SUMMARY_FILENAME
    summary_path.write_text(
        json.dumps(json_safe(summary), indent=2, ensure_ascii=False, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, output_dir)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--environment", type=Path)
    parser.add_argument("--production-hashes", type=Path)
    parser.add_argument("--orbit-site-summary", type=Path)
    parser.add_argument("--snapshot-summary", type=Path)
    parser.add_argument("--coverage-summary", type=Path)
    parser.add_argument("--night-status", type=Path)
    parser.add_argument("--frozen-row-counts", type=Path)
    parser.add_argument("--known-summary", type=Path)
    parser.add_argument(
        "--known-residuals",
        type=Path,
        help="Optional frozen Parquet with radial_residual_arcsec for percentile closure",
    )
    parser.add_argument("--random-shift-summary", type=Path)
    parser.add_argument("--unknown-stage-totals", type=Path)
    parser.add_argument("--unknown-population-summary", type=Path, required=True)
    parser.add_argument("--blinded-known-summary", type=Path)
    parser.add_argument("--review-status", type=Path, required=True)
    parser.add_argument("--retained-links", type=Path, required=True)
    parser.add_argument("--cross-night-groups", type=Path)
    parser.add_argument("--jpl-first-pass", type=Path)
    parser.add_argument("--jpl-second-pass", type=Path)
    parser.add_argument("--mpc-evidence", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--unknown-stage-scope",
        choices=("included_unknown_nights", "all_nights"),
        default="included_unknown_nights",
    )
    parser.add_argument(
        "--expected-retained-linkages", type=int, default=DEFAULT_EXPECTED_RETAINED
    )
    parser.add_argument("--emit-latex", action="store_true")
    return parser


def resolve_inputs(args: argparse.Namespace) -> tuple[dict[str, Path], Path]:
    names = (
        "config",
        "environment",
        "production_hashes",
        "orbit_site_summary",
        "snapshot_summary",
        "coverage_summary",
        "night_status",
        "frozen_row_counts",
        "known_summary",
        "known_residuals",
        "random_shift_summary",
        "unknown_stage_totals",
        "unknown_population_summary",
        "blinded_known_summary",
        "review_status",
        "retained_links",
        "cross_night_groups",
        "jpl_first_pass",
        "jpl_second_pass",
        "mpc_evidence",
    )
    inputs: dict[str, Path] = {}
    for name in names:
        value = getattr(args, name)
        if value is None:
            continue
        path = value.expanduser().resolve(strict=False)
        if not path.is_file():
            raise FileNotFoundError(f"{name}: {path}")
        inputs[name] = path
    output = args.output_dir.expanduser().resolve(strict=False)
    if output.exists() or output.with_name(output.name + ".inprogress").exists():
        raise FileExistsError(
            f"refusing to overwrite or reuse output directory: {output}"
        )
    if args.expected_retained_linkages <= 0:
        raise ValueError("--expected-retained-linkages must be positive")
    return inputs, output


def run(args: argparse.Namespace) -> None:
    inputs, output_dir = resolve_inputs(args)
    hashes = {name: sha256_file(path) for name, path in inputs.items()}
    input_hash = combined_hash(hashes)
    checks: list[dict[str, Any]] = []

    config = read_json(inputs["config"])
    environment = read_json(inputs["environment"]) if "environment" in inputs else None
    production_hashes = (
        read_csv(inputs["production_hashes"]) if "production_hashes" in inputs else None
    )
    orbit_site_summary = (
        read_json(inputs["orbit_site_summary"])
        if "orbit_site_summary" in inputs
        else None
    )
    snapshot = (
        read_json(inputs["snapshot_summary"]) if "snapshot_summary" in inputs else None
    )
    coverage = (
        read_json(inputs["coverage_summary"]) if "coverage_summary" in inputs else None
    )
    night_status = (
        read_csv(inputs["night_status"]) if "night_status" in inputs else None
    )
    row_counts = (
        read_json(inputs["frozen_row_counts"])
        if "frozen_row_counts" in inputs
        else None
    )
    known = read_json(inputs["known_summary"]) if "known_summary" in inputs else None
    known_residuals = (
        pd.read_parquet(inputs["known_residuals"], columns=["radial_residual_arcsec"])
        if "known_residuals" in inputs
        else None
    )
    random_shift = (
        read_json(inputs["random_shift_summary"])
        if "random_shift_summary" in inputs
        else None
    )
    stage_totals = (
        read_json(inputs["unknown_stage_totals"])
        if "unknown_stage_totals" in inputs
        else None
    )
    unknown = read_json(inputs["unknown_population_summary"])
    blinded_known = (
        read_json(inputs["blinded_known_summary"])
        if "blinded_known_summary" in inputs
        else None
    )
    review = read_csv(inputs["review_status"])
    retained_links = read_csv(inputs["retained_links"])
    cross_groups = (
        read_csv(inputs["cross_night_groups"])
        if "cross_night_groups" in inputs
        else None
    )
    jpl_first = (
        read_csv(inputs["jpl_first_pass"]) if "jpl_first_pass" in inputs else None
    )
    jpl_second = (
        read_csv(inputs["jpl_second_pass"]) if "jpl_second_pass" in inputs else None
    )
    mpc_evidence = (
        read_csv(inputs["mpc_evidence"]) if "mpc_evidence" in inputs else None
    )

    tables = {
        "table1": table1_configuration_environment(
            config,
            environment,
            production_hashes,
            orbit_site_summary,
            input_hash=input_hash,
            checks=checks,
        ),
        "table2": table2_data_accounting(
            snapshot,
            coverage,
            row_counts,
            night_status,
            input_hash=input_hash,
            checks=checks,
        ),
        "table3": table3_known_recovery_astrometry(
            known,
            random_shift,
            known_residuals,
            input_hash=input_hash,
            checks=checks,
        ),
        "table4": table4_unknown_funnel_retention(
            stage_totals,
            unknown,
            blinded_known,
            stage_scope=args.unknown_stage_scope,
            input_hash=input_hash,
            checks=checks,
        ),
        "table5": table5_retained_links(
            review,
            retained_links,
            cross_groups,
            jpl_first,
            jpl_second,
            mpc_evidence,
            expected_retained=args.expected_retained_linkages,
            input_hash=input_hash,
            checks=checks,
        ),
    }

    table4_retained = tables["table4"].loc[
        tables["table4"]["stage"].eq("posthoc_retained_linkage"), "value"
    ]
    if len(table4_retained) == 1 and numeric(table4_retained.iloc[0]) is not None:
        assert_close(
            checks,
            "table4_vs_table5.retained_linkage_count",
            table4_retained.iloc[0],
            len(tables["table5"]),
            left_label="Table 4 retained linkages",
            right_label="Table 5 rows",
        )

    pending_counts = {
        name: int(
            frame.apply(
                lambda column: column.astype(str).str.contains(
                    r"pending|not_available", case=False, regex=True
                )
            )
            .any(axis=1)
            .sum()
        )
        for name, frame in tables.items()
    }
    table5 = tables["table5"]
    table2_by_metric = tables["table2"].set_index("metric")["value"].to_dict()

    def table2_count(metric: str) -> int | None:
        return nonnegative_count(table2_by_metric.get(metric), f"Table 2 {metric}")

    selection_funnel = {
        "operational_review_marked_real_linkages": table2_count(
            "human_review_marked_real_linkages"
        ),
        "submission_selected_linkages": table2_count("submission_selected_linkages"),
        "audit_ledger_initial_linkages": table2_count(
            "initial_human_selected_linkages"
        ),
        "posthoc_retained_linkages": table2_count("posthoc_retained_linkages"),
    }
    review_chain = [
        selection_funnel["operational_review_marked_real_linkages"],
        selection_funnel["submission_selected_linkages"],
        selection_funnel["posthoc_retained_linkages"],
    ]
    if all(value is not None for value in review_chain):
        if not all(
            left >= right for left, right in zip(review_chain, review_chain[1:])
        ):
            raise ValueError(f"review funnel is non-monotonic: {review_chain}")
        record_check(
            checks,
            check_id="table5.review_selection_funnel",
            status="pass",
            detail=" -> ".join(str(value) for value in review_chain),
        )

    group_ids = table5["linear_motion_candidate_group_id"].map(normalize_text)
    assigned_groups = table5.loc[group_ids.ne("")].copy()
    component_sizes = (
        assigned_groups.groupby("linear_motion_candidate_group_id").size()
        if len(assigned_groups)
        else pd.Series(dtype="int64")
    )
    component_size_distribution = {
        str(int(size)): int(count)
        for size, count in component_sizes.value_counts().sort_index().items()
    }
    confirmed_jpl = table5.loc[
        table5["jpl_second_pass_numerically_confirmed_candidate"].map(truthy)
    ]
    jpl_object_counts = {
        name: int(count)
        for name, count in confirmed_jpl["jpl_second_pass_object_name"]
        .map(normalize_text)
        .value_counts()
        .items()
        if name
    }
    summary = {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_sources": {
            name: {"path": str(path), "sha256": hashes[name]}
            for name, path in inputs.items()
        },
        "combined_input_manifest_hash": input_hash,
        "data_quality_checks": checks,
        "pending_row_counts": pending_counts,
        "table5_evidence_counts": {
            "retained_link_rows": len(table5),
            "selection_funnel": selection_funnel,
            "linear_motion_candidate_group_assigned_rows": int(
                table5["linear_motion_candidate_group_id"]
                .map(normalize_text)
                .ne("")
                .sum()
            ),
            "linear_motion_candidate_component_count": int(len(component_sizes)),
            "linear_motion_candidate_component_size_distribution": component_size_distribution,
            "jpl_second_pass_numerically_confirmed_candidate_rows": int(
                table5["jpl_second_pass_numerically_confirmed_candidate"]
                .map(truthy)
                .sum()
            ),
            "jpl_second_pass_confirmed_object_counts": jpl_object_counts,
            "jpl_second_pass_confirmed_link_keys": confirmed_jpl["link_key"].tolist(),
            "authoritative_mpc_evidence_rows": int(
                table5["mpc_evidence_status"]
                .eq("verified_authoritative_mpc_evidence")
                .sum()
            ),
            "pending_mpc_evidence_rows": int(
                table5["mpc_evidence_status"].eq(MPC_PENDING).sum()
            ),
            "nonblank_provisional_designation_rows": int(
                table5["provisional_designation"].map(normalize_text).ne("").sum()
            ),
        },
        "guardrails": {
            "retained_link_grain": RETAINED_GUARDRAIL,
            "unknown_retention": UNKNOWN_RETENTION_GUARDRAIL,
            "jpl_vs_mpc": JPL_GUARDRAIL,
            "mpc_pending_rule": (
                "Without a separately supplied verified MPC evidence row, submission/"
                "identification/designation fields remain blank and ingest state is pending."
            ),
        },
        "latex_emitted": bool(args.emit_latex),
    }
    write_outputs(output_dir, tables, summary, emit_latex=args.emit_latex)
    print(
        f"[done] tables=5 retained_rows={len(table5)} "
        f"mpc_pending={summary['table5_evidence_counts']['pending_mpc_evidence_rows']} "
        f"output={output_dir}",
        flush=True,
    )


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        run(args)
    except (
        FileExistsError,
        FileNotFoundError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
