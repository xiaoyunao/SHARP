#!/usr/bin/env python3
"""Finalize the publication tables and input fingerprints for this handoff."""
import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version

import pandas as pd

from exposures import INPUT, ROOT, TABLES


def digest(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def main():
    exposures = read_json(TABLES / "observing_summary.json")
    performance = read_json(TABLES / "performance_summary.json")
    coverage = pd.read_csv(TABLES / "geometric_coverage.csv")
    area = coverage.loc[coverage.nside.eq(1024)].iloc[0]
    astrometry = pd.read_csv(TABLES / "astrometry_summary.csv")
    radial = astrometry.loc[astrometry["sample"].eq("qualified_1arcsec") &
                            astrometry.metric.eq("radial_arcsec"), "median"].iloc[0]
    results = [
        ("strict_raw_inventory", exposures["baseline_raw_file_inventory"], "files", "August frozen file count"),
        ("qualified_exposures", exposures["quality_exposures"], "exposures", "unique exposure keys, valid clocks and L2 WCS"),
        ("qualified_nights", exposures["quality_noon_nights"], "observing nights", "local noon boundary"),
        ("qualified_open_shutter", exposures["quality_shutter_hours"], "hours", "sum of exposure durations"),
        ("qualified_once", area.geometric_once_deg2, "deg2", "geometric exposure WCS sky union"),
        ("qualified_same_night_ge3", area.geometric_ge3_same_night_deg2, "deg2", "sky union of same-night>=3 exposure coverage"),
        ("recorded_field_label_plus_wcs_pointing_groups", exposures["field_nights"], "field-night pointings", "0.30 degree complete linkage"),
        ("field_night_pointings_ge3", exposures["ge3_span_n"], "field-night pointings", "at least three qualified midpoints"),
        ("adjacent_revisit_median", exposures["adjacent_min_p16_p50_p84"][1], "min", "successive qualified exposures for same pointing"),
        ("first_to_last_median", exposures["ge3_span_min_p16_p50_p84"][1], "min", "qualified groups with at least three visits"),
        ("median_depth_model_5sigma", performance["l2_quality_header_metrics"]["limiting_mag"]["median"], "mag", "catalog header MLIM"),
        ("median_1arcsec_selected_radial_oc", radial, "arcsec", "quality sample O-C conditional on 1 arcsec cut"),
        ("qualified_known_associations", performance["known"]["qualified_detected"], "associations", "1 arcsec matched detections"),
        ("qualified_unique_identities", performance["known"]["qualified_identified"], "identities", "unique known catalog identities"),
        ("conditional_link_retention", performance["retention"]["fraction"], "fraction", "frozen 1.5 arcsec detected eligible object-nights"),
        ("automatic_single_night_candidates", performance["candidates"]["auto_links"], "linkages", "frozen candidate catalog"),
        ("post_review_retained_linkages", performance["candidates"]["retained"], "linkages", "human review and frozen posthoc audit"),
    ]
    pd.DataFrame(results, columns=["quantity", "value", "unit", "definition"]).to_csv(
        TABLES / "main_results.csv", index=False)
    original_sources = read_json(INPUT / "original_hashes.json")
    original_sources["nightly_products.csv.gz"]["selected_columns"] = [
        "night", "known_match1_exists", "known_mask15_exists", "unknown_links_exists", "unknown_links_rows",
    ]
    local_input_names = {
        "de432s.bsp", "aperture_configuration_evidence.txt", "coordinate_and_calibration_source_evidence.txt",
        "instrument_and_examples.txt", "photometry_source_evidence.txt", "server_code_evidence.txt",
    }
    dependencies = {
        "distant_frozen_predictions": {
            "parent_product": "20260803 frozen_products/known_all.parquet",
            "parent_sha256": "378b29126ae85e977c4494d571f4a4903a427c5cc84beacf8411576ea776644a",
            "selected_records": 663,
            "selected_catalog_identities": 17,
            "selection_script": "scripts/select_rare_predictions.py",
        },
        "DE432s": {
            "download_url": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp",
            "sha256": "363f32e14f5255359ac32c4d38080cf28ab55564a5e16696a75f63394b666e9b",
            "local_path": "inputs/de432s.bsp",
        },
        "frozen_ASTORB": {
            "historic_path": "/pipeline/xiaoyunao/known_asteroid/astorb.dat",
            "historic_sha256": "ea22c7062e7349d9add34fb64835916d66705ad60a717cc7c86a0c47f5aa3040",
            "used_via": "byte-preserved frozen August orbit join; no dynamic identity query",
        },
        "IERS_A": {
            "historic_path": "/pipeline/xiaoyunao/survey/runtime/iers_backups/20261027T081258847508Z/finals2000A_previous.all",
            "archived_file_mtime_local": "2026-10-01",
            "included_as": "inputs/paper_iers.all",
        },
        "recorded_server_environment": {"astro-aleph": "1.0.0", "astropy": "7.1.0",
                                        "algorithm_provenance": "inputs/production_file_hashes.csv"},
    }
    manifest = {
        "analysis_generated_utc": datetime.now(timezone.utc).isoformat(),
        "history_repository": "https://github.com/xiaoyunao/SHARP",
        "history_archive_ref": "2dffff5e65081eeabd0dab0df2d1ac6ed2d2f757",
        "history_commit_utc": "2026-09-04T02:24:23Z",
        "snapshot_label": "20260803",
        "historical_algorithm_commit": "d2f0057bd4cdd4c7a3c6c2431f6449a72c327284",
        "paper_sample_archive_start": "2025-11-15",
        "paper_sample_archive_end_inclusive": "2026-07-15",
        "extension_inventory_cutoff_utc_exclusive": "2026-10-15T16:00:00Z",
        "extension_inventory_cutoff_local_exclusive": "2026-10-16T00:00:00+08:00",
        "science_time": "UTC exposure midpoint; local noon boundary (Asia/Shanghai)",
        "checked_products": "Read-only headers and preexisting evidence; no processing of post-cutoff science",
        "geometric_footprints": "40399 frozen L2 archive entries matched on file name and original size and modification time; qualified science subset used",
        "kept_summary_definitions": "Original ledger, selected WCS/sample quality, historical candidate status, historical conditional recovery",
        "dependencies": dependencies,
        "software": {"python": platform.python_version(),
                     **{name.split("==")[0]: version(name.split("==")[0])
                        for name in (ROOT / "requirements.txt").read_text().splitlines()}},
        "original_source_fingerprints": original_sources,
        "files_sha256": {str(path.relative_to(ROOT)): digest(path) for path in sorted(ROOT.rglob("*"))
                         if path.is_file() and "__pycache__" not in path.parts
                         and path.name not in {"freeze_manifest.json", ".gitignore"}
                         and not (path.parent == INPUT and path.name in local_input_names)
                         and not (path.parent == ROOT / "manuscript" and path.suffix in
                                  {".log", ".blg", ".aux", ".out", ".toc", ".fls", ".fdb_latexmk"})
                         and not any(p.startswith(".") for p in path.relative_to(ROOT).parts)},
    }
    (ROOT / "freeze_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
