# Snapshot validation report

- Overall status: **INCOMPLETE**
- Generated (UTC): `2026-08-09T13:01:51.516562+00:00`
- Snapshot label: `2026-08-03`
- Checks: 63 PASS, 0 FAIL, 1 SKIP

`SKIP` means the source artifact, author approval, or inspection backend is not yet available. It never means an observed value of zero. Expected headline values are comparison targets; observed values are independently computed.

## Headline closure

| Check | Status | Expected | Observed |
|---|---:|---:|---:|
| `inventory.strict_raw_exposures` | **PASS** | 41074 | 41074 |
| `inventory.strict_raw_nights` | **PASS** | 134 | 134 |
| `inventory.strict_raw_fields` | **PASS** | 1430 | 1430 |
| `inventory.engineering_raw_exposures` | **PASS** | 41152 | 41152 |
| `inventory.engineering_raw_nights` | **PASS** | 136 | 136 |
| `inventory.strict_l2_catalogs` | **PASS** | 40399 | 40399 |
| `inventory.strict_l2_nights` | **PASS** | 131 | 131 |
| `known.known_prediction_rows` | **PASS** | 13311871 | 13311871 |
| `known.known_matched_1arcsec_rows` | **PASS** | 534780 | 534780 |
| `known.known_mask15_rows` | **PASS** | 563612 | 563612 |
| `unknown.catalog_linkages` | **PASS** | 4762 | 4762 |
| `unknown.nonempty_nights` | **PASS** | 116 | 116 |
| `unknown.true_zero_nights` | **PASS** | 8 | 8 |
| `unknown.membership_rows` | **PASS** | 14299 | 14299 |
| `unknown.globally_unique_detections` | **PASS** | 14159 | 14159 |
| `review.initial_linkages` | **PASS** | 67 | 67 |
| `review.retained_linkages` | **PASS** | 58 | 58 |
| `review.rejected_linkages` | **PASS** | 9 | 9 |
| `review.retained_detection_memberships` | **PASS** | 179 | 179 |

## Detailed checks

### Configuration

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `config.snapshot` | **PASS** | {"schema_version": "1.0", "snapshot_label": "2026-08-03"} | {"schema_version": "1.0", "snapshot_label": "2026-08-03"} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/config/snapshot.json |

### Derived Known

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `known.derived_summary_reconciliation` | **PASS** | {"matched_1arcsec_n": 534780, "source_prediction_rows_including_duplicates": 13311871} | {"mismatches": {}, "missing_summary_keys": [], "summary_values": {"matched_1arcsec_n": 534780, "source_prediction_rows_including_duplicates": 13311871}, "unavailable_independent_values": []} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_known/known_population_summary.json |

### Figures

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `figures.fig01` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 894.24, "width_pt": 780.9420560748}], "pages": 1, "size_bytes": 31987}, "png": {"height_px": 3726, "size_bytes": 723821, "width_px": 3253}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig01_system_architecture.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig01_system_architecture.png |
| `figures.fig02` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 310.3615384615, "width_pt": 1288.9257722637}], "pages": 1, "size_bytes": 190879}, "png": {"height_px": 1293, "size_bytes": 1111489, "width_px": 5345}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig02_footprint_coverage.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig02_footprint_coverage.png |
| `figures.fig03` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 1029.168, "width_pt": 1241.9984975268}], "pages": 1, "size_bytes": 56069}, "png": {"height_px": 4288, "size_bytes": 1212714, "width_px": 5175}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig03_data_accounting.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig03_data_accounting.png |
| `figures.fig04` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 1001.6425, "width_pt": 1322.478125}], "pages": 1, "size_bytes": 148878}, "png": {"height_px": 4173, "size_bytes": 1246968, "width_px": 5509}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig04_scheduler_example.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig04_scheduler_example.png |
| `figures.fig05` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 970.919625, "width_pt": 1190.1935616477}], "pages": 1, "size_bytes": 80965}, "png": {"height_px": 4045, "size_bytes": 1050353, "width_px": 4959}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig05_known_method_and_residuals.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig05_known_method_and_residuals.png |
| `figures.fig06` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 684.745625, "width_pt": 1216.368}], "pages": 1, "size_bytes": 668747}, "png": {"height_px": 2853, "size_bytes": 1800107, "width_px": 5068}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig06_unknown_pipeline_examples.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig06_unknown_pipeline_examples.png |
| `figures.fig07` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 957.544, "width_pt": 1252.274375}], "pages": 1, "size_bytes": 69959}, "png": {"height_px": 3989, "size_bytes": 1524269, "width_px": 5220}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig07_nightly_exposure_history.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig07_nightly_exposure_history.png |
| `figures.fig08` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 1287.817985706, "width_pt": 1423.986709375}], "pages": 1, "size_bytes": 464315}, "png": {"height_px": 5364, "size_bytes": 2747265, "width_px": 5932}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig08_known_results.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig08_known_results.png |
| `figures.fig09` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 1125.1948, "width_pt": 1331.8525}], "pages": 1, "size_bytes": 90982}, "png": {"height_px": 4688, "size_bytes": 1212973, "width_px": 5549}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig09_unknown_funnel.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig09_unknown_funnel.png |
| `figures.fig10` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 971.495875, "width_pt": 1206.716875}], "pages": 1, "size_bytes": 105584}, "png": {"height_px": 4047, "size_bytes": 1156556, "width_px": 5028}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig10_high_confidence_distributions.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig10_high_confidence_distributions.png |
| `figures.fig11` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 1005.35825, "width_pt": 1267.2765136364}], "pages": 1, "size_bytes": 139273}, "png": {"height_px": 4188, "size_bytes": 1297235, "width_px": 5280}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig11_twilight_context.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig11_twilight_context.png |
| `figures.fig12` | **PASS** | {"pdf": true, "png": true, "positive_dimensions": true} | {"errors": [], "metadata": {"pdf": {"page_dimensions": [{"height_pt": 965.5405, "width_pt": 1354.0325}], "pages": 1, "size_bytes": 61707}, "png": {"height_px": 4022, "size_bytes": 1075008, "width_px": 5643}}, "present": {"pdf": true, "png": true}} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig12_operations_timeline.pdf, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures/fig12_operations_timeline.png |
| `figures.complete_set` | **PASS** | 12 | 12 | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/figures |

### Frozen Products

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `known.known_prediction_rows` | **PASS** | 13311871 | 13311871 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/known_all.parquet |
| `known.known_matched_1arcsec_rows` | **PASS** | 534780 | 534780 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/known_matched.parquet |
| `known.known_mask15_rows` | **PASS** | 563612 | 563612 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/known_mask15.parquet |
| `unknown.catalog_linkages` | **PASS** | 4762 | 4762 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/unknown_catalog.parquet |

### Inventory

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `inventory.strict_raw_exposures` | **PASS** | 41074 | 41074 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/raw_manifest.csv |
| `inventory.strict_raw_nights` | **PASS** | 134 | 134 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/raw_manifest.csv |
| `inventory.strict_raw_fields` | **PASS** | 1430 | 1430 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/raw_manifest.csv |
| `inventory.engineering_raw_exposures` | **PASS** | 41152 | 41152 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/raw_engineering_manifest.csv |
| `inventory.engineering_raw_nights` | **PASS** | 136 | 136 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/raw_engineering_manifest.csv |
| `inventory.strict_l2_catalogs` | **PASS** | 40399 | 40399 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/l2_manifest.csv |
| `inventory.strict_l2_nights` | **PASS** | 131 | 131 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/l2_manifest.csv |

### Provenance

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `provenance.raw_manifest_presence` | **PASS** | present and nonempty | {"exists": true, "size_bytes": 5955791} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/raw_manifest.csv |
| `provenance.engineering_manifest_presence` | **PASS** | present and nonempty | {"exists": true, "size_bytes": 5966277} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/raw_engineering_manifest.csv |
| `provenance.l2_manifest_presence` | **PASS** | present and nonempty | {"exists": true, "size_bytes": 25334548} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/l2_manifest.csv |
| `provenance.frozen_file_ledger_presence` | **PASS** | present and nonempty | {"exists": true, "size_bytes": 923746} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/file_status.csv |
| `provenance.review_ledger_presence` | **PASS** | present and nonempty | {"exists": true, "size_bytes": 18472} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/review_sample/review_and_mpc_status.csv |
| `provenance.frozen_product_hashes` | **PASS** | {"mismatched": 0, "missing": 0, "unrecorded": 0} | {"mismatched": [], "missing": [], "unrecorded": [], "verified": ["known_all.parquet", "known_mask15.parquet", "known_matched.parquet", "orbit_links.parquet", "orbit_obs_residuals.parquet", "unknown_catalog.parquet"]} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/row_counts.json, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/known_all.parquet, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/known_mask15.… |
| `provenance.review_hash_manifest` | **PASS** | {"malformed": 0, "mismatched": 0, "missing": 0} | {"malformed_lines": [], "mismatched": [], "missing": [], "verified": 9} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/review_sample/hashes.sha256 |
| `provenance.production_hash_manifest` | **PASS** | {"malformed": 0, "mismatched": 0, "missing": 0} | {"malformed_lines": [], "mismatched": [], "missing": [], "verified": 4} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/provenance/hashes.sha256 |
| `provenance.inventory_hash_manifest` | **PASS** | {"malformed": 0, "mismatched": 0, "missing": 0} | {"malformed_lines": [], "mismatched": [], "missing": [], "verified": 10} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/inventory/hashes.sha256 |
| `provenance.table_hash_manifest` | **PASS** | {"malformed": 0, "mismatched": 0, "missing": 0} | {"malformed_lines": [], "mismatched": [], "missing": [], "verified": 8} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/tables/hashes.sha256 |
| `provenance.derived_unknown_inputs` | **PASS** | {"mismatched": 0, "missing": 0} | {"mismatched": [], "missing": [], "verified": ["l2_manifest.csv", "unknown_catalog.parquet", "unknown_review_detections.parquet"]} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_population_summary.json |
| `provenance.derived_known_inputs` | **PASS** | {"mismatched": 0, "missing": 0} | {"mismatched": [], "missing": [], "verified": ["known_all.parquet", "known_matched.parquet", "l2_manifest.csv"]} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_known/known_population_summary.json |

### Quality Mask

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `quality.config_structure` | **PASS** | ["exclude_from_primary_science", "quarantine_unknown", "raw_without_l2"] | {"exclude_from_primary_science": ["20260111"], "quarantine_unknown": ["20251226", "20260528", "20260611"], "raw_without_l2": ["20260209", "20260312", "20260413"], "status": "provisional_pending_author_signoff"} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/config/snapshot.json |
| `quality.author_signoff` | **SKIP** | author-approved frozen mask | provisional_pending_author_signoff | A provisional mask is usable for engineering QA but is not final paper evidence. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/config/snapshot.json |
| `quality.table_reconciliation` | **PASS** | {"contradictions": 0, "night_rows": 243} | {"absent_nights": [], "duplicate_nights": 0, "extra_nights": [], "missing_excluded_codes": 0, "missing_quarantine_codes": 0, "missing_raw_without_l2_codes": 0, "night_rows": 243, "wrong_primary_flags": 0, "wrong_unknown_flags": 0} | Missing future rows are incomplete; contradictory existing flags fail. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/tables/quality_mask.csv |

### Review Sample

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `review.initial_linkages` | **PASS** | 67 | 67 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/review_sample/review_and_mpc_status.csv |
| `review.retained_linkages` | **PASS** | 58 | 58 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/review_sample/review_and_mpc_status.csv |
| `review.rejected_linkages` | **PASS** | 9 | 9 | Directly recomputed from the cited artifact. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/review_sample/review_and_mpc_status.csv |
| `review.audit_partition` | **PASS** | {"initial": 67, "rejected": 9, "retained": 58} | {"initial": 67, "ledger_rows": 67, "missing_key_rows": 0, "rejected": 9, "retained": 58, "retained_unique_keys": 58, "unrecognized_statuses": []} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/review_sample/review_and_mpc_status.csv |
| `review.retained_detection_memberships` | **PASS** | 179 | 179 | Rows are filtered by final_paper_status from the full membership table. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_all_link_memberships.csv |
| `review.retained_detection_identity` | **PASS** | {"membership_rows": 179, "missing_detection_key_rows": 0, "missing_link_key_rows": 0, "unique_detection_keys": 179, "unique_link_keys": 58} | {"membership_rows": 179, "missing_detection_key_rows": 0, "missing_link_key_rows": 0, "unique_detection_keys": 179, "unique_link_keys": 58} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_all_link_memberships.csv |
| `review.retained_key_consistency` | **PASS** | {"symmetric_difference": 0, "unique_link_keys_per_source": 58} | {"errors": [], "formal_catalog_missing_retained": 0, "missing_sources": [], "source_rows": {"derived_detections": 179, "derived_links": 58, "full_membership_retained": 179, "review_detections": 179, "review_links": 58, "review_status_retained": 58}, "source_unique_keys": {"derived_detections": 58, "derived_links": 58,… | A missing derivative is SKIP; a present source with a different key set is FAIL. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/review_sample/review_and_mpc_status.csv, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/review_sample/unknown_high_confidence_links.cs… |

### Tables

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `tables.snapshot_summary_reconciliation` | **PASS** | all available headline values equal independent counts | {"compared": {"all_mp_fits": {"independent": 41152, "summary": 41152}, "initial_review_selected_linkages": {"independent": 67, "summary": 67}, "known_matches_1arcsec": {"independent": 534780, "summary": 534780}, "known_matches_1p5arcsec": {"independent": 563612, "summary": 563612}, "known_predictions": {"independent":… | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/tables/snapshot_summary.json |

### Units

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `units.config_definitions` | **PASS** | {"known.official_match_radius_arcsec": 1.0, "known.unknown_mask_radius_arcsec": 1.5, "unknown.gaia_match_radius_arcsec": 1.5, "unknown.maximum_direction_difference_deg": 10.0, "unknown.maximum_magnitude_difference": 1.0, "unknown.maximum_speed_arcsec_per_hour": 63.0, "unknown.maximum_speed_difference_arcsec_per_hour":… | {"known.official_match_radius_arcsec": 1.0, "known.unknown_mask_radius_arcsec": 1.5, "unknown.gaia_match_radius_arcsec": 1.5, "unknown.maximum_direction_difference_deg": 10.0, "unknown.maximum_magnitude_difference": 1.0, "unknown.maximum_speed_arcsec_per_hour": 63.0, "unknown.maximum_speed_difference_arcsec_per_hour":… | The 15 s value is the standard 30 s exposure midpoint offset; it is a definition, not a correction applied here. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/config/snapshot.json |
| `units.unknown_time_reference` | **PASS** | {"orbit_and_tracklet_time_reference": "exposure_midpoint", "standard_exposure_midpoint_offset_s": 15.0, "unknown_ades_time_as_executed": "exposure_start"} | {"orbit_and_tracklet_time_reference": "exposure_midpoint", "standard_exposure_midpoint_offset_s": 15.0, "unknown_ades_time_as_executed": "exposure_start"} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/config/snapshot.json |
| `units.parquet_metadata` | **PASS** | unit metadata on all unit-bearing core columns | {"columns_checked": 98, "missing_metadata": [], "missing_products": [], "wrong": []} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/known_all.parquet, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/known_matched.parquet, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/known_m… |
| `units.unknown_membership_metadata` | **PASS** | explicit units for time, coordinates, photometry, and motion | {"missing_metadata": [], "wrong": []} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/unknown_review_detections.parquet |
| `units.unknown_speed_conversion` | **PASS** | {"bad_rows": 0, "formula": "per_day / 24"} | {"bad_rows": 0, "max_absolute_error_arcsec_per_hour": 1.4210854715202004e-14, "rows_checked": 4762} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_all_links.csv |
| `units.stage_definitions` | **PASS** | {"audit_real_n": "linkage", "gaia_survivor_n": "detection", "l2_detection_n": "detection", "link_n": "linkage", "mag_flag_prefilter_detection_n": "detection", "orbit_fit_n": "linkage", "orbit_fit_non_known_n": "linkage", "orbit_is_good_n": "linkage", "orbit_is_good_non_known_n": "linkage", "review_real_n": "linkage", … | {"definitions": {"audit_initial_n": "linkage", "audit_real_n": "linkage", "common_area_survivor_detection_n": "detection", "edge_shell_survivor_detection_n": "detection", "gaia_survivor_n": "detection", "grouped_gaia_input_detection_n": "detection", "l2_detection_n": "detection", "link_n": "linkage", "mag_flag_prefilt… | Missing definitions are unfinished; contradictory units fail. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/tables/unknown_stage_definitions.csv |

### Unknown Population

| Check | Status | Expected | Observed | Note |
|---|---:|---|---|---|
| `unknown.formal_quality_flags` | **PASS** | {"fit_ok_false_or_null": 0, "fit_ok_true": 4762, "is_good_false_or_null": 0, "is_good_true": 4762, "rows": 4762} | {"fit_ok_false_or_null": 0, "fit_ok_true": 4762, "is_good_false_or_null": 0, "is_good_true": 4762, "rows": 4762} | Flags are read from the formal Parquet rows, not from the derived summary JSON. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/unknown_catalog.parquet |
| `unknown.nonempty_nights` | **PASS** | 116 | 116 | Computed from exists/status/rows_written, so a missing file is not a true zero. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/file_status.csv |
| `unknown.true_zero_nights` | **PASS** | 8 | 8 | Computed from exists/status/rows_written, so a missing file is not a true zero. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/file_status.csv |
| `unknown.catalog_night_accounting` | **PASS** | {"nonempty_plus_zero": "116 + 8", "present_catalog_nights": 124} | {"nonempty_nights": 116, "present_catalog_nights": 124, "true_zero_nights": 8} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/frozen_products/file_status.csv |
| `unknown.membership_rows` | **PASS** | 14299 | 14299 | Every table row is one linkage-detection membership. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_all_link_memberships.csv |
| `unknown.globally_unique_detections` | **PASS** | 14159 | 14159 | Keys are independently reconstructed from night, image_name, and obj_id. Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_all_link_memberships.csv |
| `unknown.detection_key_completeness` | **PASS** | 0 | 0 | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_all_link_memberships.csv |
| `unknown.membership_unique_distinction` | **PASS** | {"difference": 140, "memberships": 14299, "unique": 14159} | {"difference": 140, "memberships": 14299, "unique": 14159} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_all_link_memberships.csv |
| `unknown.unique_table_reconciliation` | **PASS** | {"keys": 14159, "rows": 14159} | {"keys": 14159, "missing_key_rows": 0, "rows": 14159, "symmetric_difference": 0} | Sources: /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_all_link_memberships.csv, /Users/yunaoxiao/Desktop/smt_asteroid/paper_analysis_20260803/snapshot/derived_unknown/unknown_unique_detections.csv |

## Blocking failures

- None.

## Incomplete / awaiting input

- `quality.author_signoff` — The frozen quality mask has author sign-off.
