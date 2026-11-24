#!/usr/bin/env python3
"""Check frozen counts, cross-table identities, and deliverable files."""
import hashlib
import json
from datetime import datetime, timezone

import healpy as hp
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

from exposures import CUTOFF, DATA, INPUT, ROOT, TABLES


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    raw = json.loads((TABLES / "observing_summary.json").read_text())
    perf = json.loads((TABLES / "performance_summary.json").read_text())
    baseline = {
        "standard_raw_files": (41074, raw["baseline_raw_file_inventory"], "frozen raw inventory rows"),
        "raw_directory_nights": (134, raw["baseline_raw_directory_nights"], "distinct frozen raw directory labels"),
        "nominal_shutter_hours": (342.28333333333336, raw["baseline_raw_nominal_shutter_hours"], "nominal 30 s per raw inventory record"),
        "observed_fields": (1430, raw["baseline_raw_fields"], "unique field IDs in frozen raw inventory"),
        "nominal_geometry_deg2": (10387.42600221356, raw["baseline_nominal_coverage_deg2"], "frozen field polygons, NSIDE 512 inclusive"),
        "strict_l2": (40399, raw["frozen_l2_catalogs"], "frozen strict L2 rows"),
        "strict_l2_directory_nights": (131, raw["frozen_l2_directory_nights"], "frozen strict L2 directory labels"),
        "known_1arcsec": (534780, perf["known"]["baseline_detected"], "1 arcsec frozen matched rows"),
        "known_identities": (58482, perf["known"]["baseline_identified"], "distinct canonical identity keys"),
        "candidates": (4762, perf["candidates"]["auto_links"], "formal frozen candidate catalog rows"),
        "reviewed_real": (68, perf["candidates"]["marked_real"], "frozen review-true counts"),
        "selected_for_submission": (67, perf["candidates"]["selected"], "frozen selection rows"),
        "retained": (58, perf["candidates"]["retained"], "post-review retained linkages"),
        "retained_measurements": (179, perf["candidates"]["retained_measurements"], "corrected midpoint measurements"),
        "retained_nights": (34, perf["candidates"]["retained_nights"], "retained linkage nights"),
        "C2025_Y1_matches": (6, perf["candidates"]["c2025y1_linkages"], "frozen two-pass numerical ephemeris comparisons"),
        "retention_numerator": (75899, perf["retention"]["k"], "frozen eligible object-nights covered by pure link"),
        "retention_denominator": (102347, perf["retention"]["n"], "frozen >=3 detections; first-last speed 3–63 arcsec/h"),
    }
    checks = [{"metric": key, "expected": x, "recomputed": y, "definition": z, "pass": bool(np.isclose(x, y, atol=1e-10, rtol=1e-12))}
              for key, (x, y, z) in baseline.items()]
    assert all(row["pass"] for row in checks)
    pd.DataFrame(checks).to_csv(TABLES / "baseline_reproduction.csv", index=False)
    l2 = pd.read_csv(TABLES / "exposure_ledger.csv.gz", low_memory=False, dtype={"night": str, "noon_night": str})
    q = l2.loc[l2.quality]
    assert q.exposure_key.is_unique and q.l2_available.all()
    assert q.before_cutoff.all() and q.noon_night.nunique() == raw["quality_noon_nights"]
    assert not q.night.isin(["20251115", "20260111"]).any()
    assert sum(raw["exclusion_counts"].values()) == len(l2)
    assert raw["quality_exposures"]+sum(v for k,v in raw["exclusion_counts"].items() if k != "qualified") == len(l2)
    dup = pd.DataFrame(json.loads((INPUT / "duplicate_signature.json").read_text()))
    assert len(dup) == 273 and dup.a_sample_sha256.eq(dup.b_sample_sha256).all()
    assert dup.a_size.eq(dup.b_size).all() and dup.a_utc.eq(dup.b_utc).all()
    ext = pd.read_csv(TABLES / "extension_exposure_inventory.csv.gz", low_memory=False, dtype={"noon_night": str})
    assert pd.to_datetime(ext.loc[ext.entered_extension, "mid_utc"], format="mixed", utc=True).lt(CUTOFF).all()
    assert raw["extension_boundary_counts"] == {"before_cutoff": 235, "after_cutoff": 168}
    geom = pd.read_csv(TABLES / "geometric_coverage.csv")
    pix = pd.read_csv(DATA / "fig03_wcs_visits_icrs_nside512.csv.gz")
    low = geom.loc[geom.nside.eq(512)].iloc[0]
    assert np.isclose(len(pix)*hp.nside2pixarea(512, degrees=True), low.geometric_once_deg2)
    assert np.isclose(pix.ge3_same_night.sum()*hp.nside2pixarea(512, degrees=True), low.geometric_ge3_same_night_deg2)
    assert (geom.geometric_ge3_same_night_deg2 < geom.geometric_once_deg2).all()
    assert raw["qualified_with_verified_wcs"] == geom.exposures.iloc[0]
    known_sky = pd.read_csv(DATA / "fig03_known_ecliptic_nside128.csv")
    assert known_sky.known_detections.sum() == perf["known"]["qualified_detected"]
    cls = pd.read_csv(TABLES / "orbital_classes.csv")
    cls = cls.loc[cls.cohort.eq("quality")]
    assert cls.unique_objects.sum() == perf["known"]["qualified_identified"]
    assert cls.detections.sum() == perf["known"]["qualified_detected"]
    assert perf["minimum_altitude_bound_deg"] >= 0
    assert perf["spk_check_abs_difference_arcsec"]["max"] < .02
    assert sha(INPUT / "20260210_00000P2_link0035.gif") == "66c6275e59092cb4be29f36ca33772c1250e139412677bd97b91717b07392b02"
    assert sha(INPUT / "de432s.bsp") == "363f32e14f5255359ac32c4d38080cf28ab55564a5e16696a75f63394b666e9b"
    assert sha(INPUT / "survey_fov_footprints_with_visibility.fits") == "13f394877ce26e50e5cee5053a4ad845f412be4520421e46931b3e157d9e0225"
    figures = []
    for pdf in sorted((ROOT / "figures").glob("*.pdf")):
        png = pdf.with_suffix(".png")
        assert len(PdfReader(pdf).pages) == 1
        with Image.open(png) as image:
            dpi = tuple(round(k) for k in image.info["dpi"])
            assert dpi == (300, 300)
            figures.append({"stem": pdf.stem, "png_size": image.size, "png_dpi": dpi,
                            "pdf_sha256": sha(pdf), "png_sha256": sha(png)})
    assert len(figures) == 10
    qa = {"checked_utc": datetime.now(timezone.utc).isoformat(), "baseline_check_count": len(checks),
          "physical_exposure_identity_unique_in_quality_sample": True,
          "duplicate_pairs_checked_by_same_timestamp_and_four_segments": len(dup),
          "extension_boundary_partition_checked": True, "wcs_maps_and_area_consistent": True,
          "quality_known_sky_counts_close": True, "figure_pairs": figures}
    (TABLES / "validation.json").write_text(json.dumps(qa, indent=2) + "\n")
    print(f"Validated {len(checks)} baselines and {len(figures)} publication figure pairs")


if __name__ == "__main__":
    main()
