from __future__ import annotations

import json
import math
from collections import defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.table import Table
from astropy.time import Time

from .config import SchedulerConfig, SiteConfig


STATE_VERSION = 1
DEFAULT_START_NIGHT = "20260624"
DEFAULT_MAX_AGE_DAYS = 10
DEFAULT_REQUIRED_NIGHTS = 2
DEFAULT_OBS_PER_NIGHT = 5


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def split_semicolon(value: object) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def night_to_date(night: str) -> date:
    return datetime.strptime(night, "%Y%m%d").date()


def days_between(a: str, b: str) -> int:
    return (night_to_date(b) - night_to_date(a)).days


def load_state(path: Path, start_night: str = DEFAULT_START_NIGHT) -> dict[str, Any]:
    if not path.exists():
        return {
            "version": STATE_VERSION,
            "created_utc": utc_now(),
            "start_night": start_night,
            "sources": {},
            "ingested_nights": {},
        }
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Follow-up state is not a JSON object: {path}")
    data.setdefault("version", STATE_VERSION)
    data.setdefault("start_night", start_night)
    data.setdefault("sources", {})
    data.setdefault("ingested_nights", {})
    return data


def save_state(path: Path, state: dict[str, Any]) -> None:
    state["updated_utc"] = utc_now()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)


def source_key(night: str, trk_sub: str, linkage_id: object = "") -> str:
    trk = str(trk_sub or "").strip()
    if trk:
        return f"{night}:{trk}"
    return f"{night}:link{int(linkage_id):04d}"


def source_tag(source: dict[str, Any]) -> str:
    return str(source.get("trk_sub") or f"link{int(source.get('linkage_id', 0)):04d}")


def target_name_for(source: dict[str, Any], field_id: str) -> str:
    tag = source_tag(source)
    field_text = str(field_id).strip()
    if field_text.isdigit():
        field_text = f"MP{int(field_text):04d}"
    return f"MP_FU_{tag}_{field_text}"


def parse_link_detections(row: dict[str, Any], night: str) -> list[dict[str, Any]]:
    mjds = [float(x) for x in split_semicolon(row.get("mjds"))]
    ras = [float(x) for x in split_semicolon(row.get("ras_deg"))]
    decs = [float(x) for x in split_semicolon(row.get("decs_deg"))]
    image_names = split_semicolon(row.get("image_names"))
    objids = split_semicolon(row.get("objids"))
    out: list[dict[str, Any]] = []
    for idx, (mjd, ra, dec) in enumerate(zip(mjds, ras, decs)):
        det = {
            "night": night,
            "mjd": float(mjd),
            "ra_deg": float(ra) % 360.0,
            "dec_deg": float(dec),
        }
        if idx < len(image_names):
            det["image_name"] = image_names[idx]
        if idx < len(objids):
            det["objid"] = objids[idx]
        out.append(det)
    return out


def dedupe_detections(detections: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[tuple[str, str, str]] = set()
    out: list[dict[str, Any]] = []
    for det in sorted(detections, key=lambda x: float(x.get("mjd", 0.0))):
        key = (
            str(det.get("night", "")),
            str(det.get("image_name", "")),
            str(det.get("objid", "")),
        )
        fallback = ("", "", f"{float(det.get('mjd', 0.0)):.9f}:{float(det.get('ra_deg', 0.0)):.8f}:{float(det.get('dec_deg', 0.0)):.8f}")
        use_key = key if any(key) else fallback
        if use_key in seen:
            continue
        seen.add(use_key)
        out.append(det)
    return out


def fit_motion_model(detections: list[dict[str, Any]]) -> dict[str, float]:
    if len(detections) < 2:
        raise ValueError("Need at least two detections for linear follow-up prediction")
    dets = sorted(detections, key=lambda x: float(x["mjd"]))
    mjd = np.asarray([float(d["mjd"]) for d in dets], dtype=float)
    ra = np.asarray([float(d["ra_deg"]) for d in dets], dtype=float)
    dec = np.asarray([float(d["dec_deg"]) for d in dets], dtype=float)
    finite = np.isfinite(mjd) & np.isfinite(ra) & np.isfinite(dec)
    mjd, ra, dec = mjd[finite], ra[finite], dec[finite]
    if mjd.size < 2 or np.unique(mjd).size < 2:
        raise ValueError("Need at least two finite unique-MJD detections")

    ra_unwrap = np.rad2deg(np.unwrap(np.deg2rad(ra)))
    t0 = float(mjd[-1])
    ra0 = float(ra_unwrap[-1])
    dec0 = float(dec[-1])
    cos_dec = max(1.0e-6, math.cos(math.radians(dec0)))
    dt = mjd - t0
    x = (ra_unwrap - ra0) * cos_dec * 3600.0
    y = (dec - dec0) * 3600.0
    vx, x0 = np.polyfit(dt, x, 1)
    vy, y0 = np.polyfit(dt, y, 1)
    return {
        "epoch_mjd": t0,
        "ra0_deg": ra0 % 360.0,
        "dec0_deg": dec0,
        "x0_arcsec": float(x0),
        "y0_arcsec": float(y0),
        "vx_arcsec_per_day": float(vx),
        "vy_arcsec_per_day": float(vy),
        "fit_n": int(mjd.size),
        "fit_updated_utc": utc_now(),
    }


def predict_position(source: dict[str, Any], obstime: Time) -> tuple[float, float]:
    model = source.get("motion_model") or fit_motion_model(source.get("detections", []))
    dt_days = float(obstime.mjd) - float(model["epoch_mjd"])
    dec = float(model["dec0_deg"]) + (float(model["y0_arcsec"]) + float(model["vy_arcsec_per_day"]) * dt_days) / 3600.0
    cos_dec = max(1.0e-6, math.cos(math.radians(float(model["dec0_deg"]))))
    ra = float(model["ra0_deg"]) + (float(model["x0_arcsec"]) + float(model["vx_arcsec_per_day"]) * dt_days) / (3600.0 * cos_dec)
    return ra % 360.0, dec


def load_json_rows(path: Path) -> list[dict[str, Any]]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError(f"Expected JSON list: {path}")
    return [row for row in rows if isinstance(row, dict)]


def ingest_reviewed_sources(
    state: dict[str, Any],
    processed_root: Path,
    review_root: Path,
    start_night: str,
    through_night: str | None = None,
    only_night: str | None = None,
) -> dict[str, Any]:
    sources = state.setdefault("sources", {})
    ingested = state.setdefault("ingested_nights", {})
    stats = {"scanned_nights": 0, "new_sources": 0, "updated_sources": 0, "skipped_before_start": 0}
    if not review_root.exists():
        return stats
    for night_dir in sorted(review_root.iterdir()):
        if not night_dir.is_dir() or not (night_dir.name.isdigit() and len(night_dir.name) == 8):
            continue
        night = night_dir.name
        if only_night and night != only_night:
            continue
        if night < start_night:
            stats["skipped_before_start"] += 1
            continue
        if through_night and night > through_night:
            continue
        masked = processed_root / night / "L4" / f"{night}_unknown_links_submit_masked.json"
        if not masked.exists():
            continue
        signature = {
            "path": str(masked),
            "mtime_ns": masked.stat().st_mtime_ns,
            "size": masked.stat().st_size,
        }
        if ingested.get(night) == signature:
            continue
        stats["scanned_nights"] += 1
        for row in load_json_rows(masked):
            trk_sub = str(row.get("trk_sub") or row.get("trkSub") or "").strip()
            key = source_key(night, trk_sub, row.get("linkage_id", ""))
            dets = parse_link_detections(row, night)
            if len(dets) < 2:
                continue
            item = sources.get(key)
            if item is None:
                item = {
                    "source_id": key,
                    "origin_night": night,
                    "trk_sub": trk_sub,
                    "linkage_id": row.get("linkage_id"),
                    "status": "active",
                    "created_utc": utc_now(),
                    "detections": dets,
                    "associated_links": [],
                    "attempts": [],
                    "successful_nights": [],
                }
                sources[key] = item
                stats["new_sources"] += 1
            else:
                item["detections"] = dedupe_detections(list(item.get("detections", [])) + dets)
                if item.get("status") in {"abandoned", "complete"}:
                    pass
                stats["updated_sources"] += 1
            try:
                item["motion_model"] = fit_motion_model(item.get("detections", []))
            except ValueError as exc:
                item["status"] = "prediction_failed"
                item["status_reason"] = str(exc)
        ingested[night] = signature
    return stats


def footprint_radius_deg(row: Any) -> float:
    colnames = set(getattr(row, "colnames", getattr(getattr(row, "_table", None), "colnames", [])))
    center = SkyCoord(float(row["center_ra"]) * u.deg, float(row["center_dec"]) * u.deg)
    seps: list[float] = []
    for idx in range(1, 5):
        ra_key = f"corner_ra_{idx}"
        dec_key = f"corner_dec_{idx}"
        if ra_key in colnames and dec_key in colnames:
            corner = SkyCoord(float(row[ra_key]) * u.deg, float(row[dec_key]) * u.deg)
            seps.append(float(center.separation(corner).deg))
    return max(seps) if seps else 1.6


def prepare_footprint_lookup(footprints: Table) -> dict[str, Any]:
    coords = SkyCoord(np.asarray(footprints["center_ra"], dtype=float) * u.deg, np.asarray(footprints["center_dec"], dtype=float) * u.deg)
    radii = np.asarray([footprint_radius_deg(row) for row in footprints], dtype=float)
    return {"table": footprints, "coords": coords, "radii_deg": radii}


def select_followup_field(lookup: dict[str, Any], ra_deg: float, dec_deg: float, margin_deg: float = 0.05) -> tuple[dict[str, Any] | None, str]:
    target = SkyCoord(float(ra_deg) * u.deg, float(dec_deg) * u.deg)
    seps = target.separation(lookup["coords"]).deg
    inside = seps <= (lookup["radii_deg"] + margin_deg)
    if not np.any(inside):
        idx = int(np.argmin(seps))
        return None, f"outside_footprints nearest={str(lookup['table'][idx]['field_id'])} sep_deg={float(seps[idx]):.3f}"
    candidates = np.where(inside)[0]
    idx = int(candidates[np.argmin(seps[candidates])])
    row = lookup["table"][idx]
    colnames = set(getattr(row, "colnames", getattr(getattr(row, "_table", None), "colnames", [])))
    return {
        "field_id": str(row["field_id"]).strip(),
        "center_ra": float(row["center_ra"]),
        "center_dec": float(row["center_dec"]),
        "dec_band": int(row["dec_band"]) if "dec_band" in colnames else 0,
        "ra_band": int(row["ra_band"]) if "ra_band" in colnames else 0,
        "sep_deg": float(seps[idx]),
    }, ""


def is_observable(ra_deg: float, dec_deg: float, obstime: Time, site: SiteConfig, config: SchedulerConfig) -> bool:
    location = EarthLocation(lat=site.lat_deg * u.deg, lon=site.lon_deg * u.deg, height=site.height_m * u.m)
    altaz = SkyCoord(ra_deg * u.deg, dec_deg * u.deg).transform_to(AltAz(obstime=obstime, location=location))
    return bool(float(altaz.alt.deg) >= float(config.min_altitude_deg))


def plan_slots(plan: list[dict[str, Any]]) -> list[dict[str, Any]]:
    slots: list[dict[str, Any]] = []
    if not plan:
        return slots
    block_rows: list[dict[str, Any]] = []
    current_key: tuple[Any, Any] | None = None
    for idx, row in enumerate(plan):
        key = (row.get("cycle_id"), row.get("repeat"))
        if current_key is None:
            current_key = key
        if key != current_key:
            last = block_rows[-1]
            slots.append({"after_index": idx - 1, "cycle_id": current_key[0], "repeat": current_key[1], "base_start_utc": last["start_utc"]})
            block_rows = []
            current_key = key
        block_rows.append(row)
    if block_rows and current_key is not None:
        last = block_rows[-1]
        slots.append({"after_index": len(plan) - 1, "cycle_id": current_key[0], "repeat": current_key[1], "base_start_utc": last["start_utc"]})
    return slots


def source_has_attempt_for_night(source: dict[str, Any], night: str) -> bool:
    return any(str(a.get("night")) == night for a in source.get("attempts", []))


def reconcile_observed_followups(
    state: dict[str, Any],
    processed_root: Path,
    plan_night: str,
    obs_per_night: int = DEFAULT_OBS_PER_NIGHT,
) -> dict[str, Any]:
    stats = {"success": 0, "failed": 0, "checked": 0}
    for source in state.get("sources", {}).values():
        tag = source_tag(source)
        for attempt in source.get("attempts", []):
            if str(attempt.get("status")) != "planned":
                continue
            night = str(attempt.get("night"))
            if night >= plan_night:
                continue
            l2_dir = processed_root / night / "L2"
            count = 0
            if l2_dir.exists():
                count = len(list(l2_dir.glob(f"*MP_FU_{tag}*.fits*")))
            attempt["observed_l2_count"] = count
            attempt["reconciled_utc"] = utc_now()
            stats["checked"] += 1
            if count >= obs_per_night:
                attempt["status"] = "success"
                source.setdefault("successful_nights", [])
                if night not in source["successful_nights"]:
                    source["successful_nights"].append(night)
                stats["success"] += 1
            else:
                attempt["status"] = "failed"
                attempt["reason"] = "observed_less_than_required"
                stats["failed"] += 1
        if len(source.get("successful_nights", [])) >= DEFAULT_REQUIRED_NIGHTS:
            source["status"] = "complete"
            source["completed_utc"] = utc_now()
    return stats


def expire_old_sources(state: dict[str, Any], plan_night: str, max_age_days: int = DEFAULT_MAX_AGE_DAYS) -> dict[str, int]:
    stats = {"abandoned": 0}
    for source in state.get("sources", {}).values():
        if source.get("status") != "active":
            continue
        origin = str(source.get("origin_night"))
        if days_between(origin, plan_night) >= max_age_days:
            source["status"] = "abandoned"
            source["status_reason"] = f"age_days>={max_age_days}"
            source["abandoned_utc"] = utc_now()
            stats["abandoned"] += 1
    return stats


def active_sources_for_night(
    state: dict[str, Any],
    plan_night: str,
    required_nights: int = DEFAULT_REQUIRED_NIGHTS,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for source in state.get("sources", {}).values():
        if source.get("status") != "active":
            continue
        if str(source.get("origin_night")) >= plan_night:
            continue
        if len(source.get("successful_nights", [])) >= required_nights:
            continue
        if source_has_attempt_for_night(source, plan_night):
            continue
        out.append(source)
    return out


def build_followup_rows_for_slot(
    slot: dict[str, Any],
    source: dict[str, Any],
    field: dict[str, Any],
    pred_ra: float,
    pred_dec: float,
    exptime_s: float,
) -> dict[str, Any]:
    target_name = target_name_for(source, field["field_id"])
    return {
        "field_id": target_name,
        "cycle_id": slot["cycle_id"],
        "repeat": slot["repeat"],
        "mode": "followup",
        "start_utc": slot["base_start_utc"],
        "ra": float(field["center_ra"]),
        "dec": float(field["center_dec"]),
        "dec_band": int(field.get("dec_band", 0)),
        "ra_band": int(field.get("ra_band", 0)),
        "strip_score": 0.0,
        "strip_ra_min_deg": float(field["center_ra"]),
        "strip_ra_max_deg": float(field["center_ra"]),
        "strip_dec_count": 1,
        "strip_ra_count": 1,
        "exptime_s": float(exptime_s),
        "moon_phase_fraction": None,
        "moon_sep_deg": None,
        "followup_source_id": source["source_id"],
        "followup_origin_night": source.get("origin_night"),
        "followup_origin_trk_sub": source.get("trk_sub"),
        "followup_origin_linkage_id": source.get("linkage_id"),
        "followup_selected_field_id": str(field["field_id"]),
        "followup_field_center_sep_deg": float(field["sep_deg"]),
        "followup_pred_ra": float(pred_ra),
        "followup_pred_dec": float(pred_dec),
        "followup_slot_cycle_id": slot["cycle_id"],
        "followup_slot_repeat": slot["repeat"],
    }


def schedule_followups(
    plan: list[dict[str, Any]],
    state: dict[str, Any],
    footprints: Table,
    plan_night: str,
    site: SiteConfig | None = None,
    config: SchedulerConfig | None = None,
    obs_per_night: int = DEFAULT_OBS_PER_NIGHT,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    site = site or SiteConfig()
    config = config or SchedulerConfig()
    lookup = prepare_footprint_lookup(footprints)
    slots = plan_slots(plan)
    slot_rows: dict[int, list[dict[str, Any]]] = defaultdict(list)
    inserted_by_source: dict[str, list[dict[str, Any]]] = {}
    rejected: dict[str, str] = {}
    active_at_start = active_sources_for_night(state, plan_night)

    for source in active_at_start:
        candidates: list[dict[str, Any]] = []
        for slot_idx, slot in enumerate(slots):
            base_time = Time(slot["base_start_utc"], format="isot", scale="utc")
            pred_ra, pred_dec = predict_position(source, base_time)
            field, reason = select_followup_field(lookup, pred_ra, pred_dec)
            if field is None:
                rejected[source["source_id"]] = reason
                continue
            if not is_observable(float(field["center_ra"]), float(field["center_dec"]), base_time, site, config):
                rejected[source["source_id"]] = "not_observable_at_slot"
                continue
            row = build_followup_rows_for_slot(slot, source, field, pred_ra, pred_dec, config.exposure_time_s)
            row["followup_slot_index"] = slot_idx
            candidates.append(row)
            if len(candidates) >= obs_per_night:
                break
        if len(candidates) >= obs_per_night:
            inserted_by_source[source["source_id"]] = candidates[:obs_per_night]
            for row in candidates[:obs_per_night]:
                slot_rows[int(row["followup_slot_index"])].append(row)
            source.setdefault("attempts", []).append(
                {
                    "night": plan_night,
                    "status": "planned",
                    "planned_observations": obs_per_night,
                    "planned_utc": utc_now(),
                    "targets": [
                        {
                            "target_name": row["field_id"],
                            "selected_field_id": row["followup_selected_field_id"],
                            "slot_index": row["followup_slot_index"],
                            "pred_ra": row["followup_pred_ra"],
                            "pred_dec": row["followup_pred_dec"],
                        }
                        for row in candidates[:obs_per_night]
                    ],
                }
            )
        else:
            source.setdefault("attempts", []).append(
                {
                    "night": plan_night,
                    "status": "not_planned",
                    "planned_observations": 0,
                    "planned_utc": utc_now(),
                    "reason": rejected.get(source["source_id"], "insufficient_observable_slots"),
                    "candidate_slots": len(candidates),
                }
            )

    if not inserted_by_source:
        return plan, {
            "active_sources": len(active_at_start),
            "scheduled_sources": 0,
            "inserted_rows": 0,
            "slots": len(slots),
        }

    step_s = float(config.exposure_time_s + config.readout_overhead_s + config.slew_overhead_s)
    new_plan: list[dict[str, Any]] = []
    shift_s = 0.0
    slot_by_after = {slot["after_index"]: idx for idx, slot in enumerate(slots)}
    for idx, row in enumerate(plan):
        shifted = dict(row)
        shifted["start_utc"] = (Time(row["start_utc"], format="isot", scale="utc") + shift_s * u.s).isot
        new_plan.append(shifted)
        slot_idx = slot_by_after.get(idx)
        if slot_idx is None:
            continue
        inserts = slot_rows.get(slot_idx, [])
        if not inserts:
            continue
        t = Time(shifted["start_utc"], format="isot", scale="utc") + step_s * u.s
        for insert in inserts:
            out = dict(insert)
            out["start_utc"] = t.isot
            new_plan.append(out)
            t = t + step_s * u.s
            shift_s += step_s

    return new_plan, {
        "active_sources": len(active_at_start),
        "scheduled_sources": len(inserted_by_source),
        "inserted_rows": sum(len(v) for v in inserted_by_source.values()),
        "slots": len(slots),
    }


def angular_sep_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    a = SkyCoord(float(ra1) * u.deg, float(dec1) * u.deg)
    b = SkyCoord(float(ra2) * u.deg, float(dec2) * u.deg)
    return float(a.separation(b).arcsec)


def associate_unknown_links(
    state: dict[str, Any],
    night: str,
    catalog_path: Path,
    max_sep_arcsec: float = 30.0,
) -> dict[str, Any]:
    stats = {"catalog_links": 0, "matches": 0, "updated_sources": 0}
    if not catalog_path.exists():
        stats["missing_catalog"] = str(catalog_path)
        return stats
    rows = load_json_rows(catalog_path)
    stats["catalog_links"] = len(rows)
    updated_sources: set[str] = set()
    active = [s for s in state.get("sources", {}).values() if s.get("status") == "active" and str(s.get("origin_night")) < night]
    for row in rows:
        dets = parse_link_detections(row, night)
        if not dets:
            continue
        best_source: dict[str, Any] | None = None
        best_med = float("inf")
        for source in active:
            seps: list[float] = []
            for det in dets:
                pred_ra, pred_dec = predict_position(source, Time(float(det["mjd"]), format="mjd", scale="utc"))
                seps.append(angular_sep_arcsec(pred_ra, pred_dec, float(det["ra_deg"]), float(det["dec_deg"])))
            med = float(np.median(seps)) if seps else float("inf")
            if med < best_med:
                best_med = med
                best_source = source
        if best_source is None or best_med > max_sep_arcsec:
            continue
        assoc_key = f"{night}:{row.get('trk_sub') or row.get('linkage_id')}"
        existing = {str(x.get("association_key")) for x in best_source.get("associated_links", [])}
        if assoc_key in existing:
            continue
        best_source.setdefault("associated_links", []).append(
            {
                "association_key": assoc_key,
                "night": night,
                "trk_sub": row.get("trk_sub") or row.get("trkSub"),
                "linkage_id": row.get("linkage_id"),
                "median_sep_arcsec": best_med,
                "matched_utc": utc_now(),
            }
        )
        best_source["detections"] = dedupe_detections(list(best_source.get("detections", [])) + dets)
        best_source["motion_model"] = fit_motion_model(best_source["detections"])
        updated_sources.add(str(best_source["source_id"]))
        stats["matches"] += 1
    stats["updated_sources"] = len(updated_sources)
    return stats
