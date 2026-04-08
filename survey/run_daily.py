#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import shutil
from zoneinfo import ZoneInfo

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table

from .config import SchedulerConfig, ServerPaths, SiteConfig, load_ingest_state, save_ingest_state
from .history import create_history_template, load_history, save_history
from .history_update import apply_counts_to_history, count_l2_mp_exposures
from .io_utils import write_plan_json, write_plan_txt
from .scheduler import StripScheduler


def log_info(message: str) -> None:
    print(f"[INFO] {message}")


def log_section(title: str) -> None:
    print(f"[INFO] ---- {title} ----")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Daily server runner for the asteroid survey scheduler.")
    parser.add_argument("--date", default=None, help="Observing night date in YYYY-MM-DD. Defaults to today in Asia/Shanghai.")
    parser.add_argument("--root-dir", default="/pipeline/xiaoyunao/survey")
    parser.add_argument("--workspace", default="/pipeline/xiaoyunao/survey/runtime")
    parser.add_argument("--footprints", default="/pipeline/xiaoyunao/survey/footprints/survey_fov_footprints_with_visibility.fits")
    parser.add_argument("--processed-root", default="/processed1")
    parser.add_argument("--publish-dir", default="", help="Optional shared directory to receive the latest plan txt.")
    return parser.parse_args()


def _id_variants(field_id: str) -> set[str]:
    value = str(field_id).strip()
    variants = {value}
    if value.upper().startswith("MP") and value[2:].isdigit():
        variants.add(f"{int(value[2:]):04d}")
    elif value.isdigit():
        variants.add(f"MP{int(value):04d}")
    return variants


def prepare_footprints(footprints: Table) -> Table:
    out = footprints.copy()

    if "visible_days" not in out.colnames and "visible_day" in out.colnames:
        out.rename_column("visible_day", "visible_days")

    out["field_id"] = np.asarray([str(x).strip() for x in out["field_id"]]).astype("U16")

    if "field_index" not in out.colnames:
        out["field_index"] = np.arange(len(out), dtype=np.int32)

    if "dec_band" not in out.colnames:
        dec_key = np.round(np.asarray(out["center_dec"], dtype=float), 6)
        uniq_dec = np.unique(dec_key)
        dec_map = {value: idx for idx, value in enumerate(np.sort(uniq_dec))}
        out["dec_band"] = np.asarray([dec_map[value] for value in dec_key], dtype=np.int32)

    if "ra_band" not in out.colnames:
        ra_band = np.zeros(len(out), dtype=np.int32)
        center_ra = np.asarray(out["center_ra"], dtype=float)
        for dec_band in np.unique(np.asarray(out["dec_band"], dtype=int)):
            idx = np.where(np.asarray(out["dec_band"], dtype=int) == int(dec_band))[0]
            order = idx[np.argsort(center_ra[idx])]
            ra_band[order] = np.arange(order.size, dtype=np.int32)
        out["ra_band"] = ra_band

    if "ecliptic_lat_deg" not in out.colnames:
        coords = SkyCoord(
            ra=np.asarray(out["center_ra"], dtype=float) * u.deg,
            dec=np.asarray(out["center_dec"], dtype=float) * u.deg,
            frame="icrs",
        )
        out["ecliptic_lat_deg"] = coords.barycentrictrueecliptic.lat.deg.astype(float)

    return out


def remap_counts_to_history_ids(counts: dict[str, int], history_field_ids: np.ndarray) -> dict[str, int]:
    remapped: dict[str, int] = {}
    available = {str(x).strip() for x in history_field_ids}
    for field_id, count in counts.items():
        target = None
        for candidate in _id_variants(field_id):
            if candidate in available:
                target = candidate
                break
        if target is None:
            continue
        remapped[target] = remapped.get(target, 0) + int(count)
    return remapped


def resolve_night(date_arg: str | None) -> datetime.date:
    if date_arg:
        return datetime.strptime(date_arg, "%Y-%m-%d").date()
    return datetime.now(ZoneInfo("Asia/Shanghai")).date()


def main() -> None:
    args = parse_args()
    night = resolve_night(args.date)
    night_tag = night.strftime("%Y%m%d")
    publish_dir = Path(args.publish_dir) if args.publish_dir else None
    paths = ServerPaths(
        root_dir=Path(args.root_dir),
        workspace=Path(args.workspace),
        footprints_path=Path(args.footprints),
        processed_root=Path(args.processed_root),
        publish_dir=publish_dir,
    )
    paths.ensure()
    log_section("daily run start")
    log_info(f"night={night_tag}")
    log_info(f"root_dir={paths.root_dir}")
    log_info(f"workspace={paths.workspace}")
    log_info(f"footprints={paths.footprints_path}")
    log_info(f"processed_root={paths.processed_root}")
    log_info(f"publish_dir={paths.publish_dir if paths.publish_dir is not None else '<empty>'}")

    footprints = prepare_footprints(Table.read(paths.footprints_path))
    history_exists = paths.history_path.exists()
    history = load_history(paths.history_path, footprints) if history_exists else create_history_template(footprints)
    ingest_state = load_ingest_state(paths.ingest_state_path) if history_exists else {}

    updated = False
    update_counts: dict[str, int] = {}
    processed_l2_dir = paths.processed_root / night_tag / "L2"
    if history_exists:
        counts = remap_counts_to_history_ids(
            count_l2_mp_exposures(processed_l2_dir),
            np.asarray(history["field_id"]).astype("U16"),
        )
        already_ingested = ingest_state.get(night_tag, {})
        for field_id, count in counts.items():
            delta = count - int(already_ingested.get(field_id, 0))
            if delta > 0:
                update_counts[field_id] = delta
        if update_counts:
            history = apply_counts_to_history(history, update_counts, night_tag)
            merged_counts = dict(already_ingested)
            for field_id, count in counts.items():
                merged_counts[field_id] = int(count)
            ingest_state[night_tag] = merged_counts
            updated = True

    save_history(history, paths.history_path)
    save_ingest_state(paths.ingest_state_path, ingest_state)
    log_section("history")
    log_info(f"history_path={paths.history_path}")
    log_info(f"history_update={'applied' if updated else 'skipped'}")
    log_info(f"processed_l2={processed_l2_dir}")
    log_info(f"history_update_fields={len(update_counts)}")

    scheduler = StripScheduler(footprints, history, SiteConfig(), SchedulerConfig())
    plan = scheduler.build_plan(night)

    json_path = paths.plans_dir / f"{night_tag}_plan.json"
    txt_path = paths.plans_dir / f"{night_tag}_plan.txt"
    write_plan_json(plan, json_path)
    write_plan_txt(plan, txt_path)
    plots_dir = paths.nightly_plots_dir(night_tag)
    from .visualize_nightly import generate_plan_visualizations

    plot_manifest = generate_plan_visualizations(
        plan=plan,
        history=history,
        footprints=footprints,
        outdir=plots_dir,
        night_tag=night_tag,
    )
    log_section("products")
    log_info(f"plan_json={json_path} n_exp={len(plan)}")
    log_info(f"plan_txt={txt_path}")
    log_info(f"plots_dir={plots_dir}")
    log_info(f"plot_cycle_pngs={len(plot_manifest['cycle_pngs'])}")
    log_info(f"plot_summary={plot_manifest['summary_png']}")
    if plot_manifest["gif"]:
        log_info(f"plot_gif={plot_manifest['gif']}")
    published_path: Path | None = None
    latest_path: Path | None = None
    if paths.publish_dir is not None:
        published_path = paths.publish_dir / txt_path.name
        latest_path = paths.publish_dir / "current_plan.txt"
        shutil.copy2(txt_path, published_path)
        shutil.copy2(txt_path, latest_path)
    if published_path is not None and latest_path is not None:
        log_info(f"publish_txt={published_path}")
        log_info(f"publish_latest={latest_path}")
    log_section("daily run done")


if __name__ == "__main__":
    main()
