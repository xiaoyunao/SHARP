#!/usr/bin/env python3
"""Stream Xinglong evidence to stdout. Run on smtdata with its science Python."""
import argparse
import csv
import gzip
import json
import re
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS

warnings.simplefilter("ignore", AstropyWarning)
RAW = Path("/raw1")
PROCESSED = Path("/processed1")
PROJECT = Path("/pipeline/xiaoyunao")
REVIEW = PROJECT / "heliolincrr/review_packages"
BASE_END = "20260715"
END = "20261123"
CORNERS = [(0.5, 0.5), (9216.5, 0.5), (9216.5, 9232.5), (0.5, 9232.5)]
STRICT_OLD = re.compile(r"OBJ_MP_\d{4}_\d{4}\.fits(?:\.gz)?$")
RAW_NAME = re.compile(r"OBJ_.+\.fits(?:\.gz)?$")
HEADER = ("EXPSTA", "EXPEND", "OBS_DATE", "DATE-OBS", "EXPOSURE", "EXPTIME",
          "EXPREQ", "FILTER", "OBJECT", "IMAGETYP", "MJD", "N_WCS", "RA_RMS",
          "DEC_RMS", "RADESYS")


def stem(name):
    return re.sub(r"(_cat)?\.fits(?:\.gz)?$", "", name)


def read_json(path):
    return json.loads(path.read_bytes()) if path.is_file() else None


def read_metadata(path, catalog=False):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        h = fits.Header.fromfile(handle)
        if catalog:
            if h.get("NAXIS", 0):
                raise ValueError("expected empty L2 primary HDU")
            h = fits.Header.fromfile(handle)
    info = {k: h[k] for k in HEADER if k in h}
    if catalog:
        wcs = WCS(h, naxis=2)
        if not wcs.has_celestial:
            raise ValueError("missing celestial WCS")
        info["corners"] = wcs.all_pix2world(CORNERS, 1).tolist()
        info["center"] = wcs.all_pix2world([[4608.5, 4616.5]], 1)[0].tolist()
    return info


def scan(night):
    base = PROCESSED / night
    l4 = base / "L4"
    invalid = read_json(base / "ASTEROID_DATA_INVALID.json")
    plan_path = PROJECT / "survey/runtime/plans" / f"{night}_plan.json"
    plan = read_json(plan_path)
    plans = []
    if isinstance(plan, list):
        for row in plan:
            plans.append({k: row[k] for k in ("field_id", "filter", "mode", "block",
                          "start_utc", "exptime_s", "ra", "dec") if k in row})
    manifest = l4 / f"{night}_file_manifest.txt"
    science = read_json(l4 / f"{night}_science_manifest.json")
    accepted = {stem(Path(r["file"]).name): r for r in (science or {}).get("accepted", [])}
    files = sorted({stem(Path(line.strip()).name) for line in manifest.read_text().splitlines()
                    if line.strip()}) if manifest.is_file() else []
    raw_dir = RAW / night / "OBJ"
    raw = {}
    for path in sorted(raw_dir.glob("OBJ_*.fits*")):
        if RAW_NAME.fullmatch(path.name):
            raw.setdefault(stem(path.name), path)
    l2 = {}
    for path in sorted((base / "L2").glob("OBJ_*_cat.fits*")):
        l2.setdefault(stem(path.name), path)
    rows = []
    for key, path in raw.items():
        if night <= BASE_END and STRICT_OLD.fullmatch(path.name):
            continue
        cat = l2.get(key)
        row = dict(stem=key, raw=str(path), l2=str(cat) if cat else None)
        for kind, p in (("raw", path), ("l2", cat)):
            if p is None:
                continue
            try:
                row[f"{kind}_header"] = read_metadata(p, catalog=(kind == "l2"))
            except (OSError, ValueError, EOFError, TypeError) as exc:
                row[f"{kind}_error"] = str(exc)
        rows.append(row)
    status_path = l4 / f"{night}_v2_pipeline_status.json"
    status = read_json(status_path)
    return {
        "night": night, "invalid": invalid, "rows": rows, "raw_file_count": len(raw),
        "l2_only": sorted(set(l2) - set(raw)), "manifest": files, "science_accepted": accepted,
        "known_exists": (l4 / f"{night}_matched_asteroids.fits").is_file(),
        "unknown_exists": (l4 / f"{night}_unknown_links.json").is_file(),
        "manifest_exists": manifest.is_file(), "status": status, "plan": plans,
        "plan_path": str(plan_path) if plans else None,
    }


def candidate(night):
    l4 = PROCESSED / night / "L4"
    rp = REVIEW / night
    raw = read_json(l4 / f"{night}_unknown_links.json")
    screened = read_json(l4 / f"{night}_unknown_links_screened.json")
    morphology = read_json(l4 / f"{night}_unknown_morphology.json")
    package = read_json(rp / f"{night}_unknown_review_manifest.json")
    summary = read_json(PROJECT / "data/heliolincrr" / night / "analysis" /
                        f"{night}_single_night_summary.json")
    human_path = rp / f"{night}_unknown_review.csv"
    human = None
    human_mtime = None
    if human_path.is_file():
        with human_path.open(encoding="utf-8-sig") as handle:
            human = list(csv.DictReader(handle))
        human_mtime = datetime.fromtimestamp(human_path.stat().st_mtime, timezone.utc).isoformat()
    if all(item is None for item in (raw, screened, morphology, package, summary, human)):
        return None
    candidate_fields = ("trk_sub", "linkage_id", "n_obs", "image_names", "objids", "mjds",
                        "ras_deg", "decs_deg", "lin_speed_arcsec_per_day")
    screen_fields = ("trk_sub", "status", "decision", "reason", "reasons")
    replies = {}
    for p in (l4 / f"{night}_unknown_mpc_reply.txt", rp / f"{night}_unknown_mpc_reply.txt"):
        if p.is_file():
            replies[str(p)] = p.read_text(errors="replace")
    return {
        "night": night, "raw": [{k: x[k] for k in candidate_fields if k in x}
                                for x in raw] if isinstance(raw, list) else None,
        "screened_ids": [x.get("trk_sub") for x in screened] if isinstance(screened, list) else None,
        "morphology": [{k: x[k] for k in screen_fields if k in x}
                       for x in (morphology or {}).get("candidates", [])],
        "morphology_counts": (morphology or {}).get("counts"),
        "package_n": (package or {}).get("n_catalog_rows"),
        "package_obs": (package or {}).get("review_ades_rows"),
        "human": human, "human_mtime": human_mtime, "replies": replies,
        "pre_gate_count": (summary or {}).get("unknown_fit_ok_catalog", {}).get("count"),
    }


def known(night):
    path = PROCESSED / night / "L4" / f"{night}_matched_asteroids.fits"
    if not path.is_file():
        return None
    with fits.open(path, memmap=True) as hdus:
        table = hdus[1].data
        records = [
            [str(name).strip(), int(number) if number > 0 else 0, str(source).strip(),
             float(epoch), int(objid)]
            for name, number, source, epoch, objid in zip(
                *(table[k] for k in ("name", "number", "source_file", "epoch", "objID")))]
    return {"night": night, "source": str(path), "rows": records}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=["exposures", "known", "candidates"])
    parser.add_argument("--end", default=END)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    nights = sorted({p.name for root in (RAW, PROCESSED, REVIEW) for p in root.iterdir()
                     if p.is_dir() and re.fullmatch(r"\d{8}", p.name)
                     and "20251115" <= p.name <= args.end})
    fn = {"exposures": scan, "known": known, "candidates": candidate}[args.kind]
    if args.kind == "known":
        nights = [n for n in nights if BASE_END < n <= args.end and n != "20261117"]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for row in pool.map(fn, nights):
            if row is not None:
                print(json.dumps(row, ensure_ascii=False, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
