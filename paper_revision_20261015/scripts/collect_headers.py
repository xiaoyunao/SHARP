#!/usr/bin/env python3

"""Stream a read-only exposure inventory from the Xinglong filesystem."""
import argparse
import csv
import io
import zlib
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import date, timedelta
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS

START, END = date(2025, 11, 14), date(2026, 10, 17)
LAST_OLD = "20260715"
STRICT = re.compile(r"OBJ_MP_(\d{4})_(\d{4})\.fits")
WCS_CORNER = [(0.5, 0.5), (9216.5, 0.5), (9216.5, 9232.5), (0.5, 9232.5)]
WCS_KEYS = ("mjd", "obs_date", "exposure", "filter", "mlim", "seeing", "ccdzp",
            "zp_aper4", "radesys", "califref", "caliref", "cali_ref")
RAW_KEYS = ("obs_date", "expsta", "exposure", "expreq", "exptime", "filter",
            "time-sys", "timesys", "instrument", "instrume", "naxis1", "naxis2",
            "tel_ra", "tel_dec", "obj_ra", "obj_dec", "expstop")
CORNERS = [f"corner_{i}_{a}" for i in range(1, 5) for a in ("ra", "dec")]
FIELDS = ["directory_night", "kind", "source_file", "bytes", "mtime_ns", "error",
          "l1_present", "known_night_product_present", "unknown_night_product_present",
          "known_night_mtime_ns", "unknown_night_mtime_ns", "known_mask15_present",
          "l2_present", "strict_name", *RAW_KEYS,
          *[f"l2_{c}" for c in (*WCS_KEYS, *CORNERS)]]


def product_status(path):
    try:
        return path.stat().st_mtime_ns
    except FileNotFoundError:
        return ""


def header_buffer(path):
    size = 2880*24
    with path.open("rb", buffering=0) as handle:
        data = handle.read(size)
    if path.name.endswith(".gz"):
        data = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(data, size)
    return io.BytesIO(data)


def raw_header(path):
    return fits.Header.fromfile(header_buffer(path))


def scan_night(night):
    raw_dir = Path("/raw1", night, "OBJ")
    base = Path("/processed1", night)
    l1_names = {p.name for p in (base / "L1").glob("OBJ_MP_*")}
    l2_names = {p.name for p in (base / "L2").glob("OBJ_MP_*")}
    known = product_status(base / "L4" / f"{night}_matched_asteroids.fits")
    unknown = product_status(base / "L4" / f"{night}_unknown_links.fits")
    mask15 = product_status(base / "L4" / f"{night}_matched_asteroids_mask15.fits")
    rows = []
    for path in sorted(raw_dir.glob("OBJ_MP_*.fits")):
        strict = STRICT.fullmatch(path.name) is not None
        if not strict and night <= LAST_OLD:
            continue
        stem = path.stem
        l2 = next((n for n in (f"{stem}_cat.fits.gz", f"{stem}_cat.fits") if n in l2_names), None)
        row = dict(directory_night=night, kind="raw", source_file=path.name,
                   bytes=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns,
                   strict_name=strict)
        row.update(l1_present=f"{stem}.fits.gz" in l1_names or f"{stem}.fits" in l1_names,
                   l2_present=l2 is not None, known_night_product_present=bool(known),
                   unknown_night_product_present=bool(unknown), known_night_mtime_ns=known,
                   unknown_night_mtime_ns=unknown, known_mask15_present=bool(mask15))
        try:
            h = raw_header(path)
            row.update({k: h.get(k.upper(), "") for k in RAW_KEYS})
        except (OSError, ValueError, EOFError) as exc:
            row["error"] = str(exc)
        rows.append(row)
    if night <= LAST_OLD:
        for name in sorted(l2_names):
            if not STRICT.fullmatch(name.split("_cat.fits")[0] + ".fits"):
                continue
            path = base / "L2" / name
            row = dict(directory_night=night, kind="l2_wcs", source_file=name,
                       bytes=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns)
            try:
                with header_buffer(path) as handle:
                    primary = fits.Header.fromfile(handle)
                    if primary.get("NAXIS", 0):
                        raise ValueError("unexpected nonempty L2 primary HDU")
                    h = fits.Header.fromfile(handle)
                row.update({f"l2_{k}": h.get(k.upper(), "") for k in WCS_KEYS})
                corners = WCS(h, naxis=2).all_pix2world(WCS_CORNER, 1)
                row.update(dict(zip([f"l2_{c}" for c in CORNERS], corners.ravel())))
            except (OSError, ValueError, EOFError) as exc:
                row["error"] = str(exc)
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=date.fromisoformat, default=START)
    parser.add_argument("--end", type=date.fromisoformat, default=END)
    parser.add_argument("--workers", type=int, default=4)
    options = parser.parse_args()
    writer = csv.DictWriter(sys.stdout, FIELDS)
    writer.writeheader()
    nights = [(options.start + timedelta(days=i)).strftime("%Y%m%d") for i in range((options.end - options.start).days + 1)]
    with ProcessPoolExecutor(max_workers=options.workers) as pool:
        for rows in pool.map(scan_night, nights):
            writer.writerows(rows)
            sys.stdout.flush()


if __name__ == "__main__":
    main()
