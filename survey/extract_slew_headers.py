#!/usr/bin/env python3
"""Read one night's complete raw OBJ header stream; never read image pixels.

Run with the survey server's default Python and redirect stdout to a JSON file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

from astropy.io import fits


KEYS = (
    "DATE", "OBS_DATE", "DATE-OBS", "MJD-OBS", "EXPOSURE", "EXPTIME",
    "EXPREQ", "OBJ_RA", "OBJ_DEC", "TEL_RA", "TEL_DEC", "HOURANG",
    "OBLST", "TEL_AZ", "TEL_EL", "OBJECT", "CHANNEL", "FILTER", "READMODE",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("night")
    parser.add_argument("--root", type=Path, default=Path("/raw1"))
    args = parser.parse_args()
    directory = args.root / args.night / "OBJ"
    rows = []
    for path in sorted(directory.glob("*.fits*")):
        if not (path.name.endswith(".fits") or path.name.endswith(".fits.gz")):
            continue
        stat = path.stat()
        row = {"path": str(path), "name": path.name, "size_bytes": stat.st_size,
               "mtime_ns": stat.st_mtime_ns}
        match = re.search(r"_(\d+)\.fits(?:\.gz)?$", path.name)
        row["sequence"] = int(match.group(1)) if match else None
        row["standard_mp"] = bool(re.fullmatch(r"OBJ_MP_\d+_\d+\.fits(?:\.gz)?", path.name))
        try:
            header = fits.getheader(path, 0)
            row["header"] = {key: header[key] for key in KEYS if key in header}
            row["comments"] = {key: header.comments[key] for key in KEYS if key in header}
            row["header_sha256"] = hashlib.sha256(header.tostring().encode("ascii")).hexdigest()
        except Exception as exc:
            row["error"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)
    result = {"night": args.night, "directory": str(directory),
              "retrieved_utc": datetime.now(timezone.utc).isoformat(),
              "method": "astropy.io.fits.getheader(path, 0); all OBJ FITS; no image pixels",
              "rows": rows}
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
