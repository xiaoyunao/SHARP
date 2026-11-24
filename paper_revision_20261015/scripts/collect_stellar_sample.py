#!/usr/bin/env python3
"""Read the selected frozen L2 catalogs on smtdata and stream calibrated columns."""
import csv
import gzip
import hashlib
import io
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

COLS = ["objID", "RA_Win", "DEC_Win", "X_Win", "Y_Win",
        "Mag_Aper4", "MagErr_Aper4", "Flux_Aper4", "FluxErr_Aper4", "Flag", "Type", "Flag_ISO", "Class_ANN"]
HEADER = ["FILTER", "BAND", "CALI_REF", "MLIM", "FWHM", "SKY", "CCDZP", "GAIN", "RON",
          "APER_REF", "APERSIZE", "APER08", "FLUX_S"]


def main():
    manifest = csv.DictReader(open(sys.argv[1])) if len(sys.argv) > 1 else csv.DictReader(io.StringIO(MANIFEST))
    keys = ["night", "field_id", "path", "sha256", *HEADER, *COLS]
    writer = csv.writer(sys.stdout)
    writer.writerow(keys)
    for row in manifest:
        path = Path(row["path"])
        file = path.stat()
        mtime = Time(row["mtime_utc"].replace("+00:00", ""), format="isot", scale="utc").unix
        if abs(file.st_mtime-mtime) > .001 or file.st_size != int(row["size_bytes"]):
            raise ValueError(f"Frozen catalog size/mtime differ: {path}")
        with path.open("rb") as file:
            encoded = file.read()
        digest = hashlib.sha256(encoded).hexdigest()
        with fits.open(io.BytesIO(gzip.decompress(encoded)), memmap=False) as hdus:
            h, d = hdus[1].header, hdus[1].data
            values = [d[key] for key in COLS]
            keep = (d["Flag"] == 0) & (d["Flag_ISO"] == 0) & (d["Type"] == 1) & (d["Class_ANN"] > .8) & (d["Mag_Aper4"] >= 12) & (d["Mag_Aper4"] < 22)
            keep &= (d["MagErr_Aper4"] > 0) & (d["Flux_Aper4"] > 0) & (d["FluxErr_Aper4"] > 0)
            keep &= ((d["X_Win"] > 500) & (d["X_Win"] < 8716) &
                     (d["Y_Win"] > 500) & (d["Y_Win"] < 8732))
            rows = zip(*(val[keep] for val in values))
            prefix = [row["night"], row["field_id"], path, digest, *[h.get(key, "") for key in HEADER]]
            writer.writerows([*prefix, *data] for data in rows)


if __name__ == "__main__":
    main()
