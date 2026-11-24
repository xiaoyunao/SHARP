#!/usr/bin/env python3
"""Read calibrated L1 pixels on smtdata and stream the render inputs as NPZ."""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

BASE = Path("/processed1")


def midtime(header):
    start = Time(header["EXPSTA"], scale="utc")
    end = Time(header["EXPEND"], scale="utc")
    return (start + (end - start) / 2).isot


def stamps(night, internal_id):
    l4 = BASE / night / "L4" / f"{night}_unknown_links.json"
    candidate = next(row for row in json.loads(l4.read_bytes()) if row["trk_sub"] == internal_id)
    ras = np.fromstring(candidate["ras_deg"], sep=";")
    decs = np.fromstring(candidate["decs_deg"], sep=";")
    images = candidate["image_names"].split(";")
    num = 154
    scale = 1.15
    target = WCS(naxis=2)
    target.wcs.crpix = [(num+1) / 2] * 2
    target.wcs.crval = [np.mean(ras), np.mean(decs)]
    target.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    target.wcs.cdelt = [-scale / 3600, scale / 3600]
    yy, xx = np.indices((num, num))
    world = target.all_pix2world(np.stack([xx, yy], axis=-1).reshape(-1, 2), 0)
    position = target.all_world2pix(np.column_stack([ras, decs]), 0)
    frames, sources = [], []
    for filename in images:
        path = BASE / night / "L1" / filename
        with fits.open(path) as hdus:
            image = hdus["IMG"]
            original = WCS(image.header, naxis=2)
            x, y = original.all_world2pix(world, 0).T
            x0, x1 = int(x.min())-2, int(x.max())+3
            y0, y1 = int(y.min())-2, int(y.max())+3
            assert x0 >= 0 and x1 < image.shape[1] and y0 >= 0 and y1 < image.shape[0]
            cutout = image.data[y0:y1, x0:x1].astype(np.float32)
            frames.append(map_coordinates(cutout, [y-y0, x-x0], order=1).reshape(num, num))
            sources.append(dict(path=str(path), midpoint_utc=midtime(image.header),
                                exptime_s=image.header["EXPOSURE"], filter=image.header["FILTER"],
                                section=[x0, x1, y0, y1]))
    meta = dict(catalog=str(l4), candidate=candidate, projection="north-up east-left TAN",
                scale=scale, frames=sources, target_xy=position.tolist())
    return np.array(frames), meta


if __name__ == "__main__":
    asteroid, source = stamps("20260901", "00001mX")
    np.savez_compressed(sys.stdout.buffer, asteroid=asteroid,
                        metadata=np.array(json.dumps({"asteroid": source})))
