#!/usr/bin/env python3
"""Check magnitude repeatability in six historical field-night triplets."""
from itertools import combinations

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

from exposures import INPUT, TABLES, read


def directions(s):
    return SkyCoord(s.RA_Win.to_numpy()*u.deg, s.DEC_Win.to_numpy()*u.deg)


def repeat_pairs(frame):
    result = []
    for (field, night), group in frame.groupby(["field_id", "night"]):
        exposures = {name: d.reset_index(drop=True) for name, d in group.groupby("path")}
        calibrated = {}
        for path, sources in exposures.items():
            coords = directions(sources)
            _, nearest, _ = coords.match_to_catalog_sky(coords, nthneighbor=2)
            keep = nearest.arcsec > 8
            s = sources.loc[keep].reset_index(drop=True)
            calibrated[path] = s, directions(s)
        for p, q in combinations(calibrated, 2):
            a, ac = calibrated[p]
            b, bc = calibrated[q]
            other, sep, _ = ac.match_to_catalog_sky(bc)
            back, _, _ = bc.match_to_catalog_sky(ac)
            mutual = (back[other] == np.arange(len(a))) & (sep.arcsec <= .65)
            first, second = a.loc[mutual], b.iloc[other[mutual]]
            a_mag = first.Mag_Aper4.to_numpy()
            b_mag = second.Mag_Aper4.to_numpy()
            err = np.hypot(first.MagErr_Aper4.to_numpy(), second.MagErr_Aper4.to_numpy())
            result.append(pd.DataFrame(dict(
                field_id=field, night=night, exposure_a=p, exposure_b=q, g_aper_pair=(a_mag+b_mag)/2,
                dx_mag=a_mag-b_mag, expected_pair_sigma_mag=err, match_sep_arcsec=sep.arcsec[mutual],
                obj_id_a=first.objID.to_numpy(), obj_id_b=second.objID.to_numpy(),
            )))
    return pd.concat(result, ignore_index=True)


def main():
    d = read("stellar_sample.csv.gz", dtype={"field_id": str, "night": str})
    if d.path.nunique() != 18:
        raise ValueError("Expected 18 selected historical catalogs")
    headers = d[["path", "sha256", "FILTER", "BAND", "CALI_REF", "MLIM", "FWHM", "SKY",
                 "CCDZP", "GAIN", "RON", "APER_REF", "APERSIZE", "APER08", "FLUX_S"]].drop_duplicates()
    adu = 5*np.sqrt((headers.SKY*headers.GAIN + headers.RON**2)*np.pi*(headers.FWHM/2)**2)/headers.GAIN
    headers["MLIM_replayed"] = headers.CCDZP-2.5*np.log10(adu)
    headers["header_minus_replayed"] = headers.MLIM-headers.MLIM_replayed
    assert headers.header_minus_replayed.abs().max() < .006
    headers.to_csv(TABLES / "sample_photometric_header_and_depth_check.csv", index=False)
    pair = repeat_pairs(d)
    pair.to_csv(TABLES / "stellar_repeated_measurement_pairs.csv.gz", index=False)
    bins = []
    for lower in np.arange(13, 21.5, .5):
        z = pair.loc[pair.g_aper_pair.ge(lower) & pair.g_aper_pair.lt(lower+.5)]
        q25, q50, q75 = z.dx_mag.quantile([.25, .5, .75])
        bins.append({"g_aper_low": lower, "g_aper_high": lower+.5, "star_pairs": len(z),
                     "median_signed_delta_mag": q50, "empirical_per_epoch_sigma_mag": .7413*(q75-q25)/np.sqrt(2),
                     "median_catalog_per_epoch_sigma_mag": z.expected_pair_sigma_mag.median()/np.sqrt(2)})
    pd.DataFrame(bins).to_csv(TABLES / "stellar_repeatability_vs_magnitude.csv", index=False)


if __name__ == "__main__":
    main()
