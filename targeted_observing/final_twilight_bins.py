"""Nominal twilight elongation counts and HEALPix unions at exposure midpoints."""
import json
from collections import Counter

import astropy.units as u
import healpy as hp
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, get_sun
from astropy.time import Time

from .annual_coverage import pixels
from .annual_sources import dates
from .continuation_report import OUT, ROOT, sha
from .geometry import Sky, coord, setup_iers, unix
from .catalog import BASE

EDGES = [45, 55, 65, 75, 90]
LABELS = ['<45', '45–55', '55–65', '65–75', '75–90', '>90']
NSIDE = 256


def bins(elongation):
    return np.searchsorted(EDGES, np.minimum(elongation, np.nextafter(90., 0.)),
                           side='right') + (elongation > 90)


def summarize():
    spec = json.loads((BASE/'annual_defaults.json').read_text())
    setup_iers(BASE/'cache'/spec['iers_file'])
    nights = dates(spec['start_date'], spec['stop_date_exclusive'])
    results, sources = [], {}
    for strategy in ['twilight_5', 'twilight_morning10_aligned']:
        records = []
        for night in nights:
            path = ROOT/strategy/f'{night}.json'
            sources[f'{strategy}/{path.name}'] = sha(path)
            plan = json.loads(path.read_text())
            sky = Sky(plan['config'])
            for block in plan['blocks']:
                for row in block['rows']:
                    records.append((night, block['block_id'], row['field_id'],
                                    unix(row['mid_utc']), row['ra_deg'], row['dec_deg']))
        exposures = pd.DataFrame(records, columns=['date', 'block', 'field', 'mid', 'ra', 'dec'])
        exposures['side'] = exposures.block.str[:2]
        for side, group in exposures.groupby('side'):
            unique_pixels = np.zeros((len(LABELS), hp.nside2npix(NSIDE)), dtype=bool)
            for start in range(0, len(group), 256):
                batch = group.iloc[start:start+256]
                frame = AltAz(obstime=Time(batch.mid.to_numpy(), format='unix'),
                              location=sky.location, pressure=0*u.hPa)
                sun = get_sun(frame.obstime).transform_to(frame)
                centre = coord(batch.ra, batch.dec).transform_to(frame)
                exposures.loc[batch.index, 'bin'] = bins(centre.separation(sun).deg)
                ids = [pixels(ra, dec, NSIDE) for ra, dec in zip(batch.ra, batch.dec)]
                padded = np.full((max(map(len, ids)), len(ids)), -1, dtype=int)
                for column, field_pixels in enumerate(ids):
                    padded[:len(field_pixels), column] = field_pixels
                ra, dec = hp.pix2ang(NSIDE, padded.clip(min=0), lonlat=True)
                altaz = coord(ra, dec).transform_to(frame[None, :])
                pixel_bins = bins(altaz.separation(sun[None, :]).deg)
                good = padded >= 0
                unique_pixels[pixel_bins[good], padded[good]] = True
            subset = exposures.loc[group.index]
            qualified = Counter()
            for _, field in subset.groupby(['date', 'block', 'field']):
                if len(field) >= 4 and np.ptp(field.mid) >= 720:
                    middle = (field.mid.min() + field.mid.max()) / 2
                    qualified[int(field.loc[(field.mid-middle).abs().idxmin(), 'bin'])] += 1
            for index, label in enumerate(LABELS):
                results.append(dict(strategy=strategy, side=side, elongation_deg=label,
                    nominal_exposures=int((subset.bin == index).sum()),
                    nominal_exposure_denominator=len(subset),
                    nominal_qualified_field_windows=qualified[index],
                    nominal_qualified_window_denominator=sum(qualified.values()),
                    nominal_pixel_union_deg2=unique_pixels[index].sum()*hp.nside2pixarea(NSIDE, degrees=True),
                    nside=NSIDE))
    pd.DataFrame(results).to_csv(OUT/'twilight_elongation_bins.csv', index=False)
    (OUT/'twilight_elongation_inputs.json').write_text(json.dumps(sources, indent=2)+'\n')


if __name__ == '__main__':
    summarize()
