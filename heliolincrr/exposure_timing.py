"""Exposure midpoint from shutter times, never an implicit 30-second exposure."""
import math
from astropy.time import Time


def header_mid_mjd(header):
    def timestamp(value):
        try:
            return Time(str(value).strip(), format='isot', scale='utc')
        except (ValueError, TypeError):
            return None

    start = timestamp(header.get('EXPSTA'))
    end = timestamp(header.get('EXPEND'))
    if start is not None and end is not None:
        seconds = (end - start).sec
        if not math.isfinite(seconds) or seconds <= 0:
            return None
        return float((start + (end - start) / 2).mjd)
    # Explicit midpoint wins over ambiguous MJD. MJD-OBS denotes start here.
    for key in ('MJD-MID', 'MIDMJD'):
        try:
            value = float(header[key])
            if math.isfinite(value):
                return value
        except (KeyError, ValueError, TypeError):
            pass
    if start is None:
        start = timestamp(header.get('DATE-OBS'))
    if start is None:
        for key in ('MJD-OBS', 'MJD'):
            try:
                value = float(header[key])
                if math.isfinite(value):
                    start = Time(value, format='mjd', scale='utc')
                    break
            except (KeyError, ValueError, TypeError):
                pass
    for key in ('EXPOSURE', 'EXPTIME'):
        try:
            exposure = float(header[key])
            if start is not None and math.isfinite(exposure) and exposure > 0:
                from astropy.time import TimeDelta
                return float((start + TimeDelta(exposure / 2, format='sec')).mjd)
        except (KeyError, ValueError, TypeError):
            pass
    return None
