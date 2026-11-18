"""L1/L2 naming conversions, independent of observing-program prefixes."""
from pathlib import Path


def image_name(value):
    name=Path(str(value)).name
    return name.replace('_cat.fits.gz','.fits.gz').replace('_cat.fits','.fits')


def catalog_name(value):
    name=Path(str(value)).name
    if name.endswith(('_cat.fits','_cat.fits.gz')):return name
    if name.endswith('.fits.gz'):return name[:-8]+'_cat.fits.gz'
    if name.endswith('.fits'):return name[:-5]+'_cat.fits'
    raise ValueError('Expected a FITS image name')


def catalog_path(root,night,value):
    folder=Path(root)/night/'L2'
    standard=folder/catalog_name(value)
    return standard if standard.exists() else folder/Path(str(value)).name
