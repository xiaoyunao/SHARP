"""Local masked-catalog tracklet entry for the MP twilight speed profile."""
import argparse
from pathlib import Path
import subprocess
import sys
from .catalog import BASE
from .offline_twilight import load_config


def command(night, root, output, config, workers):
    if len(night)!=8 or not night.isdigit():raise ValueError('Night must be YYYYMMDD')
    return [sys.executable,str(BASE.parent/'heliolincrr/make_tracklet_linreproj.py'),night,
            '--root',str(Path(root).resolve()),'--outdir',str(Path(output).resolve()),
            '--vmin',str(config['vmin_arcsec_hour']),'--vmax',str(config['vmax_arcsec_hour']),
            '--dmag-max','1.0','--r-static','2.0','--min-repeat','2','--skip-common-area',
            '--nproc',str(workers),'--max-tracklets-per-group','100000']


def main():
    p=argparse.ArgumentParser(description='Local tracklets only. No orbit submission or server operations.')
    p.add_argument('--night',required=True);p.add_argument('--root',required=True)
    p.add_argument('--output',required=True);p.add_argument('--config');p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();source=Path(a.root)/a.night/'mask_gaia'
    if not any(source.glob('OBJ_MP_*_cat.fits.gz')):p.error('No masked local catalogs: '+str(source))
    result=subprocess.run(command(a.night,a.root,a.output,load_config(a.config),a.workers),check=False)
    raise SystemExit(result.returncode)


if __name__=='__main__':main()
