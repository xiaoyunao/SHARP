#!/usr/bin/env python3
"""Read-only SMT focus diagnostics. No production imports; outputs only to --out."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
import sep
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import median_filter, binary_dilation


def raw_tiles(data):
    """SMT 16 amplifier layout, matching inspected ccd9k.py row correction."""
    for i in range(16):
        row, col = divmod(i, 8)
        y0, x0 = row * 4700, col * 1250
        block = data[y0:y0+4700, x0:x0+1250].astype(np.float32)
        oc = block[:, -71:] if col < 4 else block[:, :71]
        bias = np.median(oc, axis=1)
        smooth = median_filter(bias, size=50)
        _, _, std = sigma_clipped_stats(bias-smooth, sigma=3)
        bad = binary_dilation(np.abs(bias-smooth) > 7*std, iterations=2)
        bias[bad] = smooth[bad]
        block -= bias[:, None]
        ys = slice(0, 4616) if row == 0 else slice(84, 4700)
        xs = slice(27, 1179) if col < 4 else slice(71, 1223)
        yield i, col*1152, row*4616, np.ascontiguousarray(block[ys, xs])


def measure_tile(data, amp, xoff, yoff, saturation=45000, snr_min=30):
    mask = np.ascontiguousarray(~np.isfinite(data) | (data > saturation))
    data = np.nan_to_num(data, copy=True)
    back = sep.Background(data, mask=mask, bw=128, bh=128, fw=3, fh=3)
    signal = np.ascontiguousarray(data-back.back(), dtype=np.float32)
    rms = np.ascontiguousarray(back.rms(), dtype=np.float32)
    # Low threshold and no restrictive roundness cut: donuts and astigmatism matter.
    objects, seg = sep.extract(signal, 3.0, err=rms, mask=mask, minarea=12,
                               deblend_nthresh=32, deblend_cont=0.01,
                               segmentation_map=True)
    rows, previews = [], []
    reasons = {}
    def reject(reason):
        reasons[reason] = reasons.get(reason, 0)+1
    for j, obj in enumerate(objects):
        if obj['flag'] != 0:
            reject('extraction_flag'); continue
        x, y, a, b = [float(obj[k]) for k in ['x','y','a','b']]
        if b < 0.65 or a > 35 or not np.isfinite(a+b):
            reject('nonstellar_scale'); continue
        # Six second-moment sigmas; validate convergence with a larger aperture.
        radius = max(12., 6*a)
        outer = int(np.ceil(radius*1.6))
        ix, iy = int(round(x)), int(round(y))
        if min(ix,iy,data.shape[1]-1-ix,data.shape[0]-1-iy) <= outer:
            reject('edge'); continue
        sl = np.s_[iy-outer:iy+outer+1, ix-outer:ix+outer+1]
        cut = signal[sl].copy()
        sg, bad = seg[sl], mask[sl]
        yy, xx = np.indices(cut.shape)
        cx, cy = x-(ix-outer), y-(iy-outer)
        rr = np.hypot(xx-cx, yy-cy)
        other = (sg != 0) & (sg != j+1)
        core = rr < radius
        if np.any(bad & core):
            reject('saturated_or_bad'); continue
        # Do not fill holes with invented stellar light.
        if np.count_nonzero(other & core) / np.count_nonzero(core) > .005:
            reject('neighbor'); continue
        sky_mask = (rr > radius*1.25) & (rr < radius*1.6) & ~other & ~bad
        if sky_mask.sum() < 60:
            reject('sky'); continue
        _, sky, noise = sigma_clipped_stats(cut[sky_mask], sigma=3, maxiters=3)
        cut -= sky
        bad = np.ascontiguousarray(bad | other)
        noise = max(float(noise),float(np.median(rms[sl])),1e-3)
        flux, error, flag = sep.sum_circle(cut,[cx],[cy],[radius],err=noise,mask=bad,subpix=5)
        # Unknown gain: empirical background plus Poisson term at conservative 1 e/ADU.
        err = np.sqrt(error[0]**2 + max(flux[0],0) + (np.pi*radius**2*noise/np.sqrt(sky_mask.sum()))**2)
        snr = flux[0]/err
        if snr < snr_min or flag[0] & 16:
            reject('snr'); continue
        radii, rf = sep.flux_radius(cut,[cx],[cy],[radius],[.5,.8],mask=bad,subpix=5)
        ext, _ = sep.flux_radius(cut,[cx],[cy],[radius*1.2],[.5,.8],mask=bad,subpix=5)
        r50,r80 = radii[0]
        convergence = float(np.max(np.abs(ext[0]/radii[0]-1)))
        if not (0.8 < r50 < r80 < .8*radius) or convergence > .1:
            reject('growth_not_converged'); continue
        # Diagnostic radial FWHM only for centrally peaked profiles; never rank donuts by it.
        bins = np.arange(0, int(radius)+1)
        profile = np.array([np.mean(cut[(rr>=k)&(rr<k+1)&~bad])
                            if np.any((rr>=k)&(rr<k+1)&~bad) else np.nan for k in bins[:-1]])
        peak_bin = int(np.nanargmax(profile))
        hole = bool(peak_bin > 1 and profile[0] < .8*np.nanmax(profile))
        fwhm = float('nan')
        if not hole:
            crossing = np.flatnonzero(profile[1:] <= profile[0]/2)+1
            if len(crossing):
                k=int(crossing[0]); fwhm=2*(k-.5+(profile[k-1]-profile[0]/2)/(profile[k-1]-profile[k]))
        rows.append(dict(x=x+xoff,y=y+yoff,amp=amp,cell=int((y+yoff)//2308)*4+int((x+xoff)//2304),
                         flux=float(flux[0]),snr=float(snr),d50=2*r50,d80=2*r80,
                         ellipticity=1-b/a,fwhm=fwhm,donut=int(hole),growth_change=convergence,
                         peak=float(obj['peak'])))
        if len(previews)<2 and snr>80:
            previews.append((cut,rows[-1]))
    return rows, previews, dict(detected=len(objects),accepted=len(rows),rejected=reasons)


def analyze(path, out, snr_min=30, saturation=45000):
    started=time.time(); path=Path(path); out=Path(out);out.mkdir(parents=True,exist_ok=True)
    before=path.stat()
    with fits.open(path,memmap=False) as hdul:
        hdu=next(h for h in hdul if h.data is not None and h.data.ndim==2)
        header=hdu.header.copy(); data=hdu.data
        band=str(header.get('FILTER',header.get('CHANNEL',''))).strip().lower()
        if band!='r': raise ValueError('Not r band: '+band)
        if data.shape==(9400,10000): tiles=raw_tiles(data)
        elif data.shape==(9232,9216):
            tiles=((i,(i%8)*1152,(i//8)*4616,
                    np.ascontiguousarray(data[(i//8)*4616:(i//8+1)*4616,(i%8)*1152:(i%8+1)*1152],dtype=np.float32)) for i in range(16))
        else: raise ValueError('Unexpected detector geometry: '+str(data.shape))
        rows=[];previews=[];stats=[]
        for amp,xoff,yoff,tile in tiles:
            r,p,s=measure_tile(tile,amp,xoff,yoff,snr_min=snr_min,saturation=saturation)
            rows.extend(r);previews.extend(p);stats.append(s)
            print(json.dumps(dict(amp=amp,**s)),flush=True)
    after=path.stat()
    if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
        raise RuntimeError('Input changed during reading')
    if not rows: raise RuntimeError('No usable sources')
    name=path.name.replace('.fits.gz','').replace('.fits','')
    with (out/(name+'_sources.csv')).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    scale=float(header.get('PIXSCAL1',1))
    summary=dict(file=str(path),name=name,band=band,exposure=header.get('EXPOSURE'),
                 utc=header.get('EXPSTA',header.get('DATE')),pixel_scale=scale,n_sources=len(rows),
                 input_size=before.st_size,input_mtime_ns=before.st_mtime_ns,
                 snr_min=snr_min,saturation_adu=saturation,amplifiers=stats,seconds=time.time()-started,
                 code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    for key in ['d50','d80','ellipticity','fwhm','growth_change','donut']:
        vals=np.array([r[key] for r in rows]); finite=vals[np.isfinite(vals)]
        summary[key+'_median']=float(np.median(finite)) if len(finite) else None
        if key in ['d50','d80','fwhm']: summary[key+'_arcsec']=summary[key+'_median']*scale if len(finite) else None
    summary['cells']={str(c):dict(n=sum(r['cell']==c for r in rows),
                                  d50=float(np.median([r['d50'] for r in rows if r['cell']==c])))
                      for c in sorted({r['cell'] for r in rows})}
    (out/(name+'_summary.json')).write_text(json.dumps(summary,indent=2,allow_nan=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(15,4.6))
    sc=axes[0].scatter([r['x'] for r in rows],[r['y'] for r in rows],c=[r['d50']*scale for r in rows],s=3,cmap='viridis',vmin=np.percentile([r['d50']*scale for r in rows],5),vmax=np.percentile([r['d50']*scale for r in rows],95))
    axes[0].set(xlim=(0,9216),ylim=(0,9232),title='All usable sources',xlabel='x [pixel]',ylabel='y [pixel]');fig.colorbar(sc,ax=axes[0],label='D50 [arcsec]')
    axes[1].hist([r['d50']*scale for r in rows],bins=50);axes[1].set(xlabel='D50 [arcsec]',ylabel='Sources',title=f'N = {len(rows)}')
    axes[2].scatter([r['snr'] for r in rows],[r['d50']*scale for r in rows],s=3);axes[2].set(xscale='log',xlabel='S/N',ylabel='D50 [arcsec]',title='Check brightness dependence')
    fig.suptitle(name);fig.tight_layout();fig.savefig(out/(name+'_diagnostic.png'),dpi=150);plt.close(fig)
    fig,axes=plt.subplots(4,8,figsize=(16,8))
    for ax,(cut,row) in zip(axes.ravel(),previews):
        lim=15;cy,cx=np.array(cut.shape)//2;stamp=cut[cy-lim:cy+lim+1,cx-lim:cx+lim+1]
        ax.imshow(stamp,origin='lower',cmap='gray',vmin=0,vmax=np.max(stamp),interpolation='nearest');ax.set_title(f"A{row['amp']} D50={row['d50']*scale:.2f}\"",fontsize=8)
    for ax in axes.ravel():ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle(name+' / representative sources (same pixel scale; full peak)');fig.tight_layout();fig.savefig(out/(name+'_stamps.png'),dpi=130);plt.close(fig)
    print('SUMMARY '+json.dumps(summary),flush=True)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('images',nargs='+');parser.add_argument('--out',required=True);parser.add_argument('--snr-min',type=float,default=30)
    parser.add_argument('--saturation',type=float,default=45000,help='Conservative ceiling after overscan; SMT r plateau observed near 51000 ADU')
    args=parser.parse_args()
    for path in args.images:analyze(path,args.out,args.snr_min,args.saturation)
