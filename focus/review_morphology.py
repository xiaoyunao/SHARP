#!/usr/bin/env python3
"""Inspect matched unsaturated stars and raw/overscan display; remote writes nothing."""
import io
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy.ndimage import shift
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from compare_focus import match

BASE=Path(__file__).resolve().parent
OUT=BASE/'results/20261111/morphology_review'
OUT.mkdir(parents=True,exist_ok=True)
names=['focus-5.84_0081','focus-5.86_0082','focus-5.88_0083','focus-5.90_0084','focus-5.92_0085']
cats=[np.genfromtxt(BASE/'results/20261111'/f'{n}_sources.csv',delimiter=',',names=True) for n in names]
ref=cats[2]
mappings=[]
keep=(ref['snr']>100)&(ref['snr']<500)&(ref['peak']<20000)
for cat in cats:
    a,b,_=match(ref,cat);mapping=np.full(len(ref),-1);mapping[a]=b;mappings.append(mapping)
    ok=np.zeros(len(ref),dtype=bool);ok[a]=(cat['snr'][b]>100)&(cat['snr'][b]<500)&(cat['peak'][b]<20000)
    keep &= ok
ids=np.flatnonzero(keep)
xy=np.array([np.column_stack((cat['x'][m[ids]],cat['y'][m[ids]])) for cat,m in zip(cats,mappings)])
payload={'names':names,'xy':xy.tolist()}
source=(BASE/'measure_focus.py').read_text()
remote='''
import types,sys,io,json,numpy as np
from astropy.io import fits
m=types.ModuleType('focus_measure');exec(SOURCE,m.__dict__)
p=PAYLOAD
allcuts=[];rawcuts=[]
for index,(name,points) in enumerate(zip(p['names'],p['xy'])):
    with fits.open('/raw1/20261111/OBJ/'+name+'.fits',memmap=False) as h:
        raw=h[0].data
        cuts=np.zeros((len(points),41,41),dtype='float32')
        for amp,xoff,yoff,tile in m.raw_tiles(raw):
            for j,(x,y) in enumerate(points):
                if int(x//1152)+int(y//4616)*8!=amp:continue
                ix,iy=int(round(x-xoff)),int(round(y-yoff))
                cuts[j]=tile[iy-20:iy+21,ix-20:ix+21]
                if index==2:
                    col=amp%8;row=amp//8
                    rx=col*1250+(27 if col<4 else 71)+ix
                    ry=(0 if row==0 else 4784)+iy
                    rawcuts.append((j,raw[ry-20:ry+21,rx-20:rx+21].astype('float32')))
        allcuts.append(cuts)
rawordered=np.stack([c for _,c in sorted(rawcuts)])
f=io.BytesIO();np.savez_compressed(f,cuts=np.stack(allcuts),raw=rawordered)
sys.stdout.buffer.write(f.getvalue())
'''.replace('SOURCE',repr(source)).replace('PAYLOAD',repr(payload))
r=subprocess.run(['ssh','-p','20093','smtpipeline@www.xinglong-naoc.cn','OPENBLAS_NUM_THREADS=1 /home/smtpipeline/Softwares/miniconda3/bin/python -B -'],input=remote.encode(),capture_output=True,check=True)
(OUT/'matched_cutouts.npz').write_bytes(r.stdout)
a=np.load(io.BytesIO(r.stdout));cuts=a['cuts'];raw=a['raw']
y,x=np.indices((41,41));radius=np.hypot(x-20,y-20);ann=(radius>15)&(radius<19)
aligned=np.zeros_like(cuts)
for f in range(len(names)):
    for j in range(len(ids)):
        cut=cuts[f,j]-np.median(cuts[f,j][ann])
        cx,cy=xy[f,j];cut=shift(cut,[-(cy-round(cy)),-(cx-round(cx))],order=3,mode='constant')
        aligned[f,j]=cut/cut[radius<12].sum()
cells=ref['cell'][ids].astype(int)
# Spatial medians on exactly the same stars, without image smoothing.
medians=np.array([[np.median(c['d50'][m[ids[cells==cell]]])*.286 for cell in range(16)] for c,m in zip(cats,mappings)])
np.savetxt(OUT/'matched_cell_d50.csv',np.column_stack((np.arange(16),np.bincount(cells,minlength=16),medians.T)),delimiter=',',header='cell,n,'+','.join(names),comments='')
fig,axes=plt.subplots(1,3,figsize=(12,4))
for ax,f in zip(axes,[1,2,3]):
    sc=ax.imshow(medians[f].reshape(4,4),origin='lower',vmin=.9,vmax=1.35,cmap='viridis',interpolation='nearest')
    for cell,v in enumerate(medians[f]):ax.text(cell%4,cell//4,f'{v:.2f}',ha='center',va='center',color='white')
    ax.set_title(names[f]);ax.set(xlabel='Detector x quarter',ylabel='Detector y quarter');ax.set_xticks(range(4));ax.set_yticks(range(4))
fig.colorbar(sc,ax=axes,label='Matched-star median D50 [arcsec]',shrink=.8);fig.suptitle(f'Identical {len(ids)} unsaturated stars in every frame');fig.savefig(OUT/'spatial_comparison.png',dpi=160,bbox_inches='tight');plt.close(fig)
# Display PSF stacks across x at mid field, identical per-panel peak normalization.
chosen=[8,9,11]
fig,axes=plt.subplots(3,5,figsize=(12,7.4))
stacks=[]
for row,cell in enumerate(chosen):
    ss=[]
    for f in range(5):
        stack=np.median(aligned[f,cells==cell],axis=0);ss.append(stack)
        core=stack[13:28,13:28];ax=axes[row,f]
        ax.imshow(core/core.max(),origin='lower',vmin=0,vmax=1,cmap='gray',interpolation='nearest')
        ax.contour(core/core.max(),levels=[.1,.3,.5,.8],colors=['#ff9830'],linewidths=.5)
        if row==0:ax.set_title(names[f].rsplit('_',1)[1])
        ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlabel(f'D50={medians[f,cell]:.3f}"')
        if f==0:ax.set_ylabel(['Left / cell 8','Inner / cell 9','Right / cell 11'][row])
    stacks.append(ss)
fig.suptitle('Same-star PSF stacks / full peak shown, no clipped white cores\n15 x 15 pixels (4.29 arcsec); linear scale; contours at 10/30/50/80%')
fig.tight_layout();fig.savefig(OUT/'psf_stacks.png',dpi=170);plt.close(fig)
# Raw vs corrected, and the former preview's percentile clipping.
fig,axes=plt.subplots(2,4,figsize=(12,6))
deltas=[]
for row,cell in enumerate([9,11]):
    candidates=np.flatnonzero(cells==cell);j=candidates[np.argsort(ref['snr'][ids[candidates]])[len(candidates)//2]]
    corr=cuts[2,j];unc=raw[j];sky0=np.median(unc[ann]);sky1=np.median(corr[ann]);peak=(corr-sky1).max()
    panels=[(unc,0,float(unc.max()),'Raw, display starts at 0'),(unc,sky0,sky0+peak,'Raw, local sky to peak'),(corr,sky1,sky1+peak,'Overscan removed, sky to peak'),(corr,sky1,float(np.percentile(corr[5:36,5:36],99)),'Previous preview: 99% clip')]
    for ax,(img,lo,hi,title) in zip(axes[row],panels):
        ax.imshow(img[10:31,10:31],origin='lower',vmin=lo,vmax=hi,cmap='gray',interpolation='nearest');ax.set_title(title,fontsize=9);ax.set_xticks([]);ax.set_yticks([])
    axes[row,0].set_ylabel(f'Cell {cell}\nx={xy[2,j,0]:.1f}, y={xy[2,j,1]:.1f}')
    delta=(unc-sky0)-(corr-sky1)
    deltas.append(dict(cell=cell,raw_sky=float(sky0),corrected_sky=float(sky1),peak_above_sky=float(peak),difference_rms=float(np.std(delta)),difference_over_peak=float(np.std(delta)/peak)))
fig.suptitle('0083: identical stars before/after overscan and display clipping');fig.tight_layout();fig.savefig(OUT/'raw_vs_corrected.png',dpi=170);plt.close(fig)
diff=(raw-np.median(raw[:,ann],axis=1)[:,None,None])-(cuts[2]-np.median(cuts[2][:,ann],axis=1)[:,None,None])
relative_rms=np.std(diff,axis=(1,2))/(cuts[2].max(axis=(1,2))-np.median(cuts[2][:,ann],axis=1))
report=dict(raw_corrected_relative_rms_percentiles=np.percentile(relative_rms,[16,50,84]).tolist(),selected_reference_indices=ids.tolist(),n_common=len(ids),cells_counts=np.bincount(cells,minlength=16).tolist(),names=names,cell_d50_arcsec=medians.tolist(),best_frame_by_cell=[names[i] for i in medians.argmin(axis=0)],display_examples=deltas,remote_writes=False)
(OUT/'review.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
