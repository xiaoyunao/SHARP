#!/usr/bin/env python3
"""Redraw SHARP image figures using the exact GOTTA cutout display recipe."""
from __future__ import annotations
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from PIL import Image, ImageSequence
from scipy.ndimage import median_filter

COLORS=["#ffcc33","#00b4d8","#fb5607","#80ed99"]

def zscale_limits(data):
    x=np.asarray(data,dtype=float);x=x[np.isfinite(x)]
    return tuple(float(v) for v in ZScaleInterval(contrast=.25,krej=2.5).get_limits(x))

def linear_stretch(data,vmin,vmax):
    x=np.clip((np.asarray(data,dtype=float)-vmin)/(vmax-vmin),0,1);x[~np.isfinite(x)]=0;return x

def flatten_background(cut,size=61):
    x=np.asarray(cut,dtype=float);finite=np.isfinite(x);fill=float(np.nanmedian(x[finite]));bg=median_filter(np.where(finite,x,fill),size=size,mode="nearest");y=x-bg;y[~finite]=np.nan;return y

def shared_limits(cuts):
    x=np.concatenate([c.ravel() for c in cuts]);x=x[np.isfinite(x)];med=np.median(x);mad=np.median(np.abs(x-med));sig=1.4826*mad
    return med-sig,med+7*sig

def cutout(data,x,y,size=180):
    half=size//2;cx=int(round(x))-1;cy=int(round(y))-1;return data[cy-half:cy+half,cx-half:cx+half]

def clear(ax):
    ax.set_xticks([]);ax.set_yticks([])
    for s in ax.spines.values():s.set_visible(False)

def label(ax,x,y,text,size=10,va="top"):
    ax.text(x,y,text,ha="left",va=va,color="white",fontsize=size,fontweight="bold",linespacing=1.15,bbox={"boxstyle":"round,pad=0.22","facecolor":"black","edgecolor":"none","alpha":.64})

def cross(ax,x,y,color,gap=11,length=17,lw=2.4):
    for x0,x1,y0,y1 in ((x-gap-length,x-gap,y,y),(x+gap,x+gap+length,y,y),(x,x,y-gap-length,y-gap),(x,x,y+gap,y+gap+length)):
        ax.plot([x0,x1],[y0,y1],color=color,lw=lw,solid_capstyle="butt")

def full_view(data,factor=6):
    ny,nx=data.shape;v=data[:ny-ny%factor,:nx-nx%factor].reshape(ny//factor,factor,nx//factor,factor).mean(axis=(1,3));return v,(nx,ny)

def fig01(root,image_path,matched_path):
    with fits.open(image_path,memmap=False) as h:data=np.asarray(h["IMG"].data,dtype=float);hdr=h["IMG"].header;mjd=float(hdr.get("MJD",hdr.get("MJD-OBS")))
    allm=pd.read_parquet(matched_path);allm=allm[(allm.night.astype(str)=="20251121")&allm.source_file.astype(str).str.contains("OBJ_MP_0940_0273")].sort_values("mag_aper4")
    src=pd.read_csv(root/"figure_data/fig01_sources.csv");selected=[]
    for s in src.itertuples(index=False):
        z=allm.iloc[np.argmin(np.hypot(allm.x_win_px-s.x_win_px,allm.y_win_px-s.y_win_px))];selected.append(z)
    view,(nx,ny)=full_view(data);fvmin,fvmax=zscale_limits(view);cuts=[flatten_background(cutout(data,r.x_win_px,r.y_win_px)) for r in selected];cvmin,cvmax=shared_limits(cuts)
    fig=plt.figure(figsize=(15.2,8.0),dpi=220);left=.015;bottom=.035;height=.93;gap=.018;cg=.018;lw=.475;cs=(1-left-lw-gap-.015-cg)/2;ch=(height-cg)/2;rl=left+lw+gap
    ax=fig.add_axes([left,bottom,lw,height]);ax.imshow(linear_stretch(view,fvmin,fvmax),cmap="gray",origin="lower",extent=[0,nx,0,ny],interpolation="nearest");clear(ax);label(ax,130,ny-150,f"OBJ_MP_0940_0273\nMJD = {mjd:.8f}",11.5)
    chosen={(round(r.x_win_px),round(r.y_win_px)) for r in selected}
    for r in allm.itertuples(index=False):
        if (round(r.x_win_px),round(r.y_win_px)) not in chosen and 58<r.x_win_px<nx-58 and 58<r.y_win_px<ny-58:
            ax.plot([r.x_win_px-58,r.x_win_px+58],[r.y_win_px,r.y_win_px],color="#2dd4bf",lw=1.2);ax.plot([r.x_win_px,r.x_win_px],[r.y_win_px-58,r.y_win_px+58],color="#2dd4bf",lw=1.2)
    for r,c in zip(selected,COLORS):
        cross(ax,r.x_win_px,r.y_win_px,c,gap=72,length=160,lw=1.5);name=str(r.asteroid_name);label(ax,min(max(r.x_win_px+155,90),nx-1450),min(max(r.y_win_px+95,170),ny-160),name,9.8,va="center")
    newrows=[]
    for i,(r,c,cut) in enumerate(zip(selected,COLORS,cuts)):
        col=i%2;row=i//2;x0=rl+col*(cs+cg);y0=bottom+(1-row)*(ch+cg);a=fig.add_axes([x0,y0,cs,ch]);a.imshow(linear_stretch(cut,cvmin,cvmax),cmap="gray",origin="lower",interpolation="nearest");clear(a);center=cut.shape[0]/2-.5;cross(a,center,center,c);num=pd.to_numeric(r.asteroid_number,errors="coerce");ident=f"({int(num)}) {r.asteroid_name}" if np.isfinite(num) else str(r.asteroid_name);dra=((r.obs_ra_win_deg-r.pred_ra_deg+180)%360-180)*np.cos(np.deg2rad(r.pred_dec_deg))*3600;ddec=(r.obs_dec_win_deg-r.pred_dec_deg)*3600;sep=np.hypot(dra,ddec);label(a,7,cut.shape[0]-7,f"{ident}\n$g_{{\\rm aper}}$ = {r.mag_aper4:.2f} mag\nO-C = {sep:.2f} arcsec",9)
        newrows.append({"image":"OBJ_MP_0940_0273.fits.gz","panel":chr(66+i),"asteroid_number":r.asteroid_number,"asteroid_name":r.asteroid_name,"epoch_mjd":r.epoch_mjd,"pred_ra_deg":r.pred_ra_deg,"pred_dec_deg":r.pred_dec_deg,"obs_ra_win_deg":r.obs_ra_win_deg,"obs_dec_win_deg":r.obs_dec_win_deg,"x_win_px":r.x_win_px,"y_win_px":r.y_win_px,"g_aper_mag":r.mag_aper4,"oc_separation_arcsec":sep})
    pd.DataFrame(newrows).to_csv(root/"figure_data/fig01_sources.csv",index=False)
    out=root/"figures/fig01_schmidt_field_and_asteroid_cutouts";fig.savefig(out.with_suffix(".png"),dpi=300,bbox_inches="tight",pad_inches=0);fig.savefig(out.with_suffix(".pdf"),bbox_inches="tight",pad_inches=0);plt.close(fig)

def fig06(root,gif_path):
    with Image.open(gif_path) as im:frames=[np.asarray(f.convert("RGB")) for f in ImageSequence.Iterator(im)]
    if len(frames)!=3:raise ValueError(f"Expected exactly 3 GIF frames, found {len(frames)}")
    fig=plt.figure(figsize=(14.55,4.85),dpi=100,facecolor="white")
    for i,frame in enumerate(frames):
        ax=fig.add_axes([i/3,0,1/3,1]);ax.imshow(frame,origin="upper",interpolation="nearest");ax.axis("off")
    out=root/"figures/fig06_candidate_three_frame_sequence";fig.savefig(out.with_suffix(".png"),dpi=300,bbox_inches=None,pad_inches=0);fig.savefig(out.with_suffix(".pdf"),bbox_inches=None,pad_inches=0);plt.close(fig)
    pd.DataFrame([{"frame":i+1,"source_gif":str(gif_path),"gif_frame_index":i,"width_px":frames[i].shape[1],"height_px":frames[i].shape[0]} for i in range(3)]).to_csv(root/"figure_data/fig06_astrometry.csv",index=False)
    (root/"figure_data/fig06_registration_check.json").write_text('{\n  "source": "author-specified three-frame review GIF",\n  "registration_check": "not re-estimated; frames retained exactly as supplied",\n  "added_overlays": false,\n  "layout": "one row by three columns"\n}\n')

def main(a):
    root=Path(a.root).resolve();fig01(root,Path(a.fig01_image),Path(a.known_matched));fig06(root,Path(a.fig06_gif))
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--root",required=True);p.add_argument("--fig01-image",required=True);p.add_argument("--known-matched",required=True);p.add_argument("--fig06-gif",required=True);main(p.parse_args())
