#!/usr/bin/env python3
"""Build real-image Figures 1 and 6 from frozen SHARP calibrated exposures."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.ndimage import map_coordinates
from common import *

def zscale(a):
    x=a[np.isfinite(a)];lo,hi=np.percentile(x,[3,99.7]);return np.arcsinh(np.clip((a-lo)/(hi-lo),0,None)*8)

def load_image(path):
    with fits.open(path,memmap=False) as h:
        return np.asarray(h["IMG"].data,dtype=np.float32),WCS(h["IMG"].header),dict(h["IMG"].header)

def reproj(img,w_in,w_out,shape):
    yy,xx=np.indices(shape);ra,dec=w_out.pixel_to_world_values(xx,yy);xi,yi=w_in.world_to_pixel_values(ra,dec)
    return map_coordinates(img,[yi,xi],order=1,mode="constant",cval=np.nan)

def make_fig01(root,img_path,matched_path):
    img,w,h=load_image(img_path); m=pd.read_parquet(matched_path);m=m[(m.night.astype(str)=="20251121")&(m.source_file.astype(str).str.contains("OBJ_MP_0940_0273"))].copy()
    sel=m.sort_values("mag_aper4").iloc[[0,2,5,12]].copy(); sel["panel"]=["B","C","D","E"]
    src=sel[["panel","asteroid_number","asteroid_name","epoch_mjd","pred_ra_deg","pred_dec_deg","obs_ra_win_deg","obs_dec_win_deg","x_win_px","y_win_px","mag_aper4"]].copy();src.insert(0,"image",img_path.name);src.to_csv(root/"figure_data/fig01_sources.csv",index=False)
    fig=plt.figure(figsize=(7.1,5.5));gs=fig.add_gridspec(2,4,height_ratios=[3.3,1.25],hspace=.12,wspace=.08);ax=fig.add_subplot(gs[0,:])
    ds=img[::16,::16];ax.imshow(zscale(ds),origin="lower",cmap="gray",vmin=0,vmax=2.7);ax.set(xticks=[],yticks=[]);ax.text(.01,.96,"A",transform=ax.transAxes,va="top",weight="bold",color="white",fontsize=11)
    ax.annotate("N",xy=(.94,.86),xytext=(.94,.70),xycoords="axes fraction",textcoords="axes fraction",color="white",ha="center",arrowprops=dict(arrowstyle="-|>",color="white",lw=1.2));ax.annotate("E",xy=(.84,.84),xytext=(.94,.84),xycoords="axes fraction",textcoords="axes fraction",color="white",va="center",arrowprops=dict(arrowstyle="-|>",color="white",lw=1.2));ax.plot([.05,.15],[.07,.07],transform=ax.transAxes,color="white",lw=2);ax.text(.10,.09,"1°",transform=ax.transAxes,ha="center",color="white",fontsize=8)
    half=55
    for j,r in enumerate(sel.itertuples(index=False)):
        a=fig.add_subplot(gs[1,j]);x,y=int(round(r.x_win_px)),int(round(r.y_win_px));cut=img[y-half:y+half,x-half:x+half];a.imshow(zscale(cut),origin="lower",cmap="gray",vmin=0,vmax=2.7);a.add_patch(Circle((half,half),9,fill=False,ec=ORANGE,lw=1.3));a.set(xticks=[],yticks=[]);a.text(.04,.92,r.panel,transform=a.transAxes,weight="bold",color="white");a.text(.5,.03,f"({int(r.asteroid_number)}) {r.asteroid_name}",transform=a.transAxes,ha="center",color="white",fontsize=7);a.plot([8,35],[10,10],color="white",lw=2);a.text(21.5,14,"30 arcsec",ha="center",color="white",fontsize=6.5)
    finish(fig,root/"figures/fig01_schmidt_field_and_asteroid_cutouts")

def make_fig06(root,imgs,cats,measurements):
    d=pd.read_parquet(measurements);d=d[(d.night.astype(str)=="20251121")&(d.trk_sub=="000004r")&(d.linkage_id==598)].sort_values("mjd").copy()
    centers=SkyCoord(d.ra_win_deg.to_numpy()*u.deg,d.dec_win_deg.to_numpy()*u.deg)
    tabs=[]
    for c in cats:
        with fits.open(c) as h: tabs.append(Table(h[1].data).to_pandas())
    t0=tabs[0];c0=SkyCoord(t0.RA_Win.to_numpy()*u.deg,t0.DEC_Win.to_numpy()*u.deg);sep=c0.separation(centers[0]).arcsec
    cand=t0[(sep>25)&(sep<120)&(t0.Mag_Aper4>14)&(t0.Mag_Aper4<19.5)&(t0.Flag==0)].copy()
    best=None
    for r in cand.itertuples(index=False):
        cc=SkyCoord(r.RA_Win*u.deg,r.DEC_Win*u.deg);ms=[]
        for tab in tabs:
            ct=SkyCoord(tab.RA_Win.to_numpy()*u.deg,tab.DEC_Win.to_numpy()*u.deg);idx,s,xx=cc.match_to_catalog_sky(ct);ms.append((tab.iloc[idx],s.arcsec))
        if max(x[1] for x in ms)<.7: best=(r,ms);break
    if best is None: raise RuntimeError("No stable reference star found")
    ref=SkyCoord(best[0].RA_Win*u.deg,best[0].DEC_Win*u.deg);mid=SkyCoord(np.mean(d.ra_win_deg)*u.deg,np.mean(d.dec_win_deg)*u.deg)
    scale=1.155/3600;n=150;wo=WCS(naxis=2);wo.wcs.crpix=[(n+1)/2,(n+1)/2];wo.wcs.crval=[mid.ra.deg,mid.dec.deg];wo.wcs.cdelt=[-scale,scale];wo.wcs.ctype=["RA---TAN","DEC--TAN"]
    cubes=[];headers=[]
    for p in imgs:
        im,wi,hh=load_image(p);cubes.append(reproj(im,wi,wo,(n,n)));headers.append(hh)
    ref_xy=np.array(wo.world_to_pixel_values([x[0].RA_Win for x in best[1]],[x[0].DEC_Win for x in best[1]])).T
    tx=np.array(wo.world_to_pixel_values(d.ra_win_deg,d.dec_win_deg)).T
    out=d[["night","trk_sub","linkage_id","image_name","mjd","ra_win_deg","dec_win_deg","x_win_px","y_win_px","mag_aper4","magerr_aper4"]].copy();out["registered_target_x_px"]=tx[:,0];out["registered_target_y_px"]=tx[:,1];out["reference_ra_deg"]=[x[0].RA_Win for x in best[1]];out["reference_dec_deg"]=[x[0].DEC_Win for x in best[1]];out["registered_reference_x_px"]=ref_xy[:,0];out["registered_reference_y_px"]=ref_xy[:,1];out.to_csv(root/"figure_data/fig06_astrometry.csv",index=False)
    scatter=np.hypot(ref_xy[:,0]-np.median(ref_xy[:,0]),ref_xy[:,1]-np.median(ref_xy[:,1]));write_json(root/"figure_data/fig06_registration_check.json",{"reference_star_selection":"fixed, unsaturated (Flag=0), isolated catalog source present in all three frames","reference_star_median_ra_deg":float(np.median(out.reference_ra_deg)),"reference_star_median_dec_deg":float(np.median(out.reference_dec_deg)),"max_registered_reference_offset_px":float(scatter.max()),"pixel_scale_arcsec":1.155,"max_registered_reference_offset_arcsec":float(scatter.max()*1.155),"pass_threshold_arcsec":0.8,"status":"PASS" if scatter.max()*1.155<.8 else "FAIL"})
    fig,axs=plt.subplots(1,3,figsize=(7.1,2.6),sharex=True,sharey=True)
    tmin=d.mjd.min()
    for k,(ax,a) in enumerate(zip(axs,cubes)):
        ax.imshow(zscale(a),origin="lower",cmap="gray",vmin=0,vmax=2.7);ax.add_patch(Circle(tx[k],7,fill=False,ec=ORANGE,lw=1.5));ax.add_patch(Rectangle(ref_xy[k]-6,12,12,fill=False,ec=BLUE,lw=1.5));ax.set(xticks=[],yticks=[]);ax.text(.04,.94,chr(65+k),transform=ax.transAxes,va="top",color="white",weight="bold");ax.text(.96,.94,f"+{(d.mjd.iloc[k]-tmin)*24*60:.1f} min",transform=ax.transAxes,ha="right",va="top",color="white",fontsize=8)
    axs[0].text(.04,.07,"moving source",transform=axs[0].transAxes,color=ORANGE,fontsize=7);axs[0].text(.04,.13,"fixed reference",transform=axs[0].transAxes,color=BLUE,fontsize=7);axs[2].annotate("N",xy=(.88,.82),xytext=(.88,.66),xycoords="axes fraction",textcoords="axes fraction",color="white",ha="center",arrowprops=dict(arrowstyle="-|>",color="white"));axs[2].annotate("E",xy=(.72,.80),xytext=(.88,.80),xycoords="axes fraction",textcoords="axes fraction",color="white",va="center",arrowprops=dict(arrowstyle="-|>",color="white"));axs[2].plot([.08,.31],[.08,.08],transform=axs[2].transAxes,color="white",lw=2);axs[2].text(.195,.11,"40 arcsec",transform=axs[2].transAxes,ha="center",color="white",fontsize=7)
    finish(fig,root/"figures/fig06_candidate_three_frame_sequence")

def main(a):
    root=Path(a.output).resolve();make_fig01(root,Path(a.fig01_image),Path(a.known_matched));make_fig06(root,[Path(x) for x in a.fig06_images],[Path(x) for x in a.fig06_catalogs],Path(a.measurements))
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True);p.add_argument("--fig01-image",required=True);p.add_argument("--known-matched",required=True);p.add_argument("--fig06-images",nargs=3,required=True);p.add_argument("--fig06-catalogs",nargs=3,required=True);p.add_argument("--measurements",required=True);main(p.parse_args())
