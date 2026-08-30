#!/usr/bin/env python3
"""Redraw SHARP science figures by adapting the bundled GOTTA source snapshots.

The ecliptic map, orbital-class plot, histogram, density-scatter, colorbar,
palette, marker, font, and save conventions are copied from the exact GOTTA
sources under ``scripts/gotta_reference_source`` and only the input-column
mapping is changed for SHARP.
"""
from __future__ import annotations
import argparse
from pathlib import Path

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, get_sun
from astropy.time import Time
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator

GROUP_ORDER = ["MBA", "OMB", "TJN", "IMB", "MCA", "NEA", "TNO", "Other"]
GROUP_STYLE = {
    "MBA": ("o", "#17becf", 0.16, 11),
    "OMB": ("s", "#f58518", 0.72, 34),
    "TJN": ("D", "#54a24b", 0.72, 34),
    "IMB": ("^", "#b279a2", 0.72, 34),
    "MCA": ("v", "#e45756", 0.78, 38),
    "NEA": ("*", "#ff9da6", 0.90, 80),
    "TNO": ("P", "#d62728", 0.88, 90),
    "Other": ("X", "#7f7f7f", 0.82, 42),
}

def setup_style():
    plt.rcParams.update({
        "font.family": "Times New Roman", "font.size": 30, "axes.linewidth": 1.4,
        "xtick.major.width": 1.2, "ytick.major.width": 1.2,
        "xtick.major.size": 7, "ytick.major.size": 7,
        "legend.frameon": False, "figure.dpi": 150,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

def save_figure(fig, root: Path, name: str, pad: float = 0.08):
    out=root/"figures"/name
    fig.savefig(out.with_suffix(".png"),dpi=300,bbox_inches="tight",pad_inches=pad)
    fig.savefig(out.with_suffix(".pdf"),bbox_inches="tight",pad_inches=pad)
    plt.close(fig)

def add_colorbar(ax,mappable,label):
    divider=make_axes_locatable(ax);cax=divider.append_axes("right",size="3.5%",pad=0.06)
    cb=ax.figure.colorbar(mappable,cax=cax);cb.set_label(label,fontsize=20);cb.ax.tick_params(labelsize=17,width=.9,length=4);cb.outline.set_linewidth(.8)

def histogram(ax,values,bins,xlabel,ylabel="Detections",logy=False,color="#4c78a8"):
    values=np.asarray(values,dtype=float);values=values[np.isfinite(values)]
    ax.hist(values,bins=bins,color=color,histtype="stepfilled",alpha=.5,edgecolor="#2b2b2b",linewidth=.55)
    ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
    if logy:ax.set_yscale("log")
    ax.set_axisbelow(True);ax.grid(alpha=.16,linewidth=.6);ax.tick_params(labelsize=22)

def density_colored_scatter(ax,x,y,xlabel,ylabel,xbins=260,ybins=190,xlim=None,ylim=None,xscale=None,point_size=3.0):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float);mask=np.isfinite(x)&np.isfinite(y)
    if xscale=="log":mask&=x>0
    if xlim is not None:mask&=(x>=xlim[0])&(x<=xlim[1])
    if ylim is not None:mask&=(y>=ylim[0])&(y<=ylim[1])
    xp,yp=x[mask],y[mask];xw=np.log10(xp) if xscale=="log" else xp
    xr=np.log10(xlim) if xscale=="log" and xlim is not None else (np.nanmin(xw),np.nanmax(xw));yr=ylim if ylim is not None else (np.nanmin(yp),np.nanmax(yp))
    counts,xe,ye=np.histogram2d(xw,yp,bins=[xbins,ybins],range=[xr,yr]);sm=gaussian_filter(counts.astype(float),sigma=2.2,mode="nearest")
    xc=.5*(xe[:-1]+xe[1:]);yc=.5*(ye[:-1]+ye[1:]);interp=RegularGridInterpolator((xc,yc),sm,bounds_error=False,fill_value=1.0)
    den=np.clip(interp(np.column_stack([xw,yp])),1,None);order=np.argsort(den)
    sc=ax.scatter(xp[order],yp[order],c=den[order],s=point_size,cmap="viridis",norm=LogNorm(vmin=1,vmax=max(1,float(np.nanmax(den)))),marker="o",linewidths=0,alpha=.82,rasterized=True)
    ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
    if xscale:ax.set_xscale(xscale)
    if xlim is not None:ax.set_xlim(*xlim)
    if ylim is not None:ax.set_ylim(*ylim)
    ax.set_axisbelow(True);ax.grid(alpha=.12,linewidth=.55);ax.tick_params(labelsize=22);add_colorbar(ax,sc,"Local density")

def running_stats(x,y,edges,min_count=40):
    d=pd.DataFrame({"x":x,"y":y}).replace([np.inf,-np.inf],np.nan).dropna();rows=[]
    for l,r in zip(edges[:-1],edges[1:]):
        z=d[(d.x>=l)&(d.x<r)].y
        if len(z)>=min_count:rows.append((np.sqrt(l*r) if l>0 and r/l>1.2 else .5*(l+r),np.percentile(z,16),np.median(z),np.percentile(z,84),len(z)))
    return pd.DataFrame(rows,columns=["x","p16","median","p84","count"])

def lon_to_mollweide_rad(lon_deg,center_deg=180.):
    centered=(np.asarray(lon_deg)-center_deg+180)%360-180;return np.deg2rad(centered)

def set_mollweide_ticks(ax):
    labels=np.arange(60,301,60);ticks=lon_to_mollweide_rad(labels);ax.set_xticks(ticks);ax.set_xticklabels([])
    for v,x in zip(labels,ticks):ax.text(x,np.deg2rad(-15),f"{v:d}°",ha="center",va="center",fontsize=24,color="black",zorder=10)

def plot_celestial_equator(ax):
    ra=np.linspace(0,360,2000);e=SkyCoord(ra=ra*u.deg,dec=np.zeros_like(ra)*u.deg,frame="icrs").barycentrictrueecliptic
    x=lon_to_mollweide_rad(e.lon.deg);order=np.argsort(x);ax.plot(x[order],np.deg2rad(e.lat.deg)[order],color="#ff8c00",linewidth=3,alpha=.95,zorder=9)

def fig03(root):
    sky=pd.read_csv(root/"figure_data/fig03_known_sky_healpix.csv");nside=int(sky.healpix_nside.iloc[0]);npix=hp.nside2npix(nside);counts=np.zeros(npix);counts[sky.pixel.astype(int)]=sky.detection_count
    theta,phi=hp.pix2ang(nside,np.arange(npix));lat=np.rad2deg(.5*np.pi-theta);lon=np.rad2deg(phi);x=lon_to_mollweide_rad(lon);y=np.deg2rad(lat)
    cmap=plt.get_cmap("rainbow").copy();cmap.set_bad(alpha=0);val=np.ma.masked_where(counts==0,np.clip(counts,1,100))
    fig=plt.figure(figsize=(14,8));ax=fig.add_subplot(111,projection="mollweide");im=ax.scatter(x,y,c=val,s=8,linewidths=0,cmap=cmap,norm=LogNorm(vmin=1,vmax=100),rasterized=True)
    ax.grid(True,alpha=.3);set_mollweide_ticks(ax);ax.tick_params(axis="y",labelsize=24);plot_celestial_equator(ax)
    cax=ax.inset_axes([.30,.09,.40,.045]);cb=fig.colorbar(im,cax=cax,orientation="horizontal");cb.set_label(r"$\mathrm{Counts\ per\ pixel}$",fontsize=18);cb.ax.xaxis.set_label_position("top");cb.ax.tick_params(labelsize=16,direction="in",length=5,width=1.1)
    fig.tight_layout();save_figure(fig,root,"fig03_known_sky_distribution",0)

def orbit_group(values):
    m={"Main belt":"MBA","Outer belt / Hilda":"OMB","Jupiter Trojan region":"TJN","Inner belt / Hungaria":"IMB","Mars-crosser":"MCA","NEO (q < 1.3 AU)":"NEA","Distant / Centaur-TNO":"TNO"}
    return np.array([m.get(x,"Other") for x in values])

def fig04(root):
    d=pd.read_csv(root/"figure_data/fig04_unique_object_orbits.csv");a=d.a_au.to_numpy();e=d.ecc.to_numpy();inc=d.inc_deg.to_numpy();groups=orbit_group(d.orbit_class);valid=np.isfinite(a)&np.isfinite(e)&np.isfinite(inc)&(a>0)&(inc>0)
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(18,9))
    for g in GROUP_ORDER:
        mask=valid&(groups==g)
        if mask.any():
            marker,color,alpha,size=GROUP_STYLE[g];ax1.scatter(a[mask],e[mask],marker=marker,color=color,s=size,alpha=alpha,linewidths=0,rasterized=True);ax2.scatter(a[mask],inc[mask],marker=marker,color=color,s=max(size,28),alpha=alpha,linewidths=0,rasterized=True)
    al=np.logspace(np.log10(1.3),2,100);ax1.plot(al,1-1.3/al,"k--",linewidth=2.3,label=r"$q=1.3\,\mathrm{AU}$");ax1.legend(fontsize=22,loc="lower right",frameon=False,handlelength=2.4);ax1.set_xscale("log");ax1.set_xlabel("Semimajor axis [AU]");ax1.set_ylabel("Eccentricity")
    ax2.set_xscale("log");ax2.set_yscale("log");ax2.set_xlabel("Semimajor axis [AU]");ax2.set_ylabel("Inclination [°]")
    handles=[Line2D([0],[0],marker=GROUP_STYLE[g][0],linestyle="None",markerfacecolor=GROUP_STYLE[g][1],markeredgecolor=GROUP_STYLE[g][1],alpha=max(GROUP_STYLE[g][2],.75),markersize=12 if g!="TNO" else 14,label=g) for g in GROUP_ORDER]
    ax2.legend(handles=handles,fontsize=22,loc="best",frameon=False,handletextpad=.4,borderpad=.2,labelspacing=.35,markerscale=1.15)
    fig.tight_layout();save_figure(fig,root,"fig04_known_orbital_distribution")

def fig04_series(root,snapshot):
    matched=pd.read_parquet(snapshot/"frozen_products/known_matched.parquet");rates=pd.read_parquet(root/"figure_data/fig04b_known_motion_rates.parquet").angular_speed_arcsec_hr
    targets=SkyCoord(matched.pred_ra_deg.to_numpy()*u.deg,matched.pred_dec_deg.to_numpy()*u.deg);sun=get_sun(Time(matched.epoch_mjd.to_numpy(),format="mjd",scale="utc")).transform_to("icrs");elong=targets.separation(sun).deg
    dist=pd.read_parquet(root/"figure_data/fig04c_known_distances.parquet")
    dra=((matched.obs_ra_win_deg-matched.pred_ra_deg+180)%360-180)*np.cos(np.deg2rad(matched.pred_dec_deg))*3600;ddec=(matched.obs_dec_win_deg-matched.pred_dec_deg)*3600;sep=np.hypot(dra,ddec)
    fig,ax=plt.subplots(1,2,figsize=(18,7));histogram(ax[0],rates,np.logspace(-1,3.2,65),r"Angular rate [arcsec h$^{-1}$]",logy=True);ax[0].set_xscale("log");med=np.nanmedian(rates);ax[0].axvline(med,color="black",ls="--",lw=1.6,label=rf"median = {med:.2f} arcsec h$^{{-1}}$");ax[0].legend(fontsize=20,loc="upper left");histogram(ax[1],elong,np.linspace(0,180,61),"Solar elongation [deg]",color="#59a14f");med=np.nanmedian(elong);ax[1].axvline(med,color="black",ls="--",lw=1.6,label=f"median = {med:.2f} deg");ax[1].legend(fontsize=20,loc="upper left");fig.tight_layout();save_figure(fig,root,"fig04b_known_motion_geometry")
    fig,ax=plt.subplots(1,2,figsize=(18,7));histogram(ax[0],dist.heliocentric_distance_au,np.linspace(0,6,61),"Heliocentric distance [AU]",color="#f28e2b");med=np.nanmedian(dist.heliocentric_distance_au);ax[0].axvline(med,color="black",ls="--",lw=1.6,label=f"median = {med:.2f} AU");ax[0].legend(fontsize=20,loc="upper right");histogram(ax[1],dist.geocentric_distance_au,np.linspace(0,5,61),"Geocentric distance [AU]",color="#e15759");med=np.nanmedian(dist.geocentric_distance_au);ax[1].axvline(med,color="black",ls="--",lw=1.6,label=f"median = {med:.2f} AU");ax[1].legend(fontsize=20,loc="upper right");fig.tight_layout();save_figure(fig,root,"fig04c_known_distances")
    mag=matched.mag_aper4;err=matched.magerr_aper4
    fig,ax=plt.subplots(1,2,figsize=(18,7));histogram(ax[0],mag,np.linspace(10,21,45),r"$g_{\rm aper}$ [mag]",color="#4c78a8");med=np.nanmedian(mag);ax[0].axvline(med,color="black",ls="--",lw=1.6,label=f"median = {med:.2f}");ax[0].legend(fontsize=20,loc="upper left");density_colored_scatter(ax[1],mag,err,r"$g_{\rm aper}$ [mag]",r"$\sigma(g_{\rm aper})$ [mag]",xlim=(10,21),ylim=(0,.6));fig.subplots_adjust(left=.08,right=.94,bottom=.16,top=.96,wspace=.34);save_figure(fig,root,"fig04d_known_photometry")
    fig,ax=plt.subplots(1,2,figsize=(19,7));density_colored_scatter(ax[0],mag,sep,r"$g_{\rm aper}$ [mag]","Separation [arcsec]",xlim=(10,21),ylim=(0,1.5));st=running_stats(mag,sep,np.linspace(10,21,23));ax[0].plot(st.x,st["median"],color="#d62728",lw=2.2,label="binned median");ax[0].fill_between(st.x,st.p16,st.p84,color="#d62728",alpha=.42,label="16--84 percentile");ax[0].legend(fontsize=18,loc="upper left")
    # Associate each finite-difference rate with the corresponding residual row by object/night ordering.
    keyed=matched.copy();n=pd.to_numeric(keyed.asteroid_number,errors="coerce");keyed["object_key"]=np.where(n.notna(),"N:"+n.fillna(-1).astype(int).astype(str),"S:"+keyed.asteroid_name.astype(str).str.lower().str.replace(r"[^a-z0-9]","",regex=True));keyed["sep"]=sep;keyed=keyed.sort_values(["object_key","night","epoch_mjd"]);keyed["rate"]=np.nan
    rv=pd.read_parquet(root/"figure_data/fig04b_known_motion_rates.parquet");rate_map={(r.object_key,str(r.night),float(r.epoch_mjd)):r.angular_speed_arcsec_hr for r in rv.itertuples()};keyed["rate"]=[rate_map.get((k,str(night),float(t)),np.nan) for k,night,t in zip(keyed.object_key,keyed.night,keyed.epoch_mjd)]
    density_colored_scatter(ax[1],keyed.rate,keyed.sep,r"Angular rate [arcsec h$^{-1}$]","Separation [arcsec]",xlim=(.7,150),ylim=(0,1.5),xscale="log");st=running_stats(keyed.rate,keyed.sep,np.logspace(np.log10(.7),np.log10(150),21),20);ax[1].plot(st.x,st["median"],color="#d62728",lw=2.2,label="binned median");ax[1].fill_between(st.x,st.p16,st.p84,color="#d62728",alpha=.42,label="16--84 percentile");ax[1].legend(fontsize=18,loc="upper left");fig.subplots_adjust(left=.07,right=.95,bottom=.17,top=.96,wspace=.34);save_figure(fig,root,"fig04e_known_astrometric_trends")
    rev=pd.read_csv(root/"figure_data/fig04f_known_revisit_statistics.csv");vmax=max(1,rev.detection_count.max());sizes=np.interp(np.log10(rev.detection_count),[0,np.log10(vmax)],[8,115]);fig,ax=plt.subplots(1,2,figsize=(18,7));sc=ax[0].scatter(rev.night_count,rev.time_baseline_days,c=rev.detection_count,s=sizes,cmap="viridis",norm=LogNorm(vmin=1,vmax=vmax),alpha=.86,linewidths=0,rasterized=True);ax[0].set_xlabel("Distinct nights per object");ax[0].set_ylabel("First-to-last baseline [days]");ax[0].set_axisbelow(True);ax[0].grid(alpha=.12,lw=.55);ax[0].tick_params(labelsize=22);add_colorbar(ax[0],sc,"Detections per object");histogram(ax[1],rev.detection_count,np.logspace(0,np.log10(vmax+1),50),"Detections per object",ylabel="Unique objects",logy=True,color="#f28e2b");ax[1].set_xscale("log");fig.subplots_adjust(left=.07,right=.95,bottom=.16,top=.96,wspace=.36);save_figure(fig,root,"fig04f_known_revisit_statistics")

def fig05(root):
    d=pd.read_csv(root/"tables/candidate_linkage_diagnostics.csv");fig,ax=plt.subplots(figsize=(10,7.4));norm=Normalize(vmin=d.n_measurements.min(),vmax=d.n_measurements.max());base=~d.is_c2025_y1.astype(bool);sc=ax.scatter(d.loc[base,"arc_length_min"],d.loc[base,"rms_arcsec"],c=d.loc[base,"n_measurements"],cmap="viridis",norm=norm,s=62,marker="o",alpha=.86,linewidths=0,rasterized=True,label="Other retained linkages");special=~base;ax.scatter(d.loc[special,"arc_length_min"],d.loc[special,"rms_arcsec"],c=d.loc[special,"n_measurements"],cmap="viridis",norm=norm,s=125,marker="*",edgecolor="#e45756",linewidth=1.3,label="C/2025 Y1 confirmations");ax.axhline(1,color="black",ls="--",lw=1.8,label="Nominal 1 arcsec threshold");ax.set_xlabel("Observed arc length [min]");ax.set_ylabel("Orbit-fit RMS [arcsec]");ax.set_ylim(0,max(1.08,float(d.rms_arcsec.max())*1.08));ax.grid(alpha=.12,lw=.55);ax.tick_params(labelsize=22);ax.legend(fontsize=17,loc="upper right");add_colorbar(ax,sc,"Measurements per linkage");fig.subplots_adjust(left=.16,right=.88,bottom=.18,top=.96);save_figure(fig,root,"fig05_candidate_orbit_diagnostics")
    tr=pd.read_csv(root/"figure_data/fig05b_representative_tracks.csv");fig,axes=plt.subplots(1,2,figsize=(18,8),constrained_layout=True)
    for ax,(label,z) in zip(axes,tr.groupby("case",sort=False)):
        ax.plot(z.model_x_arcsec,z.model_y_arcsec,color="#4c78a8",lw=2.2,label="short-arc orbit model");ax.scatter(z.observed_x_arcsec,z.observed_y_arcsec,s=85,facecolor="white",edgecolor="#f28e2b",lw=2,label="observations",zorder=4);ax.quiver(z.model_x_arcsec,z.model_y_arcsec,z.residual_dx_arcsec,z.residual_dy_arcsec,angles="xy",scale_units="xy",scale=1,color="#d62728",width=.006,zorder=3);ax.set_title(label,fontsize=25);ax.set_xlabel(r"$\Delta\alpha\cos\delta$ [arcsec]");ax.set_ylabel(r"$\Delta\delta$ [arcsec]");ax.grid(alpha=.12,lw=.55);ax.tick_params(labelsize=22);ax.set_aspect("equal",adjustable="datalim")
    axes[0].legend(fontsize=18,loc="best");save_figure(fig,root,"fig05b_representative_tracks",.14)

def main(a):
    root=Path(a.root).resolve();snapshot=Path(a.snapshot).resolve();setup_style();fig03(root);fig04(root);fig04_series(root,snapshot);fig05(root)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--root",required=True);p.add_argument("--snapshot",required=True);main(p.parse_args())
