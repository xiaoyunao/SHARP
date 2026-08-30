#!/usr/bin/env python3
"""Recompute representative SHARP distance-grid solutions with production code."""
from pathlib import Path
import argparse, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
from common import *

def predicted(module,best,mjd,loc):
    t=Time(mjd,format="mjd",scale="utc");rso=module.heliocentric_observer_xyz_au(t,loc);r=[]
    for tt in t:
        o=best["orbit"].propagate(tt-best["orbit"].epoch);r.append(o.r.to_value("AU"))
    rho=np.asarray(r)-rso;rho/=np.linalg.norm(rho,axis=1)[:,None]
    return np.degrees(np.arctan2(rho[:,1],rho[:,0]))%360,np.degrees(np.arcsin(rho[:,2]))

def case_data(d,night,trk,lid):return d[(d.night.astype(str)==str(night))&(d.trk_sub==trk)&(d.linkage_id==lid)].sort_values("mjd").drop_duplicates("detection_index")

def main(a):
    root=Path(a.output).resolve();sys.path.insert(0,str(Path(a.repo).resolve()/"heliolincrr"));import orbit_confirm_links as oc
    d=pd.read_parquet(a.measurements);cases=[("Ordinary retained linkage",20251121,"000004r",598),("C/2025 Y1 (ATLAS)",20260310,"000012H",1026)];dist=np.arange(1.3,4.1001,.2);rows=[];tracks=[]
    loc=oc.earth_location_mpc327()
    for label,night,trk,lid in cases:
        z=case_data(d,night,trk,lid);mjd=z.mjd.to_numpy();ra=z.ra_win_deg.to_numpy();dec=z.dec_win_deg.to_numpy()
        for rr in dist:
            b=oc.fit_best_hypo_lambert(mjd,ra,dec,loc,[(float(rr),0.,0.)],.02,35.)
            rows.append({"case":label,"night":night,"trk_sub":trk,"linkage_id":lid,"hypothesis_distance_au":rr,"fit_ok":b.get("ok",False),"rms_arcsec":b.get("rms_arcsec",np.nan),"median_arcsec":b.get("med_arcsec",np.nan),"max_arcsec":b.get("max_arcsec",np.nan),"site_altitude_m":oc.MPC327_ALT_M,"source_code":"heliolincrr/orbit_confirm_links.py"})
        b=oc.fit_best_hypo_lambert(mjd,ra,dec,loc,[tuple(x) for x in oc.PROFILE_DEFAULTS["single-night"]["hypos"]],.02,35.)
        pra,pdec=predicted(oc,b,mjd,loc)
        ra0=np.mean(ra);dec0=np.mean(dec);ox=(ra-ra0)*3600*np.cos(np.deg2rad(dec0));oy=(dec-dec0)*3600;px=(pra-ra0)*3600*np.cos(np.deg2rad(dec0));py=(pdec-dec0)*3600
        for k in range(len(z)):tracks.append({"case":label,"night":night,"trk_sub":trk,"linkage_id":lid,"mjd":mjd[k],"observed_x_arcsec":ox[k],"observed_y_arcsec":oy[k],"model_x_arcsec":px[k],"model_y_arcsec":py[k],"residual_dx_arcsec":ox[k]-px[k],"residual_dy_arcsec":oy[k]-py[k],"best_distance_au":b["hypo_r_au"],"rms_arcsec":b["rms_arcsec"]})
    fam=pd.DataFrame(rows);tr=pd.DataFrame(tracks);fam.to_csv(root/"figure_data/fig05c_solution_families.csv",index=False);tr.to_csv(root/"figure_data/fig05b_representative_tracks.csv",index=False)
    fig,axs=plt.subplots(1,2,figsize=(7.1,3.15))
    for ax,(label,z) in zip(axs,tr.groupby("case",sort=False)):
        ax.plot(z.model_x_arcsec,z.model_y_arcsec,color=BLUE,lw=1.2,label="distance-grid orbit model");ax.scatter(z.observed_x_arcsec,z.observed_y_arcsec,s=28,facecolor="white",edgecolor=ORANGE,lw=1.3,label="observed")
        ax.quiver(z.model_x_arcsec,z.model_y_arcsec,z.residual_dx_arcsec,z.residual_dy_arcsec,angles="xy",scale_units="xy",scale=1,color=ORANGE,width=.006)
        ax.set(title=label,xlabel=r"$\Delta\alpha\cos\delta$ (arcsec)",ylabel=r"$\Delta\delta$ (arcsec)");ax.grid(alpha=.16);ax.set_aspect("equal",adjustable="datalim")
    axs[0].legend(frameon=False);finish(fig,root/"figures/fig05b_representative_tracks")
    fig,axs=plt.subplots(1,2,figsize=(7.1,3.15),sharey=True)
    for ax,(label,z) in zip(axs,fam.groupby("case",sort=False)):
        ax.plot(z.hypothesis_distance_au,z.rms_arcsec,color=BLUE,marker="o",ms=3);ax.fill_between(z.hypothesis_distance_au,0,z.rms_arcsec,where=z.rms_arcsec<=1,color=ORANGE,alpha=.15);ax.axhline(1,color=ORANGE,ls="--",lw=1);ax.set(title=label,xlabel="Hypothesized heliocentric distance (AU)",ylabel="Orbit-fit RMS (arcsec)",ylim=(0,max(1.08,np.nanpercentile(fam.rms_arcsec,98)))) ;ax.grid(alpha=.16)
    finish(fig,root/"figures/fig05c_solution_families")
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True);p.add_argument("--measurements",required=True);p.add_argument("--repo",required=True);main(p.parse_args())
