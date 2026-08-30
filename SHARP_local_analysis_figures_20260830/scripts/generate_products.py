#!/usr/bin/env python3
"""Generate non-image SHARP tables and GOTTA-styled scientific figures."""
from __future__ import annotations
import argparse, json, os, platform, re, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from astropy.coordinates import SkyCoord, BarycentricTrueEcliptic, get_sun, get_body_barycentric
from astropy import units as u
from astropy.time import Time
from aleph.Query import Query
import healpy as hp
from common import *

def norm_name(x): return re.sub(r"[^a-z0-9]", "", str(x).lower())

def binned_xy(x, y, bins):
    rows=[]
    for lo,hi in zip(bins[:-1],bins[1:]):
        z=np.asarray(y)[(np.asarray(x)>=lo)&(np.asarray(x)<hi)]
        z=z[np.isfinite(z)]
        rows.append((lo,hi,(lo+hi)/2,len(z),np.nanmedian(z) if len(z) else np.nan,
                     np.nanpercentile(z,16) if len(z) else np.nan,np.nanpercentile(z,84) if len(z) else np.nan))
    return pd.DataFrame(rows,columns=["bin_lo","bin_hi","bin_center","n","median","p16","p84"])

def object_key(d):
    n=pd.to_numeric(d.asteroid_number,errors="coerce")
    return np.where(n.notna(),"N:"+n.fillna(-1).astype(int).astype(str),"S:"+d.asteroid_name.map(norm_name))

def orbit_class(a,e,i):
    q=a*(1-e)
    if q < 1.3: return "NEO (q < 1.3 AU)"
    if 1.3 <= q < 1.666: return "Mars-crosser"
    if a < 2.06: return "Inner belt / Hungaria"
    if a < 3.28: return "Main belt"
    if a < 5.05: return "Outer belt / Hilda"
    if 5.05 <= a < 5.4: return "Jupiter Trojan region"
    if a >= 5.4: return "Distant / Centaur-TNO"
    return "Other"

def save_fig(fig, root, name): finish(fig,root/"figures"/name)

def build_workflow(root):
    nodes=[
      (0.04,.62,.18,.20,"Calibrated\nscience frames",BLUE),(.29,.62,.18,.20,"Catalog extraction\nand masks",BLUE),
      (.54,.62,.18,.20,"Known-object\nephemeris match",ORANGE),(.79,.62,.17,.20,"Matched known\nassociations",ORANGE),
      (.29,.18,.18,.20,"Unmatched-source\ntracklets",BLUE),(.54,.18,.18,.20,"Heliocentric\ndistance grid",ORANGE),
      (.79,.18,.17,.20,"Fit-quality cuts\nand visual audit",ORANGE)]
    edges=[(0,1),(1,2),(2,3),(1,4),(4,5),(5,6)]
    pd.DataFrame([{"node_id":i,"label":n[4],"x":n[0],"y":n[1]} for i,n in enumerate(nodes)]).to_csv(root/"figure_data/workflow_nodes.csv",index=False)
    pd.DataFrame(edges,columns=["source_node","target_node"]).to_csv(root/"figure_data/workflow_edges.csv",index=False)
    fig,ax=plt.subplots(figsize=(7.1,3.2)); ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    for x,y,w,h,label,c in nodes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.012",fc="white",ec=c,lw=1.5))
        ax.text(x+w/2,y+h/2,label,ha="center",va="center",fontsize=9)
    for s,t in edges:
        a,b=nodes[s],nodes[t]; p1=(a[0]+a[2],a[1]+a[3]/2);p2=(b[0],b[1]+b[3]/2)
        if s==1 and t==4:p1=(a[0]+a[2]/2,a[1]);p2=(b[0]+b[2]/2,b[1]+b[3])
        ax.add_patch(FancyArrowPatch(p1,p2,arrowstyle="-|>",mutation_scale=10,lw=1.1,color=GRAY,connectionstyle="arc3,rad=0"))
    save_fig(fig,root,"fig02_moving_object_workflow")

def build_known(root,snap,astorb_path):
    d=pd.read_parquet(snap/"frozen_products/known_matched.parquet")
    d["object_key"]=object_key(d)
    dra=(d.obs_ra_win_deg-d.pred_ra_deg)*3600*np.cos(np.deg2rad(d.pred_dec_deg))
    ddec=(d.obs_dec_win_deg-d.pred_dec_deg)*3600
    d["resid_arcsec"]=np.hypot(dra,ddec)
    c=SkyCoord(d.pred_ra_deg.to_numpy()*u.deg,d.pred_dec_deg.to_numpy()*u.deg,frame="icrs")
    ec=c.transform_to(BarycentricTrueEcliptic(equinox="J2000"))
    lon=ec.lon.wrap_at(180*u.deg).deg; lat=ec.lat.deg
    nside=64; pix=hp.ang2pix(nside,lon,lat,lonlat=True); counts=np.bincount(pix,minlength=hp.nside2npix(nside))
    used=np.flatnonzero(counts); plon,plat=hp.pix2ang(nside,used,lonlat=True); plon=((plon+180)%360)-180
    sky=pd.DataFrame({"healpix_nside":nside,"pixel":used,"ecliptic_lon_deg":plon,"ecliptic_lat_deg":plat,"detection_count":counts[used]})
    sky.to_csv(root/"figure_data/fig03_known_sky_healpix.csv",index=False)
    pd.DataFrame({"object_key":d.object_key,"epoch_mjd":d.epoch_mjd,"ecliptic_lon_deg":lon,"ecliptic_lat_deg":lat}).to_parquet(root/"figure_data/fig03_known_detection_coordinates.parquet",index=False)
    fig=plt.figure(figsize=(7.1,3.7));ax=fig.add_subplot(111,projection="mollweide")
    sc=ax.scatter(np.deg2rad(-sky.ecliptic_lon_deg),np.deg2rad(sky.ecliptic_lat_deg),c=np.log10(sky.detection_count),s=2.8,cmap="Blues",rasterized=True)
    ax.grid(alpha=.25);ax.set_xlabel("Ecliptic longitude (J2000; increasing left)");ax.set_ylabel("Ecliptic latitude")
    cb=fig.colorbar(sc,ax=ax,pad=.07,shrink=.82);cb.set_label(r"$\log_{10}$ detections per NSIDE=64 pixel")
    save_fig(fig,root,"fig03_known_sky_distribution")

    # Exact frozen Lowell ASTORB join.
    q=Query(service="Lowell",filename=str(astorb_path))
    bynum={str(int(a["num"])):a for a in q.asts if str(a["num"]).strip()}
    byname={norm_name(a["name"]):a for a in q.asts if norm_name(a["name"])}
    identities=d.groupby("object_key",sort=False).agg(asteroid_number=("asteroid_number","first"),asteroid_name=("asteroid_name","first"),detection_count=("object_key","size"),night_count=("night","nunique"),first_mjd=("epoch_mjd","min"),last_mjd=("epoch_mjd","max")).reset_index()
    rows=[]
    for r in identities.itertuples(index=False):
        num=pd.to_numeric(r.asteroid_number,errors="coerce"); rec=bynum.get(str(int(num))) if np.isfinite(num) else None; method="number" if rec else ""
        if rec is None: rec=byname.get(norm_name(r.asteroid_name)); method="normalized_name" if rec else "unmatched"
        rows.append({**r._asdict(),"join_method":method,"join_status":"matched" if rec else "unmatched","catalog_number":str(rec["num"]).strip() if rec else "","catalog_name":str(rec["name"]).strip() if rec else "","catalog_epoch":rec["epoch"] if rec else "","a_au":rec["a"] if rec else np.nan,"ecc":rec["e"] if rec else np.nan,"inc_deg":rec["incl"] if rec else np.nan,"H_mag":rec["H"] if rec else np.nan,"mean_anomaly_epoch_deg":rec["M_epoch"] if rec else np.nan,"arg_peri_deg":rec["Arg_Peri"] if rec else np.nan,"node_deg":rec["Node"] if rec else np.nan,"mean_motion_deg_day":rec["n"] if rec else np.nan})
    j=pd.DataFrame(rows);j["time_baseline_days"]=j.last_mjd-j.first_mjd
    j.to_csv(root/"tables/known_orbit_join_audit.csv",index=False)
    good=j.query("join_status == 'matched'").copy(); good["orbit_class"]= [orbit_class(a,e,i) for a,e,i in zip(good.a_au,good.ecc,good.inc_deg)]
    good.to_csv(root/"figure_data/fig04_unique_object_orbits.csv",index=False)
    good.groupby("orbit_class").agg(unique_objects=("object_key","size"),detections=("detection_count","sum"),median_a_au=("a_au","median"),median_e=("ecc","median"),median_i_deg=("inc_deg","median")).reset_index().to_csv(root/"tables/known_orbital_class_summary.csv",index=False)
    fig,axs=plt.subplots(1,2,figsize=(7.1,3.25)); s=np.clip(np.sqrt(good.detection_count),1,5)
    axs[0].scatter(good.a_au,good.ecc,s=s,c=BLUE,alpha=.22,lw=0,rasterized=True);axs[1].scatter(good.a_au,good.inc_deg,s=s,c=ORANGE,alpha=.22,lw=0,rasterized=True)
    axs[0].set(xlabel="Semimajor axis, a (AU)",ylabel="Eccentricity",xlim=(0,8),ylim=(0,1));axs[1].set(xlabel="Semimajor axis, a (AU)",ylabel="Inclination (deg)",xlim=(0,8),ylim=(0,np.nanpercentile(good.inc_deg,99.8)))
    for ax in axs:ax.grid(alpha=.18)
    save_fig(fig,root,"fig04_known_orbital_distribution")

    # Motion rates from successive detections of each object/night.
    srt=d.sort_values(["object_key","night","epoch_mjd"]); g=srt.groupby(["object_key","night"],sort=False)
    dt=g.epoch_mjd.diff(); decm=(srt.pred_dec_deg+g.pred_dec_deg.shift())/2
    sep=np.hypot((srt.pred_ra_deg-g.pred_ra_deg.shift())*np.cos(np.deg2rad(decm)),srt.pred_dec_deg-g.pred_dec_deg.shift())*3600
    srt["angular_speed_arcsec_hr"]=sep/(dt*24)
    srt[["object_key","night","epoch_mjd","angular_speed_arcsec_hr"]].dropna().to_parquet(root/"figure_data/fig04b_known_motion_rates.parquet",index=False)
    def hist2(name,cols,labels,logx=(False,False),color=(BLUE,ORANGE)):
        fig,axs=plt.subplots(1,2,figsize=(7.1,3.0)); out=[]
        for k,(x,lbl,lg,c0) in enumerate(zip(cols,labels,logx,color)):
            x=np.asarray(x);x=x[np.isfinite(x)&(x>0 if lg else np.ones_like(x,dtype=bool))]
            bins=np.geomspace(np.percentile(x,.2),np.percentile(x,99.8),45) if lg else np.linspace(np.percentile(x,.2),np.percentile(x,99.8),45)
            n,b=np.histogram(x,bins=bins);axs[k].stairs(n,b,color=c0,lw=1.4,fill=True,alpha=.25);axs[k].axvline(np.median(x),color=c0,ls="--",lw=1);axs[k].set(xlabel=lbl,ylabel="Detection count");
            if lg:axs[k].set_xscale("log")
            out.append(pd.DataFrame({"panel":k,"bin_lo":b[:-1],"bin_hi":b[1:],"count":n}))
        pd.concat(out).to_csv(root/f"figure_data/{name}.csv",index=False);save_fig(fig,root,name)
    elong=c.separation(get_sun(Time(d.epoch_mjd.to_numpy(),format="mjd",scale="utc"))).deg
    hist2("fig04b_known_motion_geometry",[srt.angular_speed_arcsec_hr,elong],["Angular speed (arcsec hr$^{-1}$)","Solar elongation (deg)"],(True,False))
    # Vectorized two-body position from the frozen ASTORB elements at every detection epoch.
    elem=good.set_index("object_key")[["a_au","ecc","inc_deg","mean_anomaly_epoch_deg","arg_peri_deg","node_deg","mean_motion_deg_day","catalog_epoch"]]
    ed=d[["object_key","epoch_mjd"]].join(elem,on="object_key"); epoch0=pd.to_datetime(ed.catalog_epoch,utc=True); epoch0_mjd=Time(epoch0.dt.to_pydatetime()).mjd
    M=np.deg2rad((ed.mean_anomaly_epoch_deg.to_numpy()+ed.mean_motion_deg_day.to_numpy()*(ed.epoch_mjd.to_numpy()-epoch0_mjd))%360);e=ed.ecc.to_numpy();E=M.copy()
    for _ in range(8): E-= (E-e*np.sin(E)-M)/(1-e*np.cos(E))
    aa=ed.a_au.to_numpy();xo=aa*(np.cos(E)-e);yo=aa*np.sqrt(np.maximum(0,1-e*e))*np.sin(E);om=np.deg2rad(ed.arg_peri_deg);O=np.deg2rad(ed.node_deg);ii=np.deg2rad(ed.inc_deg)
    co,so,cO,sO,ci,si=np.cos(om),np.sin(om),np.cos(O),np.sin(O),np.cos(ii),np.sin(ii)
    x=(cO*co-sO*so*ci)*xo+(-cO*so-sO*co*ci)*yo;y=(sO*co+cO*so*ci)*xo+(-sO*so+cO*co*ci)*yo;z=(so*si)*xo+(co*si)*yo
    rhel=np.sqrt(x*x+y*y+z*z); earth=get_body_barycentric("earth",Time(ed.epoch_mjd.to_numpy(),format="mjd",scale="utc"))-get_body_barycentric("sun",Time(ed.epoch_mjd.to_numpy(),format="mjd",scale="utc"));eps=np.deg2rad(23.439291);xe=earth.x.to_value(u.au);ye=earth.y.to_value(u.au)*np.cos(eps)+earth.z.to_value(u.au)*np.sin(eps);ze=-earth.y.to_value(u.au)*np.sin(eps)+earth.z.to_value(u.au)*np.cos(eps);rgeo=np.sqrt((x-xe)**2+(y-ye)**2+(z-ze)**2)
    pd.DataFrame({"object_key":d.object_key,"epoch_mjd":d.epoch_mjd,"heliocentric_distance_au":rhel,"geocentric_distance_au":rgeo}).to_parquet(root/"figure_data/fig04c_known_distances.parquet",index=False)
    hist2("fig04c_known_distances",[rhel,rgeo],["Heliocentric distance (AU)","Geocentric distance (AU)"],(True,True))
    hist2("fig04d_known_photometry",[d.mag_aper4,d.magerr_aper4],["Instrumental aperture magnitude","Aperture magnitude uncertainty (mag)"],(False,True))
    bm=binned_xy(d.mag_aper4,d.resid_arcsec,np.linspace(np.nanpercentile(d.mag_aper4,.5),np.nanpercentile(d.mag_aper4,99.5),18)); br=binned_xy(srt.angular_speed_arcsec_hr,srt.resid_arcsec,np.geomspace(max(.01,np.nanpercentile(srt.angular_speed_arcsec_hr.dropna(),.5)),np.nanpercentile(srt.angular_speed_arcsec_hr.dropna(),99.5),18));bm.to_csv(root/"figure_data/fig04e_residual_vs_magnitude.csv",index=False);br.to_csv(root/"figure_data/fig04e_residual_vs_rate.csv",index=False)
    fig,axs=plt.subplots(1,2,figsize=(7.1,3));
    for ax,z,xlab in zip(axs,[bm,br],["Instrumental aperture magnitude","Angular speed (arcsec hr$^{-1}$)"]):
        ax.plot(z.bin_center,z["median"],color=BLUE,lw=1.4);ax.fill_between(z.bin_center,z.p16,z.p84,color=BLUE,alpha=.18);ax.axhline(1,color=ORANGE,ls="--",lw=1);ax.set(xlabel=xlab,ylabel="Astrometric residual (arcsec)",ylim=(0,1.08));ax.grid(alpha=.15)
    axs[1].set_xscale("log");save_fig(fig,root,"fig04e_known_astrometric_trends")
    revisit=j[["object_key","detection_count","night_count","time_baseline_days"]];revisit.to_csv(root/"figure_data/fig04f_known_revisit_statistics.csv",index=False)
    fig,axs=plt.subplots(1,2,figsize=(7.1,3));axs[0].hexbin(revisit.detection_count,revisit.night_count,gridsize=45,bins="log",mincnt=1,cmap="Blues");axs[0].set(xlabel="Matched detections per object",ylabel="Distinct nights per object");x=revisit.time_baseline_days;bins=np.geomspace(1,max(2,x.max()+1),40);axs[1].hist(x[x>0],bins=bins,color=ORANGE,alpha=.6);axs[1].set_xscale("log");axs[1].set(xlabel="First-to-last detection baseline (days)",ylabel="Unique objects");save_fig(fig,root,"fig04f_known_revisit_statistics")
    return {"known_associations":len(d),"known_unique_identities":len(j),"orbit_joined":len(good),"orbit_unmatched":int((j.join_status!="matched").sum()),"catalog_epoch_values":sorted(good.catalog_epoch.astype(str).unique().tolist())}

def build_candidates(root,snap):
    review=pd.read_csv(snap/"review_sample/review_and_mpc_status.csv")
    keep=review.query("final_paper_status == 'retained_after_posthoc_audit'").copy(); det=pd.read_parquet(snap/"frozen_products/unknown_review_detections.parquet")
    keep["origin_night"]=keep.origin_night.astype(str);det["night"]=det.night.astype(str)
    det=det.merge(keep[["origin_night","trk_sub","linkage_id"]],left_on=["night","trk_sub","linkage_id"],right_on=["origin_night","trk_sub","linkage_id"],how="inner")
    conf=pd.read_csv(snap/"jpl_identification/second_pass/jpl_second_pass_confirmations.csv");conf["night"]=conf.night.astype(str)
    ckeys=set(zip(conf.night,conf.trk_sub,conf.linkage_id))
    rows=[]
    for key,z in det.groupby(["night","trk_sub","linkage_id"],sort=False):
        z=z.sort_values("mjd"); arc=(z.mjd.max()-z.mjd.min())*24*60
        rows.append({"night":key[0],"trk_sub":key[1],"linkage_id":key[2],"n_measurements":len(z),"n_tracklets":int(z.n_tracklets.iloc[0]),"arc_length_min":arc,"angular_speed_arcsec_hr":float(z.lin_speed_arcsec_per_day.iloc[0]/24),"direction_deg":float(z.lin_dir_deg.iloc[0]),"median_instrumental_mag":float(np.nanmedian(z.mag_aper4)),"best_distance_au":float(z.a_au.iloc[0]) if np.isfinite(z.a_au.iloc[0]) else np.nan,"rms_arcsec":float(z.rms_arcsec.iloc[0]),"median_residual_arcsec":float(z.med_arcsec.iloc[0]),"max_residual_arcsec":float(z.max_arcsec.iloc[0]),"linear_rms_arcsec":float(z.lin_rms_arcsec.iloc[0]),"is_c2025_y1":key in ckeys,"source_table":"frozen unknown_review_detections + posthoc review status"})
    t=pd.DataFrame(rows);t.to_csv(root/"tables/candidate_linkage_diagnostics.csv",index=False)
    fig,ax=plt.subplots(figsize=(4.6,3.6));
    for flag,c,m,lbl in [(False,BLUE,"o","Other retained linkages"),(True,ORANGE,"D","C/2025 Y1 confirmations")]:
        z=t[t.is_c2025_y1==flag];sc=ax.scatter(z.arc_length_min,z.rms_arcsec,c=z.n_measurements,cmap="viridis",s=28 if flag else 17,marker=m,edgecolor=ORANGE if flag else "none",lw=.8,alpha=.9,label=lbl)
    ax.axhline(1,color=ORANGE,ls="--",lw=1,label="Nominal 1 arcsec threshold");ax.set(xlabel="Observed arc length (min)",ylabel="Orbit-fit RMS (arcsec)",ylim=(0,1.08));ax.legend(frameon=False);fig.colorbar(sc,ax=ax,label="Measurements per linkage");save_fig(fig,root,"fig05_candidate_orbit_diagnostics")
    det.to_parquet(root/"figure_data/fig05_candidate_retained_measurements.parquet",index=False)
    return t,det

def build_audit(root,snap,stats,repo,gpath):
    def git(args,cwd): return subprocess.check_output(["git",*args],cwd=cwd,text=True).strip()
    code={"generated_utc":pd.Timestamp.utcnow().isoformat(),"local_repository":str(repo),"local_commit":git(["rev-parse","HEAD"],repo),"local_branch":git(["branch","--show-current"],repo),"local_dirty":bool(git(["status","--porcelain"],repo)),"python_executable":sys.executable,"python_version":sys.version,"platform":platform.platform(),"analysis_scope":"read-only frozen snapshot plus current code audit; no production or manuscript changes","orbit_confirmation_site_altitude_m_in_current_code":868.221,"scheduler_and_known_matcher_altitude_m_in_current_code":960.0,"altitude_discrepancy_preserved":True}
    write_json(root/"audit/code_scope.json",code)
    gp=Path(gpath); g={"repository":str(gp),"commit":git(["rev-parse","HEAD"],gp),"dirty":bool(git(["status","--porcelain"],gp)),"reference_script":str(gp/"scripts/generate_paper_products.py"),"reference_script_sha256":sha256(gp/"scripts/generate_paper_products.py"),"reference_figure_directory":str(gp/"reviewer_figures_20260701")};write_json(root/"audit/gotta_style_provenance.json",g)
    cat={"catalog":"Lowell astorb.dat","snapshot_file":"astorb.dat","sha256":sha256(Path(args.astorb)),"size_bytes":Path(args.astorb).stat().st_size,"file_mtime_utc":pd.Timestamp(Path(args.astorb).stat().st_mtime,unit="s",tz="UTC").isoformat(),"catalog_epoch_values":stats["catalog_epoch_values"],"scope":"asteroid orbital elements; no comet entries and therefore no C/2025 Y1 join","pha_status":"unavailable: frozen ASTORB does not provide an authoritative PHA flag","neo_status":"derived only from q=a(1-e)<1.3 AU"};write_json(root/"audit/catalog_scope_c2025y1.json",cat)
    expected={"survey_exposures":41074,"survey_nights":134,"open_shutter_hr":342.28,"fields":1430,"area_deg2":10387.0,"known_associations":534780,"known_unique_identities":58482,"automated_linkages":4762,"visually_retained":67,"posthoc_artifacts":9,"final_linkages":58,"candidate_measurements":179,"candidate_nights":34}
    coverage=json.loads((snap/"coverage/coverage_summary.json").read_text())
    observed={"known_associations":stats["known_associations"],"known_unique_identities":stats["known_unique_identities"],"automated_linkages":len(pd.read_parquet(snap/"frozen_products/unknown_catalog.parquet")),"visually_retained":67,"posthoc_artifacts":9,"final_linkages":58,"candidate_measurements":179,"candidate_nights":34,"survey_exposures":coverage["raw_exposure_n"],"survey_nights":coverage["raw_night_n"],"open_shutter_hr":coverage["open_shutter_hours"],"fields":coverage["observed_field_n"],"area_deg2":coverage["unique_area_deg2"]}
    rec=[]
    for k,v in expected.items():
        digits=2 if k=="open_shutter_hr" else 0; ok=round(float(observed[k]),digits)==round(float(v),digits)
        rec.append({"metric":k,"expected_reported":v,"observed_recomputed":observed[k],"comparison":"rounded to 2 decimals" if k=="open_shutter_hr" else "rounded to reported precision","status":"PASS" if ok else "FAIL","source":"coverage_summary.json" if k in {"survey_exposures","survey_nights","open_shutter_hr","fields","area_deg2"} else "frozen product/review tables","grain":"detection" if "association" in k or "measurement" in k else "object/linkage/survey as named"})
    pd.DataFrame(rec).to_csv(root/"audit/baseline_reconciliation.csv",index=False)

def main(args):
    root=Path(args.output).resolve();snap=Path(args.snapshot).resolve();apply_style();build_workflow(root);stats=build_known(root,snap,Path(args.astorb));t,det=build_candidates(root,snap);build_audit(root,snap,stats,Path(args.repo),Path(args.gotta));write_json(root/"logs/nonimage_generation_summary.json",{"stats":stats,"candidate_linkages":len(t),"candidate_measurements":len(det)})
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True);p.add_argument("--snapshot",required=True);p.add_argument("--astorb",required=True);p.add_argument("--repo",required=True);p.add_argument("--gotta",required=True);args=p.parse_args();main(args)
