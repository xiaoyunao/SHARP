#!/usr/bin/env python3
"""Make the SHARP review atlas from the bundled, provenance-labelled inputs.

Requires numpy, pandas, matplotlib, astropy. No network or scheduler execution.
All sky positions are planning illustrations, not telescope dispatch commands.
"""
from pathlib import Path
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import astropy.units as u
from astropy.coordinates import SkyCoord, BarycentricMeanEcliptic
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
HEAD = "a5cb9454"
INK, MUTED = "#183049", "#536578"
BLUE, ORANGE, TEAL, RED, PURPLE = "#2475B5", "#D68434", "#168F86", "#C75359", "#7E61A7"
PALETTE = [BLUE, ORANGE, TEAL, PURPLE, RED, "#7B8D32", "#537F92", "#A05C85"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
    "axes.titlesize": 13, "axes.labelsize": 11, "axes.labelcolor": INK,
    "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "#FBFCFE",
    "savefig.facecolor": "white", "pdf.fonttype": 42,
    "axes.grid": True, "grid.alpha": .16, "grid.linewidth": .6})
warnings.filterwarnings("ignore", message=".*dubious year.*")

summary = pd.read_csv(DATA / "object_summary.csv", dtype={"object_id": str})
summary = summary.set_index("object_id", drop=False)
special = json.loads((DATA / "special_tracks.json").read_text())
names = {"617":"Patroclus", "17365":"Thymbraeus", "22180":"Paeon",
    "5511":"Cloanthus", "4902":"Thessandrus", "5012":"Eurymedon",
    "5284":"Orsilocus", "7543":"Prylis", "624":"Hektor",
    "11351":"Leucus", "3548":"Eurybates", "15094":"Polymele",
    "21900":"Orus", "13383":"1998 XS31", "38615":"2000 AV121",
    "228155":"2009 SF61", "2357":"Phereclos", "51364":"2000 SU333",
    "56976":"2000 SS161"}

def heading(fig, title, subtitle):
    fig.suptitle(title, x=.065, y=.978, ha="left", fontsize=20, fontweight="bold", color=INK)
    fig.text(.065, .936, subtitle, ha="left", va="top", fontsize=10.5, color=MUTED)

def footer(fig, text):
    fig.text(.065, .025, text, fontsize=8.5, color=MUTED, va="bottom", linespacing=1.4)

def save(fig, name, pdf):
    fig.savefig(OUT / (name + ".png"), dpi=165)
    pdf.savefig(fig)
    plt.close(fig)
    print(name, flush=True)

def sky(day):
    p = pd.read_csv(DATA / ("sky_" + day + ".csv.gz"), dtype={"pdes":str})
    p["q"] = p.pdes.map(summary.quality_probability)
    return p

def coord(ra, dec):
    return SkyCoord(np.asarray(ra)*u.deg, np.asarray(dec)*u.deg, frame="icrs")

def fp(ra, dec, width=3):
    h = width/2
    xx = [-h,0,h,h,h,0,-h,-h,-h]
    yy = [-h,-h,-h,0,h,h,h,0,-h]
    return SkyCoord(lon=np.asarray(xx)*u.deg, lat=np.asarray(yy)*u.deg,
                    frame=coord(ra,dec).skyoffset_frame()).icrs

def projection(ra, dec):
    return np.deg2rad((180-np.asarray(ra)+180)%360-180), np.deg2rad(dec)

def skyline(ax, ra, dec, **kwargs):
    xx, yy = projection(ra,dec)
    jumps = np.where(np.abs(np.diff(xx))>np.pi)[0]+1
    for xs,ys in zip(np.split(xx,jumps),np.split(yy,jumps)):
        ax.plot(xs,ys,**kwargs)

def setup_aitoff(ax):
    ticks = np.arange(-150,151,30)
    ax.set_xticks(np.deg2rad(ticks))
    ax.set_xticklabels([str((180-t)%360)+"°" for t in ticks], fontsize=9)
    ax.set_xlabel("Right ascension  (increases to the left)", labelpad=12)
    ax.set_ylabel("Declination", labelpad=20)
    ax.grid(True, alpha=.25)

def interp_track(key, jd):
    t = special[str(key)]
    ra = np.rad2deg(np.interp(jd,t["datetime_jd"],np.unwrap(np.deg2rad(t["RA"]))))%360
    dec = np.interp(jd,t["datetime_jd"],t["DEC"])
    v = np.interp(jd,t["datetime_jd"],t["V"])
    return ra,dec,v

def sky_overview(pdf):
    p = sky("2026-11-15")
    fig = plt.figure(figsize=(14,9))
    heading(fig,"01  |  Where the annual Trojan sample is on the sky",
        "Common-epoch RA/Dec Aitoff map  ·  15 November 2026  ·  Xinglong planning context")
    ax = fig.add_axes([.08,.24,.80,.64],projection="aitoff")
    setup_aitoff(ax)
    for group,color in [("L4",ORANGE),("L5",BLUE)]:
        g=p[p.group==group]
        xx,yy=projection(g.ra,g.dec)
        ax.scatter(xx,yy,s=2,c=color,alpha=.13,rasterized=True,linewidths=0)
    q=p[p.q.notna()]
    xx,yy=projection(q.ra,q.dec)
    sc=ax.scatter(xx,yy,s=24+30*q.q,c=q.q,cmap="viridis",vmin=0,vmax=1,
                  edgecolors="white",linewidths=.35,zorder=5)
    lon=np.linspace(0,360,1000)
    ec=SkyCoord(lon=lon*u.deg,lat=np.zeros_like(lon)*u.deg,
        frame=BarycentricMeanEcliptic(equinox=Time("J2000"))).icrs
    gc=SkyCoord(l=lon*u.deg,b=np.zeros_like(lon)*u.deg,frame="galactic").icrs
    skyline(ax,ec.ra.deg,ec.dec.deg,color=ORANGE,lw=1.15,ls="--",alpha=.8)
    skyline(ax,gc.ra.deg,gc.dec.deg,color="#89949E",lw=.8,ls=":",alpha=.9)
    skyline(ax,lon,np.full_like(lon,-19.607),color=RED,lw=1,ls="--",alpha=.75)
    for text,ra,dec,color in [("L4 cloud",219,43,ORANGE),("L5 cloud",69,53,BLUE)]:
        x,y=projection(ra,dec);ax.text(x,y,text,color=color,fontweight="bold",ha="center",fontsize=13)
    cax=fig.add_axes([.925,.35,.012,.38])
    cb=fig.colorbar(sc,cax=cax)
    cb.set_label("Generic lightcurve-quality probability",fontsize=9)
    fig.text(.08,.14,f"{len(q)} mapped IDs / {len(summary)} IDs sampled at least once across 200 simulated years",
        fontsize=11,fontweight="bold")
    fig.text(.08,.085,"Faint dots: reference orbit catalogue.  Dashed gold: ecliptic.  Dotted grey: Galactic plane.\n"
        "Red dashed line: approximate declination limit for any transit above 30° at Xinglong (an object-level bound).",fontsize=9.5,color=MUTED)
    footer(fig,"Positions: earlier two-body planning catalogue, all at one epoch; not the actual annual telescope pointings.\n"
        f"Probabilities: committed object_summary.csv at {HEAD}. No weather/precision outcome is inferred from a dot's location.")
    save(fig,"01_trojan_radec_aitoff",pdf)

def cloud_closeups(pdf):
    fig,axs=plt.subplots(1,2,figsize=(15,8.2))
    fig.subplots_adjust(left=.07,right=.915,bottom=.18,top=.83,wspace=.22)
    heading(fig,"02  |  Each cloud contains many 3° fields",
        "Seasonal cloud views  ·  background targets satisfy V ≤ 19.7 and >2 h above 30° during astronomical night")
    for ax,day,group,color in zip(axs,["2026-12-15","2027-04-15"],["L5","L4"],[BLUE,ORANGE]):
        p=sky(day);g=p[(p.group==group)&(p.V<=19.7)&(p.hours_alt30_sun18>=2)]
        ax.scatter(g.ra,g.dec,s=9,color="#B4C3D0",alpha=.7,label="Potential bright targets")
        h=g[g.q.notna()]
        sc=ax.scatter(h.ra,h.dec,c=h.q,cmap="viridis",vmin=0,vmax=1,s=25+45*h.q,
            edgecolor="white",lw=.5,zorder=4)
        xmin,xmax=g.ra.min()-3,g.ra.max()+3;ymin,ymax=g.dec.min()-3,g.dec.max()+3
        ax.set_xlim(xmax,xmin);ax.set_ylim(ymin,ymax)
        ax.set_title(f"{group}  ·  {day}\n{len(g)} geometry/brightness candidates; {len(h)} in annual sampled set",loc="left",color=color,fontsize=12)
        ax.set_xlabel("Right ascension (deg)");ax.set_ylabel("Declination (deg)")
        r0=xmin+5;d0=ymax-5
        poly=fp(r0,d0);ax.plot(poly.ra.deg,poly.dec.deg,color=INK,lw=1.5)
        ax.text(r0,d0-3,"3° × 3°",ha="center",va="top",fontsize=10,color=INK)
        keys=["617","17365","22180","5511"] if group=="L5" else ["4902","5012","5284","7543"]
        for j,k in enumerate(keys):
            rr=p[p.pdes==k]
            if rr.empty:continue
            r=rr.iloc[0]
            ax.scatter([r.ra],[r.dec],s=100,marker="*",facecolor=RED,edgecolor="white",lw=.6,zorder=7)
            dx,dy=((8,8) if j%2 else (-8,-15))
            ax.annotate(k,(r.ra,r.dec),xytext=(dx,dy),textcoords="offset points",fontsize=9,
                ha="left" if dx>0 else "right",color=RED,
                bbox=dict(fc="white",ec="none",alpha=.8,pad=1),zorder=8)
    cax=fig.add_axes([.935,.3,.012,.42]);cb=fig.colorbar(sc,cax=cax)
    cb.set_label("Generic quality probability",fontsize=9)
    fig.text(.07,.105,"Stars identify science cases discussed in the review. Highlighted annual IDs can be scheduled in a different month.",fontsize=10)
    footer(fig,"These are cloud/candidate maps, not selected annual field footprints. Positions use the earlier approximate orbit propagation.\n"
        "No Moon, field-edge, route or weather filter is applied to the background sample; a 3° field is only a small part of either cloud.")
    save(fig,"02_cloud_closeups",pdf)

def local_fields(pdf):
    jd=np.linspace(Time("2026-12-12T12:00:00").jd,Time("2026-12-19T12:00:00").jd,169)
    cases=[("17365",["17365","2357","51364"]),("22180",["22180","56976"])]
    fig,axs=plt.subplots(1,2,figsize=(15,8.5));fig.subplots_adjust(left=.07,right=.96,bottom=.18,top=.82,wspace=.27)
    heading(fig,"03  |  Science targets can anchor useful multi-object fields",
        "Proposed fixed 3° × 3° footprints for 12–19 December 2026  ·  fresh topocentric Horizons tracks")
    result=[];p=sky("2026-12-15")
    for ax,(anchor,keys) in zip(axs,cases):
        arrays=[interp_track(k,jd) for k in keys]
        allra=np.concatenate([a[0] for a in arrays]);alldec=np.concatenate([a[1] for a in arrays])
        center=coord((allra.min()+allra.max())/2,(alldec.min()+alldec.max())/2)
        # Recenter a spherical-offset bounding box, retaining a fixed PA=0 footprint.
        for _ in range(4):
            pts=coord(allra,alldec).transform_to(center.skyoffset_frame())
            lon=pts.lon.wrap_at(180*u.deg).deg;lat=pts.lat.deg
            c=SkyCoord(lon=(lon.min()+lon.max())/2*u.deg,lat=(lat.min()+lat.max())/2*u.deg,frame=center.skyoffset_frame())
            center=c.icrs
        frame=center.skyoffset_frame()
        nearby=p[(np.abs((p.ra-center.ra.deg)*np.cos(np.deg2rad(center.dec.deg)))<4.5)&(np.abs(p.dec-center.dec.deg)<4.5)&(p.V<=20.5)]
        pts=coord(nearby.ra.to_numpy(),nearby.dec.to_numpy()).transform_to(frame)
        ax.scatter(pts.lon.wrap_at(180*u.deg).deg,pts.lat.deg,s=np.clip((21-nearby.V)*12,5,55),
            c="#B5C3D0",alpha=.7,zorder=1)
        ax.add_patch(Polygon([[-1.5,-1.5],[1.5,-1.5],[1.5,1.5],[-1.5,1.5]],closed=True,
            facecolor=BLUE,alpha=.045,edgecolor="none"))
        ax.plot([-1.5,1.5,1.5,-1.5,-1.5],[-1.5,-1.5,1.5,1.5,-1.5],color=BLUE,lw=1.8,label="Proposed free centre")
        ax.plot([-1.35,1.35,1.35,-1.35,-1.35],[-1.35,-1.35,1.35,1.35,-1.35],color=BLUE,lw=.8,ls=":",label="0.15° edge margin")
        ar,ad,_=interp_track(anchor,np.array([np.mean(jd)]))
        ap=fp(ar[0],ad[0]).transform_to(frame)
        ax.plot(ap.lon.wrap_at(180*u.deg).deg,ap.lat.deg,color=RED,lw=1.1,ls="--",label="Centre on main target")
        retained={}
        for j,(k,(ra,dec,v)) in enumerate(zip(keys,arrays)):
            tr=coord(ra,dec).transform_to(frame);xx=tr.lon.wrap_at(180*u.deg).deg;yy=tr.lat.deg
            good=(np.abs(xx)<=1.35)&(np.abs(yy)<=1.35);retained[k]=float(np.mean(good))
            color=[TEAL,ORANGE,PURPLE][j%3]
            ax.plot(xx,yy,lw=2.3,color=color,zorder=6)
            ax.scatter([xx[0]],[yy[0]],marker="o",s=30,color=color,zorder=7)
            ax.scatter([xx[-1]],[yy[-1]],marker=">",s=45,color=color,zorder=7)
            ix=len(xx)//2
            ax.annotate(f"{k} {names.get(k,'')}\nV ≈ {np.median(v):.2f}", (xx[ix],yy[ix]),
                xytext=(0,13 if j%2==0 else -28),textcoords="offset points",ha="center",fontsize=9.5,
                color=color,bbox=dict(fc="white",ec="none",alpha=.85,pad=1.5),zorder=9)
        ax.scatter([0],[0],marker="+",c=BLUE,s=70,lw=1.4,zorder=8)
        ax.set_xlim(3.4,-3.4);ax.set_ylim(-3.2,3.2);ax.set_aspect("equal")
        ax.set_xlabel("Spherical tangent-frame longitude offset (deg)")
        ax.set_ylabel("Spherical tangent-frame latitude offset (deg)")
        ax.set_title(f"{anchor} {names[anchor]}\nProposed centre: RA {center.ra.deg:.3f}°, Dec {center.dec.deg:+.3f}°",loc="left",fontsize=12)
        ax.legend(loc="lower left",fontsize=8.5,framealpha=.95)
        result.append({"anchor":anchor,"companion_ids":keys,"center_ra_deg":center.ra.deg,"center_dec_deg":center.dec.deg,
            "track_start_utc":Time(jd[0],format="jd").isot,"track_stop_utc":Time(jd[-1],format="jd").isot,
            "continuous_grid_fraction_inside_0p15_margin":retained,
            "not_certified_for_execution":True,"scope":"Track-grid geometry only; no scheduled visits, Moon, stars, weather or detector mask"})
    fig.text(.07,.115,"Solid tracks: time increases from circle to triangle. Grey dots: nearby catalogue objects at the reference epoch.",fontsize=10)
    footer(fig,"A field centre need not sit on the primary target: the off-target centre can retain useful companions.\n"
        "Illustrative proposals only. Bright-star halos, actual visit times, Moon clearance, field retention and detector masks still require validation.")
    (DATA/"proposed_local_fields.json").write_text(json.dumps(result,indent=2))
    save(fig,"03_special_target_fields",pdf)

def twilight_maps(pdf):
    allsamples=json.loads((DATA/"twilight_samples.json").read_text())["samples"]
    samples=[s for s in allsamples if s["visits_per_field"]==5]
    fig=plt.figure(figsize=(15,11));gs=fig.add_gridspec(2,3,left=.065,right=.96,bottom=.15,top=.86,
        height_ratios=[1.25,1],hspace=.32,wspace=.23)
    heading(fig,"04  |  Twilight pointings: absolute sky and solar geometry",
        "Three representative nights reconstructed with the current scheduler  ·  five visits per field  ·  r / 30 s")
    ax=fig.add_subplot(gs[0,:],projection="aitoff");setup_aitoff(ax)
    for s,color in zip(samples,[BLUE,TEAL,ORANGE]):
        for b in s["blocks"]:
            for f in b["fields"]:
                rr=f["footprint_ra_deg"]+[f["footprint_ra_deg"][0]];dd=f["footprint_dec_deg"]+[f["footprint_dec_deg"][0]]
                skyline(ax,rr,dd,color=color,lw=.8,alpha=.9)
        ax.scatter([],[],marker="s",s=30,facecolor="none",edgecolor=color,label=s["date"])
    ax.legend(loc="upper center",bbox_to_anchor=(.5,1.12),ncol=3,frameon=False,fontsize=10)
    for i,s in enumerate(samples):
        ax=fig.add_subplot(gs[1,i]);ax.axvspan(-45,45,color="#DDE3EB",alpha=.45)
        ax.axhline(0,color=MUTED,lw=.8);ax.axvline(0,color=ORANGE,lw=.7)
        for b in s["blocks"]:
            color=ORANGE if b["block_id"].startswith("EV") else BLUE
            for f in b["fields"]:
                rr=f["footprint_sun_relative_lon_deg"];dd=f["footprint_ecliptic_lat_deg"]
                ax.add_patch(Polygon(np.column_stack([rr,dd]),closed=True,facecolor=color,edgecolor=color,alpha=.22,lw=.8))
                ax.scatter([f["reference_sun_relative_lon_deg"]],[f["reference_ecliptic_lat_deg"]],color=color,s=8)
        ax.set_xlim(-95,95);ax.set_ylim(-28,28);ax.set_xticks([-90,-60,0,60,90])
        ax.set_xlabel("Solar-relative ecliptic longitude (deg)")
        if i==0:ax.set_ylabel("Ecliptic latitude (deg)")
        ev=sum(len(b["fields"]) for b in s["blocks"] if b["block_id"].startswith("EV"))
        am=sum(len(b["fields"]) for b in s["blocks"] if b["block_id"].startswith("AM"))
        ax.set_title(f"{s['date']}\nEV: {ev} fields  |  AM: {am} fields",loc="left",fontsize=11.5)
        ax.text(-67,-24,"Morning",ha="center",color=BLUE,fontsize=9)
        ax.text(67,-24,"Evening",ha="center",color=ORANGE,fontsize=9)
    footer(fig,"These are new-group reconstructions (empty previous-anchor history), not exported original annual pointings; all six 5/10-visit plans passed the repository validator.\n"
        "Footprints are spherical 3° squares. Solar-relative longitude is not the same quantity as total solar elongation.")
    save(fig,"04_twilight_radec_and_solar_maps",pdf)

def twilight_altitudes(pdf):
    allsamples=json.loads((DATA/"twilight_samples.json").read_text())["samples"]
    samples=[s for s in allsamples if s["visits_per_field"]==5 and s["date"] in ["2026-11-04","2027-09-23"]]
    fig,axs=plt.subplots(2,2,figsize=(15,10.5));fig.subplots_adjust(left=.07,right=.96,bottom=.15,top=.85,hspace=.38,wspace=.22)
    heading(fig,"05  |  Twilight night cards: heights at the actual planned visits",
        "Reconstructed representative nights  ·  dots mark exposure midpoints  ·  pale bands bound the full field conservatively")
    for row,s in enumerate(samples):
        for col,side in enumerate(["EV","AM"]):
            ax=axs[row,col];g=s["environment_grids"][side];xx=np.array(g["elapsed_minutes"])
            blocks=[b for b in s["blocks"] if b["block_id"].startswith(side)]
            ax.axhspan(-40,30,color="#DDE3EB",alpha=.42)
            ax.axhline(30,color=RED,lw=1.2,ls="--");ax.axhline(-12,color=ORANGE,lw=.7,ls=":")
            for j,f in enumerate(g["field_tracks"]):
                color=PALETTE[j%len(PALETTE)];alt=np.array(f["alt_deg"])
                ax.plot(xx,alt,color=color,lw=1,alpha=.55)
                ax.fill_between(xx,alt-2.122,alt,color=color,alpha=.08)
                rows=[r for b in blocks for r in b["rows"] if r["field_id"]==f["field_id"]]
                tx=(Time([r["mid_utc"] for r in rows]).unix-Time(g["utc"][0]).unix)/60
                ax.scatter(tx,[r["alt_deg"] for r in rows],s=20,color=color,zorder=5)
            ax.plot(xx,g["sun_alt_deg"],color=ORANGE,lw=1.8,label="Sun")
            ax.plot(xx,g["moon_alt_deg"],color="#747F8A",lw=1.1,ls="--",label="Moon")
            ax.text(59,31,"30° field limit",ha="right",va="bottom",color=RED,fontsize=8.5)
            ax.text(59,-11,"Sun −12°",ha="right",va="bottom",color=ORANGE,fontsize=8.5)
            ax.set_ylim(-37,85);ax.set_xlim(0,60)
            ticks=[0,15,30,45,60];ax.set_xticks(ticks)
            ax.set_xticklabels([g["beijing_time"][i][11:16] for i in ticks])
            ax.set_xlabel("Beijing time (UTC+8)");ax.set_ylabel("Altitude (deg)")
            localdate=g["beijing_time"][0][:10]
            n=sum(len(b["rows"]) for b in blocks)
            ax.set_title(f"Night {s['date']}  |  {side} on {localdate}\n{len(g['field_tracks'])} fields · {n} exposures",loc="left",fontsize=11.5)
            if row==0 and col==0:ax.legend(loc="upper right",fontsize=9,frameon=True)
    footer(fig,"Curves are field-centre altitudes; the lower edge of each pale band is centre altitude minus 2.122°. Only exposure-time validity is decisive.\n"
        "The AM date is the following civil date. These examples have no weather losses and do not inherit the original annual anchor chain.")
    save(fig,"05_twilight_altitude_cards",pdf)

def trojan_night_cards(pdf):
    p=DATA/"special_night_tracks.json"
    if not p.exists():return
    data=json.loads(p.read_text())
    # This adapter intentionally fails on an unknown schema rather than inventing coordinates.
    nights=data["nights"]
    if isinstance(nights,dict):nights=list(nights.values())
    fig,axs=plt.subplots(2,1,figsize=(15,10.5));fig.subplots_adjust(left=.07,right=.96,bottom=.15,top=.85,hspace=.38)
    heading(fig,"06  |  Trojan night cards: targets, daylight and Moon clearance",
        "Fresh Horizons ephemerides at Xinglong  ·  geometric opportunities, not a selected or executed two-hour schedule")
    groups=[["617","17365","22180","5511"],["4902","5012","5284","7543"]]
    for ax,n,keys in zip(axs,nights,groups):
        times=n["beijing"];hour=np.arange(len(times))*5/60+16
        sun=np.asarray(n["sun_alt_deg"]);moon=np.asarray(n["moon_alt_deg"]);frac=np.asarray(n["moon_illuminated_fraction"])
        ax.fill_between(hour,0,90,where=sun>-18,color=ORANGE,alpha=.11,step="mid")
        ax.axhspan(0,30,color="#DDE3EB",alpha=.42);ax.axhline(30,color=RED,lw=1,ls="--")
        ax.plot(hour,moon,color="#87919D",ls=":",lw=1.2,label="Moon altitude")
        for j,k in enumerate(keys):
            t=n["targets"][k];a=np.asarray(t["alt_deg"]);sep=np.asarray(t["moon_sep_deg"])
            color=PALETTE[j]
            good=(sun<-18)&(a>32.8)&((moon<0)|(sep>np.where(frac>=.75,42.122,32.122)))
            ax.plot(hour,a,color=color,alpha=.35,lw=1)
            ax.plot(hour,np.where(good,a,np.nan),color=color,lw=2.5,label=f"{k} {names[k]}  (V {np.median(t['V']):.1f})")
        ax.set_xlim(16,32);ax.set_ylim(0,90);ax.set_xticks(np.arange(16,33,2));ax.set_xticklabels([f"{int(h)%24:02d}:00" for h in np.arange(16,33,2)])
        ax.set_xlabel("Beijing time  ·  right of midnight is the next civil day");ax.set_ylabel("Altitude (deg)")
        ax.set_title(n["date"]+"  |  "+("L5 cases" if keys[0]=="617" else "L4 cases"),loc="left")
        ax.legend(loc="upper right",fontsize=9,ncol=2,framealpha=.95)
    footer(fig,"Bold segments pass a review-level centre >32.8° test and 30°/40° Moon clearance plus a 2.122° field radius, with Sun <−18°.\n"
        "The operational scheduler has additional time guards, bright-star exclusions, reservations and actual route costs; these curves do not certify dispatch.")
    save(fig,"06_trojan_nightly_altitudes",pdf)

def season_targets(pdf):
    p=DATA/"special_visibility_monthly.csv"
    if not p.exists():return
    d=pd.read_csv(p,dtype={"object_id":str})
    keys=["617","17365","22180","5511","4902","5012","5284","7543","21900","15094","3548","624","11351","13383","38615","228155"]
    dates=sorted(d.date.unique());mat=np.zeros((len(keys),len(dates)));vm=mat.copy();two=np.zeros_like(mat,dtype=bool)
    for i,k in enumerate(keys):
        for j,day in enumerate(dates):
            r=d[(d.object_id==k)&(d.date==day)].iloc[0]
            mat[i,j]=r.dark_alt30_h;vm[i,j]=r.V_midnight
            two[i,j]=str(r.has_continuous_2h_window).lower()=="true"
    fig,ax=plt.subplots(figsize=(15,11));fig.subplots_adjust(left=.22,right=.925,bottom=.16,top=.85)
    heading(fig,"07  |  Which special targets are practical in this observing year?",
        "Twelve representative nights, each on the 15th  ·  colour = hours at target altitude >30° with Sun <−18°")
    im=ax.imshow(mat,origin="upper",aspect="auto",cmap="YlGnBu",vmin=0,vmax=10)
    for i in range(len(keys)):
        for j in range(len(dates)):
            if mat[i,j]>=1:
                ax.text(j,i,f"{vm[i,j]:.1f}"+("*" if two[i,j] else ""),ha="center",va="center",
                    fontsize=9,color="white" if mat[i,j]>5.2 else INK)
            else:ax.text(j,i,"—",ha="center",va="center",fontsize=9,color="#8D98A3")
    ax.set_yticks(np.arange(len(keys)));ax.set_yticklabels([k+"  "+names[k] for k in keys],fontsize=10)
    ax.set_xticks(np.arange(len(dates)));ax.set_xticklabels([x[:7] for x in dates],rotation=35,ha="right")
    ax.set_xticks(np.arange(-.5,len(dates),1),minor=True);ax.set_yticks(np.arange(-.5,len(keys),1),minor=True)
    ax.grid(which="minor",color="white",lw=.8,alpha=.75);ax.grid(which="major",visible=False);ax.tick_params(which="minor",length=0)
    ax.axhline(7.5,color=ORANGE,lw=2);ax.axhline(10.5,color=RED,lw=2)
    cax=fig.add_axes([.94,.28,.012,.44]);cb=fig.colorbar(im,cax=cax);cb.set_label("Geometrically observable hours",fontsize=9)
    footer(fig,"Cell numbers are predicted V, not calibrated r. Asterisk: a continuous 2 h review-level centre/field-radius/Moon window exists on that representative night.\n"
        "Rows 1–8: main science candidates; next 3: conditional Lucy targets; last 5: low-altitude or faint cases. A blank month is not a proof of zero availability on every date.")
    save(fig,"07_special_target_seasons",pdf)

def science_return(pdf):
    fig,axs=plt.subplots(2,2,figsize=(15,11));fig.subplots_adjust(left=.10,right=.96,bottom=.15,top=.85,wspace=.4,hspace=.44)
    heading(fig,"08  |  What the current queue rewards — and what the budget buys",
        f"Committed annual results at {HEAD}  ·  model outcomes under stated weather and photometric assumptions")
    ax=axs[0,0];roles=["control","known_other","unknown_or_weak","alias_review","literature_unchecked"]
    totals=summary.groupby("science_role").quality_probability.sum();labels=["Known control","Other known","Unknown / weak","Alias review","Unchecked"]
    values=[float(totals.get(k,0)) for k in roles]
    ax.barh(labels,values,color=[BLUE,BLUE,TEAL,ORANGE,"#ACB9C5"],height=.62)
    for i,v in enumerate(values):ax.text(v+.35,i,f"{v:.2f}",va="center",fontsize=10)
    ax.invert_yaxis();ax.set_xlim(0,32);ax.set_xlabel("Expected qualifying objects  =  sum of individual probabilities")
    ax.set_title("A. Current science-role mix",loc="left")
    ax.text(.98,.06,f"Known/control share: {(values[0]+values[1])/sum(values):.0%}\nExpected total: {sum(values):.2f}",transform=ax.transAxes,ha="right",fontsize=10,bbox=dict(fc="white",ec="none",alpha=.9))
    ax=axs[0,1];ids=["17365","22180","5511","617","4902","5012","5284","7543"]
    yy=np.arange(len(ids));probs=[float(summary.loc[k,"quality_probability"]) if k in summary.index else np.nan for k in ids]
    ax.barh(yy,np.nan_to_num(probs),height=.58,color=[RED,TEAL,BLUE,BLUE,ORANGE,BLUE,"#BBC5D0","#BBC5D0"])
    ax.set_yticks(yy);ax.set_yticklabels([k+" "+names[k] for k in ids],fontsize=9);ax.invert_yaxis();ax.set_xlim(0,1.19)
    for i,v in enumerate(probs):ax.text((v+.025) if np.isfinite(v) else .025,i,f"{v:.1%}" if np.isfinite(v) else "not sampled",va="center",fontsize=9)
    ax.set_title("B. Interesting objects can still be missed",loc="left");ax.set_xlabel("Generic quality probability; not science-case completion")
    ax=axs[1,0];b=pd.read_csv(DATA/"budget_precision.csv");xx=np.arange(len(b))
    ax.axhline(101,color=RED,lw=1.1,ls="--",label="More than 100 objects")
    for i,r in b.iterrows():
        color=BLUE if r.hours_per_block==2 else ORANGE
        ax.errorbar(i,r.P50,yerr=[[r.P50-r.P10],[r.P90-r.P50]],fmt="o",color=color,capsize=4,ms=6)
        ax.text(i,r.P90+5,f"{r.P50:g}",ha="center",fontsize=10,color=color)
    ax.set_xticks(xx);ax.set_xticklabels([f"{int(r.hours_per_block)} h\n{r.sigma_r18:.3f}" for _,r in b.iterrows()],fontsize=9)
    ax.set_ylim(0,150);ax.set_ylabel("Qualifying objects: median and P10–P90");ax.set_xlabel("Night allocation / assumed random error at r = 18 (mag)")
    ax.set_title("C. Precision and time are separate assumptions",loc="left");ax.legend(fontsize=8.5,loc="upper left")
    ax=axs[1,1];c=pd.read_csv(DATA/"cadence.csv")
    for i,r in c.iterrows():
        color=BLUE if int(r.dwell)==1 else ORANGE
        ax.plot([r.gap_s_P10/60,r.gap_s_P90/60],[i,i],color=color,lw=3,alpha=.8)
        ax.scatter([r.gap_s_P50/60],[i],color=color,s=30,zorder=3)
    ax.set_yticks(np.arange(len(c)));ax.set_yticklabels([f"{int(r.fields)} fields × {int(r.dwell)} exp/dwell" for _,r in c.iterrows()],fontsize=9)
    ax.invert_yaxis();ax.set_xlim(0,35);ax.set_xlabel("Revisit gap (min): median dot and P10–P90 interval")
    ax.set_title("D. Paired exposures hide the long return gap",loc="left")
    footer(fig,"The generic quality test is a photometry/coverage criterion; it does not by itself establish a new unique rotation period or complete a special science case.\n"
        "A 4 h allocation and 0.010 mag precision are sensitivity experiments, not demonstrated telescope performance or the adopted project time request.")
    save(fig,"08_science_yield_and_cadence",pdf)

def annual_calendar(pdf):
    d=pd.read_csv(DATA/"monthly_visibility.csv")
    fig,axs=plt.subplots(2,1,figsize=(15,9.5),gridspec_kw={"height_ratios":[1.25,1]})
    fig.subplots_adjust(left=.075,right=.96,bottom=.16,top=.84,hspace=.32)
    heading(fig,"09  |  The annual programme has two seasons and a real gap",
        "Monthly results from the submitted full-year simulation  ·  1 November 2026 to 31 October 2027")
    xx=np.arange(len(d));ax=axs[0]
    ax.bar(xx,d.L5_proxy,color=BLUE,label="L5 proxy",width=.72)
    ax.bar(xx,d.L4_proxy,bottom=d.L5_proxy,color=ORANGE,label="L4 proxy",width=.72)
    for i,r in d.iterrows():ax.text(i,r.eligible_objects+6,str(int(r.eligible_objects)),ha="center",fontsize=10)
    ax.set_ylim(0,270);ax.set_xticks(xx);ax.set_xticklabels(d.month);ax.set_ylabel("Monthly unique eligible objects")
    ax.set_title("Candidates can change strongly with season; do not force fields into every month",loc="left")
    ax.legend(ncol=2,loc="upper right",frameon=False)
    ax=axs[1]
    ax.bar(xx,d.eligible_nights,color=TEAL,width=.72,label="Nights with eligible fields")
    ax.set_ylim(0,35);ax.set_ylabel("Eligible nights in month");ax.set_xticks(xx);ax.set_xticklabels(d.month)
    for i,r in d.iterrows():ax.text(i,r.eligible_nights+.7,str(int(r.eligible_nights)),ha="center",fontsize=10)
    ax.set_title("August has no eligible Trojan fields in the current full two-hour menu",loc="left")
    footer(fig,"Source: monthly_visibility.csv. Eligibility includes the implemented full-window geometry/menu rules, not all possible shorter or redesigned visits.\n"
        "Field versions are deliberately not shown as independent sky area: a new seven-day version can cover essentially the same patch and the same objects.")
    save(fig,"09_annual_season_calendar",pdf)

def main():
    path=ROOT/"SHARP_observing_review_atlas_20261104.pdf"
    with PdfPages(path,metadata={"Title":"SHARP annual observing review — sky maps and science targets","Author":"SHARP planning review","Subject":"Read-only review of a5cb9454; explicitly labelled planning illustrations"}) as pdf:
        sky_overview(pdf);cloud_closeups(pdf);local_fields(pdf);twilight_maps(pdf);twilight_altitudes(pdf)
        trojan_night_cards(pdf);season_targets(pdf);science_return(pdf);annual_calendar(pdf)
    print(str(path),flush=True)

if __name__=="__main__":main()
