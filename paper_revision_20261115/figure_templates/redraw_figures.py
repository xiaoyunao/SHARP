#!/usr/bin/env python3
"""Redraw SHARP Figures 2, 5, and 6 using frozen products only.

Prepare compact plotting data (requires pyarrow):
  python redraw_figures.py --prepare --source /path/to/paper_revision_20261015
Redraw from the distributed compact data (no pyarrow required):
  python redraw_figures.py

This script does not import or run the production pipeline or replay linking.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

COMMIT = "2f5f71d5a074c4cfd712b0e3fbef240c40745e2a"
BLUE, ORANGE, GREEN = "#24658a", "#bf5b25", "#287460"
INK, GRAY, PALE = "#1f2930", "#78858e", "#edf0f2"
MAG = r"$g_{\rm aper}$ [mag]"
RATE = r"Predicted rate [arcsec h$^{-1}$]"
MAG_EDGES = np.round(np.arange(9.5, 22.0001, .1), 8)
ERROR_EDGES = 10 ** np.round(np.arange(-3.2, .48001, .04), 8)
RADIAL_EDGES = np.round(np.arange(0, 1.0001, .02), 8)
RET_MAG_EDGES = np.round(np.arange(10, 21.0001, .25), 8)
RET_RATE_EDGES = np.arange(3., 63.0001, 3.)
MIN_CELL_N = 30


def checksum(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_grid(path, xedges, yedges, **arrays):
    ix, iy = np.meshgrid(np.arange(len(xedges)-1), np.arange(len(yedges)-1), indexing="ij")
    d = {"x_lower": xedges[ix].ravel(), "x_upper": xedges[ix+1].ravel(),
         "y_lower": yedges[iy].ravel(), "y_upper": yedges[iy+1].ravel()}
    d.update({key: np.asarray(values).ravel() for key, values in arrays.items()})
    pd.DataFrame(d).to_csv(path, index=False)


def spherical_separation(ra1, dec1, ra2, dec2):
    """Stable great-circle separation, identical definition to frozen script."""
    a1, d1, a2, d2 = np.deg2rad([ra1, dec1, ra2, dec2])
    dl = a2-a1
    num = np.hypot(np.cos(d2)*np.sin(dl),
                   np.cos(d1)*np.sin(d2)-np.sin(d1)*np.cos(d2)*np.cos(dl))
    den = np.sin(d1)*np.sin(d2)+np.cos(d1)*np.cos(d2)*np.cos(dl)
    return np.arctan2(num, den) * (180/np.pi) * 3600


def prepare(source, out):
    source, out = Path(source), Path(out)
    data = out / "data"
    data.mkdir(parents=True, exist_ok=True)
    source_files = ["inputs/known_matched.parquet", "inputs/blinded_known_recovery_by_object.parquet",
                    "tables/exposure_ledger.csv.gz", "tables/stellar_repeatability_vs_magnitude.csv",
                    "figure_data/fig05_photometry_vs_mag.csv", "figure_data/fig05_astrometry_vs_mag.csv",
                    "figure_data/fig05_astrometry_by_night.csv", "figure_data/fig05_limiting_magnitude.csv",
                    "figure_data/fig06_retention_frozen_bins.csv"]
    provenance = {name: {"sha256": checksum(source/name), "bytes": (source/name).stat().st_size}
                  for name in source_files}
    ledger = pd.read_csv(source/"tables/exposure_ledger.csv.gz",
                         usecols=["night", "file_name", "noon_night", "quality", "limiting_mag"],
                         dtype={"night": str, "noon_night": str})
    raw = pd.read_parquet(source/"inputs/known_matched.parquet")
    raw["file_name"] = raw.source_file.astype("string").str.replace(r"_cat\.fits(?:\.gz)?$", ".fits", regex=True)
    raw = raw.merge(ledger, on=["night", "file_name"], how="left", validate="many_to_one")
    assert raw.quality.notna().all(), "Unmapped association row"
    q = raw.loc[raw.quality].copy()
    assert len(q) > 0
    q["dra_arcsec"] = (((q.obs_ra_win_deg-q.pred_ra_deg+180) % 360)-180)*np.cos(np.deg2rad(q.pred_dec_deg))*3600
    q["ddec_arcsec"] = (q.obs_dec_win_deg-q.pred_dec_deg)*3600
    q["radial_arcsec"] = spherical_separation(q.pred_ra_deg, q.pred_dec_deg, q.obs_ra_win_deg, q.obs_dec_win_deg)
    assert q.radial_arcsec.between(0, 1.0000001).all()
    phot = q.loc[q.flag.eq(0) & q.magerr_aper4.gt(0) & q.flux_aper4.gt(0) & q.fluxerr_aper4.gt(0)]
    assert len(phot) > 0
    phot_valid_n=len(phot)
    phot=phot.loc[phot.mag_aper4.between(MAG_EDGES[0],MAG_EDGES[-1],inclusive="left") &
                  phot.magerr_aper4.between(ERROR_EDGES[0],ERROR_EDGES[-1],inclusive="left")].copy()

    bcounts = np.histogram2d(phot.mag_aper4, phot.magerr_aper4, bins=[MAG_EDGES, ERROR_EDGES])[0].astype(int)
    ccounts = np.histogram2d(q.mag_aper4, q.radial_arcsec, bins=[MAG_EDGES, RADIAL_EDGES])[0].astype(int)
    assert bcounts.sum() == len(phot)
    coutside = q.loc[~q.mag_aper4.between(MAG_EDGES[0], MAG_EDGES[-1], inclusive="left")]
    assert len(coutside) + ccounts.sum() == len(q)
    coutside[["night", "noon_night", "source_file", "obj_id", "asteroid_number", "asteroid_name", "mag_aper4",
              "dra_arcsec", "ddec_arcsec", "radial_arcsec"]].to_csv(data/"fig05_c_outside_magnitude_range.csv", index=False)
    save_grid(data/"fig05_b_counts.csv.gz", MAG_EDGES, ERROR_EDGES, measurements=bcounts)
    save_grid(data/"fig05_c_counts.csv.gz", MAG_EDGES, RADIAL_EDGES, measurements=ccounts)

    photo = pd.read_csv(source/"figure_data/fig05_photometry_vs_mag.csv")
    photo = photo.loc[photo.metric.eq("magerr_aper4")].copy()
    for i,row in photo.iterrows():
        values=phot.loc[phot.mag_aper4.ge(row.g_aper_low)&phot.mag_aper4.lt(row.g_aper_high),"magerr_aper4"]
        photo.loc[i,"n"]=len(values)
        for key,prob in [("median",.5),("p16",.16),("p84",.84)]:photo.loc[i,key]=values.quantile(prob) if len(values) else np.nan
    astro = pd.read_csv(source/"figure_data/fig05_astrometry_vs_mag.csv")
    astro = astro.loc[astro.metric.eq("radial_arcsec") & astro.radius_arcsec.eq(1)].copy()
    for frame, values, column in [(photo, phot, "magerr_aper4"), (astro, q, "radial_arcsec")]:
        for row in frame.itertuples():
            s = values.loc[values.mag_aper4.ge(row.g_aper_low) & values.mag_aper4.lt(row.g_aper_high), column]
            assert len(s) == row.n
            if len(s):
                assert np.isclose(s.median(), row.median, rtol=0, atol=2e-8)
    photo.to_csv(data/"fig05_b_frozen_trends.csv", index=False)
    astro.to_csv(data/"fig05_c_frozen_trends.csv", index=False)
    stars = pd.read_csv(source/"tables/stellar_repeatability_vs_magnitude.csv")
    stars.to_csv(data/"fig05_b_stellar_repeatability.csv", index=False)

    nightly = q.groupby("noon_night").agg(n=("radial_arcsec", "size"),
                    dra_median=("dra_arcsec", "median"), ddec_median=("ddec_arcsec", "median")).reset_index(names="night")
    frozen_nights = pd.read_csv(source/"figure_data/fig05_astrometry_by_night.csv", dtype={"night": str})
    for metric, col in [("dra_arcsec", "dra_median"), ("ddec_arcsec", "ddec_median")]:
        check = nightly.merge(frozen_nights.loc[frozen_nights.metric.eq(metric), ["night", "n", "median"]],
                              on="night", suffixes=("", "_frozen"), validate="one_to_one")
        assert (check.n == check.n_frozen).all()
        assert np.allclose(check[col], check["median"], atol=1e-10, rtol=0)
    nightly.to_csv(data/"fig05_d_nightly_offsets_and_counts.csv", index=False)
    depths = ledger.loc[ledger.quality & ledger.limiting_mag.gt(0), "limiting_mag"]
    assert len(depths) > 0 and np.isfinite(depths).all()
    depth_edges = np.round(np.arange(14.4, 23.2001, .2), 8)
    depth_counts = np.histogram(depths, depth_edges)[0]
    pd.DataFrame({"lower": depth_edges[:-1], "upper": depth_edges[1:], "count": depth_counts}).to_csv(
        data/"fig05_a_model_depth_counts.csv", index=False)

    obj = pd.read_parquet(source/"inputs/blinded_known_recovery_by_object.parquet")
    e = obj.loc[obj.eligible_ge3_and_rate].copy()
    assert len(e)>0 and e.strict_linked.sum()<=len(e)
    assert not e.duplicated(["night", "object_key"]).any()
    valid_mag = e.median_measured_mag.between(RET_MAG_EDGES[0], RET_MAG_EDGES[-1], inclusive="left")
    excluded = e.loc[~valid_mag]
    assert not excluded.median_measured_mag.between(10,21,inclusive="left").any()
    excluded[["night", "object_key", "median_measured_mag", "first_last_rate_arcsec_per_hour",
              "detected_1p5arcsec_n", "strict_linked"]].to_csv(data/"fig06_magnitude_excluded.csv", index=False)
    ngrid = np.histogram2d(e.median_measured_mag, e.first_last_rate_arcsec_per_hour,
                          bins=[RET_MAG_EDGES, RET_RATE_EDGES])[0].astype(int)
    kgrid = np.histogram2d(e.loc[e.strict_linked, "median_measured_mag"],
                          e.loc[e.strict_linked, "first_last_rate_arcsec_per_hour"],
                          bins=[RET_MAG_EDGES, RET_RATE_EDGES])[0].astype(int)
    assert ngrid.sum() == valid_mag.sum() and kgrid.sum() == e.loc[valid_mag,"strict_linked"].sum()
    fraction = np.divide(kgrid, ngrid, out=np.full(ngrid.shape, np.nan), where=ngrid>0)
    save_grid(data/"fig06_mag_rate_grid.csv.gz", RET_MAG_EDGES, RET_RATE_EDGES,
              eligible_n=ngrid, retained_n=kgrid, retention=fraction, shown_retention=ngrid>=MIN_CELL_N)
    frozen_bins = pd.read_csv(source/"figure_data/fig06_retention_frozen_bins.csv")
    one = frozen_bins.loc[frozen_bins.stage.eq("strict_linked") & frozen_bins.metric.isin([
        "measured_aperture_magnitude", "first_last_rate_arcsec_per_hour"])].copy()
    for metric, col, n_expected, k_expected in [
            ("measured_aperture_magnitude", "median_measured_mag", int(valid_mag.sum()), int(e.loc[valid_mag,"strict_linked"].sum())),
            ("first_last_rate_arcsec_per_hour", "first_last_rate_arcsec_per_hour", len(e), int(e.strict_linked.sum()))]:
        sub = one.loc[one.metric.eq(metric)]
        assert sub.eligible_object_n.sum() == n_expected and sub.survived_object_n.sum() == k_expected
        for row in sub.itertuples():
            s = e.loc[e[col].ge(row.lower) & e[col].lt(row.upper)]
            assert len(s) == row.eligible_object_n and int(s.strict_linked.sum()) == row.survived_object_n
    one.to_csv(data/"fig06_frozen_1d_bins.csv", index=False)
    np.savez_compressed(data/"density_grids.npz", mag_edges=MAG_EDGES, error_edges=ERROR_EDGES,
        radial_edges=RADIAL_EDGES, photo_n=bcounts, astro_n=ccounts,
        ret_mag_edges=RET_MAG_EDGES, ret_rate_edges=RET_RATE_EDGES, ret_n=ngrid, ret_k=kgrid)
    summary = {
        "method_reference_commit": COMMIT, "scientific_interval": "2025-11-15T04:00:00Z <= midpoint < 2026-11-15T04:00:00Z", "source_files": provenance,
        "fig05": {"depth_exposures": len(depths), "depth_unavailable_exposures": int(ledger.quality.sum()-len(depths)), "depth_median": float(depths.median()),
            "depth_p25": float(depths.quantile(.25)), "depth_p75": float(depths.quantile(.75)),
            "photometry_valid_positive_n": phot_valid_n, "photometry_display_excluded_n": phot_valid_n-len(phot), "photometry_selected_n": len(phot), "photometry_displayed_n": int(bcounts.sum()),
            "astrometry_selected_n": len(q), "astrometry_displayed_n": int(ccounts.sum()),
            "astrometry_magnitude_outside_n": len(coutside), "magnitude_99_n": int(coutside.mag_aper4.eq(99).sum()),
            "other_outside_magnitudes": coutside.loc[coutside.mag_aper4.ne(99), "mag_aper4"].tolist(),
            "radial_median_arcsec": float(q.radial_arcsec.median()), "night_n": len(nightly),
            "night_association_count_min": int(nightly.n.min()), "night_association_count_max": int(nightly.n.max()),
            "density_mag_width": .1, "density_log10_error_width": .04, "density_radial_width_arcsec": .02,
            "trend_mag_width": .5, "trend_min_n": 30, "stellar_min_pairs": 50,
            "night_definition": "local-noon observing night from frozen exposure ledger",
            "night_counts_definition": "all primary 1 arcsec associations in validated exposures, including records outside magnitude panel"},
        "fig06": {"eligible_all_n": len(e), "retained_all_n": int(e.strict_linked.sum()),
            "overall_fraction": float(e.strict_linked.mean()), "magnitude_eligible_n": int(ngrid.sum()),
            "magnitude_retained_n": int(kgrid.sum()), "magnitude_excluded_n": len(excluded),
            "magnitude_excluded_retained_n": int(excluded.strict_linked.sum()),
            "excluded_mag_min": float(excluded.median_measured_mag.min()),
            "excluded_mag_max": float(excluded.median_measured_mag.max()),
            "grid_mag_width": .25, "grid_rate_width_arcsec_h": 3., "retention_min_n_per_cell": MIN_CELL_N,
            "retention_visible_cells": int((ngrid>=MIN_CELL_N).sum()),
            "retention_masked_nonempty_cells": int(((ngrid>0)&(ngrid<MIN_CELL_N)).sum()),
            "retention_visible_objects": int(ngrid[ngrid>=MIN_CELL_N].sum()),
            "denominator_definition": "annual physical-night/object identity with >=3 distinct 1.5 arcsec detections and inclusive 3--63 arcsec/h first-to-last predicted rate",
            "numerator_definition": "at least one frozen strict link whose >=3 members all map unambiguously to the same known identity",
            "one_dimensional_intervals": "frozen 95 percent Wilson intervals; descriptive binomial intervals",
            "not_measured": "image-level/L2 detection completeness or independently resolved losses at intermediate pipeline stages"}}
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({"fig05": summary["fig05"], "fig06": summary["fig06"]}, indent=2))


def style(statistical=False):
    plt.rcParams.update({"font.family": "serif" if statistical else "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
        "mathtext.fontset": "stix" if statistical else "dejavusans",
        "font.size": 9, "axes.titlesize": 9.2, "axes.labelsize": 9, "xtick.labelsize": 8.3,
        "ytick.labelsize": 8.3, "legend.fontsize": 8, "axes.linewidth": .65,
        "axes.spines.top": False, "axes.spines.right": False, "xtick.direction": "out",
        "ytick.direction": "out", "xtick.major.width": .65, "ytick.major.width": .65,
        "xtick.major.size": 3, "ytick.major.size": 3, "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK, "axes.edgecolor": INK,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "figure.facecolor": "white", "legend.frameon": False})


def title(ax, letter, text):
    ax.set_title(f"({letter}) {text}", loc="left", pad=8)


def save(fig, out, stem, svg=False):
    for ext in ["pdf", "png"] + (["svg"] if svg else []):
        fig.savefig(Path(out)/f"{stem}.{ext}", dpi=320)
    plt.close(fig)


def draw_flow(out):
    fig = plt.figure(figsize=(6.6, 2.88))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, 6.6), ylim=(.27, 3.15))
    ax.set_axis_off()
    ax.add_patch(FancyBboxPatch((.08, .38), 1.32, 2.62, boxstyle="round,pad=0.015,rounding_size=.035",
                               facecolor="#f2f4f5", edgecolor="none"))
    ax.add_patch(FancyBboxPatch((1.62, 2.14), 4.90, .71, boxstyle="round,pad=.01,rounding_size=.035",
                               facecolor="#f1f6f9", edgecolor="none"))
    ax.add_patch(FancyBboxPatch((1.62, .52), 4.90, .91, boxstyle="round,pad=.01,rounding_size=.035",
                               facecolor="#f0f6f3", edgecolor="none"))
    ax.text(.74, 2.87, "SURVEY &\nREDUCTION", ha="center", va="center", fontsize=8.2, weight="bold", color=GRAY)
    ax.text(1.67, 2.99, "KNOWN-OBJECT ASSOCIATION", fontsize=8.4, weight="bold", color=BLUE)
    ax.text(1.67, 1.61, "INDEPENDENT MOVING-SOURCE SEARCH", fontsize=8.4, weight="bold", color=GREEN)
    rows = [
        ("plan", .74, 2.36, 1.10, .46, "Field plan\n& cadence", GRAY),
        ("observe", .74, 1.73, 1.10, .46, "Repeat imaging\n30 s white light", GRAY),
        ("l2", .74, .98, 1.10, .56, "L1 calibration\nL2 catalogs", GRAY),
        ("ephemeris", 2.10, 2.49, .86, .48, "ALEPH\npredictions", BLUE),
        ("primary", 3.70, 2.49, 1.20, .48, "1″ primary\nassociation", BLUE),
        ("known", 5.45, 2.49, 1.82, .48, "Known measurements\n& ADES records", BLUE),
        ("static", 2.10, .98, .86, .66, "Gaia / static\nrejection", GREEN),
        ("link", 3.10, .98, .86, .66, "Tracklets\n& linkages", GREEN),
        ("shortarc", 4.10, .98, .86, .66, "Short-arc\nconsistency", GREEN),
        ("veto", 5.10, .98, .86, .66, "1.5″ known\nobject veto", GREEN),
        ("review", 6.10, .98, .86, .66, "Image review\nSingle-night\ncandidates", GREEN),
    ]
    nodes = {r[0]: r for r in rows}
    edges = []
    def arrow(points, dashed=False, color=GRAY):
        vertices = np.asarray(points)
        ax.add_patch(FancyArrowPatch(path=MplPath(vertices), arrowstyle="-|>", mutation_scale=7,
            linewidth=.8, color=color, linestyle=(0, (3, 2)) if dashed else "solid", zorder=2))
    def right_to_left(src, dst):
        a,b = nodes[src],nodes[dst]
        arrow([(a[1]+a[3]/2+.018,a[2]),(b[1]-b[3]/2-.025,b[2])])
        edges.append({"source":src,"target":dst,"meaning":"processing"})
    for src,dst in [("plan","observe"),("observe","l2")]:
        a,b=nodes[src],nodes[dst]
        arrow([(a[1],a[2]-a[4]/2-.025),(b[1],b[2]+b[4]/2+.025)])
        edges.append({"source":src,"target":dst,"meaning":"processing"})
    arrow([(1.31,.98),(1.50,.98),(1.50,2.49),(1.645,2.49)])
    arrow([(1.31,.98),(1.645,.98)])
    edges.extend([{"source":"l2","target":"ephemeris","meaning":"exposure and catalog input"},
                  {"source":"l2","target":"static","meaning":"catalog input"}])
    for pair in [("ephemeris","primary"),("primary","known"),("static","link"),("link","shortarc"),
                 ("shortarc","veto"),("veto","review")]:
        right_to_left(*pair)
    arrow([(2.10,2.225),(2.10,1.94),(5.10,1.94),(5.10,1.335)],dashed=True,color=BLUE)
    edges.append({"source":"ephemeris","target":"veto","meaning":"auxiliary 1.5 arcsec known associations; identity classification after moving-source search"})
    ax.text(3.68,2.035,"Auxiliary 1.5″ associations",ha="center",va="center",color=BLUE,fontsize=8.0)
    for name,x,y,w,h,label,color in rows:
        ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=.012,rounding_size=.027",
            facecolor="white",edgecolor=color,linewidth=.75,zorder=3))
        ax.text(x,y,label,ha="center",va="center",fontsize=8.25 if name!="review" else 8.0,
                linespacing=1.18,color=INK,zorder=4)
    pd.DataFrame(rows,columns=["node","x_in","y_in","width_in","height_in","label","color"]).to_csv(
        Path(out)/"data/fig02_nodes.csv",index=False)
    pd.DataFrame(edges).to_csv(Path(out)/"data/fig02_edges.csv",index=False)
    save(fig,out,"fig02_workflow",svg=True)


def draw_quality(out, summary):
    out=Path(out); data=out/"data"; grids=np.load(data/"density_grids.npz")
    fig=plt.figure(figsize=(6.6,5.85))
    gs=fig.add_gridspec(2,2,left=.105,right=.902,bottom=.09,top=.946,wspace=.34,hspace=.36)
    a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1]); c=fig.add_subplot(gs[1,0])
    dg=gs[1,1].subgridspec(2,1,height_ratios=[3.3,1],hspace=.06)
    d=fig.add_subplot(dg[0]); dn=fig.add_subplot(dg[1],sharex=d)
    stats=summary["fig05"]
    depth=pd.read_csv(data/"fig05_a_model_depth_counts.csv")
    edges=np.r_[depth.lower,depth.upper.iloc[-1]]
    a.stairs(depth["count"]/1000,edges,fill=True,color=BLUE,alpha=.18,lw=0)
    a.stairs(depth["count"]/1000,edges,color=BLUE,lw=1.0)
    a.axvspan(stats["depth_p25"],stats["depth_p75"],color=ORANGE,alpha=.10,zorder=0)
    a.axvline(stats["depth_median"],color=ORANGE,lw=1.2,zorder=3)
    a.set(xlim=(14.4,23.2),ylim=(0,float(depth["count"].max()/1000)*1.16),xticks=[15,17,19,21],
          xlabel=r"Catalog $5\sigma$ model depth [mag]",ylabel=r"Exposures [$10^3$]")
    a.text(.04,.95,f"{stats['depth_exposures']:,} exposures\nMedian {stats['depth_median']:.2f} mag\nIQR {stats['depth_p25']:.2f}–{stats['depth_p75']:.2f}",transform=a.transAxes,
           va="top",ha="left",fontsize=8.05,linespacing=1.55)
    title(a,"a","Single-exposure model depth")
    density_cmap=LinearSegmentedColormap.from_list("count_gray",["#f1f2f3","#bac1c6","#6b7882","#26343f"])
    density_cmap.set_bad("white")
    maxcount=max(grids["photo_n"].max(),grids["astro_n"].max())
    norm=LogNorm(1,10**np.ceil(np.log10(maxcount)))
    for ax,counts,ys in [(b,grids["photo_n"],grids["error_edges"]),(c,grids["astro_n"],grids["radial_edges"])]:
        mesh=ax.pcolormesh(grids["mag_edges"],ys,np.ma.masked_equal(counts.T,0),
                          cmap=density_cmap,norm=norm,shading="flat",rasterized=False)
        ax.set_xlim(9.5,22)
        ax.set_xticks([10,13,16,19,22])
        ax.set_xlabel(MAG)
    ph=pd.read_csv(data/"fig05_b_frozen_trends.csv")
    ph=ph.loc[ph.n.ge(30)]
    x=(ph.g_aper_low+ph.g_aper_high)/2
    halo=[pe.Stroke(linewidth=2.6,foreground="white",alpha=.92),pe.Normal()]
    line,=b.plot(x,ph["median"],color=ORANGE,lw=1.35,zorder=4,label="Catalog median")
    line.set_path_effects(halo)
    for col in ["p16","p84"]:
        b.plot(x,ph[col],color=ORANGE,lw=.7,ls=(0,(3,2)),alpha=.9,zorder=3)
    stars=pd.read_csv(data/"fig05_b_stellar_repeatability.csv")
    stars=stars.loc[stars.star_pairs.ge(50)]
    b.plot((stars.g_aper_low+stars.g_aper_high)/2,stars.empirical_per_epoch_sigma_mag,
           "o",ms=3.5,mec="white",mew=.45,color=BLUE,zorder=5,label="Stellar repeatability")
    b.set(yscale="log",ylim=(grids["error_edges"][0],grids["error_edges"][-1]),
          ylabel="Aperture error / scatter [mag]")
    b.yaxis.set_major_locator(LogLocator(base=10,numticks=5)); b.yaxis.set_minor_formatter(NullFormatter())
    b.legend(loc="upper left",fontsize=8.0,handlelength=1.4,borderpad=.1,labelspacing=.35)
    title(b,"b","Aperture photometry")
    ar=pd.read_csv(data/"fig05_c_frozen_trends.csv")
    ar=ar.loc[ar.n.ge(30)]; x=(ar.g_aper_low+ar.g_aper_high)/2
    line,=c.plot(x,ar["median"],color=ORANGE,lw=1.35,zorder=4)
    line.set_path_effects(halo)
    for col in ["p16","p84"]:
        c.plot(x,ar[col],color=ORANGE,lw=.8,ls=(0,(3,2)),zorder=3)
    c.set(ylim=(0,1),yticks=[0,.2,.4,.6,.8,1],ylabel="Radial O−C [arcsec]")
    c.text(.97,.96,"1″ association selection",transform=c.transAxes,ha="right",va="top",fontsize=8.0)
    c.legend(handles=[Line2D([],[],color=ORANGE,lw=1.35,label="Median"),
                      Line2D([],[],color=ORANGE,lw=.8,ls="--",label="16th–84th percentiles")],
             loc="upper left",bbox_to_anchor=(0,.88),fontsize=8.0,handlelength=1.4,borderpad=.1,labelspacing=.35)
    title(c,"c","Radial astrometric residual")
    nights=pd.read_csv(data/"fig05_d_nightly_offsets_and_counts.csv",dtype={"night":str})
    t=pd.to_datetime(nights.night,format="%Y%m%d")
    d.scatter(t,nights.dra_median,s=9,color=BLUE,linewidths=0,label=r"$\Delta\alpha\cos\delta$",zorder=3)
    d.scatter(t,nights.ddec_median,s=10,color=ORANGE,marker="x",linewidths=.7,label=r"$\Delta\delta$",zorder=3)
    d.axhline(0,color=GRAY,lw=.65,zorder=1)
    d.set(ylim=(-.15,.49),yticks=[-.1,0,.1,.2,.3,.4],ylabel="Median O−C [arcsec]")
    d.legend(loc="upper right",fontsize=8,ncol=1,handlelength=.8,borderpad=.1,labelspacing=.3)
    d.tick_params(axis="x",labelbottom=False)
    dn.bar(t,nights.n/1000,width=1.25,color="#b9c1c6",linewidth=0)
    dn.set(ylim=(0,max(14,float(nights.n.max()/1000)*1.1)),yticks=[0,10],ylabel=r"$N$ [$10^3$]")
    dn.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[12,2,4,6,8,10]))
    dn.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    dn.set_xlabel("Observing night (2025–2026)")
    dn.set_xlim(pd.Timestamp("2025-11-15"),pd.Timestamp("2026-11-15"))
    title(d,"d",f"Nightly offsets: {stats['night_n']} nights")
    cbax=fig.add_axes([.923,.28,.018,.47])
    cb=fig.colorbar(mesh,cax=cbax)
    cb.solids.set_rasterized(False)
    cb.ax.set_title(r"$N$ / cell",fontsize=8.0,pad=8)
    cb.ax.tick_params(labelsize=8,length=2,width=.6)
    cb.outline.set_linewidth(.55)
    for ax in [a,d]:
        ax.grid(axis="y",color="#cdd2d6",alpha=.4,lw=.45,zorder=0)
        ax.set_axisbelow(True)
    save(fig,out,"fig05_measurement_quality")


def draw_retention(out, summary):
    out=Path(out); data=out/"data"; grids=np.load(data/"density_grids.npz")
    stats=summary["fig06"]
    fig=plt.figure(figsize=(6.6,6.15))
    gs=fig.add_gridspec(2,2,left=.105,right=.935,bottom=.08,top=.86,wspace=.30,hspace=.44,
                        height_ratios=[1,1.05])
    a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1])
    a_counts=grids["ret_n"].T; k=grids["ret_k"].T
    count_cmap=LinearSegmentedColormap.from_list("eligible_count",["#f0f2f3","#bac4cc","#6a7e8e","#203c51"])
    count_cmap.set_bad("white")
    count_mesh=a.pcolormesh(grids["ret_mag_edges"],grids["ret_rate_edges"],np.ma.masked_equal(a_counts,0),
                           cmap=count_cmap,norm=LogNorm(1,max(10000,10**np.ceil(np.log10(a_counts.max())))),shading="flat",rasterized=False)
    fractions=np.divide(k,a_counts,out=np.full(k.shape,np.nan),where=a_counts>0)*100
    frac_cmap=plt.get_cmap("Blues").copy();frac_cmap.set_bad("#e5e8ea")
    retention_mesh=b.pcolormesh(grids["ret_mag_edges"],grids["ret_rate_edges"],
        np.ma.masked_where(a_counts<MIN_CELL_N,fractions),cmap=frac_cmap,norm=Normalize(0,100),
        shading="flat",rasterized=False)
    for ax in [a,b]:
        ax.set(xlim=(10,21),ylim=(3,63),xticks=[10,12,14,16,18,20],yticks=[3,15,30,45,60],
               xlabel=r"Median $g_{\rm aper}$ [mag]")
    a.set_ylabel(RATE)
    b.set_ylabel(RATE)
    a.set_title("(a) Eligible object-nights / cell",loc="left",pad=35)
    b.set_title("(b) Conditional retention [%]",loc="left",pad=35)
    for ax,mesh,label in [(a,count_mesh,"Eligible object-nights per cell"),(b,retention_mesh,"Retained fraction [%]")]:
        pos=ax.get_position()
        cbax=fig.add_axes([pos.x0,pos.y1+.015,pos.width,.014])
        cb=fig.colorbar(mesh,cax=cbax,orientation="horizontal")
        cb.solids.set_rasterized(False)
        cb.ax.xaxis.set_ticks_position("top")
        cb.ax.tick_params(labelsize=8,length=2,width=.55,pad=1)
        cb.outline.set_linewidth(.45)
        # Panel titles state the scale and units; extra colorbar labels would crowd them.
    b.text(.03,.96,r"Gray: $N<30$",ha="left",va="top",transform=b.transAxes,fontsize=8,
           bbox={"facecolor":"white","edgecolor":"none","alpha":.88,"pad":2})
    one=pd.read_csv(data/"fig06_frozen_1d_bins.csv")
    for j,metric,xlabel,letter,heading,xlim in [
        (0,"measured_aperture_magnitude",r"Median $g_{\rm aper}$ [mag]","c","Retention by magnitude",(10,21)),
        (1,"first_last_rate_arcsec_per_hour",RATE,"d","Retention by angular rate",(3,63))]:
        inner=gs[1,j].subgridspec(2,1,height_ratios=[3.4,1],hspace=.07)
        ax=fig.add_subplot(inner[0]); den=fig.add_subplot(inner[1],sharex=ax)
        rows=one.loc[one.metric.eq(metric)&one.eligible_object_n.gt(0)]
        mid=(rows.lower+rows.upper)/2
        y=rows.survival_fraction*100
        err=np.array([rows.survival_fraction-rows.ci95_low,rows.ci95_high-rows.survival_fraction])*100
        ax.axhline(stats["overall_fraction"]*100,color=ORANGE,lw=.8,ls=(0,(4,2)),zorder=1)
        ax.errorbar(mid,y,xerr=(rows.upper-rows.lower)/2,yerr=err,fmt="o",ms=3.6,
                    color=BLUE,mec="white",mew=.35,ecolor=BLUE,elinewidth=.9,capsize=2,zorder=3)
        ax.set(xlim=xlim,ylim=(0,100),yticks=[0,25,50,75,100])
        ax.set_ylabel("Retention [%]")
        ax.grid(axis="y",color="#cdd2d6",alpha=.4,lw=.45,zorder=0)
        ax.tick_params(axis="x",labelbottom=False)
        ax.text(.035,.955,f"All eligible: {stats['overall_fraction']*100:.2f}%",color=ORANGE,ha="left",va="top",transform=ax.transAxes,fontsize=8)
        title(ax,letter,heading)
        den.bar(mid,rows.eligible_object_n/1000,width=(rows.upper-rows.lower)*.92,
                color="#bec7cd",edgecolor="white",linewidth=.35)
        den.set(ylabel=r"$N$ [$10^3$]",xlabel=xlabel)
        den.set_yticks([0,30] if j==0 else [0,20])
        den.set_ylim(0,float(rows.eligible_object_n.max()/1000)*1.12)
        if j==0:den.set_xticks([10,12,14,16,18,20])
        else:den.set_xticks([3,15,30,45,60])
    save(fig,out,"fig06_conditional_retention")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prepare",action="store_true",help="Prepare compact data from frozen input products before plotting")
    p.add_argument("--source",type=Path,help="Frozen paper_revision_20261015 directory; used only by --prepare")
    p.add_argument("--output",type=Path,default=Path(__file__).resolve().parent)
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if args.prepare:
        if args.source is None:p.error("--prepare requires --source")
        prepare(args.source,args.output)
    summary=json.loads((args.output/"summary.json").read_text())
    style()
    draw_flow(args.output)
    style(statistical=True)
    draw_quality(args.output,summary)
    draw_retention(args.output,summary)
    print("Wrote vector PDFs and 320-dpi PNGs for Figures 2, 5, and 6.")


if __name__ == "__main__":
    main()
