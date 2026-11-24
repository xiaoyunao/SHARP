#!/usr/bin/env python3
"""Draw the SHARP overview figures from the included tabular data."""
import astropy.units as u
import healpy as hp
import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import BarycentricMeanEcliptic, SkyCoord
from astropy.time import Time
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path
from matplotlib.ticker import MaxNLocator

from exposures import DATA, ROOT, TABLES

BLUE, ORANGE, GREEN, GRAY = "#31678c", "#bd652e", "#36866c", "#64737e"
MAG = r"$g_{\rm aper}$ [mag]"
FONTS = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]
ORBIT_STYLE = {
    "Main belt": ("o", "#17becf", .18, 5),
    "Outer belt / Hilda": ("s", "#f58518", .70, 12),
    "Jupiter Trojan region": ("D", "#54a24b", .72, 12),
    "Inner belt / Hungaria": ("^", "#b279a2", .72, 14),
    "Mars-crosser": ("v", "#e45756", .78, 14),
    "NEO (q < 1.3 AU)": ("*", "#ff9da6", .90, 25),
    "Distant (a ≥ 5.4 au)": ("P", "#d62728", .88, 24),
}
plt.rcParams.update({"font.family": "serif", "font.serif": FONTS, "font.size": 10,
                     "axes.titlesize": 11, "axes.labelsize": 10, "axes.linewidth": .8,
                     "figure.dpi": 150, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "legend.frameon": False})


def write(fig, name):
    path = ROOT / "figures" / name
    for fmt in ("pdf", "png"):
        fig.savefig(path.with_suffix(f".{fmt}"), dpi=300, bbox_inches="tight", pad_inches=.08)
    plt.close(fig)


def panel(ax, label):
    ax.set_title(label, loc="left", pad=5)
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", alpha=.18, linewidth=.55)


def hist(ax, path, color=BLUE, label=None, xmax=None):
    d = pd.read_csv(path)
    if xmax is not None:
        d = d[d.iloc[:, 0] < xmax]
    edges = np.r_[d.iloc[:, 0], d.iloc[-1, 1]]
    ax.stairs(d.iloc[:, 2], edges, color=color, fill=True, alpha=.25)
    ax.stairs(d.iloc[:, 2], edges, color=color, lw=1.0, label=label)


def flow():
    nodes = [
        ("plan", 0.5, .92, .39, .07, "Ecliptic field library and nightly plan", "plan"),
        ("observe", .5, .80, .39, .07, "Repeat imaging at Xinglong", "plan"),
        ("l2", .5, .68, .44, .08, "L1 image calibration and L2 source catalogs", "common"),
        ("known", .24, .52, .39, .10, "Known-object association\nAleph predictions and 1″ matching", "known"),
        ("gaia", .76, .52, .40, .10, "Independent moving-source search\nGaia and stationary-source rejection", "unknown"),
        ("ades", .24, .36, .39, .08, "Known measurements and ADES", "known"),
        ("orbit", .76, .36, .40, .10, "Tracklets, links, and\nshort-arc orbit consistency", "unknown"),
        ("mask", .76, .22, .40, .08, "Remove links with 1.5″ known matches", "unknown"),
        ("review", .76, .09, .40, .10, "Single-night candidates\nImage review and selected ADES", "unknown"),
    ]
    cols = ["node", "x", "y", "width", "height", "label", "branch"]
    nodes = pd.DataFrame(nodes, columns=cols)
    nodes.to_csv(DATA / "fig02_workflow_nodes.csv", index=False)
    edges = [
        ("plan", "observe", "solid"), ("observe", "l2", "solid"),
        ("l2", "known", "solid"), ("l2", "gaia", "solid"), ("known", "ades", "solid"),
        ("gaia", "orbit", "solid"), ("orbit", "mask", "solid"), ("mask", "review", "solid"),
        ("known", "mask", "dashed"),
    ]
    pd.DataFrame(edges, columns=["source", "target", "linestyle"]).to_csv(DATA / "fig02_workflow_edges.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.15, 5.15))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    node = nodes.set_index("node")
    colors = {"plan": ("#eaf0f4", BLUE), "common": ("#f2f3f4", GRAY),
              "known": ("#e6f1f7", BLUE), "unknown": ("#ecf3eb", GREEN)}
    for src, dest, style in edges:
        a, b = node.loc[src], node.loc[dest]
        if style == "dashed":
            path = Path([(a.x+a.width/2, a.y), (.49, a.y), (.49, b.y), (b.x-b.width/2, b.y)])
            ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=9,
                                         color=GRAY, linewidth=1.0, linestyle="--"))
            continue
        else:
            xytext = (a.x, a.y-a.height/2)
            xy = (b.x, b.y+b.height/2)
            connection = "arc3"
        ax.annotate("", xy=xy, xytext=xytext, arrowprops=dict(
            arrowstyle="-|>", color=GRAY, lw=1.05, linestyle=style, connectionstyle=connection))
    for n in nodes.itertuples():
        fill, stroke = colors[n.branch]
        ax.add_patch(FancyBboxPatch((n.x-n.width/2, n.y-n.height/2), n.width, n.height,
                     facecolor=fill, edgecolor=stroke, linewidth=.95,
                     boxstyle="round,pad=0.008,rounding_size=0.012", zorder=3))
        ax.text(n.x, n.y, n.label, fontsize=10.2, ha="center", va="center", zorder=4)
    ax.text(.27, .61, "Known-object branch", ha="center", color=BLUE, fontsize=9)
    ax.text(.74, .61, "Moving-source branch", ha="center", color=GREEN, fontsize=9)
    write(fig, "fig02_workflow")


def maps():
    v = pd.read_csv(DATA / "fig03_wcs_visits_icrs_nside512.csv.gz")
    z = pd.read_csv(DATA / "fig03_known_ecliptic_nside128.csv")
    vmap = np.zeros(hp.nside2npix(512), dtype="uint16")
    vmap[v.ring_pixel] = v.qualified_exposures
    zmap = np.zeros(hp.nside2npix(128), dtype="uint16")
    zmap[z.ring_pixel] = z.known_detections
    ecl = BarycentricMeanEcliptic(equinox=Time("J2000"))
    x = np.linspace(-180, 180, 1441)
    y = np.linspace(-90, 90, 721)
    lon, lat = np.meshgrid((x+180) % 360, y)
    p = SkyCoord(lon*u.deg, lat*u.deg, frame=ecl).icrs
    raw = vmap[hp.ang2pix(512, p.ra.deg, p.dec.deg, lonlat=True)]
    known = zmap[hp.ang2pix(128, lon, lat, lonlat=True)]
    equator = SkyCoord(np.linspace(0, 359.99, 3600)*u.deg, np.zeros(3600)*u.deg).transform_to(ecl)
    ex = np.deg2rad(equator.lon.deg-180)
    order = np.argsort(ex)
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.45), subplot_kw={"projection": "mollweide"},
                             layout="constrained")
    cmap = plt.get_cmap("rainbow").copy()
    cmap.set_bad(alpha=0)
    corners_x, corners_y = np.meshgrid(np.deg2rad(x), np.deg2rad(y))
    for ax, values, title, unit in zip(axes, [raw, known],
             ["(a) Qualified exposure footprints", "(b) Known-object associations"],
             ["Exposures per pixel (NSIDE 512)", "Associations per pixel (NSIDE 128)"]):
        mesh = ax.pcolormesh(corners_x, corners_y, np.ma.masked_equal(values, 0),
                             shading="nearest", cmap=cmap, norm=LogNorm(vmin=1, vmax=max(10, int(np.percentile(values[values>0], 99.5)))),
                             rasterized=True)
        ax.plot(ex[order], equator.lat.rad[order], color="#ed920c", lw=.9)
        ax.grid(alpha=.25)
        ax.set_xticks(np.deg2rad(np.arange(-120, 121, 60)))
        ax.set_xticklabels([r"60°", r"120°", r"180°", r"240°", r"300°"])
        ax.set_yticks(np.deg2rad([-60, -30, 0, 30, 60]))
        ax.tick_params(labelsize=8)
        ax.set_title(title, loc="left", pad=7)
        cb = fig.colorbar(mesh, ax=ax, orientation="horizontal", shrink=.66, pad=.08)
        cb.set_label(unit, fontsize=9)
    write(fig, "fig03_coverage_and_known_sky")


def activity():
    night = pd.read_csv(DATA / "fig04_nightly_activity.csv", dtype={"night": str})
    month = pd.read_csv(DATA / "fig04_monthly_activity.csv", dtype={"month": str})
    visits = pd.read_csv(DATA / "fig04_field_night_counts.csv")
    fig = plt.figure(figsize=(8.3, 6.1), layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[1.08, 1, 1])
    ax = fig.add_subplot(gs[0, :])
    time = pd.to_datetime(night.night)
    ax.bar(time, night.raw, color="#b9c2c7", width=.9, label="Frozen raw file inventory")
    ax.bar(time, night.qualified, color=BLUE, width=.9, label="Qualified exposures")
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.set_ylabel("Exposures per night")
    panel(ax, "(a) Nightly activity")
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    ax = fig.add_subplot(gs[1, 0])
    t = np.arange(len(month))
    ax.bar(t, month.raw, color="#b9c2c7")
    ax.bar(t, month.qualified, color=BLUE)
    ax.set_xticks(t, pd.to_datetime(month.month, format="%Y%m").dt.strftime("%b"))
    ax.set_ylabel("Exposures per month")
    panel(ax, "(b) Monthly activity")
    ax = fig.add_subplot(gs[1, 1])
    ax.bar([1, 2, 3], visits.field_nights, color=[GRAY, BLUE, GREEN])
    ax.set(xticks=[1, 2, 3], xticklabels=["1", "2", "≥3"], xlabel="Exposures per field-night", ylabel="Field-nights")
    ax.set_ylim(top=visits.field_nights.max()*1.15)
    ax.bar_label(ax.containers[0], labels=[f"{n:,}" for n in visits.field_nights], padding=3, fontsize=9)
    panel(ax, "(c) Qualified visits")
    for col, (suffix, label, xlabel, lim) in enumerate([
        ("adjacent", "(d) Adjacent visits", "Midpoint separation [min]", 120),
        ("span_ge3", "(e) First to last, ≥3 visits", "Midpoint span [min]", 180),
    ]):
        ax = fig.add_subplot(gs[2, col])
        hist(ax, DATA / f"fig04_{suffix}.csv", color=[BLUE, GREEN][col], xmax=lim)
        ax.set(xlabel=xlabel, ylabel=["Adjacent pairs", "Field-nights"][col], xlim=(0, lim))
        panel(ax, label)
    write(fig, "fig04_activity_and_cadence")


def quality():
    photo = pd.read_csv(DATA / "fig05_photometry_vs_mag.csv")
    astro = pd.read_csv(DATA / "fig05_astrometry_vs_mag.csv")
    nights = pd.read_csv(DATA / "fig05_astrometry_by_night.csv", dtype={"night": str})
    fig, axs = plt.subplots(2, 2, figsize=(8.3, 5.9), layout="constrained")
    a, b, c, d = axs.ravel()
    depth = pd.read_csv(DATA / "fig05_limiting_magnitude.csv")
    depth = depth[depth.lower.ge(14) & depth.lower.lt(23.8)]
    a.stairs(depth["count"], np.r_[depth.lower, depth.upper.iloc[-1]], color=BLUE, fill=True, alpha=.28)
    a.set(xlabel="Catalog 5σ model depth [mag]", ylabel="Qualified exposures", xlim=(17, 23.2))
    panel(a, "(a) Single-exposure depth estimate")
    q = photo.loc[photo.metric.eq("magerr_aper4") & photo.n.ge(30)]
    x = (q.g_aper_low+q.g_aper_high)/2
    b.fill_between(x, q.p16, q.p84, alpha=.24, color=GREEN, label="16th–84th percentiles")
    b.plot(x, q["median"], color=GREEN, lw=1.25, label="Median")
    stars = pd.read_csv(TABLES / "stellar_repeatability_vs_magnitude.csv")
    stars = stars.loc[stars.star_pairs.ge(50)]
    b.plot((stars.g_aper_low+stars.g_aper_high)/2, stars.empirical_per_epoch_sigma_mag,
           ".", markersize=3.2, color=BLUE, label="Stellar repeatability")
    b.set(xlabel=MAG, ylabel="Aperture error or scatter [mag]", yscale="log", xlim=(13, 21.1))
    b.legend(fontsize=8, loc="upper left")
    panel(b, "(b) Aperture errors and repeatability")
    labels = [(r"$\Delta\alpha\cos\delta$", "dra_arcsec", BLUE), (r"$\Delta\delta$", "ddec_arcsec", ORANGE)]
    for label, key, color in labels:
        sub = astro.loc[astro.metric.eq(key) & astro.radius_arcsec.eq(1) & astro.n.ge(30)]
        x = (sub.g_aper_low+sub.g_aper_high)/2
        c.fill_between(x, sub["median"]-sub.sigma_iqr, sub["median"]+sub.sigma_iqr, alpha=.18, color=color)
        c.plot(x, sub["median"], color=color, lw=1.2, label=label)
        q = nights.loc[nights.metric.eq(key) & nights.qualified_n.gt(0)]
        d.plot(pd.to_datetime(q.night), q["median"], ".", color=color, markersize=3, alpha=.78, label=label)
    c.axhline(0, color=GRAY, lw=.6)
    c.set(xlabel=MAG, ylabel="Signed O−C [arcsec]", xlim=(13, 21.1))
    c.legend(fontsize=8, loc="upper left", ncol=2)
    panel(c, "(c) 1″ associations by magnitude")
    d.axhline(0, color=GRAY, lw=.6)
    d.set_ylabel("Nightly median O−C [arcsec]")
    d.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    d.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    d.legend(fontsize=8, loc="lower right", ncol=2)
    panel(d, "(d) 1″ associations by night")
    write(fig, "fig05_measurement_quality")


def retention():
    data = pd.read_csv(DATA / "fig06_retention_frozen_bins.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.3, 3.35), layout="constrained", sharey=True)
    for ax, metric, label, title in zip(axes,
            ["measured_aperture_magnitude", "first_last_rate_arcsec_per_hour"],
            [r"Median $g_{\rm aper}$ [mag]", "Predicted first-to-last rate [arcsec/h]"],
            ["(a) Calibrated aperture magnitude", "(b) Predicted angular rate"]):
        q = data.loc[data.metric.eq(metric) & data.stage.eq("strict_linked") & data.eligible_object_n.gt(0)]
        x = (q.lower+q.upper)/2
        y = 100*q.survival_fraction
        ax.errorbar(x, y, xerr=(q.upper-q.lower)/2, yerr=(np.stack([q.survival_fraction-q.ci95_low,
                     q.ci95_high-q.survival_fraction])*100), linestyle="none", marker="o", markersize=4,
                     color=BLUE, ecolor="#8297a6", capsize=2, elinewidth=.9)
        ax.axhline(74.1585, lw=.85, ls="--", color=ORANGE, label="All eligible object-nights")
        ax.set(xlabel=label, ylim=(0, 100))
        panel(ax, title)
    axes[0].set_ylabel("Conditional retention [%]")
    axes[1].legend(fontsize=8, loc="lower right")
    write(fig, "fig06_conditional_retention")


def known():
    obj = pd.read_csv(DATA / "fig07_objects_orbits_revisits.csv.gz", low_memory=False)
    fig, axs = plt.subplots(2, 2, figsize=(8.4, 6.2), layout="constrained")
    a, b, c, d = axs.ravel()
    handles = []
    for label, (marker, color, alpha, size) in ORBIT_STYLE.items():
        q = obj.loc[obj.orbit_class.eq(label)]
        a.scatter(q.a_au, q.ecc, s=size, color=color, marker=marker, alpha=alpha,
                  linewidths=0, rasterized=True)
        b.scatter(q.a_au, q.inc_deg, s=size, color=color, marker=marker, alpha=alpha,
                  linewidths=0, rasterized=True)
        handles.append(Line2D([], [], marker=marker, linestyle="none", markeredgewidth=0,
                              color=color, label=label, markersize=6))
    curve = np.geomspace(1.3, max(obj.a_au), 120)
    a.plot(curve, 1-1.3/curve, linewidth=.7, linestyle="--", color=GRAY, label="$q=1.3$ au")
    a.set(xscale="log", xlabel="Semimajor axis [au]", ylabel="Eccentricity")
    a.legend(fontsize=8, loc="lower right")
    b.set(xscale="log", yscale="log", xlabel="Semimajor axis [au]", ylabel="Inclination [deg]")
    b.legend(handles=handles, fontsize=7, loc="best")
    panel(a, "(a) Semimajor axis and eccentricity")
    panel(b, "(b) Semimajor axis and inclination")
    for ax, col, xlabel, title in [(c, "nights", "Distinct nights per object", "(c) Nightly sampling"),
                                  (d, "span_days", "First-to-last span [days]", "(d) Multi-night baseline")]:
        values = obj[col] if col == "nights" else obj.loc[obj.nights.gt(1), col]
        bins = np.arange(.5, values.max()+1.5) if col == "nights" else np.geomspace(.8, values.max()*1.001, 21)
        ax.hist(values, bins=bins, color=BLUE, alpha=.45, edgecolor=BLUE, linewidth=.4)
        ax.set(xlabel=xlabel, ylabel="Known identities", yscale="log")
        if col == "span_days":
            ax.set_xscale("log")
        else:
            ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
        panel(ax, title)
    write(fig, "fig07_orbits_and_revisits")


def appendices():
    fig, axes = plt.subplots(1, 2, figsize=(8.3, 3.25), layout="constrained")
    ax, bx = axes
    for col, color in [("predicted", BLUE), ("measured", GREEN)]:
        t = pd.read_csv(DATA / f"figA1_{col}_rate.csv").iloc[1:-1]
        ax.stairs(t["count"], np.r_[t.lower, t.upper.iloc[-1]], color=color, label=col.capitalize())
    ax.set(xlabel="Angular rate [arcsec/h]", ylabel="Successive same-night detection pairs",
           xscale="log", yscale="log", xlim=(.1, 10000))
    ax.legend(fontsize=8)
    panel(ax, "(a) Sky-plane displacement rates")
    hist(bx, DATA / "figA1_elongation.csv", color=GREEN)
    bx.set(xlabel="Topocentric solar elongation [deg]", ylabel="Known-object associations", xlim=(0, 180))
    panel(bx, "(b) Solar elongation")
    write(fig, "figA1_motion_and_solar_elongation")

    data = pd.read_csv(TABLES / "retained_linkages.csv")
    data.to_csv(DATA / "figA2_candidate_rms_arc.csv", index=False)
    fig, ax = plt.subplots(figsize=(5.4, 3.45), layout="constrained")
    other = ~data.is_c2025_y1
    for use, mark, color, label in [
        (other, "o", BLUE, "Other retained linkages"), (~other, "*", ORANGE, "Consistent with C/2025 Y1"),
    ]:
        sub = data.loc[use]
        ax.scatter(sub.arc_min_midpoint, sub.rms_arcsec, marker=mark, s=32,
                   color=color, label=label, alpha=.78, linewidths=0)
    ax.legend(fontsize=8)
    ax.set(xlabel="Observed midpoint span [min]", ylabel="Short-arc fit RMS [arcsec]")
    panel(ax, "Post-review single-night linkages")
    write(fig, "figA2_candidate_short_arc_diagnostic")


if __name__ == "__main__":
    flow()
    maps()
    activity()
    quality()
    retention()
    known()
    appendices()
