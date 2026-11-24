#!/usr/bin/env python3
"""Render the annual SHARP figures in the Pal 5 New visual style."""
import argparse
import io
import json

import healpy as hp
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator
import numpy as np
import pandas as pd
from PIL import Image
from scipy.signal import correlate2d

from summarize import EVIDENCE, ROOT

BLUE = "#1f77b4"
ORANGE = "#ff8c00"
INK = "#19252c"
MUTED = "#43535d"


def setup_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 13, "axes.labelsize": 17,
        "axes.titlesize": 18, "xtick.labelsize": 13, "ytick.labelsize": 13,
        "legend.fontsize": 12, "figure.dpi": 140, "axes.linewidth": 0.85,
        "axes.unicode_minus": False, "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK, "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    })


def save(fig, name):
    for extension in ("png", "pdf"):
        fig.savefig(ROOT / f"{name}.{extension}", dpi=300)
    plt.close(fig)


def all_sky_figure(stats):
    cov = np.load(EVIDENCE / "wcs_coverage.npz")
    nside = int(cov["nside"])
    all_counts = np.zeros(hp.nside2npix(nside), dtype=np.int32)
    all_counts[cov["ring_pixel"]] = cov["count"]
    # Mesh centres are sampled from the same equal-area grid used for the figure data.
    xedge = np.linspace(-np.pi, np.pi, 1901)
    yedge = np.linspace(-np.pi/2, np.pi/2, 951)
    x, y = np.meshgrid((xedge[1:]+xedge[:-1])/2, (yedge[1:]+yedge[:-1])/2)
    pix = hp.ang2pix(nside, np.rad2deg(-x+np.pi), np.rad2deg(y), lonlat=True)
    value = np.ma.masked_equal(all_counts[pix], 0)
    vhigh = int(all_counts.max())
    monthly = pd.read_csv(ROOT / "monthly_summary.csv").query("task == 'all_main'")
    fig = plt.figure(figsize=(12, 7.8))
    fig.text(.077, .958, "2025 Nov 16 – 2026 Nov 23", fontsize=13)
    metrics = (f"{stats['quality_nights']} nights   ·   {stats['quality_exposures']:,} exposures   ·   "
               f"{stats['quality_shutter_hours']:.1f} h   ·   {stats['area_once_deg2']:,.0f} deg²")
    fig.text(.937, .958, metrics, fontsize=13, ha="right")
    sky = fig.add_axes([.052, .354, .815, .57], projection="mollweide")
    color = plt.get_cmap("viridis").with_extremes(bad=(1, 1, 1, 0))
    m = sky.pcolormesh(xedge, yedge, value, cmap=color, norm=LogNorm(1, vhigh),
                      shading="flat", rasterized=True, zorder=1)
    sky.grid(True, linewidth=.45, alpha=.55, color="#6e7980")
    sky.set_longitude_grid(45)
    sky.set_latitude_grid(30)
    sky.set_xticklabels([])
    sky.set_yticklabels([])
    for deg in (45, 90, 135, 180, 225, 270, 315):
        sky.text(np.deg2rad(180-deg), np.deg2rad(-57), f"{deg}°", fontsize=10.5,
                 color=MUTED, ha="center", va="center", zorder=4)
    for deg in (-30, 30, 60):
        sky.annotate(f"{deg:+d}°", xy=(0, np.deg2rad(deg)), xytext=(-5, 2),
                     textcoords="offset points", ha="right", va="bottom", fontsize=9,
                     color=MUTED, zorder=4)
    eps = np.deg2rad(23.4392911)
    ra = np.linspace(0, 2*np.pi, 721)
    dec = np.arctan(np.tan(eps) * np.sin(ra))
    sky.plot(np.pi-ra, dec, lw=1.65, color=ORANGE, zorder=3)
    sky.legend(handles=[Line2D([0], [0], lw=1.8, color=ORANGE, label="Ecliptic")],
               loc="lower left", bbox_to_anchor=(.014, -.038), frameon=False, handlelength=1.9)
    fig.text(.46, .327, "R.A. (ICRS)", color=MUTED, ha="center", fontsize=12)
    fig.text(.077, .884, "(a)", fontsize=13)
    cax = fig.add_axes([.875, .492, .015, .265])
    cbar = fig.colorbar(m, cax=cax)
    ticks = [k for k in (1, 3, 10, 30, 100, 300, 1000) if k <= vhigh]
    cbar.locator = FixedLocator(ticks)
    cbar.formatter = FuncFormatter(lambda v, _: f"{v:g}")
    cbar.update_ticks()
    cbar.ax.minorticks_off()
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label("Coverage count (log)", fontsize=12, labelpad=7)
    ax = fig.add_axes([.102, .105, .831, .183])
    months = pd.to_datetime(monthly.month, format="%Y-%m")
    bars = ax.bar(np.arange(len(months)), monthly.quality_exposures, color=BLUE, width=.72, zorder=3)
    bars[-1].set_edgecolor(INK)
    bars[-1].set_linewidth(.8)
    labels = months.dt.strftime("%b").tolist()
    labels[0], labels[2] = "Nov ’25", "Jan ’26"
    ax.set_xticks(range(len(months)), labels)
    ax.set_xlim(-.62, len(months)-.38)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Exposures", fontsize=13, labelpad=8)
    fig.text(.077, .303, "(b)", fontsize=13)
    ax.tick_params(axis="both", labelsize=11)
    ax.yaxis.set_major_locator(MultipleLocator(3000))
    ax.grid(axis="y", color="#d5dce0", lw=.65)
    ax.set_axisbelow(True)
    ax.set_ylim(0, max(monthly.quality_exposures)*1.12)
    save(fig, "fig01_annual_coverage")


def standardized(frames):
    backgrounds = np.median(frames, axis=(1, 2))
    noise = 1.4826*np.median(np.abs(frames - backgrounds[:, None, None]), axis=(1, 2))
    return (frames - backgrounds[:, None, None]) / noise[:, None, None], backgrounds, noise


def stamps_qa(z, targets):
    """Use bright stationary stars to check registration after the TPV-to-TAN resampling."""
    corr = []
    first = np.clip(z[0], 0, 40)
    for frame in z[1:]:
        comparable = np.clip(frame, 0, 40)
        a = np.where(first[8:-8,8:-8] > 3, first[8:-8,8:-8], 0)
        b = np.where(comparable > 3, comparable, 0)
        cc = correlate2d(b, a, mode="valid")
        y, x = np.unravel_index(np.argmax(cc), cc.shape)
        vals = []
        for delta, m, idx in [(cc[y, :], x, 0), (cc[:, x], y, 1)]:
            if 0 < m < len(delta)-1:
                l, peak, r = delta[m-1:m+2]
                shift = .5*(l-r)/(l-2*peak+r)
            else:
                shift = 0
            vals.append(float(m + shift - 8))
        corr.append(vals)
    positions = []
    for frame, (x, y) in zip(z, targets):
        ix, iy = int(round(x)), int(round(y))
        core = frame[iy-1:iy+2, ix-1:ix+2]
        positions.append(float(core.max()))
    return corr, positions


def stamp(ax, pixels, target, timestamp, scale=True):
    top = np.arcsinh(np.clip((pixels+2.5)/7, 0, None))/np.arcsinh(22/7)
    ax.imshow(top, cmap="gray", vmin=0, vmax=1, origin="lower", interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    x, y = target
    ax.annotate("", xy=(x+2.5, y+2.5), xytext=(x+21, y+21),
                arrowprops=dict(arrowstyle="-|>", color="#ffd562", lw=1.6, mutation_scale=13))
    ax.text(.5, -.055, timestamp + " UTC", transform=ax.transAxes, ha="center", va="top", fontsize=14)
    if scale:
        ax.plot([12, 12+30/1.15], [13, 13], color="white", lw=1.5)
        ax.text(12+15/1.15, 17, "30″", color="white", ha="center", fontsize=11)
        ax.annotate("", xy=(137, 32), xytext=(137, 14),
                    arrowprops=dict(arrowstyle="-|>", color="white", lw=1.0))
        ax.annotate("", xy=(119, 14), xytext=(137, 14),
                    arrowprops=dict(arrowstyle="-|>", color="white", lw=1.0))
        ax.text(137, 34, "N", color="white", ha="center", fontsize=10)
        ax.text(116, 14, "E", color="white", ha="right", va="center", fontsize=10)


def moving_figure(data, meta):
    z, sky_levels, noise = standardized(data["asteroid"])
    targets = meta["target_xy"]
    times = [f["midpoint_utc"].split("T")[1] for f in meta["frames"]]
    corr, peaks = stamps_qa(z, targets)
    shift = np.linalg.norm(np.diff(np.asarray(targets), axis=0), axis=1)*meta["scale"]
    travel = np.linalg.norm(np.subtract(targets[-1], targets[0]))*meta["scale"]
    audit = dict(registration_static_star_offset_pixels=corr, target_max_local_sigma=peaks,
                 successive_motion_arcsec=shift.tolist(), first_to_last_arcsec=float(travel),
                 photometric_background_adu=sky_levels.tolist(), background_mad_sigma_adu=noise.tolist())
    assert np.max(np.abs(corr)) < 1.0 and min(peaks) > 5, audit
    (EVIDENCE / "image_qa.json").write_text(json.dumps(audit, indent=2) + "\n")
    fig = plt.figure(figsize=(12, 4.65))
    fig.text(.06, .934, "(282917) 2007 LA33", va="top", fontsize=20)
    fig.text(.938, .926, "2026 Sep 1", va="top", ha="right", fontsize=14, color=MUTED)
    for j in range(3):
        ax = fig.add_axes([.06 + j*.298, .143, .282, .727])
        stamp(ax, z[j], targets[j], times[j], scale=j == 0)
    save(fig, "fig02_asteroid_motion")
    gif = []
    for j in range(3):
        fig = plt.figure(figsize=(5.3, 5.7))
        fig.text(.5, .95, "(282917) 2007 LA33", ha="center", fontsize=16)
        fig.text(.5, .906, "2026 Sep 1", ha="center", color=MUTED, fontsize=12)
        ax = fig.add_axes([.10, .126, .80, .74])
        stamp(ax, z[j], targets[j], times[j])
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=150)
        plt.close(fig)
        buffer.seek(0)
        with Image.open(buffer) as img:
            gif.append(img.convert("RGB"))
    # One palette for the shared sky and the arrow, independently of the frame.
    sheet = Image.new("RGB", (gif[0].width, gif[0].height * len(gif)))
    for j, frame in enumerate(gif):
        sheet.paste(frame, (0, j * frame.height))
    palette = sheet.quantize(colors=128)
    indexed = [frame.quantize(palette=palette, dither=Image.Dither.NONE) for frame in gif]
    output = io.BytesIO()
    indexed[0].save(output, format="GIF", save_all=True, append_images=indexed[1:],
                    duration=[900, 900, 1300], loop=0, disposal=2, optimize=False)
    with Image.open(io.BytesIO(output.getvalue())) as check:
        assert check.n_frames == 3 and check.info["loop"] == 0
        for j in range(check.n_frames):
            check.seek(j)
            check.load()
    (ROOT / "fig02_asteroid_motion.gif").write_bytes(output.getvalue())


ORBIT_STYLES = {
    "Main belt": ("o", "#17becf", .18, 5),
    "Outer belt / Hilda": ("s", "#f58518", .70, 13),
    "Jupiter Trojan region": ("D", "#54a24b", .72, 13),
    "Inner belt / Hungaria": ("^", "#b279a2", .72, 15),
    "Mars-crosser": ("v", "#e45756", .78, 15),
    "NEO (q < 1.3 au)": ("*", "#ff9da6", .90, 29),
    "Distant (a ≥ 5.4 au)": ("P", "#d62728", .88, 30),
}


def orbital_figure(stats):
    orbits = pd.read_csv(ROOT / "known_object_summary.csv.gz")
    columns = ["a_au", "ecc", "inc_deg"]
    assert len(orbits) == stats["asteroids"] and orbits.identity.is_unique
    assert orbits[columns].notna().all().all()
    assert orbits.a_au.gt(0).all() and orbits.inc_deg.gt(0).all()
    counts = orbits.astorb_group.value_counts()
    assert set(counts.index) == set(ORBIT_STYLES)
    for group, count in counts.items():
        assert count == stats["asteroid_classes"][group]["objects"]

    orbit_style = {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelsize": 15.5,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "text.color": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
    }
    with plt.rc_context(orbit_style):
        fig, (left, right) = plt.subplots(1, 2, figsize=(12.4, 4.45))
        fig.subplots_adjust(left=.058, right=.987, bottom=.165, top=.905, wspace=.175)
        handles = []
        for group, (marker, color, alpha, size) in ORBIT_STYLES.items():
            members = orbits.loc[orbits.astorb_group.eq(group)]
            for ax, y in ((left, "ecc"), (right, "inc_deg")):
                points = ax.scatter(members.a_au, members[y], s=size, color=color, marker=marker,
                                    alpha=alpha, linewidths=0, rasterized=True)
                assert len(points.get_offsets()) == counts[group]
            handles.append(Line2D([], [], marker=marker, linestyle="none", markeredgewidth=0,
                                  color=color, label=group, markersize=8.2))

        xlim = (min(.55, .91*orbits.a_au.min()), max(240, 1.10*orbits.a_au.max()))
        for ax, tag, quantity in ((left, "a", "eccentricity"), (right, "b", "inclination")):
            ax.set(xscale="log", xlim=xlim, xlabel="Semimajor axis [au]")
            ax.set_title(f"({tag}) Semimajor axis and {quantity}", loc="left", fontsize=16.5, pad=7)
            ax.set_axisbelow(True)
            ax.grid(axis="y", linewidth=.65, color="#d9dde0", alpha=.7)
            ax.tick_params(which="major", length=5)
            ax.tick_params(which="minor", length=3)

        boundary = 1.3
        a = np.geomspace(boundary, xlim[1], 300)
        left.plot(a, 1-boundary/a, lw=1.05, ls="--", color="#64737e", label=r"$q = 1.3$ au")
        left.set(ylabel="Eccentricity", ylim=(-.05, 1.045))
        left.yaxis.set_major_locator(MultipleLocator(.2))
        left.legend(loc="lower right", frameon=False, fontsize=12, handlelength=2.3)

        right.set(yscale="log", ylabel="Inclination [deg]",
                  ylim=(min(.014, .8*orbits.inc_deg.min()), max(190, 1.2*orbits.inc_deg.max())))
        right.legend(handles=handles, loc="lower right", frameon=False, fontsize=10.7,
                     labelspacing=.37, handletextpad=.7, borderaxespad=.7)
        for ax, col in ((left, "ecc"), (right, "inc_deg")):
            assert orbits.a_au.between(*ax.get_xlim()).all()
            assert orbits[col].between(*ax.get_ylim()).all()
        save(fig, "fig03_orbital_distribution")

    (EVIDENCE / "orbital_figure_audit.json").write_text(json.dumps({
        "objects": len(orbits), "objects_per_panel": [len(orbits), len(orbits)],
        "classification": "Lowell ASTORB orbital-element groups",
        "objects_by_class": {name: int(counts[name]) for name in ORBIT_STYLES},
        "semimajor_axis_scale": "log", "inclination_scale": "log",
        "perihelion_boundary_au": boundary,
        "semimajor_axis_range_au": [float(orbits.a_au.min()), float(orbits.a_au.max())],
    }, ensure_ascii=False, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-only", action="store_true", help="Render figure 3 only.")
    options = parser.parse_args()
    setup_style()
    stats = json.loads((EVIDENCE / "annual_metrics.json").read_text())
    if options.orbit_only:
        orbital_figure(stats)
        return
    all_sky_figure(stats)
    with np.load(EVIDENCE / "image_pixels.npz") as data:
        meta = json.loads(str(data["metadata"]))
        (EVIDENCE / "image_sources.json").write_text(json.dumps(meta, indent=2) + "\n")
        moving_figure(data, meta["asteroid"])
    orbital_figure(stats)


if __name__ == "__main__":
    main()
