#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Create a PPT-friendly orbit-confirm diagnostic figure for one linkage."
    )
    ap.add_argument("night", help="Night name, e.g. 20260220")
    ap.add_argument("--root", default="/pipeline/xiaoyunao/data/heliolincrr")
    ap.add_argument("--linkage-id", type=int, default=None, help="Specific linkage_id to visualize.")
    ap.add_argument(
        "--selection",
        choices=["auto", "lowest-rms", "max-obs"],
        default="auto",
        help="How to auto-pick a linkage when --linkage-id is omitted.",
    )
    ap.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Defaults to <root>/<night>/rr_links/orbit_confirm/diagnostics",
    )
    return ap.parse_args()


def choose_linkage(links: Table, residuals: Table, selection: str) -> int:
    resid_counts: dict[int, int] = {}
    for row in residuals:
        lid = int(row["linkage_id"])
        resid_counts[lid] = resid_counts.get(lid, 0) + 1

    good_rows = []
    for row in links:
        if not bool(row["fit_ok"]) or not bool(row["is_good"]):
            continue
        lid = int(row["linkage_id"])
        good_rows.append(
            (
                lid,
                int(row["n_obs"]),
                resid_counts.get(lid, 0),
                float(row["rms_arcsec"]),
                float(row["max_arcsec"]),
            )
        )
    if not good_rows:
        raise RuntimeError("No fit_ok and is_good linkage found.")

    if selection == "lowest-rms":
        good_rows.sort(key=lambda item: (item[3], -item[1], item[0]))
        return good_rows[0][0]

    if selection == "max-obs":
        good_rows.sort(key=lambda item: (-item[1], item[3], item[0]))
        return good_rows[0][0]

    good_rows.sort(key=lambda item: (-item[1], item[3], item[0]))
    return good_rows[0][0]


def format_tracklets(raw: str) -> str:
    text = str(raw).strip()
    if not text or text == "--":
        return "all tracklets kept"
    return text.replace(";", "\n")


def summarize_tracklet_selection(link_row) -> str:
    rejected = str(link_row["rejected_tracklets"]).strip()
    kept = str(link_row["inlier_tracklets"]).strip()
    if rejected in {"", "--"}:
        n_tracklets = int(link_row["n_tracklets"])
        return f"All {n_tracklets} tracklets retained in the final fit."
    if kept in {"", "--"}:
        return "Tracklet filtering occurred; inspect linkage table for details."
    return f"Inlier tracklets: {kept.replace(';', ', ')}. Rejected: {rejected.replace(';', ', ')}."


def main() -> None:
    args = parse_args()

    rr_dir = Path(args.root).expanduser() / args.night / "rr_links"
    orbit_dir = rr_dir / "orbit_confirm"
    links = Table.read(orbit_dir / "orbit_links.fits")
    residuals = Table.read(orbit_dir / "orbit_obs_residuals.fits")

    linkage_id = args.linkage_id
    if linkage_id is None:
        linkage_id = choose_linkage(links, residuals, selection=args.selection)

    link_rows = links[links["linkage_id"] == linkage_id]
    if len(link_rows) != 1:
        raise RuntimeError(f"Expected one orbit_links row for linkage_id={linkage_id}, got {len(link_rows)}.")
    link = link_rows[0]

    rows = residuals[residuals["linkage_id"] == linkage_id]
    if len(rows) == 0:
        raise RuntimeError(f"No residual rows found for linkage_id={linkage_id}.")

    order = np.argsort(np.asarray(rows["mjd"], dtype=float))
    rows = rows[order]

    mjd = np.asarray(rows["mjd"], dtype=float)
    dt_min = (mjd - mjd[0]) * 24.0 * 60.0
    ra = np.asarray(rows["ra_deg"], dtype=float)
    dec = np.asarray(rows["dec_deg"], dtype=float)
    resid = np.asarray(rows["resid_arcsec"], dtype=float)
    used = np.asarray(rows["used"], dtype=bool)

    x_arcsec = (ra - np.mean(ra)) * 3600.0 * np.cos(np.deg2rad(np.mean(dec)))
    y_arcsec = (dec - np.mean(dec)) * 3600.0

    outdir = Path(args.outdir).expanduser() if args.outdir else orbit_dir / "diagnostics"
    outdir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "figure.facecolor": "#f6f1e8",
            "axes.facecolor": "#fffaf1",
            "axes.edgecolor": "#5b5247",
            "axes.labelcolor": "#352f29",
            "xtick.color": "#352f29",
            "ytick.color": "#352f29",
            "text.color": "#231f1c",
            "font.size": 12,
            "axes.titlesize": 15,
            "axes.labelsize": 12,
        }
    )

    fig = plt.figure(figsize=(15, 8.8))
    gs = fig.add_gridspec(
        2,
        3,
        width_ratios=[1.55, 1.1, 0.9],
        height_ratios=[1.0, 1.0],
        left=0.05,
        right=0.97,
        top=0.9,
        bottom=0.08,
        wspace=0.24,
        hspace=0.24,
    )
    ax_resid = fig.add_subplot(gs[:, 0])
    ax_sky = fig.add_subplot(gs[:, 1])
    ax_stats = fig.add_subplot(gs[:, 2])

    fig.text(0.05, 0.955, "Orbit Confirm Diagnostic", fontsize=24, fontweight="bold", color="#1f3b4d")
    fig.text(0.05, 0.922, f"Night {args.night}   Linkage {linkage_id}", fontsize=13, color="#6b6257")

    ax_resid.axhspan(0.0, 0.10, color="#d8efe3", alpha=0.75, zorder=0)
    ax_resid.axhspan(0.10, max(0.25, np.nanmax(resid) * 1.2), color="#fff0d9", alpha=0.7, zorder=0)
    ax_resid.plot(dt_min, resid, color="#1f3b4d", lw=2.2, alpha=0.85, zorder=2)
    ax_resid.scatter(
        dt_min[used],
        resid[used],
        s=90,
        color="#198754",
        edgecolor="white",
        linewidth=1.4,
        zorder=3,
        label="used in final fit",
    )
    if np.any(~used):
        ax_resid.scatter(
            dt_min[~used],
            resid[~used],
            s=100,
            marker="X",
            color="#c44536",
            edgecolor="white",
            linewidth=1.0,
            zorder=4,
            label="rejected",
        )
    for x0, y0, key in zip(dt_min, resid, rows["obs_key"]):
        ax_resid.text(x0, y0 + 0.008, str(key).split(":")[0].replace("_cat.fits.gz", ""), fontsize=8.5, color="#655c52")
    ax_resid.set_title("Residual Versus Time")
    ax_resid.set_xlabel("Minutes Since First Detection")
    ax_resid.set_ylabel("Angular Residual (arcsec)")
    ax_resid.grid(True, color="#d7cdbf", alpha=0.7, lw=0.8)
    ax_resid.legend(loc="upper right", frameon=False)

    ax_sky.plot(x_arcsec, y_arcsec, color="#8a8175", lw=1.6, alpha=0.9, zorder=1)
    ax_sky.scatter(
        x_arcsec,
        y_arcsec,
        c=dt_min,
        cmap="cividis",
        s=120 + 700 * np.clip(resid, 0, None),
        edgecolor="white",
        linewidth=1.2,
        zorder=3,
    )
    for idx, (x0, y0) in enumerate(zip(x_arcsec, y_arcsec), start=1):
        ax_sky.text(x0, y0, str(idx), ha="center", va="center", fontsize=8.5, color="white", fontweight="bold")
    ax_sky.set_title("Observed Sky-Plane Track")
    ax_sky.set_xlabel(r"$\Delta$RA cos(Dec) (arcsec)")
    ax_sky.set_ylabel(r"$\Delta$Dec (arcsec)")
    ax_sky.grid(True, color="#d7cdbf", alpha=0.6, lw=0.8)
    ax_sky.set_aspect("equal", adjustable="datalim")
    cbar = fig.colorbar(ax_sky.collections[0], ax=ax_sky, fraction=0.046, pad=0.04)
    cbar.set_label("Minutes Since First Detection")

    ax_stats.axis("off")
    stats_lines = [
        ("fit status", "good" if bool(link["is_good"]) else "fit_ok only"),
        ("observations", f"{int(link['n_obs'])}"),
        ("tracklets", f"{int(link['n_tracklets'])}"),
        ("RMS", f"{float(link['rms_arcsec']):.3f} arcsec"),
        ("max residual", f"{float(link['max_arcsec']):.3f} arcsec"),
        ("semi-major axis", f"{float(link['a_au']):.2f} AU"),
        ("eccentricity", f"{float(link['ecc']):.3f}"),
        ("inclination", f"{float(link['inc_deg']):.2f} deg"),
    ]
    ax_stats.text(
        0.0,
        1.0,
        "Fit Summary",
        fontsize=17,
        fontweight="bold",
        color="#1f3b4d",
        va="top",
        transform=ax_stats.transAxes,
    )
    y = 0.90
    for label, value in stats_lines:
        ax_stats.text(0.0, y, label.upper(), fontsize=9, color="#8f8274", transform=ax_stats.transAxes)
        ax_stats.text(0.0, y - 0.06, value, fontsize=14, color="#231f1c", transform=ax_stats.transAxes)
        y -= 0.12

    ax_stats.text(
        0.0,
        y - 0.02,
        "Tracklet Selection",
        fontsize=17,
        fontweight="bold",
        color="#1f3b4d",
        transform=ax_stats.transAxes,
    )
    ax_stats.text(
        0.0,
        y - 0.18,
        summarize_tracklet_selection(link),
        fontsize=12,
        color="#231f1c",
        wrap=True,
        transform=ax_stats.transAxes,
    )

    png_path = outdir / f"{args.night}_link{linkage_id}_orbit_fit_diagnostic.png"
    fig.savefig(png_path, dpi=220)
    plt.close(fig)

    print(f"[write] {png_path}")


if __name__ == "__main__":
    main()
