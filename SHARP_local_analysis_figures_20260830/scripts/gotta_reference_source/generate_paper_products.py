#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, get_sun
from astropy.time import Time
from matplotlib.colors import LogNorm
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

ADOPTED_APERTURE_INDEX = 4
ADOPTED_MAG = f"Mag_Aper{ADOPTED_APERTURE_INDEX}"
ADOPTED_MAGERR = f"MagErr_Aper{ADOPTED_APERTURE_INDEX}"
ADOPTED_FLUX = f"Flux_Aper{ADOPTED_APERTURE_INDEX}"
ADOPTED_FLUXERR = f"FluxErr_Aper{ADOPTED_APERTURE_INDEX}"
ADOPTED_MAG_LABEL = r"$g_{\rm aper}$"
ADOPTED_MAGERR_LABEL = r"$\sigma(g_{\rm aper})$"
FOUR_PANEL_FIGSIZE = (18, 14)
MEDIAN_LEGEND_FONTSIZE = 21


def as_text(values) -> pd.Series:
    return pd.Series(values).astype(str).str.strip()


def load_table(path: Path) -> pd.DataFrame:
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        columns = {}
        for name in data.names:
            arr = np.array(data[name])
            if arr.dtype.byteorder in {">", "<"} and arr.dtype.byteorder != "=":
                arr = arr.astype(arr.dtype.newbyteorder("="), copy=False)
            columns[name] = arr
        return pd.DataFrame(columns)


def parse_night(source_file: str) -> str:
    match = re.search(r"(20\d{6})", str(source_file))
    if not match:
        return "unknown"
    text = match.group(1)
    return f"{text[:4]}-{text[4:6]}-{text[6:8]}"


def parse_field(source_file: str) -> str:
    match = re.match(r"(.+?)_20\d{6}_", str(source_file))
    if match:
        return match.group(1)
    return str(source_file).split("_")[0]


def truthy(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes", "y"})


def add_derived(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in [
        "ra",
        "dec",
        "RA_Win",
        "DEC_Win",
        "Mag_Kron",
        "MagErr_Kron",
        "Flux_Kron",
        "FluxErr_Kron",
        "mag",
        "epoch",
        "ang_rate_arcsec_hour",
        "ang_rate_deg_day",
        "phase_deg",
        "r_AU",
        "delta_AU",
        ADOPTED_MAG,
        ADOPTED_MAGERR,
        ADOPTED_FLUX,
        ADOPTED_FLUXERR,
    ]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    for idx in range(1, 13):
        for prefix in ["Mag_Aper", "MagErr_Aper", "Flux_Aper", "FluxErr_Aper"]:
            col = f"{prefix}{idx}"
            if col in out:
                out[col] = pd.to_numeric(out[col], errors="coerce")
        mag_col = f"Mag_Aper{idx}"
        err_col = f"MagErr_Aper{idx}"
        if mag_col in out:
            out.loc[(out[mag_col] <= 0) | (out[mag_col] > 40), mag_col] = np.nan
        if err_col in out:
            out.loc[(out[err_col] < 0) | (out[err_col] > 10), err_col] = np.nan
    for mag_col in ["Mag_Kron", "mag"]:
        if mag_col in out:
            out.loc[(out[mag_col] <= 0) | (out[mag_col] > 40), mag_col] = np.nan

    dra_deg = ((out["RA_Win"] - out["ra"] + 180.0) % 360.0) - 180.0
    out["dra_cosdec_arcsec"] = dra_deg * np.cos(np.deg2rad(out["dec"])) * 3600.0
    out["ddec_arcsec"] = (out["DEC_Win"] - out["dec"]) * 3600.0
    out["sep_arcsec"] = np.hypot(out["dra_cosdec_arcsec"], out["ddec_arcsec"])
    out["night"] = as_text(out["source_file"]).map(parse_night)
    out["field_id"] = as_text(out["source_file"]).map(parse_field)
    out["query_id"] = as_text(out["query_id"])
    out["name"] = as_text(out["name"])
    out["object_orbit_class_code"] = as_text(out["object_orbit_class_code"])
    out["object_orbit_class_name"] = as_text(out["object_orbit_class_name"])
    out["object_neo_bool"] = truthy(out["object_neo"])
    out["object_pha_bool"] = truthy(out["object_pha"])
    out["snr_kron_proxy"] = out["Flux_Kron"] / out["FluxErr_Kron"].replace(0, np.nan)
    out.loc[out["snr_kron_proxy"] <= 0, "snr_kron_proxy"] = np.nan
    out["snr_aper_proxy"] = out[ADOPTED_FLUX] / out[ADOPTED_FLUXERR].replace(0, np.nan)
    out.loc[out["snr_aper_proxy"] <= 0, "snr_aper_proxy"] = np.nan
    out["dmag_obs_minus_pred"] = out[ADOPTED_MAG] - out["mag"]
    out["elongation_deg"] = solar_elongation_degrees(out["ra"], out["dec"], out["epoch"])
    return out


def solar_elongation_degrees(ra: pd.Series, dec: pd.Series, epoch_mjd: pd.Series) -> np.ndarray:
    ra_values = pd.to_numeric(ra, errors="coerce").to_numpy(dtype=float)
    dec_values = pd.to_numeric(dec, errors="coerce").to_numpy(dtype=float)
    epoch_values = pd.to_numeric(epoch_mjd, errors="coerce").to_numpy(dtype=float)
    elongation = np.full(ra_values.shape, np.nan, dtype=float)
    mask = np.isfinite(ra_values) & np.isfinite(dec_values) & np.isfinite(epoch_values)
    if not np.any(mask):
        return elongation
    targets = SkyCoord(ra=ra_values[mask] * u.deg, dec=dec_values[mask] * u.deg, frame="icrs")
    sun = get_sun(Time(epoch_values[mask], format="mjd", scale="utc")).transform_to("icrs")
    elongation[mask] = targets.separation(sun).deg
    return elongation


def finite(series: pd.Series) -> np.ndarray:
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    return values[np.isfinite(values)]


def q(series: pd.Series, percentile: float) -> float:
    values = finite(series)
    if values.size == 0:
        return np.nan
    return float(np.nanpercentile(values, percentile))


def rms(series: pd.Series) -> float:
    values = finite(series)
    if values.size == 0:
        return np.nan
    return float(np.sqrt(np.nanmean(values**2)))


def numeric_from_text(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.extract(r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)", expand=False)
    return pd.to_numeric(text, errors="coerce")


def fmt(value: float, ndigits: int = 3) -> str:
    if pd.isna(value):
        return "--"
    if abs(value) >= 1000 and float(value).is_integer():
        return f"{value:,.0f}"
    if abs(value) >= 1000:
        return f"{value:,.1f}"
    return f"{value:.{ndigits}f}"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.size": 30,
            "axes.linewidth": 1.4,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "legend.frameon": False,
            "figure.dpi": 150,
        }
    )


def save_figure(fig: plt.Figure, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def histogram(ax, values, bins, xlabel, ylabel="Detections", logy=False, color="#4c78a8") -> None:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    ax.hist(values, bins=bins, color=color, histtype="stepfilled", alpha=0.5, edgecolor="#2b2b2b", linewidth=0.55)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
    ax.set_axisbelow(True)
    ax.grid(alpha=0.16, linewidth=0.6)
    ax.tick_params(labelsize=22)


def add_colorbar(ax, mappable, label: str) -> None:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3.5%", pad=0.06)
    cb = ax.figure.colorbar(mappable, cax=cax)
    cb.set_label(label, fontsize=20)
    cb.ax.tick_params(labelsize=17, width=0.9, length=4)
    cb.outline.set_linewidth(0.8)


def density_map(
    ax,
    x,
    y,
    xlabel,
    ylabel,
    xbins=60,
    ybins=42,
    xscale=None,
    yscale=None,
    cmap="rainbow",
    colorbar=True,
    xlim=None,
    ylim=None,
):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if xscale == "log":
        mask &= x > 0
    if yscale == "log":
        mask &= y > 0
    if xlim is not None:
        mask &= (x >= xlim[0]) & (x <= xlim[1])
    if ylim is not None:
        mask &= (y >= ylim[0]) & (y <= ylim[1])
    x_plot = x[mask]
    y_plot = y[mask]
    if x_plot.size == 0:
        return

    x_work = np.log10(x_plot) if xscale == "log" else x_plot
    y_work = np.log10(y_plot) if yscale == "log" else y_plot
    x_edges = np.linspace(np.nanmin(x_work), np.nanmax(x_work), xbins + 1)
    y_edges = np.linspace(np.nanmin(y_work), np.nanmax(y_work), ybins + 1)
    counts, x_edges, y_edges = np.histogram2d(x_work, y_work, bins=[x_edges, y_edges])
    counts = counts.T
    masked_counts = np.ma.masked_where(counts <= 0, counts)
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad((1, 1, 1, 0))
    norm = LogNorm(vmin=1, vmax=max(1, float(np.nanmax(counts))))
    if xscale is None and yscale is None:
        mesh = ax.imshow(
            masked_counts,
            extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
            origin="lower",
            aspect="auto",
            cmap=cmap_obj,
            norm=norm,
            interpolation="bilinear",
            rasterized=True,
        )
    else:
        x_draw = 10**x_edges if xscale == "log" else x_edges
        y_draw = 10**y_edges if yscale == "log" else y_edges
        mesh = ax.pcolormesh(
            x_draw,
            y_draw,
            masked_counts,
            cmap=cmap_obj,
            norm=norm,
            shading="auto",
            rasterized=True,
        )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xscale:
        ax.set_xscale(xscale)
    if yscale:
        ax.set_yscale(yscale)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.tick_params(labelsize=22)
    if colorbar:
        add_colorbar(ax, mesh, "Detections")


def density_colored_scatter(
    ax,
    x,
    y,
    xlabel,
    ylabel,
    xbins=240,
    ybins=180,
    xlim=None,
    ylim=None,
    xscale=None,
    yscale=None,
    cmap="viridis",
    point_size=3.0,
    alpha=0.82,
    smooth_sigma=2.2,
):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if xscale == "log":
        mask &= x > 0
    if yscale == "log":
        mask &= y > 0
    if xlim is not None:
        mask &= (x >= xlim[0]) & (x <= xlim[1])
    if ylim is not None:
        mask &= (y >= ylim[0]) & (y <= ylim[1])
    x_plot = x[mask]
    y_plot = y[mask]
    if x_plot.size == 0:
        return

    x_work = np.log10(x_plot) if xscale == "log" else x_plot
    y_work = np.log10(y_plot) if yscale == "log" else y_plot
    x_range = np.log10(xlim) if xscale == "log" and xlim is not None else (np.nanmin(x_work), np.nanmax(x_work))
    y_range = np.log10(ylim) if yscale == "log" and ylim is not None else (np.nanmin(y_work), np.nanmax(y_work))
    counts, x_edges, y_edges = np.histogram2d(x_work, y_work, bins=[xbins, ybins], range=[x_range, y_range])
    smoothed = gaussian_filter(counts.astype(float), sigma=smooth_sigma, mode="nearest")
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    interpolator = RegularGridInterpolator(
        (x_centers, y_centers),
        smoothed,
        bounds_error=False,
        fill_value=1.0,
    )
    density = interpolator(np.column_stack([x_work, y_work]))
    density = np.clip(density, 1.0, None)
    order = np.argsort(density)
    sc = ax.scatter(
        x_plot[order],
        y_plot[order],
        c=density[order],
        s=point_size,
        cmap=cmap,
        norm=LogNorm(vmin=1, vmax=max(1, float(np.nanmax(density)))),
        marker="o",
        linewidths=0,
        alpha=alpha,
        rasterized=True,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xscale:
        ax.set_xscale(xscale)
    if yscale:
        ax.set_yscale(yscale)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_axisbelow(True)
    ax.grid(alpha=0.12, linewidth=0.55)
    ax.tick_params(labelsize=22)
    add_colorbar(ax, sc, "Local density")


def running_rate_statistics(rate: pd.Series, sep: pd.Series, nbins: int = 20) -> pd.DataFrame:
    data = pd.DataFrame({"rate": pd.to_numeric(rate, errors="coerce"), "sep": pd.to_numeric(sep, errors="coerce")})
    data = data[np.isfinite(data["rate"]) & np.isfinite(data["sep"]) & (data["rate"] > 0)]
    if data.empty:
        return pd.DataFrame(columns=["rate", "p16", "median", "p84", "count"])
    edges = np.logspace(np.log10(data["rate"].quantile(0.002)), np.log10(data["rate"].quantile(0.998)), nbins + 1)
    rows = []
    for left, right in zip(edges[:-1], edges[1:]):
        sub = data[(data["rate"] >= left) & (data["rate"] < right)]["sep"]
        if len(sub) < 20:
            continue
        rows.append(
            {
                "rate": np.sqrt(left * right),
                "p16": np.nanpercentile(sub, 16),
                "median": np.nanpercentile(sub, 50),
                "p84": np.nanpercentile(sub, 84),
                "std": np.nanstd(sub),
                "count": len(sub),
            }
        )
    return pd.DataFrame(rows)


def running_linear_statistics(x: pd.Series, y: pd.Series, edges: np.ndarray, min_count: int = 40) -> pd.DataFrame:
    data = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")})
    data = data[np.isfinite(data["x"]) & np.isfinite(data["y"])]
    rows = []
    for left, right in zip(edges[:-1], edges[1:]):
        sub = data[(data["x"] >= left) & (data["x"] < right)]["y"]
        if len(sub) < min_count:
            continue
        rows.append(
            {
                "x": 0.5 * (left + right),
                "p16": np.nanpercentile(sub, 16),
                "median": np.nanpercentile(sub, 50),
                "p84": np.nanpercentile(sub, 84),
                "std": np.nanstd(sub),
                "count": len(sub),
            }
        )
    return pd.DataFrame(rows, columns=["x", "p16", "median", "p84", "std", "count"])


def write_csv_and_latex(
    df: pd.DataFrame,
    csv_path: Path,
    caption: str,
    label: str,
    column_format: str | None = None,
    font_size: str = r"\footnotesize",
    tabcolsep_pt: float = 6,
    latex_columns: list[str] | None = None,
    center_over_textwidth: bool = False,
    tabular_width: str | None = None,
) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    tex_path = csv_path.with_suffix(".tex")
    if column_format is None:
        column_format = "l" * len(df.columns)
    with tex_path.open("w", encoding="utf-8") as handle:
        handle.write("\\begin{table*}\n\\centering\n")
        handle.write(f"\\caption{{{caption}}}\n")
        handle.write(f"\\label{{{label}}}\n")
        handle.write(f"{font_size}\n")
        handle.write(f"\\setlength{{\\tabcolsep}}{{{tabcolsep_pt:g}pt}}\n")
        handle.write("\\renewcommand{\\arraystretch}{1.12}\n")
        if center_over_textwidth:
            handle.write("\\makebox[\\textwidth][c]{%\n")
        if tabular_width is None:
            handle.write(f"\\begin{{tabular}}{{{column_format}}}\n")
        else:
            handle.write(f"\\begin{{tabular*}}{{{tabular_width}}}{{{column_format}}}\n")
        handle.write("\\toprule\n")
        header = latex_columns if latex_columns is not None else list(df.columns)
        if latex_columns is None:
            handle.write(" & ".join(latex_escape(col) for col in header) + " \\\\\n")
        else:
            handle.write(" & ".join(header) + " \\\\\n")
        handle.write("\\midrule\n")
        for _, row in df.iterrows():
            handle.write(" & ".join(latex_escape(x) for x in row.to_list()) + " \\\\\n")
        if tabular_width is None:
            handle.write("\\bottomrule\n\\end{tabular}")
        else:
            handle.write("\\bottomrule\n\\end{tabular*}")
        if center_over_textwidth:
            handle.write("%\n}\n")
        else:
            handle.write("\n")
        handle.write("\\end{table*}\n")


def write_overall_statistics_table(df: pd.DataFrame, tex_path: Path, caption: str, label: str) -> None:
    with tex_path.open("w", encoding="utf-8") as handle:
        handle.write("\\begin{table*}\n\\centering\n")
        handle.write(f"\\caption{{{caption}}}\n")
        handle.write(f"\\label{{{label}}}\n")
        handle.write("\\footnotesize\n")
        handle.write("\\setlength{\\tabcolsep}{5pt}\n")
        handle.write("\\renewcommand{\\arraystretch}{1.12}\n")
        handle.write("\\begin{tabularx}{0.96\\textwidth}{@{}>{\\raggedright\\arraybackslash}p{0.42\\textwidth}>{\\raggedleft\\arraybackslash}p{0.22\\textwidth}>{\\raggedright\\arraybackslash}X@{}}\n")
        handle.write("\\toprule\n")
        handle.write("Quantity & Value & Note \\\\\n")
        handle.write("\\midrule\n")
        for _, row in df.iterrows():
            handle.write(" & ".join(latex_escape(x) for x in row.to_list()) + " \\\\\n")
        handle.write("\\bottomrule\n\\end{tabularx}\n\\end{table*}\n")


def latex_escape(value) -> str:
    text = str(value)
    if "$" in text or "\\" in text:
        return text.replace("%", r"\%")
    return text.replace("_", r"\_").replace("%", r"\%")


def make_tables(df: pd.DataFrame, outdir: Path) -> dict[str, pd.DataFrame]:
    overall = pd.DataFrame(
        [
            ("Accepted source--ephemeris matches", f"{len(df):,}", "Recovered detections"),
            ("Distinct known asteroids", f"{df['query_id'].nunique():,}", "Unique minor-planet identifiers"),
            ("Exposure catalogs with accepted matches", f"{df['source_file'].nunique():,}", "Per-exposure catalog products"),
            ("UTC nights represented", f"{df['night'].nunique():,}", "Recovered sample date coverage"),
            ("Field groups represented", f"{df['field_id'].nunique():,}", "Grouped by pointing metadata"),
            ("NEO-class objects", f"{df.loc[df['object_neo_bool'], 'query_id'].nunique():,}", "Identified from small-body metadata"),
            ("PHA-class objects", f"{df.loc[df['object_pha_bool'], 'query_id'].nunique():,}", "Identified from small-body metadata"),
            ("Median ephemeris magnitude", fmt(q(df["mag"], 50), 4), "Predicted brightness"),
            (r"Median $g_{\rm aper}$", fmt(q(df[ADOPTED_MAG], 50), 4), "Adopted aperture magnitude"),
            (r"Median $\sigma(g_{\rm aper})$", fmt(q(df[ADOPTED_MAGERR], 50), 4) + " mag", "Adopted aperture uncertainty"),
            ("Median aperture S/N proxy", fmt(q(df["snr_aper_proxy"], 50), 4), "Flux/error in adopted aperture"),
            ("Median observed-minus-predicted separation", fmt(q(df["sep_arcsec"], 50), 3) + " arcsec", "Recovered detections"),
            ("84th-percentile observed-minus-predicted separation", fmt(q(df["sep_arcsec"], 84), 3) + " arcsec", "Recovered detections"),
            ("2D residual RMS", fmt(rms(df["sep_arcsec"]), 3) + " arcsec", "RMS of the 2D separation"),
            (r"RMS in $\Delta\alpha\cos\delta$", fmt(rms(df["dra_cosdec_arcsec"]), 3) + " arcsec", "Coordinate residual"),
            (r"RMS in $\Delta\delta$", fmt(rms(df["ddec_arcsec"]), 3) + " arcsec", "Coordinate residual"),
            ("Median angular rate", r"\mbox{" + fmt(q(df["ang_rate_arcsec_hour"], 50), 3) + r" arcsec h$^{-1}$}", "Derived from augmented geometry"),
            ("Median phase angle", fmt(q(df["phase_deg"], 50), 3) + " deg", "Derived from augmented geometry"),
            ("Associations without geometry augmentation", f"{int(truthy(df['horizons_failed']).sum()):,}", "Excluded from geometry statistics"),
        ],
        columns=["Quantity", "Value", "Note"],
    )

    orbit_work = df.copy()
    sparse_codes = {"AST", "CEN"}
    sparse_mask = orbit_work["object_orbit_class_code"].isin(sparse_codes)
    orbit_work.loc[sparse_mask, "object_orbit_class_code"] = "Other"
    orbit_work.loc[sparse_mask, "object_orbit_class_name"] = "Other / unclassified"
    orbit = (
        orbit_work.groupby(["object_orbit_class_code", "object_orbit_class_name"], dropna=False)
        .agg(
            Detections=("query_id", "size"),
            **{
                "Unique objects": ("query_id", "nunique"),
                r"Median $g_{\rm aper}$": (ADOPTED_MAG, "median"),
                "Median separation": ("sep_arcsec", "median"),
            },
        )
        .reset_index()
        .sort_values("Detections", ascending=False)
    )
    orbit.columns = ["Code", "Orbit class", "Detections", "Unique objects", r"Median $g_{\rm aper}$", "Median separation"]
    total_detections = orbit["Detections"].sum()
    total_objects = orbit["Unique objects"].sum()
    orbit.insert(3, "Detection fraction (%)", orbit["Detections"] / total_detections * 100.0)
    orbit.insert(5, "Object fraction (%)", orbit["Unique objects"] / total_objects * 100.0)
    orbit["Detections"] = orbit["Detections"].map(lambda x: f"{x:,}")
    orbit["Unique objects"] = orbit["Unique objects"].map(lambda x: f"{x:,}")
    orbit["Detection fraction (%)"] = orbit["Detection fraction (%)"].map(lambda x: fmt(x, 2))
    orbit["Object fraction (%)"] = orbit["Object fraction (%)"].map(lambda x: fmt(x, 2))
    orbit[r"Median $g_{\rm aper}$"] = orbit[r"Median $g_{\rm aper}$"].map(lambda x: fmt(x, 3))
    orbit["Median separation"] = orbit["Median separation"].map(lambda x: fmt(x, 3))
    orbit = orbit.rename(columns={r"Median $g_{\rm aper}$": r"Median $g_{\rm aper}$ (mag)", "Median separation": "Median separation (arcsec)"})

    work = df.copy()
    mag_bins = pd.IntervalIndex.from_tuples([(12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19), (19, 20)], closed="left")
    work["mag_bin"] = pd.cut(work[ADOPTED_MAG], mag_bins)
    mag_ast = (
        work.groupby("mag_bin", observed=True)
        .agg(
            Detections=("query_id", "size"),
            **{
                "Median separation": ("sep_arcsec", "median"),
                "RMS separation": ("sep_arcsec", rms),
                "Median mag. error": (ADOPTED_MAGERR, "median"),
            },
        )
        .reset_index()
    )
    mag_ast["Magnitude bin"] = mag_ast["mag_bin"].map(lambda x: f"$[{x.left:.0f},{x.right:.0f})$")
    mag_ast = mag_ast[["Magnitude bin", "Detections", "Median separation", "RMS separation", "Median mag. error"]]
    for col in ["Median separation", "RMS separation", "Median mag. error"]:
        mag_ast[col] = mag_ast[col].map(lambda x: fmt(x, 3))
    mag_ast["Detections"] = mag_ast["Detections"].map(lambda x: f"{x:,}")

    rate_bins = pd.IntervalIndex.from_tuples([(0, 10), (10, 20), (20, 40), (40, 80), (80, 160), (160, 320), (320, 1000)], closed="left")
    work["rate_bin"] = pd.cut(work["ang_rate_arcsec_hour"], rate_bins)
    rate_ast = (
        work.groupby("rate_bin", observed=True)
        .agg(
            Detections=("query_id", "size"),
            **{
                "Median separation": ("sep_arcsec", "median"),
                "RMS separation": ("sep_arcsec", rms),
                r"Median $g_{\rm aper}$": (ADOPTED_MAG, "median"),
            },
        )
        .reset_index()
    )
    rate_ast["Rate bin"] = rate_ast["rate_bin"].map(lambda x: f"$[{x.left:.0f},{x.right:.0f})$")
    rate_ast = rate_ast[["Rate bin", "Detections", "Median separation", "RMS separation", r"Median $g_{\rm aper}$"]]
    for col in ["Median separation", "RMS separation", r"Median $g_{\rm aper}$"]:
        rate_ast[col] = rate_ast[col].map(lambda x: fmt(x, 3))
    rate_ast["Detections"] = rate_ast["Detections"].map(lambda x: f"{x:,}")

    nightly = (
        df.groupby("night")
        .agg(
            Detections=("query_id", "size"),
            **{
                "Unique objects": ("query_id", "nunique"),
                "Exposure catalogs": ("source_file", "nunique"),
                r"Median $g_{\rm aper}$": (ADOPTED_MAG, "median"),
                "Median separation": ("sep_arcsec", "median"),
            },
        )
        .reset_index()
        .sort_values("Detections", ascending=False)
        .head(5)
    )
    nightly.columns = ["UTC date", "Detections", "Unique objects", "Exposure catalogs", r"Median $g_{\rm aper}$", "Median separation"]
    for col in ["Detections", "Unique objects", "Exposure catalogs"]:
        nightly[col] = nightly[col].map(lambda x: f"{x:,}")
    for col in [r"Median $g_{\rm aper}$", "Median separation"]:
        nightly[col] = nightly[col].map(lambda x: fmt(x, 3))
    nightly = nightly.rename(columns={r"Median $g_{\rm aper}$": r"Median $g_{\rm aper}$ (mag)", "Median separation": "Median separation (arcsec)"})

    obj = (
        df.groupby("query_id")
        .agg(
            Name=("name", "first"),
            Detections=("query_id", "size"),
            **{
                "Nights": ("night", "nunique"),
                "First epoch": ("epoch", "min"),
                "Last epoch": ("epoch", "max"),
                r"Median $g_{\rm aper}$": (ADOPTED_MAG, "median"),
                "Class": ("object_orbit_class_code", "first"),
            },
        )
        .reset_index()
        .sort_values("Detections", ascending=False)
        .head(5)
    )
    obj["Baseline days"] = obj["Last epoch"] - obj["First epoch"]
    obj = obj[["query_id", "Name", "Detections", "Nights", "First epoch", "Last epoch", "Baseline days", r"Median $g_{\rm aper}$", "Class"]]
    obj.columns = ["Object ID", "Name", "Detections", "Nights", "First MJD", "Last MJD", "Baseline (days)", r"Median $g_{\rm aper}$", "Class"]
    for col in ["First MJD", "Last MJD"]:
        obj[col] = obj[col].map(lambda x: f"{x:.1f}")
    for col in ["Baseline (days)", r"Median $g_{\rm aper}$"]:
        obj[col] = obj[col].map(lambda x: f"{x:.3f}")
    for col in ["Detections", "Nights"]:
        obj[col] = obj[col].map(lambda x: f"{x:,}")

    outdir.mkdir(parents=True, exist_ok=True)
    overall.to_csv(outdir / "overall_statistics.csv", index=False)
    write_overall_statistics_table(overall, outdir / "overall_statistics.tex", "Statistics for the recovered GOTTA known-asteroid sample.", "tab:summary")
    write_csv_and_latex(
        orbit,
        outdir / "orbit_class_statistics.csv",
        "Orbit-class composition of the recovered known-asteroid sample. Separations are in arcsec. Classes with only a few detections are grouped into ``Other/unclassified''; this group includes objects without a more specific dynamical class in the adopted small-body metadata.",
        "tab:orbit_class",
        r"@{}ll@{\hspace{10pt}}r@{\extracolsep{\fill}}rrrrr@{}",
        font_size=r"\scriptsize",
        tabcolsep_pt=3.0,
        center_over_textwidth=True,
        tabular_width=r"1.10\textwidth",
        latex_columns=[
            "Code",
            "Orbit class",
            "Detections",
            r"\begin{tabular}[c]{@{}c@{}}Detection fraction\\(\%)\end{tabular}",
            "Unique objects",
            r"\begin{tabular}[c]{@{}c@{}}Object fraction\\(\%)\end{tabular}",
            r"\begin{tabular}[c]{@{}r@{}}Median $g_{\rm aper}$\\(mag)\end{tabular}",
            r"\begin{tabular}[c]{@{}r@{}}Median separation\\(arcsec)\end{tabular}",
        ],
    )
    write_csv_and_latex(mag_ast, outdir / "astrometry_by_magnitude.csv", r"Astrometric residuals as a function of adopted aperture magnitude $g_{\rm aper}$. Separations are in arcsec.", "tab:mag_astrometry", r"lrrrr")
    write_csv_and_latex(rate_ast, outdir / "astrometry_by_rate.csv", r"Astrometric residuals as a function of sky-plane angular rate. Rate bins are in arcsec h$^{-1}$ and separations are in arcsec.", "tab:rate_astrometry", r"lrrrr")
    write_csv_and_latex(
        nightly,
        outdir / "nightly_top5.csv",
        "Five UTC nights with the largest numbers of recovered known-asteroid detections.",
        "tab:nightly_top5",
        r"lrrrrr",
        latex_columns=[
            "UTC date",
            "Detections",
            "Unique objects",
            "Exposure catalogs",
            r"\begin{tabular}[c]{@{}r@{}}Median $g_{\rm aper}$\\(mag)\end{tabular}",
            r"\begin{tabular}[c]{@{}r@{}}Median separation\\(arcsec)\end{tabular}",
        ],
    )
    write_csv_and_latex(
        obj,
        outdir / "most_observed_objects.csv",
        "Five most frequently recovered known asteroids in the GOTTA Prototype sample.",
        "tab:most_observed_objects",
        r"llrrrrrrl",
        font_size=r"\scriptsize",
        tabcolsep_pt=3.5,
    )
    return {"overall": overall, "orbit": orbit, "mag_ast": mag_ast, "rate_ast": rate_ast, "nightly": nightly, "objects": obj}


def plot_photometry(df: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=FOUR_PANEL_FIGSIZE)
    histogram(axes[0, 0], finite(df[ADOPTED_MAG]), np.linspace(10, 21, 45), f"{ADOPTED_MAG_LABEL} [mag]", color="#4c78a8")
    med = q(df[ADOPTED_MAG], 50)
    axes[0, 0].axvline(med, color="black", linestyle="--", linewidth=1.6, label=f"median = {med:.2f}")
    axes[0, 0].legend(frameon=False, fontsize=MEDIAN_LEGEND_FONTSIZE, loc="upper left")
    axes[0, 0].set_xlim(10, 21)
    histogram(axes[0, 1], finite(df["snr_aper_proxy"]), np.logspace(-0.2, 3.2, 58), "Aperture S/N proxy", logy=True, color="#f28e2b")
    axes[0, 1].set_xscale("log")
    med = q(df["snr_aper_proxy"], 50)
    axes[0, 1].axvline(med, color="black", linestyle="--", linewidth=1.6, label=f"median = {med:.1f}")
    axes[0, 1].legend(frameon=False, fontsize=MEDIAN_LEGEND_FONTSIZE, loc="upper right")
    density_colored_scatter(
        axes[1, 0],
        df[ADOPTED_MAG],
        df[ADOPTED_MAGERR],
        f"{ADOPTED_MAG_LABEL} [mag]",
        f"{ADOPTED_MAGERR_LABEL} [mag]",
        xbins=260,
        ybins=190,
        xlim=(10, 21),
        ylim=(0.0, 0.6),
        point_size=3.0,
    )
    med = q(df[ADOPTED_MAGERR], 50)
    axes[1, 0].axhline(med, color="#d62728", linewidth=2.2, label=f"median = {med:.3f}")
    axes[1, 0].legend(frameon=False, fontsize=MEDIAN_LEGEND_FONTSIZE, loc="upper left")
    axes[1, 0].set_aspect("auto")
    density_colored_scatter(
        axes[1, 1],
        df[ADOPTED_MAG],
        df["dmag_obs_minus_pred"],
        f"{ADOPTED_MAG_LABEL} [mag]",
        r"$\Delta m$ [mag]",
        xbins=260,
        ybins=190,
        xlim=(10, 21),
        ylim=(-1.5, 1.5),
        point_size=3.0,
    )
    dmag_values = finite(df["dmag_obs_minus_pred"])
    med = float(np.nanmedian(dmag_values))
    sigma = float(np.nanstd(dmag_values))
    axes[1, 1].axhline(med, color="#d62728", linewidth=2.2, label=f"median = {med:.3f}")
    axes[1, 1].axhline(
        med + sigma,
        color="#d62728",
        linestyle="--",
        linewidth=1.4,
        alpha=0.85,
        label=r"median $\pm 1\sigma$",
    )
    axes[1, 1].axhline(med - sigma, color="#d62728", linestyle="--", linewidth=1.4, alpha=0.85)
    axes[1, 1].legend(frameon=False, fontsize=MEDIAN_LEGEND_FONTSIZE, loc="upper left")
    axes[1, 1].set_aspect("auto")
    fig.tight_layout()
    save_figure(fig, outdir / "photometric_statistics")


def plot_orbit_class_pies(df: pd.DataFrame, outdir: Path) -> None:
    objects = (
        df.sort_values("epoch")
        .groupby("query_id")
        .agg(
            orbit_class=("object_orbit_class_code", "first"),
            neo=("object_neo_bool", "max"),
            pha=("object_pha_bool", "max"),
        )
    )
    main_counts = pd.Series(
        {
            "MBA": int((objects["orbit_class"] == "MBA").sum()),
            "non-MBA": int((objects["orbit_class"] != "MBA").sum()),
        }
    )
    non_mba = objects[objects["orbit_class"] != "MBA"].copy()
    detailed_counts = pd.Series(
        {
            "OMB": int((non_mba["orbit_class"] == "OMB").sum()),
            "TJN": int((non_mba["orbit_class"] == "TJN").sum()),
            "IMB": int((non_mba["orbit_class"] == "IMB").sum()),
            "MCA": int((non_mba["orbit_class"] == "MCA").sum()),
            "NEO/PHA/Other": int(
                (~non_mba["orbit_class"].isin(["OMB", "TJN", "IMB", "MCA"])).sum()
            ),
        }
    )
    main_colors = ["#4c78a8", "#f28e2b"]
    inset_colors = ["#59a14f", "#e15759", "#76b7b2", "#edc948", "#b07aa1"]

    fig = plt.figure(figsize=(10.2, 6.2))
    ax = fig.add_axes([0.04, 0.12, 0.56, 0.78])
    inset = fig.add_axes([0.58, 0.18, 0.36, 0.64])
    wedgeprops = {"edgecolor": "#2b2b2b", "linewidth": 0.65}
    main_wedges, main_texts, main_autotexts = ax.pie(
        main_counts.to_numpy(),
        labels=main_counts.index,
        autopct="%1.1f%%",
        startangle=90,
        counterclock=False,
        colors=main_colors,
        wedgeprops=wedgeprops,
        textprops={"fontsize": 17, "fontfamily": "Times New Roman"},
        pctdistance=0.62,
        labeldistance=1.06,
    )
    for text in main_autotexts:
        text.set_fontsize(15)
    for patch in main_wedges:
        patch.set_alpha(0.5)
    ax.set_aspect("equal")

    fig.patches.append(
        FancyArrowPatch(
            (0.31, 0.80),
            (0.59, 0.62),
            transform=fig.transFigure,
            arrowstyle="-|>",
            connectionstyle="arc3,rad=-0.12",
            mutation_scale=18,
            linewidth=1.4,
            color="#555555",
            alpha=0.9,
            clip_on=False,
        )
    )

    inset_wedges, _, inset_autotexts = inset.pie(
        detailed_counts.to_numpy(),
        labels=None,
        autopct="%1.1f%%",
        startangle=90,
        counterclock=False,
        colors=inset_colors,
        wedgeprops=wedgeprops,
        textprops={"fontsize": 10.5, "fontfamily": "Times New Roman"},
        pctdistance=0.68,
    )
    for patch in inset_wedges:
        patch.set_alpha(0.5)
    inset.set_aspect("equal")
    inset.legend(
        inset_wedges,
        detailed_counts.index,
        frameon=False,
        fontsize=11,
        loc="center left",
        bbox_to_anchor=(0.92, 0.5),
        handlelength=1.0,
        handletextpad=0.45,
    )
    save_figure(fig, outdir / "orbit_class_composition_pies")


def plot_astrometry(df: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=FOUR_PANEL_FIGSIZE)
    histogram(axes[0, 0], finite(df["sep_arcsec"]), np.linspace(0, 1.5, 61), "Separation [arcsec]", logy=True)
    med = q(df["sep_arcsec"], 50)
    p84 = q(df["sep_arcsec"], 84)
    axes[0, 0].axvline(med, color="black", linestyle="--", linewidth=1.4, label=f"median={med:.3f}")
    axes[0, 0].axvline(p84, color="black", linestyle=":", linewidth=1.4, label=f"84%={p84:.3f}")
    axes[0, 0].legend(fontsize=MEDIAN_LEGEND_FONTSIZE)
    histogram(axes[0, 1], finite(df["dra_cosdec_arcsec"]), np.linspace(-1.2, 1.2, 61), r"$\Delta\alpha\cos\delta$ [arcsec]", logy=True, color="#59a14f")
    histogram(axes[1, 0], finite(df["ddec_arcsec"]), np.linspace(-1.2, 1.2, 61), r"$\Delta\delta$ [arcsec]", logy=True, color="#f28e2b")
    density_colored_scatter(
        axes[1, 1],
        df[ADOPTED_MAG],
        df["sep_arcsec"],
        f"{ADOPTED_MAG_LABEL} [mag]",
        "Separation [arcsec]",
        xbins=260,
        ybins=190,
        xlim=(10, 21),
        ylim=(0, 1.5),
        point_size=3.0,
    )
    mag_stats = running_linear_statistics(df[ADOPTED_MAG], df["sep_arcsec"], np.linspace(10, 21, 23), min_count=40)
    if not mag_stats.empty:
        axes[1, 1].plot(mag_stats["x"], mag_stats["median"], color="#d62728", linewidth=2.2, label="binned median")
        axes[1, 1].fill_between(
            mag_stats["x"],
            mag_stats["p16"],
            mag_stats["p84"],
            color="#d62728",
            alpha=0.42,
            label="16--84 percentile",
        )
        axes[1, 1].legend(fontsize=MEDIAN_LEGEND_FONTSIZE, loc="upper left")
    fig.tight_layout()
    save_figure(fig, outdir / "astrometric_residuals")

    fig, ax = plt.subplots(figsize=(10, 5.3))
    data = df[np.isfinite(df["ang_rate_arcsec_hour"]) & np.isfinite(df["sep_arcsec"]) & (df["ang_rate_arcsec_hour"] > 0)]
    main = data[(data["ang_rate_arcsec_hour"] >= 0.7) & (data["ang_rate_arcsec_hour"] <= 150) & (data["sep_arcsec"] <= 1.5)]
    density_colored_scatter(
        ax,
        main["ang_rate_arcsec_hour"],
        main["sep_arcsec"],
        r"Angular rate [arcsec h$^{-1}$]",
        "Separation [arcsec]",
        xbins=260,
        ybins=190,
        xscale="log",
        xlim=(0.7, 150),
        ylim=(0, 1.5),
        point_size=3.0,
    )
    stats = running_rate_statistics(main["ang_rate_arcsec_hour"], main["sep_arcsec"])
    if not stats.empty:
        ax.plot(stats["rate"], stats["median"], color="#d62728", linewidth=2.2, label="binned median")
        ax.fill_between(stats["rate"], stats["p16"], stats["p84"], color="#d62728", alpha=0.42, label="16--84 percentile")
        ax.legend(fontsize=18, loc="upper left")
    ax.set_xscale("log")
    ax.set_xlim(0.7, 150)
    ax.set_ylim(0, 1.5)
    ax.tick_params(labelsize=22)
    fig.tight_layout()
    save_figure(fig, outdir / "separation_vs_rate")


def plot_geometry(df: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=FOUR_PANEL_FIGSIZE)
    histogram(axes[0, 0], finite(df["ang_rate_arcsec_hour"]), np.logspace(-1, 3.2, 65), r"Angular rate [arcsec h$^{-1}$]", logy=True)
    axes[0, 0].set_xscale("log")
    med = q(df["ang_rate_arcsec_hour"], 50)
    axes[0, 0].axvline(med, color="black", linestyle="--", linewidth=1.6, label=rf"median = {med:.2f} arcsec h$^{{-1}}$")
    axes[0, 0].legend(frameon=False, fontsize=MEDIAN_LEGEND_FONTSIZE, loc="upper left")
    histogram(axes[0, 1], finite(df["elongation_deg"]), np.linspace(0, 180, 61), "Solar elongation [deg]", color="#59a14f")
    med = q(df["elongation_deg"], 50)
    axes[0, 1].axvline(med, color="black", linestyle="--", linewidth=1.6, label=f"median = {med:.2f} deg")
    axes[0, 1].legend(frameon=False, fontsize=MEDIAN_LEGEND_FONTSIZE, loc="upper left")
    histogram(axes[1, 0], finite(df["r_AU"]), np.linspace(0, 6, 61), "Heliocentric distance [AU]", color="#f28e2b")
    med = q(df["r_AU"], 50)
    axes[1, 0].axvline(med, color="black", linestyle="--", linewidth=1.6, label=f"median = {med:.2f} AU")
    axes[1, 0].legend(frameon=False, fontsize=MEDIAN_LEGEND_FONTSIZE, loc="upper right")
    histogram(axes[1, 1], finite(df["delta_AU"]), np.linspace(0, 5, 61), "Topocentric distance [AU]", color="#e15759")
    med = q(df["delta_AU"], 50)
    axes[1, 1].axvline(med, color="black", linestyle="--", linewidth=1.6, label=f"median = {med:.2f} AU")
    axes[1, 1].legend(frameon=False, fontsize=MEDIAN_LEGEND_FONTSIZE, loc="upper right")
    fig.tight_layout()
    save_figure(fig, outdir / "motion_geometry")


def plot_temporal(df: pd.DataFrame, outdir: Path) -> None:
    work = df.copy()
    work["orbit_a_AU_numeric"] = numeric_from_text(work["orbit_elements_a"])
    work["orbit_e_numeric"] = numeric_from_text(work["orbit_elements_e"])
    objects = (
        work.groupby("query_id")
        .agg(
            detections=("query_id", "size"),
            a_AU=("orbit_a_AU_numeric", "first"),
            eccentricity=("orbit_e_numeric", "first"),
        )
        .reset_index()
    )
    objects = objects[np.isfinite(objects["a_AU"]) & np.isfinite(objects["eccentricity"]) & (objects["a_AU"] > 0)]
    objects = objects.sort_values("detections")
    max_detections = max(1.0, float(objects["detections"].max()))
    marker_sizes = np.interp(
        np.log10(objects["detections"].to_numpy(dtype=float)),
        [0.0, np.log10(max_detections)],
        [4.0, 115.0],
    )

    fig, ax = plt.subplots(figsize=(9, 7.2))
    sc = ax.scatter(
        objects["a_AU"],
        objects["eccentricity"],
        c=objects["detections"],
        s=marker_sizes,
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max_detections),
        alpha=0.86,
        linewidths=0,
        rasterized=True,
    )
    a_line = np.logspace(np.log10(1.3), 2, 100)
    e_line = 1 - 1.3 / a_line
    ax.plot(a_line, e_line, "k--", linewidth=2.0, label=r"$q=1.3\,\mathrm{AU}$")
    ax.legend(fontsize=18, loc="lower right", frameon=False, handlelength=2.4)
    ax.set_xscale("log")
    ax.set_xlabel("Semimajor axis [AU]")
    ax.set_ylabel("Eccentricity")
    ax.set_axisbelow(True)
    ax.grid(alpha=0.12, linewidth=0.55)
    ax.tick_params(labelsize=22)
    add_colorbar(ax, sc, "Detections per object")
    fig.tight_layout()
    save_figure(fig, outdir / "temporal_sampling")


def plot_example(df: pd.DataFrame, outdir: Path) -> None:
    counts = df.groupby("source_file").size().sort_values(ascending=False)
    source_file = counts.index[0]
    sub = df[df["source_file"] == source_file].copy()

    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    sc = axes[0, 0].scatter(sub["X_Win"], sub["Y_Win"], c=sub["sep_arcsec"], s=70, cmap="rainbow", edgecolor="black", linewidth=0.3)
    axes[0, 0].set_xlabel("X_Win [pixel]")
    axes[0, 0].set_ylabel("Y_Win [pixel]")
    axes[0, 0].set_title("Accepted matches", fontsize=26)
    add_colorbar(axes[0, 0], sc, "Separation [arcsec]")

    sample = sub.sort_values("sep_arcsec", ascending=False).head(min(45, len(sub)))
    axes[0, 1].scatter(sample["X_Win"], sample["Y_Win"], s=24, color="#b2182b", alpha=0.85, linewidth=0)
    vector_scale = 2600.0
    for _, row in sample.iterrows():
        start = (row["X_Win"], row["Y_Win"])
        end = (
            row["X_Win"] + row["dra_cosdec_arcsec"] * vector_scale,
            row["Y_Win"] + row["ddec_arcsec"] * vector_scale,
        )
        axes[0, 1].add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=19,
                linewidth=2.6,
                color="#b2182b",
                alpha=0.95,
            )
        )
    axes[0, 1].set_xlabel("X_Win [pixel]")
    axes[0, 1].set_ylabel("Y_Win [pixel]")
    axes[0, 1].set_title("Residual vectors", fontsize=26)

    histogram(axes[1, 0], sub["sep_arcsec"], np.linspace(0, max(0.2, float(sub["sep_arcsec"].max()) * 1.05), 25), "Separation [arcsec]", color="#59a14f")
    axes[1, 1].scatter(sub["mag"], sub[ADOPTED_MAG], s=55, alpha=0.8, color="#4c78a8", edgecolor="black", linewidth=0.3)
    axes[1, 1].set_xlabel("Predicted magnitude [mag]")
    axes[1, 1].set_ylabel(f"{ADOPTED_MAG_LABEL} [mag]")
    axes[1, 1].tick_params(labelsize=22)
    fig.tight_layout()
    save_figure(fig, outdir / "example_crossmatch_diagnostic")


def add_box(ax, xy, w, h, text, kind) -> None:
    x, y = xy
    colors = {
        "input": ("#dceeff", "#4f81bd"),
        "process": ("#d9ead3", "#6aa84f"),
        "decision": ("#fff2cc", "#bf9000"),
        "reject": ("#f4cccc", "#cc0000"),
        "output": ("#eadcf8", "#8e7cc3"),
        "database": ("#d9eaf7", "#6d9eeb"),
    }
    face, edge = colors[kind]
    if kind == "decision":
        patch = Polygon([(x, y + h / 2), (x + w / 2, y + h), (x + w, y + h / 2), (x + w / 2, y)], closed=True, facecolor=face, edgecolor=edge, linewidth=2.0)
    elif kind in {"input", "output"}:
        skew = 0.15 * w
        patch = Polygon([(x + skew, y), (x + w, y), (x + w - skew, y + h), (x, y + h)], closed=True, facecolor=face, edgecolor=edge, linewidth=2.0)
    elif kind == "database":
        patch = Ellipse((x + w / 2, y + h / 2), w, h, facecolor=face, edgecolor=edge, linewidth=2.0)
    else:
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03,rounding_size=0.07", facecolor=face, edgecolor=edge, linewidth=2.0)
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=17)


def arrow(ax, start, end) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=18, linewidth=1.8, color="#444444"))


def plot_flowchart(outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(15, 20))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 20)
    ax.axis("off")

    nodes = {
        "catalog": ((3.2, 18.6), 3.6, 0.7, "Photometric\ncatalog", "input"),
        "header": ((0.4, 17.1), 3.0, 0.8, "FITS header\nWCS + time", "input"),
        "orbit": ((6.6, 17.1), 3.0, 0.8, "Lowell astorb\nSBDB / Horizons", "input"),
        "read": ((3.4, 16.9), 3.2, 0.8, "Read table\nand header", "process"),
        "wcs": ((3.4, 15.6), 3.2, 0.8, "WCS footprint\nfrom CCD corners", "process"),
        "epoch": ((3.4, 14.3), 3.2, 0.8, "Mid-exposure\nepoch", "process"),
        "query": ((3.4, 13.0), 3.2, 0.8, "Ephemeris\nscreening", "process"),
        "mag": ((3.55, 11.55), 2.9, 1.0, "Field-radius\npre-filter?", "decision"),
        "faint": ((7.0, 11.7), 2.4, 0.7, "Reject outside\nfield radius", "reject"),
        "pixel": ((3.4, 10.2), 3.2, 0.8, "Predicted RA/Dec\ninto pixels", "process"),
        "inside": ((3.55, 8.75), 2.9, 1.0, "Inside CCD\nfootprint?", "decision"),
        "offccd": ((7.0, 8.9), 2.4, 0.7, "Reject off-CCD\npredictions", "reject"),
        "match": ((3.4, 7.4), 3.2, 0.8, "Nearest-neighbor\nangular match", "process"),
        "accept": ((3.55, 5.95), 2.9, 1.0, "Separation\nthreshold?", "decision"),
        "unmatched": ((7.0, 6.1), 2.4, 0.7, "Reject unmatched\nrows", "reject"),
        "gaia": ((3.4, 4.6), 3.2, 0.8, "Gaia stationary-\nsource check", "process"),
        "stellar": ((7.0, 4.75), 2.4, 0.7, "Reject likely\nstellar contaminants", "reject"),
        "merge": ((3.4, 3.25), 3.2, 0.8, "Merge ephemeris\nand catalog row", "process"),
        "augment": ((3.4, 2.0), 3.2, 0.8, "SBDB + Horizons\naugmentation", "process"),
        "final": ((3.1, 0.6), 3.8, 0.8, "Recovered sample\nstatistics + figures", "output"),
    }
    for xy, w, h, text, kind in nodes.values():
        add_box(ax, xy, w, h, text, kind)

    arrow(ax, (5.0, 18.6), (5.0, 17.7))
    arrow(ax, (3.4, 17.5), (3.4, 17.3))
    arrow(ax, (6.6, 17.5), (6.6, 17.3))
    for y0, y1 in [(16.9, 16.4), (15.6, 15.1), (14.3, 13.8), (13.0, 12.55), (11.55, 11.0), (10.2, 9.75), (8.75, 8.2), (7.4, 6.95), (5.95, 5.4), (4.6, 4.05), (3.25, 2.8), (2.0, 1.4)]:
        arrow(ax, (5.0, y0), (5.0, y1))
    arrow(ax, (6.45, 12.05), (7.0, 12.05))
    arrow(ax, (6.45, 9.25), (7.0, 9.25))
    arrow(ax, (6.45, 6.45), (7.0, 6.45))
    arrow(ax, (6.6, 5.0), (7.0, 5.1))
    ax.text(6.62, 12.28, "no", fontsize=16)
    ax.text(6.62, 9.48, "no", fontsize=16)
    ax.text(6.62, 6.68, "no", fontsize=16)
    ax.text(6.62, 5.3, "flagged", fontsize=16)
    ax.text(4.65, 11.25, "yes", fontsize=16)
    ax.text(4.65, 8.45, "yes", fontsize=16)
    ax.text(4.65, 5.65, "yes", fontsize=16)
    ax.text(4.65, 4.3, "clean", fontsize=16)
    ax.text(5.0, 19.55, "Known-object processing", ha="center", va="center", fontsize=30)
    save_figure(fig, outdir / "method_flowchart_styled")


def write_mermaid(outdir: Path) -> None:
    text = """flowchart TD
    A[/Per-exposure photometric catalog/]:::input
    B[/FITS header metadata: WCS, DATE-OBS, EXPTIME/]:::input
    C[/External orbit resources: Lowell astorb, SBDB, JPL Horizons/]:::input
    D([Read source table and header]):::process
    E([Build WCS and project CCD corners]):::process
    F([Compute mid-exposure epoch]):::process
    G([Estimate padded field search radius]):::process
    H([Ephemeris screening]):::process
    I{{Within padded field radius?}}:::decision
    J([Transform predictions to pixel coordinates]):::process
    K{{Inside CCD footprint?}}:::decision
    L([Nearest-neighbor angular cross-match]):::process
    M{{Separation threshold?}}:::decision
    G1([Gaia stationary-source check]):::process
    N([Merge prediction row with catalog source row]):::process
    O([Merge exposure/night products]):::process
    P([SBDB and Horizons augmentation]):::process
    Q[/Recovered known-asteroid sample/]:::output
    R[/Statistics, figures, light-curve inputs/]:::output
    X1([Reject predictions outside field radius]):::reject
    X2([Reject off-CCD predictions]):::reject
    X3([Reject unmatched rows]):::reject
    X4([Reject likely stellar contaminants]):::reject
    A --> D
    B --> D
    C --> H
    D --> E --> G --> H
    D --> F --> H
    H --> I
    I -- yes --> J --> K
    I -- no --> X1
    K -- yes --> L --> M
    K -- no --> X2
    M -- yes --> G1 --> N --> O --> P --> Q --> R
    M -- no --> X3
    G1 -- flagged --> X4
    classDef input fill:#dceeff,stroke:#4f81bd,stroke-width:2px;
    classDef process fill:#d9ead3,stroke:#6aa84f,stroke-width:2px;
    classDef decision fill:#fff2cc,stroke:#bf9000,stroke-width:2px;
    classDef reject fill:#f4cccc,stroke:#cc0000,stroke-width:2px;
    classDef output fill:#eadcf8,stroke:#8e7cc3,stroke-width:2px;
"""
    (outdir / "method_flowchart_styled.mmd").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate paper figures and tables from the recovered known-asteroid sample.")
    parser.add_argument("fits_path", nargs="?", default="gotta_asteroids.fits")
    parser.add_argument("--outdir", default="paper_draft")
    parser.add_argument("--paper-version", default="v3")
    args = parser.parse_args()

    root = Path(args.outdir)
    figdir = root / f"figures_{args.paper_version}"
    tabledir = root / f"tables_{args.paper_version}"
    setup_style()
    df = add_derived(load_table(Path(args.fits_path)))
    make_tables(df, tabledir)
    plot_flowchart(figdir)
    write_mermaid(figdir)
    workflow_src = Path("outputs/known_object_processing.png")
    if workflow_src.exists():
        shutil.copy2(workflow_src, figdir / "known_object_processing.png")
    previous_figdir = root / "figures_v4"
    for static_name in [
        "sitian_pilot_asteroid_statistics",
        "gotta_radec_healpix_nside64",
        "asteroid_orbits",
    ]:
        for suffix in [".png", ".pdf"]:
            src = previous_figdir / f"{static_name}{suffix}"
            if src.exists():
                shutil.copy2(src, figdir / src.name)
    plot_photometry(df, figdir)
    plot_orbit_class_pies(df, figdir)
    plot_astrometry(df, figdir)
    plot_geometry(df, figdir)
    plot_temporal(df, figdir)
    plot_example(df, figdir)
    print(f"Wrote figures to {figdir}")
    print(f"Wrote tables to {tabledir}")


if __name__ == "__main__":
    main()
