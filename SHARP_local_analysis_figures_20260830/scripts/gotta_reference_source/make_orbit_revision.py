#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.lines import Line2D


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FITS_PATH = ROOT / "gotta_asteroids.fits"
OUT_PNG = HERE / "asteroid_orbits_revised.png"
OUT_PDF = HERE / "asteroid_orbits_revised.pdf"


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


def as_text_array(values) -> np.ndarray:
    out = []
    for value in values:
        if isinstance(value, (bytes, np.bytes_)):
            out.append(value.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(value).strip())
    return np.asarray(out, dtype=str)


def str_to_float(values, unit: str | None = None) -> np.ndarray:
    out = []
    for value in values:
        if value is None:
            out.append(np.nan)
            continue
        if isinstance(value, (bytes, np.bytes_)):
            text = value.decode("utf-8", errors="ignore").strip()
        else:
            text = str(value).strip()
        if text in ("", "--", "nan", "NaN"):
            out.append(np.nan)
            continue
        if unit and text.endswith(unit):
            text = text[: -len(unit)].strip()
        try:
            out.append(float(text))
        except ValueError:
            match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text)
            out.append(float(match.group(0)) if match else np.nan)
    return np.asarray(out, dtype=float)


def load_columns(path: Path) -> dict[str, np.ndarray]:
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        return {name: np.asarray(data[name]) for name in data.names}


def orbit_group(classes: np.ndarray) -> np.ndarray:
    mapping = {
        "Main-belt Asteroid": "MBA",
        "Outer Main-belt Asteroid": "OMB",
        "Jupiter Trojan": "TJN",
        "Inner Main-belt Asteroid": "IMB",
        "Mars-crossing Asteroid": "MCA",
        "Amor": "NEA",
        "Apollo": "NEA",
        "Aten": "NEA",
        "TransNeptunian Object": "TNO",
    }
    return np.asarray([mapping.get(cls, "Other") for cls in classes], dtype=str)


def plot_orbits(columns: dict[str, np.ndarray], out_png: Path, out_pdf: Path) -> None:
    a = str_to_float(columns["orbit_elements_a"], unit="AU")
    e = str_to_float(columns["orbit_elements_e"])
    inc = str_to_float(columns["orbit_elements_i"], unit="deg")
    classes = as_text_array(columns["object_orbit_class_name"])
    groups = orbit_group(classes)
    valid = np.isfinite(a) & np.isfinite(e) & np.isfinite(inc) & (a > 0) & (inc > 0)

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.size"] = 30

    fig, axes = plt.subplots(1, 2, figsize=(18, 9))
    ax1, ax2 = axes

    for group in GROUP_ORDER:
        mask = valid & (groups == group)
        if not np.any(mask):
            continue
        marker, color, alpha, size = GROUP_STYLE[group]
        ax1.scatter(a[mask], e[mask], marker=marker, color=color, s=size, alpha=alpha, linewidths=0)
        ax2.scatter(a[mask], inc[mask], marker=marker, color=color, s=max(size, 28), alpha=alpha, linewidths=0)

    a_line = np.logspace(np.log10(1.3), 2, 100)
    e_line = 1 - 1.3 / a_line
    ax1.plot(a_line, e_line, "k--", linewidth=2.3, label=r"$q=1.3\,\mathrm{AU}$")
    ax1.legend(fontsize=22, loc="lower right", frameon=False, handlelength=2.4)
    ax1.set_xscale("log")
    ax1.set_xlabel("Semimajor axis [AU]")
    ax1.set_ylabel("Eccentricity")

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Semimajor axis [AU]")
    ax2.set_ylabel("Inclination [°]")

    handles = []
    for group in GROUP_ORDER:
        marker, color, alpha, size = GROUP_STYLE[group]
        handles.append(
            Line2D(
                [0],
                [0],
                marker=marker,
                linestyle="None",
                markerfacecolor=color,
                markeredgecolor=color,
                alpha=max(alpha, 0.75),
                markersize=12 if group != "TNO" else 14,
                label=group,
            )
        )
    ax2.legend(
        handles=handles,
        fontsize=22,
        loc="best",
        frameon=False,
        handletextpad=0.4,
        borderpad=0.2,
        labelspacing=0.35,
        markerscale=1.15,
    )

    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)


def main() -> None:
    columns = load_columns(FITS_PATH)
    plot_orbits(columns, OUT_PNG, OUT_PDF)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
