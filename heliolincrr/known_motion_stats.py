#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize nightly known-asteroid sky-motion rates for heliolincrr tuning."
    )
    parser.add_argument("night", help="Target night in YYYYMMDD format")
    parser.add_argument("--processed-root", default="/processed1")
    parser.add_argument("--matched-path", default=None, help="Override nightly matched_asteroids.fits path")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory for stats products (default: /pipeline/xiaoyunao/data/heliolincrr/<night>/analysis)",
    )
    parser.add_argument("--min-detections", type=int, default=2, help="Minimum detections per asteroid to enter rate stats")
    return parser.parse_args()


def clean_text(value: object) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="ignore").strip()
    return str(value).strip()


def wrap_delta_deg(delta_deg: np.ndarray) -> np.ndarray:
    return ((np.asarray(delta_deg, dtype=float) + 180.0) % 360.0) - 180.0


def compute_pair_rates(tbl: Table, min_detections: int) -> pd.DataFrame:
    df = tbl.to_pandas().copy()
    df["name"] = df["name"].map(clean_text)

    rows: list[dict[str, float | str]] = []
    for name, group in df.groupby("name", sort=True):
        if len(group) < min_detections:
            continue
        group = group.sort_values("MJD")
        ra = group["RA_Win"].to_numpy(dtype=float)
        dec = group["DEC_Win"].to_numpy(dtype=float)
        mjd = group["MJD"].to_numpy(dtype=float)

        dt_hr = np.diff(mjd) * 24.0
        if dt_hr.size == 0:
            continue

        dra = wrap_delta_deg(np.diff(ra)) * np.cos(np.deg2rad(0.5 * (dec[:-1] + dec[1:])))
        ddec = np.diff(dec)
        sep_arcsec = 3600.0 * np.hypot(dra, ddec)

        valid = np.isfinite(dt_hr) & (dt_hr > 0) & np.isfinite(sep_arcsec)
        if not np.any(valid):
            continue

        mjd1 = mjd[:-1][valid]
        mjd2 = mjd[1:][valid]
        dt_hr = dt_hr[valid]
        sep_arcsec = sep_arcsec[valid]
        rate_arcsec_hr = sep_arcsec / dt_hr

        for a, b, dt, sep, rate in zip(mjd1, mjd2, dt_hr, sep_arcsec, rate_arcsec_hr):
            rows.append(
                {
                    "name": name,
                    "mjd1": float(a),
                    "mjd2": float(b),
                    "dt_hr": float(dt),
                    "sep_arcsec": float(sep),
                    "rate_arcsec_hr": float(rate),
                    "rate_arcsec_min": float(rate / 60.0),
                }
            )

    if not rows:
        return pd.DataFrame(
            columns=["name", "mjd1", "mjd2", "dt_hr", "sep_arcsec", "rate_arcsec_hr", "rate_arcsec_min"]
        )
    return pd.DataFrame.from_records(rows)


def finite_stats(values: np.ndarray) -> dict[str, float | int | None]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {
            "count": 0,
            "min": None,
            "p01": None,
            "p05": None,
            "p10": None,
            "p50": None,
            "p90": None,
            "p95": None,
            "p99": None,
            "max": None,
            "mean": None,
        }
    return {
        "count": int(arr.size),
        "min": float(np.min(arr)),
        "p01": float(np.percentile(arr, 1)),
        "p05": float(np.percentile(arr, 5)),
        "p10": float(np.percentile(arr, 10)),
        "p50": float(np.percentile(arr, 50)),
        "p90": float(np.percentile(arr, 90)),
        "p95": float(np.percentile(arr, 95)),
        "p99": float(np.percentile(arr, 99)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
    }


def threshold_counts(values: np.ndarray, thresholds: list[float]) -> dict[str, int]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return {str(th): int(np.count_nonzero(arr <= th)) for th in thresholds}


def plot_rate_histogram(pair_df: pd.DataFrame, out_path: Path, night: str) -> None:
    rates = pair_df["rate_arcsec_hr"].to_numpy(dtype=float)
    rates = rates[np.isfinite(rates) & (rates >= 0)]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    if rates.size > 0:
        hi = max(5.0, float(np.percentile(rates, 99.5)))
        bins = np.linspace(0.0, hi, 50)
        ax.hist(rates, bins=bins, color="#2a9d8f", edgecolor="white", alpha=0.95)
        for q, color in [(50, "#d62828"), (90, "#1d3557"), (95, "#6a4c93"), (99, "#ff7f11")]:
            val = float(np.percentile(rates, q))
            ax.axvline(val, color=color, linestyle="--", linewidth=1.8, label=f"p{q}={val:.2f}")
        ax.legend(loc="upper right")
        ax.set_xlim(0.0, hi)
    ax.set_xlabel("Sky motion rate (arcsec / hour)")
    ax.set_ylabel("Consecutive matched pairs")
    ax.set_title("Known-asteroid Motion Histogram")

    ax = axes[1]
    if rates.size > 0:
        xs = np.sort(rates)
        ys = np.arange(1, xs.size + 1) / xs.size
        ax.plot(xs, ys, color="#457b9d", linewidth=2.0)
        for ref in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0]:
            ax.axvline(ref, color="#bbbbbb", linestyle=":", linewidth=1.0)
        ax.set_xlim(0.0, max(5.0, float(np.percentile(xs, 99.5))))
    ax.set_xlabel("Sky motion rate (arcsec / hour)")
    ax.set_ylabel("CDF")
    ax.set_title("Known-asteroid Motion CDF")
    ax.set_ylim(0.0, 1.0)

    fig.suptitle(f"Known-asteroid motion summary for {night}", fontsize=15)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()

    matched_path = (
        Path(args.matched_path)
        if args.matched_path is not None
        else Path(args.processed_root) / args.night / "L4" / f"{args.night}_matched_asteroids.fits"
    )
    outdir = (
        Path(args.outdir)
        if args.outdir is not None
        else Path("/pipeline/xiaoyunao/data/heliolincrr") / args.night / "analysis"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    if not matched_path.exists():
        raise FileNotFoundError(f"matched FITS not found: {matched_path}")

    tbl = Table.read(matched_path)
    pair_df = compute_pair_rates(tbl, min_detections=args.min_detections)

    rate_stats = finite_stats(pair_df["rate_arcsec_hr"].to_numpy(dtype=float))
    dt_stats = finite_stats(pair_df["dt_hr"].to_numpy(dtype=float))

    summary = {
        "night": args.night,
        "matched_path": str(matched_path),
        "n_matched_rows": int(len(tbl)),
        "n_unique_asteroids": int(len({clean_text(x) for x in tbl["name"]})),
        "n_rate_pairs": int(len(pair_df)),
        "min_detections": int(args.min_detections),
        "rate_arcsec_hr": rate_stats,
        "dt_hr": dt_stats,
        "count_rate_le_arcsec_hr": threshold_counts(
            pair_df["rate_arcsec_hr"].to_numpy(dtype=float),
            [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0],
        ),
    }

    json_path = outdir / f"{args.night}_known_motion_summary.json"
    csv_path = outdir / f"{args.night}_known_motion_pairs.csv"
    png_path = outdir / f"{args.night}_known_motion_summary.png"

    json_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    pair_df.to_csv(csv_path, index=False)
    plot_rate_histogram(pair_df, png_path, args.night)

    print(f"[info] matched_path={matched_path}")
    print(f"[info] n_matched_rows={summary['n_matched_rows']}")
    print(f"[info] n_unique_asteroids={summary['n_unique_asteroids']}")
    print(f"[info] n_rate_pairs={summary['n_rate_pairs']}")
    print(f"[info] rate_p50_arcsec_hr={summary['rate_arcsec_hr']['p50']}")
    print(f"[info] rate_p90_arcsec_hr={summary['rate_arcsec_hr']['p90']}")
    print(f"[info] rate_p95_arcsec_hr={summary['rate_arcsec_hr']['p95']}")
    print(f"[info] rate_p99_arcsec_hr={summary['rate_arcsec_hr']['p99']}")
    print(f"[info] json={json_path}")
    print(f"[info] csv={csv_path}")
    print(f"[info] png={png_path}")


if __name__ == "__main__":
    main()
