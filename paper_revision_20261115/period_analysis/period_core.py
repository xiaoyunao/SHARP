"""Portable, fixed-sample weighted Fourier period search (time in days).

Adapted from the numerical design in period_delivery_20261127; see README.
No reference periods, network calls, clipping, deployment or database writes.
"""
from dataclasses import dataclass

import numpy as np
from astropy.timeseries import LombScargle
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks


@dataclass
class Fit:
    frequency: float
    coefficients: np.ndarray
    residual: np.ndarray
    chi2: float
    bic: float
    order: int
    condition: float
    levels: np.ndarray


def design(t, night, frequency, order=2, mode="night"):
    levels, index = np.unique(night, return_inverse=True)
    base = np.eye(len(levels))[index] if mode == "night" else np.ones((len(t), 1))
    angle = 2 * np.pi * np.asarray(t) * frequency
    harmonics = np.column_stack([v for k in range(1, order + 1)
                                for v in (np.cos(k * angle), np.sin(k * angle))])
    return np.column_stack((base, harmonics)), levels


def fit(t, y, err, night, frequency, order=2, mode="night"):
    x, levels = design(t, night, frequency, order, mode)
    if len(y) <= x.shape[1] + 2:
        return None
    coefficients, _, rank, singular = np.linalg.lstsq(x / err[:, None], y / err, rcond=1e-10)
    condition = singular[0] / singular[-1] if singular[-1] > 0 else np.inf
    if rank != x.shape[1] or condition > 1e8:
        return None
    residual = y - x @ coefficients
    chi2 = float(np.sum((residual / err) ** 2))
    # Known (floored) heteroscedastic errors; frequency is an estimated parameter.
    bic = chi2 + (x.shape[1] + 1) * np.log(len(y))
    return Fit(float(frequency), coefficients, residual, chi2, float(bic), order, condition, levels)


def refine(t, y, err, night, seed, width, order=2, mode="night"):
    def objective(f):
        result = fit(t, y, err, night, f, order, mode)
        return result.chi2 if result else 1e100

    grid = np.linspace(max(1e-6, seed - width), seed + width, 13)
    values = np.array([objective(f) for f in grid])
    best = int(np.argmin(values))
    boundary = best in (0, len(grid) - 1)
    frequency = grid[best]
    if not boundary:
        result = minimize_scalar(objective, bounds=(grid[best - 1], grid[best + 1]),
                                 method="bounded", options={"xatol": 1e-10})
        frequency = result.x
    return fit(t, y, err, night, frequency, order, mode), boundary


def phase_gap(t, frequency):
    phase = np.sort(np.mod(np.asarray(t) * frequency, 1))
    return float(np.max(np.diff(np.r_[phase, phase[0] + 1])))


def search(t, y, err, night, min_period_h=1., max_period_h=120., oversample=8):
    """Profile all order-2 grid fits; supplement LS peaks with half/daily aliases.

    Night offsets are nuisance parameters at *every* trial frequency. LS is
    only a proposal diagnostic; ranking uses the same points and Fourier model.
    """
    t, y, err, night = map(np.asarray, (t, y, err, night))
    t = t - t.min()
    baseline = np.ptp(t)
    if baseline <= 0 or not 0 < min_period_h < max_period_h:
        raise ValueError("Need positive baseline and ordered positive period bounds")
    lo, hi = 24 / max_period_h, 24 / min_period_h
    step = 1 / (oversample * baseline)
    frequency = np.linspace(lo, hi, int(np.ceil((hi - lo) / step)) + 1)
    results = [fit(t, y, err, night, f) for f in frequency]
    chi2 = np.array([r.chi2 if r else np.inf for r in results])
    # Relative night means for LS proposal only; no observations are deleted.
    centered = y.copy()
    for label in np.unique(night):
        mask = night == label
        centered[mask] -= np.average(y[mask], weights=1 / err[mask] ** 2)
    ls = LombScargle(t, centered, err).power(frequency, method="cython")
    peaks = find_peaks(-chi2)[0]
    ls_peaks = find_peaks(ls)[0]
    seeds = [(frequency[i], "Fourier_grid") for i in peaks[np.argsort(chi2[peaks])[:12]]]
    seeds += [(frequency[i], "LS") for i in ls_peaks[np.argsort(ls[ls_peaks])[-6:]]]
    seeds += [(frequency[i], "search_boundary") for i in (0, len(frequency) - 1)]
    proposals = seeds.copy()
    for seed, origin in seeds:
        for factor in (.5, 1., 2.):
            for daily in (-1., 0., 1.):
                proposals.append((factor * seed + daily, f"{origin}:x{factor}:daily{daily:+}"))
    # Distinct local basins remain separate. Fine refinement cannot jump aliases.
    unique = []
    for f, origin in proposals:
        if lo <= f <= hi and all(abs(f - a) > .3 * step for a, _ in unique):
            unique.append((f, origin))
    candidates = []
    for seed, origin in unique:
        width = min(.5 * step, seed - lo, hi - seed)
        if width > 0:
            result, edge = refine(t, y, err, night, seed, width)
        else:
            result, edge = fit(t, y, err, night, seed), True
        if result is None:
            continue
        candidates.append((result, edge, origin, width))
    candidates.sort(key=lambda r: r[0].bic)
    merged = []
    for candidate in candidates:
        if all(abs(candidate[0].frequency - prior[0].frequency) > .75 * step for prior in merged):
            merged.append(candidate)
    window = np.abs(np.mean(np.exp(2j * np.pi * frequency[:, None] * t), axis=1)) ** 2
    spectrum = dict(frequency_cpd=frequency, fourier_chi2=chi2, ls_power=ls, spectral_window_power=window)
    return merged, spectrum


def half_cycle_test(t, y, err, night, frequency):
    """At proposed full period, compare even harmonics (2,4) to all (1..4).

    The even-only model exactly contains an order-2 model at half the period.
    This avoids unfairly rejecting 2P merely because its fit needs order 4.
    Positive delta BIC favors measured half-cycle asymmetry, not a spin claim.
    """
    full = fit(t, y, err, night, frequency, order=4)
    if full is None:
        return dict(half_cycle_delta_bic=None, half_cycle_status="unconstrained")
    x, levels = design(t, night, frequency, order=4)
    nbase = len(levels)
    keep = list(range(nbase)) + [nbase + i for i in (2, 3, 6, 7)]
    even_x = x[:, keep]
    coefficients, _, rank, _ = np.linalg.lstsq(even_x / err[:, None], y / err, rcond=1e-10)
    if rank != even_x.shape[1]:
        return dict(half_cycle_delta_bic=None, half_cycle_status="unconstrained")
    even_chi2 = np.sum(((y - even_x @ coefficients) / err) ** 2)
    scale = max(1., full.chi2 / (len(y) - len(full.coefficients)))
    delta_bic = float((even_chi2 - full.chi2) / scale - 4 * np.log(len(y)))
    return dict(half_cycle_delta_bic=delta_bic,
                half_cycle_status="asymmetry_candidate" if delta_bic > 10 else "half_double_ambiguous")


def held_night_test(t, y, err, night, frequency, width, order=2):
    """Train-only local frequency refit, conditional shape prediction on test night.

    A test-night constant is fitted from that night's data. This is explicitly
    conditional on its mean brightness, not unconditional prediction.
    """
    squared, weights, null_squared = [], [], []
    used = 0
    for label in np.unique(night):
        test = night == label
        train = ~test
        if test.sum() < 3 or len(np.unique(night[train])) < 3:
            continue
        if phase_gap(t[train], frequency) > .5:
            continue
        result, edge = refine(t[train], y[train], err[train], night[train],
                              frequency, max(width, 1e-7), order)
        if result is None or edge:
            continue
        x, _ = design(t[test], night[test], result.frequency, order)
        model = x[:, -2 * order:] @ result.coefficients[-2 * order:]
        w = 1 / err[test] ** 2
        offset = np.average(y[test] - model, weights=w)
        squared.extend((y[test] - model - offset) ** 2)
        null_squared.extend((y[test] - np.average(y[test], weights=w)) ** 2)
        weights.extend(w)
        used += 1
    if not used:
        return dict(cv_nights=0, cv_conditional_rmse_mag=None, cv_null_rmse_mag=None)
    return dict(cv_nights=used,
                cv_conditional_rmse_mag=float(np.sqrt(np.average(squared, weights=weights))),
                cv_null_rmse_mag=float(np.sqrt(np.average(null_squared, weights=weights))))
