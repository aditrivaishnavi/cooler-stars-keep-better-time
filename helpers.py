"""Everything the notebooks need, in one short file.

The method has four steps, and each is one function here:

  1. load_sector(path)              -> the brightness curve of one TESS sector
  2. fit_period(time, flux, noise, p0) -> one period (and the leftover) from that sector
  3. measure_star(...)              -> the star's wander, and its two floors from fake steady stars
  4. robust_fit(X, y)               -> a straight-line fit whose error bars do not trust equal noise

Settings are the same as the full pipeline (research/src/teff_stability_raw_fit.py)
so that the numbers can be compared one to one.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RESULTS = HERE / "results"
FIGURES = HERE / "figures"

# --- settings (identical to the full pipeline) ---------------------------------------------
SEARCH_SPAN = 0.30        # look for the period between p0/1.3 and p0/0.7 (avoids the half/double-period trap)
FINE_SPAN = 0.03          # then refine within +/-3% of the coarse peak
MINIMUM_POINTS = 200      # a sector with fewer good points is skipped
MINIMUM_SECTORS = 5       # a star with fewer usable sectors is skipped
INJECTIONS = 12           # fake steady stars per real star
SEED = 20261203           # random seed; the star's TIC id is added so every star has its own stream
FIXED_THRESHOLD = 0.015   # rad/day; the "fixed bar" used for the 83% vs 18% statement
FULLY_CONVECTIVE_TEFF = 3400.0


# --- 0. the star lists ----------------------------------------------------------------------
def read_star_list(name: str) -> list[dict]:
    """name is 'original_450' or 'replication_422'. Returns one dict per star."""
    rows = list(csv.DictReader((DATA / "star_lists" / f"{name}.csv").open()))
    for row in rows:
        row["tic_id"] = int(row["tic_id"])
        row["teff"] = float(row["teff"])
        row["tmag"] = float(row["tmag"])
        row["period"] = float(row["period"])          # the catalogue period, our starting guess
        row["color"] = float(row["color"])
        row["teff_bin_low"] = int(row["teff_bin_low"])
        row["sectors"] = int(row["sectors"])                # accepted sectors, from the catalogue
        row["sectors_allowed"] = {int(s) for s in row["sector_list"].split(",") if s}
    return rows


def sector_files(cache: Path, tic_id: int) -> list[Path]:
    return sorted(cache.glob(f"*-{tic_id:016d}_*_lc.fits"))


# --- 1. one sector ---------------------------------------------------------------------------
def load_sector(path: Path):
    """Returns (sector number, time in days, brightness relative to 1, point-to-point noise) or None."""
    with fits.open(path, memmap=False) as hdus:
        sector = int(hdus[0].header["SECTOR"])
        time = np.asarray(hdus[1].data["TIME"], float)
        flux = np.asarray(hdus[1].data["FLUX"], float)
        noise = float(hdus[1].header.get("P2PRMS", np.nan))
    keep = np.isfinite(time) & np.isfinite(flux)
    time, flux = time[keep], flux[keep]
    if len(time) < MINIMUM_POINTS:
        return None
    if not np.isfinite(noise) or noise <= 0:
        noise = 1.4826 * np.median(np.abs(np.diff(flux))) / np.sqrt(2)
    return sector, time, flux, max(noise, 1e-5)


# --- 2. one period from one sector -------------------------------------------------------------
def fit_period(time, flux, noise, p0):
    """Fit one steady wave near the catalogue period p0.

    Returns (period, amplitude, leftover) or None.
    * take out a straight-line drift
    * coarse search: which period between 0.77 and 1.43 x p0 explains the most wiggle (Lomb-Scargle)
    * fine search: +/-3% around it, with the peak's exact position found from a parabola through
      the three highest points
    * amplitude = half the 5th-to-95th percentile height of the fitted wave
    * leftover = data minus wave (this is what the fake stars are built from)
    """
    error = np.full(len(time), noise)
    weight = 1 / np.maximum(error, 0.003) ** 2
    centred = time - time.mean()
    design = np.column_stack([np.ones(len(time)), centred])
    try:
        coefficient = np.linalg.solve(design.T @ (weight[:, None] * design), design.T @ (weight * flux))
    except np.linalg.LinAlgError:
        return None
    detrended = flux - design @ coefficient
    periodogram = LombScargle(time, detrended, error, fit_mean=True, center_data=True)
    f0 = 1 / p0
    coarse = np.linspace(f0 / (1 + SEARCH_SPAN), f0 / (1 - SEARCH_SPAN), 720)
    f_coarse = coarse[int(np.nanargmax(periodogram.power(coarse, method="auto")))]
    fine = np.linspace(f_coarse * (1 - FINE_SPAN), f_coarse * (1 + FINE_SPAN), 601)
    power = periodogram.power(fine, method="auto")
    i = int(np.nanargmax(power))
    offset = 0.0
    if 0 < i < len(fine) - 1:
        y0, y1, y2 = power[i - 1], power[i], power[i + 1]
        denominator = y0 - 2 * y1 + y2
        offset = float(np.clip(0.5 * (y0 - y2) / denominator, -1, 1)) if denominator != 0 else 0.0
    f_peak = fine[i] + offset * (fine[1] - fine[0])
    wave = periodogram.model(time, f_peak)
    leftover = detrended - wave
    amplitude = 0.5 * (np.percentile(wave, 95) - np.percentile(wave, 5))
    return float(1 / f_peak), float(amplitude), leftover


def angular_mad(periods):
    """The wander: robust scatter of the angular speed 2*pi/P across sectors, in rad/day."""
    omega = 2 * np.pi / np.asarray(periods)
    return float(1.4826 * np.median(np.abs(omega - np.median(omega))))


# --- 3. one star: wander and its two floors ------------------------------------------------------
def measure_star(row: dict, cache: Path) -> dict | None:
    """Fit every accepted sector, then build fake steady stars to measure the floor two ways."""
    tic, p0 = row["tic_id"], row["period"]
    sectors = []
    for path in sector_files(cache, tic):
        loaded = load_sector(path)
        if loaded is None:
            continue
        sector, time, flux, noise = loaded
        if sector not in row["sectors_allowed"]:
            continue
        fitted = fit_period(time, flux, noise, p0)
        if fitted is None:
            continue
        period, amplitude, leftover = fitted
        sectors.append(dict(sector=sector, time=time, noise=noise, period=period, amplitude=amplitude, leftover=leftover))
    if len(sectors) < MINIMUM_SECTORS:
        return None
    sectors.sort(key=lambda s: s["sector"])
    periods = np.array([s["period"] for s in sectors])
    wander = angular_mad(periods)

    # The fake steady star, two ways of re-using the leftover as the noise.
    # Each fake star is one steady wave (random starting phase) threaded through every sector,
    # measured with the same fit_period; its angular_mad is what a star that did nothing would show.
    floors = {}
    for name in ("white", "shifted"):
        rng = np.random.default_rng(SEED + tic)        # same seeds as the full pipeline
        fake_scatters = []
        for _ in range(INJECTIONS):
            phase = rng.uniform(0, 2 * np.pi)
            fake_periods = []
            for s in sectors:
                wave = s["amplitude"] * np.sin(2 * np.pi * s["time"] / p0 + phase)
                if name == "white":
                    noise_copy = rng.permutation(s["leftover"])               # shuffled: destroys slow wiggles
                else:
                    noise_copy = np.roll(s["leftover"], int(rng.integers(1, len(s["leftover"]) - 1)))  # slid: keeps them
                fitted = fit_period(s["time"], wave + noise_copy, s["noise"], p0)
                fake_periods.append(fitted[0] if fitted else np.nan)
            fake_periods = np.array(fake_periods)
            if np.isfinite(fake_periods).sum() >= MINIMUM_SECTORS:
                fake_scatters.append(angular_mad(fake_periods[np.isfinite(fake_periods)]))
        floors[name] = float(np.median(fake_scatters)) if fake_scatters else np.nan

    early = periods[np.array([s["sector"] for s in sectors]) <= 60]
    late = periods[np.array([s["sector"] for s in sectors]) >= 61]
    out = dict(
        tic_id=tic, teff=row["teff"], color=row["color"], tmag=row["tmag"], period=p0,
        sectors=len(sectors), median_amplitude=float(np.median([s["amplitude"] for s in sectors])),
        wander=wander, floor_white=floors["white"], floor_shifted=floors["shifted"],
        early_wander=angular_mad(early) if len(early) >= 3 else np.nan,
        late_wander=angular_mad(late) if len(late) >= 3 else np.nan,
    )
    for name in ("white", "shifted"):
        floor = floors[name]
        out[f"ratio_{name}"] = wander / floor if floor > 0 else np.nan
        out[f"excess_{name}"] = float(np.sqrt(max(wander ** 2 - floor ** 2, 0.0)))
    out["sector_periods"] = ",".join(f"{s['sector']}:{s['period']:.5f}" for s in sectors)
    return out


# --- 4. statistics ---------------------------------------------------------------------------------
def robust_fit(matrix, response):
    """Least squares with heteroskedasticity-robust (HC1) standard errors.

    'Robust' means the error bars do not assume every star is equally noisy.
    Returns (coefficients, standard errors).
    """
    inverse = np.linalg.pinv(matrix.T @ matrix)
    coefficient = inverse @ (matrix.T @ response)
    residual = response - matrix @ coefficient
    n, p = matrix.shape
    meat = matrix.T @ ((residual * residual)[:, None] * matrix)
    covariance = inverse @ meat @ inverse * n / max(n - p, 1)
    return coefficient, np.sqrt(np.diag(covariance))


def design_matrix(stars, with_step=False):
    """The regression: log(excess ratio) against temperature, holding period, amplitude,
    brightness and sector count fixed. Same centring as the full pipeline."""
    columns = [
        np.ones(len(stars)),
        (stars["teff"].to_numpy() - 4500) / 1000,
        np.log2(stars["period"].to_numpy() / 3),
        np.log10(stars["median_amplitude"].to_numpy()),
        stars["tmag"].to_numpy() - 13.5,
        np.log10(stars["sectors"].to_numpy() / 10),
    ]
    names = ["intercept", "teff_per_1000K", "log2_period", "log10_amplitude", "tmag", "log10_sectors"]
    if with_step:
        columns.append((stars["teff"].to_numpy() < FULLY_CONVECTIVE_TEFF).astype(float))
        names.append("step_below_3400K")
    return np.column_stack(columns), names


def period_slope(stars, floor):
    """Slope of log10(excess in quadrature) per doubling of the period, holding the other
    covariates fixed, for stars with a positive excess. Same design as the full pipeline."""
    positive = stars[stars[f"excess_{floor}"] > 0]
    columns = [
        np.ones(len(positive)),
        np.log2(positive["period"].to_numpy() / 3),
        np.log10(positive["median_amplitude"].to_numpy()),
        positive["tmag"].to_numpy() - 13.5,
        np.log10(positive["sectors"].to_numpy() / 10),
    ]
    coefficient, se = robust_fit(np.column_stack(columns), np.log10(positive[f"excess_{floor}"].to_numpy()))
    return dict(n=len(positive), estimate=float(coefficient[1]), hc_se=float(se[1]), t=float(coefficient[1] / se[1]))
