"""Measure how much a star's rotation period wanders from one TESS sector to the next.

Both notebooks import this file. The measurement for one star is:

    load_sector    read one sector's light curve (a TARS FITS file)
    fit_period     fit one steady sine wave near the catalogue period and return its period
    angular_mad    the wander: robust scatter of the angular speed 2*pi/P across sectors
    noise_floor    the same measurement on fake, perfectly steady stars built from the star's own data
    measure_star   all of the above, for one star

Below that are the regression used in notebook 2 and the loader for the QLP light curves.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
EXAMPLES = DATA / "examples"
TARS_DIR = DATA / "lightcurves" / "tars"
QLP_DIR = DATA / "lightcurves" / "qlp"

SEARCH_SPAN = 0.30        # search for the period between 0.7 and 1.3 x p0
FINE_SPAN = 0.03          # then refine within +/-3% of the coarse peak frequency
MINIMUM_POINTS = 200      # skip a sector with fewer good points than this
MINIMUM_SECTORS = 5       # skip a star with fewer usable sectors than this
INJECTIONS = 12           # fake steady stars per real star
SEED = 20261203           # the star's TIC number is added, so every star has its own random stream
FIXED_BAR = 0.015         # rad/day
FULLY_CONVECTIVE_TEFF = 3400.0


def read_star_list(name: str) -> list[dict]:
    """name is 'original_450' or 'replication_422'."""
    rows = list(csv.DictReader((DATA / "star_lists" / f"{name}.csv").open()))
    for row in rows:
        row["tic_id"] = int(row["tic_id"])
        row["teff"] = float(row["teff"])
        row["tmag"] = float(row["tmag"])
        row["period"] = float(row["period"])      # TARS catalogue period, used as the starting guess
        row["color"] = float(row["color"])
        row["teff_bin_low"] = int(row["teff_bin_low"])
        row["sectors"] = int(row["sectors"])      # sectors TARS accepted for this star
        row["sectors_allowed"] = {int(s) for s in row["sector_list"].split(",") if s}
    return rows


def sector_files(folder: Path, tic_id: int) -> list[Path]:
    return sorted(Path(folder).glob(f"*-{tic_id:016d}_*_lc.fits"))


def sector_number(path: Path) -> int:
    return int(Path(path).name.split("_s")[1][:4])


def load_sector(path: Path):
    """(sector, time in days, relative flux, point-to-point noise), or None if too few points."""
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


def fit_period(time, flux, noise, p0):
    """Fit one steady wave near the catalogue period p0. Returns (period, amplitude, leftover) or None.

    A straight line is removed first. The coarse search covers periods of 0.7-1.3 x p0 with a
    Lomb-Scargle periodogram, so a star with two spots on opposite sides cannot pull the answer to
    p0/2. The fine search covers +/-3% in frequency around the coarse peak, and a parabola through
    the three highest points gives the peak position. Amplitude is half the 5th-95th percentile
    height of the fitted wave; the leftover is data minus wave.
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


def angular_mad(periods) -> float:
    """The wander, in rad/day: 1.4826 x median absolute deviation of 2*pi/P."""
    omega = 2 * np.pi / np.asarray(periods)
    return float(1.4826 * np.median(np.abs(omega - np.median(omega))))


def fit_sectors(row: dict, folder: Path) -> list[dict]:
    """Fit every sector TARS accepted for this star, in sector order."""
    sectors = []
    for path in sector_files(folder, row["tic_id"]):
        loaded = load_sector(path)
        if loaded is None:
            continue
        sector, time, flux, noise = loaded
        if sector not in row["sectors_allowed"]:
            continue
        fitted = fit_period(time, flux, noise, row["period"])
        if fitted is None:
            continue
        period, amplitude, leftover = fitted
        sectors.append(dict(sector=sector, time=time, flux=flux, noise=noise,
                            period=period, amplitude=amplitude, leftover=leftover))
    sectors.sort(key=lambda s: s["sector"])
    return sectors


def noise_floor(sectors: list[dict], p0: float, tic_id: int, kind: str):
    """Wander of a fake star that rotates perfectly steadily at p0.

    Each fake is a sine wave with the star's own amplitude in every sector, at the star's own
    time stamps, plus the star's own leftover as the noise. kind='white' shuffles the leftover
    (keeps its size, destroys slow wiggles); kind='shifted' slides it round in time by a random
    amount (keeps the slow wiggles). Returns (floor, list of the fake periods for each injection).
    """
    rng = np.random.default_rng(SEED + tic_id)
    scatters, runs = [], []
    for _ in range(INJECTIONS):
        phase = rng.uniform(0, 2 * np.pi)
        fake_periods = []
        for s in sectors:
            wave = s["amplitude"] * np.sin(2 * np.pi * s["time"] / p0 + phase)
            if kind == "white":
                noise = rng.permutation(s["leftover"])
            else:
                noise = np.roll(s["leftover"], int(rng.integers(1, len(s["leftover"]) - 1)))
            fitted = fit_period(s["time"], wave + noise, s["noise"], p0)
            fake_periods.append(fitted[0] if fitted else np.nan)
        fake_periods = np.array(fake_periods)
        runs.append(fake_periods)
        if np.isfinite(fake_periods).sum() >= MINIMUM_SECTORS:
            scatters.append(angular_mad(fake_periods[np.isfinite(fake_periods)]))
    return (float(np.median(scatters)) if scatters else np.nan), runs


def measure_star(row: dict, folder: Path = TARS_DIR) -> dict | None:
    """One row of results/stars_*.csv."""
    tic, p0 = row["tic_id"], row["period"]
    sectors = fit_sectors(row, folder)
    if len(sectors) < MINIMUM_SECTORS:
        return None
    periods = np.array([s["period"] for s in sectors])
    numbers = np.array([s["sector"] for s in sectors])
    wander = angular_mad(periods)
    floors = {kind: noise_floor(sectors, p0, tic, kind)[0] for kind in ("white", "shifted")}
    early, late = periods[numbers <= 60], periods[numbers >= 61]
    out = dict(
        tic_id=tic, teff=row["teff"], color=row["color"], tmag=row["tmag"], period=p0,
        sectors=len(sectors), median_amplitude=float(np.median([s["amplitude"] for s in sectors])),
        wander=wander, floor_white=floors["white"], floor_shifted=floors["shifted"],
        early_wander=angular_mad(early) if len(early) >= 3 else np.nan,
        late_wander=angular_mad(late) if len(late) >= 3 else np.nan,
    )
    for kind, floor in floors.items():
        out[f"ratio_{kind}"] = wander / floor if floor > 0 else np.nan
        out[f"excess_{kind}"] = float(np.sqrt(max(wander ** 2 - floor ** 2, 0.0)))
    out["sector_periods"] = ",".join(f"{s['sector']}:{s['period']:.5f}" for s in sectors)
    return out


def robust_fit(matrix, response):
    """Least squares with HC1 standard errors, which do not assume every star is equally noisy."""
    inverse = np.linalg.pinv(matrix.T @ matrix)
    coefficient = inverse @ (matrix.T @ response)
    residual = response - matrix @ coefficient
    n, p = matrix.shape
    meat = matrix.T @ ((residual * residual)[:, None] * matrix)
    covariance = inverse @ meat @ inverse * n / max(n - p, 1)
    return coefficient, np.sqrt(np.diag(covariance))


def design_matrix(stars, with_step: bool = False):
    """Temperature plus the things that could fake a temperature trend: period, amplitude,
    brightness and number of sectors. Optionally a step for stars below 3,400 K."""
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


def period_slope(stars, floor: str) -> dict:
    """Change in log10(excess) per doubling of the rotation period, for stars with a positive
    excess, holding amplitude, brightness and sector count fixed."""
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


def load_qlp_sector(path: Path):
    """A QLP light curve from its SAP flux: quality flags applied, a quadratic removed from each
    TESS orbit, 5-sigma outliers clipped. Same return values as load_sector.

    QLP's KSPSAP_FLUX is flattened for transit searches and takes the rotation signal out, so it
    cannot be used here.
    """
    with fits.open(path, memmap=False) as hdus:
        sector = int(hdus[0].header["SECTOR"])
        data = hdus[1].data
        names = set(hdus[1].columns.names)
        time = np.asarray(data["TIME"], float)
        flux = np.asarray(data["SAP_FLUX"], float)
        quality = np.asarray(data["QUALITY"], int)
        orbit = np.asarray(data["ORBITID"], int) if "ORBITID" in names else np.zeros(len(time), int)
    keep = np.isfinite(time) & np.isfinite(flux) & (quality == 0) & (flux > 0)
    time, flux, orbit = time[keep], flux[keep], orbit[keep]
    if len(time) < MINIMUM_POINTS:
        return None
    normalised = np.empty_like(flux)
    for value in np.unique(orbit):
        mask = orbit == value
        if mask.sum() < 50:
            normalised[mask] = flux[mask] / np.median(flux[mask])
            continue
        x = time[mask] - time[mask].mean()
        normalised[mask] = flux[mask] / np.polyval(np.polyfit(x, flux[mask], 2), x)
    flux = normalised
    residual = flux - np.median(flux)
    sigma = 1.4826 * np.median(np.abs(residual))
    keep = np.abs(residual) < 5 * max(sigma, 1e-6)
    time, flux = time[keep], flux[keep]
    if len(time) < MINIMUM_POINTS:
        return None
    noise = float(1.4826 * np.median(np.abs(np.diff(flux))) / np.sqrt(2))
    return sector, time, flux, max(noise, 1e-5)


def compare_pipelines(row: dict, tars_folder: Path = TARS_DIR, qlp_folder: Path = QLP_DIR) -> list[dict]:
    """For one star, fit the same sectors in TARS and in QLP. Returns one row per shared sector
    with each pipeline's deviation from the star's median angular speed."""
    wanted, p0 = row["overlap_sectors"], row["period"]
    periods = {"tars": {}, "qlp": {}}
    for name, folder, pattern, loader in (
        ("qlp", qlp_folder, "_llc.fits", load_qlp_sector),
        ("tars", tars_folder, "_lc.fits", load_sector),
    ):
        for path in sorted(Path(folder).glob(f"*-{row['tic_id']:016d}_*{pattern}")):
            loaded = loader(path)
            if loaded and loaded[0] in wanted:
                sector, time, flux, noise = loaded
                fitted = fit_period(time, flux, noise, p0)
                if fitted:
                    periods[name][sector] = fitted[0]
    shared = sorted(set(periods["tars"]) & set(periods["qlp"]))
    if len(shared) < MINIMUM_SECTORS:
        return []
    out = {}
    for name in ("tars", "qlp"):
        omega = 2 * np.pi / np.array([periods[name][s] for s in shared])
        out[name] = omega - np.median(omega)
    return [dict(tic_id=row["tic_id"], teff=row["teff"], sector=s, dev_tars=float(a), dev_qlp=float(b))
            for s, a, b in zip(shared, out["tars"], out["qlp"])]
