"""Fake spotted stars, measured exactly like the real ones.

Butterpy (Claytor et al. 2022) grows and decays spots on a rotating star. Each fake star copies a
real star from the sample: its period, its TESS time stamps, its amplitude and its own leftover
noise. Spot lifetimes come from the Giles, Collier Cameron & Haywood (2017) relation. Two versions:

    rigid   no differential rotation; spots emerge at 7-28 degrees latitude (like the Sun)
    shear   Collier Cameron (2007) shear for the star's temperature; spots at 30-60 degrees

The fake light curve is then measured with the same fit_period and white-noise floor as the real
star, so the model and the data go through identical code. Nothing is tuned to match the data.
"""
from __future__ import annotations

import numpy as np

import wander as w

STARS_PER_BIN = 10          # the first 10 stars (by TIC number) in each 200 K slice of the original sample
SIM_CADENCE_DAYS = 1 / 48
SEED = 20261210


def giles_lifetime_days(rms: float, teff: float) -> float:
    """Giles, Collier Cameron & Haywood 2017, eq. 8: spot lifetime in days from light-curve rms."""
    x = np.log10(max(rms, 1e-4))
    return float(10 ** (10.9252 + 3.0123 * x + 0.5062 * x**2 - 1.3606 * np.log10(teff)))


def collier_cameron_shear(teff: float) -> float:
    """Collier Cameron 2007: equator-to-pole shear in rad/day, 0.053 (Teff/5130)^8.6."""
    return 0.053 * (teff / 5130.0) ** 8.6


def model_stars(rows: list[dict]) -> list[dict]:
    chosen = []
    for low in range(3000, 6000, 200):
        in_bin = sorted((r for r in rows if low <= r["teff"] < low + 200), key=lambda r: r["tic_id"])
        chosen.extend(in_bin[:STARS_PER_BIN])
    return chosen


def simulate(row: dict, sectors: list[dict], shear_on: bool, rng: np.random.Generator) -> dict | None:
    from butterpy import Surface

    p0, teff = row["period"], row["teff"]
    observed_amplitude = float(np.median([s["amplitude"] for s in sectors]))
    lifetime_days = giles_lifetime_days(observed_amplitude / np.sqrt(2), teff)
    if shear_on:
        shear = collier_cameron_shear(teff) * p0 / (2 * np.pi)   # Butterpy wants shear relative to the rotation rate
        min_lat, max_lat = 30, 60
    else:
        shear = 0.0
        min_lat, max_lat = 7, 28
    start = min(float(s["time"].min()) for s in sectors) - 2 * lifetime_days
    stop = max(float(s["time"].max()) for s in sectors) + 1
    duration = stop - start
    np.random.seed(int(rng.integers(0, 2**31 - 1)))          # Butterpy draws from numpy's global random state
    surface = Surface()
    surface.emerge_regions(ndays=int(duration) + 5, activity_level=1, butterfly=True, cycle_period=11,
                           cycle_overlap=2, max_lat=max_lat, min_lat=min_lat)
    sim_time = np.arange(0, duration, SIM_CADENCE_DAYS)
    curve = surface.evolve_spots(time=sim_time, inclination=90, period=p0, shear=shear,
                                 tau_evol=lifetime_days / p0, alpha_med=0.005)
    sim_flux = np.asarray(curve.flux, float)
    sim_flux = sim_flux / np.nanmedian(sim_flux) - 1.0
    if not np.isfinite(sim_flux).all() or np.nanstd(sim_flux) == 0:
        return None

    # Sample the model at the real time stamps and scale it to the real star's amplitude.
    modelled = [np.interp(s["time"] - start, sim_time, sim_flux) for s in sectors]
    clean = [w.fit_period(s["time"], 1.0 + v, s["noise"], p0) for s, v in zip(sectors, modelled)]
    clean_amplitudes = [c[1] for c in clean if c is not None]
    if len(clean_amplitudes) < 3 or np.median(clean_amplitudes) <= 0:
        return None
    scale = observed_amplitude / float(np.median(clean_amplitudes))

    # Add the real star's leftover (shuffled) and measure it like a real star.
    periods, refit = [], []
    for s, v in zip(sectors, modelled):
        fitted = w.fit_period(s["time"], 1.0 + scale * v + rng.permutation(s["leftover"]), s["noise"], p0)
        if fitted is None:
            continue
        periods.append(fitted[0])
        refit.append(dict(time=s["time"], noise=s["noise"], amplitude=fitted[1], leftover=fitted[2]))
    if len(periods) < w.MINIMUM_SECTORS:
        return None
    wander = w.angular_mad(periods)

    scatters = []
    for _ in range(w.INJECTIONS):
        phase = rng.uniform(0, 2 * np.pi)
        fake = []
        for s in refit:
            wave = s["amplitude"] * np.sin(2 * np.pi * s["time"] / p0 + phase)
            fitted = w.fit_period(s["time"], wave + rng.permutation(s["leftover"]), s["noise"], p0)
            fake.append(fitted[0] if fitted else np.nan)
        fake = np.array(fake)
        if np.isfinite(fake).sum() >= 3:
            scatters.append(w.angular_mad(fake[np.isfinite(fake)]))
    floor = float(np.median(scatters)) if scatters else np.nan
    return dict(
        tic_id=row["tic_id"], teff=teff, period=p0, model="shear" if shear_on else "rigid",
        sectors=len(periods), amplitude=observed_amplitude, spot_lifetime_days=lifetime_days,
        shear_rad_per_day=collier_cameron_shear(teff) if shear_on else 0.0,
        simulated_wander=wander, simulated_floor=floor,
        simulated_excess=float(np.sqrt(max(wander**2 - floor**2, 0.0))),
    )


def run_one(job) -> dict | None:
    """job = (star row, shear_on). For a process pool."""
    row, shear_on = job
    sectors = w.fit_sectors(row, w.TARS_DIR)
    if len(sectors) < w.MINIMUM_SECTORS:
        return None
    rng = np.random.default_rng(SEED + row["tic_id"] + (1 if shear_on else 0))
    return simulate(row, sectors, shear_on, rng)
