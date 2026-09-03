"""A small forward model: fake spotted stars where WE choose the cause of the wander.

Uses Butterpy (Claytor et al. 2022) to grow and fade spots on a rotating star, samples the
result at a real star's own TESS timestamps, adds that star's own leftover noise, and measures
it with the identical fit_period. Two knobs are turned:

  shear_on         False -> rigid rotation;  True -> Collier Cameron (2007) shear for that temperature
  lifetime_factor  1 -> spot lifetimes from the Giles et al. (2017) Kepler relation; 2, 3 -> longer

Same recipe as research/src/teff_stability_forward_model.py, reduced to a small grid.
"""
from __future__ import annotations

import numpy as np

import helpers as h

SIM_CADENCE_DAYS = 1 / 48


def giles_lifetime_days(rms: float, teff: float) -> float:
    """Giles, Collier Cameron & Haywood 2017, eq. 8: spot lifetime from light-curve rms and temperature."""
    x = np.log10(max(rms, 1e-4))
    return float(10 ** (10.9252 + 3.0123 * x + 0.5062 * x**2 - 1.3606 * np.log10(teff)))


def collier_cameron_shear(teff: float) -> float:
    """Collier Cameron 2007: equator-to-pole shear in rad/day, 0.053 (Teff/5130)^8.6."""
    return 0.053 * (teff / 5130.0) ** 8.6


def real_sectors(row: dict, cache) -> list[dict]:
    """The real star's accepted sectors: timestamps, noise, fitted amplitude and leftover."""
    sectors = []
    for path in h.sector_files(cache, row["tic_id"]):
        loaded = h.load_sector(path)
        if loaded is None:
            continue
        sector, time, flux, noise = loaded
        if sector not in row["sectors_allowed"]:
            continue
        fitted = h.fit_period(time, flux, noise, row["period"])
        if fitted is None:
            continue
        sectors.append(dict(sector=sector, time=time, noise=noise, amplitude=fitted[1], leftover=fitted[2]))
    return sectors


def simulate(row: dict, sectors: list[dict], shear_on: bool, lifetime_factor: float, seed: int) -> dict | None:
    from butterpy import Surface

    rng = np.random.default_rng(seed)
    p0, teff = row["period"], row["teff"]
    observed_amplitude = float(np.median([s["amplitude"] for s in sectors]))
    lifetime_days = lifetime_factor * giles_lifetime_days(observed_amplitude / np.sqrt(2), teff)
    if shear_on:
        shear = collier_cameron_shear(teff) * p0 / (2 * np.pi)   # Butterpy wants shear relative to the rotation rate
        min_lat, max_lat = 30, 60
    else:
        shear = 0.0
        min_lat, max_lat = 7, 28
    start = min(float(s["time"].min()) for s in sectors) - 2 * lifetime_days
    stop = max(float(s["time"].max()) for s in sectors) + 1
    duration = stop - start
    np.random.seed(int(rng.integers(0, 2**31 - 1)))          # Butterpy uses numpy's global random state
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

    # Sample the model at the real timestamps; scale it to the real star's amplitude.
    modelled = [np.interp(s["time"] - start, sim_time, sim_flux) for s in sectors]
    clean = [h.fit_period(s["time"], v, s["noise"], p0) for s, v in zip(sectors, modelled)]
    clean_amplitudes = [c[1] for c in clean if c is not None]
    if len(clean_amplitudes) < 3 or np.median(clean_amplitudes) <= 0:
        return None
    scale = observed_amplitude / float(np.median(clean_amplitudes))

    # Add the real leftover (shuffled) and measure exactly as the real star was measured.
    fitted_sectors, periods = [], []
    for s, v in zip(sectors, modelled):
        fitted = h.fit_period(s["time"], scale * v + rng.permutation(s["leftover"]), s["noise"], p0)
        if fitted is None:
            continue
        periods.append(fitted[0])
        fitted_sectors.append(dict(time=s["time"], noise=s["noise"], amplitude=fitted[1], leftover=fitted[2]))
    if len(periods) < h.MINIMUM_SECTORS:
        return None
    wander = h.angular_mad(periods)

    # The same white-noise floor the real stars got (a steady wave in the simulated star's leftover).
    scatters = []
    for _ in range(h.INJECTIONS):
        phase = rng.uniform(0, 2 * np.pi)
        fake = []
        for s in fitted_sectors:
            wave = s["amplitude"] * np.sin(2 * np.pi * s["time"] / p0 + phase)
            fitted = h.fit_period(s["time"], wave + rng.permutation(s["leftover"]), s["noise"], p0)
            fake.append(fitted[0] if fitted else np.nan)
        fake = np.array(fake)
        if np.isfinite(fake).sum() >= 3:
            scatters.append(h.angular_mad(fake[np.isfinite(fake)]))
    floor = float(np.median(scatters)) if scatters else np.nan
    return dict(
        tic_id=row["tic_id"], teff=teff, period=p0, shear_on=shear_on, lifetime_factor=lifetime_factor,
        lifetime_days=lifetime_days, sectors=len(periods), simulated_wander=wander, simulated_floor=floor,
        simulated_excess=float(np.sqrt(max(wander**2 - floor**2, 0.0))),
        simulated_ratio=wander / floor if floor > 0 else np.nan,
    )


def run_config(args):
    """One (star, shear_on, lifetime_factor) job, for a process pool."""
    row, cache, shear_on, lifetime_factor = args
    sectors = real_sectors(row, cache)
    if len(sectors) < h.MINIMUM_SECTORS:
        return None
    seed = 20261210 + row["tic_id"] + (1000 if shear_on else 0) + int(10 * lifetime_factor)
    return simulate(row, sectors, shear_on, lifetime_factor, seed)
