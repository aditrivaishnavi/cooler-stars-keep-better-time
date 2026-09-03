# Cooler stars keep better time — the simplified version

*Do spotted stars' rotation readings wander less from month to month when the star is cooler? A 450-star TESS measurement, replicated on 422 more, done in six notebooks and two short Python files.*

This folder redoes the rotation-stability project step by step: one notebook per step,
plain `numpy`/`astropy` code, every number produced in front of you. It uses **exactly the same
fitting choices** as the full research pipeline in `../research/`, so the numbers are the
same — `results/reconciliation.csv` shows 44 of 44 headline quantities agreeing to the fourth decimal
(the per-star tables agree to 1e-13); every number in the v7 conference abstract (`../paper/submission/ABSTRACT.md`) comes out of notebooks 3–5. What differs is everything around the numbers: 355 lines of code
instead of 14,000, and six notebooks instead of 67 scripts.

The full pipeline is the record for the journal paper.

## What is here

| | Notebook | What it does | Time |
|---|---|---|---|
| 1 | `01_pick_stars` | Loads the two frozen star lists (450 + 422), shows the selection rule and that the slices share brightness and period windows | seconds |
| 2 | `02_get_data` | Downloads TARS light curves from MAST with `astroquery` (demonstrated live on one star), checks every star's sectors are on disk, plots one sector | seconds if cached; ~1 h to download ~900 MB |
| 3 | `03_measure` | **The method.** One period per sector; the wander; the fake steady star, two scrambles; the floor and the excess. Runs on all 872 stars | ~15 min on 16 cores |
| 4 | `04_results` | Bin medians, the fourfold fall, the 0.015 rad/day bar, the covariate regression, the step test, period dependence, early/late check, the replication against its pre-registered rules, **and the reconciliation table against the frozen abstract numbers** | seconds |
| 5 | `05_spot_model` | A small Butterpy forward model: 20 M-dwarf and 20 Sun-like twins × shear off/on × spot lifetime 1×/2×/3× | ~7 min |
| 6 | `06_independent_pipeline` | The same sectors from MIT's QLP pipeline for 20 stars: do the two pipelines see the same sectors deviate? | seconds if cached |

Supporting files:

- `helpers.py` (230 lines) — `load_sector`, `fit_period`, `angular_mad`, `measure_star`, `robust_fit`, `design_matrix`, `period_slope`. The whole method.
- `spotmodel.py` (125 lines) — the Butterpy forward model.
- `data/star_lists/` — the two frozen samples (drawn 1 Sep 2026, before any result).
- `data/frozen_numbers.json` — the abstract's numbers as produced by the full pipeline, for the reconciliation.
- `data/lightcurves/` — `original/`, `replication/`, `qlp/` (here: links to the full project's caches; notebook 2 fills them from MAST if they are empty), `demo/` (the live-download demonstration).
- `results/` — `stars_original.csv`, `stars_replication.csv` (one row per star), `reconciliation.csv`, `summary.json`, `spot_model.csv`, `qlp_sector_deviations.csv`.
- `figures/` — the five notebook figures.

## How to run it

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter lab notebooks/          # run 01 → 06 in order
```

Notebook 3 skips the 15-minute computation if `results/stars_*.csv` already exist; delete them to recompute.
Notebook 5 does the same with `results/spot_model.csv`.

## The result, as the notebooks produce it

Measured (notebooks 3–4), identical to the frozen abstract:

| | original 450 | replication 422 |
|---|---|---|
| Stars / sectors | 450 / 4,317 | 422 / 4,016 |
| Excess wander, hot slices (≥ 5,000 K) vs cool (< 3,800 K), adopted floor | 0.019–0.050 vs 0.004–0.010 rad/day | 0.023–0.034 vs 0.006–0.008 |
| Share of stars above the fixed bar 0.015 rad/day, hot ≥ 5,000 K vs cool < 3,400 K | **83% vs 18%** (white floor 88% vs 22%) | 78% vs 15% (88% vs 25%) |
| Temperature coefficient of log excess ratio, adopted floor | **+0.127 ± 0.021 dex / 1,000 K (6.1σ)** | +0.109 ± 0.023 (4.8σ) |
| Same, white floor | +0.279 ± 0.029 (9.7σ) | +0.244 ± 0.030 (8.1σ) |
| Step below 3,400 K after the trend (t) | −0.1 adopted / +0.7 white | +1.0 / +1.2 |
| Period slope, cool 3,000–3,800 K (dex per doubling, excess in quadrature) | +0.19 ± 0.06 / +0.18 ± 0.04 | +0.10 ± 0.05 / +0.18 ± 0.05 |
| Period slope, hot 4,800–6,000 K | −0.03 ± 0.05 / +0.04 ± 0.04 | +0.02 ± 0.06 / +0.08 ± 0.05 |
| Early vs late mission ratio vs temperature (Spearman) | +0.04, p = 0.44 | — |
| Replication's three pre-registered predictions | — | all three met |

Model (notebook 5; small grid, so ±30% on any entry): a rigid Sun-like model already matches the Sun-like twins (model/observed 0.95); switching Collier Cameron shear on moves it by a third. The M-dwarf twins are 3–4× steadier than the Giles-lifetime model, shear on or off; tripling the spot lifetime halves the gap. The full pipeline's 3,888-run grid closes it at ≈3× the Giles relation.

Independent pipeline (notebook 6): 20 stars, 208 shared sectors, pooled Spearman ρ = 0.62, 100% of stars positive (full pipeline, 60 stars: ρ = 0.60, 98%).

## What is simplified, what is left out

Kept exactly: the star lists; the per-sector fit (linear detrend, Lomb–Scargle in 0.77–1.43 × the catalogue period, ±3% refinement with parabolic peak); the wander; both floors with 12 fakes each and the same random seeds; the regression design and robust errors; the fixed bar; the step test; the period slopes; the early/late check; the replication rules.

Simplified: data download (one `astroquery` cell); the forward model (240 runs instead of 3,888; lifetime is the only knob besides shear); the figures (five); verification (one reconciliation cell instead of 38 tests); provenance (CSV + JSON per notebook instead of SHA-256 manifests).

Left out, with where it lives: the Doppler-imaging comparison (a literature table; `../paper/PAPER_DRAFT.md` Appendix A), the APOGEE/LAMOST temperature check, the boundary curvature and placebo scans, the 48-injection and P ≤ 5 d follow-ups, the classifier-free acceptance test, the shear sweep and full morphology grid, and all exploratory phases (`../research/EXPLORATION_LOG.md`). None of them changes a headline number.

## Honesty notes

- The numbers agree with the full pipeline because the fitting choices and random seeds are identical, not because the two were tuned to agree. `results/reconciliation.csv` is the evidence.
- The "step" t-values quoted in the paper draft (0.7 and 1.2) are the white-floor values; under the adopted floor they are −0.1 and +1.0. Both say the same thing: no detectable step.
- The forward model here turns two knobs; the conclusion "spot evolution, not shear, sets the wander" is reproduced, the exact lifetime factor is not. The full grid is the reference for that.
- What this measurement cannot say: that M dwarfs shear less. The model shows the wander barely responds to shear at realistic levels.
