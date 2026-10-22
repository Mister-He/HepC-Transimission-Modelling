# Training

Production calibration path:
1. **Nelder–Mead multi-start** finds the posterior mode (fixed prior + seroprevalence likelihood).
2. **Reference MCMC** gives the production posterior.
3. **Scenarios** run on draws from that posterior.

Cross-validation only (nothing here feeds the production path): single-round NPE trained on
50,000 Latin hypercube prior simulations, with simulation-based calibration (SBC), and a small
direct MCMC SBC.

Outputs are in `output/calibration/` (see its `README.md`). Run the steps in this order; each step
reads the previous step's output. Model design: `docs/first_arrest_multiplier.md`.

## 0. Rate fits (only after refitting the Stan rate models)

```bash
Rscript scripts/model/preprocess.R            # needs data/sps.csv (not tracked by git); ~5 min
Rscript scripts/model/update_rate_params.R    # writes lambda2, lambda3 and pi_recid into setup.R
```

- **Data:** everyone ever held for drug consumption. Screening (Dec 2014–Feb 2016) selects people
  in prison then, so in the screened subset every earlier release ends in re-arrest; the rates
  therefore use all consumption records.
- **`release_rate.stan`:** age-specific release hazard (Lexis-split stays). Both fits are saved as
  `list(summary, draws)`; `lambda2` and `lambda3` are the posterior means of `exp(log_lambda)`.
- **Also writes** `output/prison_2015_counts.csv`: the 2015 prison count of screened PWID and the
  screened seropositives by age group (aggregates only; not calibration targets).
- **`rearrest_rate.stan`:** cure-type survival model only (share `pi` who are ever re-arrested and
  age-specific re-arrest hazard `lambda3`). There is no prison steady-state block and no
  first-arrest flow. Fit 2026-10-22: `pi` 0.724, `lambda3` 0.548, 0.478, 0.469, 0.471, 0.384,
  0.265 per year, R-hat 1.00.
- **First arrest** is not estimated from the records (they cannot identify true first offenders):
  `lambda1 = c_first × lambda3`, with `c_first` calibrated (θ7).
- **Idempotent:** on the stored fits `update_rate_params.R` leaves `setup.R` unchanged (checked
  by `tests/unit/test_rate_params.R`). After a change, rerun every step below.

## Shared posterior target

- **Parameters (7):** θ1–θ6 are log contact scales (relative to `C_contact` in `setup.R`); θ7 is
  log `c`, the first-arrest multiplier. Entry into D is the CNB new-drug-abuser count, unscaled
  (`build_params()` in `scripts/calibration/likelihood.R`).
- **Bounds:** contact scales 0.01–1,000; `c` 0.5–100 (user decision: first arrest at least half
  the re-arrest rate).
- **Prior (fixed, not centred on any fit):** independent Normal, truncated to the bounds:
  θ1–θ6 ~ N(0, 1.5), θ7 ~ N(0, 0.5) (`prior_spec()`; density `nlp()`, shared by Nelder–Mead,
  MCMC and NPE).
- **Likelihood:** seroprevalence only, no counts. Each age group's observed seroprevalence is
  `Beta(p·φ, (1−p)·φ)` around the model's 2015 prison seroprevalence `p`, with fixed precision
  φ = 100 (`beta_precision` in `data/calibration_targets.json`). Prison population size is not a
  target.
- **Support:** draws outside the bounds and equilibrium failures are excluded (MCMC gives −∞;
  NPE drops the simulation and counts it).
- **Solver:** RK4 with `dt = 1/365`.
- **Mortality:** general-population `mu` × PWID SMR by state: 9.52 viraemic (a, c, t), 4.73
  otherwise (van Santen 2014). The viraemic share varies over time (2015 best fit: 0.184).
- **Details:** see `output/calibration/npe_vs_mcmc/model_audit.md`.

## 1. Deterministic calibration (posterior mode)

```bash
Rscript scripts/run_calibration.R --run-id r_20261021_agepop \
  --maxit 3000 --nm-restarts 5 --n-starts 8 --sd-perturb 0.9 --seed 20261001 --cores 8   # ~2 h
```

- **Objective:** negative log posterior (Beta likelihood + prior); every evaluation runs model
  time 0–150 and checks equilibrium over the last 10 years.
- **Starts (uninformative):** θ = 0 (contact scales 1, `c` = 1) and 7 Normal perturbations of it
  (SD 0.9), clipped to the bounds, run in parallel.
- **Restarts:** each start reruns Nelder–Mead from its own solution while the objective still
  falls by more than 0.001 (`--nm-restarts`, at most 5 restarts).
- **Result (2026-10-22):** all 8 starts reach the same mode (objective −7.83); seroprevalence
  RMSE 0.013, largest error 0.023 (50–59); `c` 1.10.
- **SMR sensitivity** (`output/calibration/smr_sensitivity/`, 3 starts each,
  `--smr-viraemic`/`--smr-nonviraemic` at the CI bounds): `c` 1.097–1.100, RMSE 0.012–0.015.

## 2. Reference MCMC (production)

```bash
Rscript scripts/run_mcmc_reference.R --cores 4   # about 4 hours
```

- **Defaults:**
  - 4 chains from independent prior draws;
  - warmup 10,000, with adaptation from proposal 5,000;
  - 50,000 evaluations per chain;
  - draws stored every 10th iteration (`--thin`);
  - seed 20261101;
  - fit `output/calibration/r_20261021_agepop/fit.rds` (used only to check that the model matches);
  - output to `output/calibration/mcmc_reference/`.
- **Sampler:** adaptive random-walk Metropolis (`arwm_budget` in `scripts/calibration/mcmc.R`). The
  initial proposal is `0.1² × diag(prior variance)`; adaptation uses the recent half of warmup
  only. Starts are truncated-prior draws.
- **Result (2026-10-22, 4.0 hours):**
  - R-hat ≤ 1.002;
  - bulk ESS 3,250–4,645 and tail ESS ≥ 3,487;
  - MCSE of the mean ≤ 0.019 posterior SD;
  - acceptance 0.18–0.22 (0.08–0.17 before the recent-half adaptation fix).
- **Prior-dominated parameters:** posterior SD above 0.7 × prior SD: θ5, θ6 (50–59 and 60+ contact
  scales) and θ7 (`c`; posterior median 1.24, 95% interval 0.57–3.13, SD ratio 0.87).
- `--targets` swaps the target file. The MCMC SBC uses it.

## 3. Cross-validation: single-round NPE with SBC

```bash
Rscript scripts/run_npe.R --step all --n-sims 50000 --seed 20261012 --train-seed 20261013 \
  --n-sbc 1000 --n-cores 10 --threads 4 --python /path/to/python
```

- **LHS design** (`scripts/calibration/prior_sampling.R`): a 50,000 × 7 design with one point per
  stratum, mapped through the truncated Normal inverse CDF of each parameter's prior.
- **Simulations:** draws outside the bounds, non-finite results and equilibrium failures are dropped
  and counted, never replaced. This run: 49,945 valid, 43 equilibrium failures and 12 non-finite
  (`simulation_summary.json`).
- **Summaries:** 6 seroprevalences, with Beta(φ = 100) observation noise, logit-transformed.
- **Network:** an equal-weight ensemble of 5 MAFs (50×5), each with batch 256, Adam lr 5e-4,
  validation fraction 0.10, patience 20, ≤ 250 epochs, gradient clipping 5.0, CPU only. Posterior
  draws outside the bounds are rejected (at most 100 rounds, then an error). Training takes
  about 6 minutes.
- **Output:** 20,000 posterior draws.
- **SBC:** 1,000 independent prior datasets (998 pass equilibrium), with ranks of the true value in
  1,000 posterior draws.
  - **Gate:** cov95 in [0.85, 1] and cov90 in [0.75, 1] for every parameter, and no systematically
    skewed rank histogram.
  - **Result:** cov95 0.933–0.960, cov90 0.872–0.911, all 7 rank shapes uniform; gate PASS
    (`NPE_OK.txt`).
- **Comparison with the reference MCMC:** NPE draws that fail the equilibrium check are dropped
  first (1,329 of 4,000, 33%). At the observed data the NPE does not match the MCMC (0 of 7
  parameters pass all checks; 95% intervals on average 28% wider; worst mean difference 0.56 SD).
  Policy reductions agree within 0.3%; joint-pass probability at 5,000/yr F25% differs by 8 pp and
  cost by 6%. See
  `output/calibration/npe_vs_mcmc/report.md`.

  ```bash
  Rscript scripts/compare_npe_mcmc.R --policy-cores 14
  ```

## 4. Optional: direct MCMC SBC

```bash
Rscript scripts/sbc_mcmc.R --n-datasets 8 --max-evals 30000 --cores 8   # about 2.5 hours
```

- The truths are the first 8 NPE SBC datasets, with synthetic seroprevalences drawn from the Beta
  likelihood. Each dataset gets one 30,000-evaluation chain.
- **Result:** pooled cov95 0.96 and cov90 0.89 over 56 ranks; mean rank 0.48; max split R-hat 1.008.
- **Limitation:** 8 datasets are a check against gross miscalibration only, not a full SBC.

## 5. 2015 predictions and calibration figure

```bash
Rscript scripts/predict_2015_posterior.R --max-draws 4000 --cores 14
Rscript scripts/plot_calibration.R --cores 14 --max-draws 2000
```

- **Prediction outputs:** `output/calibration/mcmc_reference/model_prediction_2015_strata.*`.
  Ever-PWID 50,218 (95% interval: 49,496–50,893); D 2,178 (839–4,651); J 8,451 (8,283–8,520).
  These sizes are not fitted: the observed 2015 prison count of screened PWID is 3,367
  (`output/prison_2015_counts.csv`).
- **Calibration figure:** `figures/calibration_agepop.png` compares the seroprevalence targets
  (with the 95% range of the Beta likelihood) with the reference MCMC posterior and the NPE (LHS).
- **Overlap check:** `figures/calibration_interval_overlap.csv` records whether each 95% interval
  overlaps its target (MCMC: 6 of 6).

## 6. 2026 treatment scenarios

```bash
Rscript scripts/scenarios_2026.R --root . --draws-per-seed 2000 --cores 4 --sample-seed 20261005
```

- **Draws:** 2,000 draws from the reference MCMC posterior (the default `--posterior-files`). The full
  analysis is run per draw, then summarised as the median and 95% UI.
- **Options:**
  - `--posterior-files` accepts any CSV with `theta1..theta7`;
  - `--figures true` also draws the old ggplot figures (they overwrite the `figures.py` PNGs);
  - `--min-valid-share` sets the share of draws that must pass equilibrium.
- **More detail:** see `output/scenarios_2026/README.md`.

## 7. Manuscript

```bash
python scripts/figures.py .   # redraws every figure (PNG + editable SVG)
python3 scripts/build_paper.py .
```

`figures.py` redraws all figures from the saved outputs with the `scientific-visualization` skill
(`~/.claude/skills/scientific-visualization`, or pass `--skill-dir`). It needs Python ≥ 3.11 with
`matplotlib==3.11.2 pillow==12.3.0 pandas scipy`. `--only <name> ...` redraws selected figures.
Run it after any step above that changes its inputs. The R scripts still draw their own ggplot
versions of the same files, so run `figures.py` last.
