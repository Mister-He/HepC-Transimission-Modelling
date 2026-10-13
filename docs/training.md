# Training

Current calibration workflow:
1. **Nelder–Mead** sets the prior centre.
2. **Reference MCMC** gives the production posterior.
3. **Single-round NPE**, trained on 50,000 Latin hypercube prior simulations, is the independent
   check.

Outputs are in `output/calibration/` (see its `README.md`).

## Deterministic calibration (prior centre)

```bash
Rscript scripts/run_calibration.R --run-id r_20261010_agepop \
  --maxit 3000 --n-starts 7 --sd-perturb 0.9 --seed 20261001 --cores 8
```

- **Equilibrium:** every objective evaluation runs model time 0–150 and checks equilibrium between
  years 140 and 150.
- **Starts:** the anchor, 6 perturbed starts and a population-informed start, run in parallel.
- **Result:**
  - best objective 19.01 (`perturb_1`, converged);
  - seroprevalence RMSE 0.0044 and population MAPE 2.3%.

## Shared posterior target

- **Prior:** independent Normal(`fit$best$theta`, 0.20) on the 12 log-scale parameters.
- **Likelihood:** Binomial prevalence counts and log-normal populations (σ = 0.10).
- **Equilibrium:** failure is a hard restriction (MCMC gives −∞; NPE drops the simulation).
- **Solver:** RK4 with `dt = 1/365`.
- **Details:** see `output/calibration/npe_vs_mcmc/model_audit.md`.

## Reference MCMC (production)

```bash
Rscript scripts/run_mcmc_reference.R --fit output/calibration/r_20261010_agepop/fit.rds   # about 3.8 hours
```

- **Defaults:** 4 chains from independent prior draws; warmup 10,000 with adaptation from proposal
  5,000; 50,000 evaluations per chain; draws stored every 10th iteration (`--thin`); seed 20261101;
  output to `output/calibration/mcmc_reference/`.
- **Reproducing the stored production run:** pass `--n-warmup 4999 --adapt-start 200 --thin 1`.
- **Ablation options:** `--start-file`, `--cov-file` with `--cov-scale`, `--dt` and `--adapt-every`.
- **Diagnostics written:** trace plots, per-chain summaries, adaptation drift, covariance eigenvalues,
  jump sizes and autocorrelation times.
- **Sampler:** adaptive random-walk Metropolis (`arwm_budget` in `scripts/calibration/mcmc.R`).
- **Result:** R-hat ≤ 1.002, bulk ESS 3,433–4,600, and MCSE of the mean ≤ 0.017 posterior SD.
- **Provenance:** the stored run was produced on 2026-10-10 under the script's earlier name and moved
  here unchanged.

## Reference-MCMC ablation

```bash
nohup Rscript scripts/ablation_check.R --run-id ab_20261012 > output/ablation_check/ab_20261012/driver.log 2>&1 &
```

- **Runtime:** about 6 hours on 16 cores.
- **Arms:** six covariance and initialisation arms, plus a daily-versus-monthly step check.
- **Main findings** (`output/ablation_check/ab_20261012/ablation_report.md`):
  - every arm converges (R-hat ≤ 1.005, bulk ESS ≥ 1,422) and agrees with the production posterior
    within Monte Carlo error;
  - starts do not matter;
  - posterior-shaped proposals (NPE, pilot, Hessian) give about 5–30% more ESS than the prior-scaled
    proposal;
  - the predeclared adaptation-drift tolerance (0.01) was not met, although adaptation is frozen
    before retention;
  - thinning by 10 is enough;
  - the monthly step shifts the posterior far less than its width.

## Single-round NPE (50,000 LHS prior simulations)

```bash
Rscript scripts/run_npe.R --step all --fit output/calibration/r_20261010_agepop/fit.rds \
  --out-dir output/calibration/npe_lhs --n-sims 50000 --sampling lhs --seed 20261012 \
  --train-seed 20261013 --n-cores 14 --threads 4 --python /path/to/python
```

- **LHS design** (`scripts/calibration/prior_sampling.R`): a 50,000 × 12 design with one point per
  stratum, mapped through the Normal inverse CDF of each parameter's prior.
- **Simulations:** non-finite results and equilibrium failures are dropped and counted, never
  replaced. See `simulation_summary.json` and `prior_design.csv`.
- **Network:** MAF 50×5, batch 256, Adam lr 5e-4, validation fraction 0.10, patience 20,
  ≤ 250 epochs, gradient clipping 5.0, CPU only.
- **Output:** 20,000 posterior draws.
- **Comparison with the reference MCMC:**

  ```bash
  Rscript scripts/compare_npe_mcmc.R --policy-cores 14
  ```

## 2015 predictions and calibration figure

```bash
Rscript scripts/predict_2015_posterior.R --fit output/calibration/r_20261010_agepop/fit.rds \
  --max-draws 4000 --cores 14
Rscript scripts/plot_calibration.R --cores 14 --max-draws 2000
```

- **Prediction outputs:** `output/calibration/mcmc_reference/model_prediction_2015_strata.*`. The
  ever-PWID estimate is 21,224 (95% interval: 19,226–23,437).
- **Calibration figure:** `figures/calibration_agepop.png` compares the targets with the
  reference MCMC posterior and single-round NPE (LHS).
- **Overlap check:** `figures/calibration_interval_overlap.csv` shows that all 12 intervals
  overlap the targets.

## 2026 treatment scenarios

```bash
Rscript scripts/scenarios_2026.R --root . --fit output/calibration/r_20261010_agepop/fit.rds \
  --draws-per-seed 2000 --cores 14 --sample-seed 20261005
```

- **Draws:** 2,000 draws from the reference MCMC posterior (the default `--posterior-files`). The full analysis is run per
  draw, then summarised as the median and 95% UI.
- **Options:**
  - `--posterior-files` accepts any CSV with `theta1..theta12`;
  - `--no-figures true` skips the figures;
  - `--min-valid-share` sets the share of draws that must pass equilibrium.
- **More detail:** see `output/scenarios_2026/README.md`.
