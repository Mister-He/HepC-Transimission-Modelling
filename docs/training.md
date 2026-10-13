# Training

Production calibration path:
1. **Nelder–Mead multi-start** sets the prior centre.
2. **Reference MCMC** gives the production posterior.
3. **Scenarios** run on draws from that posterior.

Cross-validation only (nothing here feeds the production path): single-round NPE trained on
50,000 Latin hypercube prior simulations, with simulation-based calibration (SBC), and a small
direct MCMC SBC.

Outputs are in `output/calibration/` (see its `README.md`). Run the steps in this order; each step
reads the previous step's output.

## 0. Rate-fit handoff (only after refitting the Stan rate models)

```bash
Rscript scripts/model/preprocess.R            # needs data/sps.csv (not in the repository)
Rscript scripts/model/update_rate_params.R    # writes lambda2, lambda3, c_composite, pi_recid into setup.R
```

- The model enforces `lambda1 = lambda3 * c_composite / beta_scale`.
- On the stored fits `update_rate_params.R` leaves `setup.R` unchanged (checked by
  `tests/unit/test_rate_params.R`). After a change, rerun every step below.

## 1. Deterministic calibration (prior centre)

```bash
Rscript scripts/run_calibration.R --run-id r_20261013_agepop \
  --maxit 3000 --n-starts 7 --sd-perturb 0.9 --seed 20261001 --cores 8
```

- **Equilibrium:** every objective evaluation runs model time 0–150 and checks equilibrium over the
  last 10 years.
- **Starts:** the anchor, 6 perturbed starts and a population-informed start, run in parallel.
- **Result:**
  - best objective 22.58 (`perturb_5`, converged; the other 7 starts reached `maxit`);
  - seroprevalence RMSE 0.033 and population MAPE 5.3%.
- **Structural misfit:** the 60+ seroprevalence is 0.508 against a target of 0.429. With the
  composite re-arrest rate, first arrests barely depend on `beta_scale`, so no parameter value
  fixes this.

## Shared posterior target

- **Prior:** independent Normal(`prior_centre(fit$best$theta)`, 0.20) on the 12 log-scale parameters.
  - `prior_centre()` (`scripts/calibration/likelihood.R`) moves each coordinate of the best fit to
    at least 5 prior SDs (1.0 on the log scale) inside the parameter bounds.
  - The best fit puts the 60+ contact scale at its lower bound (0.010), where the likelihood is
    flat. Without the move, about half the prior mass would lie outside the bounds.
  - With the move, 0 of 50,000 NPE prior draws fall outside the bounds.
- **Likelihood:** Binomial prevalence counts and log-normal populations (σ = 0.10).
- **Support:** draws outside the bounds and equilibrium failures are excluded (MCMC gives −∞; NPE
  drops the simulation and counts it).
- **Solver:** RK4 with `dt = 1/365`.
- **Details:** see `output/calibration/npe_vs_mcmc/model_audit.md`.

## 2. Reference MCMC (production)

```bash
Rscript scripts/run_mcmc_reference.R --cores 4   # about 4.3 hours
```

- **Defaults:**
  - 4 chains from independent prior draws;
  - warmup 10,000, with adaptation from proposal 5,000;
  - 50,000 evaluations per chain;
  - draws stored every 10th iteration (`--thin`);
  - seed 20261101;
  - fit `output/calibration/r_20261013_agepop/fit.rds`;
  - output to `output/calibration/mcmc_reference/`.
- **Sampler:** adaptive random-walk Metropolis (`arwm_budget` in `scripts/calibration/mcmc.R`). The
  initial proposal is `0.1² × diag(prior variance)`.
- **Diagnostics written:** trace plots, per-chain summaries, adaptation drift, covariance eigenvalues,
  jump sizes and autocorrelation times.
- **Result** (2026-10-14, these defaults):
  - R-hat ≤ 1.002;
  - bulk ESS 2,824–4,259 and tail ESS ≥ 5,119;
  - MCSE of the mean ≤ 0.019 posterior SD;
  - acceptance 0.26–0.30.
- **Prior-dominated parameters:** the 60+ contact scale (θ6) and the 60+ beta scale (θ12) have
  posterior SD ≈ the prior SD (0.20). The data do not inform them.
- `--targets` swaps the target file. The MCMC SBC uses it.

## 3. Cross-validation: single-round NPE with SBC

```bash
Rscript scripts/run_npe.R --step all --n-sims 50000 --seed 20261012 --train-seed 20261013 \
  --n-sbc 1000 --n-cores 14 --threads 4 --python /path/to/python
```

- **LHS design** (`scripts/calibration/prior_sampling.R`): a 50,000 × 12 design with one point per
  stratum, mapped through the Normal inverse CDF of each parameter's prior.
- **Simulations:** draws outside the bounds, non-finite results and equilibrium failures are dropped
  and counted, never replaced. This run: 50,000 valid, 0 excluded (`simulation_summary.json`).
- **Network:** an equal-weight ensemble of 5 MAFs (50×5), each with batch 256, Adam lr 5e-4,
  validation fraction 0.10, patience 20, ≤ 250 epochs, gradient clipping 5.0, CPU only.
  - A single network gave a rank bias of up to 0.03 on one parameter, and which parameter
    depended on the seed. The ensemble removes it.
  - Training takes about 4 minutes.
- **Output:** 20,000 posterior draws.
- **SBC:** 1,000 independent prior datasets, with ranks of the true value in 1,000 posterior draws.
  - **Gate:** cov95 in [0.85, 1] and cov90 in [0.75, 1] for every parameter, and no systematically
    skewed rank histogram (U-shaped, modal, biased or Bonferroni chi-square non-uniform).
  - **Result:** cov95 0.931–0.958, cov90 0.873–0.910, all 12 rank shapes uniform; gate PASS
    (`NPE_OK.txt`).
- **Comparison with the reference MCMC:**

  ```bash
  Rscript scripts/compare_npe_mcmc.R --policy-cores 14
  ```

## 4. Optional: direct MCMC SBC

```bash
Rscript scripts/sbc_mcmc.R --n-datasets 8 --max-evals 30000 --cores 8   # about 3.5 hours
```

- The truths are the first 8 NPE SBC datasets, with synthetic observations drawn from the MCMC
  likelihood. Each dataset gets one 30,000-evaluation chain.
- **Result:** pooled cov95 1.00 and cov90 0.95 over 96 ranks; mean rank 0.503; max split R-hat 1.021.
- **Limitation:** 8 datasets are a check against gross miscalibration only, not a full SBC.

## 5. 2015 predictions and calibration figure

```bash
Rscript scripts/predict_2015_posterior.R --max-draws 4000 --cores 14
Rscript scripts/plot_calibration.R --cores 14 --max-draws 2000
```

- **Prediction outputs:** `output/calibration/mcmc_reference/model_prediction_2015_strata.*`. The
  ever-PWID estimate is 21,169 (95% interval: 19,152–23,399).
- **Calibration figure:** `figures/calibration_agepop.png` compares the targets with the reference
  MCMC posterior and the NPE (LHS).
- **Overlap check:** `figures/calibration_interval_overlap.csv` records whether each 95% interval
  overlaps its target.

## 6. 2026 treatment scenarios

```bash
Rscript scripts/scenarios_2026.R --root . --draws-per-seed 2000 --cores 14 --sample-seed 20261005
```

- **Draws:** 2,000 draws from the reference MCMC posterior (the default `--posterior-files`). The full
  analysis is run per draw, then summarised as the median and 95% UI.
- **Options:**
  - `--posterior-files` accepts any CSV with `theta1..theta12`;
  - `--no-figures true` skips the figures;
  - `--min-valid-share` sets the share of draws that must pass equilibrium.
- **More detail:** see `output/scenarios_2026/README.md`.

## 7. Manuscript

```bash
Rscript scripts/paper_figures.R .
python3 scripts/build_paper.py .
```
