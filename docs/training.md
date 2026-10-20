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

## 0. Rate fits and handoff (only after refitting the Stan rate models)

```bash
Rscript scripts/model/preprocess.R            # needs data/sps.csv (not tracked by git); ~15 min
Rscript scripts/model/update_rate_params.R    # writes lambda2, lambda3, pi_recid and lambda1 (= lambda3) into setup.R
```

- **Population:** the modelled PWID population is everyone ever held for drug consumption with a
  known HCV screening result. It gives the 2015 prison counts (equal to the calibration targets)
  and the first offenders.
- **Rates come from all consumption records.** Screening (Dec 2014–Feb 2016) selects people in
  prison then, so in the screened subset every earlier release ends in re-arrest (re-arrest share
  0.99, against 0.72 in all records).
- **`rearrest_rate.stan`:** its 2015 prison steady state uses the ODE's ageing (0.1 per year) and
  mortality (`mu·omega`). Its first-arrest flow `f_pop` informs `lambda3` there only; it is not
  passed to the ODE (re-arrest share `pi` 0.720, R-hat 1.00).
- **First-arrest rate:** the records never observe the never-arrested, so `lambda1` (D → J) is
  assumed equal to `lambda3` and hard-coded in `setup.R` (0.546, 0.488, 0.472, 0.472, 0.386,
  0.266 per year). `update_rate_params.R` writes it.
- **Idempotent:** on the stored fits `update_rate_params.R` leaves `setup.R` unchanged (checked
  by `tests/unit/test_rate_params.R`). After a change, rerun every step below.

## 1. Deterministic calibration (prior centre)

```bash
Rscript scripts/run_calibration.R --run-id r_20261020_agepop \
  --maxit 3000 --nm-restarts 5 --n-starts 7 --sd-perturb 0.9 --seed 20261001 --cores 8   # ~2.5 h
```

- **Equilibrium:** every objective evaluation runs model time 0–150 and checks equilibrium over the
  last 10 years.
- **Starts:** the anchor, 6 perturbed starts and a population-informed start, run in parallel.
- **Restarts:** each start reruns Nelder–Mead from its own solution while the objective still falls
  by more than 0.001 (`--nm-restarts`, at most 5 restarts). Without restarts every start stopped
  early (best 25.9).
- **Result (2026-10-20):**
  - best objective 18.73 (`perturb_2`); seven of the eight starts lie within 18.734–18.747 (the
    plausible anchor stops at 35.8);
  - seroprevalence RMSE 0.00004 and population MAPE 0.05%: every target, including 60+ (0.4286
    against 0.4286), is fitted almost exactly;
  - 2015 ever-PWID 17,187, of whom 1,352 are never arrested (point estimate);
  - entry is 0.07–4.1× the CNB counts at ages < 60 and 27× at 60+ (CNB count 1).
- **Viraemic share fixed point:**
  - `run_config.json` records the chronic viraemic share of all PWID at the 2015 snapshot
    (`viraemic_share_2015`), and `setup.R` freezes it as `viraemic_share_ref`.
  - Under the current model (`lambda1 = lambda3`): 0.205 → 0.196.
  - It stops when the implied background-mortality multiplier `omega_0` changes by less than
    0.5%. This step changes it by 0.38%, so `viraemic_share_ref` stays 0.205 and no refit was needed.
- **HR sensitivity** (`output/calibration/hr_sensitivity/`, 3 starts each, `--hr-viraemic`):

  | HR | 60+ seroprevalence | Objective |
  | --- | --- | --- |
  | 1 | 0.440 | 18.91 |
  | 1.11 | 0.436 | 18.79 |
  | 1.23 | 0.431 | 18.74 |
  | 1.47 (production) | 0.429 | 18.73 |
  | 1.77 | 0.429 | 18.73 |

  - With free entry scales, 60+ fits within 0.3 SE at every HR; the HR barely changes the fit.
  - The earlier misfit was 0.508 under `c_composite`.

## Shared posterior target

- **Parameters:** θ1–θ6 are log contact scales. θ7–θ12 are log entry scales: entry into age i
  in 2015 is `exp(θ_{6+i})·β_i`, with `β` the CNB 2015 new-drug-abuser counts
  (`build_params()` in `scripts/calibration/likelihood.R`).
- **Prior:** independent Normal(`prior_centre(fit$best$theta)`, 0.20) on the 12 log-scale parameters.
  - `prior_centre()` (`scripts/calibration/likelihood.R`) moves each coordinate of the best fit to
    at least 5 prior SDs (1.0 on the log scale) inside the parameter bounds.
  - The best fit puts the 60+ contact scale at its lower bound (0.010), where the likelihood is
    flat. Without the move, about half the prior mass would lie outside the bounds.
  - With the move, 0 of 50,000 NPE prior draws fall outside the bounds.
- **Likelihood:**
  - Binomial prevalence counts.
  - Log-normal populations (σ = 0.10).
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
  - fit `output/calibration/r_20261020_agepop/fit.rds`;
  - output to `output/calibration/mcmc_reference/`.
- **Sampler:** adaptive random-walk Metropolis (`arwm_budget` in `scripts/calibration/mcmc.R`). The
  initial proposal is `0.1² × diag(prior variance)`.
- **Diagnostics written:** trace plots, per-chain summaries, adaptation drift, covariance eigenvalues,
  jump sizes and autocorrelation times.
- **Result** (2026-10-20, these defaults, 4.7 hours):
  - R-hat ≤ 1.002;
  - bulk ESS 3,053–4,511 and tail ESS ≥ 5,335;
  - MCMC standard error (MCSE) of the mean ≤ 0.018 posterior SD;
  - acceptance 0.26–0.28.
- **Prior-dominated parameters:** four parameters have posterior SD above 0.8 × the prior SD (0.20):
  θ6 (the 60+ contact scale) and θ9, θ10, θ12 (the 30–39, 40–49 and 60+ entry scales).
- `--targets` swaps the target file. The MCMC SBC uses it.

## 3. Cross-validation: single-round NPE with SBC

```bash
Rscript scripts/run_npe.R --step all --n-sims 50000 --seed 20261012 --train-seed 20261013 \
  --n-sbc 1000 --n-cores 14 --threads 4 --python /path/to/python
```

- **LHS design** (`scripts/calibration/prior_sampling.R`): a 50,000 × 12 design with one point per
  stratum, mapped through the Normal inverse CDF of each parameter's prior.
- **Simulations:** draws outside the bounds, non-finite results and equilibrium failures are dropped
  and counted, never replaced. This run: 49,278 valid and 722 equilibrium failures
  (`simulation_summary.json`).
- **Summaries:** 12: prevalence and population.
- **Network:** an equal-weight ensemble of 5 MAFs (50×5), each with batch 256, Adam lr 5e-4,
  validation fraction 0.10, patience 20, ≤ 250 epochs, gradient clipping 5.0, CPU only.
  - A single network gave a rank bias of up to 0.03 on one parameter, and which parameter
    depended on the seed. The ensemble removes it.
  - Training takes about 4 minutes.
- **Output:** 20,000 posterior draws.
- **SBC:** 1,000 independent prior datasets (980 pass equilibrium), with ranks of the true value in
  1,000 posterior draws.
  - **Gate:** cov95 in [0.85, 1] and cov90 in [0.75, 1] for every parameter, and no systematically
    skewed rank histogram (U-shaped, modal, biased or Bonferroni chi-square non-uniform).
  - **Result:** cov95 0.944–0.961, cov90 0.891–0.918, all 12 rank shapes uniform; gate PASS
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
- **Result:** pooled cov95 0.95 and cov90 0.92 over 96 ranks; mean rank 0.498; max split R-hat 1.016.
- **Limitation:** 8 datasets are a check against gross miscalibration only, not a full SBC.

## 5. 2015 predictions and calibration figure

```bash
Rscript scripts/predict_2015_posterior.R --max-draws 4000 --cores 14
Rscript scripts/plot_calibration.R --cores 14 --max-draws 2000
```

- **Prediction outputs:** `output/calibration/mcmc_reference/model_prediction_2015_strata.*`. The
  ever-PWID estimate is 17,226 (95% interval: 15,649–18,968); 1,356 (1,247–1,471) are never arrested.
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
python scripts/figures.py .   # redraws every figure (PNG + editable SVG)
python3 scripts/build_paper.py .
```

`figures.py` redraws all figures from the saved outputs with the `scientific-visualization` skill
(`~/.claude/skills/scientific-visualization`, or pass `--skill-dir`). It needs Python ≥ 3.11 with
`matplotlib==3.11.2 pillow==12.3.0 pandas scipy`. `--only <name> ...` redraws selected figures.
Run it after any step above that changes its inputs. The R scripts still draw their own ggplot
versions of the same files, so run `figures.py` last.
