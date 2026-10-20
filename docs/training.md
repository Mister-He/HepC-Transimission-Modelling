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
Rscript scripts/model/update_rate_params.R    # writes lambda2, lambda3, f_first, pi_recid into setup.R
```

- **Population:** the modelled PWID population is everyone ever held for drug consumption with a
  known HCV screening result. It gives the 2015 prison counts (equal to the calibration targets)
  and the first offenders.
- **Rates come from all consumption records.** Screening (Dec 2014–Feb 2016) selects people in
  prison then, so in the screened subset every earlier release ends in re-arrest (re-arrest share
  0.99, against 0.72 in all records).
- **`rearrest_rate.stan`:** its 2015 prison steady state uses the ODE's ageing (0.1 per year) and
  mortality (`mu·omega`). It gives the first-arrest flow `f_pop` (2026-10-15 fit: 32, 304, 102,
  75.9, 66.1 and 19.6 per year; re-arrest share `pi` 0.720, R-hat 1.00).
- **`sim.cpp` derives `lambda1` from `f_first`** through the HCV-free steady state of D
  (`first_arrest_rate()`). An entry flow too small to supply `f_first` is infeasible: the
  parameter set has zero posterior density, and NPE drops and counts it.
- **Idempotent:** on the stored fits `update_rate_params.R` leaves `setup.R` unchanged (checked
  by `tests/unit/test_rate_params.R`). After a change, rerun every step below.

## 1. Deterministic calibration (prior centre)

```bash
Rscript scripts/run_calibration.R --run-id r_20261014_agepop \
  --maxit 3000 --nm-restarts 5 --n-starts 7 --sd-perturb 0.9 --seed 20261001 --cores 8   # ~2.5 h
```

- **Equilibrium:** every objective evaluation runs model time 0–150 and checks equilibrium over the
  last 10 years.
- **Starts:** the anchor, 6 perturbed starts and a population-informed start, run in parallel.
- **Restarts:** each start reruns Nelder–Mead from its own solution while the objective still falls
  by more than 0.001 (`--nm-restarts`, at most 5 restarts). Without restarts every start stopped
  early (best 25.9).
- **Result (2026-10-15):**
  - best objective 21.44 (`population_informed`); five starts lie within 21.44–21.61;
  - seroprevalence RMSE 0.012 and population MAPE 5.7%;
  - 60+ seroprevalence 0.455 against a target of 0.429 (0.65 SE); every other age is within 0.01;
  - 2015 PWID 18,314, of whom 3,290 are never arrested;
  - entry is 0.34–1.1× the CNB counts.
- **First-arrest flow:** at the best fit, the model's 2015 first-arrest flow is within 0.04–2.7%
  of `f_first`. The shortfall comes from HCV-related deaths in D.
- **Viraemic share fixed point:**
  - `run_config.json` records the chronic viraemic share of all PWID at the 2015 snapshot
    (`viraemic_share_2015`), and `setup.R` freezes it as `viraemic_share_ref`.
  - Rounds under the final model: 0.242 → 0.205 → 0.199.
  - It stops when the implied background-mortality multiplier `omega_0` changes by less than
    0.5%. The last step changes it by 0.26%.
- **HR sensitivity** (`output/calibration/hr_sensitivity/`, 3 starts each, `--hr-viraemic`):

  | HR | 60+ seroprevalence | Objective |
  | --- | --- | --- |
  | 1 | 0.494 | 22.70 |
  | 1.11 | 0.493 | 22.16 |
  | 1.23 | 0.481 | 21.75 |
  | 1.47 (production) | 0.455 | 21.44 |
  | 1.77 | 0.439 | 21.19 |

  - With the entry anchor, the viraemic excess mortality is what brings 60+ close to its target.
  - The earlier misfit was 0.508 under `c_composite`.

## Shared posterior target

- **Parameters:** θ1–θ6 are log contact scales. θ7–θ12 are log excess-entry scales
  (`entry_scales()` in `scripts/calibration/likelihood.R`).
  - Entry into age i is `max(0, f_i − 0.1·D_{i−1}) + exp(θ_{6+i})·β_i`: the minimum that supplies
    the first-arrest flow, plus an excess.
  - So every θ is feasible. On the earlier log entry scale, 43% of prior draws could not supply
    the first arrests, mainly at 20–29, where the fit sits close to that edge.
- **Prior:** independent Normal(`prior_centre(fit$best$theta)`, 0.20) on the 12 log-scale parameters.
  - `prior_centre()` (`scripts/calibration/likelihood.R`) moves each coordinate of the best fit to
    at least 5 prior SDs (1.0 on the log scale) inside the parameter bounds.
  - The best fit puts the 60+ contact scale at its lower bound (0.010), where the likelihood is
    flat. Without the move, about half the prior mass would lie outside the bounds.
  - With the move, 0 of 50,000 NPE prior draws fall outside the bounds.
- **Likelihood:**
  - Binomial prevalence counts.
  - Log-normal populations (σ = 0.10).
  - Log-normal 2015 entry into D against the CNB new-drug-abuser counts (σ = 1.0, `nll_entry()`).
    With first arrests fixed by `f_first`, the prison data do not identify the never-arrested
    pool. Without this term, fits with about 38,000 and about 1.7 million PWID gave the same
    objective.
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
  - fit `output/calibration/r_20261014_agepop/fit.rds`;
  - output to `output/calibration/mcmc_reference/`.
- **Sampler:** adaptive random-walk Metropolis (`arwm_budget` in `scripts/calibration/mcmc.R`). The
  initial proposal is `0.1² × diag(prior variance)`.
- **Diagnostics written:** trace plots, per-chain summaries, adaptation drift, covariance eigenvalues,
  jump sizes and autocorrelation times.
- **Result** (2026-10-15, these defaults, 4.5 hours):
  - R-hat ≤ 1.003;
  - bulk ESS 2,319–3,495 and tail ESS ≥ 2,728;
  - MCMC standard error (MCSE) of the mean ≤ 0.021 posterior SD;
  - acceptance 0.25–0.27.
- **Prior-dominated parameters:** eight parameters have posterior SD above 0.8 × the prior SD (0.20):
  θ5, θ6 (the 50–59 and 60+ contact scales) and θ7–θ12 (the excess-entry scales).
- `--targets` swaps the target file. The MCMC SBC uses it.

## 3. Cross-validation: single-round NPE with SBC

```bash
Rscript scripts/run_npe.R --step all --n-sims 50000 --seed 20261012 --train-seed 20261013 \
  --n-sbc 1000 --n-cores 14 --threads 4 --python /path/to/python
```

- **LHS design** (`scripts/calibration/prior_sampling.R`): a 50,000 × 12 design with one point per
  stratum, mapped through the Normal inverse CDF of each parameter's prior.
- **Simulations:** draws outside the bounds, non-finite results and equilibrium failures are dropped
  and counted, never replaced. This run: 48,918 valid and 1,082 equilibrium failures
  (`simulation_summary.json`).
- **Summaries:** 18: prevalence, population and the 2015 entry into D (log-normal noise SD 1.0).
- **Network:** an equal-weight ensemble of 5 MAFs (50×5), each with batch 256, Adam lr 5e-4,
  validation fraction 0.10, patience 20, ≤ 250 epochs, gradient clipping 5.0, CPU only.
  - A single network gave a rank bias of up to 0.03 on one parameter, and which parameter
    depended on the seed. The ensemble removes it.
  - Training takes about 5 minutes.
- **Output:** 20,000 posterior draws.
- **SBC:** 1,000 independent prior datasets (973 pass equilibrium), with ranks of the true value in
  1,000 posterior draws.
  - **Gate:** cov95 in [0.85, 1] and cov90 in [0.75, 1] for every parameter, and no systematically
    skewed rank histogram (U-shaped, modal, biased or Bonferroni chi-square non-uniform).
  - **Result:** cov95 0.936–0.958, cov90 0.879–0.918, all 12 rank shapes uniform; gate PASS
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
- **Result:** pooled cov95 0.93 and cov90 0.83 over 96 ranks; mean rank 0.512; max split R-hat 1.016.
- **Limitation:** 8 datasets are a check against gross miscalibration only, not a full SBC.

## 5. 2015 predictions and calibration figure

```bash
Rscript scripts/predict_2015_posterior.R --max-draws 4000 --cores 14
Rscript scripts/plot_calibration.R --cores 14 --max-draws 2000
```

- **Prediction outputs:** `output/calibration/mcmc_reference/model_prediction_2015_strata.*`. The
  ever-PWID estimate is 18,453 (95% interval: 17,799–19,380).
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
