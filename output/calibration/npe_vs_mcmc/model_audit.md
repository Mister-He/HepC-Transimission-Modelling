# Statistical-target audit: single-round NPE and reference MCMC

This audit checks that the two retained inference methods target the same posterior. It was first
written for the 2026-10-10 benchmark and revised on 2026-10-14 for the current workflow: the reference
MCMC (production) and single-round NPE on 50,000 LHS prior simulations (cross-validation only).

## 1. Parameters and prior

- **Parameters.** There are 12 parameters θ on the natural-log scale (`scripts/calibration/likelihood.R`,
  `build_params`). θ1–θ6 multiply the six rows of the contact matrix and θ7–θ12 the six age-specific
  entry rates.
- **Prior.** Independent Normal(μ, 0.20²) for each θ_j, with no dependence.
  - μ is `prior_centre(fit$best$theta)` (`likelihood.R`), using
    `output/calibration/r_20261013_agepop/fit.rds`.
  - `prior_centre()` moves each coordinate to at least 5 prior SDs (1.0 on the log scale) inside the
    bounds. Only θ6, the 60+ contact scale, moves: the best fit has it at the lower bound
    log(0.01), where the likelihood is flat.
  - Both methods use the same μ. The NPE run records it in `simulation_summary.json`, the MCMC run in
    `mcmc_run_config.json`.
  - NPE implements it as a torch `Independent(Normal)` (`scripts/npe_train.py`, `build_prior`).
  - MCMC implements it as `dnorm` in `log_posterior_exact` (`scripts/calibration/mcmc.R`).
- **Support.** Both methods apply `within_bounds`: exp(θ) must lie in [0.01, 1000] for contact
  scales and [0.01, 10⁶] for entry scales.
  - The prior mass beyond these bounds is at most 2.9 × 10⁻⁷ per parameter (5 SD).
  - 0 of the 50,000 NPE design draws fell outside the bounds.
- **Transformations.** None. Both methods sample θ directly, so no Jacobian is needed.
- **LHS design for NPE** (`scripts/calibration/prior_sampling.R`):
  - a 50,000 × 12 Latin hypercube on the unit cube, with one point per stratum in every column and
    design seed 20261012;
  - each column is mapped with the Normal inverse CDF `qnorm(u, μ_j, 0.20)`.

  This targets the stated prior exactly, because the prior is independent across parameters. It is
  not uniform sampling between bounds.

## 2. Simulator

- **Solver.** Fixed-step RK4 with `dt = 1/365` over model years 0–150, with ageing and non-negativity
  steps (`scripts/model/sim.cpp`). Both methods use this solver unchanged.
- **Equilibrium check.** `check_equilibrium` (`scripts/calibration/equilibrium.R`) compares years 140
  and 150 and requires all of the following:
  - every age group's prison population changes by at most 0.01 in log terms;
  - every age group's seroprevalence changes by at most 0.005;
  - total HCV, DC and HCC each change by at most 0.01 in log terms;
  - every compartment changes by at most 0.02 in log terms.

## 3. Data and observation model

- **Data vector.** x = (p₁..p₆, N₁..N₆), the 2015 prison seroprevalence and population counts by age
  (`data/calibration_targets.json`).
- **Likelihood.**
  p(x | θ) = ∏ Binomial(xₐ | Nₐ, p̂ₐ(θ)) × ∏ LogNormal(Nₐ | log N̂ₐ(θ), 0.10)
  (`likelihood.R`: `nll_prev`, `nll_pop`).
- **NPE training data.** NPE trains on the deterministic (p̂, N̂). It adds Binomial(Nₐ, p̂ₐ)/Nₐ and
  log-normal (σ = 0.10) noise, so it learns p(θ | x) under the same likelihood and conditions on the
  same 12 observations. The network input is the bijective transform [logit p, log N].

## 4. Restrictions that change the target

- **Shared target.** posterior ∝ prior × likelihood × 1{finite and passes the equilibrium check}.
- **MCMC.** `log_posterior_exact` returns −∞ for failures. There is no soft penalty.
- **NPE.** `run_npe.R` (`simulate`) drops draws outside the bounds, non-finite simulations and
  equilibrium failures, without replacing them. All 50,000 draws were valid in the current run. It records each outcome in `prior_design.csv` and the counts in
  `simulation_summary.json`. Training on the surviving draws is training under the truncated prior.
- **Nelder–Mead.** Its plausibility and equilibrium penalties enter only the objective that sets μ.
  Neither method's target contains them.

## 5. NPE implementation

- **Package and loss.** `sbi` 0.27.0 `NPE_C`, trained in a single round on prior simulations, which
  uses the standard negative log-likelihood loss. There is no sequential proposal and no atomic
  (APT) loss.
- **Network.** An equal-weight ensemble of 5 MAFs; sampling draws n/5 from each member. Each MAF has 50 hidden features and 5 transforms, batch 256, Adam with lr 5e-4, validation
  fraction 0.10, early stopping after 20 epochs, at most 250 epochs, gradient clipping 5.0, and `sbi`
  per-feature z-scoring of θ and x.
- **Hardware.** CPU only, with 4 torch threads.
- **Seeds:**
  - prior design 20261012;
  - training seed 20261013, which gives the noise stream [20261013, 11], torch training seeds
    2026101321–2026101325 (one per member) and torch sampling seed 2026101331.
- **SBC.** 1,000 random prior datasets, ranked in 1,000 posterior draws each. cov95 is 0.931–0.958,
  cov90 0.873–0.910, and all 12 rank histograms are uniform; gate PASS.

## 6. Reference MCMC implementation

- **Sampler.** Adaptive random-walk Metropolis (`arwm_budget`, `scripts/calibration/mcmc.R`) with
  4 chains started from independent prior draws.
- **Initial proposal.** Covariance is 0.1² × the prior variance, estimated from 2,048 prior draws.
- **Adaptation.** Begins after 5,000 proposals, refreshes every 50, and freezes after the
  10,000-proposal warmup. A 5% safety mixture always uses the initial covariance.
- **Run length.** 50,000 evaluations per chain (10,000 warmup, 40,000 retained), with every 10th
  retained draw stored; seed 20261101.
- **Diagnostics.** These use every retained draw: rank-normalised R-hat ≤ 1.002, bulk ESS ≥ 2,824,
  tail ESS ≥ 5,119, and MCSE of the mean ≤ 0.019 posterior SD.
- **Script.** `scripts/run_mcmc_reference.R` with its defaults and `--cores 4`, run 2026-10-13/14.
- **Direct SBC** (`scripts/sbc_mcmc.R`, `mcmc_sbc/`): 8 datasets, one 30,000-evaluation chain each.
  Pooled cov95 is 1.00 and cov90 0.95 over 96 ranks. This is a check against gross failure only.

## 7. Conclusion

NPE and MCMC target the same posterior on θ:
- the same prior;
- the same likelihood and observation model;
- the same solver;
- the same equilibrium restriction.

They also share the same bounds and the same prior centre.
