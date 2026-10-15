# Statistical-target audit: single-round NPE and reference MCMC

This audit checks that the two retained inference methods target the same posterior. It was first
written for the 2026-10-10 benchmark and revised on 2026-10-15 for the current workflow (first-arrest flow
handoff, excess-entry parameterisation, CNB entry likelihood, viraemic excess mortality): the reference
MCMC (production) and single-round NPE on 50,000 LHS prior simulations (cross-validation only).

## 1. Parameters and prior

- **Parameters.** There are 12 parameters θ on the natural-log scale (`scripts/calibration/likelihood.R`,
  `build_params`).
  - θ1–θ6 multiply the six rows of the contact matrix.
  - θ7–θ12 are log excess-entry scales (`entry_scales`). Entry into age i is
    `max(0, f_i − 0.1·D_{i−1}) + exp(θ_{6+i})·β_i`, so every θ can supply the first-arrest flow.
- **Prior.** Independent Normal(μ, 0.20²) for each θ_j, with no dependence.
  - μ is `prior_centre(fit$best$theta)` (`likelihood.R`), using
    `output/calibration/r_20261014_agepop/fit.rds`.
  - `prior_centre()` moves each coordinate to at least 5 prior SDs (1.0 on the log scale) inside the
    bounds. θ3 and θ6 (the 30–39 and 60+ contact scales) move, because the best fit puts them within
    1.0 of the lower bound log(0.01).
  - Both methods use the same μ. The NPE run records it in `simulation_summary.json`, the MCMC run in
    `mcmc_run_config.json`.
  - NPE implements it as a torch `Independent(Normal)` (`scripts/npe_train.py`, `build_prior`).
  - MCMC implements it as `dnorm` in `log_posterior_exact` (`scripts/calibration/mcmc.R`).
- **Support.** Both methods apply `within_bounds`: exp(θ) must lie in [0.01, 1000] for contact
  scales and [10⁻⁴, 10⁶] for excess-entry scales.
  - The prior mass beyond these bounds is at most 2.9 × 10⁻⁷ per parameter (5 SD).
  - 0 of the 50,000 NPE design draws fell outside the bounds, and none could fail to supply the
    first-arrest flow.
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

- **Data vector.** x = (p₁..p₆, N₁..N₆, E₁..E₆): the 2015 prison seroprevalence and population counts
  by age, plus the 2015 CNB new-drug-abuser counts as the observed entry into D
  (`data/calibration_targets.json`).
- **Likelihood.**
  p(x | θ) = ∏ Binomial(xₐ | Nₐ, p̂ₐ(θ)) × ∏ LogNormal(Nₐ | log N̂ₐ(θ), 0.10) × ∏ LogNormal(Eₐ | log Êₐ(θ), 1.0)
  (`likelihood.R`: `nll_prev`, `nll_pop`, `nll_entry`; Ê is `pm$beta`).
- **NPE training data.** NPE trains on the deterministic (p̂, N̂, Ê). It adds Binomial(Nₐ, p̂ₐ)/Nₐ and
  log-normal noise (σ = 0.10 for N, 1.0 for E), so it learns p(θ | x) under the same likelihood and
  conditions on the same 18 observations. The network input is [logit p, log N, log E].

## 4. Restrictions that change the target

- **Shared target.** posterior ∝ prior × likelihood × 1{finite and passes the equilibrium check}.
- **MCMC.** `log_posterior_exact` returns −∞ for failures. There is no soft penalty.
- **NPE.** `run_npe.R` (`simulate`) drops draws outside the bounds, non-finite simulations and
  equilibrium failures, without replacing them. The current run has 48,918 valid draws and 1,082
  equilibrium failures. It records each outcome in `prior_design.csv` and the counts in
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
- **SBC.** 1,000 random prior datasets, of which 973 passed the equilibrium check. Each was ranked in
  1,000 posterior draws. cov95 is 0.936–0.958, cov90 0.879–0.918, and all 12 rank histograms are
  uniform; gate PASS.

## 6. Reference MCMC implementation

- **Sampler.** Adaptive random-walk Metropolis (`arwm_budget`, `scripts/calibration/mcmc.R`) with
  4 chains started from independent prior draws.
- **Initial proposal.** Covariance is 0.1² × the prior variance, estimated from 2,048 prior draws.
- **Adaptation.** Begins after 5,000 proposals, refreshes every 50, and freezes after the
  10,000-proposal warmup. A 5% safety mixture always uses the initial covariance.
- **Run length.** 50,000 evaluations per chain (10,000 warmup, 40,000 retained), with every 10th
  retained draw stored; seed 20261101.
- **Diagnostics.** These use every retained draw: rank-normalised R-hat ≤ 1.003, bulk ESS ≥ 2,319,
  tail ESS ≥ 2,728, and MCSE of the mean ≤ 0.021 posterior SD.
- **Script.** `scripts/run_mcmc_reference.R` with its defaults and `--cores 4`, run 2026-10-15
  (16,231 s).
- **Direct SBC** (`scripts/sbc_mcmc.R`, `mcmc_sbc/`): 8 datasets, one 30,000-evaluation chain each.
  - Synthetic entry observations are `E·exp(1.0 z)`.
  - Pooled over 96 ranks: cov95 0.93, cov90 0.83 and mean rank 0.51; max split R-hat 1.016.
  - This is a check against gross failure only. Ranks within a dataset are correlated, so the low
    cov90 is within noise for 8 datasets.

## 7. Conclusion

NPE and MCMC target the same posterior on θ:
- the same prior;
- the same likelihood and observation model;
- the same solver;
- the same equilibrium restriction.

They also share the same bounds and the same prior centre.
