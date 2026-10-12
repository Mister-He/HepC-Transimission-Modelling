# Statistical-target audit: single-round NPE and reference MCMC

This audit checks that the two retained inference methods target the same posterior. It was first
written for the 2026-10-10 benchmark and revised on 2026-10-11 for the current workflow: single-round
NPE on 50,000 LHS prior simulations, and the reference MCMC.

## 1. Parameters and prior

- **Parameters.** There are 12 parameters θ on the natural-log scale (`scripts/calibration/likelihood.R`,
  `build_params`). θ1–θ6 multiply the six rows of the contact matrix and θ7–θ12 the six age-specific
  entry rates.
- **Prior.** Independent Normal(μ, 0.20²) for each θ_j. μ is `fit$best$theta` from
  `output/calibration/r_20261010_agepop/fit.rds`. There is no dependence and no truncation.
  - NPE implements it as a torch `Independent(Normal)` (`scripts/npe_train.py`, `build_prior`).
  - MCMC implements it as `dnorm` in `log_posterior_exact` (`scripts/calibration/mcmc.R`).
- **Support.** MCMC also applies `within_bounds`: exp(θ) must lie in [0.01, 1000] for contact
  scales and [0.01, 10⁶] for entry scales. The prior mass beyond these bounds is below 10⁻²⁰.
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
- **NPE.** `run_npe.R` (`simulate`) drops non-finite simulations and equilibrium failures without
  replacing them. It records each outcome in `prior_design.csv` and the counts in
  `simulation_summary.json`. Training on the surviving draws is training under the truncated prior.
- **Nelder–Mead.** Its plausibility and equilibrium penalties enter only the objective that sets μ.
  Neither method's target contains them.

## 5. NPE implementation

- **Package and loss.** `sbi` 0.27.0 `NPE_C`, trained in a single round on prior simulations, which
  uses the standard negative log-likelihood loss. There is no sequential proposal and no atomic
  (APT) loss.
- **Network.** MAF with 50 hidden features and 5 transforms, batch 256, Adam with lr 5e-4, validation
  fraction 0.10, early stopping after 20 epochs, at most 250 epochs, gradient clipping 5.0, and `sbi`
  per-feature z-scoring of θ and x.
- **Hardware.** CPU only, with 4 torch threads.
- **Seeds:**
  - prior design 20261012;
  - training seed 20261013, which gives the noise stream [20261013, 11], torch training seed
    2026101321 and torch sampling seed 2026101331.

## 6. Reference MCMC implementation

- **Sampler.** Adaptive random-walk Metropolis (`arwm_budget`, `scripts/calibration/mcmc.R`) with
  4 chains started from independent prior draws.
- **Initial proposal.** Covariance is 0.1² × the prior variance, estimated from 2,048 prior draws.
- **Adaptation.** Begins after 200 proposals, refreshes every 50, and freezes after warmup. A 5%
  safety mixture always uses the initial covariance.
- **Run length.** 1 initial evaluation + 4,999 warmup + 45,000 retained proposals per chain; seed
  20261101.
- **Diagnostics.** Rank-normalised R-hat ≤ 1.002, bulk ESS ≥ 3,433, and MCSE of the mean ≤ 0.017
  posterior SD.
- **Script.** `scripts/run_mcmc_reference.R`. The results in `output/calibration/mcmc_reference/` were
  produced on 2026-10-10 by the same sampler under its earlier name `run_mcmc_validation.R`. They were
  moved here, not rerun.

## 7. Conclusion

NPE and MCMC target the same posterior on θ:
- the same prior;
- the same likelihood and observation model;
- the same solver;
- the same equilibrium restriction.

The only difference is the MCMC bounds, which have negligible prior mass.
