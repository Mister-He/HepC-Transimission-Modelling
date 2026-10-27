# Statistical-target audit: single-round NPE and reference MCMC

This audit checks that the two retained inference methods target the same posterior. It was first
written for the 2026-10-10 benchmark and revised on 2026-10-27 for the current workflow (first
arrest = `c` × re-arrest, calibrated entry scales on the CNB counts, three likelihood terms, PWID
mortality ratios by viraemic status, fixed prior; design in `docs/redesign_prompt.md`): the
reference MCMC (production) and single-round NPE on 50,000 LHS prior simulations
(cross-validation only).

## 1. Parameters and prior

- **Parameters.** There are 13 parameters θ on the natural-log scale (`scripts/calibration/likelihood.R`,
  `build_params`).
  - θ1–θ6 multiply the six rows of the contact matrix.
  - θ7–θ12 multiply the entry into D by age (`beta` = CNB 2015 new drug abusers × exp(θ)).
  - θ13 = log `c`: first arrest `λ1_i = c·λ3_i` (`c_first` in `sim.cpp`).
- **Prior.** Independent Normal for each θ_j, with no dependence, truncated to the bounds and fixed
  before fitting (`prior_spec()`): θ1–θ6 ~ N(0, 1.5), θ7–θ12 ~ N(0, 1) truncated at 0
  (half-normal), θ13 ~ N(0, 0.5). It does not use the data.
  - Both methods use the same prior. The NPE run records it in `simulation_summary.json`, the MCMC
    run in `mcmc_run_config.json`.
  - NPE implements it as a torch `Independent(Normal)` (`scripts/npe_train.py`, `build_prior`),
    with posterior draws outside the bounds rejected (`Ensemble`).
  - MCMC uses `-nlp(θ)` (the same function Nelder–Mead uses) in `log_posterior_exact`
    (`scripts/calibration/mcmc.R`); outside the bounds the target is −∞. The truncation constant does
    not depend on θ. Chain starts are truncated-prior draws.
- **Support.** Both methods apply `within_bounds`: exp(θ) in [0.01, 1000] for contact scales,
  [1, 100] for entry scales and `c` in [0.5, 100].
- **Transformations.** None. Both methods sample θ directly, so no Jacobian is needed.
- **LHS design for NPE** (`scripts/calibration/prior_sampling.R`):
  - a 50,000 × 13 Latin hypercube on the unit cube, with one point per stratum in every column and
    design seed 20261012;
  - each column is mapped with the truncated Normal inverse CDF
    `qnorm(Φ(lo) + u·(Φ(hi) − Φ(lo)), μ_j, σ_j)`.

  This targets the stated truncated prior exactly, because the prior is independent across
  parameters. 0 of the 50,000 design draws fell outside the bounds.

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

- **Data vector** (`data/calibration_targets.json`), 13 values:
  - y₁..y₆: the 2015 prison seroprevalence by age;
  - n₁..n₆: CNB 2015 new drug abusers by age (257, 644, 292, 94, 21, 1);
  - m: Teo et al. 2019 community PWID, 11,000.
- **Likelihood** (`likelihood.R`: `nll_targets`), three independent terms:
  - ∏ Beta(yₐ | p̂ₐ·φ, (1 − p̂ₐ)·φ), φ = 100 (`nll_prev`);
  - ∏ NegBin(nₐ | mean f̂ₐ, size 50), f̂ₐ = c·λ3ₐ·Dₐ in 2015 (`nll_first_arrests`);
  - Normal(log m | log(D + F), 0.2453) (`nll_community`).
  The model values come from `calib_summary()` (`model_metrics.R`), shared by every method.
- **NPE training data.** NPE trains on the deterministic summaries and adds the same noise: Beta
  with the same φ on prevalences, gamma-Poisson (negative binomial, size 50) on first arrests and
  log-normal (SD 0.2453) on D + F. So it learns p(θ | y) under the same likelihood. The network
  input is logit y (clipped to [10⁻⁴, 1 − 10⁻⁴]), log(1 + n) and log m.

## 4. Restrictions that change the target

- **Shared target.** posterior ∝ prior × likelihood × 1{within bounds, finite and passes the
  equilibrium check}.
- **MCMC.** `log_posterior_exact` returns −∞ for failures. There is no soft penalty.
- **NPE.** `run_npe.R` (`simulate`) drops draws outside the bounds, non-finite simulations and
  equilibrium failures, without replacing them. The current run (2026-10-27) has 49,941 valid
  draws, 47 equilibrium failures and 12 non-finite simulations. It records each outcome in
  `prior_design.csv` and the counts in `simulation_summary.json`.
- **Nelder–Mead.** It minimises the same negative log posterior (with an equilibrium penalty in
  place of −∞) and only reports the posterior mode; neither method's target depends on it.

## 5. NPE implementation

- **Package and loss.** `sbi` 0.27.0 `NPE_C`, trained in a single round on prior simulations, which
  uses the standard negative log-likelihood loss. There is no sequential proposal and no atomic
  (APT) loss.
- **Network.** An equal-weight ensemble of 5 MAFs; sampling draws n/5 from each member. Each MAF has
  50 hidden features and 5 transforms, batch 256, Adam with lr 5e-4, validation fraction 0.10,
  early stopping after 20 epochs, at most 250 epochs, gradient clipping 5.0, and `sbi` per-feature
  z-scoring of θ and x.
- **Hardware.** CPU only, with 4 torch threads.
- **Seeds:**
  - prior design 20261012;
  - training seed 20261013, which gives the noise stream [20261013, 11], torch training seeds
    2026101321–2026101325 (one per member) and torch sampling seed 2026101331.
- **SBC.** 1,000 random prior datasets, of which 998 passed. Each was ranked in 1,000 posterior
  draws. cov95 is 0.926–0.962 and cov90 0.871–0.917, but the rank histograms of θ2, θ10 and θ11
  are skewed; gate FAIL (`npe_lhs/NPE_NOT_OK.txt`).
- **Support at the observed data.** The NPE posterior does not respect the equilibrium restriction:
  33% of its draws fail the check. `compare_npe_mcmc.R` drops them before comparing.

## 6. Reference MCMC implementation

- **Sampler.** Adaptive random-walk Metropolis (`arwm_budget`, `scripts/calibration/mcmc.R`) with
  8 chains started from independent prior draws.
- **Initial proposal.** Covariance is 0.1² × the prior variance, estimated from 2,048 prior draws.
- **Adaptation.** Begins after 5,000 proposals, refreshes every 50 from the recent half of the
  warmup states, and freezes after the 20,000-proposal warmup. A 5% safety mixture always uses the initial covariance.
- **Run length.** 60,000 evaluations per chain (20,000 warmup, 40,000 retained), with every 10th
  retained draw stored; seed 20261101.
- **Diagnostics.** These use every retained draw: rank-normalised R-hat ≤ 1.004, bulk ESS ≥ 2,158,
  tail ESS ≥ 1,936, and MCSE of the mean ≤ 0.023 posterior SD. Acceptance 0.12–0.14.
- **Script.** `scripts/run_mcmc_reference.R` with its defaults (8 cores), run 2026-10-27 (1.9 h).
- **Direct SBC** (`scripts/sbc_mcmc.R`, `mcmc_sbc/`): 8 datasets, one 45,000-evaluation chain each
  (15,000 warmup).
  - Pooled over 104 ranks: cov95 0.952, cov90 0.942, mean rank 0.467; max split R-hat 1.015
    (one chain per dataset; three datasets are between 1.01 and 1.015). Lowest per-parameter
    cov95: θ6 6/8.
  - This is a check against gross failure only. Ranks within a dataset are correlated.

## 7. Conclusion

NPE and MCMC target the same posterior on θ:
- the same prior and bounds;
- the same likelihood and observation model;
- the same solver;
- the same equilibrium restriction.

The targets match, but the NPE approximation is inaccurate at the observed data (`report.md`).
