# Ablation and Diagnostics Plan: NPE-initialized Adaptive MCMC

## 1. Objective and scope

Evaluate and harden the reference MCMC (`scripts/run_mcmc_reference.R` +
`scripts/calibration/mcmc.R`) as an inference and validation tool. This document
specifies configuration changes, convergence/ergodicity diagnostics, and an
ablation study comparing initialization strategies. It is an executable plan,
not verified results.

Do not assume the current NPE-initialized adaptive Metropolis is well behaved.
Do not assume NPE is the best initialization. Report what is actually measured.

## 2. Working principles and guardrails

- Read `AGENTS.md`, `docs/AGENTS.md`, `docs/training.md`, and `docs/NPE_vs_MCMC_benchmark.md` first.
- Preserve the scientific model. Only change sampler configuration, diagnostics,
  and initialization. Do not change `sim.cpp`, `setup.R`, the likelihood, the
  prior, or the targets to improve results.
- Do not overwrite production posterior files. Write every experiment to an
  isolated `--out-dir` under `output/ablation_check/<run-id>/` (or an
  equivalent directory) and document which files are produced.
- Record the Git commit, and the md5 of `scripts/model/setup.R`,
  `scripts/model/sim.cpp`, `scripts/calibration/mcmc.R`,
  `scripts/run_mcmc_reference.R`, plus all CLI arguments and seeds.
- Long jobs must survive terminal closure (`nohup ... &`). Check for stale
  `mcmc` workers before starting another run.
- Preserve the target: the reference MCMC must condition on the same prior,
  likelihood, and equilibrium gate as the accepted NPE. Reconstruct the prior
  from `fit$best$theta` and pass the SDs that the NPE was trained with
  (`--contact-sd 0.20 --beta-sd 0.20`), otherwise the comparison is invalid.

## 3. Baseline to modify (current accepted reference MCMC)

| Item | Current value | Source |
| --- | --- | --- |
| Integration step | `dt = 0.0833` (monthly) | `output/calibration/mcmc_recal/mcmc_run_config.json`; `run_mcmc_reference.R:27,62` |
| Chains / iterations / burn-in / thin | 4 / 120000 / 20000 / 1 | same |
| Initial proposal covariance | `Sigma0 = cov(NPE)` + ridge | `run_mcmc_reference.R` (`Sigma0 <- stats::cov(npe)`, ~line 69-71) |
| Starting points | pooled NPE mean, `same_start = true` | `run_mcmc_reference.R` (`theta_mean`/`theta0`, ~line 72) |
| Adapt from | 200 proposals | `run_mcmc_reference.R` (`am_mcmc(... adapt_from = 200 ...)`, ~line 89) |
| Proposal scaling | `c_d = 2.38^2 / d`, `eps = 1e-8 I` | `mcmc.R:81-82,116` |
| Proposal covariance after adapt | running `M2/(n-1)` | `mcmc.R:113-115` |
| Prior | best-fit-centred Normal(0.2) on all log scales | `run_mcmc_reference.R:53-57` |
| Diagnostics currently emitted | classic R-hat (multivariate `coda::gelman.diag`), pooled ESS, mean acceptance | the `coda` diagnostics block in `run_mcmc_reference.R` |

Use this as the control arm. Every ablation arm must state how it differs.

## 4. Part A — Required configuration changes

### A1. Daily integration step (requirement 1)

- Change `integration_dt` to `1/365` (daily) for the reference MCMC
  (`run_mcmc_reference.R:27` default, then `data$dt <- DT` at `:62`).
- Keep the 150-year equilibrium horizon and model year 150 = 2015 unless a
  separate, justified experiment changes it.
- Budget the cost explicitly. Measured on the current model: one daily 150-year
  `run_sim` is about 0.36 s, so 4 chains x 120000 proposals is roughly
  48 CPU-hours (about 12 h wall on 4 cores). Confirm on the target machine and
  record it.
- Verify the daily-step target differs from the monthly-step target only by
  discretization; quantify the difference by comparing a short daily and monthly
  chain on the same prior.

### A2. Widening the thinning interval (requirement 2)

- Increase the thinning interval from the current retained `thin = 1` in the
  accepted run (`run_mcmc_reference.R` argument parsing and the
  `am_mcmc(... thin = THIN)` call; the script default is being changed from `5`
  to `1000`).
- Choose and justify the value from the measured autocorrelation time, not
  arbitrarily. Thinning does not add information; it only reduces storage and
  the autocorrelation between retained draws.
- Guard against over-thinning: with 4 chains x (120000 - 20000) iterations, a
  thin of 1000 leaves only ~400 retained draws, which is likely too few for
  reliable tail-ESS, MCSE, and 95% interval estimates. Report retained draws and
  achieved ESS; if tail ESS is insufficient, reduce the thinning interval or
  increase iterations.
- Keep one arm with no thinning for autocorrelation diagnostics.

### A3. Increasing `adapt_from` from 200 to 5000 (requirement 3)

- Pass `adapt_from = 5000` to `am_mcmc` (the `am_mcmc(...)` call in
  `run_mcmc_reference.R`, `mcmc.R:77,113`).
- Ensure `adapt_from` is exposed as a CLI argument rather than hard-coded so the
  ablation can vary it.
- Note the interaction with warm-restart: `am_mcmc` gates adaptation on the
  local iteration counter `t`, while the Welford count `n` may carry over from
  `init_state`. Document and, if warm-restart is used, fix the gate to use `n`.
- Verify adaptation completes inside burn-in, and record the iteration at which
  the running covariance stabilizes (see Part B3).

## 5. Part B — Diagnostics to add

The current code checks only classic multivariate R-hat, pooled ESS, and mean
acceptance (the `coda::gelman.diag` / `coda::effectiveSize` block in
`run_mcmc_reference.R`). It does **not** check rank-normalized
R-hat, bulk/tail ESS, traceplots, per-chain consistency, covariance
conditioning, or ergodicity of the adaptation. Add all of the following.

### B1. Diagnose whether the NPE covariance is too narrow (requirement 4)

A small-scale NPE can return an over-confident (narrow) posterior covariance,
which then seeds a proposal that under-explores. Provide evidence for or against
this, using:

- **Scale comparison:** per-parameter posterior SD from the NPE draws vs the
  MCMC posterior SD; report the ratio. A ratio near 1 with good mixing argues
  against collapse; a small NPE/MCMC ratio indicates the proposal was too narrow.
- **Proposal vs target geometry:** compare `Sigma0 = cov(NPE)` with the MAP
  Hessian (Part D) and, if available, with the chain's converged empirical
  covariance.
- **Exploration test:** run short chains from dispersed starts and measure how
  far they drift relative to `Sigma0` and the prior scale; a proposal that is too
  narrow produces high acceptance with little movement and short ESS.
- **Acceptance/step-size:** record acceptance rate and proposal jump size
  (`||theta_prop - theta||`) over time.
- **Coverage cross-check:** compare the MCMC 95% interval with the NPE 95%
  interval and the SBC-based expected coverage. Flag parameters whose MCMC
  interval is much wider than the NPE interval.
- Report a per-parameter table with NPE SD, MCMC SD, ratio, acceptance, and ESS.

### B2. Convergence diagnostics (requirement 5)

Replace/augment the current diagnostics with the modern standard:

- **Rank-normalized split R-hat** (Vehtari et al. 2021) per parameter. Do not
  rely solely on `coda::gelman.diag` classical R-hat.
- **Bulk ESS and tail ESS** per parameter (rank-normalized), not only pooled
  ESS. Tail ESS is required for interval/quantile estimates.
- **Trace plots** per parameter per chain (write to CSV/PNG).
- **Acceptance rate** per chain and pooled; flag outside roughly 0.15–0.50.
- **Chain consistency:** overlapping posteriors across the 4 chains; report
  per-chain medians/SDs and R-hat; investigate any chain that is an outlier.
- Record MCSE for the mean and for the 2.5/97.5 quantiles.

Predeclared screening thresholds (state them before inspecting results):
rank-normalized split R-hat < 1.01, bulk ESS and tail ESS each > 400 (raise to
the project's existing `--ess-min` where relevant), and no single chain
diverging.

### B3. Ergodicity / Diminishing Adaptation (requirement 6)

The current code does nothing explicit for the ergodicity of the adaptive
scheme. Add a check based on the relative change of the proposal covariance:

- Track the proposal covariance `Sigma_t` at each adaptation step (or at a fixed
  recording stride) and compute the Frobenius-norm relative change between
  consecutive values,
  `drift_t = ||Sigma_t - Sigma_{t-1}||_F / max(||Sigma_{t-1}||_F, eps)`.
- Diminishing Adaptation requires `drift_t -> 0` as `t` grows. Report the drift
  trajectory, the iteration where it falls below a stated tolerance (e.g.
  1e-3), and the drift over the retained-sampling window.
- Confirm that adaptation has effectively stopped before retained sampling
  begins; if not, the retained draws come from a non-homogeneous chain and the
  usual R-hat/ESS interpretation is weakened.
- For warm-restart blocks, compute the drift across the block boundary as well.

### B4. Proposal covariance eigenvalues and condition number (requirement 9)

- Compute the eigenvalues of the active proposal covariance (initial `Sigma0`,
  the adapted running covariance, and the MAP Hessian covariance).
- Record `condition_number = max(eig) / min(eig)` for each. Larger values mean a
  more ill-conditioned proposal; annotate the eigen-spectrum plot.
- State a flag threshold (e.g. condition number > 1e6) and whether the ridge
  (`Sigma0 + diag(ridge)` in `run_mcmc_reference.R`, `eps` in `mcmc.R`) is
  masking a near-singular covariance.
- Compare the eigen-spectrum of `cov(NPE)` to the converged chain covariance to
  see which directions the NPE proposal under- or over-covers.

## 6. Part C — Initialization alternatives (requirement 7)

A good covariance does not guarantee a good starting point. Test initialization
separately from covariance:

- **Mean start (control):** all chains at the pooled NPE mean
  (`same_start = true`), as now.
- **Dispersed NPE-draw starts:** draw each chain's start from the NPE posterior
  samples (e.g. quantile-spread or random draws), including at least one
  overdispersed start. This tests whether the chain can move away from a
  potentially unrepresentative mean.
- **Multi-modal / tail starts:** include draws from the tails of each parameter
  to probe for missed modes.
- Report posterior recovery and ESS for each start set. If dispersed starts
  converge to the same posterior, initialization is not driving the result; if
  they do not, report missed modes or non-convergence rather than averaging over
  them.

## 7. Part D — Ablation study (requirement 8)

Compare at least the following arms under the same target, same prior, same
daily integration step, and a predeclared budget:

| Arm | Proposal covariance | Starting points |
| --- | --- | --- |
| A. NPE-initialized adaptive MH (control) | `cov(NPE)` then running empirical | NPE mean (and/or NPE draws) |
| B. Pilot MCMC | short pilot chain empirical covariance, then main chain | pilot end state / dispersed |
| C. MAP + Hessian | inverse Hessian (or its negative) at the deterministic MAP from `fit.rds` | MAP |
| D. (optional) Prior-covariance / random NPE draws | fixed/prior covariance | random NPE draws |

Requirements for a fair comparison:

- Same target density, same integration step, same prior, same RNG discipline.
- Report a forward-budget-matched view (equal ODE evaluations) and a
  wall-clock-matched view (include pilot and Hessian costs). Never compare a
  budget-exhausted weak arm against a well-tuned one without saying so.
- Metrics: rank-normalized R-hat, bulk/tail ESS, ESS per CPU-hour, MCSE,
  acceptance rate, covariance condition number, adaptation drift, and agreement
  of the recovered posterior with the other arms.
- For MAP + Hessian, state how the MAP was obtained (deterministic calibration
  best fit) and how the Hessian was approximated or computed, and note that a
  Hessian is a local Gaussian approximation.

## 8. Deliverables

Write everything to an isolated directory (e.g.
`output/ablation_check/<run-id>/`), with scripts kept in the repository and free
of machine-specific hard-coded paths:

- `run_manifest.json`: commit, file hashes, seeds, CLI args, hardware, per-arm
  budgets, actual ODE evaluations, runtimes, status.
- `diagnostics.csv`: per parameter per arm: rank-normalized split R-hat, bulk
  ESS, tail ESS, MCSE, acceptance, posterior SD.
- `covariance_eigen.csv`: per arm per block: eigenvalues, condition number.
- `covariance_drift.csv`: adaptation step, Frobenius relative change.
- Traceplot data/PNGs per parameter per chain.
- `start_sensitivity.csv`: results by initialization scheme.
- `ablation_report.md`: what was run, what was not, findings, and a decision.
- Raw chains for each arm (or a documented subsample if storage is prohibitive).

## 9. Acceptance and decision criteria

State these before inspecting results:

- All arms: rank-normalized split R-hat < 1.01 and bulk/tail ESS above the
  project threshold for every parameter.
- Adaptation drift (B3) falls below the stated tolerance before retained
  sampling; otherwise the arm is reported as non-homogeneous.
- Acceptance rate within 0.15–0.50.
- Condition number recorded for every proposal; near-singularity investigated
  rather than hidden by the ridge.
- Initialization sensitivity: dispersed starts must reach the same posterior
  within Monte Carlo error; otherwise report the sensitivity explicitly.
- The ablation must conclude which initialization/covariance strategy gives the
  best ESS per unit compute at acceptable accuracy, or state that evidence is
  insufficient.

Do not choose thresholds after seeing which arm wins. Report per-parameter
worst cases, not only averages.

## 10. Reporting requirements

- Clearly separate "current code behaviour" from "proposed change" from
  "measured result".
- Explicitly answer: does the current code check rank-normalized R-hat,
  bulk/tail ESS, traceplots, acceptance, cross-chain consistency (no — only
  classic R-hat, pooled ESS, mean acceptance)? Does it enforce Diminishing
  Adaptation (no — no covariance-drift check)?
- Label everything not executed as not executed. This plan is a specification,
  not evidence.
- Keep the daily-step cost visible: it dominates all other costs, so report
  wall-clock and ODE-call counts for every arm.

## 11. Key files

- `scripts/run_mcmc_reference.R` — sampler driver, priors, start points, diagnostics.
- `scripts/calibration/mcmc.R` — `make_priors`, `log_prior`, `log_posterior`, `am_mcmc`, Welford update.
- `scripts/calibration/likelihood.R` — target parameterization and bounds.
- `scripts/calibration/equilibrium.R` — equilibrium gate.
- `scripts/model/setup.R`, `scripts/model/sim.cpp` — model and integration step.
- `output/calibration/npe_agepop_recal/` — NPE posterior samples used for `Sigma0` and starts.
- `output/calibration/r_20260929_agepop/fit.rds` — deterministic MAP and prior anchors.
- `docs/NPE_vs_MCMC_benchmark.md` — related benchmark protocol and diagnostic standards.
