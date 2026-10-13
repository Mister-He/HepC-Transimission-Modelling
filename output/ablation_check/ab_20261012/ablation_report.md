# Reference-MCMC ablation report (run `ab_20261012`)

This run implements `docs/ablation_check.md` on 2026-10-12 and 2026-10-13. Each section separates
**current code behaviour**, the **change made** and the **measured result**. All numbers come from
the CSV files in this folder.
- **Settings:** fixed in advance in `ablation_config.json`, written before any arm ran.
- **Driver:** `scripts/ablation_check.R`, run from commit `8ce1cb9`.
- **Provenance:** `run_manifest.json` records commands, seeds, file hashes and per-chain budgets.

## 0. Baseline: what the plan assumed versus what the code was

The plan's baseline (§3) describes an older workflow: NPE-initialised `am_mcmc` with a monthly step,
`coda` diagnostics, and `npe_agepop_recal` / `r_20260929_agepop` inputs. That workflow and those
folders had already been removed.

When the ablation started, the reference MCMC was `scripts/run_mcmc_reference.R` with
`arwm_budget`:
- daily step;
- 4 chains started from prior draws;
- initial proposal 0.1² × the prior variance;
- adaptation from proposal 200, refreshed every 50, and frozen after 4,999 warmup proposals;
- `thin = 1`.

**The plan's §10 questions, answered for the code at that point:**
- **Already present:** rank-normalised split R-hat, bulk and tail ESS, and MCSE of the mean and
  2.5/97.5% quantiles.
- **Missing:** trace plots, per-chain consistency summaries, per-chain acceptance flags, any check on
  covariance drift or eigenvalues, and configurable starts or covariance.

**Diminishing adaptation.** This was not checked explicitly. The kernel was frozen after warmup, so
retained draws already came from a fixed kernel.

**Inputs used here:**
- **NPE:** the current accepted posterior, `output/calibration/npe_lhs`.
- **Fit:** `output/calibration/r_20261010_agepop/fit.rds`.
- **Prior:** Normal(`fit$best$theta`, 0.20) on all 12 log-scale parameters, the same as the NPE.

## 1. Part A: configuration changes

### A1. Daily integration step

- **Current:** already `dt = 1/365`. The monthly step (0.0833) was used only by the removed
  `mcmc_recal` run.
- **Cost:** 0.29–0.32 CPU-seconds per daily evaluation, the 150-year run plus the equilibrium check.
  A monthly evaluation is about 18 times cheaper: 2.15 against 0.12 CPU-hours for the same short
  chains.
- **Measured, deterministic** (`dt_deterministic_comparison.csv`, 500 reference posterior draws):
  - **Prevalence:** the two steps differ by at most 0.18 percentage points (mean 0.02–0.13).
  - **Prison populations:** they differ by at most 6.8 people, under 0.7%.
  - **Log posterior:** the mean absolute difference is 0.11 and the maximum 0.47.
- **Measured, short chains** (`dt_chain_comparison.csv`, 4 chains × 6,000 evaluations each):
  - posterior means differ by at most 0.11 SD, within 1.3 combined MCSE;
  - these chains are short (R-hat up to 1.03, bulk ESS from 236), so this comparison is weak
    evidence.

  The discretisation shifts the posterior by much less than its width.
- **Decision:** keep the daily step, as the plan requires.

### A2. Thinning

- **Current:** `thin = 1`. The plan suggested 1,000.
- **Measured** on the 180,000 reference draws (`thinning_evidence.csv`):
  - the integrated autocorrelation time is 39–53 iterations;
  - thinning by 1,000 leaves 180 draws, with tail ESS 91 and R-hat 1.036, which **fails** the checks;
  - thinning by 10 keeps 98% of bulk ESS and 95% of tail ESS.
- **Change:** `--thin` now applies to storage only, with **default 10**, chosen by the rule fixed in
  advance: the largest interval keeping at least 90% of bulk and tail ESS. All diagnostics use every
  retained draw. The pilot and step-comparison arms store unthinned draws.

### A3. Adaptation start

- **Current:** adaptation started at proposal 200 and was not a command-line option.
- **Change:**
  - `--adapt-start` and `--adapt-every` are now options;
  - the default is `--adapt-start 5000` with warmup raised to 10,000, because with the old warmup of
    4,999 the proposal would never adapt;
  - the script warns if `--adapt-start` is at least `--n-warmup`.
- **Warm restart:** `arwm_budget` has no warm-restart path, so the plan's concern about the counter
  does not apply.
- **Measured:** the relative change of the adapted covariance first fell below 0.01 at proposals
  5,351–6,051, depending on the arm (`arm_summary.csv`).

## 2. Part B: diagnostics added

### B1. Is the NPE covariance too narrow?

From `npe_vs_mcmc_scale.csv`:
- **Spread:** NPE posterior SD / reference MCMC SD is 0.93–0.99 (mean 0.957), and the 95% width
  ratio is 0.94–0.99. So the NPE is mildly narrow but **not collapsed**.
- **Hessian:** inverse-Hessian SDs at the MAP match the MCMC SDs within about 2%, so the posterior is
  close to Gaussian.
- **Sampler behaviour with the NPE covariance** (arm A_mean):
  - acceptance 0.25–0.26;
  - mean jump size 0.30 at the start and 0.32 at the end of warmup;
  - bulk ESS 1,682–2,219.

  The NPE-based proposal does not under-explore.
- **Starts from the NPE tails** (A_tail) and **overdispersed starts** (A_disp) recover the reference
  posterior.
- **Conclusion:** the NPE covariance is slightly narrow, about 4%, which has no practical effect on
  mixing.

### B2. Convergence

`diagnostics.csv` gives per-parameter results for every arm. The thresholds, set in advance, are
R-hat < 1.01, bulk and tail ESS > 400, acceptance 0.15–0.50, and per-chain medians within 0.5 SD.

| Arm | Max R-hat | Min bulk ESS | Min tail ESS | Acceptance | Max chain median shift (SD) |
| --- | --- | --- | --- | --- | --- |
| A_mean | 1.004 | 1,682 | 3,101 | 0.25–0.26 | 0.10 |
| A_disp | 1.003 | 1,685 | 2,983 | 0.25–0.27 | 0.10 |
| A_tail | 1.005 | 1,846 | 3,197 | 0.25–0.26 | 0.09 |
| B_pilot | 1.005 | 1,494 | 2,810 | 0.26 | 0.11 |
| C_map | 1.005 | 1,797 | 2,870 | 0.24–0.26 | 0.13 |
| D_prior | 1.005 | 1,422 | 2,319 | 0.27–0.28 | 0.12 |

- **All checks:** every arm passes every one of these checks.
- **Trace plots:** in `<arm>/mcmc_traceplots.png`. They show stationary, overlapping chains.
- **Other outputs:** MCSE for the mean and the 2.5/97.5% quantiles is in `diagnostics.csv`, and
  per-chain summaries are in `<arm>/mcmc_per_chain.csv`.

### B3. Adaptation drift

Files: `covariance_drift.csv` and `<arm>/retained_block_covariance_drift.csv`.

- **Change made:** `arwm_budget` now records, at every adaptation update, the relative Frobenius
  change from the previous adapted covariance and its condition number.
- **Measured:** the change between the last two updates before freezing was **0.011–0.021**. That is
  **above the 0.01 tolerance** set in advance, so every arm formally **fails** this check
  (`pass_drift = FALSE`).
  - It fell below 0.01 at proposals 5,351–6,051 but kept fluctuating around 0.01–0.02 afterwards.
  - It never fell below 10⁻³.
  - At this point each update adds 50 draws of autocorrelated output (IAT about 40) to a window of
    about 10,000 iterations. The residual change is estimation noise in the window covariance, not a
    trend.
- **Interpretation:** the proposal is frozen at the end of warmup, so every retained draw comes from
  a single fixed Metropolis kernel. Diminishing adaptation is enforced by freezing. The 0.01
  tolerance was too strict for 50-proposal updates; the failure reflects the threshold, not
  non-stationarity.
- **Retained window:** covariances of consecutive blocks of about 2,000 retained iterations differ by
  0.20–0.27 in relative Frobenius norm. That is the sampling noise of a covariance estimated from
  about 50 effective draws per block, not adaptation. The last block is a short remainder.
- **Suggested change for a future run** (not tested): either refresh every 500 proposals, or use a
  tolerance that scales with the number of effective draws added per update.

### B4. Eigenvalues and condition numbers

File: `covariance_eigen.csv`.

| Matrix | Condition number | Smallest eigenvalue |
| --- | --- | --- |
| cov(NPE) | 103 | 3.4 × 10⁻⁴ |
| inverse Hessian at MAP | 106 | 3.7 × 10⁻⁴ |
| pilot covariance | 111 | 3.5 × 10⁻⁴ |
| reference posterior covariance | 96 | 4.0 × 10⁻⁴ |
| adapted proposals (all arms) | 105–125 | — |

- **No near-singular matrices.** Every condition number is far below the 10⁶ flag.
- **The ridge masks nothing.** The regularisation of 10⁻⁶ × prior variance (about 4 × 10⁻⁸) is four
  orders of magnitude below the smallest eigenvalue.
- **Prior-based start:** the initial proposal of arm D_prior is diagonal, with condition number 1.1.

## 3. Part C: initialisation

File: `start_sensitivity.csv`. The test is |mean_arm − mean_reference| ≤ 3 × combined MCSE, for every
parameter.

| Starts | Worst difference (SD) | Worst difference (combined MCSE) | Same posterior? |
| --- | --- | --- | --- |
| NPE mean, all chains (A_mean) | 0.045 | 1.64 | yes |
| 3 random NPE draws + 1 overdispersed, 3 SD out (A_disp) | 0.053 | 1.90 | yes |
| 4 NPE draws from the top 0.5% Mahalanobis tail (A_tail) | 0.050 | 1.80 | yes |
| 4 random NPE draws, prior-scaled proposal (D_prior) | 0.045 | 1.68 | yes |

- **All starts were valid:** every start evaluated to a finite posterior on the first attempt.
- **Conclusion:** dispersed and tail starts reach the same posterior. Initialisation does not drive
  the result, and no second mode was found.

## 4. Part D: covariance ablation

Every arm used the same target, the daily step, warmup of 10,000 with adaptation from 5,000, and a
budget of 30,000 evaluations per chain.
- **B_pilot:** the 3,000-evaluation pilot per chain was taken out of that budget, leaving 27,000.
- **C_map:** the Hessian took 314 evaluations (80 s), reported separately.

| Arm | Initial proposal | Evaluations (main + extra) | CPU-hours | Min bulk ESS | Min bulk ESS per CPU-hour |
| --- | --- | --- | --- | --- | --- |
| A_mean | (2.38²/d) cov(NPE) | 120,000 + 0 | 11.8 | 1,682 | 142 |
| A_disp | as A | 120,000 + 0 | 11.9 | 1,685 | 142 |
| A_tail | as A | 120,000 + 0 | 11.8 | 1,846 | 156 |
| B_pilot | (2.38²/d) pilot covariance | 108,000 + 12,000 | 9.7 | 1,494 | 155 |
| C_map | inverse Hessian (unscaled) | 120,000 + 314 | 9.7 | 1,797 | 185 |
| D_prior | 0.1² diag(prior variance) | 120,000 + 0 | 11.9 | 1,422 | 120 |

**Two ways to compare cost:**
- **Same number of evaluations** (all 120,000): ranked by minimum bulk ESS, A_tail > C_map > A_disp ≈
  A_mean > B_pilot > D_prior.
- **Per CPU-hour:** C_map is highest. But wave 2 (B_pilot, C_map) shared the machine with the cheap
  monthly chains and so ran with less contention: about 0.29 CPU-seconds per evaluation, against
  about 0.35 in wave 1. CPU-time comparisons between waves are therefore biased in favour of B_pilot
  and C_map.

**Early warmup behaviour** (`mcmc_jumps.csv`, first 250 proposals):
- **D_prior:** the prior-scaled proposal is small (mean jump 0.07, acceptance 0.58) until adaptation
  starts at proposal 5,000.
- **C_map:** the unscaled inverse Hessian is large (jump 0.46, acceptance 0.10). Scaling it by
  2.38²/d, as for A and B, would likely have helped; this run did not test that.

**Decision.** At the same evaluation budget, any posterior-shaped initial proposal (NPE, pilot or
Hessian) gives about 5–30% more minimum bulk ESS than the prior-scaled proposal (D_prior). Among the
posterior-shaped proposals, the differences (1,494–1,846) are within the noise of an ESS estimate
from one 4-chain run, about ±10%. Each arm ran once.

- **Strongest supported ranking:** the evidence is **insufficient to rank** NPE, pilot and Hessian
  initialisation. It is sufficient to prefer any of them over the prior-scaled start.
- **NPE initialisation:** it costs nothing extra once the NPE exists and gave the best result at the
  same evaluation budget (A_tail).
- **MAP Hessian:** it costs only 314 evaluations and is a good alternative that needs no NPE.
- **Production posterior:** it remains valid. Arm D_prior uses its covariance setting, and every arm
  agrees with it within Monte Carlo error.

## 5. What was not run

- **Classical multivariate R-hat (`coda::gelman.diag`):** replaced by rank-normalised split R-hat,
  as the plan recommends.
- **Unthinned storage:** only the pilot and the daily/monthly step chains stored all draws. The
  control arm's diagnostics use all draws in memory, but its stored draws are thinned by 10.
- **Repeated replicates of each arm:** each arm ran once (4 chains), so ESS differences between
  arms carry about ±10% Monte Carlo uncertainty.
- **A scaled-Hessian arm and a longer adaptation refresh interval:** suggested above, not run.
- **SBC-based expected coverage (B1):** not available, because the LHS NPE was not trained with SBC.
