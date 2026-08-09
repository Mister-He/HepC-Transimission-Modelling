# prompt_mcmc.md

You are the Bayesian-inference lead for the HCV PWID compartmental model of
Singapore. Starting from the Nelder-Mead calibration (final run
`run2_v1_12p_warm`, NLL 22.44, all acceptance criteria met), obtain
**posterior samples of the 12 fitted parameters** and **Bayesian credible
intervals (95%) for the 12 target summaries**, while **keeping the
Laplace-approximated 95% intervals** as reported.

Read and obey the repository-root `AGENTS.md` (including its Bayesian-phase
section) before doing anything else.

## 1. Method priority (revised)

1. **Prefer the newest practical Bayesian methods — Neural Posterior
   Estimation (NPE / SNPE, e.g. Python `sbi`'s `NPE_C` with a normalizing
   flow) is the primary method.**
2. Only if NPE cannot be made to work (installation impossible, training
   fails, posterior demonstrably unreliable) fall back to a traditional
   sampler (e.g. adaptive Metropolis-Hastings) and document why.
3. Regardless of the primary method, a **traditional MCMC validation** is
   required as an independent cross-check, and it must satisfy the strict
   convergence criteria below.
4. Language is not limited to R: the simulation stays in R (Rcpp), the
   neural training/sampling runs in Python (`sbi`/`torch`), and all figures
   remain ggplot2 (generated in R from saved samples).

## 2. Strict convergence criteria (revised)

Whenever R-hat and ESS are computed (the MCMC validation, and any MCMC
used at all), they MUST satisfy:

```text
R-hat in [0.99, 1.01]   (every parameter, univariate; report mpsrf too)
ESS > 400               (per parameter; report pooled across chains and
                         per chain)
```

If the initial run fails these, **extend the chains** (e.g. warm-restart
from the last state and double the iterations) and re-check until they
pass; do not report a run that fails them. Record the final chain length.

For NPE (no R-hat/ESS in the sampling sense), the equivalent validation is
**required**:

- **SBC** (simulation-based calibration): posterior coverage of held-out
  prior samples (report per-parameter coverage at 90%/95% bands and the
  histogram of ranks);
- **multi-start stability**: at least two independent NPE trainings
  (different seeds / weight initialisations) must give consistent
  posteriors (report per-parameter median and 95% quantile agreement);
- **cross-check against the MCMC validation posterior** (quantile
  agreement per parameter; the MCMC posterior is the reference);
- a large posterior sample (>= 50,000 draws) so Monte Carlo error on the
  95% intervals is negligible (report the MC standard error of the 2.5%
  and 97.5% quantiles).

## 3. Scope and frozen baseline

- Model structure, assumptions, constants: **frozen** (`src/setup.R` +
  `src/sim.cpp`; progression 0.027 / 0.0788 / 0.0479 / 0.0464; constant
  m = 1; 2015 mortality x omega = 14.68; tau = 0). Do not change them.
- Parameterisation unchanged: 12 log-scale parameters (6 contact row
  scales, 6 beta inflow scales). No new fitted parameters.
- Likelihood unchanged (Binomial seroprevalence + log-Normal population,
  sigma_pop = 0.10). Equilibrium gate T vs T-5 is a strong model-validity
  constraint in the log-posterior / posterior-predictive filter.

## 4. Priors (documented, model- and literature-based)

| Parameter block | Prior | Rationale |
|---|---|---|
| log contact scales `theta[1:6]` | `Normal(0, 2^2)` | Base `C_contact` is a model-specific calibrated guess; weakly informative. 95% prior interval of the scale factor ~ [0.02, 55]. |
| log beta scales `theta[7:12]` | `Student-t(3, 0, 2)` | Base `beta` is CNB-official-anchored but per-age values are placeholders; heavy tails tolerate the data-driven 60+ inflow (~110x). |

Do not construct priors from the NM fit (no double counting); NM solutions
serve only as chain starts / proposal tuning. Document any deviation.

## 5. NPE implementation (primary)

1. **Training data** (R): sample `theta ~ prior` (>= 30,000 draws),
   simulate `x = (p_1..p_6, N_1..N_6)` at the 2015 target time; keep ALL
   simulations (equilibrium is enforced at the predictive stage, not by
   deleting training data). Store `theta_train.csv` / `x_train.csv`.
2. **Train** (Python `sbi`): `NPE_C` with a masked autoregressive flow
   (or similar); 2-3 rounds (sequential proposal refinement); report
   training loss curves and the number of simulations per round.
3. **Observed summary**: `x_obs = (prev_binom, prison_total)`.
4. **Posterior sample**: >= 50,000 draws from the final posterior.
5. **SBC**: on held-out prior samples, compute posterior ranks under the
   trained flow; report coverage diagnostics.
6. Save `posterior_samples.csv`, `sbc_summary.csv`, loss curves, and the
   per-round config.

## 6. MCMC validation (reference, strict criteria)

- Sampler: adaptive Metropolis-Hastings (Haario et al. 2001); proposal
  covariance initialised from the **NPE posterior covariance** (better
  mixing than the Laplace covariance), or the Laplace covariance as
  fallback.
- **At least 3 chains from different starts** (e.g. NM best, NPE posterior
  median, another NM start). Extend iterations until
  **R-hat in [0.99, 1.01]** for every parameter and **ESS > 400**
  (pooled), then stop and record the length used.
- Report acceptance rate, trace plots (all chains overlaid), marginal
  density plots, posterior median vs NM optimum.

## 7. Posterior predictive credible intervals

- For retained NPE posterior draws and MCMC draws: simulate, keep
  equilibrium-feasible draws, compute the 12 target summaries; report
  equal-tailed 95% CrI + medians per age group.
- Compare **NPE CrI vs MCMC CrI vs Laplace 95% CI vs observed intervals**
  (Jeffreys Binomial for prevalence; log-Normal sigma_pop = 0.10 for
  population) in one table per target with explicit overlap.

## 8. Outputs

`output/calibration/npe_bayes/`:

```text
run_config.csv            priors, rounds, sims/round, seeds, chain settings
theta_train.csv, x_train.csv
posterior_samples_npe.csv     >= 50,000 draws
posterior_samples_mcmc.csv    MCMC validation draws
sbc_summary.csv               coverage diagnostics
diagnostics_mcmc.csv          R-hat (must be in [0.99, 1.01]), ESS (>400),
                              acceptance, chain length used
credible_intervals.csv        NPE + MCMC + Laplace + observed, overlaps
target_posterior.csv          predictive draws (p/N per age)
sessionInfo.txt (R), npe_config.json (Python)
ggplot2 figures: fig_trace, fig_density, fig_predictive_density,
  fig_ci_compare, fig_sbc
```

## 9. Documentation deliverables

1. `docs/calibration/bayes_methodology.md` — theory (NPE/SNPE, normalizing
   flows, SBC; Bayesian inference; relation to Laplace; why NPE + MCMC
   validation), process (data generation, training rounds, seeds, runtimes,
   chain lengths used to meet the strict criteria), and results
   interpretation (posterior medians, intervals, overlaps, 60+ and
   beta_scale implications, NPE-vs-MCMC agreement).
2. Append to `docs/calibration/DECISIONS.md` (including the deletion of
   the previous MCMC results and the method switch to NPE).
3. Update `docs/calibration/final_report.md`, `PROJECT_OVERVIEW.md`,
   `README.md`, `README.zh-CN.md`.

## 10. Final validation

- `Model schematic.pptx`, `src/sim.cpp`, `src/setup.R` unchanged;
- NPE trained and posterior sampled; SBC diagnostics reported;
- MCMC validation met R-hat in [0.99, 1.01] and ESS > 400 (recorded);
- NPE and MCMC posteriors agree (quantile agreement reported);
- credible intervals overlap the observed intervals per age group;
- Laplace intervals retained and compared;
- all figures ggplot2.
