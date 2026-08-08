# AGENTS.md

## Scope

This repository implements an age-structured HCV compartmental model for
PWID in Singapore (3 strata D/J/X, 4 liver stages, 5 HCV states, 6 age
groups, 360 compartments; see `src/sim.cpp`). The calibration goal is to fit
the equilibrium prison (`J`) age-specific HCV seroprevalence and prison
population to the six supplied targets, with **particular attention to the
60+ age group**.

Model and parameter files:

```text
src/setup.R      parameters, initial conditions, helpers
src/sim.cpp      compiled ODE solver (3-strata model)
Model schematic.pptx   visual record — never modify
```

## Baseline specification (this revision)

1. **Liver-disease progression rates (highest defensible literature
   values)** — the following rates are written into `src/setup.R` and must
   not be silently lowered:

   | Parameter | Value (/yr) | Source |
   |---|---:|---|
   | `p_CC_DC` | 0.0788 | untreated-only decompensation rate in compensated HCV cirrhosis, Alazawi et al. 2010 systematic review, Aliment Pharmacol Ther 32(3):344-55, PMID 20497143 (untreated studies 7.88%/yr) |
   | `p_CC_HCC` | 0.0479 | untreated-only HCC rate in compensated HCV cirrhosis, same review (untreated studies 4.79%/yr) |
   | `p_DC_HCC` | 0.0464 | HCC incidence in HCV-monoinfected decompensated cirrhosis, 3.31/100 py, 95% CI 2.70-4.64 — upper bound used (Rivera-Irigoin et al. 2006, AIDS Res Hum Retroviruses 22(12):1236-41, PMID 17209765) |

   `p_NC_CC = 0.027` (Thein et al. 2008) is retained. Genotype-3 relative
   risks (`r3_*`, Kanwal et al. 2014) are retained; effective
   genotype-weighted rates are `(rho*r3 + (1-rho)) * p_*`. Record the
   effective rates and their interpretation in the natural-history review.

2. **No fitted `m_min`/`m_max`** — transmission is constant over time
   (`m_min = m_max = 1` in `setup.R`); the historical transmission level is
   absorbed into the fitted contact-matrix row scales. `m_min`, `m_max`,
   `m_t0`, `m_tau` are **not** fitted parameters in this revision.

3. **2015 Singapore baseline mortality** — SingStat age-specific death
   rates (2015) × scalar PWID SMR `omega = 14.68` (Mathers et al. 2013),
   as currently in `setup.R`. `mu_DC`, `mu_HCC`, `psi_DC`, `psi_HCC` are
   disease-specific and may be revised only with explicit evidence.

4. **Fitted parameter set (baseline: 12 parameters, no additions)**:

   ```text
   6 contact-matrix row scales  theta[1:6]   contact_scale = exp(theta)
   6 beta inflow scales         theta[7:12]  beta_scale    = exp(theta)
   ```

   No new parameters are added unless the 60+ age group cannot be fitted
   with these 12 (see "Excess mortality contingency" below).

5. **Soft rule: `beta_scale > 1` for every age group** (i.e. modelled annual
   inflow per age group exceeds the `setup.R` base `beta` values, which are
   anchored to CNB official new-drug-abuser data). Preference order:
   (a) fit with all `beta_scale > 1`; (b) if that materially degrades the
   fit, allow individual groups below 1 and document which groups and why;
   (c) never force the rule by degrading other criteria.

6. **Figures: ggplot2 only, publication quality.** Every figure produced by
   this project must be generated with `ggplot2` (no base-R plotting in
   calibration outputs). Use a consistent theme, colour-blind-safe palette,
   clear labels, and a `figures/` output.

## Calibration parameterisation

```r
contact_scale <- exp(theta_contact)   # theta[1:6], row-wise
beta_scale    <- exp(theta_beta)      # theta[7:12]
pm$C_contact[i, ] <- base_params$C_contact[i, ] * contact_scale[i]
pm$beta           <- base_params$beta * beta_scale
```

`m_min = m_max = 1`, `eta_s = rep(1, 6)` unless the excess-mortality
contingency is activated. Never write fitted `C_contact` or `beta` values
back into `setup.R`; keep them in run outputs.

## Targets

```r
age_groups       <- c("<20", "20-29", "30-39", "40-49", "50-59", "60+")
target_prev      <- c(0.1118421, 0.1731044, 0.2684954, 0.4301165,
                      0.4821029, 0.3544304)          # anti-HCV serology
target_prison_total <- c(99, 1244, 1467, 1841, 1628, 409)
```

The Changi Prison universal screening (2014-2016) is anti-HCV serology, so
the primary target is **seroprevalence** `(a+c+t+s)/N`; viremic prevalence
`(a+c+t)/N` is reported as a sensitivity. Integer positives:
`x_prev = round(target_prev * prison_total)`; record the rounding.

## Equilibrium and likelihood rules

- Do not call a transient endpoint equilibrium. Compare at `T` vs `T-5`:

  ```text
  max_i |log(N_hat_i(T)/N_hat_i(T-5))| <= 0.01
  max_i |p_hat_i(T) - p_hat_i(T-5)|    <= 0.005
  ```

  Also inspect total PWID population stability. Extend the horizon only in
  a local calibration copy.
- Prevalence likelihood: `x_prev[i] ~ Binomial(n_prev[i], p_hat[i])` with
  integer successes.
- Population observation model:
  `log(N_obs[i]) ~ Normal(log(N_hat[i]), sigma_pop^2)`, `sigma_pop = 0.10`
  (model-discrepancy tolerance).
- Objective: `NLL = nll_prev + nll_pop`. No arbitrary reweighting, no
  silent priors/regularisation.

## Optimisation rules

- Endpoint: deterministic multi-start **Nelder-Mead** (`optim(..., method =
  "Nelder-Mead")`). No HMC.
- Minimum starts: 1 with all log-scales = 0 + at least 5 reproducible
  perturbed starts (+ documented informed starts if justified).
- Record per start: initial/final theta, objective, convergence, function
  evaluations, elapsed time, equilibrium status, fit metrics.
- Re-run the best solution from nearby perturbations for local stability.
- Acceptance criteria:

  ```text
  prevalence RMSE                    <= 0.02
  maximum absolute prevalence error <= 0.03
  population MAPE                    <= 0.10
  maximum population APE            <= 0.20
  equilibrium: pass
  ```

  Also report Binomial deviance, population SRSS, total NLL, and Laplace
  95% intervals with observed-interval overlap per age group.

## Uncertainty quantification (Laplace)

After the best fit: finite-difference Hessian of the combined NLL at the
optimum, invert to `Sigma_theta` (generalized inverse with documented
eigenvalue cutoff; discard draws that fail the equilibrium gate), propagate
to the 12 target summaries (linearisation or Monte Carlo), and report 95%
intervals for `p_hat[i]` and `N_hat[i]` compared with observed intervals
(Binomial/Jeffreys for prevalence; log-Normal `sigma_pop = 0.10` for
population). Record overlap per age group.

## Excess-mortality contingency (60+ age group)

Fit with the 12-parameter set first. Only if the 60+ HCV prevalence cannot
be fitted with the 12 parameters (documented failure: acceptance criteria
missed specifically in the 60+ group across reasonable multi-start runs),
may an excess-mortality mechanism be added. Any addition requires a
dedicated, self-contained file (e.g.
`docs/calibration/excess_mortality_rationale.md`) stating: the mechanism,
the biological rationale and theory, the data sources (e.g. HCV-seropositive
vs -seronegative PWID mortality HRs; age-specific SMRs), the parameter
specification, and the fit consequence. Candidates (in order of preference):
age-specific seropositive mortality multiplier `eta_s[i]` (bounded,
evidence-based), then revisions to `mu_DC`/`mu_HCC` with evidence.

## Traceability and outputs

- Every run writes to `output/calibration/<run_id>/`:
  `run_config.csv`, `targets.csv`, `initial_values.csv`,
  `optimization_history.csv`, `solutions.csv`, `predictions.csv`,
  `residuals.csv`, `equilibrium.csv`, `laplace_intervals.csv`,
  `laplace_diagnostics.csv`, `sessionInfo.txt`, `fit.rds`, and ggplot2
  figures.
- Record Git SHA, hashes of `setup.R`, `sim.cpp`, `Model schematic.pptx`,
  R/package versions, targets, progression-rate specification, equilibrium
  settings, likelihood settings, optimizer settings, initial/fitted values,
  predictions, metrics.
- Append-only decision log: `docs/calibration/DECISIONS.md`. Every material
  change to mortality, progression, likelihood, parameterisation, optimizer
  settings, or scenario acceptance is recorded with date/time, change,
  reason, evidence, runs compared, fit consequence, decision. Never delete
  failed runs from the audit trail.
- Final comprehensive file `docs/calibration/PROJECT_OVERVIEW.md` (after the
  project completes) explains the full experiment sequence: what was tried,
  which run produced which result, and what each file represents.
- README.md / README.zh-CN.md at the repository root introduce structure,
  model, sources, and reproduction (English and Chinese).

## File-edit boundaries

- `Model schematic.pptx`: never modify.
- `src/sim.cpp`: do not modify unless a structural change is explicitly
  justified, documented, and validated (3-strata structure is the current
  baseline; ageing uses `y_change = y[base+i]/10.0*dt` — verify this during
  preflight).
- `src/setup.R`: the progression rates listed above are the selected
  baseline; `m_min = m_max = 1`; mortality spec as in Section 3. Any
  deviation requires a DECISIONS.md entry.
- Candidate variants (parameters/scenarios) live in
  `src/calibration/` and `output/calibration/`; only selected, justified
  variants are promoted into `setup.R`.

## Final validation checklist

1. `Model schematic.pptx` hash unchanged;
2. `src/sim.cpp` hash unchanged (or a documented, justified change);
3. `src/setup.R` diff limited to documented parameter/comments;
4. clean-session reproduction of the selected fit;
5. equilibrium confirmed (T vs T-5);
6. acceptance criteria met; Laplace intervals overlap observed intervals;
7. multi-start stability demonstrated;
8. structural limitations documented;
9. all figures ggplot2, publication quality;
10. `PROJECT_OVERVIEW.md` and READMEs present.

Stop and report rather than forcing a fit if targets are incompatible with
the model states, evidence is missing, fits remain poor across reasonable
multi-start runs, or acceptable fits require extreme/unstable parameters.
A valid result is: "the current model cannot defensibly meet the targets
under the allowed calibration scope."
