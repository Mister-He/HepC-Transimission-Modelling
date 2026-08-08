# prompt.md

You are the lead calibration agent for an age-structured HCV compartmental
model of PWID in Singapore (3 strata D/J/X, 4 liver stages, 5 HCV states,
6 age groups, 360 compartments). Your job is to **fit the model at
equilibrium to the prison (`J`) age-specific HCV seroprevalence and
population targets**, with particular attention to the 60+ age group, and
to keep every change justifiable and documented.

Read and obey the repository-root `AGENTS.md` before doing anything else.

## Objective

```r
age_groups <- c("<20", "20-29", "30-39", "40-49", "50-59", "60+")
target_prev <- c(0.1118421, 0.1731044, 0.2684954, 0.4301165,
                 0.4821029, 0.3544304)
target_prison_total <- c(99, 1244, 1467, 1841, 1628, 409)
```

Fit every age group well, especially 60+ (both prevalence and prison
population).

## Baseline specification (this revision)

1. **Progression rates — use the highest defensible literature values
   (already selected, do not lower):**

   | Parameter | Value (/yr) | Source |
   |---|---:|---|
   | `p_CC_DC` | 0.0788 | Alazawi et al. 2010 systematic review (PMID 20497143), untreated-only decompensation in compensated HCV cirrhosis: 7.88%/yr |
   | `p_CC_HCC` | 0.0479 | Same review, untreated-only HCC: 4.79%/yr |
   | `p_DC_HCC` | 0.0464 | Rivera-Irigoin et al. 2006 (PMID 17209765): HCC in HCV-monoinfected decompensated cirrhosis 3.31/100 py, 95% CI 2.70-4.64; upper bound used |

   `p_NC_CC = 0.027` (Thein et al. 2008). GT3 relative risks (Kanwal et al.
   2014) retained. Cite sources in `setup.R` comments and in
   `docs/calibration/natural_history_review.md`.

2. **No `m_min`/`m_max` fitting** — transmission is constant
   (`m_min = m_max = 1`), merged into the contact-matrix row scales.

3. **2015 Singapore baseline mortality** (SingStat 2015 x PWID SMR
   omega = 14.68, Mathers et al. 2013).

4. **12 fitted parameters** (no additions unless the 60+ group cannot be
   fitted):

   ```text
   6 contact-matrix row scales (exp)
   6 beta inflow scales (exp)
   ```

5. **Soft rule: `beta_scale > 1` per age group.** Try to satisfy it; if it
   materially degrades the fit, drop it for the affected groups and
   document why. Never force it at the cost of the acceptance criteria.

6. **All figures with `ggplot2`, publication quality** (consistent theme,
   colour-blind-safe palette, clear labels, PNG output in `figures/`).

## Work plan

### Phase 0 — Preflight

1. Read `AGENTS.md`, `src/setup.R`, `src/sim.cpp`, `Model schematic.pptx`.
2. Record Git SHA, `git status --short`, MD5 hashes of the three files.
3. Verify the `/10.0` ageing implementation in `src/sim.cpp`.
4. Create `docs/calibration/`, `src/calibration/`, `output/calibration/`,
   `figures/`.
5. Run the baseline simulation and record its fit (no-change reference).

### Phase 1 — Model audit

Write `docs/calibration/model_audit.md`: compartment structure, ageing,
target definitions (sero vs viremic), J no-in-prison transmission, and a
diagnosis of where the age pattern comes from under constant transmission
(12-parameter fit, no time multiplier, no excess 60+ mortality yet).

### Phase 2 — Mortality and natural-history evidence

- `docs/calibration/mortality_review.md`: 2015 SingStat rates mapped to the
  six age groups; PWID SMR 14.68 (Mathers et al. 2013); source records.
- `docs/calibration/natural_history_review.md`: the adopted progression
  rates, their sources and CIs, genotype-weighted effective rates, and
  caveats (age-dependence; rates vs probabilities).

### Phase 3 — Calibration harness (`src/calibration/`)

Create:

```text
targets.R        target data + Binomial count reconstruction
model_metrics.R  J summary extraction + fit metrics
equilibrium.R    T vs T-5 stability gate
likelihood.R     parameterisation (12 params) + NLL objective
calibrate_nm.R   multi-start Nelder-Mead runner
laplace.R        Laplace approximation + interval propagation
run_calibration.R  end-to-end runner (metadata, multi-start, outputs)
plot_results.R   ggplot2 publication figures (ALL figures)
```

Requirements: reproducible seeds; failure -> large finite objective; every
evaluation logged; settings/outputs saved per run; baseline flag recording
whether the fitted model matches the committed model.

### Phase 4 — Equilibrium

`max_i |log(N_hat_i(T)/N_hat_i(T-5))| <= 0.01` and
`max_i |p_hat_i(T)-p_hat_i(T-5)| <= 0.005`, plus total-PWID stability.
Extend the horizon only in a local calibration copy (document it).

### Phase 5 — Likelihood

`x_prev[i] ~ Binomial(n_prev[i], p_hat[i])` (integer successes);
`log(N_obs[i]) ~ Normal(log(N_hat[i]), 0.10^2)`; `NLL = nll_prev +
nll_pop`. No reweighting, no priors.

### Phase 6 — Multi-start Nelder-Mead

At least 6 deterministic starts (baseline zeros + 5 reproducible
perturbations; documented informed starts optional). Record everything per
start. Re-run the best from nearby perturbations.

### Phase 7 — Iterative decision loop

In order: implementation correctness; equilibrium; initial values; mortality
scenarios; fitted-parameter magnitude and residual pattern; beta_scale>1
feasibility (document which groups, if any, remain below 1 and why); and
only then, if the 60+ group still cannot be fitted, the excess-mortality
contingency with a dedicated rationale file. Every material change is
appended to `docs/calibration/DECISIONS.md` (never delete failed attempts).

### Phase 8 — Acceptance criteria

```text
equilibrium: pass
prevalence RMSE                    <= 0.02
maximum absolute prevalence error <= 0.03
population MAPE                    <= 0.10
maximum population APE            <= 0.20
```

Report Binomial deviance, population SRSS, total NLL, and Laplace 95%
intervals with observed-interval overlap for every age group.

### Phase 9 — Uncertainty quantification (Laplace)

Hessian of the combined NLL at the optimum; covariance (documented
eigenvalue cutoff); propagate to the 12 summaries (Monte Carlo or
linearisation); discard equilibrium-infeasible draws; report 95% intervals
vs observed intervals (Binomial/Jeffreys; log-Normal sigma_pop = 0.10);
state overlap per age group.

### Phase 10 — Final deliverables

1. Promote selected parameter/structural changes into `setup.R` (documented);
2. reproduce from a clean R session;
3. write `docs/calibration/final_report.md` (specification, likelihood,
   equilibrium, fitted scales, predictions + intervals, metrics, multi-start
   stability, sensitivity, remaining limitations);
4. **`docs/calibration/PROJECT_OVERVIEW.md`** — the comprehensive file that
   walks through the whole project from start to finish: every experiment
   tried and its outcome, what each file represents, and which run produced
   which result;
5. regenerate `README.md` and `README.zh-CN.md` (structure, model, sources,
   reproduction);
6. all figures in `figures/` via ggplot2.

## Change-trace requirements

Before work: `git status --short`, `git rev-parse HEAD`. Record input
hashes and Git SHA in run metadata. Use logical commits where appropriate.
At the end show `git status --short`, `git diff -- src/setup.R`,
`git diff -- src/sim.cpp`, and confirm `Model schematic.pptx` hash
unchanged.

The final response must list:

1. files created;
2. files modified;
3. protected-file verification (`Model schematic.pptx`);
4. best calibration metrics, including Laplace interval overlap;
5. selected mortality (2015 baseline) and progression-rate specification;
6. beta_scale>1 outcome per age group;
7. excess-mortality decision (used or not, and why);
8. remaining model limitations;
9. exact reproduction commands.

## First actions

1. read `AGENTS.md`;
2. record Git state and hashes;
3. read the model files;
4. diagnose the 12-parameter baseline fit;
5. research/record 2015 mortality and the adopted progression rates;
6. implement and run calibration.
