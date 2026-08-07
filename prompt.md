# prompt.md

You are the lead calibration agent for an HCV compartmental model of PWID in Singapore.

Your job is to **research the mortality specification, audit the protected model, design and implement an equilibrium calibration likelihood, run reproducible multi-start Nelder-Mead calibration, iteratively evaluate evidence-based mortality scenarios and initial values, and leave a complete audit trail**.

Read and obey the repository-root `AGENTS.md` before doing anything else.

## Objective

Calibrate the model at equilibrium to these six prison age groups:

```r
age_groups <- c("<20", "20-29", "30-39", "40-49", "50-59", "60+")

target_prev <- c(
  0.1118421,
  0.1731044,
  0.2684954,
  0.4301165,
  0.4821029,
  0.3544304
)

target_prison_total <- c(
  99,
  1244,
  1467,
  1841,
  1628,
  409
)
```

Fit 12 positive parameters:

- 6 contact-matrix row scaling factors;
- 6 age-specific inflow (`beta`) scaling factors.

Use R:

```r
optim(..., method = "Nelder-Mead")
```

Nelder-Mead is the calibration endpoint. **Do not run HMC or create a Nelder-Mead -> HMC two-stage workflow.**

---

## Non-negotiable file constraints

### Never modify

```text
src/sim.cpp
Model schematic.pptx
```

### `src/setup.R`

Read the whole file, but modify only the final selected baseline mortality specification and its adjacent documentation:

```text
mu
omega
```

Do not change any other model parameter or design in `setup.R`.

Do not change `mu_DC`, `mu_HCC`, `psi_DC`, or `psi_HCC` in this task.

During iteration, keep candidate mortality assumptions in new calibration files rather than repeatedly editing `setup.R`.

Before finalizing, prove with `git diff` that these rules were obeyed.

---

## Important project parameterization to preserve

Use log-parameterized positive scaling factors.

For contact matrix row `i`:

```r
pm$C_contact[i, ] <- base_params$C_contact[i, ] * exp(theta_contact[i])
```

For inflow:

```r
beta_scale <- exp(theta_beta)
pm$beta <- base_params$beta * beta_scale
```

Preserve the current calibration convention linking inflow scaling to first arrest:

```r
c_true <- base_params$c_composite / beta_scale
pm$lambda1 <- base_params$lambda3 * c_true
```

Do not silently decouple this relationship.

---

## Work plan

### Phase 0 — Preflight and immutable baseline

1. Read:
   - `AGENTS.md`
   - `src/setup.R`
   - `src/sim.cpp`
   - `Model schematic.pptx`
2. Inspect repository status.
3. Record:
   - current Git SHA;
   - `git status --short`;
   - file hashes for `src/setup.R`, `src/sim.cpp`, and `Model schematic.pptx`.
4. Do not discard or overwrite unrelated changes.
5. Create the required calibration/docs directories if absent.

Write the baseline metadata into the first calibration run directory.

---

### Phase 1 — Model audit before fitting

Create:

```text
docs/calibration/model_audit.md
```

Derive the calibration target directly from the implemented model.

For prison age group `i`:

```text
N_hat[i] = total population in all J compartments for age i
```

For current/viremic HCV prevalence:

```text
I_hat[i] = J acute + J chronic + J treatment, summed across liver stages
p_hat[i] = I_hat[i] / N_hat[i]
```

Verify the interpretation of the prevalence target if a source is available.

Explicitly audit and report at least these issues:

1. `J` has no new infection force in the current C++ implementation.
2. background mortality is implemented as `mu[i] * omega`.
3. `omega` is scalar.
4. verify that the current repository contains the already-corrected ageing
   implementation:

   ```cpp
   y_change = y[base + i] / 10.0 * dt;
   ```

   Treat `/10.0` as the protected baseline. Do not revert it or further modify
   ageing. If the repository unexpectedly contains `/5.0`, stop and report the
   discrepancy rather than editing `src/sim.cpp`.
5. state `u` combines susceptible and post-SVR; therefore anti-HCV
   seroprevalence cannot be reconstructed exactly from this state partition.
6. identify any additional mismatch that could materially distort equilibrium
   age structure or HCV prevalence.

Do not patch protected structural issues.

---

### Phase 2 — Mortality evidence review

Create:

```text
docs/calibration/mortality_review.md
src/calibration/mortality_scenarios.R
```

Research two components separately.

#### A. Singapore baseline mortality

Use authoritative Singapore sources first, preferably SingStat age-specific death rates or complete life tables.

The model calibration appears intended for 2017. Therefore:

- use mortality aligned to 2017 for the primary specification when possible;
- treat newer Singapore mortality as sensitivity analysis unless there is a strong reason to remap the target year.

For each mortality source record:

- exact source and URL;
- publication/table name;
- reference year;
- resident/sex definition;
- original age bands;
- reported mortality quantity;
- conversion to annual model rate;
- mapping to the six model age groups.

Do not confuse:

- probability of death over an interval;
- age-specific death rate;
- instantaneous hazard.

Show the formula used for any conversion.

#### B. PWID excess mortality / SMR

Search peer-reviewed PWID mortality literature and systematic reviews.

The current scalar SMR is approximately the pooled magnitude reported by the Mathers et al. PWID mortality meta-analysis, but do not assume that one pooled value is appropriate for Singapore or for all ages.

Prioritize evidence in this order:

1. credible age-stratified PWID SMRs, even if they use fewer than six age bands;
2. Asian PWID cohorts with usable mortality comparisons;
3. well-described age-stratified PWID cohorts from other settings;
4. systematic reviews/meta-analyses for context.

Do **not** require one study to provide all six model age groups. Credible
two-, three-, four-, or five-band SMR evidence can be used.

When reported SMR bands are coarser than the model age groups:

- map them transparently to the six model age groups;
- default to piecewise-constant mapping when a model age group lies wholly
  inside a reported SMR band;
- do not invent smooth interpolation;
- if boundaries overlap imperfectly, document the assumption;
- use weighted aggregation only when defensible population weights are
  available;
- preserve the original reported strata and mapped six-group vector in output.

For each candidate source extract:

- country/setting;
- cohort and years;
- age groups;
- sex distribution;
- HIV/HCV context;
- deaths/person-years if available;
- SMR/CMR and 95% CI;
- reference population and standardization;
- relevance to Singapore;
- limitations.

Do not invent missing age-specific values.

#### Required mortality scenarios

Create only evidence-supported scenarios, with at least:

```text
M0_current
M1_refreshed_singapore_mu_plus_scalar_smr
```

If credible age-stratified SMR evidence is sufficient—even when it has fewer
than six reported age bands—also create:

```text
M2_age_stratified_effective_mortality
```

First map the reported SMR strata to a six-element model-age vector using the
documented piecewise/coarse-band mapping rule.

Because `sim.cpp` only accepts scalar `omega`, implement the mapped
age-specific SMR without changing C++:

```r
effective_mu <- singapore_mu * smr_model_age
omega <- 1
```

Keep all of the following separately in the scenario object and output tables:

```text
original literature SMR strata
mapping rule
mapped six-group smr_model_age
singapore_mu
effective_mu
```

Do not edit `setup.R` yet.

---

### Phase 3 — Implement calibration harness

Create modular R files under:

```text
src/calibration/
```

At minimum:

```text
targets.R
model_metrics.R
equilibrium.R
likelihood.R
calibrate_nm.R
run_calibration.R
mortality_scenarios.R
```

Requirements:

1. source the existing `src/setup.R`;
2. treat its loaded parameter list as immutable baseline input inside the calibration harness;
3. build candidate parameter lists with `modifyList()` / copies;
4. never store fitted contact/beta values back into `setup.R`;
5. use the existing Rcpp simulation unmodified;
6. make all runs reproducible with an explicit seed;
7. handle simulation failures by returning a large finite objective value;
8. save every run's settings and outputs.

---

### Phase 4 — Equilibrium check

Do not assume the final simulation row is equilibrium.

For each candidate run, compute target summaries at `T` and approximately `T - 5`.

Initial equilibrium criteria:

```text
max abs log-ratio of J age totals over final 5 years <= 0.01
max absolute change in J HCV prevalence over final 5 years <= 0.005
```

Also inspect total PWID population stability.

If the default horizon is insufficient, extend `t_end` in a local copy of the `data` object inside the calibration harness. Do not change the default `setup.R` value.

Document the horizon actually used.

---

### Phase 5 — Likelihood

Use a **Binomial likelihood** for prison HCV prevalence.

The age-specific prevalence sample size is the same age-specific prison total
population supplied as the population target:

```r
n_prev <- target_prison_total
```

The observed HCV-positive count must be integer.

First try to recover the original positive counts from the target data source.
If they are unavailable, reconstruct them explicitly:

```r
x_prev <- round(target_prev * n_prev)
```

The supplied prevalence decimals and population sizes are not exactly
integer-compatible, so do not pass `target_prev * n_prev` directly to
`dbinom()`.

Save and report:

```r
target_prev_supplied
n_prev
x_prev
target_prev_binom <- x_prev / n_prev
rounding_difference <- target_prev_binom - target_prev_supplied
```

For model-predicted prevalence `p_hat`:

```text
x_prev[i] ~ Binomial(n_prev[i], p_hat[i])
```

Use:

```r
eps <- 1e-10
p_safe <- pmin(pmax(p_hat, eps), 1 - eps)

nll_prev <- -sum(
  dbinom(x_prev, size = n_prev, prob = p_safe, log = TRUE)
)
```

For prison age-specific total population, retain the working log-scale
observation model unless better uncertainty information is available:

```text
log(N_obs[i]) ~ Normal(log(N_hat[i]), sigma_pop^2)
```

with:

```r
sigma_pop <- 0.10
```

Its objective contribution, up to constants independent of fitted parameters,
is:

```r
nll_pop <- 0.5 * sum(
  ((log(N_hat) - log(target_prison_total)) / sigma_pop)^2
)
```

Minimize:

```r
NLL <- nll_prev + nll_pop
```

Do not arbitrarily normalize or downweight the Binomial likelihood merely
because its age groups have different sample sizes: that sample-size weighting
is part of the specified observation model.

Treat `sigma_pop` as a population calibration/model-discrepancy tolerance, not
an empirical sampling standard error.

In `docs/calibration/likelihood.md`, document:

- the Binomial prevalence likelihood;
- how integer positive counts were obtained;
- any rounding from supplied prevalence to implied Binomial prevalence;
- the population likelihood and `sigma_pop`;
- the exact combined objective implemented in code.

If raw HCV-positive counts later become available, replace reconstructed
rounded counts with the raw integers and rerun calibration.

Do not silently add an informative prior or smoothness penalty to make the fit
look better. If regularization is explored, run it as a separately labelled
sensitivity analysis.

---

### Phase 6 — Multi-start Nelder-Mead

Run deterministic multi-start Nelder-Mead for each retained mortality scenario.

Use at least six starts:

- one baseline start with all log-scales = 0;
- at least five reproducible perturbations around baseline.

Record every start and result.

At minimum save:

```text
initial theta
final theta
contact row scales
beta scales
derived lambda1
objective
convergence code
function evaluations
elapsed time
equilibrium status
prevalence RMSE
maximum prevalence absolute error
population MAPE
maximum population APE
Binomial prevalence deviance
population standardized residual sum of squares
```

After finding the best candidate, run nearby perturbations to check local stability.

Do not select a result only because `optim()` reports convergence.

---

### Phase 7 — Iterative decision loop

Iterate in this order:

1. check implementation correctness;
2. check equilibrium;
3. vary Nelder-Mead initial values;
4. compare evidence-supported mortality scenarios;
5. inspect fitted parameter magnitude and residual pattern;
6. revise the working likelihood only if there is a statistical/data reason;
7. rerun.

Do **not** respond to poor fit by changing protected compartment transitions, ageing, transmission equations, incarceration equations, treatment equations, or the PowerPoint schematic.

Do **not** invent a mortality curve purely to improve objective value.

Append every material decision to:

```text
docs/calibration/DECISIONS.md
```

Each entry must include:

```text
date/time
change considered
reason
evidence
run IDs compared
fit consequence
decision
```

Never delete prior failed attempts from this log.

---

### Phase 8 — Acceptance criteria

A final candidate should satisfy all of the following unless evidence justifies a revised threshold:

```text
equilibrium: pass

prevalence RMSE                    <= 0.02
maximum absolute prevalence error <= 0.03

population MAPE                    <= 0.10
maximum population APE            <= 0.20
```

Also report, but do not use as a conventional residual-df goodness-of-fit test:

```text
Binomial prevalence deviance
population standardized residual sum of squares
total combined NLL
```

Also require:

- similar solutions from multiple starts;
- no numerical instability;
- no extreme scale factor that appears to compensate for another protected
  structural limitation;
- mortality assumptions defensible from evidence.

Because there are 12 fitted scale parameters and 12 target summaries, do not
claim conventional residual degrees of freedom for these calibration
diagnostics.

If these criteria cannot be met defensibly, the correct outcome is a documented calibration failure / structural limitation, not a forced fit.

---

### Phase 9 — Final mortality update

Only after selecting a defensible final mortality specification:

1. update `src/setup.R` only at `mu`, `omega`, and their adjacent citation comments;
2. do not change any other line for model behavior;
3. preserve the source components in `docs/calibration/mortality_review.md` and the selected run output;
4. rerun the selected calibration from a clean R session;
5. verify the result reproduces.

If the selected specification uses mapped age-stratified SMR folded into
`mu`, make the comment explicit that `mu` is now effective model-age-specific
PWID background mortality and that `omega = 1` is intentional. Preserve the
coarser original literature strata and mapping rule in the calibration
documentation.

---

## Required output structure

Use:

```text
output/calibration/<run_id>/
```

Each run should contain, where applicable:

```text
run_config.csv
initial_values.csv
optimization_history.csv
solutions.csv
predictions.csv
residuals.csv
equilibrium.csv
mortality_scenario.csv
sessionInfo.txt
fit.rds
```

Create plots for:

1. observed vs fitted prison prevalence by age;
2. observed vs fitted prison population by age;
3. residuals;
4. contact row scale factors;
5. beta scale factors;
6. objective values across starts;
7. equilibrium stability.

Put final selected report figures into `figures/`.

Create:

```text
docs/calibration/final_report.md
```

The report must state:

- final mortality specification and evidence;
- exact likelihood;
- equilibrium definition;
- best fitted scales;
- predicted vs target values by age group;
- all acceptance metrics;
- multi-start stability;
- sensitivity to mortality scenario;
- known structural limitations;
- whether the model met the target defensibly.

---

## Change-trace requirements

Before work:

```bash
git status --short
git rev-parse HEAD
```

Never destroy unrelated changes.

Record input hashes and Git SHA in run metadata.

Use Git commits as auditable milestones if the environment permits commits. Keep commits logically separated.

At the end show:

```bash
git status --short
git diff -- src/setup.R
git diff -- src/sim.cpp
```

Also verify the PowerPoint hash is unchanged.

The final response must list:

1. files created;
2. files modified;
3. protected-file verification;
4. best calibration metrics;
5. selected mortality scenario;
6. remaining model limitations;
7. exact commands needed to reproduce the final run.

---

## Multi-agent execution, if available

Use multiple agents only for safe parallel work.

### Subagent A — mortality researcher

Task:

- research Singapore baseline mortality and PWID SMR;
- build evidence table;
- propose mortality scenarios;
- write only research/scenario files;
- do not edit `setup.R`, `sim.cpp`, or the PowerPoint.

### Subagent B — model auditor / likelihood reviewer

Task:

- inspect model implementation and schematic;
- verify target extraction;
- review equilibrium and likelihood design;
- identify structural calibration limitations;
- write only audit/likelihood documentation;
- do not modify protected model files.

Run A and B in parallel if possible.

### Lead/integrator

After A and B finish:

- reconcile their findings;
- implement the calibration harness;
- run all Nelder-Mead calibration;
- own all output generation;
- be the only agent allowed to make the final permitted mortality edit to `setup.R`.

Do not run independent optimization agents that each rewrite shared source files.

If subagents are unavailable, perform the same phases sequentially yourself.

---

## First actions

Start now by:

1. reading `AGENTS.md`;
2. recording Git state and hashes;
3. reading all three protected/baseline model files;
4. producing the model audit;
5. researching mortality evidence;
6. only then implementing and running calibration.

Do not ask for permission to modify protected model design. It is not permitted in this task.
