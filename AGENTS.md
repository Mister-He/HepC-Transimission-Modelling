# AGENTS.md

## Scope

This repository implements an age-structured HCV compartmental model for PWID in Singapore.

Current calibration goal:

- fit equilibrium prison (`J`) age-specific HCV prevalence;
- fit equilibrium prison (`J`) age-specific population;
- fit 12 positive parameters:
  - 6 contact-matrix row scales;
  - 6 age-specific `beta`/inflow scales;
- use `optim(..., method = "Nelder-Mead")`;
- Nelder-Mead is the endpoint: do **not** add an HMC second stage.

Authoritative model files:

```text
src/setup.R
src/sim.cpp
Model schematic.pptx
```

Read all three before implementing calibration.

---

## Edit boundaries

### Never modify

```text
src/sim.cpp
Model schematic.pptx
```

Do not reformat, rename, regenerate, or “fix” them.

### `src/setup.R`

Read the full file.

Only the final selected baseline-mortality specification may be changed:

```text
params$mu
params$omega
```

Adjacent mortality comments/citations may also be updated.

Do not change in this task:

```text
mu_DC
mu_HCC
psi_DC
psi_HCC
```

Do not change any non-mortality quantity, including `C_contact`, `beta`,
`c_composite`, incarceration parameters, treatment, progression, initial
conditions, age structure, or simulation defaults.

Keep exploratory mortality scenarios in new files; do not repeatedly overwrite
`setup.R`.

New work belongs under:

```text
src/calibration/
docs/calibration/
output/calibration/
```

Use `figures/` only for stable final figures.

---

## Model facts to preserve

The C++ implementation has:

```text
4 strata: D, J, F, X
4 liver stages: NC, CC, DC, HCC
4 HCV states: u, a, c, t
6 age groups
384 total compartments
```

`J` is currently incarcerated.

The current implementation has no new HCV force of infection inside `J`.

Background PWID mortality enters as:

```text
mu[i] * omega
```

and `omega` is scalar.

The current repository baseline has already corrected the ageing transition to:

```cpp
y_change = y[base + i] / 10.0 * dt;
```

This is the expected protected implementation for the 10-year adult age bands.
Verify that the working tree contains this `/10.0` implementation during
preflight. Treat it as part of the current immutable baseline:

```text
do not revert it to /5.0
do not further modify ageing in this task
```

If the checked-out repository unexpectedly still contains `/5.0`, stop and
report the discrepancy rather than editing `src/sim.cpp`.

Do not force a fit by changing protected model structure.

---

## Calibration parameterization

Use unconstrained log-parameters:

```r
contact_scale <- exp(theta_contact)
beta_scale    <- exp(theta_beta)
```

Apply contact scaling row-wise:

```r
pm$C_contact[i, ] <-
  base_params$C_contact[i, ] * contact_scale[i]
```

Apply inflow scaling:

```r
pm$beta <- base_params$beta * beta_scale
```

Preserve the existing project coupling:

```r
c_true <- base_params$c_composite / beta_scale
pm$lambda1 <- base_params$lambda3 * c_true
```

Do not silently decouple `beta_scale`, `c_composite`, and `lambda1`.

Never write fitted `C_contact` or `beta` values back to `setup.R`.

---

## Target extraction

For age group `i`:

```text
N_hat[i] =
  sum of all J states across liver stage and HCV state
```

For current/viremic HCV prevalence:

```text
I_hat[i] =
  J acute + J chronic + J treatment,
  summed across liver stages

p_hat[i] = I_hat[i] / N_hat[i]
```

Verify the epidemiological definition of the observed prevalence target when
its source is available.

State `u` combines susceptible and post-SVR individuals. Therefore the current
state space cannot exactly reconstruct anti-HCV antibody prevalence. If the
target is seroprevalence rather than current/viremic infection, report the
mismatch instead of hiding it.

---

## Equilibrium and likelihood rules

Do not call a transient endpoint equilibrium.

At minimum compare target summaries at the final time and about five years
earlier. Initial working criteria:

```text
max_i |log(N_hat_i(T) / N_hat_i(T-5))| <= 0.01
max_i |p_hat_i(T) - p_hat_i(T-5)|      <= 0.005
```

Also inspect total PWID population stability.

If the default horizon is insufficient, extend it only in a local calibration
copy of the simulation data.

Use a **Binomial likelihood** for the age-specific prison HCV prevalence,
because the prison age-specific total population is also the prevalence sample
size.

Let:

```text
n_prev[i] = observed prison total population in age group i
x_prev[i] = observed HCV-positive count in age group i
p_hat[i]  = model-predicted prison HCV prevalence
```

Then:

```text
x_prev[i] ~ Binomial(n_prev[i], p_hat[i])
```

and the prevalence contribution to the objective is:

```r
nll_prev <- -sum(
  dbinom(x_prev, size = n_prev, prob = p_hat, log = TRUE)
)
```

The supplied prevalence decimals multiplied by the supplied population sizes
are not exactly integers. Therefore:

1. use original integer HCV-positive counts if they can be recovered from the
   target source;
2. otherwise reconstruct and explicitly record:

```r
x_prev <- round(target_prev * target_prison_total)
```

3. save both the supplied prevalence and the implied Binomial prevalence
   `x_prev / n_prev`;
4. do not pass non-integer successes to `dbinom()`.

For the prison population target, retain the working log-scale observation
model unless better target uncertainty is available:

```text
log(N_obs[i]) ~ Normal(log(N_hat[i]), sigma_pop^2)
sigma_pop = 0.10
```

The population `sigma_pop` is a calibration/model-discrepancy tolerance, not a
claimed sampling standard error.

The combined objective is:

```text
NLL = NLL_binomial_prevalence + NLL_population
```

Do not reweight the Binomial block merely to make optimization easier unless a
separate, explicitly documented sensitivity analysis is justified.

Do not silently add informative priors or regularization just to improve fit.

Initial target-space acceptance criteria:

```text
prevalence RMSE                    <= 0.02
maximum absolute prevalence error <= 0.03
population MAPE                    <= 0.10
maximum population APE            <= 0.20
```

Also report:

```text
Binomial deviance for prevalence
population standardized residual sum of squares
total NLL
```

Because there are 12 fitted scale parameters and 12 target summaries, do not
present any of these as a conventional goodness-of-fit test with residual
degrees of freedom.

Require stability across multiple Nelder-Mead starts.

---

## Mortality rules

Prefer authoritative Singapore mortality data, especially SingStat
age-specific death rates or complete life tables.

Use a mortality period aligned with the calibration target year when possible.

For PWID excess mortality, prioritize:

1. credible age-stratified PWID SMRs, even when the paper reports fewer than
   the model's six age groups;
2. Asian PWID cohorts where available;
3. well-described age-stratified PWID cohorts elsewhere;
4. systematic reviews/meta-analyses for context.

Do **not** require a literature source to report all six model age groups.
Two-, three-, four-, or five-band age-specific SMR evidence is acceptable when
the study is otherwise credible.

Map coarser literature age bands to the six model groups transparently. The
default mapping is piecewise constant: model age groups fully contained within
one reported SMR band inherit that band's SMR. Do not invent smooth
interpolation between reported SMRs.

If literature and model boundaries overlap imperfectly, document the mapping
assumption. Use weighted aggregation only when defensible weights are
available; otherwise retain a coarser piecewise scenario and flag the
uncertainty.

For every mortality source record:

```text
citation / DOI / PMID / URL
country/setting and study years
age bands
sex distribution
cohort size/person-years when available
HIV/HCV context when available
mortality measure and CI
reference population/standardization
reason for inclusion/exclusion
```

Do not invent missing age-specific SMRs.

Do not convert crude mortality rates into SMRs without the required reference
rates.

Because C++ accepts only scalar `omega`, do not change it to a vector.

If evidence supports age-stratified SMRs—whether six bands or fewer—first map
the reported strata to the six model age groups as documented above. Then
represent the resulting model-age-specific `SMR_i` through the existing
mortality product:

```text
effective_mu[i] = singapore_mu[i] * SMR_i
omega = 1
```

Keep `singapore_mu`, `SMR_i`, and `effective_mu` separately in calibration
files. If selected for the final baseline, explicitly document that `mu` now
represents effective age-specific PWID background mortality.

---

## Optimization rules

Use deterministic multi-start Nelder-Mead.

Minimum:

```text
1 start with all log-scales = 0
at least 5 reproducible perturbed starts
```

For every start record:

```text
initial theta
final theta
objective
convergence code
function evaluations
elapsed time
equilibrium status
fit metrics
```

Re-run the best solution from nearby perturbations to test local stability.

Do not run HMC afterward.

---

## Traceability

Before changes:

```bash
git status --short
git rev-parse HEAD
```

Do not overwrite unrelated user work.

Do not use destructive cleanup commands such as:

```bash
git reset --hard
git clean -fd
git checkout -- .
```

unless explicitly requested.

Every calibration run must record:

```text
run ID/timestamp
Git SHA
hashes of setup.R and sim.cpp
targets
mortality scenario
equilibrium settings
likelihood settings
optimizer settings
initial/fitted values
predictions
fit metrics
R/package versions
```

Maintain an append-only:

```text
docs/calibration/DECISIONS.md
```

Record every material change to mortality assumptions, likelihood, equilibrium
criteria, parameterization, optimizer settings, or scenario acceptance.

Use logical Git commits where appropriate. Do not erase failed runs from the
audit trail.

---

## Multi-agent policy

Parallelize only read-only/safe work.

Recommended roles:

```text
Mortality researcher
  -> mortality evidence and scenario proposals

Model auditor / likelihood reviewer
  -> read-only model audit and methodology review

Calibration integrator
  -> sole owner of calibration scripts, optimization,
     outputs, and final permitted setup.R mortality edit
```

The first two may run in parallel.

Do not let multiple agents edit `setup.R`.

If subagents are unavailable, perform the same roles sequentially.

---

## Final validation

Before completion:

1. prove `src/sim.cpp` is unchanged;
2. prove `Model schematic.pptx` is unchanged;
3. prove any `src/setup.R` diff is limited to permitted mortality lines/comments;
4. reproduce the selected fit from a clean R session;
5. confirm equilibrium;
6. confirm target-space metrics;
7. compare multiple starts;
8. document structural limitations.

A visually good plot is not sufficient evidence of calibration success.

Stop and report rather than forcing a fit if:

- protected changes would be required;
- the target prevalence definition is incompatible with model states;
- mortality changes lack defensible evidence;
- fit remains poor across reasonable multi-start runs;
- acceptable fit requires extreme or unstable parameters;
- the checked-out repository does not contain the expected `/10.0` ageing
  implementation;
- another protected structural limitation appears to dominate the age pattern.

A valid result is: **the current protected model cannot defensibly meet the
targets under the allowed calibration scope.**
