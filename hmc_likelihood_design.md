# HMC likelihood design for fitting an age-stratified compartmental model

## Purpose

You are modifying a local HMC-based compartmental model fitting codebase. The goal is to implement a joint likelihood that simultaneously fits:

1. age-group-specific total population size, using external population observations;
2. age-group-specific prevalence, using prevalence observations;
3. a weak-to-moderate shape regularization prior on age-group population proportions.

The important modeling principle is:

> Age-group population counts should not be treated as independent Gaussian targets. They are a compositional vector constrained by total population size. The likelihood should fit the relative age composition, while preserving the total population constraint.

The final design should combine:

```text
log posterior
= log likelihood for observed age composition
+ log likelihood for observed prevalence
+ log shape prior on age-group proportions
+ other model parameter priors
```

---

## Core model quantities

Assume there are:

```text
A = number of age groups
T = number of observation time points
```

For each time point `t` and age group `a`, the compartmental model produces:

```text
N_model[t, a] = model-predicted total population in age group a
I_model[t, a] = model-predicted number of prevalent cases, infected people, or target disease-state people in age group a
```

The model-implied age-group prevalence is:

```text
p_model[t, a] = I_model[t, a] / N_model[t, a]
```

The model-implied age composition is:

```text
pi_model[t, a] = N_model[t, a] / sum_a N_model[t, a]
```

If the model requires constant total population, enforce it deterministically:

```text
sum_a N_model[t, a] = N_total
```

or equivalently:

```text
N_model[t, a] = N_total * pi_model[t, a]
sum_a pi_model[t, a] = 1
```

Do not add a separate independent likelihood for total population unless the scientific model explicitly allows total population measurement error.

---

## Observed data inputs

The fitting code should accept the following observed data.

### Population observations

External reliable population counts:

```text
Y_pop_obs[t, a]
```

These are external source population counts for each age group.

Convert them into observed age composition:

```text
q_obs[t, a] = Y_pop_obs[t, a] / sum_a Y_pop_obs[t, a]
```

The likelihood should compare `q_obs[t, :]` with `pi_model[t, :]`, not compare `Y_pop_obs[t, a]` independently against `N_model[t, a]`.

### Prevalence observations

There are two supported cases.

#### Case A: prevalence estimate with standard error

Inputs:

```text
prev_obs[t, a]       = observed age-specific prevalence estimate
prev_se[t, a]        = standard error of observed prevalence, preferably on probability scale unless already transformed
```

Use a logit-normal likelihood.

#### Case B: numerator and denominator

Inputs:

```text
cases_obs[t, a]      = observed prevalent case count, positives, or disease-state count
sample_size[t, a]    = denominator used to estimate prevalence
```

Use a binomial or beta-binomial likelihood. Prefer beta-binomial if observations show overdispersion or survey/design heterogeneity.

---

# Recommended joint likelihood

## 1. Population composition likelihood

Use an additive log-ratio normal likelihood.

Choose one reference age group, usually the final age group:

```text
ref = A
```

For each `t` and for each `a = 1, ..., A - 1`, define:

```text
z_pop_obs[t, a]   = log(q_obs[t, a] / q_obs[t, ref])
z_pop_model[t, a] = log(pi_model[t, a] / pi_model[t, ref])
```

Then use:

```text
z_pop_obs[t, a] ~ Normal(z_pop_model[t, a], sigma_pop[a])
```

Equivalently, the log-likelihood contribution is:

```text
log L_pop = sum_t sum_{a=1}^{A-1} normal_lpdf(
    z_pop_obs[t, a] | z_pop_model[t, a], sigma_pop[a]
)
```

### Interpretation

`population observations` are reliable external sources, so `sigma_pop[a]` should usually be small, but not zero.

Suggested starting values or priors:

```text
sigma_pop[a] fixed between 0.01 and 0.05
```

or:

```text
sigma_pop[a] ~ HalfNormal(0.05)
```

If external source uncertainty is known by age group, use that information directly.

### Why not independent population-count likelihood?

Avoid this:

```text
Y_pop_obs[t, a] ~ Normal(N_model[t, a], sigma_pop_count[a])
```

because it treats age-group counts as independent. With fixed total population, age-group populations are negatively correlated components of a composition. The correct target is the relative age composition.

---

## 2. Shape regularization prior on age-group proportions

Add a second-difference shape prior on centered log proportions.

For each time point `t`, define:

```text
log_pi[t, a] = log(pi_model[t, a])
mean_log_pi[t] = mean_a log_pi[t, a]
ell[t, a] = log_pi[t, a] - mean_log_pi[t]
```

Then for each internal age group `a = 2, ..., A - 1`, compute the second difference:

```text
d2_ell[t, a] = ell[t, a + 1] - 2 * ell[t, a] + ell[t, a - 1]
```

Use a Student-t prior:

```text
d2_ell[t, a] ~ StudentT(nu_shape, 0, sigma_shape)
```

Recommended default:

```text
nu_shape = 4
sigma_shape ~ HalfNormal(0.2 to 0.5)
```

The log-prior contribution is:

```text
log p_shape = sum_t sum_{a=2}^{A-1} student_t_lpdf(
    d2_ell[t, a] | nu_shape, 0, sigma_shape
)
```

### Interpretation

This prior does not force adjacent age-group proportions to be equal. It regularizes the curvature of the log-age-composition profile.

It allows:

```text
- increasing profiles
- decreasing profiles
- hump-shaped profiles
- cohort effects
- moderate non-smoothness
```

but discourages unsupported local zig-zags.

### Why Student-t instead of Normal?

Use Student-t because age composition may have real local irregularities due to cohort structure, migration, demographic shocks, or age-bin definitions. Student-t regularizes most curvature terms but allows occasional large deviations.

Start with:

```text
nu_shape = 4
```

Avoid using extremely small `nu_shape` initially because it can make posterior geometry harder for HMC.

---

## 3. Prevalence likelihood

The prevalence likelihood depends on the format of the observed prevalence data.

---

### Option 1: logit-normal likelihood for observed prevalence estimates

Use this when observations are prevalence estimates with uncertainty.

For each `t, a`:

```text
p_model[t, a] = I_model[t, a] / N_model[t, a]
```

Apply numerical clipping before logit transforms:

```text
eps = 1e-9
p_model_clipped = min(max(p_model[t, a], eps), 1 - eps)
prev_obs_clipped = min(max(prev_obs[t, a], eps), 1 - eps)
```

Then:

```text
logit(prev_obs[t, a]) ~ Normal(
    logit(p_model[t, a]),
    sqrt(prev_logit_se[t, a]^2 + tau_prev^2)
)
```

If `prev_se[t, a]` is given on the probability scale, approximate the logit-scale standard error using the delta method:

```text
prev_logit_se[t, a] = prev_se[t, a] / (prev_obs[t, a] * (1 - prev_obs[t, a]))
```

Then add an extra discrepancy parameter:

```text
tau_prev ~ HalfNormal(0.2 to 1.0)
```

The log-likelihood contribution is:

```text
log L_prev = sum_t sum_a normal_lpdf(
    logit(prev_obs_clipped[t, a]) |
    logit(p_model_clipped[t, a]),
    sqrt(prev_logit_se[t, a]^2 + tau_prev^2)
)
```

---

### Option 2: beta-binomial likelihood for case counts

Use this when the data include numerator and denominator.

For each `t, a`:

```text
cases_obs[t, a] ~ BetaBinomial(sample_size[t, a], alpha[t, a], beta[t, a])
```

Parameterize alpha and beta using the model prevalence and an overdispersion/concentration parameter:

```text
alpha[t, a] = p_model[t, a] * phi_prev
beta[t, a]  = (1 - p_model[t, a]) * phi_prev
```

Recommended prior:

```text
phi_prev > 0
log(phi_prev) ~ Normal(log(100), 1)
```

Larger `phi_prev` means closer to binomial. Smaller `phi_prev` means more overdispersion.

If overdispersion is not needed, use:

```text
cases_obs[t, a] ~ Binomial(sample_size[t, a], p_model[t, a])
```

But beta-binomial is generally safer for real prevalence observations.

---

# Full target density

The full log posterior target should be:

```text
log_target =
    log_L_pop
  + log_L_prev
  + log_p_shape
  + log_priors_for_epidemiological_parameters
  + log_priors_for_observation_parameters
```

Where:

```text
log_L_pop    = age-composition observation likelihood
log_L_prev   = prevalence observation likelihood
log_p_shape  = second-difference shape prior on centered log proportions
```

---

# Implementation steps for the coding agent

## Step 1: Locate the current likelihood code

Find the code that currently computes the log likelihood or target density for HMC. Search for function names or blocks like:

```text
log_likelihood
likelihood
target +=
model {
observe
loss
calibration target
population target
prevalence target
```

Identify the current treatment of:

```text
- age-specific population counts
- age-specific prevalence
- total population constraint
```

If the code currently uses independent Normal or Poisson likelihoods for age-group population counts, replace that component with the composition likelihood described here.

---

## Step 2: Ensure model outputs are available

The likelihood function must have access to:

```text
N_model[t, a]
I_model[t, a]
```

Check dimensions carefully:

```text
N_model: T x A
I_model: T x A
```

If the model state has multiple compartments, compute:

```text
N_model[t, a] = sum over all living compartments in age group a
I_model[t, a] = sum over compartments that define prevalence in age group a
```

Do not accidentally include dead, removed, or cumulative-history compartments in `N_model` unless the scientific model defines them as part of the current population.

---

## Step 3: Enforce or verify total population constraint

If total population must be constant, enforce one of the following:

### Preferred deterministic construction

```text
N_model[t, a] = N_total * pi_model[t, a]
sum_a pi_model[t, a] = 1
```

### Or verification after model simulation

```text
N_total_model[t] = sum_a N_model[t, a]
```

Then check that:

```text
N_total_model[t] approximately equals N_total
```

Avoid adding a separate independent total-population likelihood unless total population is genuinely observed with measurement error.

---

## Step 4: Build observed age-composition data

From the input population observations:

```text
Y_pop_obs[t, a]
```

compute:

```text
q_obs[t, a] = Y_pop_obs[t, a] / sum_a Y_pop_obs[t, a]
```

Precompute this outside the HMC inner loop if possible.

Add a small epsilon only if zeros are possible:

```text
eps = 1e-12
q_obs_safe[t, a] = max(q_obs[t, a], eps)
q_obs_safe[t, :] = q_obs_safe[t, :] / sum_a q_obs_safe[t, a]
```

Use the safe normalized values in log-ratio transforms.

---

## Step 5: Implement population composition likelihood

For each time point `t`:

```text
pi_model[t, a] = N_model[t, a] / sum_a N_model[t, a]
```

Apply epsilon stabilization:

```text
pi_safe[t, a] = max(pi_model[t, a], eps)
pi_safe[t, :] = pi_safe[t, :] / sum_a pi_safe[t, a]
```

Let:

```text
ref = A
```

Then for `a = 1, ..., A - 1`:

```text
z_obs   = log(q_obs_safe[t, a] / q_obs_safe[t, ref])
z_model = log(pi_safe[t, a] / pi_safe[t, ref])
log_L_pop += normal_lpdf(z_obs | z_model, sigma_pop[a])
```

Use age-specific `sigma_pop[a]` if available. Otherwise use a shared `sigma_pop`.

---

## Step 6: Implement second-difference shape prior

For each time point `t`:

```text
log_pi[a] = log(pi_safe[t, a])
ell[a] = log_pi[a] - mean(log_pi)
```

For each internal age group `a = 2, ..., A - 1`:

```text
d2 = ell[a + 1] - 2 * ell[a] + ell[a - 1]
log_p_shape += student_t_lpdf(d2 | nu_shape, 0, sigma_shape)
```

Default hyperparameters:

```text
nu_shape = 4
sigma_shape ~ HalfNormal(0.2 to 0.5)
```

If the framework does not support Student-t easily, use Normal as a first-pass fallback:

```text
d2 ~ Normal(0, sigma_shape)
```

but prefer Student-t for robustness.

---

## Step 7: Implement prevalence likelihood

Compute:

```text
p_model[t, a] = I_model[t, a] / N_model[t, a]
```

Apply clipping:

```text
p_model_safe = min(max(p_model[t, a], eps), 1 - eps)
```

### If using prevalence estimates

Precompute or compute:

```text
prev_obs_safe = min(max(prev_obs[t, a], eps), 1 - eps)
prev_logit_se[t, a] = prev_se[t, a] / (prev_obs_safe * (1 - prev_obs_safe))
```

Then:

```text
sigma_prev_total[t, a] = sqrt(prev_logit_se[t, a]^2 + tau_prev^2)
log_L_prev += normal_lpdf(
    logit(prev_obs_safe) |
    logit(p_model_safe),
    sigma_prev_total[t, a]
)
```

### If using case counts

Use beta-binomial:

```text
alpha = p_model_safe * phi_prev
beta  = (1 - p_model_safe) * phi_prev
log_L_prev += beta_binomial_lpmf(cases_obs[t, a] | sample_size[t, a], alpha, beta)
```

If beta-binomial is not available, use binomial as a temporary fallback:

```text
log_L_prev += binomial_lpmf(cases_obs[t, a] | sample_size[t, a], p_model_safe)
```

---

## Step 8: Add priors for observation parameters

Recommended priors:

```text
sigma_pop[a] ~ HalfNormal(0.05)
```

or fixed values if external source uncertainty is known.

```text
sigma_shape ~ HalfNormal(0.2 to 0.5)
nu_shape = 4 fixed
```

For logit-normal prevalence likelihood:

```text
tau_prev ~ HalfNormal(0.5)
```

For beta-binomial prevalence likelihood:

```text
log(phi_prev) ~ Normal(log(100), 1)
```

Do not choose very tight priors before checking posterior predictive fit.

---

## Step 9: Numerical stability requirements

Implement these safeguards:

```text
1. Never take log(0).
2. Never take logit(0) or logit(1).
3. Use eps = 1e-9 or 1e-12 depending on framework precision.
4. Check that all N_model[t, a] > 0 before computing prevalence and proportions.
5. Check that I_model[t, a] >= 0 and I_model[t, a] <= N_model[t, a].
6. Ensure q_obs[t, :] and pi_model[t, :] each sum to 1 after stabilization.
7. Avoid using real total population size as an effective multinomial sample size unless this is intentional.
```

---

## Step 10: Diagnostics and acceptance criteria

After implementation, run short HMC tests and posterior predictive checks.

The implementation is acceptable only if:

```text
1. HMC has no or very few divergent transitions.
2. Effective sample sizes are reasonable for key parameters.
3. R-hat is close to 1 for key parameters.
4. Posterior predictive age composition matches external population composition within intended uncertainty.
5. Posterior predictive prevalence matches prevalence observations without population likelihood completely dominating.
6. Total population is constant if the model requires it.
7. The shape prior reduces unsupported local zig-zags but does not flatten real demographic structure.
```

Compare three model variants:

```text
Variant A: no shape prior
Variant B: Normal second-difference shape prior
Variant C: Student-t second-difference shape prior
```

Prefer Variant C if it improves stability and posterior predictive fit without oversmoothing.

---

# Suggested pseudocode

```pseudo
function compute_log_target(params, data):
    states = simulate_compartmental_model(params)

    N_model = extract_age_total_population(states)
    I_model = extract_age_prevalent_cases(states)

    log_L_pop = 0
    log_L_prev = 0
    log_p_shape = 0
    log_prior = compute_parameter_priors(params)

    for t in 1:T:
        N_total_t = sum_a N_model[t, a]
        pi = N_model[t, :] / N_total_t
        pi = stabilize_simplex(pi, eps)

        q = data.q_obs[t, :]
        q = stabilize_simplex(q, eps)

        ref = A

        for a in 1:(A - 1):
            z_obs = log(q[a] / q[ref])
            z_model = log(pi[a] / pi[ref])
            log_L_pop += normal_lpdf(z_obs, z_model, sigma_pop[a])

        log_pi = log(pi)
        ell = log_pi - mean(log_pi)

        for a in 2:(A - 1):
            d2 = ell[a + 1] - 2 * ell[a] + ell[a - 1]
            log_p_shape += student_t_lpdf(d2, nu_shape, 0, sigma_shape)

        for a in 1:A:
            p = I_model[t, a] / N_model[t, a]
            p = clamp(p, eps, 1 - eps)

            if data.prevalence_format == "estimate_with_se":
                y = clamp(data.prev_obs[t, a], eps, 1 - eps)
                se_logit = data.prev_se[t, a] / (y * (1 - y))
                sigma_total = sqrt(se_logit^2 + tau_prev^2)
                log_L_prev += normal_lpdf(logit(y), logit(p), sigma_total)

            if data.prevalence_format == "case_count":
                alpha = p * phi_prev
                beta = (1 - p) * phi_prev
                log_L_prev += beta_binomial_lpmf(
                    data.cases_obs[t, a],
                    data.sample_size[t, a],
                    alpha,
                    beta
                )

    return log_L_pop + log_L_prev + log_p_shape + log_prior
```

---

# Minimal Stan-style sketch

```stan
for (t in 1:T) {
    vector[A] pi_t;
    vector[A] log_pi_t;
    vector[A] ell_t;

    pi_t = to_vector(N_model[t]) / sum(to_vector(N_model[t]));

    for (a in 1:A) {
        pi_t[a] = fmax(pi_t[a], eps);
    }
    pi_t = pi_t / sum(pi_t);

    // Population composition likelihood
    for (a in 1:(A - 1)) {
        target += normal_lpdf(
            log(q_obs[t, a] / q_obs[t, A]) |
            log(pi_t[a] / pi_t[A]),
            sigma_pop[a]
        );
    }

    // Shape prior on centered log proportions
    for (a in 1:A) {
        log_pi_t[a] = log(pi_t[a]);
    }
    ell_t = log_pi_t - mean(log_pi_t);

    for (a in 2:(A - 1)) {
        target += student_t_lpdf(
            ell_t[a + 1] - 2 * ell_t[a] + ell_t[a - 1] |
            nu_shape,
            0,
            sigma_shape
        );
    }

    // Prevalence likelihood: logit-normal version
    for (a in 1:A) {
        real p_model;
        real y_obs;
        real se_logit;
        real sigma_total;

        p_model = I_model[t, a] / N_model[t, a];
        p_model = fmin(fmax(p_model, eps), 1 - eps);

        y_obs = fmin(fmax(prev_obs[t, a], eps), 1 - eps);
        se_logit = prev_se[t, a] / (y_obs * (1 - y_obs));
        sigma_total = sqrt(square(se_logit) + square(tau_prev));

        target += normal_lpdf(
            logit(y_obs) |
            logit(p_model),
            sigma_total
        );
    }
}
```

---

# Notes for refactoring existing code

When modifying the local codebase, make the likelihood components modular:

```text
compute_population_composition_loglik(...)
compute_age_shape_logprior(...)
compute_prevalence_loglik(...)
compute_total_logposterior(...)
```

This makes it easier to test each component independently.

Add unit tests for:

```text
1. q_obs conversion from population counts to proportions.
2. log-ratio likelihood dimensions: only A - 1 ratios per time point.
3. second-difference prior dimensions: only A - 2 curvature terms per time point.
4. prevalence calculation: I_model / N_model.
5. numerical clipping for zero or one prevalence values.
6. total population preservation.
```

---

# Final modeling recommendation

Use this final likelihood design unless there is strong reason to deviate:

```text
Population target:
    additive log-ratio Normal likelihood on age composition

Shape regularization:
    Student-t second-difference prior on centered log age proportions

Prevalence target:
    logit-Normal likelihood if prevalence estimates and SEs are available
    beta-binomial likelihood if numerator and denominator are available

Total population:
    deterministic constraint, not independent likelihood
```

This design respects the fact that age-group populations are compositional, preserves constant total population, uses reliable external population data without treating every age group as an independent count, and regularizes the age-profile shape without forcing an unrealistically smooth demographic curve.
