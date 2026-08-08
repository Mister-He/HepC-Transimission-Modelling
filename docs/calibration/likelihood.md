# Statistical model specification (revision v1)

## Parameters (12 fitted)

```r
theta[1:6]  contact-matrix row scales   contact_scale = exp(theta)
theta[7:12] beta inflow scales          beta_scale    = exp(theta)
pm$C_contact[i, ] <- base$C_contact[i, ] * contact_scale[i]
pm$beta           <- base$beta          * beta_scale
```

Transmission is constant (`m_min = m_max = 1`); `eta_s = rep(1, 6)`.
The excess-mortality contingency (`--excess-mortality`) adds
`theta[13]` with `eta_s[6] = 1 + 4*plogis(theta13)` (bounded [1,5]) and is
only used if the 60+ group cannot be fitted with the 12 parameters (see
AGENTS.md and the dedicated rationale file if activated).

## Observation model

```text
x_prev[i] ~ Binomial(n_prev[i], p_hat[i])     n_prev = prison_total
x_prev = round(target_prev * prison_total)    (11, 215, 394, 792, 785, 145)
p_hat  = seroprevalence (a+c+t+s)/N           (primary; anti-HCV serology)
log(N_obs[i]) ~ Normal(log(N_hat[i]), 0.10^2) (sigma_pop = 0.10)
NLL = nll_prev + nll_pop
```

No reweighting; no informative priors. The equilibrium gate
(T vs T-5, criteria 0.01 / 0.005) is enforced as a 1e6 x normalized
violation penalty during optimisation; the Laplace Hessian uses the pure
statistical NLL.

## Laplace approximation

Finite-difference Hessian (numDeriv) of the pure NLL at the optimum;
generalized inverse with relative eigenvalue cutoff 1e-4 (condition number
<= 1e4); 1000 Monte Carlo draws from N(theta_hat, Sigma) re-simulated to
the 12 target summaries; equilibrium-infeasible draws discarded; 95%
intervals per age group vs observed intervals (Jeffreys Binomial for
prevalence; log-Normal sigma_pop = 0.10 for population). Overlap is
reported per age group.

## Acceptance criteria

```text
prevalence RMSE <= 0.02; max |prev error| <= 0.03
population MAPE <= 0.10; max APE <= 0.20; equilibrium pass
```

Also reported: Binomial deviance, population SRSS, total NLL, multi-start
stability, Laplace overlap.
