# Calibration final report — HCV PWID model (3-strata), Singapore 2015

Status: COMPLETE. Final run: `run2_v1_12p_warm`
(`output/calibration/run2_v1_12p_warm/`), Git SHA `ad931f0` (working
tree), revision: highest-defensible progression rates, constant
transmission, 12 fitted parameters.

## 1. Final specification (this revision)

- **Mortality (2015 baseline)**: SingStat 2015 age-specific death rates
  (mu = 0.0002, 0.00025, 0.00045, 0.0012, 0.00335, 0.0218818 /yr for
  <20, 20-29, 30-39, 40-49, 50-59, 60+) x scalar PWID SMR
  omega = 14.68 (Mathers et al. 2013). `mu_DC = 0.130`, `mu_HCC = 0.430`
  retained.
- **Progression rates** (highest defensible values):
  `p_NC_CC = 0.027` (Thein 2008); `p_CC_DC = 0.0788`,
  `p_CC_HCC = 0.0479` (Alazawi et al. 2010, PMID 20497143,
  untreated-only); `p_DC_HCC = 0.0464` (Rivera-Irigoin et al. 2006,
  PMID 17209765, upper 95% CI). GT3 RRs (Kanwal 2014) retained. See
  [natural_history_review.md](natural_history_review.md).
- **Transmission**: constant (`m_min = m_max = 1`), no fitted multiplier;
  the historical level is absorbed by the fitted contact row scales.
- **Fitted parameters (12)**: 6 contact-matrix row scales, 6 beta inflow
  scales (`exp(theta)`). No excess-mortality parameter (not needed).

## 2. Likelihood and equilibrium

```text
x_prev[i] ~ Binomial(n_prev[i], p_hat[i])   seroprevalence (primary)
log(N_obs[i]) ~ Normal(log(N_hat[i]), 0.10^2)
NLL = nll_prev + nll_pop
```

Horizon t = -10..55 (calendar 1960-2025); targets at t = 45 (2015);
equilibrium gate at t = 55 vs 50 (2025 vs 2020). Integer positives
x_prev = (11, 215, 394, 792, 785, 145).

## 3. Results (run2_v1_12p_warm, best start `warm_12p_22p5`)

### 3.1 Fitted scales (12 parameters)

| Parameter | Value | | Parameter | Value |
|---|---:|---|---:|
| contact scale 1 | 5.51 | | beta scale 1 | 0.169 |
| contact scale 2 | 0.082 | | beta scale 2 | 0.920 |
| contact scale 3 | 0.053 | | beta scale 3 | 1.263 |
| contact scale 4 | 0.465 | | beta scale 4 | 5.544 |
| contact scale 5 | 7.449 | | beta scale 5 | 14.849 |
| contact scale 6 | 0.150 | | beta scale 6 | 109.595 |

beta_scale > 1 for age groups 3-6; groups 1-2 (<20, 20-29) remain below 1
(documented soft-rule outcome, DECISIONS.md). Modelled annual inflows:
~43, ~593, ~369, ~521, ~312, ~110/yr (total ~1948/yr, exceeding the
CNB-anchored base total 1309/yr).

### 3.2 Predicted vs target (2015)

| Age | N_obs | N_hat | p_obs (Binom) | p_hat (sero) |
|---|---:|---:|---:|---:|
| <20 | 99 | 98.8 | 0.1111 | 0.1060 |
| 20-29 | 1244 | 1255.5 | 0.1728 | 0.1731 |
| 30-39 | 1467 | 1463.6 | 0.2686 | 0.2684 |
| 40-49 | 1841 | 1880.4 | 0.4302 | 0.4303 |
| 50-59 | 1628 | 1422.4 | 0.4822 | 0.4780 |
| 60+ | 409 | 459.5 | 0.3545 | 0.3740 |

### 3.3 Acceptance metrics

| Metric | Criterion | Value | Pass |
|---|---:|---:|---|
| Prevalence RMSE | <= 0.02 | 0.00841 | PASS |
| Max absolute prevalence error | <= 0.03 | 0.01953 | PASS |
| Population MAPE | <= 0.10 | 0.04743 | PASS |
| Max population APE | <= 0.20 | 0.12631 | PASS |
| Equilibrium (T vs T-5) | pass | max_prev 0.00125 | PASS |
| Binomial deviance | report | 0.811 | - |
| Population SRSS | report | 3.231 | - |
| Total NLL | report | 22.44 (prev 20.82 + pop 1.62) | - |

### 3.4 Multi-start stability

11 deterministic starts (1 zero + 3 informed + 2 warm + 5 perturbed).
The two warm starts (run7-derived and its converged point) reproduce the
same basin: NLL 22.48 / 22.44. Other starts land in higher basins
(31.3-138.7). Local perturbations (sd 0.15) are dominated by the
equilibrium-gate penalty (objectives 34-162); the two-start convergence is
the stability evidence (same behaviour documented in earlier rounds).

## 4. Laplace uncertainty quantification

Finite-difference Hessian of the pure NLL at the optimum; generalized
inverse with relative eigenvalue cutoff 1e-4 (effective dimension 11 of
12; one near-flat contact6/beta6 trade-off direction truncated); 1000
Monte Carlo draws, all equilibrium-feasible (1000/1000); Jeffreys
intervals for observed prevalence; log-Normal (sigma_pop = 0.10) for
observed population.

| Age | p_hat | p 95% CI | p obs 95% CI | overlap | N_hat | N 95% CI | N obs 95% CI | overlap |
|---|---:|---:|---:|:---:|---:|---:|---:|:---:|
| <20 | 0.106 | 0.062-0.187 | 0.060-0.184 | YES | 98.8 | 81-120 | 81-120 | YES |
| 20-29 | 0.173 | 0.157-0.199 | 0.153-0.195 | YES | 1255.5 | 1031-1537 | 1023-1513 | YES |
| 30-39 | 0.268 | 0.251-0.296 | 0.246-0.292 | YES | 1463.6 | 1232-1832 | 1206-1785 | YES |
| 40-49 | 0.430 | 0.410-0.455 | 0.408-0.453 | YES | 1880.4 | 1585-2326 | 1513-2240 | YES |
| 50-59 | 0.478 | 0.457-0.508 | 0.458-0.506 | YES | 1422.4 | 1250-1683 | 1338-1980 | YES |
| 60+ | 0.374 | 0.323-0.417 | 0.309-0.402 | YES | 459.5 | 407-558 | 336-498 | YES |

All six age groups overlap for both prevalence and population.

## 5. Sensitivity and rejected alternatives

- **12-parameter direct multi-start (run1)**: NLL 31.33; <20 prevalence
  0.181 and 40-49 population 2597 missed the criteria. Diagnosed as a
  missed basin; fixed with the run7-derived warm start (DECISIONS.md).
- **beta_scale > 1 constraint probe**:
  `output/calibration/probe_beta_constraint_run2/`. Forcing all six
  beta scales above 1 cannot satisfy the young-age population targets
  (see DECISIONS.md); the soft rule is abandoned for groups 1-2.
- **Excess-mortality contingency (eta_s[6])**: NOT activated — the 60+
  group is fitted within criteria with 12 parameters (p6 0.374 vs 0.3545,
  max |prev err| 0.0195 < 0.03; N6 459.5 vs 409, APE 0.126 < 0.20).

## 6. Structural limitations

1. J has no in-prison transmission; prison prevalence is inherited from D
   via arrest.
2. The 60+ seroprevalence decline (0.478 at 50-59 -> 0.374 at 60+)
   relies mainly on a large net uninfected 60+ inflow (beta6 ~ 110/yr);
   the mechanism is plausible but not directly measured in Singapore.
3. The 50-59 prison population is under-fitted (-12.6%) and 60+ slightly
   over-fitted (+12.3%); both within acceptance criteria but the largest
   residuals.
4. beta_scale < 1 for <20 and 20-29 (documented).
5. Progression rates are the highest defensible literature values; the
   GT3-weighted effective rates exceed Western pooled estimates
   (intentional conservative-upper choice; see natural_history_review.md).
6. 12 fitted parameters for 12 target summaries: metrics are calibration
   acceptance criteria, not a conventional goodness-of-fit test.
7. The equilibrium gate passes with margin, but the model is a
   quasi-steady state (2015-2030 drift <= 0.006 prevalence, < 0.3% total
   population).
8. Targets are treated as anti-HCV seroprevalence; viremic prevalence is
   reported as a sensitivity in `predictions.csv`.

## 7. Verdict

The 12-parameter model (highest-defensible progression rates, constant
transmission, 2015 mortality, no excess-mortality parameter) meets every
stated acceptance criterion, and the Laplace 95% intervals overlap the
observed intervals for all six age groups. No protected files were
modified (`Model schematic.pptx` hash unchanged; `src/sim.cpp` unchanged).

## 8. Reproduction

```bash
Rscript src/calibration/run_calibration.R \
  --run-id my_run --seed 101 --maxit 3000 --n-starts 6 \
  --t-start -10 --t-end 55 --target-time 45 --target-mode sero \
  --out-dir output/calibration

Rscript src/calibration/plot_results.R \
  --fit output/calibration/<run_id>/fit.rds --out-dir figures
```

All figures in `figures/` are ggplot2-generated (`fig01_prevalence_fit.png`
... `fig07_fit_dashboard.png`).
