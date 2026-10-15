# Evaluation

Evaluate the model in four layers:

1. Age-specific seroprevalence residuals.
2. Age-specific prison-population residuals.
3. All-compartment equilibrium stability over the last 10 of the 150 model years.
4. Posterior quality:
   - **Reference MCMC:** convergence and Monte Carlo error.
   - **Single-round NPE (cross-validation only):** SBC coverage and rank uniformity, and agreement
     with the reference MCMC.

**Projection gate.** The 2015 snapshot must pass the equilibrium gate before projection.

**Scenario screen.** The primary screen is the WHO HCV impact target. In 2030, incidence must be at
or below 20% of the 2015 baseline (an 80% reduction), and HCV-related mortality at or below 35% (a
65% reduction). Ever-PWID viraemic prevalence at or below 1% is a secondary diagnostic only. It is
not evidence of national elimination.

## Current run (2026-10-15)

### Deterministic calibration

`output/calibration/r_20261014_agepop/` gives the prior centre:
- best objective 21.44 (`population_informed`, with Nelder–Mead restarts);
- seroprevalence RMSE 0.012 and population MAPE 5.7%;
- the equilibrium gate passes.

**60+ seroprevalence:** 0.455 against a target of 0.429 (SE 0.041). The earlier misfit was
0.508 under `c_composite`.
- With the CNB entry anchor, the chronic viraemic HR is what closes most of the gap: 60+ is
  0.494 at HR 1, 0.481 at 1.23, 0.455 at 1.47 and 0.439 at 1.77
  (`output/calibration/hr_sensitivity/`).

**First-arrest flow:** at the best fit, the model's 2015 first-arrest flow is within 0.04–2.7% of
Stan's `f_first`. The shortfall, from HCV-related deaths in D, is within Stan's posterior SD
(3–19%).

**Prior centre:** `prior_centre(best fit)` moves θ3 and θ6 (the 30–39 and 60+ contact scales,
near their lower bound) to 5 prior SDs inside the bound.

### Production posterior

`output/calibration/mcmc_reference/` holds the reference MCMC:
- **Sampler:** adaptive random-walk Metropolis with 4 chains × 50,000 evaluations.
  - Warmup is 10,000 proposals, with adaptation from proposal 5,000.
  - Every 10th of the 40,000 retained iterations is stored.
- **Prior:** Normal(prior centre, 0.20) on the 12 log-scale parameters, within the bounds.
- **Likelihood:** Binomial prevalence, log-normal population and log-normal 2015 entry against
  the CNB counts (σ = 1.0).
- **Equilibrium:** a hard restriction.
- **Diagnostics:** rank-normalised R-hat ≤ 1.003, bulk ESS 2,319–3,495, tail ESS ≥ 2,728, MCSE of the
  mean ≤ 0.021 posterior SD, and acceptance 0.25–0.27.
- **Prior-dominated parameters:** eight parameters have posterior SD above 0.8 × the prior SD: θ5,
  θ6 and θ7–θ12. The prison data mainly inform the contact scales of ages 20–49.
- **2015 predictions** (4,000 draws, `model_prediction_2015_strata.*`):
  - ever-PWID 18,453 (95% interval: 17,799–19,380);
  - D 3,433, J 3,253, F 2,157 and X 9,608.

### Cross-validation (does not feed the production path)

`output/calibration/npe_lhs/` holds single-round NPE:
- **Training data:** 50,000 Latin hypercube prior simulations, of which 48,918 are valid (1,082
  equilibrium failures, 0 out of bounds). There are 18 summaries, including the 2015 entry.
- **Network:** an equal-weight ensemble of 5 MAFs (50×5).
- **SBC** (973 datasets): cov95 0.936–0.958, cov90 0.879–0.918, all rank histograms uniform.
  - The gate (cov95 in [0.85, 1], cov90 in [0.75, 1], no skewed histogram) passes: `NPE_OK.txt`.
- **Accuracy against the reference** (`output/calibration/npe_vs_mcmc/report.md`):
  - all posterior means are within 0.08 SD;
  - 95% intervals match within 4% (mean ratio 0.99), and all 12 parameters pass all checks;
  - continuous policy outcomes are within 1.2%;
  - joint-pass probabilities are 0 in both posteriors.

**Direct MCMC SBC** (`output/calibration/npe_vs_mcmc/mcmc_sbc/`, 8 datasets): pooled cov95 0.93,
cov90 0.83 and mean rank 0.51. This is a gross-failure check only.

### Scenarios (`output/scenarios_2026/`, 2,000 reference draws)

- **J only:** incidence −47.6% (95% UI 36.5–57.5) and mortality −47.4% (44.7–48.8).
- **No scenario meets the incidence target.** Even the no-cap reference gives −57.3%
  (44.4–68.4). In 2030, 59% of its new infections are in D (never imprisoned), which prison and
  F/X screening do not reach.
- **The mortality target is met at the median** at 8,000/yr F25% (67.2%; 87% of draws), 8,000/yr
  F50% (65.4%; 60%) and no cap (67.4%; 89%). 5,000/yr F25% gives 64.9% (47% of draws).
- **No-treatment counterfactual (2030 vs 2015):** incidence 101.9% and deaths 100.1%.

### Calibration figure

`figures/calibration_agepop.png` and `figures/calibration_interval_overlap.csv` show that all 12 target
intervals overlap the 95% predictive intervals of both posteriors. This overlap is a calibration
diagnostic. It does not validate prospective policy predictions.

### Outcome definitions

- **Incidence:** the instantaneous annualised flow of new and repeat infections from `u` and `r` into
  `a` in D and F.
- **Mortality:** the DC/HCC HCV-attributed annual death flow.

These are event counts per year, not rates per person. Treatment starts are counted as courses.

### What these results do not show

These checks concern the fitted PWID model. They are not estimates for Singapore's general
population. The 2019-onward prison treatment assumption is not fitted to observed treatment counts.
