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

## Current run (2026-10-21)

Model: `lambda1 = lambda3` (hard-coded), plain log entry scales, no entry likelihood, viraemic
excess mortality (HR 1.47). See `docs/lambda1_revert.md`.

### Deterministic calibration

`output/calibration/r_20261020_agepop/` gives the prior centre:
- best objective 18.73 (`perturb_2`, with Nelder–Mead restarts); seven of eight starts agree to
  within 0.013;
- seroprevalence RMSE 0.00004 and population MAPE 0.05%;
- the equilibrium gate passes.

**60+ seroprevalence:** 0.4286 against a target of 0.4286 (SE 0.041). The earlier misfit was
0.508 under `c_composite`.
- With free entry scales the HR hardly matters for the fit: 60+ is 0.440 at HR 1, 0.431 at 1.23,
  0.429 at 1.47 and 0.429 at 1.77 (`output/calibration/hr_sensitivity/`).

**Prior centre:** no coordinate of the best fit lies within 5 prior SDs of a bound, so
`prior_centre(best fit)` equals the best fit.

### Production posterior

`output/calibration/mcmc_reference/` holds the reference MCMC:
- **Sampler:** adaptive random-walk Metropolis with 4 chains × 50,000 evaluations.
  - Warmup is 10,000 proposals, with adaptation from proposal 5,000.
  - Every 10th of the 40,000 retained iterations is stored.
- **Prior:** Normal(prior centre, 0.20) on the 12 log-scale parameters, within the bounds.
- **Likelihood:** Binomial prevalence and log-normal population (σ = 0.10).
- **Equilibrium:** a hard restriction.
- **Diagnostics:** rank-normalised R-hat ≤ 1.002, bulk ESS 3,053–4,511, tail ESS ≥ 5,335, MCSE of the
  mean ≤ 0.018 posterior SD, and acceptance 0.26–0.28.
- **Prior-dominated parameters:** four parameters have posterior SD above 0.8 × the prior SD: θ6,
  θ9, θ10 and θ12.
- **2015 predictions** (4,000 draws, `model_prediction_2015_strata.*`):
  - ever-PWID 17,226 (95% interval: 15,649–18,968);
  - D 1,356, J 3,374, F 2,265 and X 10,235.

### Cross-validation (does not feed the production path)

`output/calibration/npe_lhs/` holds single-round NPE:
- **Training data:** 50,000 Latin hypercube prior simulations, of which 49,278 are valid (722
  equilibrium failures, 0 out of bounds). There are 12 summaries.
- **Network:** an equal-weight ensemble of 5 MAFs (50×5).
- **SBC** (980 datasets): cov95 0.944–0.961, cov90 0.891–0.918, all rank histograms uniform.
  - The gate (cov95 in [0.85, 1], cov90 in [0.75, 1], no skewed histogram) passes: `NPE_OK.txt`.
- **Accuracy against the reference** (`output/calibration/npe_vs_mcmc/report.md`):
  - all posterior means are within 0.03 SD;
  - 95% intervals are on average 3.5% narrower (ratio 0.92–0.99); 11 of 12 parameters pass all
    checks (θ4 misses the endpoint tolerance at 0.22 SD);
  - continuous policy outcomes are within 0.3%;
  - joint-pass probability at 5,000/yr F25% is 100% in both posteriors.

**Direct MCMC SBC** (`output/calibration/npe_vs_mcmc/mcmc_sbc/`, 8 datasets): pooled cov95 0.95,
cov90 0.92 and mean rank 0.50. This is a gross-failure check only.

### Scenarios (`output/scenarios_2026/`, 2,000 reference draws)

- **J only:** incidence −81.1% (95% UI 76.5–84.3; target met in 71% of draws) and mortality
  −50.5% (49.8–51.3).
- **Every finite scenario meets the incidence target** at the median (600/yr F25%: −84.1%).
- **The mortality target is met at the median** at 5,000/yr F25% (70.4%), 5,000/yr F50% (66.9%),
  8,000/yr F25% (73.3%), 8,000/yr F50% (71.0%) and no cap (73.9%). These meet both targets in ≥ 99.9%
  of draws. 8,000/yr F75% gives 64.9% (40% of draws).
- **No-treatment counterfactual (2030 vs 2015):** incidence 132.2% and deaths 101.1%.

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
