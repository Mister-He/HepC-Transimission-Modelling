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

## Current run (2026-10-14)

### Deterministic calibration

`output/calibration/r_20261013_agepop/` gives the prior centre:
- best objective 22.58 (`perturb_5`);
- seroprevalence RMSE 0.033 and population MAPE 5.3%;
- the equilibrium gate passes.

The 60+ seroprevalence is 0.508 against a target of 0.429. This is a structural misfit: with the
composite re-arrest rate, first arrests barely depend on `beta_scale`.

The prior centre is `prior_centre(best fit)`. It moves θ6 (the 60+ contact scale, at its lower
bound in the fit) to 5 prior SDs inside the bound.

### Production posterior

`output/calibration/mcmc_reference/` holds the reference MCMC:
- **Sampler:** adaptive random-walk Metropolis with 4 chains × 50,000 evaluations.
  - Warmup is 10,000 proposals, with adaptation from proposal 5,000.
  - Every 10th of the 40,000 retained iterations is stored.
- **Prior:** Normal(prior centre, 0.20) on the 12 log-scale parameters, within the bounds.
- **Likelihood:** Binomial prevalence and log-normal population.
- **Equilibrium:** a hard restriction.
- **Diagnostics:** rank-normalised R-hat ≤ 1.002, bulk ESS 2,824–4,259, tail ESS ≥ 5,119, MCSE of the
  mean ≤ 0.019 posterior SD, and acceptance 0.26–0.30.
- **Prior-dominated parameters:** θ6 and θ12 (posterior SD ≈ prior SD).
- **2015 predictions** (4,000 draws, `model_prediction_2015_strata.*`):
  - ever-PWID 21,169 (95% interval: 19,152–23,399);
  - F 1,551 and X 15,016.

### Cross-validation (does not feed the production path)

`output/calibration/npe_lhs/` holds single-round NPE:
- **Training data:** 50,000 Latin hypercube prior simulations, all valid (0 out of bounds, 0
  equilibrium failures).
- **Network:** an equal-weight ensemble of 5 MAFs (50×5).
- **SBC** (1,000 datasets): cov95 0.931–0.958, cov90 0.873–0.910, all rank histograms uniform.
  - The gate (cov95 in [0.85, 1], cov90 in [0.75, 1], no skewed histogram) passes: `NPE_OK.txt`.
- **Accuracy against the reference** (`output/calibration/npe_vs_mcmc/report.md`):
  - all posterior means are within 0.05 SD;
  - 95% intervals are about 4% narrower, and all 12 parameters pass all checks;
  - continuous policy outcomes are within 0.6%;
  - the 5,000/yr F25% joint-pass probability is 8.7% against 11.25%.

**Direct MCMC SBC** (`output/calibration/npe_vs_mcmc/mcmc_sbc/`, 8 datasets): pooled cov95 1.00 and
cov90 0.95. This is a gross-failure check only.

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
