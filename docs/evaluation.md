# Evaluation

Evaluate the model in four layers:

1. Age-specific seroprevalence residuals.
2. Age-specific prison-population residuals.
3. All-compartment equilibrium stability between model years 140 and 150.
4. Posterior quality:
   - **Reference MCMC:** convergence and Monte Carlo error.
   - **Single-round NPE:** agreement with the reference MCMC, posterior spread, and predictive
     coverage.

**Projection gate.** The 2015 snapshot must pass the equilibrium gate before projection.

**Scenario screen.** The primary screen is the WHO HCV impact target. In 2030, incidence must be at
or below 20% of the 2015 baseline (an 80% reduction), and HCV-related mortality at or below 35% (a
65% reduction). Ever-PWID viraemic prevalence at or below 1% is a secondary diagnostic only. It is
not evidence of national elimination.

## Current run

### Deterministic calibration

`output/calibration/r_20261010_agepop/` holds the prior centre:
- best objective 19.01;
- seroprevalence RMSE 0.0044 and population MAPE 2.3%;
- the equilibrium gate passes.

### Production posterior

`output/calibration/mcmc_reference/` holds the reference MCMC:
- **Sampler:** adaptive random-walk Metropolis with 4 chains × 45,000 retained draws.
- **Prior:** Normal(best fit, 0.20) on the 12 log-scale parameters.
- **Likelihood:** Binomial prevalence and log-normal population.
- **Equilibrium:** a hard restriction.
- **Diagnostics:** rank-normalised R-hat ≤ 1.002, bulk ESS 3,433–4,600, tail ESS 6,242–8,994, and
  MCSE of the mean ≤ 0.017 posterior SD.
- **2015 predictions** (4,000 draws, `model_prediction_2015_strata.*`):
  - ever-PWID 21,224 (95% interval: 19,226–23,437);
  - F 1,532 and X 14,917.

### Independent check

`output/calibration/npe_lhs/` holds single-round NPE:
- **Training data:** 50,000 Latin hypercube prior simulations, of which 49,765 were valid; 235
  failed equilibrium.
- **Network:** MAF 50×5.
- **Quality gate:** `NPE_OK.txt` (posterior spread). SBC is not used for this check.
- **Accuracy against the reference** (`output/calibration/npe_vs_mcmc/report.md`):
  - all posterior means are within 0.1 SD;
  - 95% intervals are about 4% narrower than the reference, and 11 of 12 parameters pass all checks;
  - policy outcomes are within 0.4%.

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
