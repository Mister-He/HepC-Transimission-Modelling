# Evaluation

Evaluate the model in four layers:

1. Age-specific seroprevalence residuals (the only calibration target).
2. Modelled prison and stratum sizes, reported against the observed 2015 prison count but not
   fitted.
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

## Current run (2026-10-22)

Model (`docs/first_arrest_multiplier.md`): first arrest = `c` × re-arrest (`c` ≥ 0.5, calibrated),
CNB entry unscaled, mortality `mu` × PWID SMR (9.52 viraemic, 4.73 otherwise; van Santen 2014),
seroprevalence-only Beta likelihood (φ = 100), fixed prior (contact N(0, 1.5), log `c` N(0, 0.5)).
This run includes the audit fixes: initial state shifted to the right compartments, MCMC
adaptation from the recent half of warmup, HCV deaths including the viraemic non-liver excess, and
prison testing by the F/X rule.

### Deterministic calibration

`output/calibration/r_20261021_agepop/` gives the posterior mode:
- all 8 uninformative starts reach objective −7.83 (negative log posterior);
- seroprevalence RMSE 0.013, largest error 0.023 (50–59);
- `c` 1.10; the equilibrium gate passes.
- **SMR sensitivity** (`output/calibration/smr_sensitivity/`): at the CI bounds `c` is 1.097–1.100
  and RMSE 0.012–0.015.

### Production posterior

`output/calibration/mcmc_reference/` holds the reference MCMC:
- **Sampler:** adaptive random-walk Metropolis with 4 chains × 50,000 evaluations.
  - Warmup is 10,000 proposals, with adaptation from proposal 5,000 on the recent half of warmup.
  - Every 10th of the 40,000 retained iterations is stored.
- **Diagnostics:** rank-normalised R-hat ≤ 1.002, bulk ESS 3,250–4,645, tail ESS ≥ 3,487, MCSE of the
  mean ≤ 0.019 posterior SD, and acceptance 0.18–0.22.
- **Prior-dominated parameters:** θ5, θ6 and θ7 have posterior SD above 0.7 × the prior SD.
  `c` = 1.24 (95% interval 0.57–3.13), SD ratio 0.87: seroprevalence barely informs it.
- **2015 predictions** (4,000 draws, `model_prediction_2015_strata.*`):
  - ever-PWID 50,218 (95% interval: 49,496–50,893);
  - D 2,178, J 8,451, F 5,857 and X 33,713.
  - **Not fitted:** the modelled prison population (8,451) is 2.5 times the observed 2015 count of
    screened PWID (3,367, `output/prison_2015_counts.csv`). Absolute counts and costs are uncertain
    beyond the UIs.

### Cross-validation (does not feed the production path)

`output/calibration/npe_lhs/` holds single-round NPE:
- **Training data:** 50,000 Latin hypercube prior simulations, of which 49,945 are valid (43
  equilibrium failures, 12 non-finite). There are 6 summaries.
- **Network:** an equal-weight ensemble of 5 MAFs (50×5).
- **SBC** (998 datasets): cov95 0.933–0.960, cov90 0.872–0.911, all rank histograms uniform; gate PASS.
- **Accuracy against the reference** (`output/calibration/npe_vs_mcmc/report.md`): poor at the
  observed data. 33% of NPE draws fail the equilibrium check and are dropped; of the rest, 0 of 7
  parameters pass all checks (worst |SMD| 0.56, 95% intervals on average 28% wider). Policy
  reductions are within 0.3%; the joint-pass probability at 5,000/yr F25% differs by 11 pp.

**Direct MCMC SBC** (`output/calibration/npe_vs_mcmc/mcmc_sbc/`, 8 datasets): pooled cov95 0.96,
cov90 0.89 and mean rank 0.48. This is a gross-failure check only.

### Scenarios (`output/scenarios_2026/`, 2,000 reference draws)

- **Entry** is held at the CNB 2015 counts in every projection year (2026-10-23; previously the
  2016–2025 CNB series, then 2025 held constant).
- **J only:** incidence −84.5% (95% UI 68.7–90.5; target met in 77% of draws) and mortality
  −54.9% (52.3–57.5).
- **Every finite scenario meets the incidence target** at the median (600/yr F25%: −85.6%).
- **The mortality target is met at the median** at 5,000/yr F25% (68.4%; both targets in 95% of
  draws), 5,000/yr F50% (65.1%; 51%), 8,000/yr F25% (73.9%; 98%), 8,000/yr F50% (69.6%; 99%) and
  no cap (87.1%; 99.7%). With 75% to F it is missed (5,000/yr 61.0%, 8,000/yr 63.7%).
- **Allocation:** at every capacity, more tests for X give larger mortality falls (F holds 6%
  chronic among those eligible, X 9%; X holds 90% of chronic infections after release).
- **No-treatment counterfactual (2030 vs 2015):** incidence 100.0% and deaths 100.1% (constant entry
  keeps the untreated epidemic at its 2015 equilibrium).

### Calibration figure

`figures/calibration_agepop.png` and `figures/calibration_interval_overlap.csv` show that all 6
seroprevalence target ranges overlap the 95% intervals of the MCMC posterior. This overlap is a
calibration diagnostic. It does not validate prospective policy predictions.

### Outcome definitions

- **Incidence:** the instantaneous annualised flow of new and repeat infections from `u` and `r` into
  `a` in D and F.
- **Mortality:** the annual HCV death flow: DC/HCC liver deaths plus the non-liver excess of
  viraemic people, `mu × (smr_viraemic − smr_nonviraemic)` per person in states a, c, t.

These are event counts per year, not rates per person. Treatment starts are counted as courses.

### What these results do not show

These checks concern the fitted PWID model. They are not estimates for Singapore's general
population. The 2019-onward prison treatment assumption is not fitted to observed treatment counts.
