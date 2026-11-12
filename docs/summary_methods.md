# Methods summary (state of 2026-11-12)

Short version of the Methods and S1 Text. Commands and run times are in `docs/training.md`.

## 1. Model
- Deterministic age-structured compartmental model of PWID in Singapore (`scripts/model/sim.cpp`,
  called from R through Rcpp): 4 groups (D, J, F, X) × 4 liver stages (non-cirrhotic, compensated
  cirrhosis, DC, HCC) × 5 HCV states (susceptible, acute, chronic, on treatment, recovered) × 6 age
  groups = 480 compartments.
- RK4 with a daily step; 150-year run from a seed population to a 2015 steady state, accepted only if the
  last 10 years are stable (equilibrium check).
- Force of infection from an age contact matrix with a calibrated multiplier per age group, acting only
  in D and F.

## 2. Inputs
- **Prison rates** (release, re-arrest, π) from Singapore Prison Service records of 20,720 people,
  Bayesian survival models in Stan (cure-type model for re-arrest).
- **Entry** from CNB 2015 new drug abusers × heroin/methamphetamine share (0.861).
- Natural history, treatment efficacy, mortality from the literature (S1 Text Table A).

## 3. Calibration and validation
- 13 parameters (6 contact multipliers, 6 entry scales, first-arrest factor c), fixed priors.
- Likelihood: beta for prison seroprevalence by age (2015), negative binomial for first arrests by age,
  log-normal for the community PWID estimate (Teo 2019, 11,000).
- Multistart Nelder–Mead for the posterior mode, then adaptive random-walk Metropolis (8 chains;
  R-hat ≤ 1.005, bulk ESS ≥ 1,641). Simulation-based calibration on simulated data sets.
- Validation: all 13 posterior predictive intervals overlap their targets; the unfitted 2015 prison
  population (8,768) matches the number held (8,725); unit and integration tests, including an
  independent R implementation of the model equations. No observed series of infections or deaths
  among PWID in Singapore exists for external validation.

## 4. Policy scenarios
- Main analysis: no treatment before 2026; from 2026 prison treatment (annual start probability
  99.995%) plus 17 policies: J only (status quo), 15 finite community screening policies (600, 1,000,
  2,000, 5,000, 8,000 tests a year × 25/50/75% to F, rest to X), and uncapped screening of all F and X.
- Each year's tests find chronic infections in proportion to HCV states; those found set the year's
  treatment capacity in F or X (× linkage share in the linkage sensitivity analysis).
- Comparators: no treatment (counterfactual for outcomes averted) and J only (status quo for
  incremental costs and outcomes).
- Prison coverage grid: 0%, 50%, 90%, 99.995% (68 combinations).
- Sensitivity analyses: prison treatment from 2019; observed 2016–2025 inflow; genotype 3 share 10% /
  50%; one-way variation of 27 fixed inputs at the posterior mode (tornado); linkage 80% / 60%.
  The first four run to 2030; coverage and linkage run to 2050.
- Every scenario run for each of 2,000 posterior draws (`scripts/scenarios_2026.R`), summarised as
  median and 95% UI (2.5th–97.5th percentiles).

## 5. Outcomes
- Primary (WHO 2030 elimination targets): reduction in 2030 vs 2015 of yearly new infections (target
  80%) and yearly liver (DC/HCC) deaths (target 65%); yearly counts by trapezoid rule on daily output.
  Share of draws meeting both targets.
- Secondary: HCV deaths including the extrahepatic excess of viraemic people.
- Long horizon: projection to 2050 (`end_year`); cumulative new infections and liver deaths over
  2026–2030 and 2026–2050, averted vs no treatment and vs J only, per draw (not discounted).
- Programme budget (testing, consultations, DAA; discounted 3% to 2026), treatment courses, incremental
  budget vs J only, and cost per infection / liver death averted vs J only (NA if < 1 averted).

## 6. Strategy selection (`scripts/decision_analysis.R`)
- Strategy sets: the 17 policies (s17) and the 68 policy × prison-coverage combinations (s68).
- Criteria per horizon: infections averted, liver deaths averted (vs no treatment), P(both 2030
  targets) — benefits; programme budget, treatment courses — costs. Min–max normalised to 0–1 within the
  set, costs reverse-coded.
- Rules: (a) most benefit within a budget cap; (b) lowest cost per outcome averted vs no treatment;
  (c) TOPSIS among Pareto-efficient strategies (equal weights). Cost-effectiveness frontier per benefit
  with strong and extended dominance.
- Uncertainty: Pareto, TOPSIS and the budget rule repeated within each posterior draw → P(Pareto),
  P(ranked first), P(best within cap). Weight sensitivity: each weight doubled or zeroed, benefits only,
  costs only, 1,000 random weight vectors from the simplex.
- Main decision horizon 2026–2050.

## 7. Reporting
- Manuscript built from outputs by `scripts/build_paper.py` (PLOS Medicine, CHEERS 2022); every number
  is read from `output/` and `checks()` stops the build if a stated result changes.
- Figures by `scripts/figures.py` (matplotlib).
