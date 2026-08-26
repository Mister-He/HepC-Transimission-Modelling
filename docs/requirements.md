# Requirements

## 1. Purpose and scope

The repository implements an age-structured compartmental model of hepatitis C
virus (HCV) transmission among people who inject drugs (PWID) in Singapore,
calibrated to the age-specific HCV prevalence and population of the prison
(detained) PWID stratum. The pipeline covers: model definition, deterministic
calibration (Nelder-Mead), uncertainty quantification (Laplace approximation
and Bayesian MCMC/NPE), and 50-year treatment/screening sensitivity
projections with a practical policy recommendation.

## 2. Model requirements

- **Structure**: 4 strata `D/J/F/X` (D = never-arrested active PWID,
  J = currently detained, F = ever-arrested active PWID, X = former),
  4 liver stages (NC, CC, DC, HCC), 4 HCV states (`u/a/c/t` =
  susceptible/post-SVR, acute, chronic, treatment), 6 age groups
  (<20, 20-29, 30-39, 40-49, 50-59, 60+), 384 compartments.
- **Transmission**: only active PWID (D and F) inject; J and X do not.
  Frequency-dependent force of infection via a 6x6 contact matrix.
- **Arrest dynamics**: D -> J at the first-arrest rate `lambda1`;
  J -> F with probability `pi_recid` or J -> X with `1 - pi_recid` at the
  release rate `lambda2`; F -> J at the re-arrest rate `lambda3`; new
  susceptibles enter `D_{u,1,i}` at age-specific inflow `beta[i]`.
- **Fixed baseline parameters**: original literature progression rates
  `0.027 / 0.039 / 0.014 / 0.014` (Thein 2008; Lim et al. 2018), 2015
  SingStat background mortality, PWID SMR `omega = 14.68`, constant
  transmission multiplier (`m_min = m_max = 1`), no excess-mortality
  parameter (`eta_s = 1`).
- **Simulation horizon**: `t_start = 0`, `t_end = 150` model years;
  targets compared at `t = 45` (calendar 2015); equilibrium assessed at
  `t = 150` vs `t = 145` on all 384 compartments and total HCV/DC/HCC.

## 3. Calibration requirements

- **Fitted parameters**: exactly 12 log-parameters — 6 contact-matrix row
  scales and 6 beta inflow scales (`contact_scale = exp(theta[1:6])`,
  `beta_scale = exp(theta[7:12])`). No upper bound on beta scales.
- **Targets**: age-specific prison (`J`) HCV prevalence and prison
  population. Prevalence likelihood is Binomial with integer counts
  `x_prev = round(prev_supplied * n_prev)`; population likelihood is
  log-Normal with `sigma_pop = 0.10`.
- **Optimizer**: deterministic multi-start Nelder-Mead (1 zero start +
  >= 5 perturbed starts); all starts recorded (initial/final theta,
  objective, convergence, function evaluations, elapsed time); best
  solution re-run from nearby perturbations.
- **Acceptance criteria**:

  ```text
  prevalence RMSE                    <= 0.02
  maximum absolute prevalence error  <= 0.03
  population MAPE                    <= 0.10
  maximum population APE             <= 0.20
  equilibrium                        pass
  Laplace 95% intervals              overlap observed intervals
  ```

- **Mortality rules**: use authoritative Singapore mortality data aligned
  with the target year (2015); age-stratified PWID SMR evidence mapped
  transparently; every mortality source records citation, setting, age
  bands, sex distribution, measure, and reason for inclusion.

## 4. Bayesian requirements

- **Primary method**: strict MCMC (adaptive Metropolis) after attempting
  NPE/SNPE; NPE is used for proposals and posterior checks but MCMC is the
  reported posterior.
- **Convergence**: R-hat in `[0.995, 1.005]` and pooled ESS `> 1000` for
  every parameter (extend chains until met); 3 chains with a wide thinning
  interval (`thin = 20`).
- **Visual review**: trace plots must show overlapping chains; density
  plots must be smooth, unimodal, close to normal on the log-parameter
  scale; per-chain posterior means must be close.
- **Calibration checks**: posterior predictive 95% credible intervals per
  age group must overlap the observed intervals and the Laplace intervals.

## 5. Sensitivity-analysis requirements

- **Scenario storage**: all strategies in one machine-readable table
  `src/calibration/scenarios.csv`; each row defines per-stage treatment
  rates `tau_*`, GT3 proportion `rho`, SVR/RBV mode, post-SVR modifier
  mode, eligible strata (`elig_D/J/F/X` 0/1 flags) and a minimum eligible
  age group (`min_age_group`). Strategies must include status quo,
  community-only, prison-only (all ages and 30+/40+/50+/60+), prison
  opt-out screening, aggressive treatment, post-release linkage, and
  prison + post-release continuation, plus GT3/SVR robustness variants.
- **Projection**: run each posterior draw to full equilibrium, use the
  equilibrated state as the 2017 initial condition, project 2017-2067;
  report total HCV (a+c+t), DC, and HCC per year with median and 95%
  equal-tailed credible intervals; status quo always included.
- **Figures**: pure white background, no grid lines, 95% CrIs, legend on
  the right ordered as in `scenarios.csv`.
- **Recommendation**: justify the recommended screening/treatment strategy
  with practical considerations (prison coverage, age prioritisation,
  community reach, linkage, cost, equity) and cite literature/websites for
  every claim.

## 6. Engineering and reproducibility requirements

- Standard layout: `src/` (model + analysis library), `scripts/` (pipeline
  entry points), `tests/unit` and `tests/integration`, `docs/`
  (`requirements.md`, `architecture.md`, `test-plan.md`, calibration
  reports), `.github/workflows/ci.yml`, optional `Dockerfile` /
  `docker-compose.yml`.
- Every calibration run records run ID, Git SHA, hashes of `setup.R`,
  `sim.cpp` and the model schematic, targets, mortality scenario,
  likelihood/optimizer settings, fitted values, predictions, metrics, and
  R/package versions.
- `Model schematic.pptx` is protected and never modified; `src/sim.cpp`
  is modified only when strictly necessary and its 4-strata/384-compartment
  structure is preserved.
- CI runs unit and integration tests on every push/PR; the pipeline must
  pass locally and in CI without external secrets.
