# HCV transmission model for people who inject drugs (PWID) in Singapore

An age-structured compartmental model of hepatitis C virus (HCV) among
people who inject drugs (PWID) in Singapore, calibrated to the age-specific
HCV seroprevalence and population of the prison (detained) PWID stratum.

## 1. Project structure

```text
.
├── README.md / README.zh-CN.md  Project overview (EN / 中文)
├── AGENTS.md / prompt*.md       Working instructions for agents
├── Model schematic.pptx         Protected visual model description
├── docs/
│   ├── requirements.md          Functional & non-functional requirements
│   ├── architecture.md          Model/component architecture + data flow
│   ├── test-plan.md             Unit/integration testing strategy
│   └── calibration/             Calibration & analysis reports
│       ├── final_report.md
│       └── analysis_report.md
├── src/
│   ├── setup.R                # Parameters, initial conditions, helpers
│   ├── sim.cpp                # Compiled ODE solver (4 strata, 384 compartments)
│   └── calibration/
│       ├── targets.R          # Target data + Binomial count reconstruction
│       ├── model_metrics.R    # J summary extraction + fit metrics
│       ├── equilibrium.R      # T vs T-5 stability gate (all states)
│       ├── likelihood.R       # 12-parameterisation + NLL objective
│       ├── calibrate_nm.R     # Multi-start Nelder-Mead runner
│       ├── laplace.R          # Laplace approximation (Hessian, MC intervals)
│       ├── mcmc.R             # Adaptive Metropolis + priors
│       ├── plot_results.R     # Calibration figures (library + CLI)
│       ├── plot_mcmc.R        # Bayesian figures (library)
│       ├── scenarios.csv      # Sensitivity scenario inventory (single source)
│       └── (library modules)  # Sourced by scripts/ entry points
├── scripts/                   # Pipeline entry points
│   ├── run_calibration.R      # NM calibration end-to-end
│   ├── run_npe.R              # NPE + MCMC Bayesian pipeline
│   ├── npe_train.py           # Python NPE trainer (sbi/torch)
│   ├── run_analysis.R         # 50-year sensitivity projections (posterior CrIs)
│   ├── plot_results.R         # CLI wrapper for calibration figures
│   ├── probe_beta_constraint.R# beta_scale > 1 feasibility probe
│   └── run_tests.R            # Unit/integration test runner
├── tests/
│   ├── helper.R               # Assertion helpers
│   ├── unit/                  # Fast dependency-light unit tests
│   └── integration/           # Compile + simulate + pipeline smoke tests
├── figures/                   # Final ggplot2 figures
├── output/                    # Generated artefacts
│   ├── calibration/           # Per-run outputs (run_config, solutions, ...)
│   └── analysis/              # Sensitivity projections
├── .github/workflows/ci.yml   # CI (unit + integration tests)
├── Dockerfile                 # Optional reproducible R environment
└── docker-compose.yml         # Optional containerised workflows
```

## 2. Model structure

```text
4 strata:  D (never-arrested active), J (currently detained),
           F (ever-arrested active), X (former, non-injecting)
4 liver stages: NC, CC, DC, HCC
4 HCV states:  u (susceptible/post-SVR), a (acute), c (chronic), t (treatment)
6 age groups:  <20, 20-29, 30-39, 40-49, 50-59, 60+
384 compartments
```

- **Transmission** occurs only among active PWID (D and F inject; J and X
  do not): frequency-dependent force of infection
  `gamma_i = q * sum_j C_contact(i,j) * infectious_active_j / active_j`,
  with a **constant** transmission multiplier (m = 1; the historical level
  is absorbed by the fitted contact row scales).
- **Ageing**: 10-year bands, `y_change = y/10 per year` to the next band
  (60+ open-ended).
- **Arrest/release**: D -> J (first arrest) at `lambda1[i]`; J -> F with
  probability `pi_recid` or J -> X with `1 - pi_recid` at `lambda2[i]`;
  F -> J (re-arrest) at `lambda3[i]`.
- **Mortality**: `mu[i] * omega` (omega = 14.68 PWID SMR);
  decompensated cirrhosis adds `mu_DC = 0.130/yr`, HCC adds
  `mu_HCC = 0.430/yr`.
- **Progression rates** (original, literature-based): `p_NC_CC = 0.027`
  (Thein 2008); `p_CC_DC = 0.039`, `p_CC_HCC = 0.014`, `p_DC_HCC = 0.014`
  (Lim et al. 2018). GT3 relative risks (Kanwal et al. 2014) retained.

## 3. Data sources

- **Targets**: age-specific anti-HCV seroprevalence and prison population
  from the Changi Prison universal screening 2014-2016 (source paper:
  Park et al., medRxiv 10.1101/2025.10.24.25338708).
- **Background mortality**: Singapore SingStat age-specific death rates,
  2015, mapped to the six age groups (see `mortality_review.md`).
- **PWID excess mortality**: Mathers et al. 2013 pooled SMR 14.68
  (PMID 23554523).
- **Progression rates**: Thein et al. 2008 (NC -> CC); Lim et al. 2018
  (CC -> DC, CC -> HCC, DC -> HCC). Alazawi et al. 2010 and
  Rivera-Irigoin et al. 2006 were used in earlier higher-rate variants
  (see PROJECT_OVERVIEW.md).

## 4. Calibration

- **12 fitted parameters**: 6 contact-matrix row scales, 6 beta inflow
  scales (all `exp(theta)`). No transmission multiplier, no
  excess-mortality parameter (the 60+ group is fitted within the
  acceptance criteria with these 12).
- Objective: Binomial NLL for prevalence + log-Normal NLL for population
  (sigma_pop = 0.10); equilibrium gate (T vs T-5) as a penalty.
- Optimiser: deterministic multi-start Nelder-Mead
  (`optim`, method = "Nelder-Mead"). No HMC.
- Uncertainty: Laplace approximation at the optimum (finite-difference
  Hessian, generalized inverse with eigenvalue cutoff 1e-4, 1000 Monte
  Carlo draws filtered by equilibrium feasibility), 95% intervals compared
  with observed intervals per age group.
- Acceptance criteria: prevalence RMSE <= 0.02, max |prev err| <= 0.03,
  population MAPE <= 0.10, max APE <= 0.20, equilibrium pass.
- `beta_scale > 1` is a soft target: satisfied for age groups 3-6; the
  <20 and 20-29 groups remain below 1 because their small prison
  populations require net inflows below the CNB-anchored base `beta`
  (see DECISIONS.md).

### 4.1 Bayesian inference (MCMC)

On top of the NM fit, Bayesian posterior sampling is performed. Per
`prompt_mcmc.md` (revised), **NPE/SNPE was attempted first** (sbi `NPE_C`,
3 sequential rounds, two seeds, SBC-calibrated) but was not adopted: it
was overconfident in well-identified directions and unstable across seeds
in weakly identified directions. The **traditional MCMC (adaptive
Metropolis-Hastings) is the primary Bayesian summary**, with strict
convergence: R-hat in [0.995, 1.005] (achieved 0.9997-1.0041), pooled
ESS > 1000 (achieved 1,310-1,958), 3 chains x 30,000 iterations, burn-in
5,000, thin 20. Priors: log contact scales ~ Normal(0, 2^2); log beta
scales ~ Student-t(3, 0, 2). Posterior predictive 95% credible intervals
per age group overlap both the Laplace and observed intervals. See
`prompt_mcmc.md`.

## 5. How to regenerate

```bash
# 1. Full calibration run (defaults reproduce the final run)
Rscript scripts/run_calibration.R \
  --run-id my_run --seed 101 --maxit 3000 --n-starts 6 \
  --t-start 0 --t-end 150 --target-time 45 --target-mode sero \
  --out-dir output/calibration

# 2. Publication figures (ggplot2)
Rscript scripts/plot_results.R \
  --fit=output/calibration/<run_id>/fit.rds --out-dir=figures

# 3. Bayesian posterior (NPE 3 rounds + MCMC strict validation)
Rscript scripts/run_npe.R --step all \
  --root . --fit output/calibration/<run_id>/fit.rds \
  --out-dir output/calibration/npe_bayes \
  --n-sims 60000 --n-cores 6 --seed 2026 \
  --n-draws 60000 --n-proposal 20000 --n-rounds 3 \
  --n-iter-mcmc 20000 --burnin-mcmc 5000 --mcmc-max-iter 400000

# 4. Tests (unit + integration)
Rscript scripts/run_tests.R
```

Per-run outputs: `run_config.csv`, `targets.csv`, `initial_values.csv`,
`optimization_history.csv`, `solutions.csv`, `predictions.csv`,
`residuals.csv`, `equilibrium.csv`, `laplace_intervals.csv`,
`laplace_diagnostics.csv`, `sessionInfo.txt`, `fit.rds`, and ggplot2
figures.

### 5.1 Sensitivity analysis (50-year projections from 2017)

Scenario strategies are defined in one CSV
(`src/calibration/scenarios.csv`): treatment rates by liver stage, GT3
proportion, SVR/RBV mode, eligible strata (D/J/F/X flags), and a minimum
eligible age group. Each row is one strategy; outputs are medians and 95%
credible intervals over 300 posterior draws:

```bash
Rscript scripts/run_analysis.R \
  --root . --fit output/calibration/run1_4strata/fit.rds \
  --posterior output/calibration/npe_bayes/posterior_samples_mcmc.csv \
  --out-dir output/analysis --n-draws 300 --n-cores 4
```

Outputs: `scenario_summary.csv`, `scenario_key_years.csv`,
`fig_hcv_trajectories.png`, `fig_dc_hcc_trajectories.png`.

## 6. Results

Final run: `output/calibration/run1_4strata/`.

Acceptance criteria are met with 12 parameters and no excess mortality.
Bayesian outputs are in `output/calibration/npe_bayes/`; 50-year
sensitivity projections are in `output/analysis/`. See
`docs/calibration/final_report.md` and `docs/calibration/analysis_report.md`.
