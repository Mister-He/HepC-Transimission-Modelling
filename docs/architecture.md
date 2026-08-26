# Architecture

## 1. Directory layout

```text
.
├── README.md / README.zh-CN.md    Project overview (EN/中文)
├── AGENTS.md / prompt*.md         Working instructions for agents
├── Model schematic.pptx           Protected visual model description
├── docs/
│   ├── requirements.md            Functional & non-functional requirements
│   ├── architecture.md            This document
│   ├── test-plan.md               Testing strategy
│   └── calibration/               Calibration reports and decision log
├── src/
│   ├── setup.R                    Parameters, initial conditions, helpers
│   ├── sim.cpp                    Compiled ODE solver (RK4, 384 compartments)
│   └── calibration/               Analysis library (sourced by scripts)
│       ├── targets.R              Target data + Binomial count reconstruction
│       ├── model_metrics.R        J (prison) summary extraction + fit metrics
│       ├── equilibrium.R          T vs T-5 stability gate (all states)
│       ├── likelihood.R           12-parameter objective + penalties
│       ├── calibrate_nm.R         Multi-start Nelder-Mead runner
│       ├── laplace.R              Laplace approximation (Hessian, MC intervals)
│       ├── mcmc.R                 Adaptive Metropolis + priors
│       ├── plot_results.R         Calibration figures (library + CLI)
│       ├── plot_mcmc.R            Bayesian figures (library)
│       └── scenarios.csv          Sensitivity scenario inventory
├── scripts/                       Pipeline entry points
│   ├── run_calibration.R          NM calibration end-to-end
│   ├── run_npe.R                  NPE + MCMC Bayesian pipeline
│   ├── npe_train.py               Python NPE trainer (sbi/torch)
│   ├── run_analysis.R             50-year scenario projections
│   ├── plot_results.R             CLI wrapper for calibration figures
│   ├── probe_beta_constraint.R    Exploratory beta_scale > 1 probe
│   └── run_tests.R                Unit/integration test runner
├── tests/
│   ├── helper.R                   Assertion helpers
│   ├── unit/                      Fast, dependency-light unit tests
│   └── integration/               Compile + simulate + pipeline smoke tests
├── figures/                       Stable publication figures
├── output/                        Generated calibration/analysis artefacts
├── .github/workflows/ci.yml       CI pipeline
├── Dockerfile                     Optional reproducible R environment
└── docker-compose.yml             Optional containerised workflows
```

## 2. Component responsibilities

### `src/setup.R`
Defines `params` (all fixed model parameters: transmission, natural
history, treatment baseline `tau = 0`, arrest rates `lambda1..3`, mortality,
initial contact matrix and inflow `beta`) and `data` (horizon `0..150`,
`dt = 1/365`, initial conditions). It compiles `sim.cpp` when sourced.

### `src/sim.cpp`
Continuous-time ODE system (RK4) over 384 compartments:

```text
stratum s ∈ {D=0, J=1, F=2, X=3}
stage   k ∈ {NC=1, CC=2, DC=3, HCC=4}
state   h ∈ {u=0, a=1, c=2, t=3}
age     i ∈ {0..5}
```

Flat index: `idx(s,k,h,i) = s*96 + (k-1)*24 + h*6 + i`; the output matrix
adds a leading time column. Treatment eligibility is stratified:
`treat(s) = tau_stratum[s] * (i+1 >= tau_min_age)`, so a scenario can
target arbitrary strata/ages. J and X do not inject, so they contribute no
force of infection.

### `src/calibration/` (library)

- `targets.R` — observed age-specific prevalence and prison population;
  reconstructs integer Binomial counts.
- `model_metrics.R` — extracts J-summaries (`N_hat`, `p_sero`,
  `p_viremic`) from the output matrix and computes fit metrics.
- `equilibrium.R` — `check_equilibrium()` compares `t = T` vs `T-5` on J
  summaries, total population, total HCV/DC/HCC and all 384 compartments.
- `likelihood.R` — `build_params()` (exp parameterisation), `make_objective()`
  (Binomial prevalence + log-Normal population + soft plausibility
  penalty + equilibrium penalty), bounds and barriers.
- `calibrate_nm.R` — deterministic multi-start `optim(..., "Nelder-Mead")`
  with per-start logging.
- `laplace.R` — finite-difference Hessian, generalised inverse with
  eigenvalue cutoff, Monte-Carlo 95% intervals filtered by equilibrium
  feasibility.
- `mcmc.R` — priors, log-posterior, adaptive Metropolis sampler, posterior
  summaries and predictive intervals.
- `plot_results.R` / `plot_mcmc.R` — publication-quality ggplot2 figures.
- `scenarios.csv` — the single source of truth for sensitivity strategies.

### `scripts/` (entry points)

- `run_calibration.R` — sources `setup.R` + library modules, runs
  multi-start NM, Laplace intervals, fit metrics, writes `fit.rds` and
  figures under `output/calibration/<run_id>/`.
- `run_npe.R` — generates prior training data (R, parallel), trains the
  neural posterior estimator (`npe_train.py`), runs strict MCMC validation
  and posterior predictive checks, writes `output/calibration/npe_bayes/`.
- `run_analysis.R` — reads the MCMC posterior and `scenarios.csv`, runs
  every scenario x posterior-draw projection, writes
  `output/analysis/scenario_summary.csv`, key-year table and figures.
- `run_tests.R` — executes `tests/unit` and `tests/integration`.

## 3. Data flow

```text
setup.R ──> sim.cpp ──> likelihood.R (objective)
                           │
        calibrate_nm.R <───┘   multi-start NM  ──> fit.rds (best theta,
                                                       metrics, equilibrium)
                           │
        laplace.R ────────> Laplace 95% intervals (overlap check)
        run_npe.R ────────> NPE + MCMC posterior ──> posterior_samples_mcmc.csv
                           │
        run_analysis.R <───┘  scenarios.csv
                           └──> output/analysis/* (summary, key years,
                                                    figures)
```

## 4. Time conventions

- Model year `t = 0` corresponds to calendar 1970.
- Calibration runs `t = 0..150`; targets are read at `t = 45` (2015);
  equilibrium is checked at `t = 150` vs `t = 145`.
- Analysis uses the equilibrated state as the 2017 initial condition
  (`t = 47`) and projects to `t = 97` (2067).
