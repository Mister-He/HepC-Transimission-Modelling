# AGENTS.md

## Scope

Calibrate the Singapore HCV PWID model on branch `dev_params_reasonable`.

## Baseline

- Original progression rates: `0.027 / 0.039 / 0.014 / 0.014`.
- Constant transmission (`m_min = m_max = 1`).
- 2015 Singapore baseline mortality and `omega = 14.68`.
- **No excess mortality**: `eta_s = 1`; do not fit or add `eta_s`.
- 12 fitted parameters: 6 contact scaling factors, 6 beta scaling factors.
- **Simulation horizon**: at least 150 model years; verify all-state
  convergence (including DC/HCC totals) before reporting equilibrium.
- **Beta scaling factors have no upper bound**.

## Soft preferences

- Scaled contact-matrix **diagonal** should be unimodal, peak near 30-39
  (or 20-29), with 50+/60+ diagonal entries small. Enforced as a small
  penalty in `likelihood.R`.
- Beta scales should be monotone non-decreasing and preferably `>1`.
- These preferences may be relaxed if they materially degrade fit.

## Editing rules

- `Model schematic.pptx`: never modify.
- `src/sim.cpp`: do not modify unless strictly necessary; currently unchanged.
- `src/setup.R`: original progression rates are required.
- Keep generated files minimal.
- All ggplot figures: pure white background, no grid lines; analysis legend
  on the right, scenario order per `analysis_report.md`, 95% CrIs.

## Calibration rules

- Multi-start Nelder-Mead; record starts, objectives, convergence, metrics.
- Acceptance criteria:

  ```text
  prevalence RMSE                    <= 0.02
  maximum absolute prevalence error <= 0.03
  population MAPE                    <= 0.10
  maximum population APE            <= 0.20
  equilibrium                       pass
  ```

- Laplace intervals must overlap observed intervals.

## Bayesian rules

- NPE/SNPE primary; strict MCMC validation required.
- Use a **wide thinning interval** (e.g. `thin = 20`).
- **Strict convergence**: R-hat in `[0.995, 1.005]` and pooled ESS `> 1000`
  for every parameter. Extend chains until both are met.
- Before reporting, numerically and visually review trace plots and density
  plots: chains must overlap, densities must be smooth and unimodal, and
  per-chain posterior means must be close.

## Sensitivity analysis

After a passing MCMC fit, execute `prompt_analysis.md`:

- treatment-rate configurations by liver stage;
- GT3 proportion configurations;
- other documented scenarios;
- 50-year projection from 2017;
- report total HCV, DC, and HCC counts over time;
- propagate MCMC posterior samples into credible intervals;
- report tables and figures.
