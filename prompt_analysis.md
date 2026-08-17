# prompt_analysis.md

Run this phase only after the NM and MCMC results pass all checks.

## Goal

Project the calibrated compartmental model 50 years forward from 2017 under
different intervention and epidemiological configurations, and report the
counts over time with Bayesian credible intervals propagated from the MCMC
posterior samples.

## Scenarios

Explore at least the following dimensions:

1. **Treatment rate by liver stage** (`tau = (NC, CC, DC, HCC)`):

   - `no treatment`: `(0, 0, 0, 0)` (reference);
   - `early treatment`: `(0.5, 0.5, 0, 0)`;
   - `broad treatment`: `(0.5, 0.5, 0.3, 0.2)`;
   - `aggressive treatment`: `(0.8, 0.8, 0.6, 0.4)`.

2. **GT3 proportion** (`rho`):

   - baseline `0.78`;
   - lower `0.6`;
   - higher `0.9`.

3. **Other configurations to consider**:

   - SVR efficacy for DC without RBV (`alpha_DC_neg`) instead of the
     RBV-containing regimen;
   - post-SVR progression reduction disabled (natural-history sensitivity);
   - a combined scenario (broad treatment + `rho = 0.6`).

## Projection setup

- First run the model to **at least 150 model years** with each posterior
  parameter draw; all states (including DC/HCC) must be converged.
- Take the equilibrated state as the initial condition for calendar year
  2017.
- Simulate 50 years forward from 2017 to 2067 under each scenario.
- Record annual counts of:

  ```text
  total HCV (acute + chronic + treatment) across all strata/stages
  decompensated cirrhosis (DC)
  hepatocellular carcinoma (HCC)
  ```

## Posterior propagation

- Draw a representative subsample of the MCMC posterior (e.g. 300-500
  equilibrium-feasible draws).
- For each draw and scenario, run the projection.
- Report median and 95% equal-tailed credible interval at key years:
  `2017, 2027, 2037, 2047, 2057, 2067`.
- The **status quo (no treatment)** scenario must always appear, with its
  own 95% credible interval.

## Figure requirements

- All credible intervals: 95%.
- Pure white background, no grid lines.
- Legend placed on the right.
- Scenario order in legends matches `analysis_report.md`:

  ```text
  no_treatment, early_treatment, broad_treatment, aggressive_treatment,
  rho_low, rho_high, dc_without_rbv, no_postsvr_modifiers,
  broad_treatment_rho60
  ```

## Outputs

```text
output/analysis/
  scenario_summary.csv     median/CrI by year and scenario
  fig_hcv_trajectories.png count vs year, total HCV
  fig_dc_hcc_trajectories.png DC/HCC count vs year
  analysis_report.md       concise interpretation
```

All figures must use `ggplot2`.
