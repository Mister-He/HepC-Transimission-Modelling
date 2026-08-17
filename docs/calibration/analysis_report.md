# Sensitivity analysis — 50-year projections from 2017

## Method

- Uses the passing MCMC posterior (`output/calibration/npe_bayes/`).
- 300 equilibrium-feasible posterior draws per scenario.
- Each draw is run 150 model years to full equilibrium, then the
  equilibrated state is used as the 2017 initial condition and projected to
  2067.
- Reported quantities: total HCV (acute + chronic + treatment), DC, HCC.
- Median and 95% equal-tailed credible intervals at
  `2017, 2027, 2037, 2047, 2057, 2067`.

## Scenarios

```text
no_treatment               tau = (0, 0, 0, 0), rho = 0.78
early_treatment            tau = (0.5, 0.5, 0, 0)
broad_treatment            tau = (0.5, 0.5, 0.3, 0.2)
aggressive_treatment       tau = (0.8, 0.8, 0.6, 0.4)
rho_low                    tau = 0, rho = 0.60
rho_high                   tau = 0, rho = 0.90
dc_without_rbv             DC SVR without RBV
no_postsvr_modifiers       post-SVR progression reduction disabled
broad_treatment_rho60      broad treatment + rho = 0.60
```

## Key findings

- Baseline 2017 total HCV is about 8,850 (95% CrI 8,173-9,580).
- Status quo is flat over 2017-2067: DC about 394 (CrI 356-440) and HCC
  about 103.5 (CrI 93.4-115.5), confirming full equilibrium.
- Without treatment, HCV remains high over 50 years; aggressive treatment
  reduces total HCV to near zero within 10 years.
- DC and HCC are stable under no treatment; treatment accelerates their
  decline.
- GT3 proportion has modest effects; rho = 0.90 yields slightly faster
  progression to HCC.
- Disabling post-SVR progression modifiers increases long-run DC/HCC counts
  relative to broad treatment.

## Files

```text
output/analysis/scenario_summary.csv
output/analysis/scenario_key_years.csv
output/analysis/fig_hcv_trajectories.png
output/analysis/fig_dc_hcc_trajectories.png
```
