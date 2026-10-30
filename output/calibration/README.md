# Calibration output guide

All folders share one posterior target:
- **Parameters (13):** six log contact scales, six log entry scales (entry into D = scale × CNB
  2015 new drug abusers, scale 1–100) and log `c` (first arrest = `c` × re-arrest, 0.5–100).
- **Prior:** fixed independent Normal, truncated to the bounds: contact scales N(0, 1.5), entry
  scales N(0, 1) truncated at 0 (half-normal), log `c` N(0, 0.5) (`prior_spec()` in
  `scripts/calibration/likelihood.R`).
- **Likelihood (`nll_targets()`):** 2015 prison seroprevalence, Beta(p·φ, (1−p)·φ) with φ = 100;
  CNB 2015 new drug abusers against the model first arrests `c·λ3·D` by age, negative binomial
  (size 50); Teo et al. 2019 community PWID (11,000) against model D + F, log-normal (SD 0.245).
- **Restriction:** parameter bounds and the 10-year equilibrium check.
- **Solver:** daily RK4.

See `docs/training.md` for the commands.

| Folder | Role | Contents | Written by |
| --- | --- | --- | --- |
| `r_20261030_agepop/` | production | Nelder–Mead multi-start fit (posterior mode) | `scripts/run_calibration.R` |
| `mcmc_reference/` | **production posterior** | 8 chains × 40,000 retained iterations, stored every 10th (32,000 draws, `mcmc_draws.csv`); diagnostics, trace plots, adaptation and jump records, run config; 2015 stratum predictions (`model_prediction_2015_*`) and a predictive cache for the calibration figure | `scripts/run_mcmc_reference.R`, `scripts/predict_2015_posterior.R` |
| `smr_sensitivity/` | sensitivity only (old 7-parameter design, 2026-10-22; not rerun) | Nelder–Mead fits with the PWID mortality ratios at their CI bounds (predictions, solutions, config) | `scripts/run_calibration.R --smr-viraemic --smr-nonviraemic` |
| `mcmc_sbc/` | check only | Direct SBC of the reference MCMC: true prior draws (`theta_true.csv`, `x_true.csv`), synthetic targets, one chain per dataset, ranks and coverage (`mcmc_sbc_summary.csv`) | `scripts/sbc_mcmc.R` |
