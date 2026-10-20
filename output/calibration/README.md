# Calibration output guide

All folders share one posterior target:
- **Prior:** independent Normal(`prior_centre(best fit)`, 0.20) on the 12 log-scale parameters.
  `prior_centre()` moves each coordinate of the Nelder–Mead best fit to at least 5 prior SDs inside
  the parameter bounds (`scripts/calibration/likelihood.R`).
- **Likelihood:** Binomial prevalence and log-normal population (σ = 0.10).
- **Restriction:** parameter bounds and the 10-year equilibrium check.
- **Solver:** daily RK4.

See `npe_vs_mcmc/model_audit.md` for the full audit and `docs/training.md` for the commands.

| Folder | Role | Contents | Written by |
| --- | --- | --- | --- |
| `r_20261020_agepop/` | production | Nelder–Mead multi-start fit; `fit.rds` gives the prior centre | `scripts/run_calibration.R` |
| `mcmc_reference/` | **production posterior** | 4 chains × 40,000 retained iterations, stored every 10th (16,000 draws, `mcmc_draws.csv`); diagnostics, trace plots, adaptation and jump records, run config; 2015 stratum predictions (`model_prediction_2015_*`) and a predictive cache for the calibration figure | `scripts/run_mcmc_reference.R`, `scripts/predict_2015_posterior.R` |
| `npe_lhs/` | cross-validation only | Single-round NPE (5-network ensemble) trained on 50,000 Latin hypercube prior simulations: design and outcomes (`prior_design.csv`), counts (`simulation_summary.json`), training pairs, posterior draws, config, timings, SBC set (`theta_sbc.csv`, `x_sbc.csv`), SBC results (`sbc_summary.csv`, `sbc_ranks.npy`, `sbc_rank_histograms.png`), gate (`NPE_OK.txt`) | `scripts/run_npe.R`, `scripts/npe_train.py` |
| `hr_sensitivity/` | sensitivity only | Nelder–Mead fits at viraemic mortality HR 1, 1.11, 1.23 and 1.77 (predictions, solutions, config) | `scripts/run_calibration.R --hr-viraemic` |
| `npe_vs_mcmc/` | cross-validation only | NPE versus reference MCMC: parameter, joint and policy comparisons, figures, `report.md`, `model_audit.md`; `policy/npe_lhs/` holds the NPE policy draws (the reference side is `output/scenarios_2026/`); `mcmc_sbc/` holds the 8-dataset direct MCMC SBC | `scripts/compare_npe_mcmc.R`, `scripts/sbc_mcmc.R` |
