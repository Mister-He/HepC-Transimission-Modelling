# Calibration output guide

All folders share one posterior target:
- **Parameters:** six log contact scales and log `c` (first arrest = `c` × re-arrest), 7 in all.
- **Prior:** fixed independent Normal, truncated to the bounds: contact scales N(0, 1.5), log `c`
  N(0, 0.5), with `c` in 0.5–100 (`prior_spec()` in `scripts/calibration/likelihood.R`).
- **Likelihood:** seroprevalence only: Beta(p·φ, (1−p)·φ) per age group with φ = 100. No counts
  and no population target.
- **Restriction:** parameter bounds and the 10-year equilibrium check.
- **Solver:** daily RK4.

See `npe_vs_mcmc/model_audit.md` for the full audit and `docs/training.md` for the commands.

| Folder | Role | Contents | Written by |
| --- | --- | --- | --- |
| `r_20261021_agepop/` | production | Nelder–Mead multi-start fit (posterior mode) | `scripts/run_calibration.R` |
| `mcmc_reference/` | **production posterior** | 4 chains × 40,000 retained iterations, stored every 10th (16,000 draws, `mcmc_draws.csv`); diagnostics, trace plots, adaptation and jump records, run config; 2015 stratum predictions (`model_prediction_2015_*`) and a predictive cache for the calibration figure | `scripts/run_mcmc_reference.R`, `scripts/predict_2015_posterior.R` |
| `npe_lhs/` | cross-validation only | Single-round NPE (5-network ensemble) trained on 50,000 Latin hypercube prior simulations: design and outcomes (`prior_design.csv`), counts (`simulation_summary.json`), training pairs, posterior draws, config, timings, SBC set (`theta_sbc.csv`, `x_sbc.csv`), SBC results (`sbc_summary.csv`, `sbc_ranks.npy`, `sbc_rank_histograms.png`), gate (`NPE_OK.txt`) | `scripts/run_npe.R`, `scripts/npe_train.py` |
| `smr_sensitivity/` | sensitivity only | Nelder–Mead fits with the PWID mortality ratios at their CI bounds (predictions, solutions, config) | `scripts/run_calibration.R --smr-viraemic --smr-nonviraemic` |
| `npe_vs_mcmc/` | cross-validation only | NPE versus reference MCMC: parameter, joint and policy comparisons, figures, `report.md`, `model_audit.md`; `policy/npe_lhs/` holds the NPE policy draws (the reference side is `output/scenarios_2026/`); `mcmc_sbc/` holds the 8-dataset direct MCMC SBC | `scripts/compare_npe_mcmc.R`, `scripts/sbc_mcmc.R` |
