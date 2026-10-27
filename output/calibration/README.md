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

See `npe_vs_mcmc/model_audit.md` for the full audit and `docs/training.md` for the commands.

| Folder | Role | Contents | Written by |
| --- | --- | --- | --- |
| `r_20261027_agepop/` | production | Nelder–Mead multi-start fit (posterior mode) | `scripts/run_calibration.R` |
| `mcmc_reference/` | **production posterior** | 8 chains × 40,000 retained iterations, stored every 10th (32,000 draws, `mcmc_draws.csv`); diagnostics, trace plots, adaptation and jump records, run config; 2015 stratum predictions (`model_prediction_2015_*`) and a predictive cache for the calibration figure | `scripts/run_mcmc_reference.R`, `scripts/predict_2015_posterior.R` |
| `npe_lhs/` | cross-validation only | Single-round NPE (5-network ensemble) trained on 50,000 Latin hypercube prior simulations: design and outcomes (`prior_design.csv`), counts (`simulation_summary.json`), training pairs, posterior draws, config, timings, SBC set (`theta_sbc.csv`, `x_sbc.csv`), SBC results (`sbc_summary.csv`, `sbc_ranks.npy`, `sbc_rank_histograms.png`), gate (`NPE_OK.txt`, or `NPE_NOT_OK.txt` on failure; the 2026-10-27 run failed on rank shape) | `scripts/run_npe.R`, `scripts/npe_train.py` |
| `smr_sensitivity/` | sensitivity only (old 7-parameter design, 2026-10-22; not rerun) | Nelder–Mead fits with the PWID mortality ratios at their CI bounds (predictions, solutions, config) | `scripts/run_calibration.R --smr-viraemic --smr-nonviraemic` |
| `npe_vs_mcmc/` | cross-validation only | NPE versus reference MCMC: parameter, joint and policy comparisons, figures, `report.md`, `model_audit.md`; `policy/npe_lhs/` holds the NPE policy draws (the reference side is `output/scenarios_2026/`); `mcmc_sbc/` holds the 8-dataset direct MCMC SBC | `scripts/compare_npe_mcmc.R`, `scripts/sbc_mcmc.R` |
