# Calibration output guide

All folders share one posterior target:
- **Prior:** independent Normal(Nelder–Mead best fit, 0.20) on the 12 log-scale parameters.
- **Likelihood:** Binomial prevalence and log-normal population.
- **Restriction:** the 10-year equilibrium check.
- **Solver:** daily RK4.

See `npe_vs_mcmc/model_audit.md` for the full audit.

| Folder | Contents | Written by |
| --- | --- | --- |
| `r_20261010_agepop/` | Nelder–Mead multi-start fit; `fit.rds` is the prior centre for NPE and MCMC | `scripts/run_calibration.R` |
| `mcmc_reference/` | **Production posterior**: 4 chains × 45,000 retained draws (`mcmc_draws.csv`), diagnostics, run config, 2015 stratum predictions (`model_prediction_2015_*`) and a predictive cache for the calibration figure | `scripts/run_mcmc_reference.R`, `scripts/predict_2015_posterior.R` |
| `npe_lhs/` | Single-round NPE trained on 50,000 Latin hypercube prior simulations: design and outcomes (`prior_design.csv`), counts (`simulation_summary.json`), training pairs, posterior draws, config, timings, gate (`NPE_OK.txt`) | `scripts/run_npe.R`, `scripts/npe_train.py` |
| `npe_vs_mcmc/` | Comparison of NPE with the reference MCMC: parameter, joint and policy comparisons, figures, `report.md`, `model_audit.md`; `policy/` holds the 1,000-state scenario runs (reference and NPE-LHS); `previous_random_npe/` holds the earlier random-sampling single-round NPE replicates (posterior draws and policy draws) for context | `scripts/compare_npe_mcmc.R` |

The MCMC draws file has two columns from its original benchmark run, `forward_budget` (all TRUE) and
`time_matched` (all NA). Downstream scripts read only `theta1..theta12`.
