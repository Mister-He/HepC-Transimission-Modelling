# Output guide

Generated results; reruns overwrite these folders. Do not create new run folders unless versioned
parallel results are requested.

- `calibration/`: calibration results. See `calibration/README.md`.
  - `r_20261030_agepop/`: Nelder–Mead fit (posterior mode).
  - `smr_sensitivity/`: Nelder–Mead fits with the PWID mortality ratios at their CI bounds (sensitivity only).
  - `mcmc_reference/`: reference MCMC posterior (production).
  - `npe_lhs/`: single-round NPE from 50,000 LHS prior simulations, with SBC (cross-validation only).
  - `npe_vs_mcmc/`: NPE versus reference MCMC comparison and the direct MCMC SBC (cross-validation only).
- `scenarios_2026/`: 2026–2030 policy scenarios and manuscript figures, run from the reference MCMC
  posterior. See `scenarios_2026/README.md`.
- `scenarios_2026_jcov/`: the same scenarios at prison (J) treatment coverage 0%, 50% and 90% (main
  analysis has 99.995%), with the coverage heatmap.
- `scenarios_2026_j2019/`: sensitivity analysis with prison treatment from 2019.
- `release_rate_fit.rds`, `rearrest_rate_fit.rds`: Stan fits of prison release and re-arrest
  (survival models only), written by `scripts/model/preprocess.R`. `scripts/model/update_rate_params.R` copies their summaries into
  `scripts/model/setup.R`; the simulator itself does not read them.
