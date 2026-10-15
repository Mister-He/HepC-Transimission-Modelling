# Output guide

Generated results; reruns overwrite these folders. Do not create new run folders unless versioned
parallel results are requested.

- `calibration/`: calibration results. See `calibration/README.md`.
  - `r_20261014_agepop/`: Nelder–Mead fit (prior centre).
  - `hr_sensitivity/`: Nelder–Mead fits at other viraemic mortality HRs (sensitivity only).
  - `mcmc_reference/`: reference MCMC posterior (production).
  - `npe_lhs/`: single-round NPE from 50,000 LHS prior simulations, with SBC (cross-validation only).
  - `npe_vs_mcmc/`: NPE versus reference MCMC comparison and the direct MCMC SBC (cross-validation only).
- `scenarios_2026/`: 2026–2030 policy scenarios and manuscript figures, run from the reference MCMC
  posterior. See `scenarios_2026/README.md`.
- `release_rate_fit.rds`, `rearrest_rate_fit.rds`: Stan fits of prison release and re-arrest
  (with the 2015 first-arrest flow `f_pop`), written by `scripts/model/preprocess.R`. `scripts/model/update_rate_params.R` copies their summaries into
  `scripts/model/setup.R`; the simulator itself does not read them.
