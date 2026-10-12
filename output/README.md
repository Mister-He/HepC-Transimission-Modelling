# Output guide

Generated results; reruns overwrite these folders. Do not create new run folders unless versioned
parallel results are requested.

- `calibration/`: calibration results. See `calibration/README.md`.
  - `r_20261010_agepop/`: Nelder–Mead fit (prior centre).
  - `mcmc_reference/`: reference MCMC posterior (production).
  - `npe_lhs/`: single-round NPE from 50,000 LHS prior simulations.
  - `npe_vs_mcmc/`: NPE versus reference MCMC comparison.
- `scenarios_2026/`: 2026–2030 policy scenarios and manuscript figures, run from the reference MCMC
  posterior. See `scenarios_2026/README.md`.
- `release_rate_fit.rds`, `rearrest_rate_fit.rds`: Stan fits of prison release and re-arrest written by
  `scripts/model/preprocess.R`. The simulator does not read them; they are kept for provenance.
