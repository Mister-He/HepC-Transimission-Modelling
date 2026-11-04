# Output guide

Generated results; reruns overwrite these folders. Do not create new run folders unless versioned
parallel results are requested.

- `calibration/`: calibration results. See `calibration/README.md`.
  - `r_20261103_agepop/`: Nelder–Mead fit (posterior mode) after the 2026-11-03 mortality fix (one route per
    cause of death; 60+ background rate for ages 60–74). Current.
  - `r_20261030_agepop/`: previous fit (old mortality). Superseded; kept until the authors remove it.
  - `smr_sensitivity/`: Nelder–Mead fits with the PWID mortality ratios at their CI bounds, from the
    2026-10-22 7-parameter design (historical; not rerun).
  - `mcmc_reference/`: reference MCMC posterior (production).
  - `mcmc_sbc/`: direct MCMC simulation-based calibration on 8 synthetic datasets (check only). NPE was
    removed on 2026-10-30.
- `scenarios_2026/`: 2026–2030 policy scenarios and manuscript figures, run from the reference MCMC
  posterior. See `scenarios_2026/README.md`.
- `scenarios_2026_jcov/`: the same scenarios at prison (J) treatment coverage 0%, 50% and 90% (main
  analysis has 99.995%), with the coverage heatmap.
- `scenarios_2026_j2019/`: sensitivity analysis with prison treatment from 2019.
- `scenarios_2026_inflow/`: sensitivity analysis with the observed 2016–2025 CNB inflow series
  (`--inflow observed`; `REVIEW_REPORT.md` R-02).
- `scenarios_2026_gt3/`: sensitivity analysis with a genotype 3 share of 10% (`rho010/`) and 50% (`rho050/`)
  instead of 78% (`--param-file`; posterior not refitted).
- `tornado/`: one-way sensitivity of the 2030 results to each fixed input at the posterior mode
  (`scripts/tornado.R`). See `tornado/README.md`.
- `prison_2015_counts.csv`: aggregate 2015 prison counts by age (`n_held`: everyone held for drug
  consumption on 31 Dec 2015; `n_in_prison`, `positive`, `screened`: the screened subset), written by
  `scripts/model/preprocess.R`.
- `release_rate_fit.rds`, `rearrest_rate_fit.rds`: Stan fits of prison release and re-arrest
  (survival models only), written by `scripts/model/preprocess.R`. `scripts/model/update_rate_params.R` copies their summaries into
  `scripts/model/setup.R`; the simulator itself does not read them.
