# Two-step calibration changes

This repository now supports a two-step calibration workflow:

1. `nelder_mead.R` estimates MAP spline coefficients for the two age-varying
   input vectors.
2. `hmc.R` reads `two-steps-calibration/nelder_mead_fit.rds`, centers the HMC
   spline-coefficient priors on the Nelder-Mead estimates, and samples from the
   posterior.

Run the full workflow with:

```sh
./calib.sh
```

Optional runtime controls:

```sh
CHAINS=2 CORES=2 WARMUP=100 ITER=150 PPC=100 ./calib.sh
```

All generated RDS files, plots, and logs are written under
`two-steps-calibration/`.

## Files modified

- `HMC_core.r`
  - Keeps the binomial prevalence likelihood:
    `obs_pos[a] ~ Binomial(obs_tot[a], q_age[a])`.
  - Keeps the input transform as a B-spline basis with 5 internal knots for
    both `C_contact_scale` and `tot_in_scaling_fct`.
  - Adds `configure_spline_priors()` so HMC can use Nelder-Mead coefficients as
    valid prior centers while preserving positivity through the existing
    exponential transform.
  - Strengthens spline smoothness with `SPLINE_RW_SD = 0.22` and loosens the
    first two coefficient anchors to `SPLINE_ANCHOR_SD = 0.75`.

- `nelder_mead.R`
  - Uses the same observation vectors and binomial prevalence target as `hmc.R`.
  - Uses shared population/shape regularization settings:
    `sigma_pop = c(0.20, rep(0.12, 9))`, `sigma_shape = 0.20`,
    `nu_shape = 7`.
  - Adds deterministic early-age starts to reduce the chance of the local
    solution with a near-zero first-age prevalence.
  - Saves all outputs to `two-steps-calibration/` by default.

- `hmc.R`
  - Refactored into `run_hmc_calibration()` plus a CLI entry point.
  - Reads the Nelder-Mead fit via `--nm-fit=...`, converts `theta_hat` into
    HMC spline-prior centers, and initializes chains near the MAP.
  - Saves HMC outputs and plots under `two-steps-calibration/`.

- `HMC_conv.R`
  - Fixes convergence diagnostics so log-scale diagnostics report the sampled
    spline coefficients only.
  - Adds `plot_posterior_fitted_curves()`, which plots smooth posterior fitted
    prevalence and population-count curves from model means instead of noisy
    posterior predictive replicates.

- `calib.sh`
  - Runs the complete two-step workflow and writes step logs to
    `two-steps-calibration/nelder_mead.log` and
    `two-steps-calibration/hmc.log`.
