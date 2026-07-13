# Two-Step Calibration Workflow Changes

This repository now supports a two-step calibration workflow:

1. Run Nelder-Mead MAP calibration to obtain reasonable log-scale parameter estimates.
2. Use the spline-smoothed Nelder-Mead estimates as the prior center for HMC posterior sampling.

Run the full workflow with:

```bash
./calib.sh
```

All generated outputs are written under `two-steps-calibration/`.

## Files Modified

- `HMC_core.r`
  - Replaced independent age-specific Normal priors with a shared 5-basis natural-spline Gaussian prior for both log-parameter families.
  - Added helpers for constructing, validating, smoothing, and evaluating spline-basis priors.
  - Kept positivity constraints through the existing `exp(theta)` transform and added finite log-scale bounds to reject numerically invalid proposals.
  - Added prevalence-curve shape regularization alongside the existing population-composition shape regularization.

- `nelder_mead.R`
  - Uses the shared 5-basis spline prior during MAP calibration.
  - Saves both raw and spline-smoothed Nelder-Mead estimates.
  - Accepts an output directory argument and writes fit objects and plots into that directory.

- `hmc.R`
  - Accepts an output directory and Nelder-Mead fit path.
  - Constructs the HMC spline prior from the spline-smoothed Nelder-Mead estimates.
  - Initializes chains near the smoothed Nelder-Mead prior center.
  - Writes HMC outputs and plots into the requested output directory.

- `calib.sh`
  - Runs the complete two-step workflow and stores all outputs in `two-steps-calibration/`.

## Smoothing Choices

Both calibration steps use spline-basis priors with `n_knots = 5`. The prior covariance is built from the spline basis plus a small residual variance, so age-specific parameters can deviate from a smooth spline curve but jagged profiles are penalized. The HMC prior is centered on a spline projection of the Nelder-Mead estimate rather than the raw MAP vector to avoid carrying over overfitted age-to-age noise.
