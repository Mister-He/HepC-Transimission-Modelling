# prompt_mcmc.md

Run the Bayesian phase only after the NM fit passes acceptance.

## Method

- NPE/SNPE primary; strict MCMC validation required.
- Priors in log scaling-factor space, centred on the NM anchors.
- At least three chains from different starts.
- **Wide thinning**: `thin = 20` (or larger).
- **Strict convergence**:

  ```text
  R-hat in [0.995, 1.005]   (every parameter)
  pooled ESS > 1000         (every parameter)
  ```

- Extend warm-restarted chains until both criteria pass.

## Review before reporting

Inspect:

- trace plots: chains mix, no long-term drift, three chains overlap;
- density plots: smooth, unimodal, chains overlap;
- per-chain posterior means and quantiles are close.

Do not report an MCMC result that fails visual or numerical review.

## Outputs

```text
output/calibration/npe_bayes/
  posterior_samples_mcmc.csv
  diagnostics_mcmc.csv
  credible_intervals.csv
  fig_trace.png
  fig_density_mcmc.png
  fig_ci_compare.png
```
