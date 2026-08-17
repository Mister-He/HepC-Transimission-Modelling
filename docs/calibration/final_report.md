# Final report

## Run

```text
output/calibration/run5_redo150/
```

## Baseline

- Original progression rates: `0.027 / 0.039 / 0.014 / 0.014`.
- No excess mortality (`eta_s = 1`).
- 2015 Singapore mortality, `omega = 14.68`.
- Simulation horizon: `t = -10 .. 140` (150 model years); all states,
  including DC/HCC totals, are converged at the equilibrium gate.
- Beta scaling factors: no upper bound.

## Fitted scales

```text
contact scales: 5.558, 0.0811, 0.0517, 0.4570, 7.120, 0.539
beta scales:    0.1688, 0.9140, 1.2559, 5.4211, 14.090, 118.15
```

## Fit metrics

```text
best NLL                          22.80
prevalence RMSE                   0.0088
maximum prevalence error          0.0209
population MAPE                   0.0493
maximum population APE            0.1414
equilibrium (all states)          pass
```

Laplace intervals overlap all observed intervals.

## Bayesian

```text
output/calibration/npe_bayes/
```

- NPE: 30,000 simulations, 2 rounds, 2 seeds.
- MCMC validation: 3 chains, 30,000 total iterations, burn-in 5,000,
  thin 20.
- R-hat range: `1.0002 - 1.0035`.
- Pooled ESS range: `1432.0 - 2026.2`.

## Sensitivity analysis

```text
output/analysis/
docs/calibration/analysis_report.md
```

Status quo (no treatment) is now flat over 2017-2067:
DC remains about 394 and HCC about 103.5, confirming equilibrium.
