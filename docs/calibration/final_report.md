# Final report

## Run

```text
output/calibration/run1_4strata/
```

## Model

- 4 strata: `D/J/F/X` (`lambda1` = first arrest, `lambda3` = re-arrest).
- 4 liver stages: NC, CC, DC, HCC.
- 4 HCV states: `u/a/c/t` (no separate cleared state; `u` includes
  susceptible and post-SVR).
- 6 age groups.
- Simulation: `t_start = 0`, `t_end = 150`.
- Original progression rates; no excess mortality; beta scales unbounded.

## Fit metrics

```text
best NLL                          21.42
prevalence RMSE                   0.0051
maximum prevalence error          0.0123
population MAPE                   0.0360
maximum population APE            0.0935
equilibrium (all states)          pass
```

Laplace intervals overlap all observed intervals.

## Bayesian

```text
output/calibration/npe_bayes/
```

- thin 20; R-hat range `0.9997 - 1.0041`; pooled ESS `1309.6 - 1958.2`.
- Density plots are reported on the log-parameter scale; skewness is
  moderate (max absolute skew about 0.9) and marginals are close to normal.

## Analysis

```text
output/analysis/
```

Status quo is flat over 2017-2067: DC about 479 (CrI 429-535) and HCC about
125.6 (CrI 112.7-140.1).
