# Evaluation

Evaluate the model in four layers:

1. age-specific seroprevalence residuals;
2. age-specific prison-population residuals;
3. all-compartment equilibrium stability between model years 145 and 150; and
4. NPE simulation-based calibration, posterior spread, and predictive coverage.

The 2015 snapshot must pass the equilibrium gate before projection. The primary
scenario screen is the WHO HCV impact target: incidence at or below 20% of the
2015 baseline and HCV-related mortality at or below 35% of the 2015 baseline
in 2030. Ever-PWID viremic prevalence at or below 1% is retained only as a
secondary diagnostic and is not evidence of national elimination.

Run-specific results are stored under `output/calibration/` and
`output/scenarios_2026/`. NPE results are primary only when `NPE_OK.txt` is
present.

## Current run

The corrected 2015 age-specific prevalence targets and population observation
sigma (`0.10`) are in `data/calibration_targets.json`. The deterministic
multi-start calibration in `output/calibration/r_20260929_agepop/` achieved an
objective of `19.56661`, prevalence RMSE `0.00872`, population MAPE `0.037718`,
and a passing equilibrium gate. See `predictions.csv`, `solutions.csv`, and
`run_config.json` there.

The primary posterior is `output/calibration/npe_agepop_recal/`. It uses one
shared set of 3,000 training simulations, one-round MAF NPE-C, two independently
trained seeds, 5,000 posterior draws per seed, 128 hidden features, 10
transforms, 500 maximum epochs, and a
best-fit-centred log-scale Normal prior with SD `0.20` for both parameter
blocks. Prevalence observation noise is Binomial and population noise is
log-normal with sigma `0.10`. Both seeds pass the 95% and 90% SBC coverage
checks and the minimum posterior SD `0.05` check. The authoritative gate is
`NPE_OK.txt`; the seed-specific results are in `npe_evaluation*.csv` and
`sbc_summary*.csv`.

`output/scenarios_2026/` now contains a J-only comparator, 15 finite annual
F/X screening-quota policies and an unconstrained-screening stress test. All
policies share J treatment from 2019 and the same 2026 branch state. The
WHO-style model-wide target is the joint 80% incidence and 65% HCV-related
mortality reduction versus 2015; F/X is also evaluated separately. Three of
15 finite policies meet both components at the NPE median: annual screening
5,000 with F 25%, and 8,000 with F 25% or 50%. The 5,000/F25 mortality
margin is only about 0.13 percentage points. The ideal screening case meets
both within this model, requiring about 80,908
F/X screenings and 1,877 F/X treatment starts over 2026--2030. It is a
feasibility stress test rather than an observed program. Cost accounting now
includes only 2026--2030 F/J/X services, discounted 3% per year to 2026.
The SOF/VEL model price is USD 2,577/week for 12 weeks, converted with the
2025 average SGD/USD rate; this and the service charges are not verified
Singapore payer net costs.
See `output/scenarios_2026/README.md` for methods, evidence and limitations.

Incidence is the instantaneous annualized flow of new and repeat infections
from `u` and `r` into `a` in D/F; mortality is the DC/HCC HCV-attributed
annual death flow. These are event counts per year, not per-person rates.
Treatment starts, completions and cures are integrated model events, counted
as courses and not distinct people. Re-treatment consumes a new start quota.
The calibration was not repeated after changing the 2019-onward historical J
assumption, so 2019 observed treatment counts are a separate validation check.

## Posterior and predictive checks

The 2015 stratum-level posterior prediction in
`output/calibration/npe_agepop_recal/` uses all 10,000 combined NPE draws;
every draw was finite and passed the equilibrium check. Its median model
populations are approximately 2,900 current PWID (D+F) and 20,915 ever-PWID
(D+J+F+X). The NPE/MCMC age calibration plot is in `figures/`. These checks
concern the fitted PWID model; they are not estimates for Singapore's general
population.

`figures/calibration_interval_overlap.csv` shows that all six age-specific
seroprevalence and all six prison-population 95% predictive intervals from
both NPE and MCMC overlap their respective target intervals. In ages 30--39
and 40--49, the NPE prevalence interval widths are about 4.45 and 4.83
percentage points, respectively, below the target observation interval widths
of about 6.32 and 7.12 points. Overlap is a calibration diagnostic, not
validation of prospective treatment-policy predictions.

MCMC is initialized using the combined NPE covariance and targets the same
likelihood over the full 150-year equilibrium horizon with a monthly
integration step. Its run settings are in
`output/calibration/mcmc_recal/mcmc_run_config.json`. The strict gate requires
every parameter to have ESS above 4,000 and R-hat below 1.001; the authoritative
result is `MCMC_STRICT_GATE.txt`. This four-chain run passed: the minimum pooled
ESS is 10,302 and the maximum R-hat is 1.000552. MCMC provides an independent
posterior sensitivity check; its parameter medians differ from the combined
NPE medians by at most 0.0195 on the log scale. The 11-policy screen uses the quality-gated NPE posterior median.
