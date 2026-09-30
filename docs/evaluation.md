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

The deterministic age/population run is
`output/calibration/r_20260929_agepop/`. Its selected objective is `33.1808`,
with prevalence RMSE `0.0649`, population MAPE `0.0723`, and a passing
equilibrium gate. The multi-start table and age-specific residuals are in
`solutions.csv` and `predictions.csv`.

The accepted NPE result is in `output/calibration/npe_agepop/`. It used 600
simulations, one-round NPE-C, two seeds, 1,000 posterior draws, and 250 maximum
epochs. Both seeds passed the 95%/90% SBC coverage gate and the posterior-spread
gate. Two-round sensitivity runs are retained separately because each had one
parameter just below the coverage threshold.

The 2026 scenario run is in `output/scenarios_2026/`. The scenario script uses
the calibrated `fit$best$theta` to construct the 2015 snapshot and evaluates
75 D-free combinations. Under the WHO screen, 6 of 75 pass; all six are in
the all-age setting. No 20--50 or 30--60-only combination passes both targets.
The endpoint and treatment-demand results are in `scenario_2030_summary.csv`
and `target_scenario_treatment_summary.csv`. The primary trade-off output is
`scenario_tradeoff.csv`, with target-reaching policies in
`target_scenario_tradeoff.csv` and the combined figure in
`scenario_tradeoff.png`.

The 2015 stratum-level posterior predictions are in
`output/calibration/npe_agepop/model_prediction_2015_strata.csv`, with a
report-ready table in `model_prediction_2015_strata_pretty.csv` and a figure
in `model_prediction_2015_strata.png`. They use the combined accepted NPE
draws from both seeds; 1,991 of 2,000 finite draws pass the equilibrium gate.
Total population includes `u+a+c+t+r`, while active HCV cases are `a+c+t`.
These are conditional model diagnostics, not national elimination estimates.
The CSV also retains
`current_viremic_prev` for the D+F subgroup.

Incidence is the model flow of new infections into acute infection from `u`
and reinfected `r` in D/F. HCV-related mortality is the model-specific
disease mortality component from DC/HCC; background PWID mortality is not
counted as HCV-related. The reported 2015 and 2030 rates are reference-year
annualized rates evaluated at the corresponding model snapshots.

Treatment starts are numerically integrated from `tau * c` for each stratum
over 2026--2030. `treatment_courses_per_2026_person_*` is a secondary
descriptive ratio. It is cumulative courses divided by the population present
in 2026, not the fraction of unique people treated and may include repeat
courses after reinfection. The primary efficiency metric is treatment starts
per 1,000 cumulative viremic person-years averted.

The scenario code also contains an explicit regression fix for treatment-rate
indexing: C++ reads `tau_by_stratum` in stratum-major order, so the R matrix is
transposed before flattening. Earlier scenario files produced before this fix
must not be used for interpretation.

## Current NPE and MCMC sensitivity

The current NPE result is in `output/calibration/npe_agepop_recal/`. It uses
3,000 simulations per seed, one-round MAF NPE-C, 128 hidden features, 10
transforms, 500 maximum epochs, and a best-fit-centred log-scale Normal prior
with SD 0.15. Prevalence observation noise is Binomial and population noise is
log-normal with sigma 0.10. Both seeds pass the SBC 95%/90% coverage gate.

The predictive interval comparison is in
`figures/calibration_interval_overlap.csv`; the figure is
`figures/calibration_agepop.png`. The narrower prior reduces the seroprevalence
BCI width, but the 30--49 age groups remain wider than their target intervals,
indicating residual parameter non-identifiability. The 60+ posterior remains
centred above its target, so overlap is not claimed for that age group.

MCMC is initialized from the combined NPE covariance and uses the same target
likelihood. The current validation uses the full 150-year horizon with an
explicit monthly integration step recorded in
`output/calibration/mcmc_recal/mcmc_run_config.json`. It is accepted only when
every parameter has ESS > 4,000 and R-hat < 1.0001; the result of that gate is
in `output/calibration/mcmc_recal/MCMC_STRICT_GATE.txt`. The current gate is
FAIL: ESS passes, but several R-hat values remain above 1.0001. Therefore NPE
remains the usable posterior and MCMC remains a diagnostic sensitivity check.

The NPE-median scenario rerun overwrote `output/scenarios_2026_recal_npe/`:
75 combinations were evaluated and 6 all-age combinations passed both WHO
impact thresholds; no 20--50 or 30--60-only combination passed.
