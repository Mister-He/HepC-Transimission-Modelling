# Problems and follow-up

## 2026-09-29: age/population calibration and 2026 scenarios

- The 47% halfway-house/release-cohort observation is no longer a calibration
  target. The active likelihood contains only six age-specific
  seroprevalence targets and six age-specific population counts from
  `data/calibration_targets.json`; the age targets are not multiplied by
  `0.74`.
- The deterministic fit used 12 log-scale parameters, five recorded
  multi-starts, `maxit = 1000`, and a 150-year simulation for every objective
  evaluation. All recorded starts passed the equilibrium gate. The selected
  start has objective `33.1808`, prevalence RMSE `0.0649`, and population
  MAPE `0.0723`. The fit is usable as the NPE simulator centre, but residuals
  remain notable for `<20` and `60+`; this is a limitation of the current
  12-parameter structure, not evidence of exact identifiability.
- After correcting the post-2015 calendar-row index, a provenance rerun was
  started with the same full `maxit = 1000` setting. It was stopped during a
  later start because the corrected index is inactive during the 0--150
  calibration horizon; no partial output was used. The completed five-start
  fit above remains the calibration result.
- The model now uses the official SingStat/CNB age-specific new-drug-abuser
  series for 2016--2025, holds 2025 constant after the latest observation,
  and applies the same calibrated age-specific beta scales through 2030.
- An implementation issue was found in the first scenario run: it used
  uncalibrated `fit$base_params` instead of applying `fit$best$theta`. The
  scenario script now constructs the snapshot from calibrated parameters and
  all figures/tables were regenerated. Under the corrected model, status quo
  current-PWID viremic prevalence is about `0.284` in 2015 and `0.181` in
  2030.
- The plots and endpoint target were subsequently changed to Ever-PWID
  (D+J+F+X) viremic prevalence. A second scenario implementation issue was
  then found: R flattened the stratum-by-stage treatment matrix in a different
  order from the C++ index. This made the first scenario outputs invalid; the
  `tau_by_stratum` construction now uses stratum-major order and the outputs
  were regenerated. In that previous low-intensity scenario set, the best tested prison/community
  intensity-grid point is about `0.0050` at `(prison_rate, community_rate) =
  (2, 2)`. The current J/F/X-only `mixed_all_ages` scenario is about `0.0115`,
  just above the <=1% operational screen. A counterfactual that also treats D
  at 1.0/y is about `0.00455`. The CSV retains current-PWID prevalence as a
  secondary diagnostic. Those values belong to the previous scenario set and
  are not part of the current D-free policy comparison. These are conditional
  model results, not a national elimination conclusion.
- The scenarios were redesigned around treatment initiation probability, with D
  treatment prohibited in every policy. The current run has 75 combinations:
  J and F/X at 0%, 25%, 50%, 75%, and 100%, crossed with age settings 20--50,
  30--60, and all ages. Six all-age combinations meet both WHO impact targets;
  no age-restricted combination does. The previous 1% Ever-PWID prevalence
  screen is now secondary.
  Treatment-demand tables integrate `tau*c` and report treatment
  courses, not unique treated people. Unique-person treatment counts cannot be
  recovered from the aggregate compartment states without adding individual
  treatment-history compartments.
- Average courses per 2026 person is retained only as a secondary descriptive
  ratio. The primary trade-off uses treatment starts versus cumulative viremic
  person-years averted relative to status quo, plus a Pareto-efficient flag and
  treatment starts per 1,000 viremic person-years averted.
- The WHO mortality target uses the model's DC/HCC disease-specific mortality
  component, excluding background PWID mortality. Incidence and mortality are
  reference-year annualized model rates at the 2015 and 2030 snapshots, not
  national surveillance estimates.
- During implementation of the coverage parameterization, an intermediate
  draft accidentally multiplied the converted J hazard twice through the
  stage and stratum factors. This was caught by inspecting the generated
  endpoints before accepting the run; the final implementation uses a neutral
  stage multiplier and one converted hazard per stratum.

## NPE status

- NPE-C was rerun with `/tmp/bayes-venv/bin/python` (`torch 2.14.0`, `sbi
  0.27.0`), JSON targets, binomial prevalence noise, log-normal population
  noise, two seeds, and a best-fit-centred prior. No Laplace approximation
  was used.
- A 300-simulation, two-round run failed narrowly: one seed had `theta8`
  95% SBC coverage `0.83`. A 600-simulation, two-round run also failed
  narrowly for `theta2` in the second seed (`0.84`). These diagnostics are
  retained under `output/calibration/npe_agepop_two_round_300/` and
  `output/calibration/npe_agepop_two_round_600/`.
- The accepted result is a 600-simulation, two-seed, one-round NPE run in
  `output/calibration/npe_agepop/`. Both seeds pass the 95%/90% SBC coverage
  ranges and all posterior standard deviations exceed the collapse threshold.
  One-round training was selected because it retains prior support and was
  more stable than the small-sample sequential proposal.
- `tests/unit/test_npe_coverage.R` now executes rather than skips and passes.
  The accepted posterior is still based on a moderate simulation budget;
  before publication, increase simulations and epochs and re-run the same
  gate as a sensitivity check.

## Structural limitations to review

- The prison age/population likelihood alone does not identify all contact
  and inflow scales strongly. The multi-start range and age-specific residuals
  should be reported with any posterior result.
- The WHO impact screen is applied only to the modelled PWID network. It cannot
  support a whole-Singapore claim because non-PWID populations and national
  incidence and mortality are not modelled. The 1% prevalence result remains a
  secondary operational diagnostic.
- The SingStat new-abuser series is an observed enforcement-based inflow
  proxy, not a complete census of new PWID. Scenario conclusions are
  conditional on that interpretation and on holding the latest observed
  inflow constant through 2030.

## Historical attempts

- Earlier runs included the 47% release-cohort target and produced a
  persistent model-data conflict. That target was removed as requested rather
  than resolved by artificially widening uncertainty.
- Earlier Laplace and release-cohort outputs are no longer active inputs and
  were removed from the current calibration workflow.

## 2026-09-29: widened NPE and NPE-initialized MCMC validation

- A deterministic `maxit = 3000` multi-start recalibration was attempted in
  `output/calibration/r_20260929_agepop_recal/`, but the first start was still
  running after an extended wall-clock interval and was stopped before any
  fit, solutions, or predictions were written. The previous accepted
  `r_20260929_agepop` fit was therefore retained as the simulator centre.
- A wider one-round NPE was then trained in
  `output/calibration/npe_agepop_recal/` using the same target observation
  assumptions: Binomial age-specific prevalence noise and log-normal
  population noise with sigma 0.10. The log-scale contact and beta prior SD
  was increased to 1.2, with 1,200 training simulations, 3,000 posterior
  draws per seed, two seeds, and 300 maximum epochs. Both seeds passed the
  existing SBC/spread gate.
- The NPE age-specific predictive intervals were recomputed from 6,000 draws
  and 1,991 equilibrium-pass draws. All six prevalence intervals overlap the
  corresponding 95% Binomial target interval, including 60+; all six
  population intervals overlap the corresponding 95% log-normal target
  interval. The overlap table is
  `figures/calibration_interval_overlap.csv`.
- Adaptive Metropolis was initialized with the empirical covariance of the
  combined NPE seeds and the same best-fit-centred Normal(1.2) log-scale prior.
  Three chains were run for 1,000 iterations each, with 200 burn-in and thin
  2. The mean acceptance rate was 0.270, but R-hat ranged from 1.06 to 3.03
  and pooled ESS ranged from 62 to 104. The chains are not converged; this
  MCMC run cannot validate NPE. Its 60+ seroprevalence predictive interval
  also does not overlap the 95% target interval. A longer exact MCMC run is
  still computationally expensive because every proposal integrates the full
  150-year model.
- The calibration figure now shows Target, NPE posterior, and MCMC posterior
  medians with their respective 95% target/credible intervals. MCMC is shown
  for transparency but must not be treated as an accepted posterior until
  convergence is improved.
- Scenarios were rerun to `output/scenarios_2026_recal/` using the retained
  deterministic fit. There are 75 policy combinations and 6 WHO-passing
  combinations; all 6 are in the all-ages setting. Because the long
  deterministic recalibration did not complete, these scenarios are not
  conditional on a newly optimized point estimate.
- A second scenario rerun using the combined NPE posterior median is in
  `output/scenarios_2026_recal_npe/`. It has 7 WHO-passing combinations, all
  in the all-ages setting. This is a posterior-median sensitivity, not a full
  uncertainty propagation over NPE/MCMC draws.

## 2026-09-30: tightened NPE and long-chain MCMC gate

- The NPE result in `output/calibration/npe_agepop_recal/` was replaced with a
  one-round, two-seed run using 3,000 simulations per seed, 5,000 draws per
  seed, 128 hidden features, 10 MAF transforms, 500 epochs, and centred
  log-scale Normal prior SD 0.15. Both SBC runs pass the 95%/90% coverage
  gate. The predictive seroprevalence BCI widths are approximately 0.127,
  0.077, 0.140, 0.148, 0.112, and 0.100 for `<20` through `60+`; the target
  widths are 0.129, 0.049, 0.069, 0.071, 0.082, and 0.106. Thus the prior
  tightening helped, but 30--49 remains wider because the current summaries
  do not identify those parameters sufficiently, and 60+ remains centred
  above target.
- An exact daily-step MCMC attempt was stopped after approximately 37 minutes
  without completing because every proposal integrates the full 150-year
  horizon at 1/365-year steps. No partial exact chain was accepted as a
  result. Subsequent validation retained the 150-year horizon but used an
  explicit monthly integration step (`dt = 0.0833`), recorded in
  `mcmc_run_config.json`.
- Monthly adaptive MCMC runs were attempted with 4 parallel chains and
  NPE-covariance initialization: 12,000 iterations (ESS 756--1,076), 60,000
  iterations (ESS 5,053--5,565), 120,000 iterations after 20,000 burn-in
  (ESS 10,305--10,848), and a 180,000-iteration run after 60,000 burn-in
  (ESS 12,558--12,998). All long runs passed ESS > 4,000, but none passed
  the strict Rhat < 1.0001 gate; the final same-start run has Rhat
  1.0000898--1.0006926. `MCMC_STRICT_GATE.txt` therefore correctly reports
  FAIL. The MCMC output is diagnostic and must not be described as strict
  convergence validation.
