# Data guide

Input files, not generated results. `calibration_targets.json` holds the three calibration targets:
2015 age-specific prison seroprevalence with the fixed precision `beta_precision` of its Beta
likelihood (values rebuilt from an earlier NPE run's `x_obs.csv`); `first_arrests`, the CNB 2015
new drug abusers by age (257, 644, 292, 94, 21, 1) times the share whose main drug was heroin or
methamphetamine (1,127 / 1,309 = 0.861; `singstat_drug_abusers_by_drug.csv`), rounded to 221, 554,
251, 81, 18, 1, with the negative binomial size `first_arrests_nb_size`; and `community_pwid`, the
Teo et al. 2019 network scale-up estimate of injecting drug users in Singapore (11,000, 95% CI
6,500–17,000, survey year `community_pwid_year` 2017), with log-normal SD `community_pwid_log_sd`
taken from the CI.

- `singstat_new_drug_abusers.csv`: SingStat/CNB new drug abusers by year and age. Only the 2015 row
  is used: `beta` in `setup.R` is these counts times 0.861 (unrounded). Entry into D is `beta` times
  the calibrated entry scales θ7–θ12 (≥ 1), held constant in every model year including 2016–2030.
- `singstat_drug_abusers_by_drug.csv`: SingStat table M891521 (CNB, number of drug abusers by main
  drug abused), 2015 values for all and new abusers, downloaded 2026-10-30 from
  `https://tablebuilder.singstat.gov.sg/api/table/tabledata/M891521`. Heroin and methamphetamine
  are taken as the injected drugs; the table has no age split, so one share applies to all ages.
- `scenario_policy_config.json`: 2026 policy and cost parameters, reconstructed from
  `docs/structured_prompt.md` and `output/scenarios_2026/cost_evidence.csv`. `end_year` (2050) is the
  last projection year; the WHO targets are still checked in 2030 (`WHO_YEAR` in `scripts/scenarios_2026.R`).
- `sps.csv` (Singapore Prison Service records; not tracked by git) is read only by
  `scripts/model/preprocess.R` for the rate fits.
