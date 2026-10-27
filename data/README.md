# Data guide

Input files, not generated results. `calibration_targets.json` holds the three calibration targets:
2015 age-specific prison seroprevalence with the fixed precision `beta_precision` of its Beta
likelihood (values rebuilt from an earlier NPE run's `x_obs.csv`); `first_arrests`, the CNB 2015
new drug abusers by age (257, 644, 292, 94, 21, 1), with the negative binomial size
`first_arrests_nb_size`; and `community_pwid`, the Teo et al. 2019 network scale-up estimate of
injecting drug users in Singapore (11,000, 95% CI 6,500–17,000, survey year
`community_pwid_year` 2017), with log-normal SD `community_pwid_log_sd` taken from the CI.
`output/calibration/npe_lhs/x_obs.csv` is regenerated from this file. `singstat_new_drug_abusers.csv`
(SingStat/CNB new drug abusers by year; only the 2015 row is used: `beta` in `setup.R`; entry into
D is `beta` times the calibrated entry scales θ7–θ12 (≥ 1), held constant in every model year
including 2016–2030), and `scenario_policy_config.json` (2026 policy and cost
parameters, reconstructed from `docs/structured_prompt.md` and
`output/scenarios_2026/cost_evidence.csv`). `sps.csv` (Singapore Prison Service records; not tracked by git) is read only by
`scripts/model/preprocess.R` for the rate fits.
