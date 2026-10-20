# Data guide

Input files, not generated results. `calibration_targets.json` (2015 age-specific
seroprevalence and prison-population targets; rebuilt from the observation vector
in an earlier NPE run's `x_obs.csv`; `output/calibration/npe_lhs/x_obs.csv` is regenerated from this file), `singstat_new_drug_abusers.csv`
(SingStat/CNB new drug abusers: the 2015 row is the base entry into D, `beta`, scaled by the
calibrated entry scales; later rows are the 2016–2025 inflow schedule), and `scenario_policy_config.json` (2026 policy and cost
parameters, reconstructed from `docs/structured_prompt.md` and
`output/scenarios_2026/cost_evidence.csv`). `sps.csv` (Singapore Prison Service records; not tracked by git) is read only by
`scripts/model/preprocess.R` for the rate fits.
