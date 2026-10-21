# Data guide

Input files, not generated results. `calibration_targets.json` (2015 age-specific
seroprevalence targets and the fixed precision `beta_precision` of their Beta likelihood; no
counts or population sizes; values rebuilt from the observation vector
in an earlier NPE run's `x_obs.csv`; `output/calibration/npe_lhs/x_obs.csv` is regenerated from this file), `singstat_new_drug_abusers.csv`
(SingStat/CNB new drug abusers: the 2015 row is the entry into D, `beta`, used without
scaling; later rows are the 2016–2025 inflow schedule), and `scenario_policy_config.json` (2026 policy and cost
parameters, reconstructed from `docs/structured_prompt.md` and
`output/scenarios_2026/cost_evidence.csv`). `sps.csv` (Singapore Prison Service records; not tracked by git) is read only by
`scripts/model/preprocess.R` for the rate fits.
