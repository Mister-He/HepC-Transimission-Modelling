# First-arrest multiplier c, fixed CNB entry, PWID mortality ratios

Structured prompt (2026-10-21). It supersedes `docs/lambda1_revert.md` (λ1 = λ3, free entry
scales) and the viraemic-share fixed point of `docs/first_arrest_mortality_update.md` §3.2.
Everything else in `docs/AGENTS.md` and `docs/training.md` still applies.

## 1. Objective
- Rates: the SPS records estimate only the release rate (λ2), the re-arrest rate (λ3) and
  `pi_recid`. Refit them and hard-code them in `setup.R`.
- Calibrate 7 parameters: six contact scales and one first-arrest multiplier `c`, with
  `λ1_i = c · λ3_i` for every age.
- Entry into D is the CNB 2015 new drug abuser count, with no scaling (2016–2025: the CNB rows,
  unscaled).
- Mortality: the general-population rate times a fixed PWID mortality ratio, one for viraemic
  and one for non-viraemic states. Drop the SMR 14.68 and the frozen viraemic share.
- Rerun the whole chain, then the scenarios, figures, manuscript and SI. Then run the audit,
  refactor and review in `docs/instructions.md` §8.
- Use as few commits as possible (`commit_rule.md`). Do not push.

## 2. Decisions
1. **Rates (user):**
   - `rearrest_rate.stan` keeps only the cure-type survival model for `pi` and `λ3`. Remove
     the 2015 equilibrium block: `J`, `F0`, `mu`, `a`, `j_total`, `j_props`, `r_pop`, `f_pop`.
   - `preprocess.R` drops the inputs that only that block needed (`J_pop`, `f_yr_dist`, `mu`).
   - The release model is unchanged. Release and re-arrest still use all consumption records.
2. **Parameters (user):** θ1–θ6 = log contact scales; θ7 = log `c`. `N_THETA` = 7.
3. **Prior (user, 2026-10-21: fixed and more informative; no data-dependent centre):**
   - Independent Normal on θ, truncated to the bounds, shared by Nelder–Mead, the MCMC and NPE
     (`prior_spec()` in `likelihood.R`).
   - θ1–θ6 (log contact scales relative to `setup.R`'s `C_contact`): `Normal(0, 1.5)`; 95%
     range of a scale 0.05–19.
   - θ7 = log `c`: `Normal(0, 0.5)`, centred on λ1 = λ3. Bounds: **c in 0.5–100** (user: c must
     be above 50%).
   - Nelder–Mead minimises the negative log posterior, so the best fit is the posterior mode.
     `prior_centre()` and the contact `pattern_penalty()` are removed.
4. **Likelihood (user, 2026-10-21): seroprevalence only, no counts at all.**
   - `y_i ~ Beta(p_i·φ, (1 − p_i)·φ)` for each age group, with `p_i` the model's 2015 prison
     seroprevalence and `φ` = 100 fixed (`beta_precision` in `calibration_targets.json`; SD
     about 0.045 at p = 0.3).
   - `calibration_targets.json` holds only `age_groups`, `seroprevalence`, `beta_precision`
     and the equilibrium settings. `population` and `population_sigma_log` are removed;
     modelled prison size is reported only as a model output.
   - `c` is informed by seroprevalence through the D/J mix, and otherwise by its prior.
   - NPE uses the six prevalence summaries with Beta noise; the MCMC SBC simulates Beta data.
   - History: the fit with the population target converged at c = 0.034 (objective 119,
     prison 20–29 381 against 1,001, 60+ 318 against 147). A binomial with the screened
     counts and a fractional Bernoulli were considered and rejected by the user.

5. **Mortality ratios (searched 2026-10-21):**
   - Source: van Santen et al. 2014, *AIDS* 28:2589–99 (doi:10.1097/QAD.0000000000000450),
     Amsterdam Cohort Studies, 1985–2012, Table 3. These are all-cause SMRs against the
     age-, sex- and period-matched general Dutch population.
   - **Non-viraemic:** HCV/HIV-uninfected **4.73 (3.47–6.45)**.
   - **Viraemic:** HCV-monoinfected **9.52 (8.06–11.26)**.
   - No cohort reports SMRs by HCV RNA status. The anti-HCV-positive group is used for
     viraemic as the closest available (about 75% are chronic).
   - **States:** viraemic = acute, chronic and on treatment (h = a, c, t). Non-viraemic =
     uninfected and cleared or SVR (u, r), as Lee 2012 found RNA-negative people at negative
     risk.
   - **Formula:** `background = mu_i × SMR_state`, with DC/HCC liver deaths added on top as
     before.
   - **Share:** the viraemic share is no longer a parameter. It follows the model state in
     every year, so treatment lowers mortality as people clear.
   - **Sensitivity:** refit at both CI bounds (low 3.47/8.06, high 6.45/11.26), replacing
     the earlier HR sensitivity.
6. **Scenarios:** use the same per-state mortality, so the viraemic share changes year by
   year with treatment.

## 3. Requirements
### 3.1 Model and rates
- **`sim.cpp`:**
  - read `c_first` and set `λ1 = c_first · λ3`;
  - read `smr_viraemic` and `smr_nonviraemic`;
  - remove `omega`, `omega_neg`, `hr_viraemic`, `viraemic_share_ref` and the `lambda1`
    input.
- **`setup.R`:**
  - remove `omega`, `hr_viraemic`, `viraemic_share_ref` and `lambda1`;
  - add `smr_viraemic = 9.52`, `smr_nonviraemic = 4.73` and `c_first = 1` (a placeholder,
    since it is calibrated);
  - write the new rates.
- **`update_rate_params.R`:** write `lambda2`, `lambda3` and `pi_recid` only.
- **`preprocess.R` and `rearrest_rate.stan`:** as in §2.1. Refit with
  `Rscript scripts/model/preprocess.R`, then `update_rate_params.R`.

### 3.2 Calibration
- **`likelihood.R`:**
  - `N_THETA` = 7, with bounds `C_LO`/`C_HI`;
  - `build_params` sets `c_first = exp(θ7)` and leaves `beta` at CNB;
  - the pattern penalty keeps only the contact terms;
  - add `log_prior_c()`, used by Nelder–Mead, the MCMC and NPE.
- **`mcmc.R` and `run_mcmc_reference.R`:** the prior mean and SD vectors become
  (centre, 0.2) for θ1–θ6 and (0, 1) for θ7.
- **`calibrate_nm.R`:**
  - `informed_theta0` picks `c` by a 1-D search that matches the J-implied first-arrest flow
    from the D steady state with CNB entry;
  - the starts are 7-dimensional.
- **`run_calibration.R`:**
  - `solutions.csv` reports `c` in place of `beta_scale`;
  - the HR flag becomes `--smr-viraemic`/`--smr-nonviraemic`;
  - drop the viraemic-share fixed point.
- **`run_npe.R`, `npe_train.py`, `sbc_mcmc.R`, `compare_npe_mcmc.R`,
  `predict_2015_posterior.R`, `plot_calibration.R`, `scenarios_2026.R`, `figures.py`:**
  7 parameters; θ7 named `c`. Pair plots: (θ1, θ7), (θ6, θ7), (θ2, θ3), (θ4, θ5).
- **`model_metrics.R`:** `viraemic_share` counts a, c and t; it is reported, not frozen.
- **New fit:** `output/calibration/r_20261021_agepop`. Delete `r_20261020_agepop`. Move the
  SMR sensitivity fits to `output/calibration/smr_sensitivity/`, replacing `hr_sensitivity/`.

### 3.3 Tests
- Update the `N_THETA`, `build_params`, prior and mortality tests:
  - with both SMRs at 1, every state dies at `mu`;
  - a viraemic person dies faster than a non-viraemic one;
  - `λ1 = c·λ3`.
- `Rscript scripts/run_tests.R --scope all` passes.

### 3.4 Rerun chain (launch with `nohup` from frozen copies)
1. Stan rate refit (`preprocess.R`), then `update_rate_params.R`.
2. **Pilot:** a short Nelder–Mead run to check that the fit is feasible. Fit quality is
   reported as it is; the model is not changed to rescue a poor fit.
3. Nelder–Mead (8 starts, maxit 3000, `--nm-restarts 5`). Starts are uninformative (user):
   θ = 0 (contact scales 1, c = 1) and 7 Normal perturbations of it (SD 0.9), clipped to the
   bounds.
4. These three in parallel:
   - the reference MCMC;
   - NPE + SBC;
   - the MCMC SBC and the SMR-bound sensitivity fits.
5. 2015 predictions, the calibration figure, the scenarios (2,000 draws), NPE vs MCMC.
6. `figures.py`, then `build_paper.py`. Rewrite the narrative and assertions to match the
   results.

### 3.5 Documents
- Update:
  - `docs/training.md`, `docs/evaluation.md`;
  - `scripts/model/README.md`, `scripts/README.md`, `data/README.md`;
  - `output/README.md`, `output/calibration/README.md`, `output/scenarios_2026/README.md`,
    `README.zh-CN.md`;
  - `output/calibration/npe_vs_mcmc/{report.md,model_audit.md}`;
  - `claude.md`.
- Add `vansanten2014` to `references.bib`. Drop the `omega`/Mathers row from the manuscript's
  parameter table (Mathers may still be cited in the discussion).
- Historical docs stay as they are.

## 4. Verification
- The tests pass.
- Nelder–Mead converges and passes the equilibrium check.
- MCMC R-hat ≤ 1.01.
- The NPE SBC gate is reported as it is.
- `build_paper.py` runs with no failed assertion.
- No file outside the historical docs still mentions `omega`, `hr_viraemic`,
  `viraemic_share_ref`, `beta_scale`, `f_pop` or `lambda1 =`.

## 5. Audit and commits
- Run `docs/instructions.md` §8: ponytail-audit, then ponytail-full refactor (preserving
  numerical behaviour; tests after each batch), then ponytail-review, and fix its findings.
- Commits:
  1. `feat(model)!: estimate first arrest as a multiple of re-arrest, fix CNB entry and use
     PWID mortality ratios`: code, tests and code-level docs, with a `BREAKING CHANGE` footer.
  2. `feat(calibration): refit rates and regenerate calibration, scenarios and manuscript`.
  3. `refactor: apply ponytail audit simplifications`.
  4. `fix: address ponytail review findings` (only if there are findings).
- Leave alone: the user's deletion of `docs/final_review.md`, and the untracked
  `docs/instructions.md`.
