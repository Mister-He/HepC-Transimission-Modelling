# First-arrest rate back to lambda1 = lambda3

Structured prompt (2026-10-20). It supersedes §3.1 (the first-arrest handoff) and the entry
anchor of `docs/first_arrest_mortality_update.md`. Its §3.2 (viraemic excess mortality) and
everything else in `docs/final_review.md` and `docs/AGENTS.md` still apply.

## 1. Objective

- Return the first-arrest rate to the pre-`7c3dad3` form: `lambda1` is hard-coded in `setup.R`
  with the same values as `lambda3`.
- Remove everything that existed only because `lambda1` was derived from `f_first`.
- Change every script, test and document to match, then regenerate all results with the current
  pipeline (`docs/training.md`).
- Use as few commits as possible, each following `commit_rule.md`. Do not push.

## 2. Decisions (user, 2026-10-20)

1. **First-arrest rate:** `lambda1_i = lambda3_i` for every age, as a hard-coded vector in
   `setup.R` (current SPS fit: 0.546, 0.488, 0.472, 0.472, 0.386, 0.266).
2. **Entry anchor removed too.** θ7–θ12 return to plain log entry scales,
   `entry_i = exp(θ_{6+i}) · beta_i` (`beta` = CNB 2015 counts), with bounds 0.01–10⁶. No entry
   term in the likelihood, the MCMC target, the NPE summaries or the MCMC SBC.
3. **Viraemic excess mortality kept:** `hr_viraemic` = 1.47 and the frozen `viraemic_share_ref`,
   recomputed by the fixed point after the refit.
4. **Full chain rerun**, including NPE + SBC, NPE vs MCMC, the MCMC SBC and the HR sensitivity
   fits.

## 3. Requirements

### 3.1 Model and rates
- **`scripts/model/sim.cpp`:**
  - read `lambda1` from the parameters (length 6, finite, non-negative);
  - remove `first_arrest_rates()` and the exported `first_arrest_rate()`;
  - keep `AGE_RATE` and the viraemic mortality.
- **`scripts/model/setup.R`:** replace `f_first` with
  `lambda1 = c(0.546, 0.488, 0.472, 0.472, 0.386, 0.266)`, with a comment saying it equals
  `lambda3` by assumption. Fix the stale `pi_recid` 0.59 comment. Leave the user's whitespace
  edit in the index comment alone.
- **`scripts/model/update_rate_params.R`:** stop writing `f_first`; write `lambda1` with the same
  values as `lambda3`. It must still leave `setup.R` unchanged on the stored fits.
- **Stan and `preprocess.R`: unchanged.** `f_pop` stays inside `rearrest_rate.stan`, where its
  multinomial term on the first-offender age shape informs `lambda3`; it is no longer handed to
  the ODE. The rates do not change, so no Stan refit is needed. `preprocess.R` holds an
  uncommitted user edit and is not touched.

### 3.2 Calibration
- **`scripts/calibration/likelihood.R`:**
  - remove `entry_scales()`, `excess_theta()`, `nll_entry()` and `AGE_RATE`;
  - `build_params()` sets `beta = beta_base · exp(θ7–12)`;
  - `BETA_LO` = 0.01.
- **`scripts/calibration/mcmc.R`:** drop `nll_entry` from `log_posterior_exact()`.
- **`scripts/calibration/targets.R` and `data/calibration_targets.json`:** remove `entry` and
  `entry_sigma_log`.
- **`scripts/calibration/calibrate_nm.R`:** rewrite `informed_theta0()` as a steady-state
  balance:
  - J targets with `lambda2`, `lambda3`, `pi` and ageing give F and the first-arrest flow;
  - divide by `lambda1` for D;
  - the D balance gives entry.
- **`scripts/run_calibration.R`:**
  - starts are entry scales directly (no `excess_theta` conversion);
  - `solutions.csv` reports `beta_scale = exp(θ7–12)`.
- **NPE back to 12 summaries:**
  - `scripts/run_npe.R`: `simulate()` returns `c(p_sero, N_hat)`; `x_obs` has 12 columns; drop
    `--entry-sigma-log`;
  - `scripts/npe_train.py`: drop the entry noise and the `--entry-sigma-log` option;
    `transform_x` uses columns 7–12.
- **`scripts/sbc_mcmc.R`:** drop the entry observation.
- **`scripts/scenarios_2026.R`:** `j_positive_entries()` uses `pm$lambda1`.
- **Fit path:** the new Nelder–Mead run is `output/calibration/r_20261020_agepop`. Update the
  default `--fit` in every script and in `build_paper.py`, and delete the superseded
  `r_20261014_agepop/`.

### 3.3 Tests
- Remove the `f_first`, `first_arrest_rate`, `entry_scales`/`excess_theta` and `nll_entry`
  checks, and the CNB entry check in `test_targets.R`.
- Add:
  - `setup.R` `lambda1` equals `lambda3`;
  - `run_sim` uses `lambda1`: doubling it lowers D and raises J in a short run without
    infection;
  - `build_params()` multiplies `beta` by `exp(θ7–12)`.
- `Rscript scripts/run_tests.R --scope all` must pass.

### 3.4 Rerun chain (`docs/training.md` order)
Launch each long run detached (`nohup`) from a frozen copy of the scripts (Rscript reads
lazily).

| Step | Command (summary) | Cores | Approx. time |
| --- | --- | --- | --- |
| 1. Nelder–Mead | `run_calibration.R --run-id r_20261020_agepop --maxit 3000 --nm-restarts 5 --n-starts 7 --sd-perturb 0.9 --seed 20261001` | 8 | 2.5 h |
| 1b. Viraemic share fixed point | Compare `viraemic_share_2015` with `viraemic_share_ref`. If the implied `omega_0` moves ≥ 0.5%, set the new value in `setup.R` and refit. | 8 | 0–2.5 h |
| 2. Reference MCMC | `run_mcmc_reference.R --cores 4` (defaults) | 4 | 4.5 h |
| 2b. HR sensitivity (parallel with 2) | HR 1, 1.11, 1.23, 1.77; `--n-starts 2 --maxit 3000 --nm-restarts 5`; keep 4 files each | 12 | 2.5 h |
| 3. 2015 predictions, calibration figure | `predict_2015_posterior.R --max-draws 4000`, `plot_calibration.R --max-draws 2000` | 14 | 0.5 h |
| 4. Scenarios | `scenarios_2026.R --draws-per-seed 2000 --cores 14 --sample-seed 20261005` | 14 | 0.5 h |
| 5. NPE + SBC | `run_npe.R --step all --n-sims 50000 --seed 20261012 --train-seed 20261013 --n-sbc 1000` | 14 | 1.5 h |
| 6. NPE vs MCMC | `compare_npe_mcmc.R --policy-cores 14` | 14 | 0.3 h |
| 7. MCMC SBC | `sbc_mcmc.R --n-datasets 8 --max-evals 30000 --cores 8` | 8 | 3.5 h |
| 8. Figures and manuscript | `python scripts/figures.py .`, then `python3 scripts/build_paper.py .` | 1 | minutes |

### 3.5 Documents
- **Update to the new model and numbers:**
  - `docs/training.md`;
  - `scripts/model/README.md`, `scripts/README.md`;
  - `output/calibration/README.md`, `output/calibration/hr_sensitivity/README.md`;
  - `output/calibration/npe_vs_mcmc/{report.md,model_audit.md}`;
  - `output/scenarios_2026/README.md` (Chinese);
  - `data/README.md`, `docs/evaluation.md` (where it describes the handoff);
  - `claude.md`.
- **`scripts/build_paper.py`:** methods describe `lambda1 = lambda3` and free entry scales, with
  no anchor. Every result sentence and assertion must follow the new outputs; revise the
  narrative if the assertions fail.
- **Historical records stay as they are:** `docs/first_arrest_mortality_update.md` and
  `docs/final_review.md`. Their supersession is noted here and in `claude.md`.

## 4. Verification
- Tests pass.
- Nelder–Mead converges; the best fit passes the equilibrium check.
- MCMC R-hat ≤ 1.01; the NPE SBC gate passes (or failures are reported as they are).
- Every number in the manuscript comes from the new outputs, and `build_paper.py` runs without
  failed assertions.
- No file still refers to `f_first`, `first_arrest_rate`, `entry_scales`, `nll_entry` or
  `entry_sigma_log`, except the historical docs.

## 5. Commits (no push)
1. `feat(figures): …`: the pending scientific-visualization redraw work (`scripts/figures.py`,
   the two cost CSV writers, the SVGs and regenerated PNGs, its docs).
2. `feat(model)!: …`: `lambda1 = lambda3`, no entry anchor; code, tests and code-level docs, with
   a `BREAKING CHANGE` footer (θ7–θ12 meaning, targets file, old fits incompatible).
3. `feat(calibration): …`: all regenerated outputs, documents and the manuscript.

Unrelated user changes stay unstaged:
- `Hepc_modelling_technical_details.pdf`;
- the deletions of `codex.md`, `route.md` and `old_trajectory.png`;
- `validity_audit_2026-10-11.md`;
- the whitespace edits in `preprocess.R` and `setup.R`;
- `r_20261014_agepop/targets.json`, which is deleted with its folder in commit 3.
