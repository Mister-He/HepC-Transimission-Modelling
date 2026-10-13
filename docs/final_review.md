# Final Review and Simplification: NM → Reference MCMC → Scenarios

## 1. Objective

Reduce the project to one production inference path, one retained result set,
and exactly one independent cross-validation:

```
Nelder–Mead (multi-start)  →  prior mean  →  reference MCMC (prior-scaled initial
covariance, long adaptive random-walk chain)  →  posterior + diagnostics
                                                    →  scenarios_2026

independent cross-validation:  single-round NPE  vs.  reference MCMC
                               + Simulation-Based Calibration (SBC)
```

NPE is no longer primary and no longer part of the production path. Everything
else around it — multi-round NPE, MAP+Hessian, Pilot, and the ablation — is
legacy and must be removed. This document is a specification, not a record of
results. Every step must be executed and reported; steps that cannot be carried
out must be reported as blocked, never silently skipped.

## 2. Working principles and guardrails

- Read `AGENTS.md`, `docs/AGENTS.md`, `docs/training.md`, and each directory's
  Markdown guide before changing or deleting anything.
- Do not change the scientific model, the likelihood, the prior form, the
  equilibrium rule, or the numerical solver settings of the retained pipeline.
- Preserve all existing functionality, numerical behavior, and public
  interfaces of the retained production path. Simplification must not change
  results.
- Before deleting, inventory the candidates and justify each removal from its
  generating script, configuration, and current dependencies — not from age.
  Retain anything whose purpose is unclear and flag it in the report.
- Do not delete or overwrite user changes without identifying their exact scope.
  Run `git status` first.
- **Edit `docs/AGENTS.md` directly.** Remove the "NPE-C/SNPE is primary"
  calibration policy and make the reference MCMC the production posterior. Keep
  the model contract lines intact (strata, states, 480 compartments, 2015
  anchor, equilibrium rule).
- **Keep the new reference-MCMC defaults**: warmup 10,000, adaptation from
  proposal 5,000, storage thinning 10.
- Use `commit_rule.md` (Conventional Commits) and create one focused commit per
  step below. Do not push.
- Capture a baseline before deleting: record the current `git status`, the list
  of files to remove, and the retained posterior/scenario hashes.

## 3. Target end state

- **Single production path.** Nelder–Mead multi-start sets the prior centre; the
  reference MCMC produces the production posterior and its diagnostics.
  Multi-round NPE, MAP+Hessian, Pilot, and the ablation no longer exist in code,
  config, docs, or outputs.
- **One independent cross-validation, made comprehensive with SBC.** The
  single-round NPE posterior is compared with the reference MCMC posterior for
  parameters and policy outcomes, and Simulation-Based Calibration (SBC) checks
  whether the reported uncertainties are calibrated (empirical coverage of
  nominal intervals). These exist only to validate the inference: they must
  never feed the prior centre, the proposal covariance, the starting points, or
  the scenarios, and the production scripts must not import them.
- **Single scenario result set.** Only scenario results generated from the
  reference MCMC posterior (the prior-scaled initial-covariance chain) are
  retained.
- **Reference MCMC scope.** Keep the long adaptive random-walk Metropolis
  sampler with a prior-scaled initial covariance and the diagnostic outputs
  (trace plots, per-chain summaries, adaptation drift, covariance eigenvalues,
  jump sizes, autocorrelation times). Remove only the NPE/Hessian/Pilot hooks
  added for the ablation.
- **Documentation.** Every retained substantive directory has an accurate guide;
  no document references a deleted method, script, or output.
- **Verification.** Tests pass; the retained entry points run; the retained
  results are reproducible and are the ones described in the docs.

## 4. Part 1 — Reduce to one production path plus one cross-validation

### 4.1 Retain the single-round NPE comparison

Keep only the single-round NPE and its comparison with the reference MCMC:

- `scripts/run_npe.R` — single-round only; remove any sequential-proposal logic.
- `scripts/npe_train.py` — single-round training.
- `scripts/compare_npe_mcmc.R` — the parameter and policy comparison.
- `scripts/calibration/prior_sampling.R` — LHS prior sampling for the NPE.
- `tests/unit/test_npe_coverage.R`.
- `output/calibration/npe_lhs/` — the single-round NPE posterior.
- `output/calibration/npe_vs_mcmc/` — keep the LHS comparison; delete the
  `previous_random_npe/` subfolder and any other non-LHS comparison material.
- `figures/calibration_agepop.png`, `figures/calibration_interval_overlap.csv` —
  the NPE-versus-MCMC calibration comparison.

Confirm that `run_mcmc_reference.R`, `scenarios_2026.R`, and the manuscript build
do not import any of it.

### 4.2 Add SBC to the cross-validation

Make SBC the calibration-coverage arm of the cross-validation, using the
amortized single-round NPE (it can be re-run cheaply per synthetic dataset,
whereas SBC with the reference MCMC would require a full chain per dataset).

Procedure:

1. Generate `N_SBC` synthetic datasets from the prior: `theta ~ Normal(NM best,
   0.20)`, then `x` from the same observation model used for calibration
   (Binomial age-specific prevalence and log-normal population counts with
   `sigma = 0.10`).
2. For each synthetic dataset, draw from the single-round NPE posterior (the
   amortized network conditions on the synthetic `x`).
3. For every parameter, compute the rank of the true `theta` within the posterior
   draws. Under a well-calibrated posterior the ranks are Uniform(0, 1).
4. Report per-parameter `coverage_95` and `coverage_90`, rank histograms, and a
   uniformity test (for example Kolmogorov–Smirnov or chi-square).

Requirements and acceptance:

- Run through `Rscript scripts/run_npe.R --step all ... --n-sbc <N>` with `N` large
  enough for stable coverage (target at least several hundred; state the achieved
  number and its binomial uncertainty).
- Keep the existing gate: `coverage_95` in `[0.85, 1]` and `coverage_90` in
  `[0.75, 1]` for every parameter, and flag any parameter whose rank histogram is
  systematically skewed (U-shaped or modal), not merely borderline coverage.
- Keep `sbc_summary.csv`, `sbc_ranks.npy`, and `tests/unit/test_npe_coverage.R` as
  the retained SBC artifacts.
- Record the limitation explicitly: SBC on the amortized NPE certifies the NPE's
  uncertainty; because NPE and the reference MCMC share the prior, likelihood,
  simulator, and equilibrium rule, NPE coverage together with NPE–MCMC agreement
  is evidence for the reference MCMC's calibration, but it does not validate the
  model or the data. If budget allows, add a small direct MCMC-based SBC (re-run
  the reference MCMC on a handful of synthetic datasets) as a stronger check, and
  report clearly if it is not run.

### 4.3 Remove all other NPE / MAP+Hessian / Pilot / ablation artifacts

Inventory the following candidates, confirm each against its generator and
current dependencies, then remove the obsolete ones and update every reference.

**Code to remove (verify first):**

- `scripts/ablation_check.R`.
- The ablation-only options in `scripts/run_mcmc_reference.R` and
  `scripts/calibration/mcmc.R` that exist to inject NPE / inverse-Hessian /
  pilot inputs or starts (for example `--cov-file`, `--cov-scale`,
  `--start-file`, `--save-unthinned`). Keep the sampler, the prior-scaled initial
  covariance, and the diagnostics. Confirm no retained path uses the removed
  options before deleting them.
- Any multi-round / APT-loss NPE code, proposals, or configurations.

**Documents to remove or rewrite:**

- Remove `docs/ablation_check.md` and `docs/calibration_update.md` (the latter
  mandates retaining NPE as a production method, which is now obsolete).
- Review `docs/NPE_vs_MCMC_benchmark.md`: reduce it to the retained single-round
  cross-validation or delete it if it no longer matches, and document the
  retained validation procedure in `docs/training.md` and
  `output/calibration/npe_vs_mcmc/report.md`.
- Rewrite the NPE sections and the "NPE is primary" policy in `docs/AGENTS.md`,
  `docs/training.md`, `docs/model_card.md`, `docs/evaluation.md`,
  `docs/structured_prompt.md`, `README.md`, `scripts/README.md`,
  `output/calibration/README.md`, `output/scenarios_2026/README.md`, and
  `tests/unit/README.md` so that the reference MCMC is the production posterior
  and the single-round NPE is only the cross-validation.
- **Manuscript.** Modify every NPE reference in the manuscript pipeline
  (`scripts/build_paper.py`, `scripts/paper_figures.R`, `docs/main.docx`,
  `docs/Supplementary Information.docx`) to match the new production path, and
  rebuild the documents. `scripts/build_paper.py` is the source of the prose, so
  edit it there, not in Word.

**Outputs to remove (verify first):**

- `output/ablation_check/`.
- `output/calibration/npe_vs_mcmc/previous_random_npe/` and any other
  non-LHS comparison material.

**Outputs to keep:** `output/calibration/r_20261010_agepop/`,
`output/calibration/mcmc_reference/`, `output/calibration/npe_lhs/`,
`output/calibration/npe_vs_mcmc/` (LHS comparison only),
`output/scenarios_2026/`.

**Reference sweep:** after deletion, grep the repository (code, docs, configs,
JSON, tests) for `multi-round`, `apt`, `ablation`, `hessian`, and `pilot`, and
fix every remaining reference. No retained file may point at a deleted path.

## 5. Part 2 — Keep only reference-MCMC scenario results

- Confirm that every retained scenario output under `output/scenarios_2026/` was
  generated from `output/calibration/mcmc_reference/mcmc_draws.csv` (the
  prior-scaled initial-covariance reference MCMC posterior) and from
  `output/calibration/r_20261010_agepop/fit.rds`.
- Remove any scenario result, figure, or table generated from the NPE posterior
  or from any other posterior. Only the `output/scenarios_2026/` set may remain.
- Rewrite `output/scenarios_2026/README.md` so its provenance, method text, and
  any quoted numbers describe the reference MCMC posterior. Remove the
  comparison with the old NPE posterior and any `npe_vs_mcmc` link as a
  *production* justification; the NPE comparison may be cited only as the
  independent cross-validation.
- If the retained results and the documentation disagree, state which is
  authoritative and, if a rerun is needed, run
  `Rscript scripts/scenarios_2026.R --root .` with the retained fit and posterior
  before treating the docs as current.

## 6. Part 3 — Enforce the calibration flow

The required production flow is exactly:

1. **Nelder–Mead multi-start fit.** `scripts/run_calibration.R` runs multiple
   starts (anchor, perturbed, population-informed) in parallel and selects the
   best objective; it writes `fit.rds`, whose `best$theta` is the prior centre.
2. **Reference MCMC.** `scripts/run_mcmc_reference.R` builds the prior
   `Normal(fit$best$theta, 0.20)` on the 12 log-scale parameters, initialises the
   chains from prior draws with a prior-scaled initial proposal covariance
   (`0.1^2 × prior variance`), adapts during warmup, freezes the kernel, and
   retains draws until the evaluation budget. It writes the production posterior
   (`mcmc_draws.csv`) and the diagnostics.

Requirements:

- Verify the retained code implements exactly this flow and nothing else.
- Keep the new defaults: warmup 10,000, adaptation from proposal 5,000, storage
  thinning 10. Document them in `docs/training.md` and in the reference MCMC run
  config. Note the provenance of the stored posterior (produced with the earlier
  settings `--n-warmup 4999 --adapt-start 200 --thin 1`).
- Keep the reference MCMC diagnostics and the check that all chains produced
  retained output (R-hat, bulk/tail ESS, MCSE, acceptance, trace plots,
  adaptation drift, covariance condition numbers).
- Remove any code path that lets NPE, MAP+Hessian, or Pilot seed the prior mean,
  the initial covariance, or the starting points.

### 6.1 Harmonize rearrest-rate parameters and refresh dependent results

This parameterization update is an explicit, scoped change to the scientific
model. It supersedes the no-model-change and unchanged-results guardrails above
for rearrest rates and every calibration or scenario result that depends on
them. Keep the remaining model contract and inference workflow intact.

1. **Use one name for the composite parameter.** In both the active files under
   `scripts/model/` and the legacy top-level copies, rename the Stan/preprocessing
   parameter `c` to `c_composite`. It is the same age-specific baseline
   calibration quantity represented by `c_composite` in `setup.R` and consumed
   by `sim.cpp`, not a second independent parameter. Keep Stan initialization,
   summaries, saved draws, and downstream extraction on the same name.
2. **Automate the rate-fit handoff.** Add a script that reads the relevant
   `output/*_rate_fit.rds` files and updates only the corresponding rate and
   composite-parameter entries in `scripts/model/setup.R`: map the release fit's
   log-rate summary to `lambda2`, and the rearrest fit summaries to `lambda3`,
   `pi_recid`, and `c_composite`. Preserve the current summary conventions
   (including exponentiating the posterior mean of the log-rate where used),
   validate vector lengths, parameter domains, and finite values, and fail
   without writing if an expected parameter or unique setup block is missing.
   Write the updated file only after validation; this replaces manual copying
   of fitted values.
3. **Enforce the composite relation in calibration.** Define
   `c_composite = c * beta_scale`, so the effective age-specific scaling is
   `c = c_composite / beta_scale`. Derive the first-arrest rate in the simulator
   as `lambda1 = lambda3 * c`; every calibration proposal must satisfy the
   composite relation by construction. Treat the Stan estimate as the
   `beta_scale = 1` baseline composite, and apply the calibrated beta scale
   consistently to inflows and the derived first-arrest rate.
4. **Regenerate model-dependent results after cleanup.** Once the obsolete
   methods and outputs in Part 1 are removed and the parameter handoff is in
   place, rerun the release/rearrest rate fits, update `setup.R` through the new
   script, then rerun the deterministic calibration, reference MCMC, retained
   single-round NPE/SBC cross-validation, and scenarios in dependency order.
   Replace stale calibration and scenario outputs with results from this model
   version, record their provenance, and update the reports. Do not describe
   old results as unchanged or reproducible under the new parameterization.

## 7. Part 4 — Verify calibration and scenarios

After Parts 1–3, run and report:

For outputs affected by §6.1, verification means successfully regenerating and
validating the new model-version outputs; equality with the old retained fit,
posterior, or scenarios is not expected. The existing reproducibility checks
below apply to the newly generated artifacts.

- `Rscript scripts/run_tests.R --scope all` — all retained tests pass.
- Rerun deterministic calibration with `Rscript scripts/run_calibration.R
  --run-id <id> ...`; retain the new `fit.rds` and its equilibrium diagnostics.
- Rerun reference MCMC with `Rscript scripts/run_mcmc_reference.R --fit ...`
  and the new defaults; retain the new posterior and confirm every chain passes
  its diagnostics.
- Rerun scenarios with `Rscript scripts/scenarios_2026.R --root .` from the new
  fit and posterior, and verify their provenance.
- Rerun the single-round NPE cross-validation with `Rscript
  scripts/run_npe.R ...` and `Rscript scripts/compare_npe_mcmc.R`; update the
  comparison report for the new model version.
- Rerun SBC with `Rscript scripts/run_npe.R --step all ... --n-sbc <N>`,
  with every parameter passing the coverage gate and no systematically skewed
  rank histogram. Report the achieved `N` and the per-parameter coverage.
- Static checks: no retained file references a deleted method, script, or
  output; every retained substantive directory has an accurate guide; the
  scripts, docs, and results describe only NM → reference MCMC → scenarios plus
  the one NPE cross-validation.
- Report any mismatch between the documented flow and the actual code as a
  finding, not as a silent fix.

## 8. Part 5 — Ponytail audit, refactor, and review

Run these in order.

1. **Audit.**
   ```
   ponytail-audit
   ```
   Have the agent analyse the whole project and produce a list of
   optimisation/simplification opportunities for the retained pipeline.

2. **Refactor.** Run `ponytail-full` and give the agent exactly this prompt:

   ```
   Based on the previous Ponytail audit, implement the
   recommended code simplifications.

   Prioritize high-impact, low-risk improvements.
   Preserve all existing functionality, numerical behavior,
   and public interfaces.

   Implement changes incrementally and run the relevant
   tests after each batch.
   ```

3. **Review.** After the refactor completes, run
   ```
   ponytail-review
   ```
   to check the changes are correct, then have the agent fix every issue found.
   Re-run `Rscript scripts/run_tests.R --scope all` after the fixes.

## 9. Commits (per `commit_rule.md`, no push)

Create one focused Conventional Commit per step, for example:

- `chore: remove ablation, Hessian and pilot inference paths and outputs`
- `refactor(npe): keep a single-round cross-validation only`
- `test(npe): enable and document SBC coverage in the cross-validation`
- `docs: make the reference MCMC the sole production posterior`
- `refactor(mcmc): drop ablation-only proposal and start hooks`
- `docs: rewrite scenarios provenance for the reference MCMC`
- `docs(manuscript): rebuild without NPE production references`
- `refactor: apply ponytail audit simplifications`
- `fix: address ponytail review findings`

Follow the format `<type>(scope): <description>`, mark any breaking change with
`!` and a `BREAKING CHANGE:` footer, review the staged diff before committing,
and exclude unrelated changes, temporary files, and secrets. **Do not push.**

## 10. Deliverables

- The cleaned repository: one production path, one retained cross-validation,
  one scenario result set, updated docs and manuscript.
- A short final report stating:
  - what was removed and what was retained, with reasons;
  - the finalized calibration flow and the exact commands for NM → reference MCMC
    → scenarios, and for the single-round NPE cross-validation and its SBC;
  - the verification results (tests, fit, MCMC diagnostics, scenarios, NPE
    comparison, SBC coverage);
  - the Ponytail audit findings, the simplifications applied, and the review
    outcome;
  - anything retained but flagged as uncertain, and any blocked step;
  - the commit list.

## 11. Key files

- `scripts/run_calibration.R`, `scripts/calibration/calibrate_nm.R` — NM prior centre.
- `scripts/run_mcmc_reference.R`, `scripts/calibration/mcmc.R` — reference MCMC and diagnostics.
- `scripts/scenarios_2026.R`, `scripts/post_calibration_analysis.R` — scenarios.
- `scripts/run_npe.R`, `scripts/npe_train.py`, `scripts/compare_npe_mcmc.R`, `scripts/calibration/prior_sampling.R` — the single independent cross-validation (single-round NPE vs reference MCMC).
- `output/calibration/npe_lhs/sbc_summary.csv`, `output/calibration/npe_lhs/sbc_ranks.npy`, `tests/unit/test_npe_coverage.R` — the retained SBC coverage check.
- `scripts/predict_2015_posterior.R`, `scripts/plot_calibration.R` — 2015 predictions and figure.
- `output/calibration/r_20261010_agepop/`, `output/calibration/mcmc_reference/`, `output/calibration/npe_lhs/`, `output/calibration/npe_vs_mcmc/`, `output/scenarios_2026/` — retained results.
- `docs/AGENTS.md`, `docs/training.md`, `docs/model_card.md`, `docs/evaluation.md` — calibration policy and methods.
- `commit_rule.md` — commit format.
