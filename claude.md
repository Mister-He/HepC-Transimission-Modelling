# claude.md — handoff notes (2026-10-05)

## Task done
The 17 policy scenarios (J-only, 15 finite F/X screening policies, no-cap reference) were regenerated so that every
result is the **median over posterior uncertainty with a 95% uncertainty interval (UI = 2.5th–97.5th percentile)**,
aggregating (1) posterior parameter draws and (2) the two NPE seeds. Before this, `codex.md` describes a single
parameter point (pooled NPE median) with no intervals.

## How it works (`scripts/scenarios_2026.R`)
- Draws: `--draws-per-seed` (default 500; the delivered run used **1000**) random rows from each of
  `posterior_samples_npe.csv` and `posterior_samples_npe_seed2.csv` in `output/calibration/npe_agepop_recal/`,
  pooled (2000 draws, `--sample-seed 20261005`). "Seeds" = the two NPE training seeds; the ODE model itself is deterministic.
- `evaluate_draw(theta)` runs the whole pipeline per draw: build params, 150-yr equilibrium check, 2016–2025 common
  history (J from 2019), all 17 scenarios 2026–2030, costs, per-draw consistency `stopifnot`s. Run in parallel with
  `parallel::mclapply` (`--cores`, 14 used). ~18 min for 2000 draws on 14 cores.
- `agg_draws()` aggregates per-draw tables to median + `_lo` + `_hi`; logical columns become `*_prob`.
  Derived quantities (reductions, costs, incremental cost vs J-only) are computed **per draw first**, then summarised.
- PASS/FAIL (80% incidence, 65% mortality reduction vs 2015) is decided from the **median** reductions (what the
  figures show); `joint_pass_probability` = share of draws passing both.
- Outputs (all in `output/scenarios_2026/`, overwritten): `scenario_2030_summary.csv` (median/_lo/_hi, `*_text`
  columns formatted "xx (95% UI: xx–xx)", sensitivity costs), `scenario_annual_services_costs.csv`,
  `scenario_trajectories.csv`, `scenario_treatment_summary.csv`, `target_scenario_treatment_summary.csv`,
  `scenario_results.xlsx` (sheet `BySeed` = per-seed medians), `calibration_source.txt`, all PNGs (medians only,
  no UI drawn), and the updated `README.md` (table in UI format, Chinese).
- Rerun: `Rscript scripts/scenarios_2026.R --root . --draws-per-seed 1000 --cores 14`.

## Missing `data/` reconstructed
The worktree had no `data/`. Created: `data/scenario_policy_config.json` (from `docs/structured_prompt.md` +
`cost_evidence.csv`), `data/singstat_new_drug_abusers.csv` (copied from the sibling `HepC-Transimission-Modelling`
checkout), `data/calibration_targets.json` (rebuilt from `npe_agepop_recal/x_obs.csv`; the sibling checkout's copy is
an older target set), `data/README.md`, and a `.gitignore` whitelist for the config. **Validation:** the old
single-point script run with these files reproduced the previous summary exactly (max difference 0). The original
outputs were backed up to the session scratchpad only (not in the repo).

## Key results (median, 95% UI)
- J-only: incidence −75.9% (70.5–79.8), mortality −42.6% (41.8–43.5); FAIL.
- Median PASS: 5000/yr F25% (93.3%, 65.1% [63.7–66.4]; only **55.9%** of draws pass → not robust), 8000/yr F25%
  (100% of draws), 8000/yr F50% (91.9%), and the no-cap reference. The other 12 finite policies fail the mortality
  target. Incidence target is met by most finite policies from ≥600/yr with F share ≥25%.
- Seeds agree closely (e.g. 5000 F25% mortality 65.06% vs 65.14%). Full table: README / summary CSV.

## Caveats / next steps
- UIs reflect **posterior parameter uncertainty only**; not drug price, screening uptake, linkage, or capacity.
- Stacked cost-bar segments are per-stratum medians and pie shares are normalised item medians; labelled totals are
  median totals, so they may differ slightly. Cascade/annual figures sum medians.
- Tests: `Rscript scripts/run_tests.R --scope all` was rerun after the refactor: 6 files, 62 checks, 0 failures.
  None of them exercise `scripts/scenarios_2026.R`; that was validated by the exact match to the old single-point output.
- Two old PNGs (`scenario_2025_nature_style_*.png`) show as deleted in git; they were already gone before this session's run.
- Stale text still based on the single-point analysis: `结论.md` (Chinese conclusions), `codex.md`, and the
  "NPE median" wording in `docs/structured_prompt.md`. The README paragraph on why low-capacity starts are small
  quotes the old single-point numbers (flagged as such). Update these if the user wants consistency.
- `docs/AGENTS.md` says not to run post-calibration analysis unless asked; the user asked here.

## Manuscript (2026-10-05, second task)
- Deliverables: `docs/main.docx` (PLOS Medicine layout: structured abstract, author summary, IMRaD,
  statements, figure legends, Table 1) and `docs/Supplementary Information.docx` (S1 Text: model,
  calibration, scenarios, costs, all results with 95% UIs, by-seed table, cost sensitivity, CHEERS 2022).
  References: `docs/references.bib` (43 entries; journal metadata pulled from Crossref by DOI, web
  sources added by hand) rendered with `docs/plos.csl` (PLOS Vancouver style).
- Rebuild: `Rscript scripts/scenarios_2026.R --root . --draws-per-seed 1000 --cores 14` (~18 min) →
  `Rscript scripts/paper_figures.R .` → `python3 scripts/build_paper.py .`. Prose lives in
  `scripts/build_paper.py`; every number is read from `output/scenarios_2026/*.csv`, so edit text there
  (edits made directly in Word are lost on rebuild).
- New analysis output: `scenarios_2026.R` now also writes `no_treatment_trajectory.csv` (no treatment in
  any stratum 2015–2030). 2030 incidence 130.3% (121.6–140.0) and HCV deaths 101.2% (100.6–102.0) of 2015.
  Scenario results are unchanged (same sampling seed).
- Core argument: J treatment from 2019 cuts incidence 75.9% by 2030 but deaths only 42.6%; in 2026 X holds
  94% of released chronic infections and F/X account for 97% of HCV deaths; the F pool (~1,560 eligible/yr)
  saturates, so extra F quota is wasted while mortality gains come from X. Allocation rule: give F enough
  tests to cover its eligible pool, the rest to X.
- Placeholders to fill: affiliation, corresponding email, acknowledgments, funding, competing interests,
  author contributions, repository link. Verify the WHO target wording (80%/65% vs 2015) against the cited
  WHO documents before submission. PLOS formatting extras (double spacing, line numbers) are not applied.

## Learning route (2026-10-06)
- `route.md` (repo root): staged learning path and full reproduction guide, with commands, code examples,
  `file:line` locations, inputs, outputs, runtimes and expected values for every stage. Examples 1 and the
  treatment exercise were executed and verified.
- Fixed in S1 Text: `pi_recid` is 0.59 in `setup.R:154` (the value the model runs with), not 0.65 as in the
  technical PDF and the `setup.R:34` comment. Arrest/release rates in `setup.R:150-153` are labelled "GUESS";
  confirm their provenance with the author.

## Recalibration + NPE vs MCMC benchmark (2026-10-10/11) — superseded by the cleanup below; `benchmark/` was removed
- Equilibrium check is now 10 years (user commit `5e3bf80`). Old outputs still pass it (checked), but were superseded.
- Nelder–Mead rerun: `output/calibration/r_20261010_agepop` (8 parallel starts, `--sd-perturb 0.9`, maxit 3000,
  `--cores`). Best objective 19.01 (`perturb_1`, converged), RMSE 0.0044, MAPE 2.3%. It is only the prior centre.
- Benchmark per `docs/NPE_vs_MCMC_benchmark.md` in `benchmark/npe_vs_mcmc/` (read `README.md`, `report.md`,
  `model_audit.md`): 3 replicates × {2-round SNPE-C/APT, 1-round NPE, budgeted adaptive RWM MCMC} at ~26,000
  forward runs / ~1,800 s on 4 cores each, plus a 200,000-evaluation reference MCMC. Policy runs: 1,000 states each.
- Verdict: comparable at matched cost on policy outcomes (all within 1.2% of reference; NPE/MCMC within 0.7%).
  Single-round NPE has the best means but 95% intervals ~4.5% too narrow (systematic); SNPE worst (47% of cells
  pass); budgeted MCMC unbiased but noisy (R-hat up to 1.07). **MCMC chosen**; production posterior =
  `output/calibration/mcmc_reference/mcmc_draws.csv` (R-hat ≤ 1.002, bulk ESS ≥ 3,433, MCSE ≤ 0.017 SD).
- Code changes: `npe_train.py` rewritten (fixed biased `force_first_round_loss=True` in round 2; APT atomic loss;
  MAF 50×5; CPU); `run_npe.R` (equilibrium filter, no padding of failed sims, timings, SBC optional);
  `mcmc.R` (`log_posterior_exact`, `arwm_budget`); `run_mcmc_validation.R` rewritten (prior-started chains, budget,
  deadline, posterior-package diagnostics); `run_calibration.R` (`--cores`, perturb SD 0.9);
  `scenarios_2026.R` (`--posterior-files`, `--no-figures`, `--min-valid-share`, `scenario_2030_draws.csv`, UI error
  bars/ribbons in figures). Commands: `docs/training.md`. Tests: 6 files, 63 checks pass.
- Scenarios rerun from 2,000 reference-MCMC draws (`output/scenarios_2026/`, README table updated):
  J-only incidence −75.6% (68.5–80.4), mortality −42.6% (41.6–43.7). **5000/yr F25% now FAILS at the median**
  (mortality 64.9%, 63.2–66.5; 45.4% of draws pass). Median PASS: 8000/yr F25% (100% of draws), 8000/yr F50% (82.0%),
  no cap. No-treatment 2030 incidence 132.5% (120.9–147.3) of 2015.
- Updated 2026-10-11: manuscript rebuilt (`python3 scripts/build_paper.py .`); methods now describe Nelder–Mead prior centre +
  adaptive Metropolis MCMC + the benchmark; pass lists, capacity multiples, drug share and 2015 populations are derived from
  outputs (2/15 finite scenarios pass). SI Sections 2.3–2.6 replaced (MCMC diagnostics Table S3, benchmark Table S4, posterior
  Table S5, fit Table S6); old seed table removed, later tables renumbered S7–S12. New refs `haario2001`, `vehtari2021`.
  Calibration figure `figures/calibration_agepop.png` + `figures/calibration_interval_overlap.csv` regenerated (targets vs
  reference MCMC vs single-round NPE; all 12 overlap). 2015 stratum predictions from 4,000 reference draws are in
  `output/calibration/mcmc_reference/model_prediction_2015_strata.*` (ever-PWID 21,224, 19,226–23,437).
  `plot_calibration.R` now takes `--mcmc-file/--npe-file`; `predict_2015_posterior.R` takes `--posterior-files/--max-draws`.
  Incidence threshold is 80% (user confirmed 2026-10-11; the earlier 90% note in `review.md`/`codex.md` is withdrawn).

## Inference cleanup + LHS NPE (2026-10-11, `docs/calibration_update.md`)
- Kept only **single-round NPE** and the **reference MCMC**. Removed multi-round APT/SNPE (round-2/proposal/atomic-loss
  code in `npe_train.py`/`run_npe.R`), short-chain MCMC (deadline/time-matched logic; dead `am_mcmc`, `make_priors`,
  `log_prior`, `log_posterior`, `predictive_intervals` in `mcmc.R`), the `benchmark/` folder (driver, config, report,
  smoke tests, SNPE/short-MCMC runs and policy runs) and their sbi logs. `run_mcmc_validation.R` → `run_mcmc_reference.R`
  (defaults reproduce the production run). Reference MCMC results moved unchanged, not rerun.
- `output/calibration/` now: `r_20261010_agepop/` (prior centre), `mcmc_reference/` (production posterior + 2015
  predictions), `npe_lhs/` (new NPE), `npe_vs_mcmc/` (comparison, `report.md`, `model_audit.md`, `policy/{reference,npe_lhs}`,
  `previous_random_npe/rep1..3` compact record). Deleted superseded `npe_agepop*`, `mcmc_recal`, `r_20260929_agepop*`,
  two aborted timestamped runs.
- LHS NPE: `scripts/calibration/prior_sampling.R` (`lhs_unit` + Normal inverse CDF per parameter; prior is independent
  Normal, no truncation). 50,000 attempted, 49,765 valid, 235 equilibrium failures, 0 non-finite (not replaced).
  Seeds: design 20261012, training 20261013. Same network/training (MAF 50×5, batch 256, lr 5e-4, patience 20, ≤250 epochs).
  Trained 71 epochs in 37.5 s; gate PASS.
- Comparison (`scripts/compare_npe_mcmc.R`): all means within 0.1 SD; 95% widths 0.94–0.99 (mean 0.961); 11/12 parameters
  pass (θ10 lower endpoint 0.26 SD); policy medians within 0.4%; joint-pass 44.5% vs 43.5% (within MC error). Improvement
  over the earlier 26,000-random runs is small and confounded with sample size.
- Every project script now has a short header in the `tests/helper.R` style. Calibration figure regenerated
  (MCMC vs NPE-LHS). Manuscript rebuilt (SI §2.5/Table S4 now NPE-LHS vs MCMC). Tests pass.
- Flags: protected `docs/AGENTS.md` still says SNPE is primary and MCMC a fallback; `docs/NPE_vs_MCMC_benchmark.md`
  (user's protocol) kept; `problems.md`/`回答.md` are historical logs that name deleted folders;
  `tests/integration/.ipynb_checkpoints/` is a generated copy; `output/*_rate_fit.rds` are not read by the model.

## Reference-MCMC ablation (2026-10-12/13, `docs/ablation_check.md`)
- The plan's §3 baseline (NPE-initialised `am_mcmc`, monthly step, coda R-hat) described already-removed code;
  the current reference MCMC already used the daily step and rank-normalised R-hat / bulk-tail ESS / MCSE.
- Code (commits `094c8af`, `8ce1cb9`): `arwm_budget` records adaptation drift (relative Frobenius change),
  condition numbers, jump sizes, CPU time and accepts a fixed start; `run_mcmc_reference.R` gains
  `--adapt-start/--adapt-every/--thin/--start-file/--cov-file/--cov-scale/--save-unthinned` and writes trace plots,
  per-chain summaries, drift/eigen/jump CSVs and autocorrelation times. **New defaults:** warmup 10,000,
  adapt-start 5,000, thin 10 (storage only). Stored production run = `--n-warmup 4999 --adapt-start 200 --thin 1`.
- Driver `scripts/ablation_check.R` → `output/ablation_check/ab_20261012/` (report: `ablation_report.md`).
  Arms (same target, daily dt, 4 chains × 30,000 evaluations): A_mean/A_disp/A_tail (NPE covariance; mean,
  dispersed, tail starts), B_pilot (pilot covariance), C_map (inverse Hessian at fit best; 314 evals), D_prior
  (production prior-scaled covariance, random NPE starts); plus dt_daily vs dt_monthly short chains and a
  500-draw deterministic step comparison.
- Results: all arms R-hat ≤ 1.005, bulk ESS 1,422–1,846, tail ESS ≥ 2,319, acceptance 0.24–0.28, starts
  irrelevant (all within 1.9 combined MCSE of production). Posterior-shaped proposals give ~5–30% more ESS than
  prior-scaled; NPE vs pilot vs Hessian not distinguishable (±10% ESS noise, one run each). All condition numbers
  ~100 (no near-singularity). NPE SD/MCMC SD 0.93–0.99 (slightly narrow, not collapsed). Predeclared drift
  tolerance 0.01 **not met** (last update change 0.011–0.021) although the kernel is frozen before retention.
  Thinning 1000 fails (tail ESS 91); thin 10 keeps ≥95% ESS. Monthly vs daily: prevalence ≤0.18 pp, populations
  <0.7%, posterior means within 0.11 SD; monthly is ~18× cheaper.
- Suggested follow-ups (not run): scale the inverse Hessian by 2.38²/d; refresh adaptation every 500 proposals or
  use an ESS-scaled drift tolerance; replicate arms to rank NPE/pilot/Hessian.
