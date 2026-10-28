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

## Recalibration + NPE vs MCMC benchmark (2026-10-10/11) — superseded; see "Final review" at the end
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

## Inference cleanup + LHS NPE (2026-10-11) — superseded; see "Final review" at the end
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

## Reference-MCMC ablation (2026-10-12/13) — removed per `docs/final_review.md`; code, outputs and doc deleted
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

## Final review (2026-10-13/14, `docs/final_review.md`) — superseded; see "First-arrest handoff" below
- Production path: Nelder–Mead (`output/calibration/r_20261013_agepop`, objective 22.58) → `prior_centre()` (moves
  any coordinate to ≥5 prior SD inside the bounds; only θ6, the 60+ contact scale, moves) → reference MCMC
  (`output/calibration/mcmc_reference`, defaults warmup 10,000 / adapt-start 5,000 / thin 10; R-hat ≤ 1.002, bulk
  ESS ≥ 2,824) → scenarios (`output/scenarios_2026`, 2,000 draws). Commands in `docs/training.md`.
- Cross-validation only: `npe_lhs` (50,000 LHS, all valid, 5-network MAF ensemble; SBC n=1,000 gate PASS),
  `npe_vs_mcmc` (12/12 parameters pass, policy within 0.6%), `npe_vs_mcmc/mcmc_sbc` (8 datasets).
- Removed: ablation, MAP+Hessian, pilot, multi-round NPE, `--sampling`, `--start-file/--cov-file/--cov-scale/
  --save-unthinned`, `previous_random_npe/`, `r_20261010_agepop/`, `probe_beta_constraint.R`, generated binaries.
- Model change in user commit `7c3dad3` (`lambda1 = lambda3*c_composite/beta_scale`, `psi_DC` 0.5) moved results:
  J-only incidence −59.2% (46.4–69.0), mortality −41.8%; only 8,000/yr F25% (100% of draws) and no-cap pass.
- `scripts/model/update_rate_params.R` copies Stan rate-fit summaries into `setup.R` (idempotent on stored fits;
  refitting needs `data/sps.csv`, not in the repo).
- Flags: 60+ seroprevalence misfit is structural (0.508 vs 0.429); θ6/θ12 prior-dominated; the two
  `preprocess.R` copies differ scientifically (left as is); long runs must be launched from a copy of the
  script (Rscript reads lazily).

## First-arrest handoff + viraemic mortality (2026-10-14/15, `docs/first_arrest_mortality_update.md`) — superseded; see "lambda1 = lambda3" at the end
- **Rates:**
  - `scripts/model/preprocess.R` is the single preprocessing script; the top-level copies were removed.
  - Release, re-arrest and `pi` (0.720) come from all consumption records.
  - The 2015 prison counts and first offenders come from the screened subset, which the user
    defined as the PWID. The counts equal the targets.
  - Stan's prison steady state uses ageing 0.1 and `mu*omega`. It gives `f_first` (32, 304, 102,
    75.9, 66.1, 19.6 per year), and `update_rate_params.R` writes it.
- **Model:**
  - `sim.cpp` derives `lambda1` from `f_first` through the HCV-free D steady state
    (`first_arrest_rate()`).
  - Chronic viraemic people (state c) have non-liver mortality × `hr_viraemic` = 1.47
    (Lee 2012, Table 3).
  - `omega` is rescaled with the frozen `viraemic_share_ref` = 0.205 (fixed point).
- **Calibration:**
  - θ7–θ12 are log excess-entry scales (`entry_scales()`), so every θ is feasible.
  - New likelihood term `nll_entry`: the 2015 entry against the CNB counts, log-normal SD 1.0 (user
    decision). It is needed because D was unidentified (38k–1.7M PWID at the same objective).
  - Nelder–Mead uses `--nm-restarts 5`. Result: objective 21.44; 60+ seroprevalence 0.455 against
    0.429 (0.65 SE); 2015 PWID 18,314.
  - HR sensitivity (`output/calibration/hr_sensitivity/`): 60+ is 0.494 at HR 1 and 0.439 at HR 1.77.
- **Cross-validation:** NPE has 48,918 of 50,000 valid (1,082 equilibrium failures), with 18
  summaries (adds the 2015 entry). SBC n = 973: cov95 0.936–0.958, all ranks uniform; gate PASS.
- **Reference MCMC:** R-hat ≤ 1.003, bulk ESS ≥ 2,319, 4.5 h. Eight parameters are
  prior-dominated (θ5, θ6, θ7–θ12).
- **2015 PWID:** ever-PWID 18,453 (17,799–19,380). All 12 predictive intervals overlap the
  targets. NPE agrees: 12/12 parameters, policy within 1.2%. MCMC SBC (8 datasets): cov95 0.93.
- **Scenarios (2,000 draws):**
  - J only: incidence −47.6% (36.5–57.5) and mortality −47.4% (44.7–48.8).
  - No scenario meets the incidence target. No cap gives −57.3%, because 59% of 2030 infections
    are in D (never imprisoned).
  - Mortality is met at the median only at 8,000/yr F25% (67.2%; 87% of draws), F50% (65.4%;
    60%) and no cap.
  - The joint-pass probability is 0 everywhere.
- **Manuscript:** the narrative was rewritten for these results. `build_paper.py` asserts the
  claims it relies on and stops if the results change.
- **Flags:**
  - The 2015 prison counts both inform `f` and are targets.
  - The HR comes from a general-population cohort.
  - The first-arrest flow is realised within 2.7% (HCV deaths in D).
  - The NPE cross-validation counts equilibrium failures as exclusions.

## Figures redrawn with the scientific-visualization skill (2026-10-19)
- `scripts/figures.py` redraws all 21 figures from saved outputs (no re-simulation) with the user-level skill
  `~/.claude/skills/scientific-visualization` (style "default", PLOS widths 190.5/132 mm, PNG 300 dpi + editable SVG).
  Needs Python ≥ 3.11 with `matplotlib==3.11.2 pillow==12.3.0 pandas scipy`. Run it last (`docs/training.md` §7).
- `scenarios_2026.R` now also writes `scenario_stratum_totals.csv` and `scenario_item_totals.csv` (cost figures).
  Rerun with the production seed reproduced every other scenario CSV byte for byte.
- The old ggplot/base/matplotlib plotting code in the R scripts and `npe_train.py` was left in place (removing it was
  blocked by the auto-mode classifier); it writes the same PNG names, so running it after `figures.py` reverts them.
- Data notes surfaced while redrawing: F/X treatment quota equals actual starts in every year (dashed quota lines are
  hidden; capacity-utilisation figure is 100% everywhere); no scenario passes jointly, so Fig S6 shows only J only.

## lambda1 = lambda3, no entry anchor (2026-10-20/21, `docs/lambda1_revert.md`) — superseded; see "First-arrest multiplier" at the end
- **Model (user decision):** `lambda1` is hard-coded in `setup.R` equal to `lambda3` (0.546, 0.488, 0.472,
  0.472, 0.386, 0.266). `f_first`, `first_arrest_rate()`, `entry_scales()`, `excess_theta()`, `nll_entry` and
  the targets' `entry` are gone. θ7–θ12 are plain log entry scales on the CNB 2015 counts (bounds 0.01–10⁶).
  Viraemic HR 1.47 kept; fixed point 0.205 → 0.196 moves `omega_0` by 0.38% (< 0.5%), so no refit.
- **Fit** `output/calibration/r_20261020_agepop` (8 starts): objective 18.73, every target fitted almost exactly
  (60+ 0.4286). Prior centre = best fit (nothing near a bound). HR sensitivity: 60+ 0.440 (HR 1) to 0.429 (1.77).
- **MCMC:** R-hat ≤ 1.002, bulk ESS ≥ 3,053, 4.7 h. Prior-dominated: θ6, θ9, θ10, θ12. 2015 ever-PWID 17,226
  (15,649–18,968); D only 1,356.
- **Cross-validation:** NPE 49,278/50,000 valid, SBC (980) cov95 0.944–0.961, gate PASS; 11/12 parameters pass
  (θ4 endpoint 0.22 SD), policy within 0.3%. MCMC SBC (8) cov95 0.95.
- **Scenarios (2,000 draws):** J only incidence −81.1% (76.5–84.3; 71% of draws pass), mortality −50.5%.
  Every finite scenario meets incidence. Mortality and joint pass (≥ 99.9% of draws): 5,000/yr F25%/F50%,
  8,000/yr F25%/F50%, no cap. 8,000/yr F75% 64.9% (40%). F eligible pool ~2,330/yr.
- **Manuscript** rewritten for these results (`build_paper.py` asserts them). Main caveat now stated: the
  result rests on `lambda1 = lambda3` (the records cannot identify it); a lower first-arrest rate would leave
  more injectors outside prison.
- `docs/first_arrest_mortality_update.md` and `docs/final_review.md` are historical.

## First-arrest multiplier, seroprevalence-only fit, audit fixes (2026-10-21/22, `docs/first_arrest_multiplier.md`) — superseded; see "Redesign" at the end
- **Rates:** `rearrest_rate.stan` is survival-only (no prison steady state, no `f_pop`); `pi` 0.724,
  `lambda3` 0.548, 0.478, 0.469, 0.471, 0.384, 0.265. Both fits saved as `list(summary, draws)`;
  `lambda2`/`lambda3` = posterior mean of `exp(log_lambda)`. `preprocess.R` is vectorised (bit-identical
  Stan inputs) and writes `output/prison_2015_counts.csv` (aggregates; 3,367 in prison).
- **Model:** first arrest = `c_first` × `lambda3`; CNB entry unscaled; background mortality
  `mu` × SMR by state (9.52 viraemic a/c/t, 4.73 u/r; van Santen 2014, antibody-based). `y0` now
  starts in the right compartments (was shifted by one).
- **Calibration (user decisions):** 7 parameters (6 contact scales + log `c`, `c` 0.5–100);
  seroprevalence only, Beta(p·φ, (1−p)·φ), φ = 100, no counts; fixed prior contact N(0, 1.5),
  log `c` N(0, 0.5), one density `nlp()` for Nelder–Mead, MCMC and NPE; uninformative NM starts.
  Discarded attempts: CNB entry with the population target (c = 0.034, objective 119); binomial
  with screened counts (rejected by the user).
- **Fit** `output/calibration/r_20261021_agepop`: all 8 starts agree, RMSE 0.013, `c` 1.10.
  SMR-bound fits in `smr_sensitivity/` (`c` 1.097–1.100).
- **MCMC** (adaptation on the recent half of warmup): R-hat ≤ 1.002, bulk ESS ≥ 3,250, acceptance
  0.18–0.22, 4.0 h; `c` 1.24 (0.57–3.13), prior-dominated (as are θ5, θ6). 2015 J 8,451 vs observed
  3,367 (not fitted); ever-PWID 50,218.
- **Cross-validation:** NPE SBC gate PASS, but 33% of NPE draws fail equilibrium (dropped) and 0/7
  parameters pass (intervals 28% wider); policy reductions within 0.3%, joint-pass 8 pp lower.
  MCMC SBC cov95 0.96. NPE is not validated under this prior.
- **Outcomes (audit fixes):** HCV deaths = DC/HCC liver deaths + viraemic non-liver excess
  `mu(9.52 − 4.73)`; prison tests follow the F/X rule (antibody for entrants u/a/c/r, RNA for a/c/r).
- **Projection entry (2026-10-23):** entry into D is held at the CNB 2015 counts in every year
  (the 2016–2025 schedule was removed from `setup.R`/`sim.cpp`); calibration is unaffected.
- **Scenarios:** J only incidence −84.5% (68.7–90.5; 77% of draws), mortality −54.9%. All finite
  scenarios meet incidence; mortality at the median at 5,000/8,000 per year with F 25%/50% (joint
  95%, 51%, 98%, 99%); F 75% misses. No cap 87.1%. More tests to X always lower deaths more.
  No-treatment 2030 = 100% of 2015.
- **Code:** `scenarios_2026.R` draws its old ggplot figures only with `--figures true`; prediction,
  plot and scenario scripts stop if the fit does not match `setup.R`.
- **Manuscript** rewritten; `build_paper.py` asserts every result claim it makes.

## Redesign: entry scales, first-arrest and community targets (2026-10-27, `docs/redesign_prompt.md`) — current state
- **Rates:** rerun, unchanged; `setup.R` not modified.
- **Model fixes:** ψ mortality fix in `sim.cpp`; `split_by_age` fix in `preprocess.R`.
- **Calibration (user decisions):** 13 parameters: θ1–θ6 contact N(0, 1.5); θ7–θ12 log entry scales on the CNB 2015
  counts, N⁺(0, 1), scale 1–100; θ13 log `c` N(0, 0.5), 0.5–100. Likelihood: Beta seroprevalence (φ 100) + negative
  binomial first arrests vs CNB (257, 644, 292, 94, 21, 1; size 50) + log-normal D + F vs Teo 2019 (11,000,
  6,500–17,000; SD 0.2453). Community = D + F.
- **Fit** `output/calibration/r_20261027_agepop` (maxit 5000): 7/8 starts agree (objective 18.8–19.5; one local optimum
  121.5); RMSE 0.014, `c` 2.28, D + F 7,720. Entry scales 30+ at the bound 1; first arrests over target at 30+.
- **MCMC** (8 chains × 60,000, 1.9 h): R-hat ≤ 1.004, bulk ESS ≥ 2,158, acceptance 0.12–0.14. `c` 2.78 (1.55–5.49).
  Prior-dominated (SD ratio > 0.7): θ5, θ6, θ13. 2015: D + F 7,957, J 9,741 (observed 3,367, not fitted), ever-PWID 57,742.
  MCMC SBC (8 datasets, 104 ranks): cov95 0.952, max split R-hat 1.015.
- **NPE:** gate **FAIL** (θ2, θ11 biased, θ10 U-shaped; coverage fine), 24% of draws fail equilibrium, 2/13 parameters
  agree with MCMC. User decision: accept and document; NPE will likely be deleted later. The test now checks the
  marker matches the SBC result instead of requiring a pass.
- **Scenarios (2,000 draws):** J only incidence −84.4% (67.9–90.4), mortality −54.0%. All finite scenarios meet
  incidence; mortality at the median only at S5000_F25 (66.2%, joint 80%), S8000_F25 (71.6%, 99.85%), S8000_F50 (67.4%).
  X holds 89% of released chronic infections. Manuscript rebuilt; `build_paper.py` asserts these claims.
- Removed: superseded `output/calibration/r_20261021_agepop/`. Next per `docs/instructions.md`: ponytail audit →
  refactor → review, one commit per step, no push.

## Ponytail pass, DC fix and new manuscript (2026-10-28)
- **Ponytail** (`docs/instructions.md`): audit → refactor (`equilibrium_data()`, `rate_components()` moved to
  `model_metrics.R`, fit.rds without NM history, figure grids 3×5, read_targets everywhere) → review fix.
- **Model bug fixed** (`sim.cpp`): recovered people in DC were also added back to DC when progressing to HCC
  (people created). Found by the new R-vs-C++ HCV-death test. Calibration unchanged (no r in DC before 2019);
  scenarios rerun: mortality +0.4–0.5 pp (J only 54.4%), same pass set, S5000_F25 joint 87%.
- **Manuscript rewritten from scratch** with the scientific-writing skill (PLOS Medicine): new
  `scripts/build_paper.py` + `scripts/paper/`. Registries in `docs/paper/` (276 claims, 85 sources, all
  unverified). Placeholders: authors/affiliation, ethics, funding, conflicts, CRediT, data availability,
  AI statement, registration. Open items in `docs/paper/unresolved_issues.md` (seroprevalence targets differ
  from current prison processing by up to 0.5 pp; death-rate and contact-matrix sources; CHEERS wording).
- NPE: gate FAIL accepted (user); user plans to delete NPE later (`sbc_mcmc.R` depends on NPE SBC files).


## Yearly counts for the WHO outcomes (2026-10-28)
- Trajectories and 2030 results are now events **summed over each calendar year** (`annual_events()` in
  `model_metrics.R`, trapezoid rule on daily output), not the 31 December rate, which overstated 2030 reductions.
  2015 = one year at the equilibrium rate (2016 count within 0.1%, checked in `scenarios_2026.R`).
- Rerun: reductions about 1.3 pp lower. J only incidence 83.2% (66.5–89.4), mortality 53.1%. Same median pass set,
  but S5000_F25 is narrow (65.4%, joint 60%); S8000_F50 66.6% (87%); S8000_F25 70.8% (99.85%).
- Trajectory figures: yearly counts at mid-year joined by a monotone PCHIP curve, ticks at 1 January.
- Manuscript, S1 Text and the scenarios README updated; `build_paper.py` asserts S5000_F25 stays narrow.
