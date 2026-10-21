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

## First-arrest multiplier, seroprevalence-only fit (2026-10-21/22, `docs/first_arrest_multiplier.md`) — current state
- **Rates:** `rearrest_rate.stan` is survival-only (no prison steady state, no `f_pop`); refit
  `pi` 0.724, `lambda3` 0.548, 0.478, 0.469, 0.471, 0.384, 0.265; release unchanged.
- **Model:** first arrest = `c_first` × `lambda3`; CNB entry unscaled; background mortality
  `mu` × SMR by state (9.52 viraemic a/c/t, 4.73 u/r; van Santen 2014 Table 3, Amsterdam,
  antibody-based); no `omega`, `hr_viraemic` or `viraemic_share_ref`.
- **Calibration (user decisions):** 7 parameters (6 contact scales + log `c`, `c` 0.5–100);
  seroprevalence only, Beta(p·φ, (1−p)·φ) with φ = 100, no counts (`population` removed from
  the targets JSON); fixed prior contact N(0, 1.5), log `c` N(0, 0.5); Nelder–Mead = posterior
  mode from uninformative starts. Two earlier attempts were discarded: CNB entry with the
  population target fit at c = 0.034 (objective 119; prison 20–29 381 vs 1,001), then a binomial
  with screened counts was rejected by the user.
- **Fit** `output/calibration/r_20261021_agepop`: all 8 starts agree, RMSE 0.013, `c` 1.10.
  SMR-bound fits in `smr_sensitivity/` (`c` 1.097–1.100).
- **MCMC:** R-hat ≤ 1.005, bulk ESS ≥ 1,798, acceptance 0.08–0.17, 4.0 h; `c` 1.21 (0.55–3.15),
  SD ratio 0.90 (prior-dominated, as are θ5, θ6). 2015 J 8,451 vs observed 3,367 (not fitted);
  ever-PWID 50,223.
- **Cross-validation:** NPE SBC gate PASS but NPE vs MCMC 0/7 parameters pass (intervals 23%
  wider); policy reductions within 2.4%. MCMC SBC cov95 0.96.
- **Scenarios:** J only incidence −85.6% (68.9–91.3), mortality −45.8%; all finite scenarios meet
  incidence, none meet mortality (best 8,000/yr F25% 60.5%); only no cap meets both (99.8%).
  More tests to X always lower deaths more (F eligible ~5,600/yr, 6% chronic; X 9%).
- **Manuscript** rewritten for these results; `build_paper.py` asserts them.

