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
