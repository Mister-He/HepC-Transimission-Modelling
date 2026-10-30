# PWID share of the inflow, and treatment from 2026 as the main analysis

Status: done 2026-10-30 (Section 6 holds the user decisions; they override Sections 2–3). Results: Section 7.

## 0. Goal and boundary

- Two changes, then rerun the whole chain: calibration → 2015 predictions → calibration figure →
  scenarios → figures → manuscript.
- **Change only what Sections 2–3 specify.** Everything else stays as at HEAD `f20eed1`:
  - compartments D/J/F/X, six age groups, HCV states u/a/c/t/r, liver stages;
  - rate fits (`lambda2`, `lambda3`, `pi`), `sim.cpp` logic, mortality, HCV-death definition;
  - 13 parameters, priors, likelihood forms, seroprevalence and community targets;
  - 150-year equilibrium to 2015, RK4 `dt = 1/365`, yearly-count outcomes (`annual_events()`);
  - F/X screening quotas, costs, discounting, WHO targets (80% incidence, 65% mortality vs 2015).
- Constraints: no push; commits per `commit_rule.md`; `data/sps.csv` aggregates only, never
  sent to Ollama; long runs with `nohup` from frozen script copies; old R plotting code kept.

## 1. Source: SingStat M891521 (CNB, "Number of drug abusers by main drug abused, annual")

- Downloaded 2026-10-30 from `https://tablebuilder.singstat.gov.sg/api/table/tabledata/M891521`.
- 2015 new drug abusers: total 1,309; heroin 113; methamphetamine 1,014.
- PWID share = (113 + 1,014) / 1,309 = **0.8610**.
- The table has no age split, so the same share applies to every age group.

| Age | CNB 2015 | × 0.8610 |
|---|---|---|
| <20 | 257 | 221.3 |
| 20–29 | 644 | 554.5 |
| 30–39 | 292 | 251.4 |
| 40–49 | 94 | 80.9 |
| 50–59 | 21 | 18.1 |
| 60+ | 1 | 0.86 |

## 2. Change 1 — inflow restricted to heroin and methamphetamine users

1. Store the 2015 drug counts in `data/singstat_new_drug_abusers_by_drug.csv` (aggregate, public).
2. `setup.R` `beta` = CNB 2015 × 0.8610 (entry into D before the entry scales θ7–θ12).
3. `data/calibration_targets.json` `first_arrests` = CNB 2015 × 0.8610.
4. Recalibrate: Nelder–Mead (8 starts) → prior centre → reference MCMC (8 × 60,000, ~2 h) →
   2015 predictions → calibration figure.
5. Update `data/README.md`, `setup.R` comments, tests that hard-code 257…1.

## 3. Change 2 — main analysis: no treatment before 2026

### 3.1 Main analysis
- 2015–2025: no treatment in any stratum. The model stays at the 2015 equilibrium, so the 2026
  state equals the 2015 state.
- 2026–2030: scenarios start; prison (J) treatment coverage becomes a policy lever (Section 6, Q4).

### 3.2 Sensitivity analysis
- The current design: J treatment from 2019 (annual start probability 0.99995), F/X screening
  from 2026. Same 17 scenarios as now, reported in S1 Text.

### 3.3 Code
- `scenarios_2026.R`: a `--j-start-year` option (main 2026, sensitivity 2019) and an output
  directory per analysis (`output/scenarios_2026/` main, `output/scenarios_2026_j2019/`
  sensitivity).
- `data/scenario_policy_config.json`: `j_start_year` 2026 (main).
- `figures.py`, `build_paper.py`, `scripts/paper/*`: main results from the main run, sensitivity
  table/figure in S1 Text; asserts rewritten for the new results.

## 4. Verification
- `Rscript scripts/run_tests.R --scope all` passes.
- Nelder–Mead starts agree; MCMC R-hat ≤ 1.01, bulk ESS ≥ 400.
- Main run: the 2016–2025 counts equal 2015 within 1e-3 (no treatment, equilibrium).
- Sensitivity run with the new posterior uses the same code path as the old main run.

## 5. Then follow `docs/instructions.md`
- ponytail audit → refactor → review and fixes → tests; one commit per step, no push.

## 6. User decisions (2026-10-30)
1. The 0.861 share applies to both `beta` and the first-arrest target, same share at every age.
2. The target is rounded (221, 554, 251, 81, 18, 1); `beta` is not rounded.
3. Rerun Nelder–Mead, reference MCMC and the MCMC SBC (8 datasets, own prior draws via
   `prior_sampling.R`, output `output/calibration/mcmc_sbc/`).
4. Main analysis has both:
   - (a) J treatment from 2026 at the current coverage (0.99995/yr) in all 17 policies, plus a
     "no treatment" reference row;
   - (b) J coverage grid 0%, 50%, 90%, 99.995% × all 17 policies (68 scenarios), main-text table
     and heatmap.
5. Sensitivity: J from 2019 (current design) on the new posterior, 2,000 draws, S1 Text table +
   trajectory figure, `output/scenarios_2026_j2019/`.
6. Delete NPE and NPE SBC (code, tests, outputs, docs, manuscript text) in one commit before the
   recalibration.

## 7. Results (2026-10-30)
- Fit `output/calibration/r_20261030_agepop`: 7/8 starts agree (objective 18.77–18.84), `c` 2.11,
  D + F 6,786. MCMC: R-hat ≤ 1.004, bulk ESS ≥ 2,137, 1.3 h. MCMC SBC (8): cov95 0.92.
- Main run: 2016–2025 within 0.05% of 2015 (check passed). J only −49.8% incidence, −32.6% deaths.
  Five finite policies meet incidence; none meets mortality; only no cap meets both.
- J coverage grid (`output/scenarios_2026_jcov/`): every reduction rises with coverage; no finite
  policy meets mortality at any coverage.
- Sensitivity (`output/scenarios_2026_j2019/`): J only −83.4% / −53.2%; mortality met at
  5,000/yr F25%, 8,000/yr F25% and F50% (as in the previous main analysis).
- Manuscript: main text Table 2 and Fig 5 (coverage), S1 Text Section 9 (Table J, Fig F).
