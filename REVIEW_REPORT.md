# Independent adversarial review — Singapore PWID HCV transmission model

Reviewer: Claude (Opus 5.5) acting as an independent reviewer, 2026-11-03.
Reviewed commit: `5fc2b95` on `feature/fixinflow`, plus the uncommitted manuscript edits that were in the
working tree at the start of the review (`docs/paper/*`, `scripts/paper/*`, `scripts/build_paper.py`).
The reviewer did not edit or rebuild those manuscript files.

Status legend: **OPEN** (not yet addressed) · **FIXED** (changed and verified) · **DOCUMENTED** (needs an author
or research decision; evidence and options recorded here) · **PARTIAL**.

---

## 0. What the project claims

- **Question.** Can prison HCV treatment (stratum J) from 2026, alone or with capped community screening of
  released people (F = still injecting, X = stopped), meet the WHO 2030 targets: an 80% fall in incidence and a
  65% fall in HCV deaths versus 2015?
- **Model.** Deterministic ODE with 480 compartments: 4 strata (D, J, F, X) × 4 liver stages × 5 HCV states
  (u, a, c, t, r) × 6 age groups. The C++ solver (`scripts/model/sim.cpp`) uses daily RK4 with
  operator-split ageing.
- **Calibration.**
  - 13 log-scale parameters: 6 contact scales, 6 entry scales and the first-arrest multiplier `c`.
  - Likelihood: 2015 prison seroprevalence (Beta, φ = 100), CNB first arrests (negative binomial, size 50)
    and the Teo 2019 community-PWID estimate (log-normal).
  - A Nelder–Mead posterior mode, then an adaptive random-walk Metropolis "reference MCMC".
- **Scenarios.**
  - 17 policies (J only, 15 capped F/X screening policies, no cap).
  - Each policy is run over 2,000 posterior draws and reported as the median with a 95% UI.
- **Headline results** (`output/scenarios_2026/scenario_2030_summary.csv`):
  - J only: incidence −49.8% (37.7–56.9) and HCV deaths −32.6% (30.6–34.8).
  - Five finite policies meet the incidence target.
  - No finite policy meets the mortality target. The best is 8,000/yr with 25% to F (62.6%, 59.0–65.6;
    joint pass 5.7% of draws).
  - Only the uncapped policy meets both targets.

## 1. What was executed (evidence base)

All commands were run on macOS (Darwin 25.6), R 4.5.1, Rcpp 1.1.0, RcppArmadillo 14.6.3.1, posterior 1.6.1,
cmdstan 2.39.0, 16 cores. "Fresh clone" means `git clone` of commit `5fc2b95` into the session scratchpad.

| # | Check | Result |
|---|---|---|
| E1 | `Rscript scripts/run_tests.R --scope all` (working tree) | 6 files, 105 checks, 0 failures |
| E2 | **Independent re-implementation of the ODE right-hand side** (pure R, written from the model contract, not from `sim.cpp`), one RK4 day from a random state with treatment, age limits and the RBV-free DC branch switched on | max relative difference from `sim.cpp` **1.8e-16** (identical to machine precision) |
| E3 | Mass balance of the independent RHS: d(total)/dt − (entry − deaths) | −1.0e-12 |
| E4 | Numerical convergence at the best fit, 150-year burn-in: `sim.cpp` dt = 1/365 vs dt = 1/1460 vs `deSolve::lsoda` (rtol 1e-10) on the independent RHS with continuous ageing | all 14 calibration summaries agree within 4.3e-4 relative (most < 2e-4) |
| E5 | `event_rates()` HCV deaths vs an independent liver + excess decomposition | identical (245.5079) |
| E6 | Scenario pipeline from a fresh clone, `scenarios_2026.R --draws-per-seed 2000 --cores 8 --sample-seed 20261005` (~16 min) | **all 10 result CSVs and `draw_validity.txt` byte-identical** to the committed outputs; `calibration_source.txt` differs only in absolute paths (R-16) |
| E7 | Headline numbers recomputed from `scenario_2030_draws.csv` (medians, 2.5/97.5% quantiles, joint-pass share) | match the summary CSV and the manuscript abstract |
| E8 | MCMC convergence diagnostics recomputed with `posterior` from the stored thinned draws | R-hat 1.001–1.004, bulk ESS 2,124–3,476, tail ESS ≥ 2,230 (matches the claims in `docs/training.md`) |
| E9 | Nelder–Mead calibration rerun from a fresh clone (8 starts, documented command) | **bit-identical** best fit and `solutions.csv` (§5) |
| E10 | Reference MCMC rerun from a fresh clone (documented defaults) | **bit-identical** draws and diagnostics, except the wall-clock column (§5) |
| E11 | Prison counts from `output/cleaned_data.rds` (aggregates only) | 8,725 people held for drug consumption on 31 Dec 2015; 3,367 of them screened |
| E12 | GitHub Actions status of HEAD (`gh run view 41399261724`) | **CI fails**: job "Python NPE syntax check" can't find `scripts/npe_train.py` |
| E13 | Citation check of the 16 fixed-parameter sources (literature search by a sub-agent; see §4) | see §4 |

---

## 2. Prioritised issues

### R-01 · HCV-death outcome is dominated by an SMR-based non-liver excess concentrated in the 60+ group — **High** (scientific validity of the mortality conclusions)
- **Location:**
  - `scripts/calibration/model_metrics.R:86-94` (`event_rates`)
  - `scripts/model/sim.cpp:74-85` (`mortality`)
  - `scripts/model/setup.R:118-125`
- **Problem:** "HCV deaths" are the WHO-target outcome. They are defined as liver deaths in DC/HCC plus
  `mu_i × (9.52 − 4.73)` for every viraemic person. Three problems follow:
  - The SMRs come from an all-cause comparison (van Santen 2014, antibody-based, Amsterdam 1985–2012).
  - They are multiplied by the Singapore 60+ general-population death rate, 21.9/1000, which includes people
    aged 85+. The result is an all-cause death rate of **20.8% per year** for viraemic 60+ PWID.
  - Liver deaths `mu_DC`/`mu_HCC` are then added on top. If the SMR already contains liver deaths, they are
    counted twice.
- **Evidence (E5):**
  - At the best fit, 2015 HCV deaths are 245.5. Of these, **133.3 (54%)** are the viraemic non-liver excess
    and 112.3 are liver deaths.
  - 84.3 of the 133.3 excess deaths are in the 60+ group.
  - Cure moves a person straight from SMR 9.52 to 4.73, so more than half of the "HCV deaths" respond to
    treatment at once, with no time for liver disease to regress.
- **Impact:**
  - The size of the mortality reduction, and so the main conclusion ("no finite policy meets the mortality
    target"), depends on a definition that is not the WHO's. The WHO counts deaths from HCV-related cirrhosis
    and liver cancer.
  - The direction of the bias on the pass/fail verdict is not obvious. It has to be quantified.
- **Recommended fix:**
  - Report a liver-only HCV-death outcome (DC/HCC liver deaths) next to the current definition.
  - State in the manuscript which definition the WHO target uses.
  - Choosing the primary definition is an author decision.
- **Verification:** the scenarios output liver-only death counts and 2030 reductions, with UIs, for every
  policy. Their sum with the excess reproduces the current mortality column exactly.

### R-02 · Inflow held at the 2015 level for 2016–2030, contradicting the model contract and the observed series — **High**
- **Location:**
  - `scripts/model/setup.R:159-172, 201-202`
  - `scripts/scenarios_2026.R:55-62`
  - Contract: `docs/AGENTS.md` ("2016–2025 inflow uses the official SingStat/CNB age-specific series …").
  - `docs/equations.tex:5-7`
- **Problem:**
  - Entry into D is held at the 2015 CNB counts in every projection year.
  - The repository's own `data/singstat_new_drug_abusers.csv` shows new drug abusers falling from 1,309
    (2015) to 802 (2022) and 1,165 (2025).
  - The protected model contract and the equation reference both say that the observed series is used.
- **Impact:**
  - Incidence (and, more slowly, deaths) in 2026–2030 depends on how many susceptibles enter. A constant
    inflow overstates entry by 11–39% in 2020–2025 and biases the 2030 incidence ratio.
  - The code contradicts its own specification.
- **Recommended fix:**
  - Add the observed series as a scenario sensitivity analysis.
  - Leave the main analysis unchanged: the user chose constant inflow on 2026-10-30, per `claude.md`.
  - Ask the authors to reconcile `docs/AGENTS.md` and `equations.tex` with that decision.
- **Verification:** a sensitivity run with the observed 2016–2025 series (2025 held to 2030), with every
  2030 reduction reported against the main analysis.

### R-03 · Population definitions are inconsistent; the manuscript's "8,493 vs 3,367" comparison is misleading — **High**
- **Location:**
  - `scripts/model/preprocess.R:40-55, 153-167`
  - Manuscript abstract limitation ("a modelled prison population (8,493) far above the 3,367 screened
    people recorded")
  - `docs/training.md` §4
- **Problem:**
  - Release, re-arrest and π come from **all** 20,720 people held for drug consumption.
  - Seroprevalence comes from the screened subset (3,367 at the end of 2015).
  - Entry is new drug abusers × 0.861 (heroin + methamphetamine).
  - J is calibrated to none of these counts. Its posterior value is 8,493 (7,687–9,783).
- **Evidence (E11):** the SPS records hold **8,725** people for drug consumption on 31 Dec 2015, so the
  model's J matches the population the rates were fitted to. It is 2.5× the **screened** subset only.
- **Impact:**
  - The manuscript presents a definitional mismatch as a large model misfit.
  - The real open question is the opposite one: the seroprevalence target (screened, possibly selected
    people) is applied to a J that represents all consumption inmates, and J represents all consumption
    inmates while the model calls them PWID.
  - Prison treatment counts and costs scale with J.
- **Recommended fix:**
  - Authors decide whether J represents all drug-consumption inmates (then compare with 8,725) or only
    PWID (then the rates and entry need restricting).
  - Revise the limitation sentence accordingly.
- **Verification:** the manuscript compares J with the population its rates describe; `docs/training.md` §4
  records both counts.

### R-04 · The PWID proxy (heroin + methamphetamine = 86.1% of new abusers) is unverified — **High** (documented by the authors as a limitation)
- **Location:** `setup.R:159-171`; `data/calibration_targets.json` (`first_arrests`)
- **Problem:**
  - Methamphetamine users are 1,014 of the 1,127 counted. The repository has no evidence on what share of
    them inject.
  - The proxy fixes both the entry lower bound and the first-arrest target, and the posterior sits at that
    bound (R-12). The proxy therefore sets the size of the injecting population almost directly.
- **Recommended fix:** a source for the injecting share of methamphetamine users in Singapore, or a
  sensitivity analysis on the share. This is a research decision, so it is documented, not changed.

### R-05 · Uncertainty intervals cover only 13 calibrated parameters, and their width is driven by fixed, arbitrary likelihood precisions — **High**
- **Location:**
  - `data/calibration_targets.json` (`beta_precision` 100, `first_arrests_nb_size` 50)
  - Every fixed parameter in `setup.R`
- **Problem:**
  - Fixed at point values: transmission probability q, clearance, SVR, progression, liver mortality, SMRs,
    and the prison rates λ2, λ3 and π (which have Stan posteriors that are not propagated).
  - The seroprevalence likelihood uses φ = 100 whatever the sample size (screened n = 122–995 per age). The
    user chose this over a binomial likelihood.
- **Impact:** the 95% UIs are conditional on all of the above. "Joint pass probability" values such as
  5.7% are therefore not probabilities in any population sense.
- **Recommended fix:**
  - Label the UIs as conditional.
  - Propagate the λ2/λ3/π Stan draws, and run a one-way or probabilistic sensitivity analysis on the most
    influential fixed parameters (SMR, `mu_DC`/`mu_HCC`, q, p_CC_DC).
  - These are research decisions, so they are documented here.

### R-06 · Calibration seroprevalence targets cannot be reproduced from the repository — **Medium** (provenance)
- **Location:**
  - `data/calibration_targets.json`
  - `data/README.md` ("values rebuilt from an earlier NPE run's `x_obs.csv`")
  - `scripts/model/preprocess.R:162-167`
- **Evidence:**
  - The targets equal 13/102, 192/1001, 242/806, 364/754, 296/557 and 63/147.
  - The denominators are the current-age prison counts.
  - The current preprocessing gives 15/122, 191/995, 241/803, 361/747, 296/557 and 61/143, using age at
    admission.
  - The largest difference is 0.45 percentage points. That is small next to the Beta(φ = 100) noise (SD
    about 3–5 pp), but the numerators come from code that is not in the repository.
- **Recommended fix:** the authors choose the age definition and regenerate the targets from
  `preprocess.R`, or document the original derivation. This is not changed here, because changing the
  targets is a calibration decision.

### R-07 · Continuous integration fails on HEAD — **Medium** (reproducibility; engineering)
- **Location:** `.github/workflows/ci.yml:68-80`
- **Evidence (E12):** the last run on this branch failed. `python -m py_compile scripts/npe_train.py`
  can't find the file, because NPE was removed on 2026-10-30.
- **Fix:** replace the job with a syntax check of the Python scripts that exist.
- **Verification:** the job passes locally (`py_compile` on every tracked `.py` file).

### R-08 · Docker image cannot run the test suite; dependencies are unpinned and undocumented — **Medium**
- **Location:** `Dockerfile:25`; `docker-compose.yml`; no `renv.lock` or `requirements.txt`
- **Problem:**
  - The image installs `Rcpp RcppArmadillo dplyr ggplot2 MASS patchwork` but not **jsonlite**.
    `read_targets()` stops without it, so the image's default `CMD` (the tests) fails.
  - It also lacks `posterior` and `writexl`, which the MCMC and scenario scripts need.
  - The header still refers to the removed `scripts/run_npe.R`.
  - The figure and manuscript scripts need Python packages that are listed only in prose
    (`docs/training.md`).
- **Fix:**
  - Add the missing R packages and remove the stale NPE text.
  - Add `requirements.txt` with the versions already stated in `docs/training.md`.
- **Verification:** the Docker image builds and runs the tests, where Docker is available. Otherwise, a
  script reports any package the R scripts load that the Dockerfile does not install.

### R-09 · Figures and manuscript build depend on files outside the repository — **Medium**
- **Location:** `scripts/figures.py:24-32`; `scripts/build_paper.py:21`
- **Problem:** both scripts import code from `~/.claude/skills/scientific-visualization` and
  `~/.claude/skills/scientific-writing`, which are user-level tool files, not part of the repository. An
  independent researcher can't redraw the figures.
- **Recommendation:** vendor the few style and export helpers into `scripts/`, or document the exact files
  and versions needed. This is a packaging decision for the author, so it is documented here.

### R-10 · Stan rate fits are not seeded — **Medium**
- **Location:** `scripts/model/preprocess.R:202, 218-225, 257, 273-280`
- **Problem:** `set.seed(114514)` seeds R, but `cmdstanr::sample()` gets no `seed`, so every refit draws a
  new random seed. The rates in `setup.R` are rounded to 3 significant figures, so the effect is small, but
  a refit cannot be repeated exactly.
- **Fix:** pass `seed = 114514` to both `model$sample()` calls.
- **Verification:** two refits give identical draws. This needs `data/sps.csv`, which is local and not in
  git.

### R-11 · Protected model documentation contradicts the code — **Medium**
- **Location:** `docs/AGENTS.md`, `docs/equations.tex`
- **Examples:**
  - `docs/AGENTS.md` lists 12 parameters, prison-population targets with `sigma_pop = 0.10`, NPE as the
    cross-validation, and a 2016–2025 SingStat inflow schedule. The code has 13 parameters and none of the
    others.
  - `equations.tex` still uses ω-scaled mortality.
  - It applies ψ to background mortality after SVR; the code applies ψ to the liver excess only.
  - It has a separate λ⁽¹⁾; the code uses λ1 = c·λ3.
- **Impact:** a reader who checks the equations against the code finds them inconsistent. `AGENTS.md`
  requires that the two stay consistent.
- **Recommendation:** these files are marked protected, so they are not edited here. The author should
  update them.

### R-12 · Entry scales are identified mainly by their lower bound — **Medium**
- **Location:** `scripts/calibration/likelihood.R:13-15`
- **Evidence:**
  - Best-fit entry scales θ7–θ12 are 0 to 0.13 on the log scale, i.e. at or near the bound of 1.
  - In the posterior, 7–33% of draws per age lie within 0.05 of the bound.
  - The model's first arrests sit above target at ages 30+ (`docs/training.md`).
- **Impact:** the bound (a user decision) carries much of the information about entry, and so about
  population size and incidence. The tension is a sign that the model structure is misspecified rather
  than a parameter estimate.
- **Recommendation:** report a sensitivity analysis with the bound relaxed. Documented, not changed.

### R-13 · Structural simplifications that need justification — **Low to Medium**
Each is documented, not changed, because each is a model-design decision:
- **(a) Age bands.** Ageing is 0.1/yr out of every band (`sim.cpp:18-19`). The `<20` band of PWID is
  narrower (about 15–19), which makes the stay in `<20` too long.
- **(b) 60+ deaths.** The 60+ background rate includes people aged 85+ (`setup.R:110-116`).
- **(c) ψ during treatment.** The SVR liver-mortality modifier ψ applies to people still on treatment
  (state t) as well as to those cured (`sim.cpp:80-83`).
- **(d) Meaning of X.** X means "never re-arrested", which is not the same as "stopped injecting". People
  who keep injecting but are not re-arrested get no infection risk and no onward transmission.
- **(e) No treatment before 2026.** In the main analysis nobody is treated anywhere from 2015 to 2026, and
  the 2015 state is an equilibrium.
- **(f) Repeated screening.** Screening is a fresh year-start snapshot each year with no memory. Everyone,
  including people already known to be antibody-positive or cured, is re-tested every year when the quota
  allows. This inflates test costs.
- **(g) Drug price.** SOF/VEL is priced at a US list price (`cost_evidence.csv`); the authors acknowledge
  this.

### R-14 · Scenario logic has no unit tests; no test checks the full right-hand side or the sampler — **Medium**
- **Location:** `scripts/scenarios_2026.R` (functions defined inline in the script); `tests/`
- **Problem:**
  - The integration test checks individual flows. Nothing checks the complete RHS against an independent
    derivation.
  - Nothing checks that `arwm_budget` samples a known target correctly.
  - `screen_plan`, `agg_draws`, `j_entries` and the cost code can't be sourced on their own.
- **Fix:**
  - Add (i) an independent-RHS equivalence test (E2 as a regression test) and (ii) an ARWM test on a
    Gaussian target with known moments.
  - Leave the scenario script's structure as is (changing it is a refactor, not a fix) and document that.

### R-15 · Fragile step count in `run_sim` — **Low**
- **Location:** `sim.cpp:408`
- **Problem:** `n_steps = (int)((t_end - t_start)/dt) + 1` truncates. For dt values that aren't exactly
  representable, the run can stop one step early. It is exact for every interval used here (checked: 54,750
  and 366 steps). `run_year` checks the endpoint, but the burn-in does not.
- **Fix:** round to the nearest step and stop if `t_end − t_start` is not a whole number of steps
  (relative tolerance 1e-9).
- **Verification:** results are unchanged for dt = 1/365 and 1/12, and an error is raised for a
  non-dividing dt.

### R-16 · Output provenance is machine-specific or incomplete — **Low**
- `calibration_source.txt` records absolute paths.
- `mcmc_run_config.json` has no git SHA and no hashes of `setup.R` and `sim.cpp`; `run_config.json` does.
- **Fix:** write paths relative to the repository root and add the git SHA and file hashes to the MCMC
  config.

### R-17 · Negative values are clamped silently — **Low**
- **Location:** `sim.cpp:484-486`
- **Problem:** negative compartments are set to 0 after each step without any record, which can hide
  instability or break mass conservation.
- **Recommendation:** count clamping events or warn above a threshold. Not changed here: the RHS was
  verified (E2–E4), and no clamping-sensitive behaviour was seen.

---

## 3. Items checked and found correct
- The RHS, treatment flows, prison flows, reinfection of r and SVR branching all match an independent
  derivation (E2).
- Mass balance holds (E3).
- RK4 at a daily step has converged (E4).
- `event_rates` matches the ODE death terms (E5).
- The 2015 baseline is a true equilibrium for the posterior draws: 2016 differs from 2015 by less than
  1e-3, asserted per draw in the code.
- The Metropolis acceptance rule is valid for the symmetric mixture proposal. Adaptation is frozen before
  retained draws, so the retained chain is a valid MCMC. Out-of-bounds proposals are rejected (truncated
  prior).
- The WHO pass logic in `add_pass` (incidence ≤ 0.20 and deaths ≤ 0.35 of 2015) matches the stated
  80%/65% targets.
- Discounting: (1.03)^−(year−2026), asserted.
- Courses × S$40,405.30 = the DAA cost, asserted.
- Per-draw derivation of incremental costs and pass flags is correct.

## 4. Citation check (E13)

**How this was done:**
- A literature sub-agent searched the web for the 16 sources behind the fixed parameters in
  `docs/paper/si.md` Table A and the WHO targets. Access ranged from abstract only to full text, as stated
  per row.
- The reviewer did not re-read every full text. Treat the verdicts as **verified by literature search**,
  not as an author-level audit.
- Under that caveat, the main finding (Kanwal 2014 reports adjusted HRs of 1.31 and 1.80) matches the
  published abstract.

| # | Parameter (value) — source | Verdict | Note |
|---|---|---|---|
| 1 | q = 0.0057 — Boelen 2014 | Partial | Per injecting episode **with sharing, source status unknown**. The model multiplies q by the infectious share of contacts. The q×C product is calibrated, so the effect is absorbed into the contact scales, but the definitions differ. |
| 2 | κ = 0.26 — Micallef 2006 | Supported | Pooled 0.26 (0.22–0.29); few injecting-drug-use cohorts in the pool. |
| 3 | ι₁ = 180 days — Fasano 2024 | Partial | Six months is the *definition* of the acute phase, not a measured mean duration. |
| 4 | Genotype 3 share 0.78 — Zhao 2019 | Unverifiable | The abstract has no genotype split; the full text is paywalled. |
| 5 | SVR 0.976 / 0.988 — EASL 2020 | Partial | Derivable from the real-world per-protocol SOF/VEL data EASL quotes (a secondary source). Per-protocol figures overstate intention-to-treat effectiveness. |
| 6 | DC SVR without RBV, GT3 0.45 / other 0.90 — Xiao 2025 | Partial | GT3 values are in the abstract (0.45 rests on 7/15 patients). "Other" 0.90/0.92 are not found. |
| 7 | NC→CC 0.027 — Thein 2008 | Supported (derived) | Mean time to F4 from the stage-specific rates is about 37 years. |
| 8 | CC→DC 0.039, CC→HCC 0.014, DC→HCC 0.014, μ_DC 0.13, μ_HCC 0.43 — Lim 2018 | Partial (secondary) | In Lim's supplement (UK cohorts). **The footnote says they were calibrated for Pakistan's ~80% GT3 share.** |
| 9 | GT3 relative risks 1.36 / 1.93 — Kanwal 2014 | **Not supported** | Kanwal's adjusted GT3 HRs are **1.31 (1.22–1.39)** for cirrhosis and **1.80 (1.61–2.03)** for HCC. 1.36 is a calendar-period HR; 1.93 is a White-subgroup HR. Applying them to CC→DC and DC→HCC is not supported either. |
| 10 | Post-SVR CC→DC 0.07 — van der Meer 2012 | Supported | Liver-failure HR; interferon era; Ishak 4–6. |
| 11 | Post-SVR CC→HCC 0.23 — Morgan 2013 | Supported | HCC HR in advanced disease. |
| 12 | ψ_DC 0.5 — Premkumar 2024 | Unverifiable / probably unsupported | The abstract has no SVR vs non-SVR mortality ratio. |
| 13 | ψ_HCC 0.37 — Dang 2020 | Supported | 5-year **all-cause** mortality HR (the model applies it to liver deaths). |
| 14 | SMR 9.52 / 4.73 — van Santen 2014 | Partial | Table 3 values confirmed. **Antibody-based**, **all-cause (includes liver deaths)**, 1985–2012, people who *use* drugs (72% ever injected), ages 20–64. |
| 15 | 11,000 (6,500–17,000) — Teo 2019 | Partial | Network scale-up of "intravenous drug users", fielded Jul–Aug 2017, with no injecting recall period. |
| 16 | WHO 80% / 65% | Supported, with caveats | Headline incidence target 90% (HBV + HCV). Mortality 65% refers to HBV + HCV deaths from acute hepatitis, cirrhosis and liver cancer. The 2021 validation guidance adds absolute targets (not opened). |

New issues raised by this check are listed below. They are **DOCUMENTED** because each would change
calibrated inputs and require a full recalibration and an author decision.

### R-18 · Genotype 3 relative risks mis-cited — **High → Low (resolved 2026-11-03: derived values)**
- **Location:** `setup.R:90-94`; SI Table A
- **Problem:** the model uses 1.36 and 1.93. Kanwal 2014 reports 1.31 and 1.80. The model also applies the
  cirrhosis HR to CC→DC and the HCC HR to DC→HCC, which Kanwal does not estimate.
- **Impact:**
  - Progression is overstated. At ρ = 0.78, the CC→HCC multiplier is 1.73 instead of 1.62, and NC→CC is
    1.28 instead of 1.24.
  - This feeds liver deaths and so the mortality outcome.
- **Fix:** correct the values (or cite the true source of 1.36/1.93) and recalibrate. Not changed here,
  because a full refit and new scenarios take hours and change the published results; the author should
  make that call together with R-19.

### R-19 · Genotype 3 is probably counted twice in progression — **REFUTED (2026-11-03)**
- **Original claim:** Lim 2018's base rates are already for an ~80% GT3 population, so multiplying them by
  the GT3 relative risk counts GT3 twice.
- **Follow-up check (Lim 2018 full text and supplement, PMC5913612):**
  - 0.039 (CC→DC) and 0.014 (CC/DC→HCC, one parameter ξ) are in Table S2.
  - Their sources are UK cohort and HTA studies (Grieve 2006, Shepherd 2007, Wright 2003), where about
    50% of patients had genotype 1.
  - Lim applies a separate GT3 adjustment (1.30/1.80, Kanwal 2014) for Pakistan's 80% GT3, in its own row
    of Table S2. The model does the same, so GT3 is not counted twice.
  - One footnote ("calibrated to reflect the higher proportion of HCV genotype 3") is ambiguous. The listed
    values equal the UK values, so the evidence favours "no double count".
- **Residual differences (author decision: keep as is):**
  - Lim's values are annual probabilities; the model uses them as rates. −ln(1−p) differs by 1–2%.
  - Lim's GT3 adjustment is 1.30/1.80 at an 80% share; the model's is 1.36/1.93 at 78%.
- **R-18 resolution:** 1.36/1.93 are not Kanwal's reported values. They were derived from Kanwal's
  genotype-specific HRs relative to a prevalence-weighted composite of the other genotypes
  (`Hepc_modelling_technical_details.pdf` §2.4). The SI must cite them as derived from Kanwal 2014.

### R-20 · SMR misuse — strengthens R-01 — **High**
- The 9.52 SMR is all-cause and includes liver deaths. Adding μ_DC/μ_HCC on top double-counts liver
  mortality.
- HCV status in the cohort is antibody-based, yet the model gives the "non-viraemic" 4.73 to cured people
  and treats 9.52 as a viraemia effect.
- See R-01 for the size of the effect on the results.

### R-21 · Other source mismatches — **Low to Medium**
- R-21a: ψ_DC 0.5 has no visible source.
- R-21b: ψ_HCC is an all-cause HR but is applied to liver deaths.
- R-21c: SVR values are per-protocol.
- R-21d: "Other" DC SVR values were not found in Xiao 2025.
- R-21e: the WHO wording should say "65% reduction in HBV + HCV deaths", and the incidence target needs
  care (80% for HCV within a 90% combined headline).
- R-21f: ι₁ = 180 days is a definition, not a measured mean.
- R-21g: the genotype 3 share 0.78 could not be verified.

### R-22 · Nelder–Mead modes not tightly converged — **Low**
- **Evidence:** the 7 starts that "agree" (objective 18.77–18.84) differ by up to 0.64 posterior SD (θ2),
  0.38 (θ1) and 0.43 (θ3). Range across starts ÷ MCMC posterior SD.
- **Impact:** low. The mode only serves `check_fit_current()` and the documentation; the MCMC starts
  from prior draws. Calling it "the posterior mode" overstates its precision.

## 5. Reproduction from scratch

**Environment:** fresh `git clone` of `5fc2b95`; R 4.5.1 on macOS, 16 cores. The two long runs shared
the machine, so wall times are longer than documented.

| Stage | Command (as documented in `docs/training.md`) | Outcome |
|---|---|---|
| Stan rate fits | `Rscript scripts/model/preprocess.R` (needs the untracked `data/sps.csv`), twice, with the R-10 seed fix | Draws **bit-identical** between the two runs. λ2 = 0.462, 0.636, 0.454, 0.337, 0.325, 0.394; π = 0.724: identical to `setup.R`. λ3 = 0.549, 0.478, 0.469, 0.471, 0.384, 0.265: `setup.R` has 0.548 for `<20` (third-figure Monte Carlo difference from the earlier unseeded fit). `prison_2015_counts.csv` reproduced exactly. The stored fits were not replaced. |
| Nelder–Mead | `run_calibration.R --run-id r_repro --maxit 5000 --nm-restarts 5 --n-starts 8 --sd-perturb 0.9 --seed 20261001 --cores 8` (about 2 h under load) | **Bit-identical:** best start `perturb_5`, objective 18.7658, θ identical (max difference 0), `solutions.csv` identical; all 8 start objectives identical. |
| Reference MCMC | `run_mcmc_reference.R --out-dir output/calibration/mcmc_repro --cores 8` (6,098 s under load) | **Bit-identical:** every θ, lp, acceptance and evaluation column of `mcmc_draws.csv` matches; only `elapsed_sec` (wall clock) differs. `mcmc_diagnostics.csv`, `mcmc_parameter_summary.csv`, `mcmc_warmup.csv`, `mcmc_per_chain.csv` and `covariance_drift.csv` are byte-identical. |
| Scenarios (main) | `scenarios_2026.R --root . --draws-per-seed 2000 --cores 8 --sample-seed 20261005` | **Bit-identical** (E6). |
| Tests | `Rscript scripts/run_tests.R --scope all` | Pass (before: 6 files, 105 checks; after the fixes: 8 files, 119 checks). |
| Figures / manuscript | `python scripts/figures.py .`, `python3 scripts/build_paper.py .` | **Not reproduced.** The figures need the external skill files (R-09) and matplotlib 3.11.2, which is not installed here. The manuscript sources were mid-edit by the author and were not rebuilt. |
| Docker | `docker build .` | **Not run:** no Docker daemon was available. The Dockerfile fix (R-08) is unverified in a container. |

## 6. Fixes implemented and verification

| ID | Status | What changed | Verification |
|---|---|---|---|
| R-01 | **FIXED (reporting) / DOCUMENTED (definition)** | `event_rates()` also returns `liver_mortality` (DC/HCC liver deaths only). The scenarios report `liver_mortality*` and `liver_mortality_reduction*` (median, UI, text) in the summary, trajectories, no-treatment and per-draw files. Pass/fail still uses the main definition. | New test: deaths = liver + viraemic excess exactly. Main rerun: **every pre-existing column of every result CSV is identical**. Results below. |
| R-02 | **FIXED (sensitivity) / DOCUMENTED (main design)** | `scenarios_2026.R --inflow observed`: 2016–2025 entry follows the observed CNB series by age relative to 2015; 2025 held to 2030. Output in `output/scenarios_2026_inflow/` with a README. | Run complete: 2,000/2,000 valid draws. Results below. |
| R-03, R-04, R-05, R-06, R-09, R-11, R-12, R-13, R-18, R-19, R-20, R-21 | **DOCUMENTED** | Each needs an author or research decision, a recalibration, or an edit to a protected document; details are in each issue. | — |
| R-07 | **FIXED** | CI job replaced with a syntax check of the existing `scripts/*.py` and `scripts/paper/*.py`; stale comments fixed. | `py_compile` passes under Python 3.11.16 (the CI version), on both HEAD and the working tree. |
| R-08 | **FIXED (unverified in container)** | Dockerfile installs `jsonlite posterior writexl dplyr lubridate ggplot2 patchwork scales tidyr viridisLite`; stale NPE text removed; new `requirements.txt` lists the versions from `docs/training.md`. | The package list was checked against every `library()`/`::` use in the R scripts. No Docker daemon was available to build the image. |
| R-10 | **FIXED** | `seed = 114514` passed to both `model$sample()` calls. | Two refits give bit-identical draws (table above). |
| R-14 | **FIXED (partly)** | New `tests/integration/test_independent_rhs.R` (6 checks: full-RHS equivalence, mass balance, death split, step count) and `tests/unit/test_mcmc.R` (8 checks: the sampler recovers a correlated Gaussian and a half-normal truncated by −Inf). The scenario script's inline functions still have no unit tests. | All 8 test files pass. |
| R-15 | **FIXED** | `n_steps = floor((t_end − t_start)/dt + 1e-9) + 1`. | Old vs new `sim.cpp`: 150-year burn-in and annual runs **bit-identical**. For dt = 0.1, t_end = 0.3 the old code stopped at 0.2 and the new code reaches 0.3 (regression test added). |
| R-16 | **FIXED** | `calibration_source.txt` is repo-relative; `mcmc_run_config.json` records the git SHA and the md5 of `setup.R`/`sim.cpp`. | Main rerun writes `Sources: output/calibration/mcmc_reference/mcmc_draws.csv`. |
| R-17, R-22 | **DOCUMENTED** | — | — |

Other documentation fixes:
- `output/README.md` no longer lists the deleted NPE folders and now lists `mcmc_sbc/` and
  `scenarios_2026_inflow/`.
- `docs/training.md` §5 documents `--inflow observed` and the liver-only outcome.
- The test READMEs describe the new tests.
- `output/scenarios_2026/README.md` reports the liver-only outcome.
- The manuscript (`docs/paper`, `scripts/paper`, `docs/*.docx`) was **not** changed: it held the author's
  uncommitted edits.

### 6.1 Effect of R-01: liver-only HCV deaths (main analysis, 2,000 draws)

| Policy | 2030 fall in HCV deaths, current definition | Liver deaths only |
|---|---|---|
| J only | 32.6% (30.6–34.8) | 24.0% (22.9–25.2) |
| 5,000/yr, F 25% | 53.8% (50.7–56.7) | 40.3% (38.1–42.0) |
| 8,000/yr, F 25% (best finite) | 62.6% (59.0–65.6) | 47.3% (44.7–49.3) |
| 8,000/yr, F 75% | 49.0% (46.1–51.9) | 35.9% (34.3–37.5) |
| No cap | 85.1% (84.1–86.1) | 67.4% (67.2–67.6) |

- **The qualitative conclusion is robust:** no finite policy reaches 65% on either definition, and only
  the uncapped policy passes (100% of draws on both).
- The ranking (more tests to X lowers deaths more) is unchanged.
- **The magnitudes are not robust:** every reduction is about a quarter smaller on the liver-only
  definition. The 5.7% joint-pass probability of 8,000/yr F 25% falls to 0% on that definition.
- About half of the "HCV deaths" averted in the main results are SMR-based non-liver excess deaths.

### 6.2 Effect of R-02: observed 2016–2025 inflow (sensitivity, 2,000 draws)

- J only: incidence −54.2% (41.4–61.2), against −49.8% in the main analysis; deaths −32.5%, against
  −32.6%.
- Incidence reductions rise by up to 4.4 points (less where they are already close to 100%).
  - Finite policies meeting the incidence target at the median: **6** (adds 2,000/yr F 75%, 82.9%),
    against 5.
- Mortality is unchanged in substance. Only the uncapped policy passes. 8,000/yr F 25% reaches 63.2%,
  with joint pass 12.4% of draws, against 5.7%.
- The headline conclusions hold, but the incidence pass set is sensitive to the inflow assumption.

## 7. Final assessment

**Decision: MAJOR REVISION REQUIRED**

**What holds up (verified by execution):**
- **Computational reproducibility is excellent.** Rate fits (after the seed fix), Nelder–Mead, the
  reference MCMC and all 17 scenarios × 2,000 draws reproduce **bit for bit** from a fresh clone with the
  documented commands.
- **The simulator is mathematically correct for the model it states.** It matches an independent
  derivation to machine precision, conserves mass, and has converged numerically (about 1e-4 relative
  error at a daily step).
- The MCMC sampler recovers known targets. Its convergence diagnostics are what the docs claim.
- The headline numbers in the manuscript abstract match the outputs.

**Why it is not PASS WITH CONDITIONS:** the unresolved problems are scientific, and some of them alter
the inputs behind the main conclusions.
- **Mis-cited inputs:**
  - The genotype 3 relative risks (1.36/1.93) are not the cited Kanwal estimates (1.31/1.80) (R-18).
  - Genotype 3 is probably applied twice to progression rates that were already calibrated for a GT3-heavy
    population (R-19).
  - Both inflate liver progression and liver deaths. Their effect is **not quantified**, because fixing
    them needs a full recalibration, which is an author decision.
- **The mortality outcome is not the WHO's:**
  - It is half SMR-based non-liver excess from an antibody-based, all-cause SMR, with probable
    double-counting of liver deaths (R-01, R-20).
  - On a liver-only definition every mortality reduction is about a quarter smaller. The conclusion "no
    finite policy meets the mortality target" survives, and gets stronger.
- **Population definitions conflict** (R-03), and the PWID proxy is unverified (R-04). The manuscript
  misreads a definitional mismatch as a 2.5× misfit.
- **The UIs are conditional** on 13 calibrated parameters and on arbitrary likelihood precisions (R-05),
  so the "probability of meeting the targets" figures should not be read as such.
- **The design deviates from its own contract:** constant inflow despite a falling observed series
  (R-02). In the sensitivity run this adds one finite policy to the incidence pass set.
- The protected equation and contract documents contradict the code (R-11). The seroprevalence targets
  cannot be regenerated from the repository (R-06).

**Direction of the results:** the qualitative message (starting in 2026 is too late for capped
screening to meet the mortality target, and incidence needs thousands of community tests) survived both
sensitivity analyses run here. It has **not** been tested against R-18/R-19 or the fixed-parameter
uncertainty.

**To reach "pass with conditions":**
1. Correct R-18/R-19 and recalibrate.
2. Choose and justify the HCV-death definition (R-01/R-20) and report the other one as a sensitivity
   analysis.
3. Reconcile the J population definition and revise the manuscript limitation (R-03).
4. Propagate the λ2/λ3/π posteriors and run a sensitivity analysis on the key fixed parameters (R-05).
5. Regenerate the seroprevalence targets from code (R-06).
6. Update `docs/AGENTS.md` and `equations.tex` (R-11).
7. Vendor or document the figure dependencies (R-09).

**Not verified in this review:**
- The Docker build (no daemon).
- Figure and manuscript regeneration (external skill files; matplotlib 3.11.2 not installed; manuscript
  mid-edit).
- Full-text checks of the paywalled sources (rows 4 and 12 of §4).
- Any re-analysis of the raw SPS data beyond rerunning `preprocess.R`.

Nothing was committed. All changes are in the working tree next to the author's earlier uncommitted
manuscript edits, which were left untouched.

## 8. Follow-up: author decisions and recalibration (2026-11-03/04, branch `review/liver-deaths-recalibration`)

Sections 0–7 record the review as run. This section records what the authors decided afterwards and
what was changed, rerun and checked. Nothing is committed.

| ID | New status | What changed |
|---|---|---|
| R-01, R-20 | **FIXED (model + outcome)** | One route per cause of death: non-liver deaths = `mu` × 4.73 (`smr_nonviraemic`, van Santen 2014, HCV/HIV-uninfected PWID) for u, r and × 4.73 × 1.47 (`smr_viraemic`, Lee 2012 extrahepatic HR) for a, c, t; liver deaths only from `mu_DC`/`mu_HCC`. Primary outcome (`mortality`) = liver deaths (the WHO definition); `hcv_mortality` = liver + extrahepatic excess is the secondary outcome (S1 Text Section 10). The §6.1 column `liver_mortality` no longer exists (the main column is now liver-only). Background `mu` for 60+ uses ages 60–74 (10.14/1000, SingStat M810141/M810011; author decision) instead of 21.88. |
| R-18 | **RESOLVED** | 1.36/1.93 are derived in §2.4 of the technical PDF from Kanwal's genotype-specific HRs; kept (author decision). |
| R-19 | **REFUTED** | See the R-19 section: Lim 2018 base rates are non-GT3. |
| R-02 | **DOCUMENTED + SENSITIVITY** | Constant inflow kept for the main analysis (author decision); observed series in S1 Text Section 11; inflow passages of `docs/AGENTS.md` and `docs/equations.tex` updated (authorised). |
| R-03 | **FIXED (manuscript)** | `preprocess.R` writes `n_held` (8,725 held for drug consumption at end-2015) beside the screened subset (3,367). The manuscript compares modelled J (8,768, 8,000–10,110) with 8,725; the limitation is now the representativeness of the screened subset. |
| R-05 | **PARTLY FIXED** | Conditional UI wording in the Methods; one-way sensitivity of 27 fixed inputs at the posterior mode (`scripts/tornado.R`, `output/tornado/`, S1 Text Section 12, Fig G). λ2/λ3/π are varied over their Stan 95% intervals there, not propagated through the posterior. |
| R-11 | **PARTLY FIXED** | Only the inflow passages were authorised. The mortality descriptions in `docs/AGENTS.md`/`docs/equations.tex` still describe the old SMR design and now contradict the code. |

**Reruns (all from the new fit):** Nelder–Mead `r_20261103_agepop` (8 starts within 1, objective 19.25;
seroprevalence RMSE 0.023, against 0.014 before; 60+ 0.468 vs 0.429). Reference MCMC: R-hat ≤ 1.006,
bulk ESS ≥ 1,641, 2.0 h. MCMC SBC: cov95 0.952. Calibration overlap 13/13. All six scenario runs
(2,000/2,000 valid draws each). Tornado: 55 runs, 2 without a 2015 equilibrium (`q` low: HCV dies out;
`smr` low: not converged at 150 years). Tests: 8 files, 119 checks, 0 failures. Figures redrawn with
matplotlib 3.11.2; manuscript rebuilt, every `checks()` assertion rewritten for the new results.

**Results (main analysis, liver deaths):**

| Policy | New infections | Liver deaths | Liver + extrahepatic |
|---|---|---|---|
| J only | 51.5% (40.5–57.2) | 21.2% (20.3–22.1) | 25.2% (23.8–26.7) |
| 8,000/yr, F 25% (best finite) | 85.4% (78.3–88.9) | 42.2% (39.7–44.0) | 50.3% (47.2–52.8) |
| No cap | 99.8% (98.7–99.9) | 67.1% (66.9–67.3) | 77.7% (77.0–78.5) |

- Same incidence pass set (5 finite policies). No finite policy meets the mortality target or both
  targets in any draw, in the main run, at any J coverage, with observed inflow or with J from 2019
  (best 55.0%). This conclusion is now stronger than in §7.
- The uncapped policy passes narrowly. It fails without J treatment (61.6%), and in the one-way analysis
  it fails with ψ_DC, ψ_HCC or φ_CC_DC high, or progression to HCC low.
- The 2019 sensitivity no longer has finite policies meeting the mortality target (§7's "several thousand
  tests" result under the old definition is gone).

**Still open:** R-04, R-06, R-09, R-12, R-13, R-21 as in §6; Premkumar 2024 (ψ_DC = 0.5) could not be
verified (full text inaccessible; range 0.17–0.96 used in the tornado); Lim 2018 probabilities are used as
rates (1–2% difference, left as is by the authors); `output/calibration/smr_sensitivity/` is from the
2026-10-22 design and was not rerun.

## 9. Follow-up: ψ for cured people only, new ψ values, genotype 3 sensitivity (2026-11-04)

**Author decisions:** ψ multiplies the liver death rate of cured people (r) only; people on treatment
(t) keep the untreated rate (was t and r; the technical PDF says t only). R-21a/R-21b resolved by
switching to liver-related death hazard ratios, SVR vs no SVR, read from the abstracts:
ψ_DC 0.40 (0.17–0.96, Pageaux 2022, HEPATHER) and ψ_HCC 0.34 (0.13–0.88, Dang 2020). Premkumar 2024 is
no longer used. Genotype 3 share 10% and 50% added as a sensitivity analysis (main 78%).

**Code:** `sim.cpp` `mortality()`, `model_metrics.R` `event_rates()`, `tests/integration/test_independent_rhs.R`
(h == r), `setup.R`, `scripts/tornado.R` (ψ_HCC range 0.13–0.88). The R-vs-C++ death-rate test
(`test_sim_and_pipeline.R`) now compares instantaneous rates (one 1e-6-year step, tolerance 1e-5); with a
one-day step its within-day depletion (0.14%) exceeded the 1e-3 tolerance once t no longer had ψ.

**Calibration unaffected:** nobody is cured before treatment starts, so ψ does not enter the 2015
equilibrium. At the NM best fit and at 200 MCMC draws the 150-year trajectories under old and new ψ
were identical (maximum difference 0). NM rerun into `r_20261103_agepop`: all 8 starts identical
(objective 19.24979). MCMC, SBC and 2015 predictions were not rerun (target unchanged).

**Reruns:** tornado (55 runs; same 2 failures), all six scenario runs, two genotype 3 runs
(`output/scenarios_2026_gt3/`; `--param-file`, posterior not refitted: at the mode ρ 0.10 raises the
target negative log-likelihood by 0.27). Figures redrawn; manuscript rebuilt (new S1 Text Section 12,
Table M; tornado now Section 13/Table N; assumptions Section 14/Table O); tests 8 files, 119 checks, 0 failures.

| Liver-death reduction, 2030 | §8 | Now |
|---|---|---|
| J only | 21.2% | 22.1% (21.2–23.0) |
| 8,000/yr, F 25% (best finite) | 42.2% | 44.0% (41.4–45.9) |
| No cap | 67.1% | 71.4% (71.3–71.6) |
| No cap, J coverage 0 | 61.6% (fail) | 65.4% (65.1–65.7, pass) |
| No cap, J from 2019 | 73.3% | 76.4% |
| Best finite, J from 2019 | 55.0% | 56.1% |

- Incidence results unchanged. No finite policy meets the mortality target anywhere.
- No cap now meets both targets at every J coverage. In the one-way analysis only ψ_DC high (0.96) makes it
  fail (59.4%); ψ_HCC, φ_CC_DC and p_HCC no longer flip it.
- Genotype 3 share 10% / 50%: identical pass sets; J only liver deaths −19.7% / −21.2%, no cap 71.1% / 71.3%.

**Still open:** the SVR-vs-no-SVR hazard ratios may favour cure (healthier patients reach SVR); Pageaux's
treated-vs-untreated comparison found no significant effect on liver deaths. `docs/tex_consistency_instructions.md`
lists the changes needed in the protected `docs/equations.tex` and the PDF Notation table (not applied).
