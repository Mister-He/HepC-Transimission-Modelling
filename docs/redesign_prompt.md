# Redesign: entry scales, first-arrest target and community PWID target

Status: implemented 2026-10-27. User decisions on Section 5 override the draft below:
- community target is Teo et al. 2019 (about 11,000 inject, 95% CI 6,500–17,000; log SD 0.2453),
  not Quaye; community = D + F (X excluded);
- entry scales are bounded at ≥ 1; first arrests compared in 2015;
- both audit fixes (ψ in `sim.cpp`, `split_by_age` in `preprocess.R`) are included.
Results: `CLAUDE.md` ("Redesign") and `docs/training.md`.

## 0. Goal and boundary

- Redesign the calibration and rerun the whole chain: rate estimation → calibration → 2015
  predictions → calibration figure → treatment scenarios → manuscript.
- **Change only what Sections 2–3 specify.** Everything else stays as at HEAD `152e1bc`:
  - compartments D/J/F/X, six age groups, HCV states u/a/c/t/r, liver stages;
  - `sim.cpp` simulation logic (first arrest `lambda1 = c × lambda3`, mortality `mu × SMR`,
    HCV deaths definition, entry held constant in every projection year);
  - 150-year equilibrium run to 2015, RK4 with `dt = 1/365`, equilibrium check;
  - treatment scenario design, costs and WHO targets (`docs/structured_prompt.md`).
- Constraints: no push; commits per `commit_rule.md`; `data/sps.csv` aggregates only and never
  sent to Ollama; long runs with `nohup` from frozen script copies.

## 1. Estimation (unchanged code, rerun)

1. Run `scripts/model/preprocess.R` (needs `data/sps.csv`).
2. Fit `release_rate.stan` → age-specific release rate `lambda2` (6).
3. Fit `rearrest_rate.stan` → age-specific re-arrest rate `lambda3` (6) and `pi_recid`.
4. Run `scripts/model/update_rate_params.R` to write them into `setup.R`.
5. Check: values match the stored fits (`pi_recid` 0.724; `lambda3` 0.548 … 0.265). If they
   match, `setup.R` is unchanged.

## 2. Calibration

### 2.1 Simulator
- Current `sim.cpp` and `setup.R`, fed the Section 1 rates.
- **Change:** entry into D by age = `s_i × CNB_2015_i` (CNB 2015 new drug abusers: 257, 644, 292,
  94, 21, 1), held constant in every year.

### 2.2 Parameters (13)
| # | Parameter | Scale | Bounds | Role |
|---|---|---|---|---|
| θ1–θ6 | contact-matrix row scales | log | 0.01–1,000 | seroprevalence by age |
| θ7–θ12 | entry scales `s_i` (initial value 1 = CNB counts) | log | `s_i` 1–100 (log ≥ 0) | first arrests by age |
| θ13 | first-arrest multiplier `c` | log | 0.5–100 | community PWID size |

### 2.3 Calibration targets
1. **Seroprevalence:** 2015 prison HCV seroprevalence by age (6), as now.
2. **First arrests:** model annual flow D → J by age in 2015 (`c × lambda3_i × D_i`) against
   CNB 2015 new drug abusers by age (6).
3. **Community PWID:** model D + F in 2019 against Quaye et al. 2023 (JIAS 26:e25973):
   3,470 (95% CI 1,540–7,830).

### 2.4 Prior, likelihood, algorithm (proposed; user asked me to choose)
- **Prior** (independent, truncated to bounds, fixed before fitting):
  - θ1–θ6 ~ N(0, 1.5) (unchanged);
  - θ7–θ12 ~ half-normal on log scale, N⁺(0, 1.0) (prior median entry ≈ 2.0 × CNB);
  - θ13 ~ N(0, 0.5) (unchanged).
- **Likelihood:**
  - seroprevalence: Beta(p·φ, (1−p)·φ), φ = 100 (unchanged);
  - first arrests: negative binomial on the observed counts, size 50 (CV ≈ 0.15 at large counts;
    handles the 60+ count of 1). Sensitivity: Poisson and size 20;
  - community PWID: log-normal, centre log 3,470, SD 0.415 (from the 95% CI).
- **Algorithm:** Nelder–Mead multi-start (mode) → reference adaptive RWM MCMC (production
  posterior), same sampler and diagnostics as now (R-hat ≤ 1.01, bulk ESS ≥ 1,000). NPE + SBC as
  cross-validation only (summaries: 6 seroprevalences + 6 first arrests + 1 community count).
- **Code touch points:** `data/calibration_targets.json`, `build_params()`/`prior_spec()`/`nlp()`
  in `scripts/calibration/likelihood.R`, `mcmc.R`, `model_metrics.R` (first arrests, D + F),
  `run_npe.R`/`npe_train.py`, tests.

### 2.5 Calibration figure
- Update `plot_calibration.R` and `figures.py`: three panels — seroprevalence by age, first
  arrests by age, community PWID (model D + F vs Quaye interval); MCMC posterior and NPE.
- Update `calibration_interval_overlap.csv` to cover all 13 targets.

## 3. Treatment scenarios
- Same design: J only, 15 finite F/X policies, no cap; 2,000 posterior draws, median and 95% UI.
- Rerun `scenarios_2026.R`, `figures.py`, `build_paper.py`; update README, `docs/training.md`,
  claude.md.

## 4. Acceptance checks
1. `Rscript scripts/run_tests.R --scope all` passes (new tests for the two new targets).
2. All 8 Nelder–Mead starts reach the same mode; no parameter at a bound.
3. MCMC R-hat ≤ 1.01; report which parameters stay prior-dominated.
4. Every 95% posterior predictive interval overlaps its target (13 targets).
5. Scenario and manuscript numbers regenerate with no stale text (`build_paper.py` asserts).

## 5. Open questions (decide before coding)
1. **Blocker — F is too large:** model F ≈ 5,900 at every `c`, already above Quaye's 3,470.
   D + F cannot reach the target, so `c` would sit at its bound. Check why F is large first
   (time in F before re-arrest, `pi_recid`, mortality). ~15 min.
2. **Entry scale ≥ 1 may bind:** almost everyone entering D is eventually arrested, so first
   arrests ≈ entry, and the CNB target pins `s_i` ≈ 1. Scales above 1 are possible only through
   deaths or exits from D before arrest.
3. **Which population does Quaye measure:** assumed D + F (injecting, not in prison). Confirm
   whether X belongs in it.
4. **Target year:** first arrests compared in 2015 only (equal to every pre-2019 year at
   equilibrium). Confirm.
5. **Open audit fixes:** include the ψ mortality fix and the `split_by_age` fix in this rerun,
   or keep them out (the brief says keep everything else the same).

## 6. Time estimate
- Rates: ~10 min. Code changes and tests: ~3 h. Nelder–Mead: ~2 h. MCMC (13 parameters): ~6–8 h.
  NPE + SBC: ~3 h. Scenarios, figures, manuscript: ~1.5 h. Total ≈ 1.5 days wall time.
