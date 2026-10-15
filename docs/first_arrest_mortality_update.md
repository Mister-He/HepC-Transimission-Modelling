# First-arrest flow handoff and viraemic excess mortality

Structured prompt (2026-10-14). It supersedes §6.1 item 3 of `docs/final_review.md`
(`c_composite = c * beta_scale`). Everything else in `docs/final_review.md` and `docs/AGENTS.md`
still applies.

## 1. Objective

- Make the hand-off from the Stan rate fits to the ODE mathematically consistent.
- Remove the 60+ seroprevalence misfit with an evidence-based mortality assumption.
- Refit the rates, rerun the whole calibration and scenario chain, and update the code, docs and
  manuscript.
- Check the misfit again, run the audit and review steps in `docs/final_review.md` §8, and commit
  without pushing.

## 2. Findings that justify the change

1. **`c_composite` is redundant and mixes units.**
   - In `rearrest_rate.stan`, `c` is a free parameter for each age, fitted to a single count. It
     ends up at about `c_i ≈ f_i / (λ3_i · tot_in_i)`. Here `f_i` is the steady-state first-arrest
     flow and `tot_in_i` is the CNB new-abuser count per year.
   - So `λ1 = λ3·c_composite/beta_scale = f_i / (beta_scale_i · tot_in_i)`. This divides one flow
     by another, but the result is used as a rate per year.
   - In Stan, `beta_scale` scales the never-arrested **stock**. In `sim.cpp` it scales the entry
     **flow**. The two agree only when people spend exactly one year in D.
   - At the 2026-10-13 fit, the model's first arrests were 0.29 to 4.5 times Stan's. At 60+,
     `beta_scale` cancels from the first-arrest flow, so θ12 is unidentified.
2. **Stan's prison balance did not match the ODE.**
   - Stan used ageing `a = 0.2` for every band. The ODE uses 1/10 for every band.
   - Stan used general-population mortality `mu` (60+: 0.0065). The ODE uses `mu·omega` (60+: 0.32).
3. **The 60+ misfit is structural.**
   - The ODE moves 1/10 of each band up an age band every year. So the 60+ prison population is
     fed mainly by people ageing in from 50–59, and it inherits their seroprevalence.
   - Background mortality is the same for every HCV state. So nothing can lower 60+ seroprevalence
     except dilution by uninfected first arrests.
   - Before 7c3dad3, `beta_scale_6` (44.9 new 60+ users a year, against 1 recorded by CNB) bought
     that dilution. After 7c3dad3 the fit gave 0.508 against a target of 0.429.
4. **The evidence for excess mortality** (sources checked against the original tables):
   - **Lee et al. 2012, J Infect Dis 206:469–477 (REVEAL-HCV), Table 3:** deaths from causes other
     than the liver, compared with anti-HCV negative people.
     - Detectable HCV RNA: **HR 1.47 (1.23–1.77)**.
     - Anti-HCV positive but RNA undetectable: **0.90 (0.64–1.28)**.
     - All causes, detectable RNA: 2.20 (1.90–2.55). This includes liver deaths, which the model
       already represents in its cirrhosis and liver-cancer stages.
   - **Piselli et al. 2021, BMC Infect Dis 21:667 (Naples):** deaths from causes other than the
     liver, RNA positive vs anti-HCV negative: HR 1.11 (0.82–1.49). This is the lower sensitivity
     bound.
   - **Cepeda et al. 2022, Ann Intern Med (ALIVE, PWID):** all-cause mortality, undetectable vs
     detectable RNA: HR 0.54 (0.38–0.77). This supports the direction of the effect in PWID, but
     is not used as the HR because the comparison group is self-selected (cleared or treated).
   - **Mathers et al. 2013, Bull WHO 91:102–123:** pooled SMR for PWID **14.68 (13.01–16.35)**,
     from cohorts that mix HCV-positive and HCV-negative people.

## 3. Requirements

### 3.1 Hand over the first-arrest flow `f_i`

- **`rearrest_rate.stan` (both copies):**
  - Pass the ageing rate `a` and the death rate `mu` as data. `preprocess.R` sets them to the
    ODE's values: `a = 0.1`, and `mu·omega` from `setup.R`.
  - Remove `c_composite`, `tot_in` and `i_pop`. With one free value per age they absorb their own
    data term exactly, so they carry no information.
  - Keep `f_pop`, and save its draws.
- **`preprocess.R` (both copies):** run from the repository root, with Stan file paths resolved
  inside `scripts/model/`.
- **`update_rate_params.R`:** write `f_first = posterior mean of f_pop` into `setup.R` in place
  of `c_composite`. Keep the validation and the rule that nothing is written unless every check
  passes.
- **`sim.cpp`:** derive `λ1` from the ODE's own 2015 steady-state balance for D, before any
  infection:
  - `D_i = (β_i + a·D_{i−1} − f_i) / (mu_i·omega + a_i)`, with `a_6 = 0`;
  - `λ1_i = f_i / D_i`.
  - Here `β_i` is the scaled entry rate. A parameter set whose numerator is ≤ 0 is infeasible:
    there are not enough never-arrested people to supply the first arrests. `run_sim` stops, so
    the parameter set gets zero posterior density.
- **One implementation:** export `first_arrest_rate(params)` from `sim.cpp` for R callers
  (scenarios, tests). Drop `beta_scale` and `c_composite` from the simulator inputs.
- **`informed_theta0`:** seed `D_i = f_i/λ3_i` (first arrest as fast as re-arrest). Mark it as a
  shortcut.

### 3.2 Excess non-liver mortality for chronic viraemic people

- **Rates:**
  - Background mortality becomes `mu_i·omega_neg` for people who are not chronically viraemic,
    and `hr_viraemic·mu_i·omega_neg` for chronic HCV (state c, all strata and stages).
  - `omega_neg = omega / (1 + (hr_viraemic − 1)·v)`, so the PWID average stays at Mathers'
    `omega`.
  - The HR multiplies only the background part of cirrhosis and liver-cancer mortality, never
    `mu_DC` or `mu_HCC`. Liver deaths are already explicit, so this does not double count them.
- **Only state c.** Lee's HR describes chronic carriers. Acute infection (6 months) and treatment
  (12 weeks) are short, and cleared or SVR people (r) had HR 0.90.
- **`hr_viraemic = 1.47`.**
- **`viraemic_share_ref = v`** is frozen as a constant. It is the chronic viraemic share among all
  PWID at the 2015 snapshot of the Nelder–Mead best fit. Fixed point: start at 0.20, refit while
  the implied `omega_0` still moves by more than 0.5%, and record `v` in the fit config.
- **Why `v` is frozen:** if `v(t)` changed over time, treating one person would raise everyone
  else's death rate.
- **Remove** the dead `eta_s` entry.
- **The HCV-death outcome is unchanged** (liver deaths only). The manuscript says the extra deaths
  are not counted.
- **Sensitivity:** Nelder–Mead refits at HR 1.11, 1.23 and 1.77, reporting the 60+ fit. Their
  scenario effects are left as future work unless they are cheap.

### 3.3 Regenerate everything in dependency order

Launch every long run from a frozen copy of its script (`scratchpad/frozen/`).

1. `preprocess.R`, then `update_rate_params.R`, then the rate tests.
2. `run_calibration.R` (Nelder–Mead, new run id `r_20261014_agepop`), then the `v` fixed point.
3. `run_npe.R --step all` (50,000 LHS draws, SBC n = 1,000). Report draws excluded as infeasible
   and do not replace them.
4. `run_mcmc_reference.R --cores 4`, in parallel with `sbc_mcmc.R --cores 8`.
5. `predict_2015_posterior.R`, `plot_calibration.R`, `scenarios_2026.R`, `compare_npe_mcmc.R`,
   `paper_figures.R`, `build_paper.py`.
6. Delete superseded outputs, e.g. `r_20261013_agepop/`.

### 3.4 Docs and code

- Update every reference to `c_composite`, `beta_scale` in `λ1`, `eta_s` and the 60+ misfit text:
  README files, `docs/training.md`, `evaluation.md`, `model_card.md`, the `scripts/model`
  README, `claude.md` and `build_paper.py` (methods, parameter table, limitations, references).

### 3.5 Checks

- **60+ seroprevalence:** report the best fit and the posterior predictive interval against the
  target (SE ≈ 0.041). Also report the fit without the HR, to separate the two effects.
- **First-arrest flow:** the model's 2015 first-arrest flow `λ1·D` must be within Stan's posterior SD of `f_i`
  (the closed form leaves out HCV-related deaths in D).
- **Gates:** MCMC R-hat ≤ 1.01; NPE SBC gate; all 12 predictive intervals overlap the targets.

### 3.6 Audit and review

Follow `docs/final_review.md` §8: ponytail-audit, refactor (preserving behaviour), ponytail-review,
then `Rscript scripts/run_tests.R --scope all`.

### 3.7 Commits (`commit_rule.md`, as few as possible, no push)

1. `feat(model)!: hand over the first-arrest flow and add viraemic excess mortality`: code, tests
   and the rate refit.
2. `feat(calibration): regenerate calibration, cross-validation and scenarios`: outputs.
3. `docs: …` / `docs(manuscript): …`, merged into one if possible.
4. Ponytail fixes, only if there are any.

## 4. Known limits

- The 2015 prison counts inform `f` in Stan and are also calibration targets. With `f` fixed,
  population fit is almost automatic, so the population term carries little information about θ.
- The HR comes from a general-population cohort in Taiwan, not from PWID. Its PWID value is
  uncertain, so the sensitivity range is reported.
- `λ1` is solved from the HCV-free balance. Excess deaths in D shift the realised flow slightly;
  §3.5 checks the size of that shift.

## 5. Decisions made during execution (2026-10-14/15)

1. **One preprocessing script.**
   - The top-level `scripts/preprocess.R`, `rearrest_rate.stan` and `release_rate.stan` were
     removed.
   - The active `scripts/model/preprocess.R` adopted their data choices (held for drug
     consumption; 2015 prison counts with a known screening result). Those choices reproduce
     `data/calibration_targets.json` exactly, and the stored release fit.
2. **Population (user decision).**
   - The modelled PWID population is the people with a known HCV screening result (unscreened
     people are assumed not to be PWID).
   - It defines the 2015 prison counts in Stan (= the calibration targets) and the first
     offenders.
   - Release and re-arrest rates use all consumption records. In the screened subset, every
     earlier release ends in re-arrest by selection (re-arrest share 0.993 against 0.720).
3. **Age bands.** `cut(..., right = FALSE)` so that exact ages 20, 30, … fall in the older band,
   as in the targets and `sim.cpp`.
4. **Excess-entry parameterisation.**
   - θ7–θ12 = log(excess entry / β_i). Entry is `max(0, f_i − 0.1·D_{i−1})` plus the excess.
   - On the log entry scale, 43% of prior draws could not supply the first arrests, mostly at
     20–29.
   - Under this parameterisation every θ is feasible.
5. **Entry anchor (user decision).**
   - The 2015 CNB new-drug-abuser counts are one year's entry into D: log-normal likelihood,
     σ = 1.0.
   - Without it, the never-arrested pool is unidentified (about 38,000 to 1.7 million PWID
     at the same objective).
6. **Nelder–Mead `maxit` 10,000.** With the anchor, 3,000 iterations left 7 of 8 starts
   unconverged (objective 23.2). Restarts reached 21.05.
