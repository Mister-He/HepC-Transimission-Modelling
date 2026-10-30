# Training

Production calibration path:
1. **Nelder–Mead multi-start** finds the posterior mode (fixed prior + seroprevalence likelihood).
2. **Reference MCMC** gives the production posterior.
3. **Scenarios** run on draws from that posterior.

Check only (nothing here feeds the production path): a small direct simulation-based
calibration (SBC) of the MCMC. The NPE cross-validation was removed on 2026-10-30.

Outputs are in `output/calibration/` (see its `README.md`). Run the steps in this order; each step
reads the previous step's output. Model design: `docs/first_arrest_multiplier.md`, calibration
redesign (entry scales, first-arrest and community-PWID targets): `docs/redesign_prompt.md`; PWID share
of the inflow and treatment from 2026 as the main analysis: `docs/pwid_inflow_main_analysis.md`.

## 0. Rate fits (only after refitting the Stan rate models)

```bash
Rscript scripts/model/preprocess.R            # needs data/sps.csv (not tracked by git); ~5 min
Rscript scripts/model/update_rate_params.R    # writes lambda2, lambda3 and pi_recid into setup.R
```

- **Data:** everyone ever held for drug consumption. Stays are Lexis-split by age; a stay that
  crosses an age boundary carries its release event only in the age band where it ends. Screening (Dec 2014–Feb 2016) selects people
  in prison then, so in the screened subset every earlier release ends in re-arrest; the rates
  therefore use all consumption records.
- **`release_rate.stan`:** age-specific release hazard (Lexis-split stays). Both fits are saved as
  `list(summary, draws)`; `lambda2` and `lambda3` are the posterior means of `exp(log_lambda)`.
- **Also writes** `output/prison_2015_counts.csv`: the 2015 prison count of screened PWID and the
  screened seropositives by age group (aggregates only; not calibration targets).
- **`rearrest_rate.stan`:** cure-type survival model only (share `pi` who are ever re-arrested and
  age-specific re-arrest hazard `lambda3`). There is no prison steady-state block and no
  first-arrest flow. Fit 2026-10-22: `pi` 0.724, `lambda3` 0.548, 0.478, 0.469, 0.471, 0.384,
  0.265 per year, R-hat 1.00. Refit 2026-10-27 (after the `split_by_age` event fix): same values
  to 3 significant figures, so `setup.R` is unchanged.
- **First arrest** is not estimated from the records (they cannot identify true first offenders):
  `lambda1 = c_first × lambda3`, with `c_first` calibrated (θ13).
- **Idempotent:** on the stored fits `update_rate_params.R` leaves `setup.R` unchanged (checked
  by `tests/unit/test_rate_params.R`). After a change, rerun every step below.

## Shared posterior target

- **Parameters (13):** θ1–θ6 are log contact scales (relative to `C_contact` in `setup.R`);
  θ7–θ12 are log entry scales `s_i` (entry into D_i = `s_i` × `beta_i`, held constant in every
  year; `beta` in `setup.R` = CNB 2015 new drug abusers × 0.861, unrounded); θ13 is log `c`, the first-arrest multiplier (`build_params()` in
  `scripts/calibration/likelihood.R`). Design: `docs/redesign_prompt.md`.
- **Bounds:** contact scales 0.01–1,000; entry scales 1–100 (user decision: the CNB counts are a
  lower bound); `c` 0.5–100 (user decision: first arrest at least half the re-arrest rate).
- **Prior (fixed, not centred on any fit):** independent Normal, truncated to the bounds:
  θ1–θ6 ~ N(0, 1.5), θ7–θ12 ~ N(0, 1) truncated at 0 (half-normal), θ13 ~ N(0, 0.5)
  (`prior_spec()`; density `nlp()`, shared by Nelder–Mead and MCMC).
- **Likelihood (`nll_targets()`), three independent terms** (`data/calibration_targets.json`):
  - seroprevalence: each age group's observed 2015 prison seroprevalence is
    `Beta(p·φ, (1−p)·φ)` around the model value `p`, φ = 100 (`beta_precision`);
  - first arrests: CNB 2015 new drug abusers by age (257, 644, 292, 94, 21, 1) × 0.861, rounded
    (221, 554, 251, 81, 18, 1), are negative binomial around the model D → J flow `c·λ3_i·D_i`,
    size 50 (`first_arrests_nb_size`);
  - PWID share 0.861 = (heroin 113 + methamphetamine 1,014) / 1,309: the 2015 new drug abusers by
    main drug (SingStat M891521, `data/singstat_drug_abusers_by_drug.csv`); no age split, so the same
    share at every age (user decision 2026-10-30);
  - community PWID: Teo et al. 2019 (network scale-up, 2017 survey) 11,000 (6,500–17,000)
    injecting drug users in Singapore; log-normal around the model D + F, log SD 0.245
    (`community_pwid_log_sd`). The model value is the 2015 equilibrium, which also stands for
    2017 (constant entry, no treatment before 2019 in any analysis).
  - Seroprevalence is unchanged when all entry scales rise together (transmission depends on the
    infectious share), so entry is pinned by the two count targets. `c` is pinned by the tension
    between them: a low `c` keeps more people in D (raising D + F) but shifts first arrests to
    older ages.
- **Support:** draws outside the bounds and equilibrium failures are excluded (MCMC gives −∞).
- **Solver:** RK4 with `dt = 1/365`.
- **Mortality:** general-population `mu` × PWID SMR by state: 9.52 viraemic (a, c, t), 4.73
  otherwise (van Santen 2014). In DC and HCC the liver excess `mu_DC`/`mu_HCC` is added, and after
  SVR only that liver excess is multiplied by ψ (0.5 DC, 0.37 HCC).

## 1. Deterministic calibration (posterior mode)

```bash
Rscript scripts/run_calibration.R --run-id r_20261030_agepop \
  --maxit 5000 --nm-restarts 5 --n-starts 8 --sd-perturb 0.9 --seed 20261001 --cores 8   # ~2 h
```

- **Objective:** negative log posterior (`nll_targets()` + prior); every evaluation runs model
  time 0–150 and checks equilibrium over the last 10 years.
- **Starts (uninformative):** θ = 0 (contact scales 1, entry scales 1, `c` = 1) and 7 Normal perturbations of it
  (SD 0.9), clipped to the bounds, run in parallel.
- **Restarts:** each start reruns Nelder–Mead from its own solution while the objective still
  falls by more than 0.001 (`--nm-restarts`, at most 5 restarts).
- **Result (2026-10-30, 1.2 h on 8 cores):** 7 of 8 starts reach the same mode (objective 18.77–18.84;
  `perturb_2` stops in a local optimum at 73.1). Best (`perturb_5`): seroprevalence RMSE 0.014, largest
  error 0.024; `c` 2.11; D + F 6,786 (Teo 11,000, 6,500–17,000).
  - Entry scales sit at about 1.0–1.15 (lower bound 1).
  - First arrests are too high at older ages (30–39 281 vs 251; 40–49 99 vs 81; 50–59 25 vs 18;
    60+ 3.2 vs 1): people who start injecting young age into the older groups before their first
    arrest, and entry cannot fall below the scaled CNB counts.
- **SMR sensitivity** (`output/calibration/smr_sensitivity/`) is from the 2026-10-22 7-parameter
  design and was not rerun.

## 2. Reference MCMC (production)

```bash
Rscript scripts/run_mcmc_reference.R   # about 5 hours on 8 cores
```

- **Defaults:**
  - 8 chains from independent prior draws, 8 cores;
  - warmup 20,000, with adaptation from proposal 5,000;
  - 60,000 evaluations per chain;
  - draws stored every 10th iteration (`--thin`);
  - seed 20261101;
  - fit `output/calibration/r_20261030_agepop/fit.rds` (used only to check that the model matches);
  - output to `output/calibration/mcmc_reference/`.
- **Sampler:** adaptive random-walk Metropolis (`arwm_budget` in `scripts/calibration/mcmc.R`). The
  initial proposal is `0.1² × diag(prior variance)`; adaptation uses the recent half of warmup
  only. Starts are truncated-prior draws.
- **Result (2026-10-30, 1.3 hours on 8 cores):**
  - R-hat ≤ 1.004;
  - bulk ESS 2,137–3,474 and tail ESS ≥ 2,219;
  - MCSE of the mean ≤ 0.022 posterior SD;
  - acceptance 0.12–0.13 (below the 0.15 flag; mixing is still adequate by ESS).
- **Posterior:** `c` median 2.59 (95% interval 1.45–5.24) against a prior median of about 1.05;
  posterior SD 0.78 × the truncated prior SD. Entry scales <20 to 50–59 have medians 1.09–1.17
  (upper 95% limits 1.41–1.61); 60+ is wider (median 1.50, up to 3.97).
- **Prior-dominated parameters** (posterior SD above 0.7 × the truncated prior SD): θ5 and θ6
  (50–59 and 60+ contact scales) and θ13 (`c`, 0.78; its location still moves well above the
  prior). θ12 (60+ entry) is at 0.62.
- `--targets` swaps the target file. The MCMC SBC uses it.

## 3. Optional: direct MCMC SBC

```bash
Rscript scripts/sbc_mcmc.R --n-datasets 8 --cores 8   # 15,000 warmup, 45,000 evaluations; about 1.8 hours
```

- The truths are the first 8 random prior draws (`prior_sampling.R`) that pass the equilibrium
  check (`theta_true.csv`), with synthetic targets drawn from the three likelihood terms. Each dataset gets one 45,000-evaluation chain.
- **Result (2026-10-30, 1.8 h):** pooled cov95 0.923 and cov90 0.894 over 104 ranks; mean rank 0.490; max split
  R-hat 1.011 (one chain per dataset). Lowest per-parameter cov95: θ13 (`c`), 6 of 8.
- **Limitation:** 8 datasets are a check against gross miscalibration only, not a full SBC.

## 4. 2015 predictions and calibration figure

```bash
Rscript scripts/predict_2015_posterior.R --max-draws 4000 --cores 14
Rscript scripts/plot_calibration.R --cores 14 --max-draws 2000
```

- **Prediction outputs:** `output/calibration/mcmc_reference/model_prediction_2015_strata.*`.
  Ever-PWID 50,494 (95% interval: 45,344–58,894); D 1,040 (490–1,903); F 5,934 (5,353–6,886);
  D + F 7,005 (6,125–8,319) against the Teo target 11,000; J 8,493 (7,687–9,783). The prison
  size is not fitted: the observed 2015 prison count of screened PWID is 3,367
  (`output/prison_2015_counts.csv`).
- **Calibration figure:** `figures/calibration_agepop.png` has three panels (seroprevalence,
  first arrests, community PWID). Each target is shown with the 95% range of its likelihood,
  against 95% posterior intervals of the model values (no observation noise) from the reference MCMC.
- **Overlap check:** `figures/calibration_interval_overlap.csv` records whether each 95% interval
  overlaps its target (MCMC: 13 of 13). The posterior medians of first
  arrests sit above the CNB counts at ages 30+, but within the negative binomial range.

## 5. 2026 treatment scenarios

```bash
SC="--root . --draws-per-seed 2000 --cores 8 --sample-seed 20261005"   # ~7.5 min per run on 8 cores
Rscript scripts/scenarios_2026.R $SC                                    # main: J from 2026, coverage 99.995%
Rscript scripts/scenarios_2026.R $SC --j-probability 0   --out-dir output/scenarios_2026_jcov/j000
Rscript scripts/scenarios_2026.R $SC --j-probability 0.5 --out-dir output/scenarios_2026_jcov/j050
Rscript scripts/scenarios_2026.R $SC --j-probability 0.9 --out-dir output/scenarios_2026_jcov/j090
Rscript scripts/scenarios_2026.R $SC --j-start-year 2019 --out-dir output/scenarios_2026_j2019   # sensitivity
```

- **Draws:** 2,000 draws from the reference MCMC posterior (the default `--posterior-files`). The full
  analysis is run per draw, then summarised as the median and 95% UI.
- **Design (user decision 2026-10-30):** main analysis = no treatment before 2026 (the model stays at
  the 2015 equilibrium; 2016–2025 within 0.05% of 2015), J treatment and F/X screening from 2026.
  `--j-start-year` (default `j_start_year` in `data/scenario_policy_config.json`, 2026) and
  `--j-probability` (default 0.99995) set the start and annual start probability of J treatment.
- **Options:**
  - `--posterior-files` accepts any CSV with `theta1..theta13`;
  - `--figures true` also draws the old ggplot figures, into `output/scenarios_2026/legacy/`
    (`scripts/paper_figures.R` writes its legacy versions there too);
  - `--min-valid-share` sets the share of draws that must pass equilibrium.
- **Results (2026-10-30):** main J only: incidence −49.8% (37.7–56.9), mortality −32.6%. Five finite
  policies meet incidence (8,000/yr any F share; 5,000/yr F50%/F75%); none meets mortality (best
  8,000/yr F25% 62.6%, joint 5.7%); only no cap meets both. Coverage grid: every reduction rises with
  J coverage; no finite policy meets mortality at any coverage. Sensitivity (J from 2019): J only
  −83.4% / −53.2%; mortality met at 5,000/yr F25% (66.9%, joint 90%), 8,000/yr F25% and F50%.
- **Figures:** after the runs, `python scripts/figures.py .` (main, includes `jcov_heatmap`), then
  `python scripts/figures.py . --scen-dir output/scenarios_2026_j2019 --only scenario_target_trajectories paper_fig2_trajectories`.
- **More detail:** see `output/scenarios_2026/README.md`.

## 6. Manuscript

```bash
python scripts/figures.py .   # redraws every figure (PNG + editable SVG)
python3 scripts/build_paper.py . --python /path/to/python-with-pillow
```

- **`build_paper.py`** writes the PLOS Medicine files: `docs/main.docx`, `docs/Supplementary
  Information.docx` (S1 Text), `docs/S1_Checklist.docx` (CHEERS 2022) and `docs/paper/figures/Fig1–4.tif`.
  The docx files are double spaced with continuous line numbers and page numbers.
- **Prose and numbers:** text is in `scripts/paper/main_text.py` and `si_text.py`; every number comes
  from `scripts/paper/facts.py`, which reads `output/` and `figures/` (no simulation). The build stops
  if a result the text relies on changes (`main_text.checks`).
- **Evidence registries** (scientific-writing skill format) in `docs/paper/`: each factual line has a
  claim ID bound to evidence IDs (`claims_*.csv`, `source_manifest.json`, `*_audit.md`). All are
  `unverified` until an author checks them. If the skill is installed, the build runs its local audits;
  the expected failures are the unverified claims and the missing statements.
  `docs/paper/unresolved_issues.md` lists what the authors must supply or decide.

`figures.py` redraws all figures from the saved outputs with the `scientific-visualization` skill
(`~/.claude/skills/scientific-visualization`, or pass `--skill-dir`). It needs Python ≥ 3.11 with
`matplotlib==3.11.2 pillow==12.3.0 pandas scipy`. `--only <name> ...` redraws selected figures.
Run it after any step above that changes its inputs. The R scripts still draw their own ggplot
versions of the same files, so run `figures.py` last.
