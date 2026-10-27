# Training

Production calibration path:
1. **Nelder–Mead multi-start** finds the posterior mode (fixed prior + seroprevalence likelihood).
2. **Reference MCMC** gives the production posterior.
3. **Scenarios** run on draws from that posterior.

Cross-validation only (nothing here feeds the production path): single-round NPE trained on
50,000 Latin hypercube prior simulations, with simulation-based calibration (SBC), and a small
direct MCMC SBC.

Outputs are in `output/calibration/` (see its `README.md`). Run the steps in this order; each step
reads the previous step's output. Model design: `docs/first_arrest_multiplier.md`, calibration
redesign (entry scales, first-arrest and community-PWID targets): `docs/redesign_prompt.md`.

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
  θ7–θ12 are log entry scales `s_i` (entry into D_i = `s_i` × CNB 2015 new drug abusers, held
  constant in every year); θ13 is log `c`, the first-arrest multiplier (`build_params()` in
  `scripts/calibration/likelihood.R`). Design: `docs/redesign_prompt.md`.
- **Bounds:** contact scales 0.01–1,000; entry scales 1–100 (user decision: the CNB counts are a
  lower bound); `c` 0.5–100 (user decision: first arrest at least half the re-arrest rate).
- **Prior (fixed, not centred on any fit):** independent Normal, truncated to the bounds:
  θ1–θ6 ~ N(0, 1.5), θ7–θ12 ~ N(0, 1) truncated at 0 (half-normal), θ13 ~ N(0, 0.5)
  (`prior_spec()`; density `nlp()`, shared by Nelder–Mead, MCMC and NPE).
- **Likelihood (`nll_targets()`), three independent terms** (`data/calibration_targets.json`):
  - seroprevalence: each age group's observed 2015 prison seroprevalence is
    `Beta(p·φ, (1−p)·φ)` around the model value `p`, φ = 100 (`beta_precision`);
  - first arrests: CNB 2015 new drug abusers by age (257, 644, 292, 94, 21, 1) are negative
    binomial around the model D → J flow `c·λ3_i·D_i`, size 50 (`first_arrests_nb_size`);
  - community PWID: Teo et al. 2019 (network scale-up, 2017 survey) 11,000 (6,500–17,000)
    injecting drug users in Singapore; log-normal around the model D + F, log SD 0.245
    (`community_pwid_log_sd`). The model value is the 2015 equilibrium, which also stands for
    2017 (constant entry, no treatment before 2019).
  - Seroprevalence is unchanged when all entry scales rise together (transmission depends on the
    infectious share), so entry is pinned by the two count targets. `c` is pinned by the tension
    between them: a low `c` keeps more people in D (raising D + F) but shifts first arrests to
    older ages.
- **Support:** draws outside the bounds and equilibrium failures are excluded (MCMC gives −∞;
  NPE drops the simulation and counts it).
- **Solver:** RK4 with `dt = 1/365`.
- **Mortality:** general-population `mu` × PWID SMR by state: 9.52 viraemic (a, c, t), 4.73
  otherwise (van Santen 2014). In DC and HCC the liver excess `mu_DC`/`mu_HCC` is added, and after
  SVR only that liver excess is multiplied by ψ (0.5 DC, 0.37 HCC).
- **Details:** see `output/calibration/npe_vs_mcmc/model_audit.md`.

## 1. Deterministic calibration (posterior mode)

```bash
Rscript scripts/run_calibration.R --run-id r_20261027_agepop \
  --maxit 5000 --nm-restarts 5 --n-starts 8 --sd-perturb 0.9 --seed 20261001 --cores 8   # ~2 h
```

- **Objective:** negative log posterior (`nll_targets()` + prior); every evaluation runs model
  time 0–150 and checks equilibrium over the last 10 years.
- **Starts (uninformative):** θ = 0 (contact scales 1, entry scales 1, `c` = 1) and 7 Normal perturbations of it
  (SD 0.9), clipped to the bounds, run in parallel.
- **Restarts:** each start reruns Nelder–Mead from its own solution while the objective still
  falls by more than 0.001 (`--nm-restarts`, at most 5 restarts).
- **Result (2026-10-27):** 7 of 8 starts reach the same mode (objective 18.80–19.53; `perturb_4`
  stops in a local optimum at 121.5). Best: seroprevalence RMSE 0.014, largest error 0.023;
  `c` 2.28; D + F 7,720 (Teo 11,000, 6,500–17,000).
  - Entry scales sit at their lower bound 1 for ages 30+ and at about 1.1 for <20 and 20–29.
  - First arrests are too high at older ages (30–39 324 vs 292; 40–49 113 vs 94; 50–59 28 vs 21;
    60+ 3.5 vs 1): people who start injecting young age into the older groups before their first
    arrest, and entry cannot fall below the CNB counts.
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
  - fit `output/calibration/r_20261027_agepop/fit.rds` (used only to check that the model matches);
  - output to `output/calibration/mcmc_reference/`.
- **Sampler:** adaptive random-walk Metropolis (`arwm_budget` in `scripts/calibration/mcmc.R`). The
  initial proposal is `0.1² × diag(prior variance)`; adaptation uses the recent half of warmup
  only. Starts are truncated-prior draws.
- **Result (2026-10-27, 1.9 hours on 8 cores):**
  - R-hat ≤ 1.004;
  - bulk ESS 2,158–3,425 and tail ESS ≥ 1,936;
  - MCSE of the mean ≤ 0.023 posterior SD;
  - acceptance 0.12–0.14 (below the 0.15 flag; mixing is still adequate by ESS).
- **Posterior:** `c` median 2.78 (95% interval 1.55–5.49) against a prior median of about 1.05;
  posterior SD 0.75 × the truncated prior SD. Entry
  scales <20 to 50–59 have medians 1.09–1.15 (upper 95% limits 1.41–1.58); 60+ is wider
  (median 1.47, up to 3.77).
- **Prior-dominated parameters** (posterior SD above 0.7 × the truncated prior SD): θ5 and θ6
  (50–59 and 60+ contact scales) and θ13 (`c`, 0.75; its location still moves well above the
  prior). θ12 (60+ entry) is at 0.60.
- `--targets` swaps the target file. The MCMC SBC uses it.

## 3. Cross-validation: single-round NPE with SBC

```bash
Rscript scripts/run_npe.R --step all --n-sims 50000 --seed 20261012 --train-seed 20261013 \
  --n-sbc 1000 --n-cores 10 --threads 4 --python /path/to/python
```

- **LHS design** (`scripts/calibration/prior_sampling.R`): a 50,000 × 13 design with one point per
  stratum, mapped through the truncated Normal inverse CDF of each parameter's prior.
- **Simulations:** draws outside the bounds, non-finite results and equilibrium failures are dropped
  and counted, never replaced. This run (2026-10-27): 49,941 valid, 47 equilibrium failures and 12
  non-finite (`simulation_summary.json`).
- **Summaries (13):** 6 seroprevalences with Beta(φ = 100) noise (logit-transformed), 6 first
  arrests with negative binomial (size 50) noise (`log1p`), and D + F with log-normal (SD 0.245)
  noise (log).
- **Network:** an equal-weight ensemble of 5 MAFs (50×5), each with batch 256, Adam lr 5e-4,
  validation fraction 0.10, patience 20, ≤ 250 epochs, gradient clipping 5.0, CPU only. Posterior
  draws outside the bounds are rejected (at most 100 rounds, then an error). Training takes
  about 10 minutes (152–236 epochs per member).
- **Output:** 20,000 posterior draws.
- **SBC:** 1,000 independent prior datasets (998 pass equilibrium), with ranks of the true value in
  1,000 posterior draws.
  - **Gate:** cov95 in [0.85, 1] and cov90 in [0.75, 1] for every parameter, and no systematically
    skewed rank histogram.
  - **Result (2026-10-27):** cov95 0.926–0.962 and cov90 0.871–0.917 (all pass), but three
    rank shapes fail: θ2 (20–29 contact, mean rank 0.463, 4.0 SE low), θ11 (50–59 entry, 0.456,
    4.8 SE low) and θ10 (40–49 entry, outer-10% share 0.129, 3.1 SE high). Gate FAIL
    (`NPE_NOT_OK.txt`). The entry scales have posteriors piled at their lower bound, which a
    flow with rejection sampling fits poorly. The NPE is cross-validation only, so it was not
    retrained.
- **Comparison with the reference MCMC (2026-10-27):** NPE draws that fail the equilibrium check
  are dropped first (945 of 4,000, 24%). At the observed data the NPE does not match the MCMC
  (2 of 13 parameters pass all checks; 95% intervals on average 28% wider, range 0.92–2.25×;
  worst mean difference 1.08 SD). Policy medians agree within 5.4% (J-only incidence) and 1.3%
  for mortality; the joint-pass probability at 5,000/yr F25% differs by 34 pp and cost by 9%. See
  `output/calibration/npe_vs_mcmc/report.md`.

  ```bash
  Rscript scripts/compare_npe_mcmc.R --policy-cores 7
  ```

## 4. Optional: direct MCMC SBC

```bash
Rscript scripts/sbc_mcmc.R --n-datasets 8 --cores 8   # 15,000 warmup, 45,000 evaluations; about 1.8 hours
```

- The truths are the first 8 NPE SBC datasets, with synthetic targets drawn from the three
  likelihood terms. Each dataset gets one 45,000-evaluation chain.
- **Result (2026-10-27):** pooled cov95 0.952 and cov90 0.942 over 104 ranks; mean rank 0.467; max split
  R-hat 1.015 (one chain per dataset). Lowest per-parameter cov95: θ6, 6 of 8.
- **Limitation:** 8 datasets are a check against gross miscalibration only, not a full SBC.

## 5. 2015 predictions and calibration figure

```bash
Rscript scripts/predict_2015_posterior.R --max-draws 4000 --cores 14
Rscript scripts/plot_calibration.R --cores 14 --max-draws 2000
```

- **Prediction outputs:** `output/calibration/mcmc_reference/model_prediction_2015_strata.*`.
  Ever-PWID 57,742 (95% interval: 52,350–66,850); D 1,130 (577–2,004); F 6,802 (6,189–7,827);
  D + F 7,957 (7,062–9,298) against the Teo target 11,000; J 9,741 (8,889–11,104). The prison
  size is not fitted: the observed 2015 prison count of screened PWID is 3,367
  (`output/prison_2015_counts.csv`).
- **Calibration figure:** `figures/calibration_agepop.png` has three panels (seroprevalence,
  first arrests, community PWID). Each target is shown with the 95% range of its likelihood,
  against 95% posterior intervals of the model values (no observation noise) from the reference MCMC and
  the NPE (LHS; 1,527 of 2,000 NPE draws pass equilibrium).
- **Overlap check:** `figures/calibration_interval_overlap.csv` records whether each 95% interval
  overlaps its target (MCMC: 13 of 13; NPE: 13 of 13). The posterior medians of first
  arrests sit above the CNB counts at ages 30+, but within the negative binomial range.

## 6. 2026 treatment scenarios

```bash
Rscript scripts/scenarios_2026.R --root . --draws-per-seed 2000 --cores 4 --sample-seed 20261005
```

- **Draws:** 2,000 draws from the reference MCMC posterior (the default `--posterior-files`). The full
  analysis is run per draw, then summarised as the median and 95% UI.
- **Options:**
  - `--posterior-files` accepts any CSV with `theta1..theta13`;
  - `--figures true` also draws the old ggplot figures, into `output/scenarios_2026/legacy/`
    (`scripts/paper_figures.R` writes its legacy versions there too);
  - `--min-valid-share` sets the share of draws that must pass equilibrium.
- **More detail:** see `output/scenarios_2026/README.md`.

## 7. Manuscript

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
