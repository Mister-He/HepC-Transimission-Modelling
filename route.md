# Learning route: Singapore PWID HCV model

This guide walks you through the repository in the order the analysis was built:
concepts → simulator → calibration → posterior → scenarios → manuscript. Every
stage gives the command that reproduces its results, a short code example, the
location of the original code (`file:line`), the inputs it reads, the outputs it
writes and how long it takes. Line numbers refer to the files as of 2026-10-06.

> **Rule before you run anything:** most scripts overwrite their result folders
> (see `AGENTS.md`). To experiment, pass a different `--out-dir` / `--run-id`,
> and delete it afterwards.

---

## 0. Map of the repository

| Path | What it is | Guide |
| --- | --- | --- |
| `AGENTS.md`, `docs/AGENTS.md` | Working rules and the model contract (strata, states, targets, gates) | — |
| `Hepc_modelling_technical_details.pdf` | Assumptions, equations, parameter tables, schematic | — |
| `docs/equations.tex` | Full equation reference (480 ODEs) | — |
| `docs/model_card.md`, `docs/training.md`, `docs/evaluation.md` | Purpose, run commands, current calibration results | — |
| `data/` | Inputs: targets, inflow series, policy and cost config | `data/README.md` |
| `scripts/model/` | `sim.cpp` (simulator), `setup.R` (parameters, initial state), rate-fit code | `scripts/model/README.md` |
| `scripts/calibration/` | Targets loader, summaries, likelihood, equilibrium check, Nelder–Mead, MCMC | `scripts/calibration/README.md` |
| `scripts/*.R`, `scripts/*.py` | Entry points for each pipeline stage | `scripts/README.md` |
| `output/calibration/` | Fits, NPE posteriors, MCMC samples (active: `r_20260929_agepop`, `npe_agepop_recal`, `mcmc_recal`) | `output/calibration/README.md` |
| `output/scenarios_2026/` | 2026–2030 policy results and figures | `output/scenarios_2026/README.md` |
| `figures/` | Calibration figures | `figures/README.md` |
| `docs/main.docx`, `docs/Supplementary Information.docx` | Manuscript and S1 Text | — |
| `tests/` | Unit and integration tests | `tests/README.md` |
| `claude.md`, `codex.md` | Hand-off notes from previous agent sessions | — |

The pipeline and its dependencies:

```text
data/sps.csv ──► preprocess.R ──► rate fits (informational; values hard-coded in setup.R)
                                        │
data/calibration_targets.json ──┐       ▼
data/singstat_new_drug_abusers ─┼─► run_calibration.R ──► output/calibration/r_20260929_agepop/fit.rds
                                │                               │
                                │        ┌──────────────────────┤
                                │        ▼                      ▼
                                │   run_npe.R (+npe_train.py)   run_mcmc_validation.R
                                │        │  npe_agepop_recal/   │  mcmc_recal/
                                │        ├──► predict_2015_posterior.R
                                │        ├──► plot_calibration.R ──► figures/calibration_agepop.png
                                │        ▼
data/scenario_policy_config ────┴─► scenarios_2026.R ──► output/scenarios_2026/*.csv, *.png
                                         ▼
                                    paper_figures.R ──► output/scenarios_2026/paper_fig*.png
                                         ▼
                                    build_paper.py ──► docs/main.docx, docs/Supplementary Information.docx
```

---

## Stage 1. Concepts (≈ 2 hours, reading only)

Read in this order:

1. `docs/model_card.md`: purpose and limits, one page.
2. `docs/AGENTS.md`: the model contract. Memorise these points:
   - Groups `D` (active, never imprisoned), `J` (in prison), `F` (released, still injecting) and `X` (released, stopped).
   - Liver stages `NC, CC, DC, HCC`.
   - HCV states `u, a, c, t, r`.
   - Six age groups.
   - 4 × 4 × 5 × 6 = **480 compartments**.
   - Transmission only in D and F.
   - Model year 150 = 2015.
3. `Hepc_modelling_technical_details.pdf`: sections 1 to 3 (assumptions, equations) and tables 1 and 2 (parameters).
4. `docs/evaluation.md`: what "good calibration" means here.
5. `output/scenarios_2026/README.md` (Chinese): the policy questions and results.

**Key idea: one flat index.** Every compartment is one element of a 480-vector:

```cpp
// scripts/model/sim.cpp:15
inline int idx(int s, int k, int h, int i) {   // s: D,J,F,X = 0..3; k: stage 1..4
    return s*4*5*6 + (k-1)*5*6 + h*6 + i;      // h: u,a,c,t,r = 0..4; i: age 0..5
}
```

R uses the same formula plus 2, because column 1 of the output matrix is time. See
`scripts/model/setup.R:18` and `scripts/calibration/model_metrics.R:2`.

---

## Stage 2. Environment (≈ 1 hour)

Versions used for the current results:

| Tool | Version | Needed for |
| --- | --- | --- |
| R | 4.5.1 (CI uses 4.4.1) | everything |
| R packages | Rcpp, RcppArmadillo, dplyr, jsonlite, ggplot2, tidyr, patchwork, scales, viridisLite, ragg, writexl, readxl, MASS | model, calibration, scenarios, figures |
| R packages (optional) | lubridate, cmdstanr + CmdStan 2.39 | rate fitting (Stage 4) |
| Python | conda env `bayes-venv`: sbi 0.27.0, torch 2.14.0, numpy | NPE (Stage 6) |
| pandoc | 3.8 | manuscript (Stage 10) |

```bash
Rscript -e 'install.packages(c("Rcpp","RcppArmadillo","dplyr","jsonlite","ggplot2","tidyr",
  "patchwork","scales","viridisLite","ragg","writexl","readxl","MASS"))'
conda create -n bayes-venv python=3.11 && conda activate bayes-venv && pip install sbi torch numpy
brew install pandoc
```

An optional container for the R part is in `Dockerfile` and `docker-compose.yml`:
`docker compose run --rm r-test`. CI is in `.github/workflows/ci.yml`.

**Check the setup** (≈ 5 s). The `sim.cpp` integration test compiles the simulator:

```bash
Rscript scripts/run_tests.R --scope all      # expect: All 6 test file(s) passed.
```

---

## Stage 3. Inputs in `data/`

| File | Content | Used by |
| --- | --- | --- |
| `calibration_targets.json` | 2015 prison seroprevalence and population in six age groups, `population_sigma_log = 0.10`, horizon 150, lag 5 | calibration, NPE, MCMC, scenarios |
| `singstat_new_drug_abusers.csv` | New drug users by age, 2015–2025 ([data.gov.sg](https://data.gov.sg/datasets/d_f174679381e9085e32869c32bdb01b87/view)) | `setup.R:205` (inflow schedule) |
| `scenario_policy_config.json` | Screening totals, F shares, start probabilities, unit costs, discounting | `scenarios_2026.R:10` |
| `sps.csv` (**not in this worktree**) | Individual Singapore Prison Service records, sensitive and gitignored | `scripts/model/preprocess.R:9` only |

> **Pitfall.** `data/` is gitignored except the whitelisted files.
> - **Missing folder:** a fresh worktree may have no `data/` at all. In this
>   worktree the three JSON/CSV files were rebuilt on 2026-10-05 (see
>   `data/README.md`). The rebuild was verified by reproducing the earlier
>   scenario results exactly.
> - **Where `sps.csv` is:** only in the main checkout
>   (`../HepC-Transimission-Modelling/data/`).
> - **Stale targets:** that main checkout holds an *older* `calibration_targets.json`.
>   Do not copy it here.

Reading the targets:

```r
# scripts/calibration/targets.R:2
source("scripts/calibration/targets.R")
targets <- read_targets("data/calibration_targets.json")
targets$x_prev        # binomial counts = round(seroprevalence * population)
```

---

## Stage 4 (optional). Arrest, release and recidivism rates

- **Code:** `scripts/model/preprocess.R`, run from the repo root. `scripts/preprocess.R` is an
  older copy that uses `../` paths. Stan models are in `scripts/model/release_rate.stan`
  and `scripts/model/rearrest_rate.stan`. The binaries `scripts/release_rate` and
  `scripts/rearrest_rate` are compiled CmdStan executables for macOS arm64. Rebuild
  them; do not reuse them.
- **What it does:**
  - Cleans `sps.csv` into `output/cleaned_data.rds` (line 138).
  - Prints the age-specific prison counts and seroprevalence used as targets (lines 150–175).
  - Fits release rates (`release_rate.stan`, saved at line 253) and re-arrest and
    recidivism (`rearrest_rate.stan`).
- **Run:** the file is an interactive script. Run it section by section in RStudio with
  the working directory at the repo root. It needs `data/sps.csv` and `cmdstanr`.
- **Outputs:** `output/cleaned_data.rds`, `output/release_rate_fit.rds`, `output/rearrest_rate_fit.rds`.

> **Important.** The simulator does **not** read these fits. The rates are hard-coded:
> - **Arrest and release rates:** `lambda1`, `lambda2` and `lambda3` at `setup.R:150-153`,
>   commented "GUESS".
> - **Recidivism:** `pi_recid = 0.59` at `setup.R:154`. The header comment
>   (`setup.R:34`) and the technical PDF say 0.65, but **0.59 is what runs**.
>
> If you change these values, every later stage must be rerun.

---

## Stage 5. The simulator (≈ 3 hours with exercises)

### 5.1 Code tour

| Location | What to read |
| --- | --- |
| `scripts/model/sim.cpp:20` | `struct Params`, all parameters the C++ side expects |
| `sim.cpp:43` | `inflow_rate`: before model year 150 uses `beta`; after it uses the SingStat schedule × calibrated scale |
| `sim.cpp:54` | `force_of_infection`: infectious (a + c) over active population, D and F only |
| `sim.cpp:76` | `rhs`: the 480 derivatives (match against `docs/equations.tex`) |
| `sim.cpp:266` | `treatment_events`: integrates starts, completions and cures per stratum |
| `sim.cpp:300` | `run_sim(params, data)`: exported to R; returns time + 480 states (481 columns). If `data$treatment_capacity` (length 4) is supplied, 12 cumulative event columns are added (493 columns: starts, completions and cures for D, J, F, X) and the per-stratum caps are enforced (`sim.cpp:366-406`) |
| `scripts/model/setup.R:41` | `params` list, with literature sources in comments (lines 20–38) |
| `setup.R:116-122` | background mortality `mu` and the PWID multiplier `omega` |
| `setup.R:158`, `setup.R:172` | contact matrix `C_contact` and inflow `beta` (scaled by calibration) |
| `setup.R:190` | initial state `y0` (D/u and D/a seeds; all `r` start at zero) |
| `setup.R:215` | `data` list: `t_start`, `t_end`, `dt = 1/365`, inflow schedule starting at model year 150 |

### 5.2 Example 1: run one 150-year simulation (0.24 s)

Save as `ex1.R` and run `Rscript ex1.R` from the repo root. This is the pattern every entry
script uses. `setup.R` must be sourced from inside `scripts/model/` because it reads the
inflow CSV by a relative path.

```r
setup_env <- new.env()
old <- setwd("scripts/model"); sys.source("setup.R", setup_env); setwd(old)
Rcpp::sourceCpp("scripts/model/sim.cpp")                 # compiles once per session
source("scripts/calibration/targets.R")
targets <- read_targets("data/calibration_targets.json")  # needed by the helpers below
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R"))
  source(file.path("scripts/calibration", f))

fit <- readRDS("output/calibration/r_20260929_agepop/fit.rds")
pm  <- build_params(fit$best$theta, fit$base_params)     # likelihood.R:32
d   <- setup_env$data; d$t_start <- 0; d$t_end <- 150
out <- run_sim(pm, d)                                    # 54,751 × 481 matrix
J_summary_final(out)                                     # model_metrics.R:24
check_equilibrium(out)$pass                              # equilibrium.R:8 → TRUE
```

Expected output: the J-stratum summary equals `output/calibration/r_20260929_agepop/predictions.csv`.
For example, age `<20` gives `N_hat = 103.3` and `p_sero = 0.131`.

### 5.3 Exercises

1. **Run without inflow.** Set `d$inflow_schedule <- NULL` and run to year 165.
   Compare total population after 2015.
2. **Transmission and prevalence.** Double `pm$q` and watch J seroprevalence change.
   This is why contact scales and q are not separately identifiable.
3. **Treatment in prison.** Treatment rates come from `pm$tau_by_stratum`, a length-16
   stratum-major vector (`sim.cpp:109-112`), for ages `tau_min_age` to `tau_max_age`. Start
   from the 2015 state (`d$y0 <- out[nrow(out), 2:481]`, `d$t_start <- 150`, `d$t_end <- 160`),
   set `pm$tau_by_stratum <- rep(0, 16)`, `pm$tau_by_stratum[5:8] <- -log1p(-0.99995)` (J) and
   `d$treatment_capacity <- rep(Inf, 4)`. Then read cumulative J starts in column 483.
   `scenarios_2026.R:84` (`run_year`) does exactly this.

---

## Stage 6. Deterministic calibration (≈ 30–60 min)

- **Code:**
  - `scripts/run_calibration.R:32-64` (starts, objective, best run).
  - `scripts/calibration/likelihood.R:32` (`build_params`: θ1–6 scale the contact rows,
    θ7–12 scale `beta`).
  - `likelihood.R:42-50` (binomial and log-normal NLL).
  - `likelihood.R:52` (`make_objective`, which adds the plausibility and equilibrium penalties).
  - `scripts/calibration/calibrate_nm.R:24,54,61` (Nelder–Mead wrapper, anchors, start set).
- **Command** (reproduces the active fit; overwrites it):

```bash
Rscript scripts/run_calibration.R \
  --run-id r_20260929_agepop --maxit 1000 --n-starts 5 --seed 20261001
```

- **Smoke test first** (≈ 1 min, separate folder):

```bash
Rscript scripts/run_calibration.R --run-id smoke_test --n-starts 1 --maxit 50
rm -r output/calibration/smoke_test
```

- **Outputs:** in `output/calibration/<run-id>/`: `fit.rds` (the input to every later
  stage), `solutions.csv`, `predictions.csv`, `equilibrium.csv`,
  `optimization_history.csv`, `initial_values.csv`, `run_config.json` and a copy of the targets.
- **Expected:** best objective **19.5666**, prevalence RMSE **0.0087**, population MAPE
  **3.77%**, equilibrium pass (`run_config.json`, `docs/evaluation.md`).
- **Runtime:** each objective call is a 150-year run (~0.25 s). There are 7 starts × ≤1000 iterations.

---

## Stage 7. Neural posterior estimation (NPE) (≈ 30–60 min on CPU)

- **Code:**
  - `scripts/run_npe.R:33` `sample_prior`: Normal on log scale, centred on the best fit.
  - `run_npe.R:42` `simulate`: one model run returns 6 seroprevalences + 6 population sizes.
  - `run_npe.R:67` `step_data`: 3,000 training + 300 SBC simulations.
  - `run_npe.R:79` `py_round`: calls Python.
  - `run_npe.R:109` `step_train`: seed `SEED` and `SEED+1`, suffix `seed2`.
  - `run_npe.R:118` `step_evaluate`: the quality gate. SBC 95% coverage must be 0.85–1,
    90% coverage 0.75–1, and every posterior SD at least 0.05.
  - `scripts/npe_train.py:56` (observation noise), `:80` (prior), `:178` (`NPE_C` with a MAF),
    `:203` (writes posterior samples).
- **Command** (`docs/training.md`):

```bash
Rscript scripts/run_npe.R --step all \
  --fit output/calibration/r_20260929_agepop/fit.rds \
  --out-dir output/calibration/npe_agepop_recal \
  --n-sims 3000 --n-draws 5000 --n-cores 4 --n-rounds 1 \
  --max-epochs 500 --seed 20261001 \
  --npe-contact-sd 0.20 --npe-beta-sd 0.20 \
  --npe-hidden-features 128 --npe-num-transforms 10 \
  --python /opt/homebrew/Caskroom/miniconda/base/envs/bayes-venv/bin/python
```

  Steps can be run separately: `--step data`, then `--step train`, then `--step evaluate`.
- **Outputs** (`output/calibration/npe_agepop_recal/`):
  - Training and test data: `theta_train.csv`, `x_train.csv`, `theta_sbc.csv`, `x_sbc.csv`, `x_obs.csv`.
  - Posterior draws: `posterior_samples_npe.csv` and `posterior_samples_npe_seed2.csv`, 5,000 each.
  - Diagnostics and model: `sbc_summary*.csv`, `npe_evaluation*.csv`, `npe_inference_round1*.pkl`, `npe_config*.json`.
  - **The gate file:** `NPE_OK.txt` or `NPE_NOT_OK.txt`.
- **Expected:** `NPE_OK.txt`. Minimum 95% SBC coverage is 0.927 in both seeds, and minimum
  posterior SD is 0.0549 and 0.0544. `tests/unit/test_npe_coverage.R` checks this.
- **Reproducibility note:** neural network training is not bit-for-bit reproducible across
  machines or library versions. A rerun gives a statistically similar posterior but
  different draws, and every later number shifts slightly. To reproduce the published
  numbers exactly, keep the existing posterior CSVs.
- **Historical runs** (not active): `npe_agepop/`, `npe_agepop_two_round_300/`,
  `npe_agepop_two_round_600/`. Each has its own README.

Example 2: inspect the posterior.

```r
p <- rbind(read.csv("output/calibration/npe_agepop_recal/posterior_samples_npe.csv"),
           read.csv("output/calibration/npe_agepop_recal/posterior_samples_npe_seed2.csv"))
t(apply(exp(p), 2, quantile, c(.025, .5, .975)))   # multiplicative scales on contact / inflow
```

---

## Stage 8. MCMC cross-check (several hours)

- **Code:** `scripts/run_mcmc_validation.R:67-69` (NPE covariance as the proposal start),
  `:84` (`am_mcmc` calls) and `:132` (strict gate). `scripts/calibration/mcmc.R:22`
  (priors), `:39` (`log_posterior`), `:76` (adaptive Metropolis), `:151` (summaries).
- **Command:**

```bash
Rscript scripts/run_mcmc_validation.R \
  --fit output/calibration/r_20260929_agepop/fit.rds \
  --npe-dir output/calibration/npe_agepop_recal \
  --out-dir output/calibration/mcmc_recal \
  --n-iter 120000 --burnin 20000 --thin 1 \
  --n-chains 4 --parallel-chains 4 --seed 20261002 \
  --contact-sd 0.20 --beta-sd 0.20 --dt 0.0833 \
  --ess-min 4000 --rhat-max 1.001 --same-start true
```

  `--dt 0.0833` uses a monthly step for speed. For a quick look use
  `--n-iter 4000 --burnin 1000` with a separate `--out-dir`.
- **Outputs:** `posterior_samples_mcmc.csv` (≈ 90 MB), `mcmc_parameter_summary.csv`,
  `mcmc_diagnostics.csv`, `npe_covariance_used.csv`, `MCMC_STRICT_GATE.txt` and `mcmc_run_config.json`.
- **Expected:** `PASS`. Minimum ESS is 10,302 and maximum R-hat 1.000552. Medians are
  within 0.02 (log scale) of the pooled NPE medians.

---

## Stage 9. 2015 predictions and calibration figures (≈ 10–20 min)

```bash
Rscript scripts/predict_2015_posterior.R \
  --fit output/calibration/r_20260929_agepop/fit.rds \
  --out-dir output/calibration/npe_agepop_recal --year 2015 --cores 4
Rscript scripts/plot_calibration.R \
  --fit output/calibration/r_20260929_agepop/fit.rds \
  --npe-dir output/calibration/npe_agepop_recal \
  --mcmc-dir output/calibration/mcmc_recal --cores 4 --max-draws 2000
```

| Script | What it does | Outputs |
| --- | --- | --- |
| `predict_2015_posterior.R` | Runs all 10,000 NPE draws to 2015 and summarises by stratum (`:82-114`) | `model_prediction_2015_*.csv/.md/.png` |
| `plot_calibration.R` | Draws 2,000 draws per method and compares intervals with targets (`:173`, `:194`) | `figures/calibration_agepop.png`, `figures/calibration_interval_overlap.csv`, predictive caches |

Expected values:
- **2015 populations** (`model_prediction_2015_strata.md`): D 1,386, J 3,340, F 1,513, X 14,680.
- **Ever-PWID total:** 20,915.
- **Calibration check:** all 12 interval overlaps are `TRUE`.

Run `plot_calibration.R` only after the MCMC run finishes.

---

## Stage 10. 2026–2030 policy scenarios (≈ 18 min on 14 cores)

### 10.1 What the script does (`scripts/scenarios_2026.R`)

| Line | Step |
| --- | --- |
| `:10` | Reads and checks `data/scenario_policy_config.json` |
| `:26-38` | Samples `--draws-per-seed` posterior rows from each NPE seed (`--sample-seed 20261005`) |
| `:43` | `rate_components`: annual new infections and HCV deaths, total and for F + X |
| `:84` | `run_year`: one calendar year (model time `150 + year − 2016` onward) with a treatment-capacity vector |
| `:94` | `j_positive_entries`: antibody tests on prison entry (cost accounting) |
| `:106` | `screen_plan`: the screening quota algorithm. Screened = min(quota, eligible); treatment cap = screened × c / N |
| `:120`, `:152` | `process_row` and `add_cost`: services, unit costs, 3% discounting |
| `:169-191` | `scenario_table`: J only + 15 finite (5 totals × 3 F shares) + no cap = 17 |
| `:195` | `evaluate_draw(theta)`: **the whole analysis for one parameter set** (equilibrium, history 2016–2025 with J from 2019, no-treatment counterfactual at `:253`, all 17 policies at `:233`) |
| `:349` | `agg_draws`: median and 2.5/97.5 percentiles across draws, plus `_prob` for logical columns |
| `:364` | `parallel::mclapply` over draws (forks; not available on Windows) |
| `:384-450` | Pass flags from medians, per-seed table, "xx (95% UI: xx–xx)" text, CSV/XLSX writing |
| `:455` onward | Scenario figures |

### 10.2 Commands

```bash
# quick test, separate folder (≈ 1 min); then delete it
Rscript scripts/scenarios_2026.R --root . --draws-per-seed 4 --cores 8 --out-dir _tmp_test
rm -r _tmp_test

# full reproduction (overwrites output/scenarios_2026/)
Rscript scripts/scenarios_2026.R --root . --draws-per-seed 1000 --cores 14
```

Run long jobs in a way that survives the terminal closing, for example `nohup ... &`. Check
that no older copy is still running (`pgrep -f scenarios_2026.R`). Two copies write to the
same folder.

### 10.3 Outputs (`output/scenarios_2026/`)

| File | Content |
| --- | --- |
| `scenario_2030_summary.csv` | One row per scenario. Every metric has median, `_lo`, `_hi`, plus `*_text` ("xx (95% UI: xx–xx)"), `joint_pass_probability` and cost sensitivity |
| `scenario_trajectories.csv` | 2015–2030 incidence and mortality ratios by scenario (median and UI) |
| `no_treatment_trajectory.csv` | No-treatment counterfactual 2015–2030 |
| `scenario_annual_services_costs.csv` | Per scenario × year × stratum services and costs, with cumulative columns |
| `scenario_treatment_summary.csv`, `target_scenario_treatment_summary.csv` | Starts, completions and cures. The second file has scenarios passing both targets |
| `scenario_results.xlsx` | All tables, including sheet `BySeed` |
| `scenario_*.png` | Scenario figures (medians only) |
| `calibration_source.txt` | Which posterior files, draws and sampling seed were used |

### 10.4 Expected key results (median, 95% UI)

| Quantity | Value |
| --- | --- |
| No treatment, 2030 incidence vs 2015 | 130.3% (121.6%–140.0%) |
| J only, incidence reduction 2030 | 75.9% (70.5%–79.8%) |
| J only, mortality reduction 2030 | 42.6% (41.8%–43.5%) |
| 5,000/yr, F 25%: mortality reduction / share of draws passing both | 65.1% (63.7%–66.4%) / 56% |
| 8,000/yr, F 25% vs F 75%: mortality reduction | 70.7% vs 57.4% (incidence 94.2% in both) |
| Finite scenarios meeting both targets at the median | 3 of 15 |

### 10.5 Example 3: read and query results

```r
s <- read.csv("output/scenarios_2026/scenario_2030_summary.csv")
s[, c("scenario_id", "incidence_reduction_text", "mortality_reduction_text", "joint_pass_probability")]
a <- read.csv("output/scenarios_2026/scenario_annual_services_costs.csv")
subset(a, scenario_id == "J_only" & year == 2026 & stratum %in% c("F","X"),
       c(stratum, screening_eligible, chronic_stock_start))   # F ≈ 103, X ≈ 1,620 chronic
```

### 10.6 Exercise: add a scenario

Edit `data/scenario_policy_config.json`. Note that `scenarios_2026.R:11-14` asserts the
current totals (600 to 8,000) and that `:191` asserts 17 scenarios. Relax those checks
first. Then run with a small `--draws-per-seed` and a separate `--out-dir`.

---

## Stage 11. Manuscript (≈ 1 min)

```bash
Rscript scripts/paper_figures.R .     # Fig 1–3 and S5 → output/scenarios_2026/paper_fig*.png
python3 scripts/build_paper.py .      # → docs/main.docx, docs/Supplementary Information.docx
```

- **`scripts/paper_figures.R`:**
  - Fig 1 schematic: `:14`.
  - Fig 2 trajectories (no treatment, J only, 8,000/yr at F 25% and F 75%): `:40`.
  - Fig 3 reductions by capacity and F share: `:69`.
  - Fig S5 quota use: `:95`.
- **`scripts/build_paper.py`:**
  - Helpers: `:24-70`.
  - Main-text template: `main = f"""…` at `:75`.
  - S1 Text template: `si = f"""…` at `:316`.
  - pandoc call: `:609`.
  - Every number is read from the CSVs, so edit prose **in this file**. Edits made in Word
    are lost on rebuild.
- **Citations:** `docs/references.bib` (BibTeX keys such as `[@chaillon2021]`), rendered with
  `docs/plos.csl` (PLOS Vancouver style).

---

## Stage 12. Tests

```bash
Rscript scripts/run_tests.R --scope unit          # targets, metrics, likelihood, equilibrium, NPE gate
Rscript scripts/run_tests.R --scope integration   # compiles sim.cpp; flows, reinfection, treatment, inflow
Rscript scripts/run_tests.R --scope all           # expect 6 files, 62 checks, 0 failures
```

`tests/helper.R` holds the shared assertion helpers. No test exercises `scenarios_2026.R`.
Its checks are the `stopifnot` calls inside `evaluate_draw` (capacity never exceeded, costs
reconcile, discount factors correct).

---

## Full reproduction in one block

```bash
cd HepC-Transimission-Modelling-feature-AddRecovery
PY=/opt/homebrew/Caskroom/miniconda/base/envs/bayes-venv/bin/python

Rscript scripts/run_tests.R --scope all                                    # 0. sanity
# (optional) scripts/model/preprocess.R interactively                      # 1. rate fits
Rscript scripts/run_calibration.R --run-id r_20260929_agepop \
  --maxit 1000 --n-starts 5 --seed 20261001                                # 2. best fit
Rscript scripts/run_npe.R --step all --fit output/calibration/r_20260929_agepop/fit.rds \
  --out-dir output/calibration/npe_agepop_recal --n-sims 3000 --n-draws 5000 --n-cores 4 \
  --n-rounds 1 --max-epochs 500 --seed 20261001 --npe-contact-sd 0.20 --npe-beta-sd 0.20 \
  --npe-hidden-features 128 --npe-num-transforms 10 --python $PY             # 3. posterior
Rscript scripts/run_mcmc_validation.R --fit output/calibration/r_20260929_agepop/fit.rds \
  --npe-dir output/calibration/npe_agepop_recal --out-dir output/calibration/mcmc_recal \
  --n-iter 120000 --burnin 20000 --thin 1 --n-chains 4 --parallel-chains 4 --seed 20261002 \
  --contact-sd 0.20 --beta-sd 0.20 --dt 0.0833 --ess-min 4000 --rhat-max 1.001 --same-start true  # 4.
Rscript scripts/predict_2015_posterior.R --fit output/calibration/r_20260929_agepop/fit.rds \
  --out-dir output/calibration/npe_agepop_recal --year 2015 --cores 4       # 5.
Rscript scripts/plot_calibration.R --fit output/calibration/r_20260929_agepop/fit.rds \
  --npe-dir output/calibration/npe_agepop_recal --mcmc-dir output/calibration/mcmc_recal \
  --cores 4 --max-draws 2000                                               # 6.
Rscript scripts/scenarios_2026.R --root . --draws-per-seed 1000 --cores 14 # 7. scenarios
Rscript scripts/paper_figures.R .                                          # 8. paper figures
python3 scripts/build_paper.py .                                           # 9. manuscript
Rscript scripts/run_tests.R --scope all                                    # 10. recheck
```

**To reproduce the published numbers exactly,** keep the existing `fit.rds` and NPE
posterior CSVs and start at step 7. Steps 2 to 4 regenerate them, and the NPE step in
particular will not be bit-identical.

---

## Known pitfalls

1. **Missing `data/` in new worktrees.** Restore the whitelisted files (Stage 3) before anything else.
2. **Parameter mismatches.** `pi_recid` is 0.59 in code and 0.65 in the PDF, and the arrest
   and release rates are marked "GUESS" in `setup.R`. The manuscript S1 Text reports the
   values in code.
3. **The 2019 prison treatment rate is an assumption.** It is
   `j_annual_start_probability = 0.99995` in the config and is not fitted to observed
   treatment counts.
4. **Outcomes are events per year, not rates per person.** Treatment starts count courses, not people.
5. **Unused files.** `scripts/model/Rplots.pdf` and `scripts/__pycache__/` are artefacts.
   `scripts/probe_beta_constraint.R` is a one-off diagnostic of a beta constraint, not part of the pipeline.
6. **Historical documents.** `回答.md` describes an earlier 75-scenario design and
   `结论.md` the earlier single-point results. The current results are in
   `output/scenarios_2026/README.md` and the manuscript.
7. **Background runs.** `mclapply` workers can outlive a closed shell. Kill stale runs
   before starting a new one.
