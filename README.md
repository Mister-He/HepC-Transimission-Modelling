# HCV transmission model for people who inject drugs (PWID) in Singapore

An age-structured compartmental model of hepatitis C virus (HCV) among
people who inject drugs (PWID) in Singapore, calibrated to the age-specific
HCV seroprevalence and population of the prison (detained) PWID stratum.

## 1. Project structure

```text
.
├── src/
│   ├── setup.R                # Parameters, initial conditions, helpers
│   ├── sim.cpp                # Compiled ODE solver (3-strata model)
│   └── calibration/
│       ├── targets.R          # Calibration target data + Binomial counts
│       ├── model_metrics.R    # J summary extraction + fit metrics
│       ├── equilibrium.R      # T vs T-5 stability gate
│       ├── likelihood.R       # 12-parameterisation + NLL objective
│       ├── calibrate_nm.R     # Multi-start Nelder-Mead runner + informed/warm starts
│       ├── laplace.R          # Laplace approximation (Hessian, MC intervals)
│       ├── run_calibration.R  # End-to-end calibration runner
│       ├── plot_results.R     # ggplot2 publication figures (all figures)
│       ├── probe_beta_constraint.R  # beta_scale > 1 feasibility probe
│       └── (analysis probes)  # Variant/sensitivity probes
├── docs/calibration/
│   ├── preflight.md           # Pre-run record (Git SHA, hashes)
│   ├── model_audit.md         # Model audit + failure diagnosis
│   ├── mortality_review.md    # 2015 mortality + PWID SMR evidence
│   ├── natural_history_review.md  # Liver progression rate evidence
│   ├── likelihood.md          # Statistical model specification
│   ├── DECISIONS.md           # Append-only decision log
│   ├── final_report.md        # Final calibration report
│   └── PROJECT_OVERVIEW.md    # Full experiment walk-through (what was tried, file provenance)
├── output/calibration/        # Per-run outputs (run_config, solutions, predictions, ...)
├── figures/                   # Final ggplot2 figures
├── README.md / README.zh-CN.md
├── AGENTS.md / prompt.md      # Working instructions for agents
└── Model schematic.pptx       # Visual record (protected, never modified)
```

## 2. Model structure

```text
3 strata:  D (community PWID, inject), J (detained), X (former PWID, non-injecting)
4 liver stages: NC, CC, DC, HCC
5 HCV states:  u (susceptible), a (acute), c (chronic), t (treatment), s (seropositive cleared/post-SVR)
6 age groups:  <20, 20-29, 30-39, 40-49, 50-59, 60+
360 compartments
```

- **Transmission** occurs only in D (detained PWID do not inject; X do not
  re-inject): frequency-dependent force of infection
  `gamma_i = q * sum_j C_contact(i,j) * infectious_D_j / active_D_j`,
  with a **constant** transmission multiplier (m = 1; the historical level
  is absorbed by the fitted contact row scales).
- **Ageing**: 10-year bands, `y_change = y/10 per year` to the next band
  (60+ open-ended).
- **Arrest/release**: D -> J at `lambda1[i]`; J -> D with probability
  `pi_recid` or J -> X with `1 - pi_recid` at `lambda2[i]`.
- **Mortality**: `mu[i] * omega` (omega = 14.68 PWID SMR);
  decompensated cirrhosis adds `mu_DC = 0.130/yr`, HCC adds
  `mu_HCC = 0.430/yr`.
- **Progression rates** (highest defensible literature values):
  `p_NC_CC = 0.027` (Thein 2008); `p_CC_DC = 0.0788` and
  `p_CC_HCC = 0.0479` (Alazawi et al. 2010, PMID 20497143, untreated-only
  rates); `p_DC_HCC = 0.0464` (Rivera-Irigoin et al. 2006, PMID 17209765,
  upper 95% CI). GT3 relative risks (Kanwal et al. 2014) retained.

## 3. Data sources

- **Targets**: age-specific anti-HCV seroprevalence and prison population
  from the Changi Prison universal screening 2014-2016 (source paper:
  Park et al., medRxiv 10.1101/2025.10.24.25338708).
- **Background mortality**: Singapore SingStat age-specific death rates,
  2015, mapped to the six age groups (see `mortality_review.md`).
- **PWID excess mortality**: Mathers et al. 2013 pooled SMR 14.68
  (PMID 23554523).
- **Progression rates**: Alazawi et al. 2010 (PMID 20497143);
  Rivera-Irigoin et al. 2006 (PMID 17209765); Thein et al. 2008.

## 4. Calibration

- **12 fitted parameters**: 6 contact-matrix row scales, 6 beta inflow
  scales (all `exp(theta)`). No transmission multiplier, no
  excess-mortality parameter (the 60+ group is fitted within the
  acceptance criteria with these 12).
- Objective: Binomial NLL for prevalence + log-Normal NLL for population
  (sigma_pop = 0.10); equilibrium gate (T vs T-5) as a penalty.
- Optimiser: deterministic multi-start Nelder-Mead
  (`optim`, method = "Nelder-Mead"). No HMC.
- Uncertainty: Laplace approximation at the optimum (finite-difference
  Hessian, generalized inverse with eigenvalue cutoff 1e-4, 1000 Monte
  Carlo draws filtered by equilibrium feasibility), 95% intervals compared
  with observed intervals per age group.
- Acceptance criteria: prevalence RMSE <= 0.02, max |prev err| <= 0.03,
  population MAPE <= 0.10, max APE <= 0.20, equilibrium pass.
- `beta_scale > 1` is a soft target: satisfied for age groups 3-6; the
  <20 and 20-29 groups remain below 1 because their small prison
  populations require net inflows below the CNB-anchored base `beta`
  (see DECISIONS.md).

### 4.1 Bayesian inference (MCMC)

On top of the NM fit, Bayesian posterior sampling is performed. Per
`prompt_mcmc.md` (revised), **NPE/SNPE was attempted first** (sbi `NPE_C`,
3 sequential rounds, two seeds, SBC-calibrated) but was not adopted: it
was overconfident in well-identified directions and unstable across seeds
in weakly identified directions. The **traditional MCMC (adaptive
Metropolis-Hastings) is the primary Bayesian summary**, with strict
convergence: R-hat in [0.99, 1.01] (achieved 1.0008-1.0096), ESS > 400
(achieved 1,374-2,457), 3 chains x 40,000 iterations. Priors: log contact
scales ~ Normal(0, 2^2); log beta scales ~ Student-t(3, 0, 2). Posterior
predictive 95% credible intervals per age group overlap both the Laplace
and observed intervals. See `prompt_mcmc.md` and
`docs/calibration/bayes_methodology.md`.

## 5. How to regenerate

```bash
# 1. Full calibration run (defaults reproduce the final run)
Rscript src/calibration/run_calibration.R \
  --run-id my_run --seed 101 --maxit 3000 --n-starts 6 \
  --t-start -10 --t-end 55 --target-time 45 --target-mode sero \
  --out-dir output/calibration

# 2. Publication figures (ggplot2)
Rscript src/calibration/plot_results.R \
  --fit output/calibration/<run_id>/fit.rds --out-dir figures

# 3. Bayesian posterior (NPE 3 rounds + MCMC strict validation)
Rscript src/calibration/run_npe.R --step all \
  --root . --fit output/calibration/<run_id>/fit.rds \
  --out-dir output/calibration/npe_bayes \
  --n-sims 60000 --n-cores 6 --seed 2026 \
  --n-draws 60000 --n-proposal 20000 --n-rounds 3 \
  --n-iter-mcmc 20000 --burnin-mcmc 5000 --mcmc-max-iter 400000
```

Per-run outputs: `run_config.csv`, `targets.csv`, `initial_values.csv`,
`optimization_history.csv`, `solutions.csv`, `predictions.csv`,
`residuals.csv`, `equilibrium.csv`, `laplace_intervals.csv`,
`laplace_diagnostics.csv`, `sessionInfo.txt`, `fit.rds`, and ggplot2
figures.

## 6. Results

Final run: `output/calibration/run5_redo150/`.

Acceptance criteria are met with 12 parameters and no excess mortality.
Bayesian outputs are in `output/calibration/npe_bayes/`; 50-year
sensitivity projections are in `output/analysis/`. See
`docs/calibration/final_report.md` and `docs/calibration/analysis_report.md`.
