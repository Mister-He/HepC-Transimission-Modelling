# Single-round NPE (50,000 LHS) versus reference MCMC

Written 2026-10-22 (after the audit fixes) from the outputs of `scripts/compare_npe_mcmc.R` and `scripts/sbc_mcmc.R`.
**Cross-validation only:** the reference MCMC is the production posterior, and nothing here feeds
the production path.

Both methods target the same posterior (`model_audit.md`):
- **Parameters:** six log contact scales and log `c` (first arrest = `c` × re-arrest).
- **Prior and bounds:** the same fixed truncated Normal prior (contact N(0, 1.5), log `c`
  N(0, 0.5), `c` 0.5–100).
- **Likelihood:** seroprevalence only, Beta(p·φ, (1−p)·φ) with φ = 100; no counts.
- **Simulator:** CNB entry unscaled; mortality `mu` × SMR (9.52 viraemic, 4.73 otherwise).

## Runs compared

| Run | Cost | Notes |
| --- | --- | --- |
| Reference MCMC (`../mcmc_reference/`) | 4 × 50,000 likelihood evaluations, 4.0 h on 4 cores | Production posterior. Warmup 10,000, adaptation from 5,000 on the recent half of warmup, thin 10 (16,000 stored draws). R-hat ≤ 1.002, bulk ESS ≥ 3,250, tail ESS ≥ 3,487, MCSE ≤ 0.019 SD |
| **NPE, LHS** (`../npe_lhs/`) | **50,000 prior simulations (49,945 valid)**, 1,733 s of simulation on 10 workers, 345 s of training | 43 equilibrium failures, 12 non-finite, 0 out of bounds. 6 summaries (logit seroprevalence with Beta noise). Design seed 20261012, training seed 20261013. Equal-weight ensemble of 5 MAFs (50×5), single-round `NPE_C` with the standard NLL loss; draws outside the bounds rejected |

**Training settings:** batch 256, Adam with lr 5e-4, validation fraction 0.10, patience 20, at most
250 epochs and gradient clipping 5.0. The members trained for 90–117 epochs. The run drew 20,000
posterior draws.

## Simulation-based calibration

| Check | Datasets | 95% coverage | 90% coverage | Rank shape |
| --- | --- | --- | --- | --- |
| NPE (`../npe_lhs/sbc_summary.csv`) | 998 random prior draws (2 failed) | 0.933–0.960 per parameter | 0.872–0.911 per parameter | all 7 uniform |
| Reference MCMC (`mcmc_sbc/`) | 8 (first 8 NPE SBC truths) | 0.96 pooled over 56 ranks | 0.89 pooled | mean rank 0.48 |

- **Gate:** cov95 in [0.85, 1], cov90 in [0.75, 1] and no systematically skewed rank histogram.
  The NPE passes (`NPE_OK.txt`).
- **MCMC SBC:** each of the 8 datasets got one chain of 30,000 evaluations (max split R-hat
  1.008). Ranks within a dataset are correlated, so these runs check for gross failure only.

## Parameter accuracy

Tolerances: |SMD| ≤ 0.1, |95% endpoint| ≤ 0.2 SD and width ratio 0.9–1.1 (`accuracy_summary.csv`).
The NPE draws are first restricted to the MCMC support: of 4,000 random NPE draws, 1,329 (33%) fail
the 2015 equilibrium check and are dropped (`npe_support.txt`).

| Parameter | SMD | 95% width ratio | Passes all |
| --- | --- | --- | --- |
| θ1 (<20 contact) | −0.56 | 1.29 | no |
| θ2 (20–29) | 0.08 | 1.20 | no |
| θ3 (30–39) | −0.04 | 1.37 | no |
| θ4 (40–49) | −0.14 | 1.89 | no |
| θ5 (50–59) | −0.26 | 1.04 | no |
| θ6 (60+) | −0.13 | 1.00 | no |
| θ7 (log `c`) | 0.35 | 1.15 | no |

- **Summary:** worst |SMD| 0.56, worst endpoint 2.12 SD, mean width ratio 1.28 (1.00–1.89),
  0 of 7 parameters pass all checks; max |correlation difference| 0.30.
- **Reading:** the NPE is calibrated on average over prior datasets (SBC) but not accurate at the
  observed data: a third of its mass lies outside the model's equilibrium support, and the rest is
  wider than the MCMC posterior. With six noisy summaries and a wide prior, almost no training
  simulation lies near the observed data, so the flow extrapolates there.
- **Prior-dominated parameters (MCMC):** θ5, θ6 and θ7 have posterior SD above 0.7 × the prior SD.

## Policy outcomes

Seven predeclared outcomes compare 1,000 NPE posterior states with the production scenario run
(`output/scenarios_2026/`, 2,000 reference MCMC states). Source: `policy_comparison.csv`.

| Outcome | NPE: median (95% interval) | Reference | Error |
| --- | --- | --- | --- |
| J only, incidence reduction | 85.9% (40.8–93.7) | 85.7% (70.5–91.3) | 0.28% |
| J only, mortality reduction | 54.7% (50.0–61.0) | 54.9% (52.3–57.4) | 0.28% |
| 5,000/yr F 25%, mortality reduction | 68.4% (65.2–72.3) | 68.5% (66.6–70.4) | 0.23% |
| 5,000/yr F 25%, joint pass probability | 88.8% | 97.0% | 8.2 pp |
| 8,000/yr F 25%, mortality reduction | 74.0% (71.3–77.3) | 74.1% (72.4–75.7) | 0.19% |
| 8,000/yr F 75%, mortality reduction | 63.7% (59.6–68.7) | 63.8% (61.4–66.2) | 0.15% |
| 8,000/yr F 25%, discounted cost | S$96.5m (18.4–163.1) | S$102.6m (86.1–122.6) | 6.0% |

- **Reductions:** medians within 0.3% of the reference, with wider intervals; the wider intervals
  lower the joint-pass probability by 8 percentage points.
- **Cost:** 6% lower with a much wider interval, because cost depends on absolute population size,
  which the NPE posterior spreads more widely.

## Conclusion

- **Agreement:** poor for parameters (0 of 7 pass, a third of draws outside the support),
  adequate for policy reduction medians (≤ 0.3%).
- **Production posterior:** the reference MCMC samples the exact target, passes its own SBC
  check and fits all six targets, so it remains the production posterior.
- **Role of NPE:** independent cross-validation, not validated under this prior. It costs about 29
  minutes of simulation on 10 workers plus 6 minutes of training, against 4.0 hours for the MCMC on
  4 cores.
