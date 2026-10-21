# Single-round NPE (50,000 LHS) versus reference MCMC

Written 2026-10-22 from the outputs of `scripts/compare_npe_mcmc.R` and `scripts/sbc_mcmc.R`.
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
| Reference MCMC (`../mcmc_reference/`) | 4 × 50,000 likelihood evaluations, 4.0 h on 4 cores | Production posterior. Warmup 10,000, adaptation from 5,000, thin 10 (16,000 stored draws). R-hat ≤ 1.005, bulk ESS ≥ 1,798, tail ESS ≥ 2,459, MCSE ≤ 0.024 SD |
| **NPE, LHS** (`../npe_lhs/`) | **50,000 prior simulations (49,945 valid)**, 1,784 s of simulation on 10 workers, 380 s of training | 42 equilibrium failures, 13 non-finite, 0 out of bounds. 6 summaries (logit seroprevalence with Beta noise). Design seed 20261012, training seed 20261013. Equal-weight ensemble of 5 MAFs (50×5), single-round `NPE_C` with the standard NLL loss; draws outside the bounds rejected |

**Training settings:** batch 256, Adam with lr 5e-4, validation fraction 0.10, patience 20, at most
250 epochs and gradient clipping 5.0. The members trained for 95–138 epochs. The run drew 20,000
posterior draws.

## Simulation-based calibration

| Check | Datasets | 95% coverage | 90% coverage | Rank shape |
| --- | --- | --- | --- | --- |
| NPE (`../npe_lhs/sbc_summary.csv`) | 998 random prior draws (2 failed) | 0.935–0.959 per parameter | 0.876–0.911 per parameter | all 7 uniform |
| Reference MCMC (`mcmc_sbc/`) | 8 (first 8 NPE SBC truths) | 0.96 pooled over 56 ranks | 0.89 pooled | mean rank 0.47 |

- **Gate:** cov95 in [0.85, 1], cov90 in [0.75, 1] and no systematically skewed rank histogram.
  The NPE passes (`NPE_OK.txt`).
- **MCMC SBC:** each of the 8 datasets got one chain of 30,000 evaluations (max split R-hat
  1.008). Ranks within a dataset are correlated, so these runs check for gross failure only.

## Parameter accuracy

Tolerances: |SMD| ≤ 0.1, |95% endpoint| ≤ 0.2 SD and width ratio 0.9–1.1 (`accuracy_summary.csv`).

| Parameter | SMD | 95% width ratio | Passes all |
| --- | --- | --- | --- |
| θ1 (<20 contact) | −0.83 | 1.19 | no |
| θ2 (20–29) | −0.10 | 1.13 | no |
| θ3 (30–39) | −0.38 | 1.33 | no |
| θ4 (40–49) | −0.63 | 1.79 | no |
| θ5 (50–59) | −0.07 | 1.06 | no |
| θ6 (60+) | −0.20 | 1.03 | no |
| θ7 (log `c`) | −0.10 | 1.08 | no |

- **Summary:** worst |SMD| 0.83, worst endpoint 2.10 SD (θ4), mean width ratio 1.23 (1.03–1.79),
  0 of 7 parameters pass all checks; max |correlation difference| 0.33.
- **Reading:** the NPE is calibrated on average over prior datasets (SBC) but not accurate at the
  observed data: its posterior is wider and shifted for the young and middle-age contact scales,
  and its predictive seroprevalence is biased low (`figures/calibration_agepop.png`). With six
  noisy summaries and a wide prior, the flow spends its capacity across the whole prior range.
- **Prior-dominated parameters (MCMC):** θ5, θ6 and θ7 have posterior SD above 0.7 × the prior SD.

## Policy outcomes

Seven predeclared outcomes compare 1,000 NPE posterior states with the production scenario run
(`output/scenarios_2026/`, 2,000 reference MCMC states). Source: `policy_comparison.csv`.

| Outcome | NPE: median (95% interval) | Reference | Error |
| --- | --- | --- | --- |
| J only, incidence reduction | 87.6% (44.1–95.0) | 85.6% (68.9–91.3) | 2.4% |
| J only, mortality reduction | 46.1% (43.6–48.9) | 45.8% (44.2–47.4) | 0.68% |
| 5,000/yr F 25%, mortality reduction | 56.2% (54.5–58.0) | 56.1% (54.9–57.2) | 0.15% |
| 5,000/yr F 25%, joint pass probability | 0% | 0% | 0 pp |
| 8,000/yr F 25%, mortality reduction | 60.5% (59.2–62.0) | 60.5% (59.5–61.3) | 0.01% |
| 8,000/yr F 75%, mortality reduction | 52.4% (50.3–54.8) | 52.2% (50.8–53.6) | 0.38% |
| 8,000/yr F 25%, discounted cost | S$85.0m (23.5–151.9) | S$101.5m (85.9–122.7) | 16% |

- **Reductions:** medians within 2.4% of the reference, with wider intervals.
- **Cost:** 16% lower with a much wider interval, because cost depends on absolute population
  size, which the NPE posterior spreads more widely.

## Conclusion

- **Agreement:** poor for parameters (0 of 7 pass), adequate for policy reductions (≤ 2.4%).
- **Production posterior:** the reference MCMC samples the exact target, passes its own SBC
  check and fits all six targets, so it remains the production posterior.
- **Role of NPE:** independent cross-validation. It costs about 30 minutes of simulation on 10
  workers plus 6 minutes of training, against 4.0 hours for the MCMC on 4 cores.
