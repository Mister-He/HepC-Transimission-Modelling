# Single-round NPE (50,000 LHS) versus reference MCMC

Written 2026-10-21 from the outputs of `scripts/compare_npe_mcmc.R` and `scripts/sbc_mcmc.R`.
**Cross-validation only:** the reference MCMC is the production posterior, and nothing here feeds
the production path.

Both methods target the same posterior (`model_audit.md`):
- **Prior and bounds:** the same prior centre and bounds.
- **Likelihood:** prevalence and population (no entry term).
- **Simulator:** `lambda1 = lambda3` and the viraemic excess mortality.

## Runs compared

| Run | Cost | Notes |
| --- | --- | --- |
| Reference MCMC (`../mcmc_reference/`) | 4 × 50,000 likelihood evaluations, 16,787 s on 4 cores | Production posterior. Warmup 10,000, adaptation from 5,000, thin 10 (16,000 stored draws). R-hat ≤ 1.002, bulk ESS ≥ 3,053, tail ESS ≥ 5,335, MCSE ≤ 0.018 SD |
| **NPE, LHS** (`../npe_lhs/`) | **50,000 prior simulations (49,278 valid)**, 1,616 s of simulation on 12 workers, 236 s of training | 722 equilibrium failures, 0 out of bounds, 0 non-finite. 12 summaries (prevalence, population). Design seed 20261012, training seed 20261013. Equal-weight ensemble of 5 MAFs (50×5), single-round `NPE_C` with the standard NLL loss |

**Training settings:** batch 256, Adam with lr 5e-4, validation fraction 0.10, patience 20, at most
250 epochs and gradient clipping 5.0. The members trained for 54–93 epochs. The run drew 20,000
posterior draws.

## Simulation-based calibration

| Check | Datasets | 95% coverage | 90% coverage | Rank shape |
| --- | --- | --- | --- | --- |
| NPE (`../npe_lhs/sbc_summary.csv`) | 980 random prior draws (20 failed equilibrium) | 0.944–0.961 per parameter | 0.891–0.918 per parameter | all 12 uniform |
| Reference MCMC (`mcmc_sbc/`) | 8 (first 8 NPE SBC truths) | 0.95 pooled over 96 ranks | 0.92 pooled | mean rank 0.498 |

- **Gate:** cov95 in [0.85, 1], cov90 in [0.75, 1] and no systematically skewed rank histogram.
  The NPE passes (`NPE_OK.txt`).
- **MCMC SBC:**
  - Each of the 8 datasets got one chain of 30,000 evaluations (max split R-hat 1.016).
  - Ranks within a dataset are correlated, so these runs check for gross failure only.

## Parameter accuracy

Tolerances: |SMD| ≤ 0.1, |95% endpoint| ≤ 0.2 SD and width ratio 0.9–1.1 (`accuracy_summary.csv`).

| Worst \|SMD\| | Worst endpoint (SD) | 95% width ratio, mean (range) | Parameters passing all checks | Max \|corr difference\| |
| --- | --- | --- | --- | --- |
| 0.030 | 0.218 | 0.965 (0.92–0.99) | 11/12 | 0.034 |

- **Means:** every mean is within 0.03 posterior SD of the reference.
- **Intervals:** NPE 95% intervals are on average 3.5% narrower. θ4 (the 40–49 contact scale)
  fails only the endpoint check (0.22 SD, against a tolerance of 0.2); its width ratio is 0.92.
- **Prior-dominated parameters:** four parameters have posterior SD above 0.8 × the prior SD: θ6
  (the 60+ contact scale) and θ9, θ10, θ12 (the 30–39, 40–49 and 60+ entry scales).

## Policy outcomes

Seven predeclared outcomes compare 1,000 NPE posterior states with the production scenario run
(`output/scenarios_2026/`, 2,000 reference MCMC states). Source: `policy_comparison.csv`.

| Outcome | NPE: median (95% interval) | Reference | Error |
| --- | --- | --- | --- |
| J only, incidence reduction | 81.1% (76.8–84.3) | 81.1% (76.5–84.3) | 0.01% |
| J only, mortality reduction | 50.5% (49.8–51.3) | 50.5% (49.8–51.3) | 0.00% |
| 5,000/yr F 25%, mortality reduction | 70.4% (69.5–71.3) | 70.4% (69.4–71.3) | 0.06% |
| 5,000/yr F 25%, joint pass probability | 100% | 100% | 0 pp |
| 8,000/yr F 25%, mortality reduction | 73.3% (72.8–73.8) | 73.3% (72.8–73.8) | 0.04% |
| 8,000/yr F 75%, mortality reduction | 64.8% (63.8–65.9) | 64.9% (63.7–66.0) | 0.10% |
| 8,000/yr F 25%, discounted cost | S$55.9m (51.0–60.5) | S$55.8m (50.9–61.0) | 0.23% |

- **Continuous outcomes:** all medians are within 0.3% of the reference.

## Conclusion

- **Agreement:** the NPE agrees with the reference MCMC on 11 of 12 parameter checks; its 95%
  intervals are slightly narrow (mean ratio 0.965). Policy medians are within 0.3%.
- **Production posterior:** the reference MCMC samples the exact target, so it remains the
  production posterior.
- **Role of NPE:** independent cross-validation. It costs about 27 minutes of simulation on 12
  workers plus 4 minutes of training, against 4.7 hours for the MCMC on 4 cores.
