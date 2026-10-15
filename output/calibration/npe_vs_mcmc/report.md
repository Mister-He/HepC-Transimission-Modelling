# Single-round NPE (50,000 LHS) versus reference MCMC

Written 2026-10-15 from the outputs of `scripts/compare_npe_mcmc.R` and `scripts/sbc_mcmc.R`.
**Cross-validation only:** the reference MCMC is the production posterior, and nothing here feeds
the production path.

Both methods target the same posterior (`model_audit.md`):
- **Prior and bounds:** the same prior centre and bounds.
- **Likelihood:** prevalence, population and the CNB entry anchor.
- **Simulator:** the first-arrest flow handoff and the viraemic excess mortality.

## Runs compared

| Run | Cost | Notes |
| --- | --- | --- |
| Reference MCMC (`../mcmc_reference/`) | 4 × 50,000 likelihood evaluations, 16,231 s on 4 cores | Production posterior. Warmup 10,000, adaptation from 5,000, thin 10 (16,000 stored draws). R-hat ≤ 1.003, bulk ESS ≥ 2,319, tail ESS ≥ 2,728, MCSE ≤ 0.021 SD |
| **NPE, LHS** (`../npe_lhs/`) | **50,000 prior simulations (48,918 valid)**, 1,612 s of simulation on 12 workers, 299 s of training | 1,082 equilibrium failures, 0 out of bounds, 0 non-finite. 18 summaries (prevalence, population, 2015 entry). Design seed 20261012, training seed 20261013. Equal-weight ensemble of 5 MAFs (50×5), single-round `NPE_C` with the standard NLL loss |

**Training settings:** batch 256, Adam with lr 5e-4, validation fraction 0.10, patience 20, at most
250 epochs and gradient clipping 5.0. The members trained for 89–102 epochs. The run drew 20,000
posterior draws.

## Simulation-based calibration

| Check | Datasets | 95% coverage | 90% coverage | Rank shape |
| --- | --- | --- | --- | --- |
| NPE (`../npe_lhs/sbc_summary.csv`) | 973 random prior draws (27 failed equilibrium) | 0.936–0.958 per parameter | 0.879–0.918 per parameter | all 12 uniform |
| Reference MCMC (`mcmc_sbc/`) | 8 (first 8 NPE SBC truths) | 0.93 pooled over 96 ranks | 0.83 pooled | mean rank 0.512 |

- **Gate:** cov95 in [0.85, 1], cov90 in [0.75, 1] and no systematically skewed rank histogram.
  The NPE passes (`NPE_OK.txt`).
- **MCMC SBC:**
  - Each of the 8 datasets got one chain of 30,000 evaluations (max split R-hat 1.016).
  - Ranks within a dataset are correlated, so the pooled cov90 of 0.83 is within noise for 8
    datasets.
  - These runs check for gross failure only.

## Parameter accuracy

Tolerances: |SMD| ≤ 0.1, |95% endpoint| ≤ 0.2 SD and width ratio 0.9–1.1 (`accuracy_summary.csv`).

| Worst \|SMD\| | Worst endpoint (SD) | 95% width ratio, mean (range) | Parameters passing all checks | Max \|corr difference\| |
| --- | --- | --- | --- | --- |
| 0.075 | 0.121 | 0.992 (0.96–1.01) | 12/12 | 0.042 |

- **Means:** every mean is within 0.08 posterior SD of the reference.
- **Intervals:** 95% intervals match to within 4%.
- **Prior-dominated parameters:** eight parameters have posterior SD above 0.8 × the prior SD: θ5,
  θ6 (the 50–59 and 60+ contact scales) and θ7–θ12 (the excess-entry scales).
  - The prison data inform the contact scales of ages 20–49 most. The CNB anchor (SD 1.0) informs
    entry only weakly.

## Policy outcomes

Seven predeclared outcomes compare 1,000 NPE posterior states with the production scenario run
(`output/scenarios_2026/`, 2,000 reference MCMC states). Source: `policy_comparison.csv`.

| Outcome | NPE: median (95% interval) | Reference | Error |
| --- | --- | --- | --- |
| J only, incidence reduction | 47.1% (36.3–56.9) | 47.6% (36.5–57.5) | 1.17% |
| J only, mortality reduction | 47.2% (44.6–48.8) | 47.4% (44.7–48.8) | 0.28% |
| 5,000/yr F 25%, mortality reduction | 64.7% (60.9–67.0) | 64.9% (60.9–67.0) | 0.24% |
| 5,000/yr F 25%, joint pass probability | 0% | 0% | 0 pp |
| 8,000/yr F 25%, mortality reduction | 67.0% (63.1–69.4) | 67.2% (63.1–69.4) | 0.25% |
| 8,000/yr F 75%, mortality reduction | 60.1% (56.6–62.2) | 60.3% (56.7–62.2) | 0.23% |
| 8,000/yr F 25%, discounted cost | S$69.0m (64.3–74.0) | S$68.8m (63.7–73.8) | 0.19% |

- **Continuous outcomes:** all medians are within 1.2% of the reference.
  - The largest error is J-only incidence (0.56 percentage points). That is about 2 combined Monte
    Carlo SE of the two medians (0.0022 and 0.0017).
- **Joint-pass probability:** this is no longer informative, because no draw meets the incidence
  target in either posterior.

## Conclusion

- **Agreement:** the NPE agrees with the reference MCMC on every parameter check, with 95%
  intervals within 4%. Policy medians are within 1.2%.
- **Production posterior:** the reference MCMC samples the exact target, so it remains the
  production posterior.
- **Role of NPE:** independent cross-validation. It costs about 27 minutes of simulation on 12
  workers plus 5 minutes of training, against 4.5 hours for the MCMC on 4 cores.
