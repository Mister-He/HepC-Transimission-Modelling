# Single-round NPE (50,000 LHS) versus reference MCMC

Written 2026-10-14 from the outputs of `scripts/compare_npe_mcmc.R` and `scripts/sbc_mcmc.R`.
**Cross-validation only:** the reference MCMC is the production posterior, and nothing here feeds
the production path. Both methods target the same posterior, with the same prior centre and bounds
(`model_audit.md`).

## Runs compared

| Run | Cost | Notes |
| --- | --- | --- |
| Reference MCMC (`../mcmc_reference/`) | 4 × 50,000 likelihood evaluations, 15,548 s on 4 cores | Production posterior. Warmup 10,000, adaptation from 5,000, thin 10 (16,000 stored draws). R-hat ≤ 1.002, bulk ESS ≥ 2,824, tail ESS ≥ 5,119, MCSE ≤ 0.019 SD |
| **NPE, LHS** (`../npe_lhs/`) | **50,000 prior simulations (all valid)**, 1,618 s of simulation, 243 s of training | 0 out-of-bounds, 0 equilibrium failures and 0 non-finite runs. Design seed 20261012, training seed 20261013. Equal-weight ensemble of 5 MAFs (50×5), single-round `NPE_C` with the standard NLL loss |

**Training settings:** batch 256, Adam with lr 5e-4, validation fraction 0.10, patience 20, at most
250 epochs and gradient clipping 5.0. The members trained for 54–78 epochs. The run drew 20,000
posterior draws.

## Simulation-based calibration

| Check | Datasets | 95% coverage | 90% coverage | Rank shape |
| --- | --- | --- | --- | --- |
| NPE (`../npe_lhs/sbc_summary.csv`) | 1,000 random prior draws | 0.931–0.958 per parameter | 0.873–0.910 per parameter | all 12 uniform |
| Reference MCMC (`mcmc_sbc/`) | 8 (first 8 NPE SBC truths) | 1.00 pooled over 96 ranks | 0.95 pooled | mean rank 0.503 |

- **Gate:** cov95 in [0.85, 1], cov90 in [0.75, 1] and no systematically skewed rank histogram.
  The NPE passes (`NPE_OK.txt`).
- **Rank shape rule:** the histogram fails if the outer-10% share or the mean rank departs by more
  than 3 SE, or if the 20-bin chi-square p is below 0.001/12.
- **Single network:** one network failed the rule on θ1 (mean rank 0.530 on 1,000 datasets, 0.519 on
  an independent 3,000). Other seeds moved the offset to other parameters. The 5-network ensemble
  removed it (max |z| 2.1 on the independent 3,000).
- **MCMC SBC:** 8 datasets check only for gross failure; each chain had 30,000 evaluations (max split
  R-hat 1.021).

## Parameter accuracy

Tolerances: |SMD| ≤ 0.1, |95% endpoint| ≤ 0.2 SD and width ratio 0.9–1.1 (`accuracy_summary.csv`).

| Worst \|SMD\| | Worst endpoint (SD) | 95% width ratio, mean (range) | Parameters passing all checks | Max \|corr difference\| |
| --- | --- | --- | --- | --- |
| 0.047 | 0.200 | 0.963 (0.93–1.01) | 12/12 | 0.034 |

- **Means:** every mean is within 0.05 posterior SD of the reference.
- **Intervals:** 95% intervals are about 4% narrower than the reference on average. The narrowest
  are θ8 and θ10, at 0.93.
  - Only θ6 (prior-dominated) is wider, at 1.01.
  - This mild narrowing matches the NPE SBC coverage, whose cov95 mean is 0.949.
- **Shape:** the marginals and pair plots (`figures/marginals.png`, `figures/pairs.png`) show one
  mode and no missing mass.
- **Prior-dominated parameters:** θ6 and θ12 (the 60+ contact and entry scales) have posterior
  SD ≈ prior SD in both methods.

## Policy outcomes

Seven predeclared outcomes compare 1,000 NPE posterior states with the production scenario run
(`output/scenarios_2026/`, 2,000 reference MCMC states). Source: `policy_comparison.csv`.

| Outcome | NPE: median (95% interval) | Reference | Error |
| --- | --- | --- | --- |
| J only, incidence reduction | 59.5% (46.4–69.1) | 59.2% (46.4–69.0) | 0.57% |
| J only, mortality reduction | 41.7% (40.9–42.7) | 41.8% (40.9–42.8) | 0.08% |
| 5,000/yr F 25%, mortality reduction | 63.9% (62.5–65.4) | 64.0% (62.5–65.6) | 0.17% |
| 5,000/yr F 25%, joint pass probability | 8.7% | 11.25% | 2.55 pp (combined MCSE ≈ 1.1 pp) |
| 8,000/yr F 25%, mortality reduction | 69.5% (68.3–70.6) | 69.6% (68.3–70.8) | 0.15% |
| 8,000/yr F 75%, mortality reduction | 56.4% (55.0–57.8) | 56.5% (55.0–58.1) | 0.20% |
| 8,000/yr F 25%, discounted cost | S$78.2m (72.2–84.9) | S$77.9m (71.3–84.7) | 0.35% |

- **Continuous outcomes:** all medians are within 0.6% of the reference, and the 95% endpoints are
  within 1.4%.
- **Joint-pass probability:** the gap is about 2.2 combined MCSE. This is a tail probability near
  the 65% threshold, and the slightly narrower NPE posterior puts fewer draws above the threshold.
- **Equilibrium leakage:** all 1,000 NPE states passed the equilibrium check.

## Conclusion

- **Agreement:** the NPE agrees with the reference MCMC on every parameter check and on the
  continuous policy outcomes within 1%.
- **Spread:** its intervals are about 4% narrower, which is consistent with its SBC coverage.
- **Joint-pass probability:** this tail probability near a threshold is the most sensitive quantity.
  NPE gives 8.7% against 11.25%.
- **Production posterior:** the reference MCMC is the production posterior because it samples the
  exact target.
- **Role of NPE:** independent cross-validation. It costs about 27 minutes of simulation on 11 workers
  plus 4 minutes of training, against 4.3 hours for the MCMC on 4 cores.
