# Single-round NPE (50,000 LHS) versus reference MCMC

Generated 2026-10-11 by `scripts/compare_npe_mcmc.R`. Both methods target the same posterior
(`model_audit.md`).

## Runs compared

| Run | Prior simulations (valid) | Notes |
| --- | --- | --- |
| Reference MCMC (`../mcmc_reference/`) | 200,000 likelihood evaluations | Adaptive random-walk Metropolis, 4 × 45,000 retained draws; R-hat ≤ 1.002, bulk ESS ≥ 3,433, MCSE ≤ 0.017 SD. Production posterior |
| **NPE, LHS** (`../npe_lhs/`) | **50,000 (49,765)** | 235 equilibrium failures and 0 non-finite runs, dropped and not replaced. Design seed 20261012, training seed 20261013. Simulation took 1,389 s on 14 workers and training 37.5 s (71 epochs) |
| NPE, random, replicates 1–3 (`previous_random_npe/`) | 26,000 each (25,865–25,885) | Earlier runs with plain random prior draws and the same network and training settings. Shown for context |

Network and training are the same for all NPE runs:
- **Network:** MAF with 50 hidden features and 5 transforms, using single-round `NPE_C` with the
  standard NLL loss.
- **Training:** batch 256, Adam with lr 5e-4, validation fraction 0.10, patience 20, at most 250
  epochs, gradient clipping 5.0.
- **Output:** 20,000 posterior draws.

## Parameter accuracy

Tolerances are |SMD| ≤ 0.1, |95% endpoint| ≤ 0.2 SD and width ratio 0.9–1.1. The table uses
`accuracy_summary.csv`.

| Run | Worst \|SMD\| | Worst endpoint (SD) | 95% width ratio, mean (range) | Parameters passing all checks | Max \|corr difference\| |
| --- | --- | --- | --- | --- | --- |
| **NPE, LHS** | **0.091** | **0.26** | **0.961 (0.94–0.99)** | **11/12** | **0.052** |
| NPE, random, replicate 1 | 0.089 | 0.28 | 0.961 (0.92–1.05) | 10/12 | 0.047 |
| NPE, random, replicate 2 | 0.091 | 0.27 | 0.949 (0.90–1.00) | 10/12 | 0.048 |
| NPE, random, replicate 3 | 0.088 | 0.24 | 0.953 (0.91–0.99) | 11/12 | 0.059 |

**LHS run in detail.**
- **Means:** every parameter passes the mean check.
- **Endpoint failure:** the one failure is θ10, the entry scale for ages 40–49, whose lower 95%
  endpoint sits 0.26 SD above the reference.
- **Systematic narrowing:** the 95% intervals are 1–6% narrower than the reference for every
  parameter. The lower endpoints lie above the reference by +0.08 to +0.26 SD for most parameters.
  The narrowing pattern is the same as in the random-sampling runs, but the lowest width ratio rose
  from 0.90–0.92 to 0.94.

**What can be attributed to LHS.** The LHS run differs from the earlier runs in both sampling
method and size (50,000 versus 26,000). The improvement can't be attributed to LHS alone, and it is
within the spread of the three random replicates.

Marginals and pair plots (`figures/marginals.png`, `figures/pairs.png`) show one mode and no missing
mass.

## Policy outcomes

Seven predeclared outcomes, each from 1,000 posterior states, compared with the 1,000-state reference
run (`policy_comparison.csv`):

| Outcome | NPE, LHS: median (95% interval) | Reference | Error |
| --- | --- | --- | --- |
| J only, incidence reduction | 75.7% (70.2–80.5) | 75.4% (68.8–80.3) | 0.38% |
| J only, mortality reduction | 42.6% (41.6–43.7) | 42.6% (41.7–43.6) | 0.13% |
| 5,000/yr F 25%, mortality reduction | 64.9% (63.3–66.3) | 64.9% (63.4–66.5) | 0.03% |
| 5,000/yr F 25%, joint pass probability | 44.5% | 43.5% | 1.0 pp (MCSE ≈ 1.6 pp) |
| 8,000/yr F 25%, mortality reduction | 70.6% (69.3–71.7) | 70.6% (69.3–71.8) | 0.03% |
| 8,000/yr F 75%, mortality reduction | 57.2% (55.7–58.7) | 57.2% (55.7–58.8) | 0.03% |
| 8,000/yr F 25%, discounted cost | S$75.1m (69.7–81.7) | S$75.3m (69.2–81.6) | 0.21% |

- **Continuous outcomes:** all medians are within 0.4% of the reference, and the 95% endpoints are
  within 2.1%.
- **Joint-pass probability:** the difference is within its Monte Carlo error.
- **Equilibrium leakage:** all 1,000 NPE states passed the equilibrium check.
- **Random-sampling runs:** worst median errors were 0.54% (replicate 1), 0.22% (replicate 2) and
  0.16% (replicate 3).

## Conclusion

- **Policy outcomes:** single-round NPE on 50,000 LHS prior simulations reproduces the reference
  results well within 1%.
- **Parameter intervals:** about 4% too narrow on average, which is a mild, consistent
  overconfidence.
- **Production posterior:** the reference MCMC remains the production posterior because it samples
  the exact target.
- **Role of NPE:** a fast, independent check. It takes about 23 minutes of simulation on 14 cores
  plus under 1 minute of training, against about 3.8 hours for the reference MCMC on 4 cores.
