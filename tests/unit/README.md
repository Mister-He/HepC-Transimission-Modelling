# Unit tests guide

Unit tests cover metrics, equilibrium checks, likelihood bounds and the prior
centre, target loading, the rate-fit handoff into `setup.R`, and the adaptive Metropolis sampler on known
targets (`test_mcmc.R`), and the decision rules of `scripts/decision_analysis.R` (`test_decision.R`:
Pareto efficiency, normalisation, TOPSIS, efficiency frontier, budget rule), without requiring the
full pipeline.
