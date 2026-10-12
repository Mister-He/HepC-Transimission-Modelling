# =============================================================================
# prior_sampling.R — draws from the independent Normal prior on the 12 log-scale
# calibration parameters, by Latin hypercube sampling (LHS) or plain random sampling.
# =============================================================================

# LHS design on the unit hypercube: each column has exactly one point in each of
# the n equal-probability strata (random permutation plus uniform jitter).
lhs_unit <- function(n, d, seed) {
  set.seed(seed)
  vapply(seq_len(d), function(j) (sample.int(n) - runif(n)) / n, numeric(n))
}

# Map unit-cube points to the prior with the Normal inverse CDF, column by column.
# The prior is independent Normal(mean, sd) on theta (log scale), so no Jacobian or
# dependence structure is needed and there is no truncation.
prior_from_unit <- function(u, mean, sd) {
  stopifnot(ncol(u) == length(mean), length(sd) %in% c(1L, length(mean)))
  sd <- rep_len(sd, length(mean))
  vapply(seq_along(mean), function(j) qnorm(u[, j], mean[j], sd[j]), numeric(nrow(u)))
}

sample_prior_lhs <- function(n, mean, sd, seed) prior_from_unit(lhs_unit(n, length(mean), seed), mean, sd)

sample_prior_random <- function(n, mean, sd, seed) {
  set.seed(seed)
  sd <- rep_len(sd, length(mean))
  vapply(seq_along(mean), function(j) rnorm(n, mean[j], sd[j]), numeric(n))
}
