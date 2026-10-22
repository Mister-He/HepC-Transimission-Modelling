# =============================================================================
# prior_sampling.R — draws from the independent (truncated) Normal prior on the log-scale
# calibration parameters, by Latin hypercube sampling (LHS) or plain random sampling.
# =============================================================================

# LHS design on the unit hypercube: each column has exactly one point in each of
# the n equal-probability strata (random permutation plus uniform jitter).
lhs_unit <- function(n, d, seed) {
  set.seed(seed)
  vapply(seq_len(d), function(j) (sample.int(n) - runif(n)) / n, numeric(n))
}

# Map unit-cube points to the prior with the Normal inverse CDF, column by column. The prior is
# independent Normal(mean, sd) on theta (log scale), truncated to [lo, hi] (default: none), so
# no Jacobian or dependence structure is needed.
prior_from_unit <- function(u, mean, sd, lo = -Inf, hi = Inf) {
  stopifnot(ncol(u) == length(mean), length(sd) %in% c(1L, length(mean)))
  sd <- rep_len(sd, length(mean)); lo <- rep_len(lo, length(mean)); hi <- rep_len(hi, length(mean))
  matrix(vapply(seq_along(mean), function(j) {  # matrix() keeps one-row input a matrix
    a <- pnorm(lo[j], mean[j], sd[j]); b <- pnorm(hi[j], mean[j], sd[j])
    qnorm(a + u[, j] * (b - a), mean[j], sd[j])
  }, numeric(nrow(u))), nrow(u))
}

sample_prior_lhs <- function(n, mean, sd, seed, lo = -Inf, hi = Inf)
  prior_from_unit(lhs_unit(n, length(mean), seed), mean, sd, lo, hi)

sample_prior_random <- function(n, mean, sd, seed, lo = -Inf, hi = Inf) {
  set.seed(seed)
  prior_from_unit(matrix(runif(n * length(mean)), n), mean, sd, lo, hi)
}
