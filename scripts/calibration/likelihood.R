# =============================================================================
# likelihood.R — likelihood, parameter bounds and parameterisation for the
# six contact and six inflow scales used by every calibration method.
# =============================================================================

N_THETA <- 12L
TARGET_MODE <- "sero"
EPS <- 1e-10
EQ_PENALTY_FACTOR <- 1e6
BOUND_PENALTY_FACTOR <- 1e6
PLAUSIBILITY_FACTOR <- 1e-2
CONTACT_LO <- rep(0.01, 6); CONTACT_HI <- rep(1000, 6)
BETA_LO <- rep(0.01, 6); BETA_HI <- rep(1e6, 6)

# Prior centre shared by the reference MCMC and NPE: the Nelder-Mead best fit, moved to at
# least PRIOR_MARGIN_SD prior SDs inside the bounds. The fit can end on a bound where the
# likelihood is flat (the 60+ contact scale), which would put half the prior out of bounds.
PRIOR_MARGIN_SD <- 5
prior_centre <- function(theta, sd = 0.20) {
  m <- PRIOR_MARGIN_SD * sd
  pmin(pmax(theta, log(c(CONTACT_LO, BETA_LO)) + m), log(c(CONTACT_HI, BETA_HI)) - m)
}

within_bounds <- function(theta) {
  cs <- exp(theta[1:6]); bs <- exp(theta[7:12])
  all(cs >= CONTACT_LO & cs <= CONTACT_HI & bs >= BETA_LO & bs <= BETA_HI)
}

bound_barrier <- function(theta) {
  cs <- exp(theta[1:6]); bs <- exp(theta[7:12])
  sum(pmax(log(CONTACT_LO / cs), 0)^2 + pmax(log(cs / CONTACT_HI), 0)^2) +
    sum(pmax(log(BETA_LO / bs), 0)^2 + pmax(log(bs / BETA_HI), 0)^2)
}

pattern_penalty <- function(pm, theta) {
  d <- diag(pm$C_contact); bs <- exp(theta[7:12])
  peak <- pmax(log(d[2] / d[3]), 0) + pmax(log(d[1] / d[2]), 0) +
    pmax(log(d[4] / d[3]), 0) + pmax(log(d[5] / d[4]), 0) +
    pmax(log(d[6] / d[5]), 0)
  peak + pmax(log(d[5] / (0.5 * d[3])), 0) +
    pmax(log(d[6] / (0.5 * d[3])), 0) +
    sum(pmax(log(bs[1:5] / bs[2:6]), 0)) + 0.01 * sum(pmax(1 - bs, 0)^2)
}

build_params <- function(theta, base_params) {
  stopifnot(length(theta) == N_THETA)
  cs <- exp(theta[1:6]); bs <- exp(theta[7:12])
  pm <- base_params
  pm$C_contact <- base_params$C_contact * cs
  pm$beta <- base_params$beta * bs
  pm$beta_base <- base_params$beta
  pm$beta_scale <- bs
  pm
}

nll_prev <- function(p_hat) {
  p <- pmin(pmax(p_hat, EPS), 1 - EPS)
  -sum(dbinom(targets$x_prev, targets$population, p, log = TRUE))
}

nll_pop <- function(N_hat) {
  0.5 * sum(((log(N_hat) - log(targets$population)) /
               targets$population_sigma_log)^2)
}

make_objective <- function(base_params, data, equilibrium_penalty = TRUE) {
  function(theta) {
    if (!within_bounds(theta)) return(1e12 + BOUND_PENALTY_FACTOR * bound_barrier(theta))
    pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
    out <- if (!is.null(pm)) tryCatch(run_sim(pm, data), error = function(e) NULL) else NULL
    if (is.null(out)) return(1e12)
    s <- J_summary_final(out)
    if (any(!is.finite(s$N_hat)) || any(s$N_hat <= 0) || any(!is.finite(s$p_sero)))
      return(1e12)
    value <- nll_prev(s$p_sero) + nll_pop(s$N_hat) +
      PLAUSIBILITY_FACTOR * pattern_penalty(pm, theta)
    if (equilibrium_penalty) {
      eq <- tryCatch(check_equilibrium(out), error = function(e) NULL)
      if (is.null(eq)) return(1e12)
      if (!eq$pass) value <- value + EQ_PENALTY_FACTOR * eq$violation
    }
    value
  }
}
