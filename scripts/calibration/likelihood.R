# =============================================================================
# likelihood.R — likelihood, parameter bounds and parameterisation for the
# six contact scales and six excess-entry scales used by every calibration method.
# =============================================================================

N_THETA <- 12L
TARGET_MODE <- "sero"
EPS <- 1e-10
EQ_PENALTY_FACTOR <- 1e6
BOUND_PENALTY_FACTOR <- 1e6
PLAUSIBILITY_FACTOR <- 1e-2
CONTACT_LO <- rep(0.01, 6); CONTACT_HI <- rep(1000, 6)
# Bounds on the excess-entry scale exp(theta[7:12]) (see entry_scales()).
BETA_LO <- rep(1e-4, 6); BETA_HI <- rep(1e6, 6)
AGE_RATE <- 0.1  # sim.cpp AGE_RATE

# theta[7:12] = log(excess entry / beta_i). Entry into age i is the minimum that supplies the
# first-arrest flow, max(0, f_i - 0.1*D_{i-1}), plus that excess, so every theta is feasible
# for first_arrest_rate() in sim.cpp. D follows the same HCV-free steady state as sim.cpp.
# Returns the entry scale entry_i / beta_i applied to setup.R's beta.
entry_scales <- function(theta, base_params) {
  m <- base_params$mu * base_params$omega; f <- base_params$f_first; b <- base_params$beta
  bs <- numeric(6); d_prev <- 0
  for (i in 1:6) {
    entry <- max(0, f[i] - AGE_RATE * d_prev) + exp(theta[6 + i]) * b[i]
    d_prev <- (entry + AGE_RATE * d_prev - f[i]) / (m[i] + if (i < 6) AGE_RATE else 0)
    bs[i] <- entry / b[i]
  }
  bs
}

# Inverse of entry_scales() for starts given as entry scales; an entry below the minimum
# becomes the smallest excess (BETA_LO).
excess_theta <- function(beta_scale, base_params) {
  m <- base_params$mu * base_params$omega; f <- base_params$f_first; b <- base_params$beta
  th <- numeric(6); d_prev <- 0
  for (i in 1:6) {
    need <- max(0, f[i] - AGE_RATE * d_prev)
    th[i] <- log(max((beta_scale[i] * b[i] - need) / b[i], BETA_LO[i]))
    entry <- need + exp(th[i]) * b[i]
    d_prev <- (entry + AGE_RATE * d_prev - f[i]) / (m[i] + if (i < 6) AGE_RATE else 0)
  }
  th
}

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
  d <- diag(pm$C_contact); bs <- pm$beta / pm$beta_base
  peak <- pmax(log(d[2] / d[3]), 0) + pmax(log(d[1] / d[2]), 0) +
    pmax(log(d[4] / d[3]), 0) + pmax(log(d[5] / d[4]), 0) +
    pmax(log(d[6] / d[5]), 0)
  peak + pmax(log(d[5] / (0.5 * d[3])), 0) +
    pmax(log(d[6] / (0.5 * d[3])), 0) +
    sum(pmax(log(bs[1:5] / bs[2:6]), 0)) + 0.01 * sum(pmax(1 - bs, 0)^2)
}

build_params <- function(theta, base_params) {
  stopifnot(length(theta) == N_THETA)
  cs <- exp(theta[1:6])
  pm <- base_params
  pm$C_contact <- base_params$C_contact * cs
  pm$beta <- base_params$beta * entry_scales(theta, base_params)
  pm$beta_base <- base_params$beta
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

# 2015 entry into D (pm$beta) against the CNB count: log-normal, SD entry_sigma_log. Without it
# the never-arrested pool is unidentified (first arrests are fixed by f_first).
nll_entry <- function(entry_hat) {
  0.5 * sum(((log(entry_hat) - log(targets$entry)) / targets$entry_sigma_log)^2)
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
    value <- nll_prev(s$p_sero) + nll_pop(s$N_hat) + nll_entry(pm$beta) +
      PLAUSIBILITY_FACTOR * pattern_penalty(pm, theta)
    if (equilibrium_penalty) {
      eq <- tryCatch(check_equilibrium(out), error = function(e) NULL)
      if (is.null(eq)) return(1e12)
      if (!eq$pass) value <- value + EQ_PENALTY_FACTOR * eq$violation
    }
    value
  }
}
