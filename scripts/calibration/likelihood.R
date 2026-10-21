# =============================================================================
# likelihood.R — likelihood, prior, parameter bounds and parameterisation for the six
# contact scales and the first-arrest multiplier c used by every calibration method.
# =============================================================================

N_THETA <- 7L
TARGET_MODE <- "sero"
EPS <- 1e-10
EQ_PENALTY_FACTOR <- 1e6
BOUND_PENALTY_FACTOR <- 1e6
CONTACT_LO <- rep(0.01, 6); CONTACT_HI <- rep(1000, 6)
# theta[7] = log c: first arrest lambda1_i = c * lambda3_i. Entry into D is the CNB count,
# unscaled. c >= 0.5 (user decision: first arrest at least half the re-arrest rate).
C_LO <- 0.5; C_HI <- 100
THETA_LO <- log(c(CONTACT_LO, C_LO)); THETA_HI <- log(c(CONTACT_HI, C_HI))

# Fixed prior, shared by Nelder-Mead (posterior mode), the reference MCMC and NPE. It does not
# depend on the data: independent Normal on theta, truncated to the bounds.
# - theta1-6 (log contact scales relative to setup.R's C_contact): Normal(0, 1.5); 95% range
#   of a scale 0.05-19.
# - theta7 (log c): Normal(0, 0.5), centred on first arrest as fast as re-arrest. Seroprevalence
#   barely informs c, so this prior largely sets it.
PRIOR_MEAN <- rep(0, N_THETA)
PRIOR_SD <- c(rep(1.5, 6), 0.5)
prior_spec <- function() list(mean = PRIOR_MEAN, sd = PRIOR_SD, lo = THETA_LO, hi = THETA_HI)
nlp <- function(theta) 0.5 * sum(((theta - PRIOR_MEAN) / PRIOR_SD)^2)  # -log prior + const

within_bounds <- function(theta) all(theta >= THETA_LO & theta <= THETA_HI)

bound_barrier <- function(theta) sum(pmax(THETA_LO - theta, 0)^2 + pmax(theta - THETA_HI, 0)^2)

build_params <- function(theta, base_params) {
  stopifnot(length(theta) == N_THETA)
  pm <- base_params
  pm$C_contact <- base_params$C_contact * exp(theta[1:6])
  pm$c_first <- exp(theta[7])
  pm
}

# Seroprevalence is the only target, with no counts: each age group's observed value y_i is
# Beta(p_i * phi, (1 - p_i) * phi) around the model value p_i, with fixed precision phi
# (targets$beta_precision; variance p(1 - p) / (phi + 1)).
nll_prev <- function(p_hat) {
  p <- pmin(pmax(p_hat, EPS), 1 - EPS); phi <- targets$beta_precision
  -sum(dbeta(targets$seroprevalence, p * phi, (1 - p) * phi, log = TRUE))
}

# Nelder-Mead objective: negative log posterior (likelihood + prior), so the best fit is the
# posterior mode; failed equilibrium adds a large penalty.
make_objective <- function(base_params, data, equilibrium_penalty = TRUE) {
  function(theta) {
    if (!within_bounds(theta)) return(1e12 + BOUND_PENALTY_FACTOR * bound_barrier(theta))
    pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
    out <- if (!is.null(pm)) tryCatch(run_sim(pm, data), error = function(e) NULL) else NULL
    if (is.null(out)) return(1e12)
    s <- J_summary_final(out)
    if (any(!is.finite(s$N_hat)) || any(s$N_hat <= 0) || any(!is.finite(s$p_sero)))
      return(1e12)
    value <- nll_prev(s$p_sero) + nlp(theta)
    if (equilibrium_penalty) {
      eq <- tryCatch(check_equilibrium(out), error = function(e) NULL)
      if (is.null(eq)) return(1e12)
      if (!eq$pass) value <- value + EQ_PENALTY_FACTOR * eq$violation
    }
    value
  }
}
