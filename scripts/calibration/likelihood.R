# =============================================================================
# likelihood.R — likelihood, prior, parameter bounds and parameterisation for the six
# contact scales, the six entry scales and the first-arrest multiplier c used by every
# calibration method.
# =============================================================================

N_THETA <- 13L
IDX_CONTACT <- 1:6; IDX_ENTRY <- 7:12; IDX_C <- 13L
EPS <- 1e-10
EQ_PENALTY_FACTOR <- 1e6
BOUND_PENALTY_FACTOR <- 1e6
CONTACT_LO <- rep(0.01, 6); CONTACT_HI <- rep(1000, 6)
# theta7-12 = log entry scales: entry into D_i = s_i * CNB 2015 new drug abusers_i (setup.R beta).
# s_i >= 1 (user decision: the CNB counts are a lower bound on entry into injecting).
ENTRY_LO <- rep(1, 6); ENTRY_HI <- rep(100, 6)
# theta13 = log c: first arrest lambda1_i = c * lambda3_i. c >= 0.5 (user decision: first arrest
# at least half the re-arrest rate).
C_LO <- 0.5; C_HI <- 100
THETA_LO <- log(c(CONTACT_LO, ENTRY_LO, C_LO)); THETA_HI <- log(c(CONTACT_HI, ENTRY_HI, C_HI))

# Fixed prior, shared by Nelder-Mead (posterior mode), the reference MCMC and NPE. It does not
# depend on the data: independent Normal on theta, truncated to the bounds.
# - theta1-6 (log contact scales relative to setup.R's C_contact): Normal(0, 1.5); 95% range
#   of a scale 0.05-19.
# - theta7-12 (log entry scales): Normal(0, 1) truncated at 0, i.e. half-normal; median scale
#   1.4, 95% below 7.1.
# - theta13 (log c): Normal(0, 0.5), centred on first arrest as fast as re-arrest.
PRIOR_MEAN <- rep(0, N_THETA)
PRIOR_SD <- c(rep(1.5, 6), rep(1, 6), 0.5)
prior_spec <- function() list(mean = PRIOR_MEAN, sd = PRIOR_SD, lo = THETA_LO, hi = THETA_HI)
nlp <- function(theta) 0.5 * sum(((theta - PRIOR_MEAN) / PRIOR_SD)^2)  # -log prior + const

within_bounds <- function(theta) all(theta >= THETA_LO & theta <= THETA_HI)

bound_barrier <- function(theta) sum(pmax(THETA_LO - theta, 0)^2 + pmax(theta - THETA_HI, 0)^2)

build_params <- function(theta, base_params) {
  stopifnot(length(theta) == N_THETA)
  pm <- base_params
  pm$C_contact <- base_params$C_contact * exp(theta[IDX_CONTACT])
  pm$beta <- base_params$beta * exp(theta[IDX_ENTRY])
  pm$c_first <- exp(theta[IDX_C])
  pm
}

# Three independent likelihood terms (targets in data/calibration_targets.json):
# - seroprevalence: each age group's observed value y_i is Beta(p_i * phi, (1 - p_i) * phi)
#   around the model value p_i, with fixed precision phi (variance p(1 - p) / (phi + 1));
# - first arrests: the CNB 2015 new drug abusers by age are negative binomial around the model
#   D -> J flow, with fixed size (CV about 1/sqrt(size) for large counts; handles the 60+ count
#   of 1);
# - community PWID: log(estimate) is Normal around log(model D + F), SD from the 95% interval.
nll_prev <- function(p_hat) {
  p <- pmin(pmax(p_hat, EPS), 1 - EPS); phi <- targets$beta_precision
  -sum(dbeta(targets$seroprevalence, p * phi, (1 - p) * phi, log = TRUE))
}
nll_first_arrests <- function(m)
  -sum(dnbinom(targets$first_arrests, size = targets$first_arrests_nb_size, mu = pmax(m, EPS), log = TRUE))
nll_community <- function(m)
  -dnorm(log(targets$community_pwid), log(max(m, EPS)), targets$community_pwid_log_sd, log = TRUE)
nll_targets <- function(cs)
  nll_prev(cs$p_sero) + nll_first_arrests(cs$first_arrests) + nll_community(cs$community_pwid)

# Finite, positive model summaries (otherwise the parameter set is outside the support).
summary_ok <- function(out, cs) {
  s <- J_summary_final(out)
  all(is.finite(s$N_hat)) && all(s$N_hat > 0) && all(is.finite(cs$p_sero)) &&
    all(is.finite(cs$first_arrests)) && is.finite(cs$community_pwid) && cs$community_pwid > 0
}

# Nelder-Mead objective: negative log posterior (likelihood + prior), so the best fit is the
# posterior mode; failed equilibrium adds a large penalty.
make_objective <- function(base_params, data, equilibrium_penalty = TRUE) {
  function(theta) {
    if (!within_bounds(theta)) return(1e12 + BOUND_PENALTY_FACTOR * bound_barrier(theta))
    pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
    out <- if (!is.null(pm)) tryCatch(run_sim(pm, data), error = function(e) NULL) else NULL
    if (is.null(out)) return(1e12)
    cs <- calib_summary(out, pm)
    if (!summary_ok(out, cs)) return(1e12)
    value <- nll_targets(cs) + nlp(theta)
    if (equilibrium_penalty) {
      eq <- tryCatch(check_equilibrium(out), error = function(e) NULL)
      if (is.null(eq)) return(1e12)
      if (!eq$pass) value <- value + EQ_PENALTY_FACTOR * eq$violation
    }
    value
  }
}

# Stop if a stored Nelder-Mead fit no longer matches the current model: the parameters must equal
# setup.R's and re-simulating the best fit must reproduce its stored seroprevalence, which catches
# changes to sim.cpp as well as to setup.R.
check_fit_current <- function(fit, base_params, data) {
  stopifnot("fit made with different setup.R parameters" = isTRUE(all.equal(fit$base_params, base_params)))
  d <- data; d$t_start <- 0; d$t_end <- targets$equilibrium_horizon_years
  d$dt <- fit$data$dt; d$save_every <- round(1 / d$dt)  # the fit's own step, whatever the caller uses
  p <- J_summary_final(run_sim(build_params(fit$best$theta, base_params), d))$p_sero
  stopifnot("fit made with a different model (sim.cpp)" =
              isTRUE(all.equal(p, fit$predictions$seroprevalence_hat, tolerance = 1e-8)))
  invisible(TRUE)
}
