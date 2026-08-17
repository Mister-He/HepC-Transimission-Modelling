# =============================================================================
# likelihood.R — 12-parameter calibration with soft plausibility preferences
#
#   theta[1:6]  contact-matrix row scales (exp)
#   theta[7:12] beta inflow scales (exp)
#
# No excess-mortality parameter is fitted or used (eta_s = 1).
#
# Soft preferences (small penalty, not hard constraints):
#   - scaled contact-matrix diagonal should be unimodal, peak at 30-39
#     (or 20-29), and 50+/60+ diagonal entries should be relatively small;
#   - beta scales should be monotone non-decreasing with age and preferably
#     > 1.
# =============================================================================

sigma_pop <- 0.10
eps_prev  <- 1e-10
EQ_PENALTY_FACTOR <- 1e6
BOUND_PENALTY_FACTOR <- 1e6
PLAUSIBILITY_FACTOR <- 1e-2

N_THETA <- 12
TARGET_MODE <- "sero"

CONTACT_LO <- rep(0.01, 6)
CONTACT_HI <- rep(1000, 6)
BETA_LO <- rep(0.01, 6)
BETA_HI <- rep(1e6, 6)

within_bounds <- function(theta, tol = 1e-8) {
  contact_scale <- exp(theta[1:6])
  beta_scale <- exp(theta[7:12])
  all(contact_scale >= CONTACT_LO - tol & contact_scale <= CONTACT_HI + tol) &&
    all(beta_scale >= BETA_LO - tol & beta_scale <= BETA_HI + tol)
}

bound_barrier <- function(theta) {
  contact_scale <- exp(theta[1:6])
  beta_scale <- exp(theta[7:12])
  sum(pmax(log(CONTACT_LO) - log(contact_scale), 0)^2) +
    sum(pmax(log(contact_scale) - log(CONTACT_HI), 0)^2) +
    sum(pmax(log(BETA_LO) - log(beta_scale), 0)^2) +
    sum(pmax(log(beta_scale) - log(BETA_HI), 0)^2)
}

pattern_penalty <- function(pm, theta) {
  diag <- diag(pm$C_contact)
  pen <- 0
  # unimodal diagonal, peak at 30-39
  pen <- pen + pmax(log(diag[2]) - log(diag[3]), 0) + pmax(log(diag[1]) - log(diag[2]), 0)
  pen <- pen + pmax(log(diag[4]) - log(diag[3]), 0) + pmax(log(diag[5]) - log(diag[4]), 0) +
    pmax(log(diag[6]) - log(diag[5]), 0)
  # 50+/60+ should stay relatively small: below half of the 30-39 peak
  pen <- pen + pmax(log(diag[5]) - log(0.5 * diag[3]), 0) +
    pmax(log(diag[6]) - log(0.5 * diag[3]), 0)
  # beta monotone non-decreasing with age, values above 1 preferred
  beta_scale <- exp(theta[7:12])
  for (i in 1:5) pen <- pen + pmax(log(beta_scale[i]) - log(beta_scale[i + 1]), 0)
  pen + 0.01 * sum(pmax(1 - beta_scale, 0)^2)
}

build_params <- function(theta, base_params) {
  stopifnot(length(theta) == N_THETA)
  contact_scale <- exp(theta[1:6])
  beta_scale <- exp(theta[7:12])

  pm <- base_params
  for (i in 1:6) {
    pm$C_contact[i, ] <- base_params$C_contact[i, ] * contact_scale[i]
  }
  pm$beta <- base_params$beta * beta_scale
  pm
}

nll_prev <- function(p_hat, x = cal_targets$x_prev, n = cal_targets$n_prev) {
  p_safe <- pmin(pmax(p_hat, eps_prev), 1 - eps_prev)
  -sum(dbinom(x, size = n, prob = p_safe, log = TRUE))
}

nll_pop <- function(N_hat, N_obs = cal_targets$prison_total, sigma = sigma_pop) {
  0.5 * sum(((log(N_hat) - log(N_obs)) / sigma)^2)
}

make_objective <- function(base_params, data, target_mode = TARGET_MODE,
                           equilibrium_penalty = TRUE) {
  force(base_params); force(data); force(target_mode); force(equilibrium_penalty)
  function(theta) {
    if (!within_bounds(theta)) {
      return(1e12 + BOUND_PENALTY_FACTOR * bound_barrier(theta))
    }

    pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
    if (is.null(pm)) return(1e12)

    out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
    if (is.null(out) || nrow(out) < 2) return(1e12)

    s <- J_summary_final(out)
    p_hat <- if (target_mode == "sero") s$p_sero else s$p_viremic
    if (!all(is.finite(s$N_hat)) || !all(is.finite(p_hat)) ||
        any(s$N_hat <= 0)) return(1e12)

    nll <- nll_prev(p_hat) + nll_pop(s$N_hat) +
      PLAUSIBILITY_FACTOR * pattern_penalty(pm, theta)
    if (!is.finite(nll)) return(1e12)

    if (equilibrium_penalty) {
      eq <- tryCatch(check_equilibrium(out, target_mode = target_mode),
                     error = function(e) NULL)
      if (!is.null(eq) && !eq$pass) {
        pen_log  <- pmax(eq$max_log_pop_ratio / 0.01, 0)
        pen_prev <- pmax(eq$max_prev_change   / 0.005, 0)
        pen_tot  <- pmax(eq$total_log_ratio   / 0.01, 0)
        pen_state <- pmax(eq$state_log_ratio   / 0.01, 0)
        pen_comp  <- pmax(eq$max_comp_log_ratio / 0.02, 0)
        nll <- nll + EQ_PENALTY_FACTOR * (pen_log + pen_prev + pen_tot +
                                          pen_state + pen_comp)
      }
    }
    nll
  }
}
