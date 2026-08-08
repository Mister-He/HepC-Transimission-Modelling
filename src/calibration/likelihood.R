# =============================================================================
# likelihood.R — parameterisation (12 fitted parameters) and combined NLL
#
#   theta[1:6]  contact-matrix row scales     (exp)
#   theta[7:12] beta inflow scales            (exp)
#
# Constant transmission (m_min = m_max = 1, NOT fitted; merged into the
# contact scales). No excess-mortality parameter by default: FIT_ETA_S6 is
# FALSE. Only if the 60+ age group cannot be fitted with the 12 parameters
# may FIT_ETA_S6 be switched to TRUE (adds theta[13] = eta_s[6] via
# 1 + 4*plogis(theta13), bounded to [1,5]) — see AGENTS.md "Excess-mortality
# contingency" and the dedicated rationale file.
#
# Likelihood:
#   x_prev[i] ~ Binomial(n_prev[i], p_hat[i])    p_hat = seroprevalence
#   log(N_obs[i]) ~ Normal(log(N_hat[i]), 0.10^2)
#   NLL = nll_prev + nll_pop
# =============================================================================

sigma_pop <- 0.10
eps_prev  <- 1e-10
EQ_PENALTY_FACTOR <- 1e6

N_THETA <- 12
FIT_ETA_S6 <- FALSE
TARGET_MODE <- "sero"

build_params <- function(theta, base_params) {
  stopifnot(length(theta) == N_THETA)
  contact_scale <- exp(theta[1:6])
  beta_scale    <- exp(theta[7:12])

  pm <- base_params
  for (i in 1:6) {
    pm$C_contact[i, ] <- base_params$C_contact[i, ] * contact_scale[i]
  }
  pm$beta <- base_params$beta * beta_scale

  # Constant transmission: m_min/m_max = 1 from setup.R; not fitted here.
  # Contact row scales absorb the historical transmission level.
  if (FIT_ETA_S6) {
    pm$eta_s[6] <- 1 + 4 * plogis(theta[13])
  }
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
    pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
    if (is.null(pm)) return(1e12)

    out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
    if (is.null(out) || nrow(out) < 2) return(1e12)

    s <- J_summary_final(out)
    p_hat <- if (target_mode == "sero") s$p_sero else s$p_viremic
    if (!all(is.finite(s$N_hat)) || !all(is.finite(p_hat)) ||
        any(s$N_hat <= 0)) return(1e12)

    nll <- nll_prev(p_hat) + nll_pop(s$N_hat)
    if (!is.finite(nll)) return(1e12)

    # Equilibrium feasibility gate (during optimisation only; the Laplace
    # Hessian uses the pure statistical likelihood).
    if (equilibrium_penalty) {
      eq <- tryCatch(check_equilibrium(out, target_mode = target_mode),
                     error = function(e) NULL)
      if (!is.null(eq) && !eq$pass) {
        pen_log  <- pmax(eq$max_log_pop_ratio / 0.01, 0)
        pen_prev <- pmax(eq$max_prev_change   / 0.005, 0)
        pen_tot  <- pmax(eq$total_log_ratio   / 0.01, 0)
        nll <- nll + EQ_PENALTY_FACTOR * (pen_log + pen_prev + pen_tot)
      }
    }
    nll
  }
}
