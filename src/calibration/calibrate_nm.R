# =============================================================================
# calibrate_nm.R — deterministic multi-start Nelder-Mead runner
# =============================================================================

make_logged_objective <- function(fn) {
  env <- new.env(parent = emptyenv())
  env$n_eval <- 0L
  env$history <- vector("list", 0L)
  logged <- function(theta) {
    val <- fn(theta)
    env$n_eval <- env$n_eval + 1L
    env$history[[env$n_eval]] <- c(iteration = env$n_eval, theta, objective = val)
    val
  }
  attr(logged, "get_history") <- function() {
    if (env$n_eval == 0L) return(NULL)
    h <- do.call(rbind, env$history)
    colnames(h) <- c("iteration", paste0("theta", 1:N_THETA), "objective")
    h
  }
  logged
}

run_nm_start <- function(theta0, objective_fn, start_id = "start",
                        maxit = 3000, reltol = 1e-8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  logged <- make_logged_objective(objective_fn)
  t0 <- proc.time()
  fit <- optim(par = theta0, fn = logged, method = "Nelder-Mead",
               control = list(maxit = maxit, reltol = reltol, trace = 0))
  elapsed <- as.numeric((proc.time() - t0)["elapsed"])
  list(
    start_id       = start_id,
    seed           = seed,
    theta0         = theta0,
    theta          = fit$par,
    objective      = fit$value,
    convergence    = fit$convergence,
    fn_evaluations = fit$counts[["function"]],
    elapsed_sec    = elapsed,
    message        = fit$message,
    history        = attr(logged, "get_history")()
  )
}

make_start_set <- function(n_perturbed = 5, sd_perturb = 0.8, seed_base = 101) {
  starts <- list(list(id = "baseline_zeros", theta0 = rep(0, N_THETA), seed = seed_base))
  for (j in seq_len(n_perturbed)) {
    set.seed(seed_base + j)
    starts[[j + 1]] <- list(
      id = paste0("perturb_", j),
      theta0 = rnorm(N_THETA, mean = 0, sd = sd_perturb),
      seed = seed_base + j
    )
  }
  starts
}

# Deterministic population-informed start: sequential steady-state balance
# for D_i and J_i given the prison population targets J_i*, including ageing
# flows (10-year bands -> 1/10 per year). Returns a length-12 log-scale theta.
informed_theta0 <- function(base_params, targets) {
  lambda1 <- base_params$lambda1
  lambda2 <- base_params$lambda2
  mu_eff  <- base_params$mu * base_params$omega
  pi      <- base_params$pi_recid
  J_target <- targets$prison_total

  D <- numeric(6)
  J <- numeric(6)
  beta_need <- numeric(6)
  for (i in 1:6) {
    J[i] <- J_target[i]
    j_in <- if (i == 1) 0 else J[i - 1] / 10
    j_out_rate <- if (i == 6) 0 else 1 / 10
    d_in <- if (i == 1) 0 else D[i - 1] / 10
    d_out_rate <- if (i == 6) 0 else 1 / 10

    D[i] <- ((lambda2[i] + mu_eff[i] + j_out_rate) * J[i] - j_in) / lambda1[i]
    D[i] <- max(D[i], 1e-3)
    beta_need[i] <- (lambda1[i] + mu_eff[i] + d_out_rate) * D[i] -
                    pi * lambda2[i] * J[i] - d_in
  }
  beta_scale <- beta_need / base_params$beta
  beta_scale <- pmax(beta_scale, 1e-3)
  c(log(rep(1, 6)), log(beta_scale))
}

# Informed variant with the 60+ inflow scale lifted (the baseline informed
# start puts beta6 ~ 0, which traps the optimizer; lifting beta6 to ~90
# anchors the 60+ prison share).
informed_theta0_b6 <- function(base_params, targets, beta6_anchor = 90) {
  th <- informed_theta0(base_params, targets)
  th[12] <- log(beta6_anchor)
  th
}

# Informed variant with both the <20 contact row and the 60+ inflow lifted.
informed_theta0_b6_c1 <- function(base_params, targets,
                                  beta6_anchor = 90, c1_anchor = 8) {
  th <- informed_theta0(base_params, targets)
  th[1]  <- log(c1_anchor)
  th[12] <- log(beta6_anchor)
  th
}

# Warm start derived from the previous round's converged 15-parameter
# solution (run7_final_y1960, best start warm_run3): the fitted
# transmission multiplier was nearly flat (m ~ 0.57), so under constant
# transmission the contact row scales are rescaled by 0.57 and the beta
# scales are kept; the excess-mortality parameter is dropped. This basin
# (NLL ~22.5, all acceptance criteria met with the adopted progression
# rates) is not reachable from the population-informed or random starts
# (see DECISIONS.md "12-parameter warm start").
warm_run7_derived_12p <- c(
  log(5.777), log(0.082), log(0.053), log(0.465),
  log(7.5),   log(0.23),
  log(0.169), log(0.911), log(1.283), log(5.435),
  log(15.135), log(112.264)
)

# Converged 12-parameter solution reached from warm_run7_derived_12p with
# Nelder-Mead (maxit 6000, seed 101): objective 22.4675, equilibrium pass
# (max |prev change| 0.00124). Hardcoded for reproducible multi-start
# stability checking.
warm_12p_22p5 <- c(
  1.7529, -2.5014, -2.9376, -0.7657, 2.0149, -1.4697,
  -1.7778, -0.0932, 0.2492, 1.6931, 2.7173, 4.7208
)

run_stability_check <- function(theta_best, objective_fn,
                                n = 5, sd = 0.15, seed = 9091) {
  out <- vector("list", n)
  for (j in seq_len(n)) {
    set.seed(seed + j)
    theta_j <- theta_best + rnorm(length(theta_best), 0, sd)
    out[[j]] <- data.frame(
      perturbation = j,
      objective = objective_fn(theta_j),
      max_abs_theta_shift = max(abs(theta_j - theta_best))
    )
  }
  do.call(rbind, out)
}
