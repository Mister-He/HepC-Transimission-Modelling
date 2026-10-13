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

clip_to_bounds <- function(theta) {
  contact_scale <- exp(theta[1:6])
  beta_scale <- exp(theta[7:12])
  contact_scale <- pmin(pmax(contact_scale, CONTACT_LO), CONTACT_HI)
  beta_scale <- pmin(pmax(beta_scale, BETA_LO), BETA_HI)
  c(log(contact_scale), log(beta_scale))
}

plausible_anchors <- function() {
  list(
    contact = c(6.3713, 0.0857, 0.0589, 0.5584, 9.5232, 0.2558),
    beta    = c(0.1670, 0.9443, 1.2297, 5.2385, 16.9018, 110.5976)
  )
}

make_start_set <- function(n_perturbed = 5, sd_perturb = 0.4, seed_base = 101) {
  a <- plausible_anchors()
  anchor <- c(log(a$contact), log(a$beta))
  starts <- list(list(id = "plausible_anchor", theta0 = anchor, seed = seed_base))
  set.seed(seed_base)
  for (j in seq_len(n_perturbed)) {
    theta_j <- anchor + rnorm(N_THETA, mean = 0, sd = sd_perturb)
    starts[[j + 1]] <- list(
      id = paste0("perturb_", j),
      theta0 = clip_to_bounds(theta_j),
      seed = seed_base + j
    )
  }

  starts
}

# Deterministic population-informed start: sequential steady-state balance
# for D_i and J_i given the prison population targets, including ageing
# flows (10-year bands -> 1/10 per year). Returns a length-12 log-scale theta.
informed_theta0 <- function(base_params, targets) {
  # shortcut: assumes beta_scale = 1 in lambda1 = lambda3 * c_composite / beta_scale; it only
  # seeds one start, so Nelder-Mead corrects it. Upgrade to a fixed-point pass if that start matters.
  lambda1 <- base_params$lambda3 * base_params$c_composite
  lambda2 <- base_params$lambda2
  mu_eff  <- base_params$mu * base_params$omega
  pi      <- base_params$pi_recid
  J_target <- targets$population

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
