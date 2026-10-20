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

# Restarts rerun Nelder-Mead from its own solution (a fresh simplex) while the objective still
# falls by more than 1e-3: the simplex can collapse early on this objective.
run_nm_start <- function(theta0, objective_fn, start_id = "start",
                        maxit = 3000, reltol = 1e-8, seed = NULL, restarts = 0) {
  if (!is.null(seed)) set.seed(seed)
  logged <- make_logged_objective(objective_fn)
  t0 <- proc.time()
  nm <- function(par) optim(par = par, fn = logged, method = "Nelder-Mead",
                            control = list(maxit = maxit, reltol = reltol, trace = 0))
  fit <- nm(theta0); evals <- fit$counts[["function"]]
  for (k in seq_len(restarts)) {
    nxt <- nm(fit$par); evals <- evals + nxt$counts[["function"]]
    better <- nxt$value < fit$value - 1e-3
    if (nxt$value <= fit$value) fit <- nxt
    if (!better) break
  }
  fit$counts[["function"]] <- evals
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

# Population-informed start: contact scales 1 and the entry scales of the HCV-free steady state
# that holds the prison at its 2015 targets. J_i gives F_i (pi*lambda2 releases plus ageing), the
# first-arrest flow f_i = J outflow - ageing in - lambda3*F_i, D_i = f_i / lambda1_i, and the D
# balance gives entry. Returns a length-12 theta; Nelder-Mead corrects the shortcut.
informed_theta0 <- function(base_params, targets) {
  l1 <- base_params$lambda1; l2 <- base_params$lambda2; l3 <- base_params$lambda3
  m <- base_params$mu * base_params$omega; pi <- base_params$pi_recid; J <- targets$population
  a <- c(rep(0.1, 5), 0)  # ageing out of each band (none out of 60+)
  D <- Fp <- entry <- numeric(6)
  for (i in 1:6) {
    J_in <- if (i > 1) 0.1 * J[i - 1] else 0
    F_in <- if (i > 1) 0.1 * Fp[i - 1] else 0
    D_in <- if (i > 1) 0.1 * D[i - 1] else 0
    Fp[i] <- (pi * l2[i] * J[i] + F_in) / (l3[i] + m[i] + a[i])
    f <- max((l2[i] + m[i] + a[i]) * J[i] - J_in - l3[i] * Fp[i], 1e-3)
    D[i] <- f / l1[i]
    entry[i] <- (l1[i] + m[i] + a[i]) * D[i] - D_in
  }
  c(rep(0, 6), log(pmax(entry, 1e-3) / base_params$beta))
}
