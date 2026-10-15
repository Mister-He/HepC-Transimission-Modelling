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

# Population-informed start: contact scales 1 and the entry that keeps each never-arrested
# pool at D_i = f_first_i / lambda3_i in steady state. Returns a length-12 theta.
informed_theta0 <- function(base_params) {
  # shortcut: assumes first arrest is as fast as re-arrest (lambda1 = lambda3); it only seeds
  # one start, so Nelder-Mead corrects it.
  D <- base_params$f_first / base_params$lambda3
  out <- base_params$mu * base_params$omega + c(rep(AGE_RATE, 5), 0)
  entry <- out * D + base_params$f_first - AGE_RATE * c(0, D[-6])
  c(rep(0, 6), excess_theta(pmax(entry, 0) / base_params$beta, base_params))
}
