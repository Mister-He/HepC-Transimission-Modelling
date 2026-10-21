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

clip_to_bounds <- function(theta) pmin(pmax(theta, THETA_LO), THETA_HI)

# Uninformative starts: theta = 0 (every contact scale 1, c = 1) and n_perturbed draws
# around it (log-scale SD sd_perturb), clipped to the bounds. Nothing comes from earlier fits.
make_start_set <- function(n_perturbed = 5, sd_perturb = 0.4, seed_base = 101) {
  anchor <- rep(0, N_THETA)
  starts <- list(list(id = "uninformative", theta0 = anchor, seed = seed_base))
  set.seed(seed_base)
  for (j in seq_len(n_perturbed)) {
    starts[[j + 1]] <- list(
      id = paste0("perturb_", j),
      theta0 = clip_to_bounds(anchor + rnorm(N_THETA, mean = 0, sd = sd_perturb)),
      seed = seed_base + j
    )
  }
  starts
}
