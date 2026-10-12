#!/usr/bin/env Rscript
# =============================================================================
# run_mcmc_reference.R — reference MCMC calibration: adaptive random-walk Metropolis
# on the exact posterior (4 prior-started chains), with convergence diagnostics.
#
# Target (shared with single-round NPE): independent Normal(fit$best$theta, prior_sd)
# prior on the 12 log-scale parameters, Binomial prevalence + log-normal population
# likelihood, and the 2015 equilibrium check as a hard restriction. Defaults reproduce
# the production run: 1 init slot + 4,999 warmup + 45,000 retained proposals per chain.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) default else args[i + 1]
}
T_SCRIPT0 <- proc.time()[["elapsed"]]
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
FIT <- normalizePath(arg("--fit", file.path(ROOT,
  "output/calibration/r_20261010_agepop/fit.rds")), mustWork = TRUE)
OUT <- file.path(ROOT, arg("--out-dir", "output/calibration/mcmc_reference"))
N_CHAINS <- as.integer(arg("--n-chains", "4"))
CORES <- as.integer(arg("--cores", "4"))
N_WARMUP <- as.integer(arg("--n-warmup", "4999"))
MAX_EVALS <- as.integer(arg("--max-evals", "50000"))
SEED <- as.integer(arg("--seed", "20261101"))
PRIOR_SD <- as.numeric(arg("--prior-sd", "0.20"))
INIT_SCALE <- as.numeric(arg("--init-scale", "0.1"))
N_PRIOR_VAR <- as.integer(arg("--n-prior-var", "2048"))
DT <- as.numeric(arg("--dt", as.character(1 / 365)))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

source(file.path(ROOT, "scripts/calibration/targets.R"))
targets <- read_targets(file.path(ROOT, "data/calibration_targets.json"))
old <- setwd(file.path(ROOT, "scripts/model"))
setup_env <- new.env(parent = globalenv())
sys.source("setup.R", setup_env)
setwd(old)
Rcpp::sourceCpp(file.path(ROOT, "scripts/model/sim.cpp"))
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R", "mcmc.R"))
  sys.source(file.path(ROOT, "scripts/calibration", f), envir = environment())
fit <- readRDS(FIT)
base_params <- setup_env$params
data <- setup_env$data
data$t_start <- 0
data$t_end <- targets$equilibrium_horizon_years
data$dt <- DT
startup_sec <- proc.time()[["elapsed"]] - T_SCRIPT0

# ---- sampling ----
T_METHOD0 <- proc.time()[["elapsed"]]
prior_mean <- as.numeric(fit$best$theta)
prior_sd <- rep(PRIOR_SD, 12)
logpost <- function(theta) log_posterior_exact(theta, base_params, data, prior_mean, prior_sd)

# Transformed-prior variances from prior draws (no ODE solves). theta is already
# unconstrained and the prior is defined on theta, so no Jacobian is needed.
set.seed(SEED + 3000)
prior_draws <- sweep(matrix(rnorm(N_PRIOR_VAR * 12), N_PRIOR_VAR, 12), 2, prior_sd, "*") +
  matrix(prior_mean, N_PRIOR_VAR, 12, byrow = TRUE)
prior_var <- apply(prior_draws, 2, stats::var)
Sigma_init <- diag(INIT_SCALE^2 * prior_var)
reg <- diag(1e-6 * prior_var)
pre_chain_sec <- proc.time()[["elapsed"]] - T_METHOD0

run_chain <- function(i) {
  set.seed(SEED + 2000 + i)
  draw_prior <- function() rnorm(12, prior_mean, prior_sd)
  res <- arwm_budget(logpost, draw_prior, Sigma_init, reg, n_warmup = N_WARMUP,
                     max_evals = MAX_EVALS,
                     seed = SEED + 1000 + i, t_offset = pre_chain_sec)
  res$chain <- i
  res
}
chains <- if (CORES > 1 && .Platform$OS.type != "windows") {
  parallel::mclapply(seq_len(N_CHAINS), run_chain, mc.cores = min(CORES, N_CHAINS),
                     mc.set.seed = FALSE)
} else lapply(seq_len(N_CHAINS), run_chain)
method_sec <- proc.time()[["elapsed"]] - T_METHOD0

bad <- vapply(chains, function(ch) !is.list(ch) || !identical(ch$status, "ok"), logical(1))
chain_summary <- do.call(rbind, lapply(seq_along(chains), function(i) {
  ch <- chains[[i]]
  if (!is.list(ch)) return(data.frame(chain = i, status = "error", stringsAsFactors = FALSE))
  data.frame(chain = i, status = ch$status, init_tries = ch$init_tries %||% NA,
             warmup_acceptance = ch$warmup_acceptance %||% NA,
             retained_total = if (!is.null(ch$info)) nrow(ch$info) else 0L,
             acceptance = if (!is.null(ch$info)) mean(ch$info[, "accepted"]) else NA,
             total_evals = ch$total_evals %||% ch$evals, total_ode = ch$total_ode %||% ch$ode,
             warmup_sec = ch$warmup_sec %||% NA, total_sec = ch$total_sec %||% NA,
             stringsAsFactors = FALSE)
}))
write.csv(chain_summary, file.path(OUT, "mcmc_chain_summary.csv"), row.names = FALSE)
if (any(bad)) {
  writeLines(paste("FAIL: chain(s)", paste(which(bad), collapse = ","), "did not produce retained output"),
             file.path(OUT, "MCMC_STATUS.txt"))
  stop("MCMC chain failure; see mcmc_chain_summary.csv")
}

draws <- do.call(rbind, lapply(chains, function(ch) {
  x <- as.data.frame(ch$samples); names(x) <- paste0("theta", 1:12)
  cbind(data.frame(chain = ch$chain, iteration = seq_len(nrow(x))), x, as.data.frame(ch$info))
}))
write.csv(draws, file.path(OUT, "mcmc_draws.csv"), row.names = FALSE)
warm <- do.call(rbind, lapply(chains, function(ch) {
  x <- as.data.frame(ch$warmup); names(x) <- paste0("theta", 1:12)
  cbind(data.frame(chain = ch$chain, warmup_iteration = seq_len(nrow(x)) - 1L), x)
}))
write.csv(warm, file.path(OUT, "mcmc_warmup.csv"), row.names = FALSE)

# Rank-normalised split R-hat, bulk/tail ESS and MCSE (posterior package).
diag_for <- function(sel, label) {
  d <- draws[sel, ]
  n_min <- min(table(d$chain))
  arr <- array(NA_real_, c(n_min, N_CHAINS, 12))
  for (k in seq_len(N_CHAINS)) arr[, k, ] <- as.matrix(d[d$chain == k, paste0("theta", 1:12)][seq_len(n_min), ])
  do.call(rbind, lapply(1:12, function(j) {
    x <- arr[, , j]
    data.frame(subset = label, parameter = paste0("theta", j), draws_per_chain = n_min,
               mean = mean(x), sd = sd(as.vector(x)),
               rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x),
               ess_tail = posterior::ess_tail(x), mcse_mean = posterior::mcse_mean(x),
               mcse_q025 = posterior::mcse_quantile(x, 0.025),
               mcse_q975 = posterior::mcse_quantile(x, 0.975))
  }))
}
diagnostics <- diag_for(rep(TRUE, nrow(draws)), "all_retained")
write.csv(diagnostics, file.path(OUT, "mcmc_diagnostics.csv"), row.names = FALSE)
write.csv(posterior_summaries(as.matrix(draws[, paste0("theta", 1:12)])),
          file.path(OUT, "mcmc_parameter_summary.csv"), row.names = FALSE)

jsonlite::write_json(list(
  fit = FIT, sampler = "adaptive random-walk Metropolis (arwm_budget, scripts/calibration/mcmc.R)",
  n_chains = N_CHAINS, cores = CORES, n_warmup = N_WARMUP, max_evals_per_chain = MAX_EVALS,
  adapt_start = 200, adapt_every = 50, safety_p = 0.05,
  adaptive_scale = "2.38^2/d", init_scale = INIT_SCALE, n_prior_var_draws = N_PRIOR_VAR,
  regularisation = "1e-6 * diag(prior variance)", thinning = 1,
  prior = sprintf("independent Normal(fit$best$theta, %s) on theta", PRIOR_SD),
  target = "Binomial prevalence + log-normal population likelihood; equilibrium failure = -Inf",
  integration_dt = DT, seed = SEED,
  seeds = list(prior_variance = SEED + 3000, chain_init = "SEED+2000+i", chain = "SEED+1000+i"),
  startup_sec = startup_sec, pre_chain_sec = pre_chain_sec, method_sec = method_sec,
  r_version = R.version.string
), file.path(OUT, "mcmc_run_config.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines("OK: all chains produced retained output; see mcmc_diagnostics.csv",
           file.path(OUT, "MCMC_STATUS.txt"))
cat("MCMC done in", round(method_sec), "s; outputs in", OUT, "\n")
