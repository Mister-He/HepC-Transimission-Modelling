#!/usr/bin/env Rscript
# =============================================================================
# run_mcmc_reference.R — reference MCMC calibration: adaptive random-walk Metropolis
# on the exact posterior, with convergence, adaptation and proposal diagnostics.
#
# Target (shared with single-round NPE): the fixed truncated Normal prior of likelihood.R
# (prior_spec()) on the 13 log-scale parameters, the seroprevalence, first-arrest and
# community-PWID likelihood (nll_targets), and the 2015 equilibrium check as a hard
# restriction. Defaults: 8 prior-started chains, 20,000 warmup proposals with adaptation from
# proposal 5,000, 60,000 evaluations per chain, draws stored every 10th iteration (--thin; all
# diagnostics use every retained draw). The stored production run
# (output/calibration/mcmc_reference) used these defaults.
# Starts and the initial proposal covariance come from the prior only.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) default else args[i + 1]
}
T_SCRIPT0 <- proc.time()[["elapsed"]]
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
FIT <- normalizePath(arg("--fit", file.path(ROOT,
  "output/calibration/r_20261027_agepop/fit.rds")), mustWork = TRUE)
OUT <- file.path(ROOT, arg("--out-dir", "output/calibration/mcmc_reference"))
N_CHAINS <- as.integer(arg("--n-chains", "8"))
CORES <- as.integer(arg("--cores", "8"))
N_WARMUP <- as.integer(arg("--n-warmup", "20000"))
MAX_EVALS <- as.integer(arg("--max-evals", "60000"))
ADAPT_START <- as.integer(arg("--adapt-start", "5000"))
ADAPT_EVERY <- as.integer(arg("--adapt-every", "50"))
THIN <- as.integer(arg("--thin", "10"))
SEED <- as.integer(arg("--seed", "20261101"))
INIT_SCALE <- as.numeric(arg("--init-scale", "0.1"))
N_PRIOR_VAR <- as.integer(arg("--n-prior-var", "2048"))
DT <- as.numeric(arg("--dt", as.character(1 / 365)))
if (ADAPT_START >= N_WARMUP) warning("--adapt-start >= --n-warmup: the proposal never adapts")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

source(file.path(ROOT, "scripts/calibration/targets.R"))
TARGETS <- normalizePath(arg("--targets", file.path(ROOT, "data/calibration_targets.json")), mustWork = TRUE)
targets <- read_targets(TARGETS)  # --targets: synthetic data for the MCMC-based SBC
old <- setwd(file.path(ROOT, "scripts/model"))
setup_env <- new.env(parent = globalenv())
sys.source("setup.R", setup_env)
setwd(old)
Rcpp::sourceCpp(file.path(ROOT, "scripts/model/sim.cpp"))
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R", "mcmc.R", "prior_sampling.R"))
  sys.source(file.path(ROOT, "scripts/calibration", f), envir = environment())
fit <- readRDS(FIT)
base_params <- setup_env$params
data <- setup_env$data
check_fit_current(fit, base_params, data)
data$dt <- DT; data <- equilibrium_data(data)
startup_sec <- proc.time()[["elapsed"]] - T_SCRIPT0
TH <- paste0("theta", 1:N_THETA)

# ---- sampling ----
T_METHOD0 <- proc.time()[["elapsed"]]
prior <- prior_spec()
prior_mean <- prior$mean; prior_sd <- prior$sd
logpost <- function(theta) log_posterior_exact(theta, base_params, data)

# Transformed-prior variances from prior draws (no ODE solves). theta is already
# unconstrained and the prior is defined on theta, so no Jacobian is needed.
prior_draws <- sample_prior_random(N_PRIOR_VAR, prior_mean, prior_sd, SEED + 3000, THETA_LO, THETA_HI)
prior_var <- apply(prior_draws, 2, stats::var)
reg <- diag(1e-6 * prior_var)
Sigma_init <- diag(INIT_SCALE^2 * prior_var)
pre_chain_sec <- proc.time()[["elapsed"]] - T_METHOD0

run_chain <- function(i) {
  draw_prior <- function() prior_from_unit(matrix(runif(N_THETA), 1), prior_mean, prior_sd, THETA_LO, THETA_HI)[1, ]
  res <- arwm_budget(logpost, draw_prior, Sigma_init, reg, n_warmup = N_WARMUP,
                     max_evals = MAX_EVALS, adapt_start = ADAPT_START, adapt_every = ADAPT_EVERY,
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
  acc <- if (!is.null(ch$info)) mean(ch$info[, "accepted"]) else NA
  data.frame(chain = i, status = ch$status, start_source = ch$start_source %||% NA,
             init_tries = ch$init_tries %||% NA, warmup_acceptance = ch$warmup_acceptance %||% NA,
             retained_total = if (!is.null(ch$info)) nrow(ch$info) else 0L, acceptance = acc,
             acceptance_flag = !is.na(acc) && (acc < 0.15 || acc > 0.50),
             mean_jump_retained = if (!is.null(ch$info)) mean(ch$info[, "jump"]) else NA,
             total_evals = ch$total_evals %||% ch$evals, total_ode = ch$total_ode %||% ch$ode,
             warmup_sec = ch$warmup_sec %||% NA, total_sec = ch$total_sec %||% NA,
             cpu_sec = ch$cpu_sec %||% NA, stringsAsFactors = FALSE)
}))
write.csv(chain_summary, file.path(OUT, "mcmc_chain_summary.csv"), row.names = FALSE)
if (any(bad)) {
  writeLines(paste("FAIL: chain(s)", paste(which(bad), collapse = ","), "did not produce retained output"),
             file.path(OUT, "MCMC_STATUS.txt"))
  stop("MCMC chain failure; see mcmc_chain_summary.csv")
}

draws <- do.call(rbind, lapply(chains, function(ch) {
  x <- as.data.frame(ch$samples); names(x) <- TH
  cbind(data.frame(chain = ch$chain, iteration = seq_len(nrow(x))), x, as.data.frame(ch$info))
}))
keep_rows <- (draws$iteration - 1L) %% THIN == 0L
write.csv(draws[keep_rows, ], file.path(OUT, "mcmc_draws.csv"), row.names = FALSE)
warm <- do.call(rbind, lapply(chains, function(ch) {
  x <- as.data.frame(ch$warmup); names(x) <- TH
  cbind(data.frame(chain = ch$chain, warmup_iteration = seq_len(nrow(x)) - 1L), x)
}))
write.csv(warm[(warm$warmup_iteration %% THIN) == 0L, ], file.path(OUT, "mcmc_warmup.csv"), row.names = FALSE)

# ---- convergence: rank-normalised split R-hat, bulk/tail ESS, MCSE (all retained draws) ----
n_min <- min(table(draws$chain))
arr <- array(NA_real_, c(n_min, N_CHAINS, N_THETA))
for (k in seq_len(N_CHAINS)) arr[, k, ] <- as.matrix(draws[draws$chain == k, TH][seq_len(n_min), ])
diagnostics <- do.call(rbind, lapply(1:N_THETA, function(j) {
  x <- arr[, , j]
  data.frame(subset = "all_retained", parameter = TH[j], draws_per_chain = n_min,
             mean = mean(x), sd = sd(as.vector(x)),
             rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x),
             ess_tail = posterior::ess_tail(x), mcse_mean = posterior::mcse_mean(x),
             mcse_q025 = posterior::mcse_quantile(x, 0.025),
             mcse_q975 = posterior::mcse_quantile(x, 0.975),
             iat = n_min * N_CHAINS / posterior::ess_basic(x))
}))
write.csv(diagnostics, file.path(OUT, "mcmc_diagnostics.csv"), row.names = FALSE)
write.csv(posterior_summaries(as.matrix(draws[, TH])), file.path(OUT, "mcmc_parameter_summary.csv"), row.names = FALSE)
per_chain <- do.call(rbind, lapply(seq_len(N_CHAINS), function(k) {
  x <- as.matrix(draws[draws$chain == k, TH])
  data.frame(chain = k, parameter = TH, median = apply(x, 2, median), sd = apply(x, 2, sd),
             q025 = apply(x, 2, quantile, 0.025), q975 = apply(x, 2, quantile, 0.975))
}))
pooled_sd <- setNames(diagnostics$sd, TH); pooled_mean <- setNames(diagnostics$mean, TH)
per_chain$median_shift_sd <- (per_chain$median - pooled_mean[per_chain$parameter]) / pooled_sd[per_chain$parameter]
write.csv(per_chain, file.path(OUT, "mcmc_per_chain.csv"), row.names = FALSE)

# ---- adaptation (diminishing adaptation) and proposal geometry ----
drift <- do.call(rbind, lapply(chains, function(ch) if (is.null(ch$adapt_trace)) NULL else
  data.frame(chain = ch$chain, ch$adapt_trace)))
if (!is.null(drift)) write.csv(drift, file.path(OUT, "covariance_drift.csv"), row.names = FALSE)
eig_row <- function(label, chain, S) {
  e <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  data.frame(matrix_label = label, chain = chain, t(setNames(e, paste0("eigen", 1:N_THETA))),
             condition_number = max(e) / min(e), min_eigen = min(e))
}
eig <- rbind(eig_row("initial_proposal", NA, Sigma_init),
             do.call(rbind, lapply(chains, function(ch) if (is.null(ch$proposal_cov_adapt)) NULL else
               eig_row("adapted_proposal", ch$chain, ch$proposal_cov_adapt))),
             eig_row("retained_draws_covariance", NA, stats::cov(as.matrix(draws[, TH]))))
write.csv(eig, file.path(OUT, "covariance_eigen.csv"), row.names = FALSE)
# Retained-window stability: covariance of consecutive blocks of retained draws (pooled chains).
block <- max(10L, n_min %/% 10L)
blocks <- split(seq_len(n_min), ceiling(seq_len(n_min) / block))
bcov <- lapply(blocks, function(ix) stats::cov(do.call(rbind, lapply(seq_len(N_CHAINS), function(k) matrix(arr[ix, k, ], nrow = length(ix))))))
nb <- length(bcov)
retained_drift <- data.frame(block = seq_len(max(0L, nb - 1L)) + 1L, block_size = rep(block, max(0L, nb - 1L)),
  drift = vapply(seq_len(max(0L, nb - 1L)) + 1L, function(b) norm(bcov[[b]] - bcov[[b - 1]], "F") /
                   norm(bcov[[b - 1]], "F"), numeric(1)))
write.csv(retained_drift, file.path(OUT, "retained_block_covariance_drift.csv"), row.names = FALSE)
jumps <- do.call(rbind, lapply(chains, function(ch) {
  w <- ch$warmup_jump; bin <- ceiling(seq_len(nrow(w)) / 250)
  data.frame(chain = ch$chain, phase = "warmup", bin_end = tapply(seq_len(nrow(w)), bin, max),
             mean_jump = tapply(w[, "jump"], bin, mean), acceptance = tapply(w[, "accepted"], bin, mean))
}))
write.csv(jumps, file.path(OUT, "mcmc_jumps.csv"), row.names = FALSE)

# ---- trace plots (stored draws) ----
if (requireNamespace("ggplot2", quietly = TRUE)) {
  tp <- draws[keep_rows, c("chain", "iteration", TH)]
  tl <- data.frame(chain = factor(rep(tp$chain, N_THETA)), iteration = rep(tp$iteration, N_THETA),
                   parameter = factor(rep(TH, each = nrow(tp)), TH), value = unlist(tp[, TH]))
  p <- ggplot2::ggplot(tl, ggplot2::aes(iteration, value, colour = chain)) +
    ggplot2::geom_line(linewidth = 0.15) + ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 3) +
    ggplot2::labs(x = "Retained iteration", y = "theta", colour = "Chain") + ggplot2::theme_classic(base_size = 8)
  ggplot2::ggsave(file.path(OUT, "mcmc_traceplots.png"), p, width = 11, height = 9, dpi = 150)
}

jsonlite::write_json(list(
  fit = FIT, targets = TARGETS, sampler = "adaptive random-walk Metropolis (arwm_budget, scripts/calibration/mcmc.R)",
  n_chains = N_CHAINS, cores = CORES, n_warmup = N_WARMUP, max_evals_per_chain = MAX_EVALS,
  adapt_start = ADAPT_START, adapt_every = ADAPT_EVERY, safety_p = 0.05, adaptive_scale = "2.38^2/d",
  thin_storage = THIN, init_scale = INIT_SCALE, n_prior_var_draws = N_PRIOR_VAR,
  initial_covariance = "init_scale^2 * diag(prior variance)", starts = "independent prior draws",
  regularisation = "1e-6 * diag(prior variance)",
  prior = "independent Normal(prior_mean, prior_sd) on theta, truncated to the bounds (likelihood.R)",
  prior_mean = prior_mean, prior_sd = prior_sd, theta_lo = THETA_LO, theta_hi = THETA_HI,
  target = "Beta seroprevalence + negative binomial first arrests + log-normal community PWID; equilibrium failure = -Inf",
  integration_dt = DT, seed = SEED,
  seeds = list(prior_variance = SEED + 3000, chain = "SEED+1000+i (initial prior draws and sampling)"),
  startup_sec = startup_sec, pre_chain_sec = pre_chain_sec, method_sec = method_sec,
  cpu_sec_total = sum(chain_summary$cpu_sec), r_version = R.version.string
), file.path(OUT, "mcmc_run_config.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines("OK: all chains produced retained output; see mcmc_diagnostics.csv",
           file.path(OUT, "MCMC_STATUS.txt"))
cat("MCMC done in", round(method_sec), "s; outputs in", OUT, "\n")
