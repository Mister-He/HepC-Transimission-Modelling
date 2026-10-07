#!/usr/bin/env Rscript
# Adaptive-Metropolis validation initialized from the NPE posterior covariance.
# The target likelihood is the same Binomial prevalence + log-normal population
# likelihood used by deterministic calibration and NPE observation augmentation.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) default else args[i + 1]
}
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
FIT <- normalizePath(arg("--fit", file.path(ROOT,
  "output/calibration/r_20260929_agepop/fit.rds")), mustWork = TRUE)
NPE_DIR <- normalizePath(arg("--npe-dir", file.path(ROOT,
  "output/calibration/npe_agepop_recal")), mustWork = TRUE)
OUT <- file.path(ROOT, arg("--out-dir", "output/calibration/mcmc_recal"))
N_ITER <- as.integer(arg("--n-iter", "4000"))
BURNIN <- as.integer(arg("--burnin", "1000"))
THIN <- as.integer(arg("--thin", "5"))
N_CHAINS <- as.integer(arg("--n-chains", "3"))
PARALLEL_CHAINS <- as.integer(arg("--parallel-chains", "1"))
SEED <- as.integer(arg("--seed", "20260930"))
CONTACT_SD <- as.numeric(arg("--contact-sd", "0.35"))
BETA_SD <- as.numeric(arg("--beta-sd", "0.35"))
ESS_MIN <- as.numeric(arg("--ess-min", "4000"))
RHAT_MAX <- as.numeric(arg("--rhat-max", "1.001"))
DT <- as.numeric(arg("--dt", as.character(1 / 365)))
SAME_START <- identical(arg("--same-start", "false"), "true")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
gate_file <- file.path(OUT, "MCMC_STRICT_GATE.txt")
if (file.exists(gate_file)) file.remove(gate_file)

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

files <- c(file.path(NPE_DIR, "posterior_samples_npe.csv"),
           file.path(NPE_DIR, "posterior_samples_npe_seed2.csv"))
files <- files[file.exists(files)]
if (!length(files)) stop("No NPE posterior samples found")
npe <- do.call(rbind, lapply(files, function(f) read.csv(f)[, paste0("theta", 1:12)]))
npe <- as.matrix(npe)
if (any(!is.finite(npe))) stop("NPE posterior contains non-finite values")

# Match the NPE regularized prior: best-fit-centred log-normal scales.
priors <- make_priors(
  contact_anchor = exp(fit$best$theta[1:6]),
  beta_anchor = exp(fit$best$theta[7:12]),
  contact_sd = CONTACT_SD, beta_sd = BETA_SD, beta_df = Inf
)
base_params <- setup_env$params
data <- setup_env$data
data$t_start <- 0
data$t_end <- targets$equilibrium_horizon_years
data$dt <- DT
logpost <- function(theta)
  log_posterior(theta, base_params, data, priors, target_mode = "sero",
                equilibrium_penalty = TRUE)

Sigma0 <- stats::cov(npe)
ridge <- max(1e-6, median(diag(Sigma0)) * 1e-4)
Sigma0 <- Sigma0 + diag(ridge, ncol(Sigma0))
theta_mean <- colMeans(npe)
set.seed(SEED)
theta0 <- if (SAME_START) {
  matrix(theta_mean, nrow = N_CHAINS, ncol = length(theta_mean), byrow = TRUE)
} else rbind(theta_mean,
  MASS::mvrnorm(N_CHAINS - 1L, theta_mean, Sigma0))
rownames(theta0) <- paste0("chain", seq_len(N_CHAINS))

if (!all(vapply(seq_len(N_CHAINS), function(i) is.finite(logpost(theta0[i, ])), logical(1))))
  stop("At least one NPE-initialized MCMC chain has non-finite log posterior")

cat("Running", N_CHAINS, "adaptive MCMC chains; NPE draws:", nrow(npe),
    "iterations per chain:", N_ITER, "\n")
run_chain <- function(i) {
  am_mcmc(theta0[i, ], logpost, n_iter = N_ITER, burnin = BURNIN,
          Sigma0 = Sigma0, seed = SEED + i, adapt_from = 200, thin = THIN)
}
chains <- if (PARALLEL_CHAINS > 1 && .Platform$OS.type != "windows") {
  parallel::mclapply(seq_len(N_CHAINS), run_chain,
                     mc.cores = min(PARALLEL_CHAINS, N_CHAINS))
} else lapply(seq_len(N_CHAINS), run_chain)

chain_rows <- do.call(rbind, lapply(seq_along(chains), function(i) {
  x <- as.data.frame(chains[[i]]$samples)
  names(x) <- paste0("theta", 1:12)
  x$chain <- paste0("chain", i)
  x$iteration <- chains[[i]]$iteration
  x
}))
write.csv(chain_rows, file.path(OUT, "posterior_samples_mcmc.csv"), row.names = FALSE)
write.csv(posterior_summaries(as.matrix(chain_rows[, paste0("theta", 1:12)])),
          file.path(OUT, "mcmc_parameter_summary.csv"), row.names = FALSE)

if (requireNamespace("coda", quietly = TRUE)) {
  mlist <- coda::mcmc.list(lapply(chains, function(x) coda::mcmc(x$samples)))
  gd <- coda::gelman.diag(mlist, autoburnin = FALSE, multivariate = TRUE)
  ess <- coda::effectiveSize(mlist)
  diagnostics <- data.frame(
    parameter = paste0("theta", 1:12),
    rhat = gd$psrf[, "Point est."],
    ess_pooled = as.numeric(ess),
    acceptance_mean = mean(vapply(chains, `[[`, numeric(1), "acceptance")),
    n_iter = N_ITER, burnin = BURNIN, thin = THIN
  )
  write.csv(diagnostics, file.path(OUT, "mcmc_diagnostics.csv"), row.names = FALSE)
} else {
  diagnostics <- data.frame(parameter = paste0("theta", 1:12),
    rhat = NA_real_, ess_pooled = NA_real_,
    acceptance_mean = mean(vapply(chains, `[[`, numeric(1), "acceptance")),
    n_iter = N_ITER, burnin = BURNIN, thin = THIN)
  write.csv(diagnostics, file.path(OUT, "mcmc_diagnostics.csv"), row.names = FALSE)
}

strict_pass <- all(is.finite(diagnostics$ess_pooled) &
  diagnostics$ess_pooled > ESS_MIN &
  is.finite(diagnostics$rhat) & diagnostics$rhat < RHAT_MAX)
writeLines(if (strict_pass) {
  sprintf("PASS: ESS > %.0f and Rhat < %.4f for every parameter.",
          ESS_MIN, RHAT_MAX)
} else {
  sprintf("FAIL: strict gate requires ESS > %.0f and Rhat < %.4f for every parameter.",
          ESS_MIN, RHAT_MAX)
}, file.path(OUT, "MCMC_STRICT_GATE.txt"))

write.csv(Sigma0, file.path(OUT, "npe_covariance_used.csv"), row.names = FALSE)
jsonlite::write_json(list(
  fit = FIT, npe_dir = NPE_DIR, n_npe_draws = nrow(npe),
  n_chains = N_CHAINS, n_iter = N_ITER, burnin = BURNIN, thin = THIN,
  parallel_chains = PARALLEL_CHAINS,
  same_start = SAME_START,
  prior = sprintf("best-fit-centred Normal(%s) on all log scales; same Binomial/log-normal target likelihood",
                  CONTACT_SD),
  beta_sd = BETA_SD, ess_min = ESS_MIN, rhat_max = RHAT_MAX,
  integration_dt = DT,
  strict_pass = strict_pass,
  covariance = "empirical covariance of combined NPE seed1 + seed2 posterior"
), file.path(OUT, "mcmc_run_config.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Wrote MCMC validation outputs to", OUT, "\n")
