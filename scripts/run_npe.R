#!/usr/bin/env Rscript
# NPE-C calibration. This script uses the same 2015 equilibrium snapshot and
# age-specific observations as run_calibration.R; no Laplace approximation is used.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(name, args); if (is.na(i)) default else args[i + 1]
}
STEP <- arg("--step", "all")
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
FIT <- normalizePath(arg("--fit", file.path(ROOT, "output/calibration/fit.rds")), mustWork = TRUE)
OUT <- file.path(ROOT, arg("--out-dir", "output/calibration/npe"))
PYTHON <- arg("--python", Sys.which("python3"))
N_SIMS <- as.integer(arg("--n-sims", "20000")); N_PROPOSAL <- as.integer(arg("--n-proposal", "6000"))
N_DRAWS <- as.integer(arg("--n-draws", "20000")); N_CORES <- as.integer(arg("--n-cores", "4"))
N_ROUNDS <- as.integer(arg("--n-rounds", "2")); MAX_EPOCHS <- as.integer(arg("--max-epochs", "250"))
SEED <- as.integer(arg("--seed", "20260929"))
CONTACT_SD <- as.numeric(arg("--npe-contact-sd", "0.50")); BETA_SD <- as.numeric(arg("--npe-beta-sd", "0.50"))
HIDDEN <- as.integer(arg("--npe-hidden-features", "64")); TRANSFORMS <- as.integer(arg("--npe-num-transforms", "8"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts/calibration/targets.R"))
targets <- read_targets(file.path(ROOT, "data/calibration_targets.json"))
old <- setwd(file.path(ROOT, "scripts/model")); setup_env <- new.env(parent = globalenv())
sys.source("setup.R", setup_env); setwd(old)
Rcpp::sourceCpp(file.path(ROOT, "scripts/model/sim.cpp"))
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R", "calibrate_nm.R"))
  sys.source(file.path(ROOT, "scripts/calibration", f), envir = environment())
fit <- readRDS(FIT); base_params <- setup_env$params; data <- setup_env$data
data$t_start <- 0; data$t_end <- targets$equilibrium_horizon_years
center <- as.numeric(fit$best$theta); N_PARAM <- 12L; N_SUMMARY <- 12L

sample_prior <- function(n, seed) {
  set.seed(seed); z <- matrix(rnorm(n * N_PARAM), n, N_PARAM)
  z[, 1:6] <- sweep(sweep(z[, 1:6, drop = FALSE], 2, CONTACT_SD, "*"),
                    2, center[1:6], "+")
  z[, 7:12] <- sweep(sweep(z[, 7:12, drop = FALSE], 2, BETA_SD, "*"),
                     2, center[7:12], "+")
  z
}

simulate <- function(theta) {
  pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
  if (is.null(pm)) return(NULL)
  out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  s <- J_summary_final(out)
  if (any(!is.finite(s$p_sero)) || any(!is.finite(s$N_hat)) || any(s$N_hat <= 0)) return(NULL)
  c(s$p_sero, s$N_hat)
}

simulate_many <- function(theta) {
  one <- function(i) simulate(theta[i, ])
  ans <- if (N_CORES > 1 && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(nrow(theta)), one, mc.cores = N_CORES) else
    lapply(seq_len(nrow(theta)), one)
  ok <- !vapply(ans, is.null, logical(1))
  list(x = do.call(rbind, ans[ok]), ok = ok)
}

write_obs <- function() {
  x <- c(targets$p_prev_binom, targets$population)
  write.csv(setNames(as.data.frame(t(x)), c(paste0("p", 1:6), paste0("N", 1:6))),
            file.path(OUT, "x_obs.csv"), row.names = FALSE)
}

step_data <- function() {
  old_gates <- file.path(OUT, c("NPE_OK.txt", "NPE_NOT_OK.txt"))
  file.remove(old_gates[file.exists(old_gates)])
  th <- sample_prior(N_SIMS, SEED); sim <- simulate_many(th)
  write.csv(th[sim$ok, , drop = FALSE], file.path(OUT, "theta_train.csv"), row.names = FALSE)
  write.csv(sim$x, file.path(OUT, "x_train.csv"), row.names = FALSE)
  th_sbc <- sample_prior(min(300, N_SIMS), SEED + 1); sbc <- simulate_many(th_sbc)
  write.csv(th_sbc[sbc$ok, , drop = FALSE], file.path(OUT, "theta_sbc.csv"), row.names = FALSE)
  write.csv(sbc$x, file.path(OUT, "x_sbc.csv"), row.names = FALSE); write_obs()
  cat("finite simulations:", sum(sim$ok), "/", N_SIMS, "\n")
}

py_round <- function(seed, round, suffix = "") {
  tag <- if (nzchar(suffix)) paste0("_", suffix) else ""
  py <- file.path(ROOT, "scripts/npe_train.py")
  a <- c(py, "--theta", file.path(OUT, "theta_train.csv"), "--x", file.path(OUT, "x_train.csv"),
         "--x-obs", file.path(OUT, "x_obs.csv"), "--theta-sbc", file.path(OUT, "theta_sbc.csv"),
         "--x-sbc", file.path(OUT, "x_sbc.csv"), "--out", OUT, "--seed", seed,
         "--round", round, "--n-rounds", N_ROUNDS, "--n-draws", N_DRAWS,
         "--n-proposal", N_PROPOSAL, "--contact-anchor", paste(exp(center[1:6]), collapse = ","),
         "--beta-anchor", paste(exp(center[7:12]), collapse = ","), "--contact-sd", CONTACT_SD,
         "--beta-sd", BETA_SD, "--hidden-features", HIDDEN, "--num-transforms", TRANSFORMS,
         "--population-sigma-log", targets$population_sigma_log,
         "--max-epochs", MAX_EPOCHS)
  if (nzchar(suffix)) a <- c(a, "--suffix", suffix)
  if (round > 1) a <- c(a, "--theta2", file.path(OUT, paste0("proposal_theta_round", round, tag, ".csv")),
                         "--x2", file.path(OUT, paste0("x_round", round, tag, ".csv")))
  status <- system2(PYTHON, a, env = "MPLCONFIGDIR=/tmp/mplcache")
  if (status != 0) stop("NPE failed at round ", round, ", seed ", seed)
}

simulate_proposal <- function(round, suffix = "") {
  tag <- if (nzchar(suffix)) paste0("_", suffix) else ""
  proposal <- read.csv(file.path(OUT, paste0("proposal_theta_round", round, tag, ".csv")))
  sim <- simulate_many(as.matrix(proposal[, paste0("theta", 1:12)]))
  x <- matrix(NA_real_, nrow(proposal), N_SUMMARY)
  x[sim$ok, ] <- sim$x
  x[!sim$ok, ] <- matrix(rep(c(rep(1e-6, 6), rep(1e12, 6)), sum(!sim$ok)),
                          nrow = sum(!sim$ok), byrow = TRUE)
  write.csv(x, file.path(OUT, paste0("x_round", round, tag, ".csv")), row.names = FALSE)
}

step_train <- function() {
  for (s in list(c(SEED, ""), c(SEED + 1, "seed2"))) {
    for (r in seq_len(N_ROUNDS)) {
      py_round(as.integer(s[1]), r, s[2])
      if (r < N_ROUNDS) simulate_proposal(r + 1, s[2])
    }
  }
}

step_evaluate <- function() {
  read_quality <- function(tag = "") {
    sbc <- read.csv(file.path(OUT, paste0("sbc_summary", tag, ".csv")))
    post <- read.csv(file.path(OUT, paste0("posterior_samples_npe", tag, ".csv")))
    sd <- apply(post[paste0("theta", 1:12)], 2, stats::sd)
    within <- sbc$coverage_95 >= 0.85 & sbc$coverage_95 <= 1 &
      sbc$coverage_90 >= 0.75 & sbc$coverage_90 <= 1
    data.frame(parameter = sbc$parameter, coverage_95 = sbc$coverage_95,
               coverage_90 = sbc$coverage_90, posterior_sd = sd,
               coverage_ok = within, posterior_broad_ok = sd >= 0.05)
  }
  q <- read_quality(); write.csv(q, file.path(OUT, "npe_evaluation.csv"), row.names = FALSE)
  ok <- all(q$coverage_ok & q$posterior_broad_ok)
  if (file.exists(file.path(OUT, "sbc_summary_seed2.csv"))) {
    q2 <- read_quality("_seed2"); write.csv(q2, file.path(OUT, "npe_evaluation_seed2.csv"), row.names = FALSE)
    ok <- ok && all(q2$coverage_ok & q2$posterior_broad_ok)
  }
  writeLines(if (ok) "PASS: NPE posterior accepted." else
               "FAIL: NPE quality gate failed; inspect SBC before using posterior.",
             file.path(OUT, if (ok) "NPE_OK.txt" else "NPE_NOT_OK.txt"))
  stale_gate <- file.path(OUT, if (ok) "NPE_NOT_OK.txt" else "NPE_OK.txt")
  if (file.exists(stale_gate)) file.remove(stale_gate)
  cat("NPE quality gate:", if (ok) "PASS" else "FAIL", "\n")
  invisible(ok)
}

if (STEP %in% c("data", "all")) step_data()
if (STEP %in% c("train", "all")) step_train()
if (STEP %in% c("evaluate", "all")) step_evaluate()
