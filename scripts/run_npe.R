#!/usr/bin/env Rscript
# =============================================================================
# run_npe.R — single-round NPE cross-validation: draws Latin hypercube prior parameter
# sets, simulates the 2015 targets, and trains NPE via npe_train.py. Not used by the
# production path (Nelder-Mead -> reference MCMC -> scenarios).
#
# Steps: data (prior simulations), train (Python), evaluate (posterior-spread gate,
# plus, when --n-sbc > 0, SBC coverage and a rank-histogram skew check). A simulation
# counts only if it is finite and passes check_equilibrium() (the same restriction as
# the reference MCMC target); failures are dropped and counted by reason, never replaced.
# =============================================================================

args       <- commandArgs(trailingOnly = TRUE)
arg        <- function(name, default = NULL) {
  i        <- match(name, args); if (is.na(i)) default else args[i + 1]
}
STEP       <- arg("--step", "all")
ROOT       <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
FIT        <- normalizePath(arg("--fit", file.path(ROOT, "output/calibration/r_20261010_agepop/fit.rds")), mustWork = TRUE)
OUT        <- file.path(ROOT, arg("--out-dir", "output/calibration/npe_lhs"))
PYTHON     <- arg("--python", Sys.which("python3"))
N_SIMS     <- as.integer(arg("--n-sims", "50000"))
N_DRAWS    <- as.integer(arg("--n-draws", "20000")); N_CORES <- as.integer(arg("--n-cores", "4"))
MAX_EPOCHS <- as.integer(arg("--max-epochs", "250"))
SEED       <- as.integer(arg("--seed", "20261012"))          # prior design; training uses --train-seed
TRAIN_SEED <- as.integer(arg("--train-seed", SEED + 1L))
N_SBC      <- as.integer(arg("--n-sbc", "0"))
THREADS    <- as.integer(arg("--threads", "4"))
CONTACT_SD <- as.numeric(arg("--npe-contact-sd", "0.20")); BETA_SD <- as.numeric(arg("--npe-beta-sd", "0.20"))
HIDDEN     <- as.integer(arg("--npe-hidden-features", "50")); TRANSFORMS <- as.integer(arg("--npe-num-transforms", "5"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

source(file.path(ROOT, "scripts/calibration/targets.R"))
targets <- read_targets(file.path(ROOT, "data/calibration_targets.json"))
old <- setwd(file.path(ROOT, "scripts/model")); setup_env <- new.env(parent = globalenv())
sys.source("setup.R", setup_env); setwd(old)
Rcpp::sourceCpp(file.path(ROOT, "scripts/model/sim.cpp"))
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R", "prior_sampling.R"))
  sys.source(file.path(ROOT, "scripts/calibration", f), envir = environment())
fit <- readRDS(FIT); base_params <- setup_env$params; data <- setup_env$data
data$t_start <- 0; data$t_end <- targets$equilibrium_horizon_years
center <- as.numeric(fit$best$theta); N_PARAM <- 12L
prior_sd <- c(rep(CONTACT_SD, 6), rep(BETA_SD, 6))

log_timing <- function(stage, seconds, n_attempted = NA, n_ok = NA, n_nonfinite = NA,
                       n_equilibrium_fail = NA) {
  f <- file.path(OUT, "r_timings.csv")
  row <- data.frame(stage = stage, seconds = seconds, n_attempted = n_attempted, n_ok = n_ok,
                    n_nonfinite = n_nonfinite, n_equilibrium_fail = n_equilibrium_fail,
                    workers = N_CORES, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
  write.table(row, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
}

# Returns c(p_sero, N_hat) or a failure code: "nonfinite" or "equilibrium".
simulate <- function(theta) {
  pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
  if (is.null(pm)) return("nonfinite")
  out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
  if (is.null(out)) return("nonfinite")
  s <- J_summary_final(out)
  if (any(!is.finite(s$p_sero)) || any(!is.finite(s$N_hat)) || any(s$N_hat <= 0)) return("nonfinite")
  eq <- tryCatch(check_equilibrium(out), error = function(e) NULL)
  if (is.null(eq) || !eq$pass) return("equilibrium")
  c(s$p_sero, s$N_hat)
}

simulate_many <- function(theta, stage) {
  t0 <- proc.time()[["elapsed"]]
  one <- function(i) simulate(theta[i, ])
  ans <- if (N_CORES > 1 && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(nrow(theta)), one, mc.cores = N_CORES) else
    lapply(seq_len(nrow(theta)), one)
  ok <- vapply(ans, is.numeric, logical(1))
  code <- vapply(ans, function(a) if (is.numeric(a)) "ok" else as.character(a)[1], character(1))
  log_timing(stage, proc.time()[["elapsed"]] - t0, nrow(theta), sum(ok),
             sum(code == "nonfinite"), sum(code == "equilibrium"))
  list(x = do.call(rbind, ans[ok]), ok = ok, code = code)
}

write_pairs <- function(theta, sim, theta_file, x_file) {
  write.csv(setNames(as.data.frame(theta[sim$ok, , drop = FALSE]), paste0("theta", 1:12)),
            file.path(OUT, theta_file), row.names = FALSE)
  write.csv(setNames(as.data.frame(sim$x), c(paste0("p", 1:6), paste0("N", 1:6))),
            file.path(OUT, x_file), row.names = FALSE)
}

step_data <- function() {
  old_gates <- file.path(OUT, c("NPE_OK.txt", "NPE_NOT_OK.txt"))
  file.remove(old_gates[file.exists(old_gates)])
  th <- sample_prior_lhs(N_SIMS, center, prior_sd, SEED)
  sim <- simulate_many(th, "prior_simulation")
  write_pairs(th, sim, "theta_train.csv", "x_train.csv")
  # Every attempted parameter set with its outcome, so exclusions can be audited.
  write.csv(data.frame(setNames(as.data.frame(th), paste0("theta", 1:12)), outcome = sim$code),
            file.path(OUT, "prior_design.csv"), row.names = FALSE)
  if (N_SBC > 0) {
    th_sbc <- sample_prior_random(N_SBC, center, prior_sd, SEED + 1L)
    sbc <- simulate_many(th_sbc, "sbc_simulation")
    write_pairs(th_sbc, sbc, "theta_sbc.csv", "x_sbc.csv")
  }
  x <- c(targets$p_prev_binom, targets$population)
  write.csv(setNames(as.data.frame(t(x)), c(paste0("p", 1:6), paste0("N", 1:6))),
            file.path(OUT, "x_obs.csv"), row.names = FALSE)
  jsonlite::write_json(list(
    sampling = "Latin hypercube on the unit cube, mapped by the Normal inverse CDF per parameter",
    prior = sprintf("independent Normal(fit$best$theta, sd = %s/%s) on log-scale theta", CONTACT_SD, BETA_SD),
    prior_centre = center, design_seed = SEED, n_attempted = N_SIMS, n_valid = sum(sim$ok),
    n_excluded_nonfinite = sum(sim$code == "nonfinite"), n_excluded_equilibrium = sum(sim$code == "equilibrium"),
    exclusion_rule = "drop simulations that are non-finite or fail check_equilibrium (10-year lag); no replacement",
    fit = FIT), file.path(OUT, "simulation_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
  cat("valid prior simulations:", sum(sim$ok), "/", N_SIMS,
      "| non-finite:", sum(sim$code == "nonfinite"), "| equilibrium failures:", sum(sim$code == "equilibrium"), "\n")
}

step_train <- function() {
  py <- file.path(ROOT, "scripts/npe_train.py")
  a <- c(py, "--theta", file.path(OUT, "theta_train.csv"), "--x", file.path(OUT, "x_train.csv"),
         "--x-obs", file.path(OUT, "x_obs.csv"), "--out", OUT, "--seed", TRAIN_SEED, "--n-draws", N_DRAWS,
         "--contact-anchor", paste(exp(center[1:6]), collapse = ","),
         "--beta-anchor", paste(exp(center[7:12]), collapse = ","), "--contact-sd", CONTACT_SD,
         "--beta-sd", BETA_SD, "--hidden-features", HIDDEN, "--num-transforms", TRANSFORMS,
         "--population-sigma-log", targets$population_sigma_log, "--max-epochs", MAX_EPOCHS,
         "--threads", THREADS, "--sampling-label", "lhs")
  if (N_SBC > 0) a <- c(a, "--n-sbc", N_SBC, "--theta-sbc", file.path(OUT, "theta_sbc.csv"),
                        "--x-sbc", file.path(OUT, "x_sbc.csv"))
  t0 <- proc.time()[["elapsed"]]
  status <- system2(PYTHON, a, env = c("MPLCONFIGDIR=/tmp/mplcache", paste0("OMP_NUM_THREADS=", THREADS)))
  log_timing("python_training", proc.time()[["elapsed"]] - t0)
  if (status != 0) stop("NPE training failed")
}

step_evaluate <- function() {
  post <- read.csv(file.path(OUT, "posterior_samples_npe.csv"))
  sd <- apply(post[paste0("theta", 1:12)], 2, stats::sd)
  q <- data.frame(parameter = paste0("theta", 1:12), posterior_sd = sd, posterior_broad_ok = sd >= 0.05,
                  coverage_95 = NA_real_, coverage_90 = NA_real_, coverage_ok = NA,
                  rank_shape = NA_character_, rank_shape_ok = NA)
  sbc_file <- file.path(OUT, "sbc_summary.csv")
  if (file.exists(sbc_file)) {
    sbc <- read.csv(sbc_file)
    q$coverage_95 <- sbc$coverage_95; q$coverage_90 <- sbc$coverage_90
    q$coverage_ok <- sbc$coverage_95 >= 0.85 & sbc$coverage_95 <= 1 & sbc$coverage_90 >= 0.75 & sbc$coverage_90 <= 1
    q$rank_shape <- sbc$rank_shape; q$rank_shape_ok <- sbc$rank_shape == "uniform"
  }
  write.csv(q, file.path(OUT, "npe_evaluation.csv"), row.names = FALSE)
  ok <- all(q$posterior_broad_ok) && all(is.na(q$coverage_ok) | q$coverage_ok) &&
    all(is.na(q$rank_shape_ok) | q$rank_shape_ok) && all(is.finite(as.matrix(post)))
  note <- if (file.exists(sbc_file)) "" else
    " (SBC not run; accuracy is checked against the reference MCMC in output/calibration/npe_vs_mcmc/)"
  writeLines(paste0(if (ok) "PASS: NPE posterior accepted." else
               "FAIL: NPE quality gate failed; inspect diagnostics before using posterior.", note),
             file.path(OUT, if (ok) "NPE_OK.txt" else "NPE_NOT_OK.txt"))
  stale <- file.path(OUT, if (ok) "NPE_NOT_OK.txt" else "NPE_OK.txt")
  if (file.exists(stale)) file.remove(stale)
  cat("NPE quality gate:", if (ok) "PASS" else "FAIL", "\n")
  invisible(ok)
}

if (STEP %in% c("data", "all")) step_data()
if (STEP %in% c("train", "all")) step_train()
if (STEP %in% c("evaluate", "all")) step_evaluate()
