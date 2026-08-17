# =============================================================================
# run_calibration.R — deterministic multi-start Nelder-Mead calibration
#
# Usage:
#   Rscript src/calibration/run_calibration.R \
#     [--root <repo root>] [--run-id <id>] [--seed 101] [--maxit 3000] \
#     [--n-starts 6] [--sd-perturb 0.8] [--t-start -10] [--t-end 55] \
#     [--target-time 45] [--target-mode sero] [--excess-mortality] \
#     [--out-dir output/calibration]
#
# Model-time convention: t = 0 <-> calendar 1970. Simulation runs at least
# 150 model years (t = -10..140 by default); targets compared at t = 45
# (calendar 2015); equilibrium gate at t = 140 vs t = 135.
#
# 12 fitted log-parameters: 6 contact row scales, 6 beta inflow scales.
# Transmission is constant (m_min = m_max = 1, merged into contact scales).
# --excess-mortality adds the 13th parameter eta_s[6] (bounded 1..5) only
# when the 60+ group cannot otherwise be fitted (AGENTS.md contingency).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(name, default = NULL) {
  pos <- match(name, args)
  if (is.na(pos)) return(default)
  args[pos + 1]
}

ROOT       <- normalizePath(arg_val("--root", getwd()), mustWork = TRUE)
RUN_ID     <- arg_val("--run-id", format(Sys.time(), "%Y%m%dT%H%M%S"))
SEED       <- as.integer(arg_val("--seed", "101"))
MAXIT      <- as.integer(arg_val("--maxit", "3000"))
N_STARTS   <- as.integer(arg_val("--n-starts", "6"))
SD_PERTURB <- as.numeric(arg_val("--sd-perturb", "0.8"))
T_START    <- as.numeric(arg_val("--t-start", "-10"))
T_END      <- as.numeric(arg_val("--t-end", "140"))
TARGET_TIME <- as.numeric(arg_val("--target-time", "45"))
TARGET_MODE <- arg_val("--target-mode", "sero")
EXCESS     <- "--excess-mortality" %in% args
OUT_DIR    <- file.path(normalizePath(arg_val("--out-dir",
                     file.path(ROOT, "output", "calibration")), mustWork = FALSE), RUN_ID)
stopifnot(TARGET_MODE %in% c("viremic", "sero"))

cat("== Calibration run:", RUN_ID, "==\n")
cat("Root:", ROOT, "\n")
cat("Excess-mortality parameter (eta_s[6]):", EXCESS, "\n")

# ---------------------------------------------------------------------------
# 1. Baseline model (current tree)
# ---------------------------------------------------------------------------
old_wd <- setwd(file.path(ROOT, "src"))
setup_env <- new.env(parent = globalenv())
sys.source("setup.R", envir = setup_env)
setwd(old_wd)

if (!exists("run_sim", envir = globalenv(), inherits = FALSE)) {
  Rcpp::sourceCpp(file.path(ROOT, "src", "sim.cpp"))
}
run_sim <- get("run_sim", envir = globalenv(), inherits = FALSE)
base_params <- setup_env$params
data_base   <- setup_env$data

data_local <- data_base
data_local$t_start <- T_START
data_local$t_end   <- T_END
cat("Simulation horizon (model years):", T_START, "to", T_END, "\n")
cat("Target time (model year of the 2014-2016 prison screening):", TARGET_TIME, "\n")
cat("Progression rates (p_CC_DC, p_CC_HCC, p_DC_HCC):",
    base_params$p_CC_DC, base_params$p_CC_HCC, base_params$p_DC_HCC, "\n")

# ---------------------------------------------------------------------------
# 2. Calibration modules
# ---------------------------------------------------------------------------
for (f in c("targets.R", "model_metrics.R", "equilibrium.R",
            "likelihood.R", "calibrate_nm.R", "laplace.R", "plot_results.R")) {
  sys.source(file.path(ROOT, "src", "calibration", f), envir = environment())
}
if (EXCESS) {
  FIT_ETA_S6 <- TRUE
  N_THETA <- 13
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package MASS is required for the Laplace approximation (mvrnorm).")
}

# ---------------------------------------------------------------------------
# 3. Metadata
# ---------------------------------------------------------------------------
git_sha <- system("git rev-parse HEAD", intern = TRUE)
setup_sha  <- tools::md5sum(file.path(ROOT, "src", "setup.R"))
simcpp_sha <- tools::md5sum(file.path(ROOT, "src", "sim.cpp"))
pptx_sha   <- tools::md5sum(file.path(ROOT, "Model schematic.pptx"))

run_config <- data.frame(
  field = c(
    "run_id", "git_sha", "setup_R_md5", "sim_cpp_md5", "pptx_md5",
    "R_version", "Rcpp_version", "RcppArmadillo_version", "dplyr_version",
    "MASS_version",
    "target_year", "horizon_start", "horizon_end", "target_time",
    "age_groups", "target_prev_supplied", "target_prison_total",
    "x_prev_rounded", "sigma_pop", "eps_prev", "equilibrium_t_lag",
    "equilibrium_crit_log_pop", "equilibrium_crit_prev",
    "progression_p_NC_CC", "progression_p_CC_DC", "progression_p_CC_HCC",
    "progression_p_DC_HCC", "progression_sources",
    "transmission_multiplier", "excess_mortality_eta_s6_fitted",
    "optimizer", "maxit", "reltol", "n_starts", "perturb_sd", "seed",
    "fitted_parameters",
    "likelihood", "equilibrium_policy", "laplace_n_draws", "target_mode"
  ),
  value = c(
    RUN_ID, git_sha, unname(setup_sha), unname(simcpp_sha), unname(pptx_sha),
    R.version.string,
    as.character(packageVersion("Rcpp")),
    as.character(packageVersion("RcppArmadillo")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("MASS")),
    "2015", as.character(T_START), as.character(T_END), as.character(TARGET_TIME),
    paste(cal_targets$age_groups, collapse = "|"),
    paste(cal_targets$prev_supplied, collapse = "|"),
    paste(cal_targets$prison_total, collapse = "|"),
    paste(cal_targets$x_prev, collapse = "|"),
    as.character(sigma_pop), as.character(eps_prev), "5", "0.01", "0.005",
    as.character(base_params$p_NC_CC),
    as.character(base_params$p_CC_DC),
    as.character(base_params$p_CC_HCC),
    as.character(base_params$p_DC_HCC),
    "Alazawi et al. 2010 (PMID 20497143); Rivera-Irigoin et al. 2006 (PMID 17209765)",
    "constant (m_min = m_max = 1, merged into contact scales)",
    as.character(EXCESS),
    "optim Nelder-Mead", as.character(MAXIT), "1e-8",
    as.character(N_STARTS), as.character(SD_PERTURB), as.character(SEED),
    if (EXCESS) "13: 6 contact scales, 6 beta scales, eta_s[6] (1+4*plogis)"
      else "12: 6 contact scales, 6 beta scales",
    paste0("nll_prev(Binomial ", TARGET_MODE, " prevalence) + nll_pop(log-Normal)"),
    "penalty gate 1e6 * normalized violation",
    "1000",
    TARGET_MODE
  ),
  stringsAsFactors = FALSE
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(run_config, file.path(OUT_DIR, "run_config.csv"), row.names = FALSE)
write.csv(targets_table, file.path(OUT_DIR, "targets.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

# ---------------------------------------------------------------------------
# 4. Fit
# ---------------------------------------------------------------------------
obj_fn <- make_objective(base_params, data_local, target_mode = TARGET_MODE)
starts <- make_start_set(n_perturbed = max(0, N_STARTS - 1),
                         sd_perturb = SD_PERTURB, seed_base = SEED)
a <- plausible_anchors()
cat("  plausible anchors:", paste(round(a$contact, 3), collapse = ","),
    "|", paste(round(a$beta, 3), collapse = ","), "\n")

runs <- lapply(seq_along(starts), function(j) {
  st <- starts[[j]]
  cat("  start", j, "/", length(starts), ":", st$id, "\n")
  run_nm_start(st$theta0, obj_fn, start_id = st$id, maxit = MAXIT, seed = st$seed)
})

best <- runs[[which.min(sapply(runs, `[[`, "objective"))]]
cat("  best objective:", best$objective, "(start:", best$start_id, ")\n")

stability <- run_stability_check(best$theta, obj_fn, seed = SEED + 1000)
write.csv(stability, file.path(OUT_DIR, "stability.csv"), row.names = FALSE)

solution_rows <- lapply(runs, function(r) {
  pm <- build_params(r$theta, base_params)
  out <- run_sim(pm, data_local)
  s <- J_summary_final(out)
  eq <- check_equilibrium(out, target_mode = TARGET_MODE)
  m <- fit_metrics(s, target_mode = TARGET_MODE)
  row <- data.frame(
    start_id = r$start_id, seed = r$seed,
    objective = r$objective, convergence = r$convergence,
    fn_evaluations = r$fn_evaluations, elapsed_sec = r$elapsed_sec,
    setNames(as.list(r$theta), paste0("theta", 1:N_THETA)),
    setNames(as.list(exp(r$theta[1:6])), paste0("contact_scale", 1:6)),
    setNames(as.list(exp(r$theta[7:12])), paste0("beta_scale", 1:6)),
    equilibrium_pass = eq$pass,
    equilibrium_max_log_pop_ratio = eq$max_log_pop_ratio,
    equilibrium_max_prev_change = eq$max_prev_change,
    prevalence_rmse = m$prevalence_rmse,
    prevalence_max_abs_err = m$prevalence_max_abs_err,
    population_mape = m$population_mape,
    population_max_ape = m$population_max_ape,
    binomial_deviance = m$binomial_prevalence_deviance,
    population_srss = m$population_srss,
    nll_prev = m$nll_prev, nll_pop = m$nll_pop
  )
  if (EXCESS) row$eta_s6 <- 1 + 4 * plogis(r$theta[13])
  row
})
solutions <- do.call(rbind, solution_rows)
write.csv(solutions, file.path(OUT_DIR, "solutions.csv"), row.names = FALSE)

initial_rows <- lapply(seq_along(starts), function(j) {
  data.frame(start_id = starts[[j]]$id, seed = starts[[j]]$seed,
             setNames(as.list(starts[[j]]$theta0), paste0("theta0_", 1:N_THETA)))
})
write.csv(do.call(rbind, initial_rows), file.path(OUT_DIR, "initial_values.csv"), row.names = FALSE)

hist_rows <- lapply(seq_along(runs), function(j) {
  h <- runs[[j]]$history
  if (is.null(h)) return(NULL)
  cbind(data.frame(start_id = runs[[j]]$start_id), h)
})
write.csv(do.call(rbind, hist_rows), file.path(OUT_DIR, "optimization_history.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Best-fit predictions, residuals, equilibrium, Laplace intervals
# ---------------------------------------------------------------------------
pm_best <- build_params(best$theta, base_params)
out_best <- run_sim(pm_best, data_local)
s_best <- J_summary_final(out_best)
eq_best <- check_equilibrium(out_best, target_mode = TARGET_MODE)
m_best <- fit_metrics(s_best, target_mode = TARGET_MODE)

predictions <- data.frame(
  age_group = s_best$age_group,
  N_obs = cal_targets$prison_total,
  N_hat = s_best$N_hat,
  p_obs_supplied = cal_targets$prev_supplied,
  p_obs_binom = cal_targets$prev_binom,
  p_hat_sero = s_best$p_sero,
  p_hat_viremic = s_best$p_viremic,
  p_hat_target = if (TARGET_MODE == "sero") s_best$p_sero else s_best$p_viremic
)
residuals <- data.frame(
  age_group = s_best$age_group,
  prev_residual = (if (TARGET_MODE == "sero") s_best$p_sero else s_best$p_viremic) -
    cal_targets$prev_binom,
  log_pop_residual = log(s_best$N_hat) - log(cal_targets$prison_total),
  pop_ape = abs(s_best$N_hat - cal_targets$prison_total) / cal_targets$prison_total
)
eq_table <- eq_best$by_age

write.csv(predictions, file.path(OUT_DIR, "predictions.csv"), row.names = FALSE)
write.csv(residuals, file.path(OUT_DIR, "residuals.csv"), row.names = FALSE)
write.csv(eq_table, file.path(OUT_DIR, "equilibrium.csv"), row.names = FALSE)

cat("Computing Laplace intervals (this re-simulates ~1000 draws)...\n")
lap <- laplace_intervals(best$theta, base_params, data_local,
                         target_mode = TARGET_MODE, n_draws = 1000,
                         rel_thresh = 1e-4, seed = SEED + 500)
if (lap$success) {
  write.csv(lap$intervals, file.path(OUT_DIR, "laplace_intervals.csv"), row.names = FALSE)
  lap_diag <- data.frame(
    field = c("hessian_method", "rank", "effective_dimension", "n_draws_used",
              "n_parameters", "rel_thresh", "n_draws_total"),
    value = c(lap$hessian_method, lap$rank, lap$effective_dimension,
              lap$n_draws_used, N_THETA, "1e-4", lap$n_draws_total)
  )
  write.csv(lap_diag, file.path(OUT_DIR, "laplace_diagnostics.csv"), row.names = FALSE)
  write.csv(data.frame(eigenvalue = lap$eigenvalues),
            file.path(OUT_DIR, "laplace_eigenvalues.csv"), row.names = FALSE)
  write.csv(data.frame(lap$eigenvectors),
            file.path(OUT_DIR, "laplace_eigenvectors.csv"), row.names = TRUE)
} else {
  writeLines(lap$reason, file.path(OUT_DIR, "laplace_FAILED.txt"))
}

# Mortality / progression scenario record
scenario_tbl <- data.frame(
  field = c("name", "baseline_year", "mu", "omega", "smr_source",
            "p_CC_DC", "p_CC_HCC", "p_DC_HCC", "progression_sources",
            "m_min", "m_max", "eta_s"),
  value = c(
    "M1_2015_singapore_mu_plus_scalar_smr",
    "2015",
    paste(base_params$mu, collapse = "|"),
    as.character(base_params$omega),
    "Mathers et al. 2013 pooled PWID SMR (PMID 23554523)",
    as.character(base_params$p_CC_DC),
    as.character(base_params$p_CC_HCC),
    as.character(base_params$p_DC_HCC),
    "Alazawi et al. 2010 (PMID 20497143); Rivera-Irigoin et al. 2006 (PMID 17209765)",
    "1", "1", paste(base_params$eta_s, collapse = "|")
  )
)
write.csv(scenario_tbl, file.path(OUT_DIR, "mortality_scenario.csv"), row.names = FALSE)

fit_obj <- list(
  run_id = RUN_ID, scenario = "M1_2015", best = best, runs = runs,
  stability = stability, solutions = solutions, predictions = predictions,
  residuals = residuals, equilibrium = eq_best, metrics = m_best,
  laplace = lap, targets = cal_targets,
  base_params_md5 = unname(setup_sha), simcpp_md5 = unname(simcpp_sha),
  git_sha = git_sha, t_start = data_local$t_start, t_end = data_local$t_end,
  excess_mortality = EXCESS
)
saveRDS(fit_obj, file.path(OUT_DIR, "fit.rds"))

# ---------------------------------------------------------------------------
# 6. ggplot2 figures (all figures are ggplot2, publication quality)
# ---------------------------------------------------------------------------
plot_all(fit_obj, OUT_DIR)
cat("Figures written to", OUT_DIR, "\n")
cat("Done. Best NLL:", best$objective, "\n")
