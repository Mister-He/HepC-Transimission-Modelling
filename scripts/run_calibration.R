#!/usr/bin/env Rscript
# =============================================================================
# run_calibration.R — deterministic multi-start Nelder-Mead calibration of the 13
# log-scale parameters (six contact scales, six entry scales, first-arrest multiplier c); the
# best fit is the posterior mode (seroprevalence, first-arrest and community-PWID likelihood
# + fixed prior, likelihood.R).
#
# Each objective evaluation runs model years 0–150; the final row is the equilibrium-checked
# 2015 snapshot. Starts: uninformative (theta = 0: contact scales 1, entry = CNB counts, c = 1)
# and perturb_1..n
# (theta = 0 plus Normal noise, SD --sd-perturb). --smr-viraemic and
# --smr-nonviraemic override the PWID mortality ratios for sensitivity fits.
# =============================================================================

# Read arguments from CommandLine
args <- commandArgs(trailingOnly = TRUE)
arg  <- function(name, default = NULL) {
  i <- match(name, args); if (is.na(i)) default else args[i + 1]
}
ROOT       <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
RUN_ID     <- arg("--run-id", format(Sys.time(), "%Y%m%dT%H%M%S"))
OUT        <- file.path(ROOT, arg("--out-dir", "output/calibration"), RUN_ID)
SEED       <- as.integer(arg("--seed", "20260929"))
MAXIT      <- as.integer(arg("--maxit", "1000"))
RESTARTS   <- as.integer(arg("--nm-restarts", "0"))     # Nelder-Mead restarts per start
N_STARTS   <- as.integer(arg("--n-starts", "8"))
PERTURB_SD <- as.numeric(arg("--sd-perturb", "0.9"))   # log-scale SD of perturbed starts (was 0.7)
N_CORES    <- as.integer(arg("--cores", "1"))         # starts run in parallel when > 1

# Read from repo path
source(file.path(ROOT, "scripts", "calibration", "targets.R"))
targets    <- read_targets(file.path(ROOT, "data", "calibration_targets.json"))
old        <- setwd(file.path(ROOT, "scripts", "model"))
setup_env  <- new.env(parent = globalenv()); sys.source("setup.R", setup_env); setwd(old)

# Read simulation logic from cpp
Rcpp::sourceCpp(file.path(ROOT, "scripts", "model", "sim.cpp"))
base_params  <- setup_env$params
for (nm in c("smr_viraemic", "smr_nonviraemic")) {  # sensitivity only: override setup.R
  v <- arg(paste0("--", gsub("_", "-", nm)))
  if (!is.null(v)) base_params[[nm]] <- as.numeric(v)
}
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R", "calibrate_nm.R"))
  sys.source(file.path(ROOT, "scripts", "calibration", f), envir = environment())
data         <- equilibrium_data(setup_env$data)

# Save calibration results insdie output/
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(ROOT, "data", "calibration_targets.json"),
          file.path(OUT, "targets.json"), overwrite = TRUE)

starts <- make_start_set(max(0, N_STARTS - 1), PERTURB_SD, SEED)
objective <- make_objective(base_params, data)

cat("== calibration", RUN_ID, "==\n")
cat("snapshot: model year", targets$equilibrium_horizon_years,
    "| maxit:", MAXIT, "| starts:", length(starts), "\n")
# Starts are independent, so they can run in parallel (forked workers; not on Windows).
run_start <- function(st) {
  cat("start", st$id, "\n")
  run_nm_start(st$theta0, objective, st$id, MAXIT, seed = st$seed, restarts = RESTARTS)
}
runs <- if (N_CORES > 1 && .Platform$OS.type != "windows") {
  parallel::mclapply(starts, run_start, mc.cores = min(N_CORES, length(starts)),
                     mc.preschedule = FALSE)
} else lapply(starts, run_start)
failed <- vapply(runs, function(r) inherits(r, "try-error") || is.null(r$objective), logical(1))
if (any(failed)) stop("Nelder-Mead start failed: ", paste(which(failed), collapse = ", "))

# Record the best parammeter set with lowest objective values
best         <- runs[[which.min(vapply(runs, `[[`, numeric(1), "objective"))]]
pm_best      <- build_params(best$theta, base_params)
out_best     <- run_sim(pm_best, data)
summary_best <- J_summary_final(out_best)
cs_best      <- calib_summary(out_best, pm_best)
eq_best      <- check_equilibrium(out_best)

# solutions: Summary of model fitting performance
summarise_run <- function(r) {
  pm <- build_params(r$theta, base_params); out <- run_sim(pm, data)
  cs <- calib_summary(out, pm); eq <- check_equilibrium(out); m <- fit_metrics(cs)
  data.frame(start_id = r$start_id, seed = r$seed, objective = r$objective,
             convergence = r$convergence, fn_evaluations = r$fn_evaluations,
             equilibrium_pass = eq$pass, equilibrium_violation = eq$violation, m,
             setNames(as.list(r$theta), paste0("theta", 1:N_THETA)),
             setNames(as.list(exp(r$theta[IDX_CONTACT])), paste0("contact_scale", 1:6)),
             setNames(as.list(exp(r$theta[IDX_ENTRY])), paste0("entry_scale", 1:6)),
             c_first = exp(r$theta[IDX_C]))
}
solutions <- do.call(rbind, lapply(runs, summarise_run))
write.csv(solutions, file.path(OUT, "solutions.csv"), row.names = FALSE)

# history: Parameter values in each iteration per start points
history <- do.call(rbind, lapply(runs, function(r) {
  					 if (is.null(r$history)) return(NULL)
  					 cbind(data.frame(start_id = r$start_id), r$history)
}))
write.csv(history, file.path(OUT, "optimization_history.csv"), row.names = FALSE)
write.csv(data.frame(start_id = vapply(starts, `[[`, character(1), "id"),
                     seed = vapply(starts, `[[`, numeric(1), "seed")),
          file.path(OUT, "initial_values.csv"), row.names = FALSE)

# predictions: targets against the best fit (and model viraemic prevalence); the community
# PWID target is one number, repeated on every row.
predictions <- data.frame(
  age_group              = summary_best$age_group,
  seroprevalence_target  = targets$seroprevalence,
  seroprevalence_hat     = summary_best$p_sero,
  viremic_prevalence_hat = summary_best$p_viremic,
  first_arrests_target   = targets$first_arrests,
  first_arrests_hat      = cs_best$first_arrests,
  entry_hat              = pm_best$beta,
  community_pwid_target  = targets$community_pwid,
  community_pwid_hat     = cs_best$community_pwid
)
write.csv(predictions, file.path(OUT, "predictions.csv"), row.names = FALSE)
write.csv(eq_best$by_age, file.path(OUT, "equilibrium.csv"), row.names = FALSE)

# config: Config used in current calibration run 
config <- list(
  run_id = RUN_ID, 
  seed = SEED, 
  maxit = MAXIT, 
  nm_restarts = RESTARTS,
  sd_perturb = PERTURB_SD,
  convergence = vapply(runs, `[[`, numeric(1), "convergence"),
  n_starts = length(starts),
  horizon = targets$equilibrium_horizon_years, 
  target_year = targets$target_year,
  fitted_parameters = N_THETA, 
  optimizer = "Nelder-Mead",
  target_source = "data/calibration_targets.json",
  inflow_rule = "entry scale (theta7-12) x CNB 2015 age-specific counts (setup.R beta), held constant in every year",
  best_objective = best$objective, 
  equilibrium_pass = eq_best$pass,
  viraemic_share_2015 = viraemic_share(out_best),  # reported only; it varies over time
  smr_viraemic = base_params$smr_viraemic,
  smr_nonviraemic = base_params$smr_nonviraemic,
  c_first = exp(best$theta[IDX_C]),
  entry_scales = exp(best$theta[IDX_ENTRY]),
  community_pwid = cs_best$community_pwid,
  git_sha = system("git rev-parse HEAD", intern = TRUE),
  setup_md5 = unname(tools::md5sum(file.path(ROOT, "scripts/model/setup.R"))),
  sim_cpp_md5 = unname(tools::md5sum(file.path(ROOT, "scripts/model/sim.cpp")))
)
jsonlite::write_json(config, file.path(OUT, "run_config.json"), pretty = TRUE,
                     auto_unbox = TRUE)

# fit: Overall summary of solutions, predictions and etc in fit.rds
fit <- list(run_id = RUN_ID, 
	    best = within(best, rm(history)), 
	    runs = lapply(runs, function(r) { r$history <- NULL; r }),  # history: optimization_history.csv
	    solutions = solutions,
            predictions = predictions, 
            equilibrium = eq_best, 
            targets = targets,
            base_params = base_params, 
            data = data)
saveRDS(fit, file.path(OUT, "fit.rds"))
cat("best objective:", best$objective, "| equilibrium:", eq_best$pass,
    "| c:", round(exp(best$theta[IDX_C]), 4), "| D+F:", round(cs_best$community_pwid), "| viraemic share 2015:", round(viraemic_share(out_best), 4), "\n")
