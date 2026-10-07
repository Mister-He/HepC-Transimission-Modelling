#!/usr/bin/env Rscript
# Deterministic calibration. Every objective evaluation runs model years 0--150;
# the final row is the equilibrium-checked 2015 snapshot.

# Read arguments from CommandLine
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(name, args); if (is.na(i)) default else args[i + 1]
}
ROOT       <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
RUN_ID     <- arg("--run-id", format(Sys.time(), "%Y%m%dT%H%M%S"))
OUT        <- file.path(ROOT, arg("--out-dir", "output/calibration"), RUN_ID)
SEED       <- as.integer(arg("--seed", "20260929"))
MAXIT      <- as.integer(arg("--maxit", "1000"))
N_STARTS   <- as.integer(arg("--n-starts", "8"))
PERTURB_SD <- as.numeric(arg("--sd-perturb", "0.7"))

# Read from repo path
source(file.path(ROOT, "scripts", "calibration", "targets.R"))
targets    <- read_targets(file.path(ROOT, "data", "calibration_targets.json"))
old        <- setwd(file.path(ROOT, "scripts", "model"))
setup_env  <- new.env(parent = globalenv()); sys.source("setup.R", setup_env); setwd(old)

# Read simulation logic from cpp
Rcpp::sourceCpp(file.path(ROOT, "scripts", "model", "sim.cpp"))
base_params  <- setup_env$params
data         <- setup_env$data
data$t_start <- 0
data$t_end   <- targets$equilibrium_horizon_years
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R", "calibrate_nm.R"))
  sys.source(file.path(ROOT, "scripts", "calibration", f), envir = environment())

# Save calibration results insdie output/
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(ROOT, "data", "calibration_targets.json"),
          file.path(OUT, "targets.json"), overwrite = TRUE)

starts <- make_start_set(max(0, N_STARTS - 1), PERTURB_SD, SEED)
starts[[length(starts) + 1L]] <- list(
  id = "population_informed", theta0 = clip_to_bounds(informed_theta0(base_params, targets)),
  seed = SEED + 2000
)
objective <- make_objective(base_params, data)

cat("== calibration", RUN_ID, "==\n")
cat("snapshot: model year", targets$equilibrium_horizon_years,
    "| maxit:", MAXIT, "| starts:", length(starts), "\n")
runs <- lapply(starts, function(st) {
  cat("start", st$id, "\n")
  run_nm_start(st$theta0, objective, st$id, MAXIT, seed = st$seed)
})
best <- runs[[which.min(vapply(runs, `[[`, numeric(1), "objective"))]]
pm_best <- build_params(best$theta, base_params)
out_best <- run_sim(pm_best, data)
summary_best <- J_summary_final(out_best)
eq_best <- check_equilibrium(out_best)

summarise_run <- function(r) {
  pm <- build_params(r$theta, base_params); out <- run_sim(pm, data)
  s <- J_summary_final(out); eq <- check_equilibrium(out); m <- fit_metrics(s)
  data.frame(start_id = r$start_id, seed = r$seed, objective = r$objective,
             convergence = r$convergence, fn_evaluations = r$fn_evaluations,
             equilibrium_pass = eq$pass, equilibrium_violation = eq$violation,
             prevalence_rmse = m$prevalence_rmse, population_mape = m$population_mape,
             setNames(as.list(r$theta), paste0("theta", 1:N_THETA)),
             setNames(as.list(exp(r$theta[1:6])), paste0("contact_scale", 1:6)),
             setNames(as.list(exp(r$theta[7:12])), paste0("beta_scale", 1:6)))
}
solutions <- do.call(rbind, lapply(runs, summarise_run))
write.csv(solutions, file.path(OUT, "solutions.csv"), row.names = FALSE)

predictions <- data.frame(
  age_group = summary_best$age_group,
  population_target = targets$population,
  population_hat = summary_best$N_hat,
  seroprevalence_target = targets$p_prev_binom,
  seroprevalence_hat = summary_best$p_sero,
  viremic_prevalence_hat = summary_best$p_viremic
)
write.csv(predictions, file.path(OUT, "predictions.csv"), row.names = FALSE)
write.csv(eq_best$by_age, file.path(OUT, "equilibrium.csv"), row.names = FALSE)

history <- do.call(rbind, lapply(runs, function(r) {
  if (is.null(r$history)) return(NULL)
  cbind(data.frame(start_id = r$start_id), r$history)
}))
write.csv(history, file.path(OUT, "optimization_history.csv"), row.names = FALSE)
write.csv(data.frame(start_id = vapply(starts, `[[`, character(1), "id"),
                     seed = vapply(starts, `[[`, numeric(1), "seed")),
          file.path(OUT, "initial_values.csv"), row.names = FALSE)

config <- list(
  run_id = RUN_ID, seed = SEED, maxit = MAXIT, n_starts = length(starts),
  horizon = targets$equilibrium_horizon_years, target_year = targets$target_year,
  fitted_parameters = 12, optimizer = "Nelder-Mead",
  target_source = "data/calibration_targets.json",
  inflow_source = "data/singstat_new_drug_abusers.csv",
  inflow_rule = "2015 baseline; 2016--2025 SingStat age-specific counts; latest year held constant through 2030; all multiplied by calibrated beta scale",
  best_objective = best$objective, equilibrium_pass = eq_best$pass,
  git_sha = system("git rev-parse HEAD", intern = TRUE),
  setup_md5 = unname(tools::md5sum(file.path(ROOT, "scripts/model/setup.R"))),
  sim_cpp_md5 = unname(tools::md5sum(file.path(ROOT, "scripts/model/sim.cpp")))
)
jsonlite::write_json(config, file.path(OUT, "run_config.json"), pretty = TRUE,
                     auto_unbox = TRUE)

fit <- list(run_id = RUN_ID, best = best, runs = runs, solutions = solutions,
            predictions = predictions, equilibrium = eq_best, targets = targets,
            base_params = base_params, data = data)
saveRDS(fit, file.path(OUT, "fit.rds"))
cat("best objective:", best$objective, "| equilibrium:", eq_best$pass, "\n")
