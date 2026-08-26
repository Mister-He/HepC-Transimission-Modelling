# =============================================================================
# Integration tests: compile sim.cpp, simulate, and smoke-test the
# calibration + scenario pipeline end-to-end.
#
# Run: Rscript tests/integration/test_sim_and_pipeline.R
# Requires: Rcpp, RcppArmadillo, dplyr (analysis path)
# =============================================================================

.args0 <- commandArgs(FALSE)
.file_arg <- sub("^--file=", "", .args0[grepl("^--file=", .args0)])
root <- normalizePath(file.path(dirname(normalizePath(.file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))

suppressMessages({
  library(Rcpp)
  library(RcppArmadillo)
})

old_wd <- setwd(file.path(root, "src"))
setup_env <- new.env(parent = globalenv())
sys.source("setup.R", envir = setup_env)
setwd(old_wd)
if (!exists("run_sim", envir = globalenv(), inherits = FALSE)) {
  Rcpp::sourceCpp(file.path(root, "src", "sim.cpp"))
}
run_sim <- get("run_sim", envir = globalenv(), inherits = FALSE)
base_params <- setup_env$params
data_base <- setup_env$data

for (f in c("targets.R", "model_metrics.R", "equilibrium.R",
            "likelihood.R", "calibrate_nm.R")) {
  sys.source(file.path(root, "src", "calibration", f), envir = environment())
}

data_short <- data_base
data_short$t_start <- 0
data_short$t_end <- 1
data_short$dt <- 1 / 52

# ---------------------------------------------------------------------------
# 1. Simulation invariants
# ---------------------------------------------------------------------------
out <- run_sim(base_params, data_short)
check_identical(dim(out)[2], 385L, "output has 385 columns (time + 384 states)")
check_true(all(out[, -1] >= -1e-6), "no negative compartment values")
check_true(all(is.finite(out)), "no NaN/Inf in the output")

tot <- rowSums(out[, -1])
check_true(all(tot > 0) && all(is.finite(tot)), "total population finite > 0")
check_true(max(tot) / min(tot) < 2,
           "total population stays within a factor of 2 over 1 year")

s <- J_summary_at(out, nrow(out))
check_true(all(s$N_hat >= 0), "J age populations are non-negative")
check_true(all(s$p_sero >= 0 & s$p_sero <= 1),
           "J sero prevalence is in [0, 1]")

# Determinism: identical inputs give identical outputs.
out2 <- run_sim(base_params, data_short)
check_identical(out, out2, "simulation is deterministic")

# ---------------------------------------------------------------------------
# 2. Calibration pipeline smoke test (theta = 0, short horizon)
# ---------------------------------------------------------------------------
data_cal <- data_base
data_cal$t_end <- 10
obj <- make_objective(base_params, data_cal)
nll0 <- obj(rep(0, 12))
check_true(is.finite(nll0), "objective is finite at theta = 0")
check_true(nll0 > 0, "objective is positive at theta = 0")

# ---------------------------------------------------------------------------
# 3. Scenario wiring smoke test (regression for the tau_stratum bug)
#    Treating J (prison) should reduce total HCV relative to no treatment.
# ---------------------------------------------------------------------------
scen <- read.csv(file.path(root, "src", "calibration", "scenarios.csv"),
                 stringsAsFactors = FALSE)
check_true(nrow(scen) >= 14, "scenario inventory is populated")

idx_hcv <- as.vector(sapply(0:3, function(s) sapply(1:4, function(k)
  sapply(1:3, function(h) sapply(0:5, function(i)
    s * 96 + (k - 1) * 24 + h * 6 + i + 2L)))))

project_scenario <- function(pm, t_end = 5) {
  data_p <- data_base
  data_p$t_end <- t_end
  run_sim(pm, data_p)
}

# Status quo: no treatment
out0 <- project_scenario(base_params)
hcv0 <- sum(out0[nrow(out0), idx_hcv])

# Prison-only treatment (all stages, all ages)
pm_tx <- base_params
pm_tx$tau <- as.numeric(scen[scen$scenario == "prison_only_treatment",
                            c("tau_NC", "tau_CC", "tau_DC", "tau_HCC")])
pm_tx$tau_stratum <- c(0, 1, 0, 0)   # J only
pm_tx$tau_min_age <- 1L
out_tx <- project_scenario(pm_tx)
hcv_tx <- sum(out_tx[nrow(out_tx), idx_hcv])

check_true(hcv_tx < hcv0, "prison-only treatment reduces total HCV")

# Age-restricted prison treatment must differ from all-age prison treatment.
pm_tx40 <- pm_tx
pm_tx40$tau_min_age <- 4L
out_tx40 <- project_scenario(pm_tx40)
hcv_tx40 <- sum(out_tx40[nrow(out_tx40), idx_hcv])
check_true(hcv_tx < hcv_tx40, "all-age prison treatment beats 40+-only treatment")

# Community-only treatment (D + F) also reduces HCV.
pm_com <- base_params
pm_com$tau <- as.numeric(scen[scen$scenario == "broad_treatment_community",
                              c("tau_NC", "tau_CC", "tau_DC", "tau_HCC")])
pm_com$tau_stratum <- c(1, 0, 1, 0)
pm_com$tau_min_age <- 1L
out_com <- project_scenario(pm_com)
check_true(sum(out_com[nrow(out_com), idx_hcv]) < hcv0,
           "community treatment reduces total HCV")

n_fail <- test_report("tests/integration/test_sim_and_pipeline.R")
if (n_fail > 0) quit(status = 1)
