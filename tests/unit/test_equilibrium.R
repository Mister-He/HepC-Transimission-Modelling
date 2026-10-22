# =============================================================================
# tests/unit/test_equilibrium.R — unit tests for the equilibrium (stability) gate with five HCV states.
# =============================================================================

args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))
source(file.path(root, "scripts", "calibration", "targets.R"))
targets <- read_targets(file.path(root, "data", "calibration_targets.json"))
source(file.path(root, "scripts", "calibration", "model_metrics.R"))
source(file.path(root, "scripts", "calibration", "equilibrium.R"))

stable <- matrix(1, nrow = 6, ncol = 481); stable[, 1] <- 0:5
eq <- check_equilibrium(stable)
check_true(eq$pass, "unchanged five-state trajectory passes")
check_equal(eq$max_log_pop_ratio, 0, "stable population ratio is zero")
check_equal(eq$max_prev_change, 0, "stable seroprevalence change is zero")

perturbed <- stable; perturbed[6, -1] <- 1.5
eq2 <- check_equilibrium(perturbed)
check_true(!eq2$pass, "a final-state jump fails the gate")
check_true(eq2$max_comp_log_ratio > 0.02, "five-state compartment jump is detected")
check_equal(eq$T_minus_lag, 0, "the 10-year lag reaches back to the first row")

n_fail <- test_report("tests/unit/test_equilibrium.R")
if (n_fail > 0) quit(status = 1)
