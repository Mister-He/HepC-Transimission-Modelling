# =============================================================================
# Unit tests: equilibrium gate (equilibrium.R)
# Run: Rscript tests/unit/test_equilibrium.R
# =============================================================================

.args0 <- commandArgs(FALSE)
.file_arg <- sub("^--file=", "", .args0[grepl("^--file=", .args0)])
root <- normalizePath(file.path(dirname(normalizePath(.file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))
source(file.path(root, "src", "calibration", "targets.R"))
source(file.path(root, "src", "calibration", "model_metrics.R"))
source(file.path(root, "src", "calibration", "equilibrium.R"))

# Stable trajectory: identical states at T and T-5 -> gate passes.
stable <- matrix(1.0, nrow = 6, ncol = 385)
stable[, 1] <- 0:5
eq <- check_equilibrium(stable)
check_true(eq$pass, "identical end states pass the equilibrium gate")
check_equal(eq$max_log_pop_ratio, 0, "population log-ratio is 0")
check_equal(eq$max_prev_change, 0, "prevalence change is 0")

# Perturbed final row (all compartments x1.5) -> gate fails on state/compartments.
perturbed <- stable
perturbed[6, -1] <- 1.5
eq2 <- check_equilibrium(perturbed)
check_true(!eq2$pass, "a 50% state jump at T fails the gate")
check_true(eq2$max_comp_log_ratio > 0.02,
           "compartment ratio exceeds the 0.02 criterion")

# Default criteria match the documented acceptance rules.
check_equal(check_equilibrium(stable)$T, 5, "T is the final time")
check_equal(check_equilibrium(stable)$T_minus_lag, 0,
            "T-5 is the comparison time")

n_fail <- test_report("tests/unit/test_equilibrium.R")
if (n_fail > 0) quit(status = 1)
