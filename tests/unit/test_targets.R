# =============================================================================
# tests/unit/test_targets.R — unit tests for loading the 2015 calibration targets and inflow data.
# =============================================================================

args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests/helper.R"))
source(file.path(root, "scripts/calibration/targets.R"))
targets <- read_targets(file.path(root, "data/calibration_targets.json"))

check_equal(targets$target_year, 2015, "targets describe the 2015 snapshot")
check_equal(length(targets$age_groups), 6, "six age groups are loaded")
check_true(all(abs(targets$seroprevalence - targets$p_prev_binom) <=
                 0.5 / targets$population),
           "age-specific seroprevalence is preserved up to count rounding")
check_true(all(targets$x_prev >= 0 & targets$x_prev <= targets$population),
           "reconstructed binomial counts are valid")
check_true(!"release_sero_2017_target" %in% names(targets),
           "the former 47 percent target is absent")
check_true(file.exists(file.path(root, "data/singstat_new_drug_abusers.csv")),
           "SingStat inflow data are versioned in the repository")

n_fail <- test_report("tests/unit/test_targets.R")
if (n_fail > 0) quit(status = 1)
