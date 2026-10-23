# =============================================================================
# tests/unit/test_targets.R — unit tests for loading the 2015 calibration targets and the CNB data.
# =============================================================================

args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests/helper.R"))
source(file.path(root, "scripts/calibration/targets.R"))
targets <- read_targets(file.path(root, "data/calibration_targets.json"))

check_equal(targets$target_year, 2015, "targets describe the 2015 snapshot")
check_equal(length(targets$age_groups), 6, "six age groups are loaded")
check_true(all(targets$seroprevalence > 0 & targets$seroprevalence < 1),
           "six seroprevalence targets in (0, 1)")
check_true(targets$beta_precision > 0, "the Beta likelihood has a positive fixed precision")
check_true(is.null(targets$population) && is.null(targets$x_prev),
           "no population or count fields: seroprevalence is the only target")
check_true(!"release_sero_2017_target" %in% names(targets),
           "the former 47 percent target is absent")
check_true(is.null(targets$entry) && is.null(targets$population_sigma_log),
           "no entry or population likelihood settings")
check_true(file.exists(file.path(root, "data/singstat_new_drug_abusers.csv")),
           "SingStat new drug abuser data (Table S2 source) are versioned in the repository")

n_fail <- test_report("tests/unit/test_targets.R")
if (n_fail > 0) quit(status = 1)
