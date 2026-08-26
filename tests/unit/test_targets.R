# =============================================================================
# Unit tests: calibration target data (targets.R)
# Run: Rscript tests/unit/test_targets.R
# =============================================================================

.args0 <- commandArgs(FALSE)
.file_arg <- sub("^--file=", "", .args0[grepl("^--file=", .args0)])
root <- normalizePath(file.path(dirname(normalizePath(.file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))
source(file.path(root, "src", "calibration", "targets.R"))

check_true(exists("cal_targets"), "cal_targets object is defined")
check_identical(length(cal_targets$age_groups), 6L,
                "six age groups")
check_identical(cal_targets$age_groups,
                c("<20", "20-29", "30-39", "40-49", "50-59", "60+"),
                "age-group labels match the model")

check_true(all(cal_targets$prev_supplied > 0 & cal_targets$prev_supplied < 1),
           "supplied prevalence is strictly inside (0, 1)")
check_true(all(cal_targets$prison_total > 0), "prison totals are positive")
check_true(all(cal_targets$prison_total ==
                 as.integer(cal_targets$prison_total)),
           "prison totals are integers")

# Binomial count reconstruction: x = round(prev * n), so the implied
# Binomial prevalence differs from the supplied decimal by at most 0.5/n.
max_rounding <- max(abs(cal_targets$prev_binom - cal_targets$prev_supplied))
check_true(max_rounding <= 0.5 / min(cal_targets$prison_total),
           "implied Binomial prevalence is within rounding tolerance")

check_identical(cal_targets$x_prev,
                round(cal_targets$prev_supplied * cal_targets$prison_total),
                "x_prev equals round(prev * n)")
check_true(all(cal_targets$x_prev >= 0 & cal_targets$x_prev <= cal_targets$n_prev),
           "x_prev is a valid Binomial success count")

check_true(all(c("age_group", "n_prev", "x_prev",
                 "target_prev_supplied", "target_prev_binom") %in%
                 names(targets_table)),
           "targets_table exposes the documented columns")

n_fail <- test_report("tests/unit/test_targets.R")
if (n_fail > 0) quit(status = 1)
