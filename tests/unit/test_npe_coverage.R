# Simulation-based coverage checks for the current NPE output.
args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests/helper.R"))
npe_dir <- file.path(root, "output", "calibration", "npe_agepop_recal")
eval_path <- file.path(npe_dir, "npe_evaluation.csv")
post_path <- file.path(npe_dir, "posterior_samples_npe.csv")
if (!file.exists(eval_path) || !file.exists(post_path)) {
  cat("SKIP: current NPE output is not present\n")
  quit(status = 0)
}

eval <- read.csv(eval_path); post <- read.csv(post_path)
theta <- post[, paste0("theta", 1:12), drop = FALSE]
check_identical(nrow(eval), 12L, "NPE evaluation covers all 12 parameters")
check_true(all(eval$coverage_95 >= 0.85 & eval$coverage_95 <= 1),
           "NPE SBC 95% coverage is within the acceptance range")
check_true(all(eval$coverage_90 >= 0.75 & eval$coverage_90 <= 1),
           "NPE SBC 90% coverage is within the acceptance range")
check_true(all(eval$posterior_sd >= 0.05),
           "NPE posterior standard deviations are not collapsed")
check_true(nrow(post) >= 500 && all(is.finite(as.matrix(theta))),
           "NPE posterior draws are finite and sufficiently large")
check_true(file.exists(file.path(npe_dir, "NPE_OK.txt")) &&
             !file.exists(file.path(npe_dir, "NPE_NOT_OK.txt")),
           "active NPE quality gate is accepted without a stale failure marker")

eval2_path <- file.path(npe_dir, "npe_evaluation_seed2.csv")
if (file.exists(eval2_path)) {
  eval2 <- read.csv(eval2_path)
  check_true(all(eval2$coverage_95 >= 0.85 & eval2$coverage_95 <= 1),
             "second-seed SBC 95% coverage is within the acceptance range")
  check_true(all(eval2$coverage_90 >= 0.75 & eval2$coverage_90 <= 1),
             "second-seed SBC 90% coverage is within the acceptance range")
}

n_fail <- test_report("tests/unit/test_npe_coverage.R")
if (n_fail > 0) quit(status = 1)
