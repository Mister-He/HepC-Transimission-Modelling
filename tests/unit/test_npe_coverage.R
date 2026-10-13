# =============================================================================
# tests/unit/test_npe_coverage.R — checks the LHS prior sampler and the current
# single-round NPE output (simulation accounting, posterior spread, quality gate, SBC).
# =============================================================================
args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests/helper.R"))
source(file.path(root, "scripts/calibration/prior_sampling.R"))

# LHS: one point per stratum in every column, and Normal prior marginals on theta.
n <- 2000; mu <- seq(-1, 1, length.out = 12)
u <- lhs_unit(n, 12, 1)
check_true(all(u > 0 & u < 1), "LHS unit design lies inside the unit cube")
check_true(all(apply(u, 2, function(v) all(sort(floor(v * n)) == 0:(n - 1)))),
           "LHS places exactly one point in each stratum of every column")
th <- prior_from_unit(u, mu, 0.2)
check_identical(dim(th), c(as.integer(n), 12L), "LHS prior draws have 12 parameter columns")
check_true(max(abs(colMeans(th) - mu)) < 0.01 && all(abs(apply(th, 2, sd) - 0.2) < 0.01),
           "LHS draws reproduce the Normal(mean, 0.20) prior marginals")
check_true(all(sapply(1:12, function(j) ks.test(th[, j], "pnorm", mu[j], 0.2)$p.value) > 0.01),
           "LHS marginals pass a Kolmogorov-Smirnov test against the prior")

npe_dir <- file.path(root, "output", "calibration", "npe_lhs")
post_path <- file.path(npe_dir, "posterior_samples_npe.csv")
if (file.exists(post_path)) {
  post <- read.csv(post_path); eval <- read.csv(file.path(npe_dir, "npe_evaluation.csv"))
  ss <- jsonlite::fromJSON(file.path(npe_dir, "simulation_summary.json"))
  theta <- as.matrix(post[, paste0("theta", 1:12)])
  check_identical(as.integer(ss$n_attempted), 50000L, "NPE uses 50,000 attempted LHS prior simulations")
  check_identical(as.integer(ss$n_valid + ss$n_excluded_nonfinite + ss$n_excluded_equilibrium +
                             ss$n_excluded_bounds),
                  as.integer(ss$n_attempted), "valid plus excluded simulations equal attempted simulations")
  check_identical(nrow(read.csv(file.path(npe_dir, "theta_train.csv"))), as.integer(ss$n_valid),
                  "training set contains exactly the valid simulations (no replacements)")
  check_true(nrow(post) >= 500 && all(is.finite(theta)), "NPE posterior draws are finite and sufficiently large")
  check_true(all(eval$posterior_sd >= 0.05), "NPE posterior standard deviations are not collapsed")
  check_true(file.exists(file.path(npe_dir, "NPE_OK.txt")) && !file.exists(file.path(npe_dir, "NPE_NOT_OK.txt")),
             "NPE quality gate is accepted without a stale failure marker")
  sbc_path <- file.path(npe_dir, "sbc_summary.csv")
  check_true(file.exists(sbc_path) && file.exists(file.path(npe_dir, "sbc_ranks.npy")),
             "SBC summary and ranks are retained")
  if (file.exists(sbc_path)) {
    sbc <- read.csv(sbc_path)
    check_identical(sbc$parameter, paste0("theta", 1:12), "SBC covers all 12 parameters")
    check_true(all(sbc$n_sbc >= 500), "SBC uses at least 500 synthetic datasets")
    check_true(all(sbc$coverage_95 >= 0.85 & sbc$coverage_90 >= 0.75),
               "SBC coverage passes the 95% and 90% gates for every parameter")
    check_true(all(sbc$rank_shape == "uniform"), "no SBC rank histogram is systematically skewed")
  }
} else cat("SKIP: single-round NPE output is not present\n")

n_fail <- test_report("tests/unit/test_npe_coverage.R")
if (n_fail > 0) quit(status = 1)
