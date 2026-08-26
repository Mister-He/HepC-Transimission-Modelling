# =============================================================================
# Unit tests: J (prison) summary extraction and fit metrics (model_metrics.R)
# Run: Rscript tests/unit/test_metrics.R
# =============================================================================

.args0 <- commandArgs(FALSE)
.file_arg <- sub("^--file=", "", .args0[grepl("^--file=", .args0)])
root <- normalizePath(file.path(dirname(normalizePath(.file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))
source(file.path(root, "src", "calibration", "targets.R"))
source(file.path(root, "src", "calibration", "model_metrics.R"))

# Compartment indexing: 384 unique columns, inside the output matrix
all_idx <- as.vector(sapply(0:3, function(s) sapply(1:4, function(k)
  sapply(0:3, function(h) sapply(0:5, function(i) idx(s, k, h, i))))))
check_identical(length(all_idx), 384L, "384 compartment columns")
check_identical(length(unique(all_idx)), 384L, "all compartment indices are unique")
check_true(all(all_idx >= 2 & all_idx <= 385),
           "compartment indices fall in the output columns (2..385)")

for (i in 1:6) {
  check_identical(length(j_pop_idx(i)), 16L,
                  sprintf("j_pop_idx(%d) returns 16 J compartments", i))
  check_true(all(j_pop_idx(i) %in% all_idx),
             sprintf("j_pop_idx(%d) indices are valid", i))
}

# Synthetic output: age group 4 (40-49, i = 3): J{NC, acute} = 10,
# J{NC, u} = 90, all other J compartments zero -> N_hat = 100, p = 0.10.
out <- matrix(0, nrow = 1, ncol = 385)
out[1, 1] <- 45
out[1, idx(1, 1, 1, 3)] <- 10   # J acute
out[1, idx(1, 1, 0, 3)] <- 90   # J susceptible/post-SVR

s <- J_summary_at(out, 1)
check_equal(s$N_hat[4], 100, "synthetic J population for 40-49 is 100")
check_equal(s$p_sero[4], 0.10, "synthetic sero prevalence for 40-49 is 0.10")
check_equal(s$p_viremic[4], 0.10, "synthetic viremic prevalence is 0.10")
check_equal(s$N_hat[-4], rep(0, 5), "other age groups have zero population")

# fit_metrics at the observed values: all errors zero
obs <- data.frame(
  age_group = cal_targets$age_groups,
  i = 1:6,
  N_hat = cal_targets$prison_total,
  p_sero = cal_targets$prev_binom,
  p_viremic = cal_targets$prev_binom
)
m <- fit_metrics(obs, target_mode = "sero")
check_equal(m$prevalence_rmse, 0, "RMSE is 0 at the observed prevalence")
check_equal(m$prevalence_max_abs_err, 0, "max prevalence error is 0")
check_equal(m$population_mape, 0, "population MAPE is 0")
check_equal(m$population_max_ape, 0, "population max APE is 0")
check_equal(m$binomial_prevalence_deviance, 0,
            "binomial deviance is 0 at the observed prevalence")
check_equal(m$population_srss, 0,
            "population SRSS is 0 at the observed population")

n_fail <- test_report("tests/unit/test_metrics.R")
if (n_fail > 0) quit(status = 1)
