# =============================================================================
# tests/unit/test_metrics.R — unit tests for five-state compartment indexing and seroprevalence summaries.
# =============================================================================

args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests/helper.R"))
source(file.path(root, "scripts/calibration/targets.R"))
targets <- read_targets(file.path(root, "data/calibration_targets.json"))
source(file.path(root, "scripts/calibration/model_metrics.R"))

all_idx <- as.vector(sapply(0:3, function(s) sapply(1:4, function(k)
  sapply(0:4, function(h) sapply(0:5, function(i) idx(s, k, h, i))))))
check_identical(length(all_idx), 480L, "five-state model has 480 compartments")
check_identical(length(unique(all_idx)), 480L, "compartment indices are unique")
check_true(all(all_idx >= 2 & all(all_idx <= 481)), "indices fit output")
check_equal(length(j_pop_idx(1)), 20L, "J population includes r")
check_equal(length(j_sero_idx(1)), 16L, "J seropositive states are a,c,t,r")

out <- matrix(0, nrow = 1, ncol = 481); out[1, 1] <- 150
out[1, idx(1, 1, 0, 3)] <- 50
out[1, idx(1, 1, 1, 3)] <- 10
out[1, idx(1, 1, 4, 3)] <- 40
s <- J_summary_at(out, 1)
check_equal(s$N_hat[4], 100, "J denominator includes u and r")
check_equal(s$p_sero[4], 0.50, "J seroprevalence includes recovered")
check_equal(s$p_viremic[4], 0.10, "viremic prevalence excludes recovered")

# C++ reads treatment rates as stratum-major: all four stages for D, then J,
# F, and X. Keep this explicit because R's default matrix flattening is
# column-major and silently produces the wrong scenario policy.
tau_stratum <- c(10, 20, 30, 40)
tau_stage <- c(1, 2, 3, 4)
tau_by_stratum <- as.vector(t(outer(tau_stratum, tau_stage)))
check_equal(tau_by_stratum,
            c(10, 20, 30, 40, 20, 40, 60, 80,
              30, 60, 90, 120, 40, 80, 120, 160),
            "treatment-rate vector follows C++ stratum-major indexing")

vs <- matrix(0, 1, 481)
vs[1, c(idx(0, 1, 2, 0), idx(3, 4, 3, 5), idx(1, 1, 1, 2))] <- 1  # c, t, a: viraemic
vs[1, c(idx(2, 2, 4, 3), idx(1, 3, 0, 1))] <- 1                    # r, u: not viraemic
check_equal(viraemic_share(vs), 0.6, "viraemic share counts a, c and t states over all PWID")

n_fail <- test_report("tests/unit/test_metrics.R")
if (n_fail > 0) quit(status = 1)
