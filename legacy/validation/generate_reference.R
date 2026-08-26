# =============================================================================
# generate_reference.R — generate reference outputs from the legacy R/C++ model
#
# Purpose: produce machine-readable reference data used by tests/ to verify
# that the NumPyro/JAX port (src/simulator.py) reproduces the legacy simulator.
#
# Usage (from repo root):
#   Rscript legacy/validation/generate_reference.R <out_dir>
# Requires: R, Rcpp, RcppArmadillo. Writes:
#   <out_dir>/summary_nm.csv          final J-stratum summary at t=45
#   <out_dir>/state_T.csv             full 384-compartment state at t=150
#   <out_dir>/state_T5.csv            full 384-compartment state at t=145
#   <out_dir>/traj_subset.csv         every 1000th row of the full trajectory
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- normalizePath(args[1], mustWork = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

library(Rcpp)
root <- getwd()
setwd(file.path(root, "src"))
env <- new.env(parent = globalenv())
sys.source("setup.R", envir = env)
setwd(root)
for (f in c("targets.R", "model_metrics.R")) {
  sys.source(file.path(root, "src", "calibration", f), envir = env)
}

# Final Nelder-Mead point estimate (legacy/nelder_mead_point_estimate.csv)
theta_nm <- c(
  1.85180416132217, -2.456917410171, -2.83108349195722, -0.582694482047231,
  2.25372719933782, -1.36343925736614, -1.78977034422545, -0.0573350799585561,
  0.206762639702967, 1.65603708644149, 2.82741816698456, 4.70589814858433
)

base <- env$params
pm <- base
for (i in 1:6) pm$C_contact[i, ] <- base$C_contact[i, ] * exp(theta_nm[i])
pm$beta <- base$beta * exp(theta_nm[7:12])

data_local <- env$data
out <- run_sim(pm, data_local)
times <- out[, 1]
row_T  <- nrow(out)
T      <- times[row_T]
row_T5 <- which.min(abs(times - (T - 5)))
row45  <- which.min(abs(times - 45))

summary_df <- env$J_summary_at(out, row45)
summary_df$p_hat_target <- summary_df$p_sero
write.csv(summary_df, file.path(OUT_DIR, "summary_nm.csv"), row.names = FALSE)
write.csv(out[row_T, , drop = FALSE], file.path(OUT_DIR, "state_T.csv"),
          row.names = FALSE)
write.csv(out[row_T5, , drop = FALSE], file.path(OUT_DIR, "state_T5.csv"),
          row.names = FALSE)

idx_sub <- seq(1, row_T, by = 1000)
write.csv(out[idx_sub, , drop = FALSE],
          file.path(OUT_DIR, "traj_subset.csv"), row.names = FALSE)

cat("T =", T, "row_T5 =", row_T5, "row45 =", row45, "\n")
cat("Reference files written to", OUT_DIR, "\n")
