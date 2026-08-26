# =============================================================================
# probe_beta_constraint.R — beta_scale > 1 soft-rule probe
#
# Reparameterises the six beta inflow scales so that beta_scale = 1 + exp(phi)
# (i.e. beta_scale in (1, Inf)), starts from the best fit of a reference run,
# and re-optimises with Nelder-Mead. Compares acceptance metrics against the
# reference. If the constrained fit fails the acceptance criteria, the soft
# rule is abandoned for the affected age groups and this probe documents why.
#
# Usage:
#   Rscript scripts/probe_beta_constraint.R \
#     --root . --fit output/calibration/<run_id>/fit.rds \
#     --out-dir output/calibration/probe_beta_constraint_<run_id> \
#     [--maxit 2000]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(name, default = NULL) {
  pos <- match(name, args)
  if (is.na(pos)) return(default)
  args[pos + 1]
}

ROOT     <- normalizePath(arg_val("--root", getwd()), mustWork = TRUE)
FIT_PATH <- normalizePath(arg_val("--fit"), mustWork = TRUE)
OUT_DIR  <- normalizePath(arg_val("--out-dir",
              file.path(ROOT, "output", "calibration", "probe_beta_constraint")),
            mustWork = FALSE)
MAXIT    <- as.integer(arg_val("--maxit", "2000"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

old_wd <- setwd(file.path(ROOT, "src"))
setup_env <- new.env(parent = globalenv())
sys.source("setup.R", envir = setup_env)
setwd(old_wd)

if (!exists("run_sim", envir = globalenv(), inherits = FALSE)) {
  Rcpp::sourceCpp(file.path(ROOT, "src", "sim.cpp"))
}
run_sim <- get("run_sim", envir = globalenv(), inherits = FALSE)
base_params <- setup_env$params
data_base   <- setup_env$data

fit <- readRDS(FIT_PATH)
theta0 <- fit$best$theta
data_local <- data_base
data_local$t_start <- fit$t_start
data_local$t_end   <- fit$t_end

for (f in c("targets.R", "model_metrics.R", "equilibrium.R",
            "likelihood.R", "calibrate_nm.R")) {
  sys.source(file.path(ROOT, "src", "calibration", f), envir = environment())
}
TARGET_TIME <- 45

# --- constrained parameterisation: beta_scale = 1 + exp(phi) ---------------
constrained_objective <- function(base_params, data, target_mode = TARGET_MODE) {
  force(base_params); force(data); force(target_mode)
  function(phi) {
    if (length(phi) != N_THETA) return(1e12)
    theta <- phi
    theta[7:12] <- log1p(exp(phi[7:12]))   # beta_scale = 1 + exp(phi)
    pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
    if (is.null(pm)) return(1e12)
    out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
    if (is.null(out)) return(1e12)
    s <- J_summary_final(out)
    p_hat <- if (target_mode == "sero") s$p_sero else s$p_viremic
    if (!all(is.finite(s$N_hat)) || !all(is.finite(p_hat)) || any(s$N_hat <= 0)) return(1e12)
    nll <- nll_prev(p_hat) + nll_pop(s$N_hat)
    if (!is.finite(nll)) return(1e12)
    eq <- tryCatch(check_equilibrium(out, target_mode = target_mode), error = function(e) NULL)
    if (!is.null(eq) && !eq$pass) {
      nll <- nll + EQ_PENALTY_FACTOR * (pmax(eq$max_log_pop_ratio / 0.01, 0) +
                                        pmax(eq$max_prev_change / 0.005, 0))
    }
    nll
  }
}

cat("== beta_scale > 1 probe ==\n")
cat("Reference fit:", FIT_PATH, "\n")

# Map the reference theta to phi space (phi[7:12] = log(beta_scale - 1))
ref_beta_scale <- exp(theta0[7:12])
if (any(ref_beta_scale <= 1)) {
  cat("Reference has beta scales <= 1:",
      paste(round(ref_beta_scale, 3), collapse = ", "), "\n")
  cat("Initialising those at beta_scale = 2 (phi = 0)...\n")
}
phi0 <- theta0
phi0[7:12] <- log(pmax(ref_beta_scale, 2) - 1)

obj_c <- constrained_objective(base_params, data_local)
set.seed(101)
res <- run_nm_start(phi0, obj_c, start_id = "constrained_beta>1",
                    maxit = MAXIT, seed = 101)

theta_c <- res$theta
theta_c[7:12] <- log1p(exp(res$theta[7:12]))
beta_scale_c <- exp(theta_c[7:12])

pm_c <- build_params(theta_c, base_params)
out_c <- run_sim(pm_c, data_local)
s_c <- J_summary_final(out_c)
eq_c <- check_equilibrium(out_c, target_mode = TARGET_MODE)
m_c <- fit_metrics(s_c, target_mode = TARGET_MODE)

# Reference metrics
pm_ref <- build_params(theta0, base_params)
out_ref <- run_sim(pm_ref, data_local)
s_ref <- J_summary_final(out_ref)
eq_ref <- check_equilibrium(out_ref, target_mode = TARGET_MODE)
m_ref <- fit_metrics(s_ref, target_mode = TARGET_MODE)

compare <- rbind(
  data.frame(version = "reference", objective = fit$best$objective,
             equilibrium_pass = eq_ref$pass,
             eq_max_prev = eq_ref$max_prev_change,
             prev_rmse = m_ref$prevalence_rmse,
             prev_max_abs = m_ref$prevalence_max_abs_err,
             pop_mape = m_ref$population_mape,
             pop_max_ape = m_ref$population_max_ape,
             N6 = s_ref$N_hat[6], p6 = s_ref$p_sero[6],
             beta_scale = paste(round(exp(theta0[7:12]), 3), collapse = "|")),
  data.frame(version = "beta>1 constrained", objective = res$objective,
             equilibrium_pass = eq_c$pass,
             eq_max_prev = eq_c$max_prev_change,
             prev_rmse = m_c$prevalence_rmse,
             prev_max_abs = m_c$prevalence_max_abs_err,
             pop_mape = m_c$population_mape,
             pop_max_ape = m_c$population_max_ape,
             N6 = s_c$N_hat[6], p6 = s_c$p_sero[6],
             beta_scale = paste(round(beta_scale_c, 3), collapse = "|"))
)

write.csv(compare, file.path(OUT_DIR, "beta_constraint_compare.csv"), row.names = FALSE)
saveRDS(list(compare = compare, res = res, theta_constrained = theta_c,
             ref_theta = theta0),
        file.path(OUT_DIR, "beta_constraint.rds"))

cat("\nReference:\n")
print(compare[1, ])
cat("\nConstrained (beta_scale > 1):\n")
print(compare[2, ])
cat("\nSaved:", file.path(OUT_DIR, "beta_constraint_compare.csv"), "\n")
