#!/usr/bin/env Rscript
# =============================================================================
# tornado.R — one-way sensitivity of the 2030 results to the fixed inputs. The 13 calibrated
# parameters stay at the Nelder–Mead best fit (no refit); each input is set to its low and then
# its high value, one at a time, and scripts/scenarios_2026.R runs all 17 policies for that
# single parameter set (--param-file). Ranges: the Notation table of
# Hepc_modelling_technical_details.pdf, otherwise the cited literature (REVIEW_REPORT.md R-05).
# Run: Rscript scripts/tornado.R --root . --cores 14
# Output: output/tornado/ (tornado_results.csv, pass_flips.csv, failed_runs.csv, ranges.csv, runs/).
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) { i <- match(name, args); if (is.na(i)) default else args[i + 1] }
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
FIT <- arg("--fit", "output/calibration/r_20261103_agepop/fit.rds")
OUT <- file.path(ROOT, arg("--out-dir", "output/tornado"))
N_CORES <- as.integer(arg("--cores", 8))
dir.create(file.path(OUT, "runs"), recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(file.path(ROOT, FIT)); bp <- fit$base_params
theta_file <- file.path(OUT, "best_theta.csv")
write.csv(setNames(as.data.frame(t(fit$best$theta)), paste0("theta", seq_along(fit$best$theta))),
          theta_file, row.names = FALSE)

# Prison rates: 2.5% and 97.5% posterior quantiles of the Stan fits (all ages move together).
as_draws <- function(f) { d <- readRDS(file.path(ROOT, f))[[2]]; matrix(d, ncol = dim(d)[3], dimnames = list(NULL, dimnames(d)[[3]])) }
rate_q <- function(f, p) { d <- as_draws(f); unname(apply(exp(d[, grep("^log_lambda\\[", colnames(d))]), 2, quantile, p)) }
pi_q <- function(p) unname(quantile(plogis(as_draws("output/rearrest_rate_fit.rds")[, "pi_raw"]), p))

# SVR composites in setup.R are rho * GT3 + (1 - rho) * other; one genotype's SVR moves the
# composite by its weight (RBV-free DC is what the scenarios use for DC; RBV+ is used for HCC).
gt <- list(a3NC = 0.976, aoNC = 0.988, a3DCm = 0.45, aoDCm = 0.90, a3DCp = 0.87, aoDCp = 0.92)
svr <- function(name, v) {
  w <- if (startsWith(name, "a3")) bp$rho else 1 - bp$rho; d <- w * (v - gt[[name]])
  switch(substr(name, 3, 5), NC = list(alpha_NC = bp$alpha_NC + d), DCm = list(alpha_DC_neg = bp$alpha_DC_neg + d),
         DCp = list(alpha_DC_pos = bp$alpha_DC_pos + d, alpha_HCC = bp$alpha_HCC + d))
}
rho_set <- function(v) list(rho = v,
  alpha_NC = bp$alpha_NC + (v - bp$rho) * (gt$a3NC - gt$aoNC),
  alpha_DC_neg = bp$alpha_DC_neg + (v - bp$rho) * (gt$a3DCm - gt$aoDCm),
  alpha_DC_pos = bp$alpha_DC_pos + (v - bp$rho) * (gt$a3DCp - gt$aoDCp),
  alpha_HCC = bp$alpha_HCC + (v - bp$rho) * (gt$a3DCp - gt$aoDCp))
hr0 <- bp$smr_viraemic / bp$smr_nonviraemic  # extrahepatic HR (Lee 2012), 1.47

# name, label, low, high, source, function(value) -> named list of overrides
inputs <- list(
  list("q", "Transmission probability q", 0.0032, 0.0105, "PDF Table 1 (Boelen 2014)", function(v) list(q = v)),
  list("kappa", "Spontaneous clearance", 0.22, 0.29, "PDF Table 1 (Micallef 2006)", function(v) list(kappa = v)),
  list("iota1", "Acute duration (years)", 0.25, 1, "Grebely 2014 (3-12 months)", function(v) list(iota1 = v)),
  list("rho", "GT3 share", 0.74, 0.82, "PDF Table 1 (Zhao 2019)", rho_set),
  list("p_NC_CC", "NC to CC", 0.027 - 1.96 * 0.0008, 0.027 + 1.96 * 0.0008, "PDF Table 2 (Thein 2008)", function(v) list(p_NC_CC = v)),
  list("p_CC_DC", "CC to DC", qbeta(.025, 14.6, 360.2), qbeta(.975, 14.6, 360.2), "PDF Table 2 (Lim 2018)", function(v) list(p_CC_DC = v)),
  list("p_HCC", "CC/DC to HCC", qbeta(.025, 1.9, 136.1), qbeta(.975, 1.9, 136.1), "PDF Table 2 (Lim 2018, one parameter)",
       function(v) list(p_CC_HCC = v, p_DC_HCC = v)),
  list("r3_cirr", "GT3 RR, cirrhosis", 1.27, 1.46, "PDF Table 2 (derived from Kanwal 2014)",
       function(v) list(r3_NC_CC = v, r3_CC_DC = v)),
  list("r3_hcc", "GT3 RR, HCC", 1.71, 2.17, "PDF Table 2 (derived from Kanwal 2014)",
       function(v) list(r3_CC_HCC = v, r3_DC_HCC = v)),
  list("phi_CC_DC", "Post-SVR RR, CC to DC", 0.03, 0.20, "PDF Table 2 (van der Meer 2012)", function(v) list(phi_CC_DC = v)),
  list("phi_CC_HCC", "Post-SVR RR, CC to HCC", 0.16, 0.35, "PDF Table 2 (Morgan 2013)", function(v) list(phi_CC_HCC = v)),
  list("phi_DC_HCC", "Post-SVR RR, DC to HCC", 0.52, 1.00, "Pageaux 2022 (HR 0.52); base 1 assumed", function(v) list(phi_DC_HCC = v)),
  list("mu_DC", "Liver mortality, DC", qbeta(.025, 147, 984), qbeta(.975, 147, 984), "PDF Table 2 (Lim 2018)", function(v) list(mu_DC = v)),
  list("mu_HCC", "Liver mortality, HCC", qbeta(.025, 117.1, 155.2), qbeta(.975, 117.1, 155.2), "PDF Table 2 (Lim 2018)", function(v) list(mu_HCC = v)),
  list("psi_DC", "Post-SVR mortality RR, DC", 0.17, 0.96, "Pageaux 2022 (liver death HR 0.40, 0.17-0.96)", function(v) list(psi_DC = v)),
  list("psi_HCC", "Post-SVR mortality RR, HCC", 0.13, 0.88, "Dang 2020 (liver death HR 0.34, 0.13-0.88)", function(v) list(psi_HCC = v)),
  list("smr", "PWID SMR", 3.47, 6.45, "van Santen 2014 Table 3", function(v) list(smr_nonviraemic = v, smr_viraemic = v * hr0)),
  list("hr_extrahepatic", "Extrahepatic HR, viraemic", 1.23, 1.77, "Lee 2012 Table 3", function(v) list(smr_viraemic = bp$smr_nonviraemic * v)),
  list("a3NC", "SVR GT3, NC/CC", 0.968, 0.982, "PDF Table 1 (EASL 2020)", function(v) svr("a3NC", v)),
  list("aoNC", "SVR other GT, NC/CC", 0.984, 0.991, "PDF Table 1 (EASL 2020)", function(v) svr("aoNC", v)),
  list("a3DCm", "SVR GT3, DC (RBV-free)", 0.13, 0.79, "PDF Table 1 (EASL 2020; Xiao 2025)", function(v) svr("a3DCm", v)),
  list("aoDCm", "SVR other GT, DC (RBV-free)", 0.82, 0.96, "PDF Table 1 (EASL 2020; Xiao 2025)", function(v) svr("aoDCm", v)),
  list("a3DCp", "SVR GT3, HCC (RBV+)", 0.71, 0.98, "PDF Table 1 (EASL 2020; Xiao 2025)", function(v) svr("a3DCp", v)),
  list("aoDCp", "SVR other GT, HCC (RBV+)", 0.87, 0.95, "PDF Table 1 (EASL 2020; Xiao 2025)", function(v) svr("aoDCp", v)),
  list("lambda2", "Release rates", NA, NA, "Stan fit 95% interval", NULL),
  list("lambda3", "Re-arrest rates", NA, NA, "Stan fit 95% interval", NULL),
  list("pi_recid", "Return to injecting after release", pi_q(.025), pi_q(.975), "Stan fit 95% interval", function(v) list(pi_recid = v)))
names(inputs) <- vapply(inputs, `[[`, "", 1)
runs <- list(base = list())
for (x in inputs) for (side in c("low", "high")) {
  ov <- if (x[[1]] == "lambda2") list(lambda2 = rate_q("output/release_rate_fit.rds", if (side == "low") .025 else .975)) else
    if (x[[1]] == "lambda3") list(lambda3 = rate_q("output/rearrest_rate_fit.rds", if (side == "low") .025 else .975)) else
    x[[6]](if (side == "low") x[[3]] else x[[4]])
  runs[[paste(x[[1]], side, sep = "_")]] <- ov
}
ranges <- do.call(rbind, lapply(inputs, function(x) data.frame(input = x[[1]], label = x[[2]],
  base = paste(signif(unlist(bp[names(runs[[paste0(x[[1]], "_low")]])]), 4), collapse = "; "),
  low = paste(signif(unlist(runs[[paste0(x[[1]], "_low")]]), 4), collapse = "; "),
  high = paste(signif(unlist(runs[[paste0(x[[1]], "_high")]]), 4), collapse = "; "), source = x[[5]])))
write.csv(ranges, file.path(OUT, "ranges.csv"), row.names = FALSE)

run_one <- function(nm) {
  dir <- file.path(OUT, "runs", nm); dir.create(dir, showWarnings = FALSE)
  a <- c(file.path(ROOT, "scripts/scenarios_2026.R"), "--root", ROOT, "--fit", FIT, "--posterior-files", theta_file,
         "--draws-per-seed", "1", "--cores", "1", "--min-valid-share", "1", "--out-dir", file.path("output/tornado/runs", nm))
  if (length(runs[[nm]])) { f <- file.path(dir, "override.rds"); saveRDS(runs[[nm]], f); a <- c(a, "--param-file", f) }
  status <- system2("Rscript", a, stdout = file.path(dir, "log.txt"), stderr = file.path(dir, "log.txt"))
  if (status != 0) return(NULL)
  s <- read.csv(file.path(dir, "scenario_2030_summary.csv"))
  data.frame(run = nm, s[c("scenario_id", "incidence_reduction", "mortality_reduction", "hcv_mortality_reduction",
                           "who_incidence_pass", "who_mortality_pass", "who_joint_pass")])
}
res <- parallel::mclapply(names(runs), run_one, mc.cores = N_CORES, mc.preschedule = FALSE)
# A run fails when its 2015 snapshot is not at equilibrium after the 150-year burn-in (for example,
# HCV dies out at a low transmission probability when the contact scales are not refitted). Failed
# runs are listed, not dropped silently; the base run must succeed.
failed <- names(runs)[vapply(res, is.null, TRUE)]
stopifnot("base run failed" = !"base" %in% failed)
write.csv(data.frame(run = failed), file.path(OUT, "failed_runs.csv"), row.names = FALSE)
if (length(failed)) cat("runs without a 2015 equilibrium (see output/tornado/runs/<run>/log.txt):", failed, "\n")
res <- do.call(rbind, res)
res$input <- sub("_(low|high)$", "", res$run); res$side <- ifelse(res$run == "base", "base", sub("^.*_", "", res$run))
write.csv(res, file.path(OUT, "tornado_results.csv"), row.names = FALSE)

# Policies whose pass/fail (incidence, mortality, joint) differs from the base run.
base <- res[res$run == "base", ]
m <- merge(res[res$run != "base", ], base[c("scenario_id", "who_incidence_pass", "who_mortality_pass", "who_joint_pass")],
           by = "scenario_id", suffixes = c("", "_base"))
flip <- with(m, who_incidence_pass != who_incidence_pass_base | who_mortality_pass != who_mortality_pass_base |
                who_joint_pass != who_joint_pass_base)
write.csv(m[flip, c("input", "side", "scenario_id", "incidence_reduction", "mortality_reduction",
                    "who_incidence_pass_base", "who_incidence_pass", "who_mortality_pass_base", "who_mortality_pass")],
          file.path(OUT, "pass_flips.csv"), row.names = FALSE)
cat("runs:", length(runs), " pass/fail changes:", sum(flip), "\n")
