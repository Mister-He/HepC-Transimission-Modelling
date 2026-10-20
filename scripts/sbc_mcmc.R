#!/usr/bin/env Rscript
# =============================================================================
# sbc_mcmc.R — small direct simulation-based calibration of the reference MCMC: reruns
# run_mcmc_reference.R on K synthetic datasets and ranks the true parameters in each
# posterior. Cross-validation only; nothing here feeds the production path.
#
# True parameters and noise-free outputs are the first K rows of the NPE SBC set
# (output/calibration/npe_lhs/theta_sbc.csv, x_sbc.csv: prior draws that pass the
# equilibrium check). Observations follow the MCMC likelihood: N_obs = N * exp(0.10 z),
# n = round(N_obs) and seropositives ~ Binomial(n, p). Each dataset gets one chain with
# the production sampler settings and --max-evals evaluations.
# Usage: Rscript scripts/sbc_mcmc.R [--root .] [--n-datasets 8] [--max-evals 30000] [--cores 8]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) { i <- match(name, args); if (is.na(i)) default else args[i + 1] }
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE); setwd(ROOT)
K <- as.integer(arg("--n-datasets", "8")); MAX_EVALS <- as.integer(arg("--max-evals", "30000"))
CORES <- as.integer(arg("--cores", "8")); SEED <- as.integer(arg("--seed", "20261201"))
FIT <- arg("--fit", "output/calibration/r_20261020_agepop/fit.rds")
NPE <- "output/calibration/npe_lhs"; OUT <- "output/calibration/npe_vs_mcmc/mcmc_sbc"
TH <- paste0("theta", 1:12)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

theta <- as.matrix(read.csv(file.path(NPE, "theta_sbc.csv"))[seq_len(K), TH])
x <- as.matrix(read.csv(file.path(NPE, "x_sbc.csv"))[seq_len(K), ])
stopifnot(nrow(theta) == K, nrow(x) == K)
source("scripts/calibration/likelihood.R")
stopifnot(isTRUE(all.equal(jsonlite::fromJSON(file.path(NPE, "simulation_summary.json"))$prior_centre,
                           prior_centre(readRDS(FIT)$best$theta))))  # same prior as the NPE SBC set
base <- jsonlite::fromJSON("data/calibration_targets.json")
set.seed(SEED)
for (k in seq_len(K)) {
  n <- pmax(1L, as.integer(round(x[k, 7:12] * exp(rnorm(6, 0, base$population_sigma_log)))))
  tg <- base; tg$population <- n; tg$seroprevalence <- rbinom(6, n, x[k, 1:6]) / n
  jsonlite::write_json(tg, file.path(OUT, sprintf("targets_ds%02d.json", k)), auto_unbox = TRUE, digits = NA,
                       pretty = TRUE)
}

run_one <- function(k) {
  d <- file.path(OUT, sprintf("ds%02d", k))
  st <- system2("Rscript", c("scripts/run_mcmc_reference.R", "--root", ".", "--fit", FIT,
                             "--targets", file.path(OUT, sprintf("targets_ds%02d.json", k)), "--out-dir", d,
                             "--n-chains", 1, "--cores", 1, "--max-evals", MAX_EVALS, "--seed", SEED + k),
                stdout = paste0(d, ".log"), stderr = paste0(d, ".log"))
  unlink(file.path(d, setdiff(list.files(d), c("mcmc_draws.csv", "mcmc_diagnostics.csv", "mcmc_run_config.json",
                                               "MCMC_STATUS.txt"))))
  st
}
status <- unlist(parallel::mclapply(seq_len(K), run_one, mc.cores = min(CORES, K), mc.preschedule = FALSE))
if (any(status != 0)) stop("MCMC SBC run(s) failed: ", paste(which(status != 0), collapse = ", "))

ranks <- t(vapply(seq_len(K), function(k) {
  dr <- as.matrix(read.csv(file.path(OUT, sprintf("ds%02d", k), "mcmc_draws.csv"))[, TH])
  colMeans(sweep(dr, 2, theta[k, ], "<"))
}, numeric(12)))
rhat <- vapply(seq_len(K), function(k)
  max(read.csv(file.path(OUT, sprintf("ds%02d", k), "mcmc_diagnostics.csv"))$rhat), numeric(1))
write.csv(data.frame(dataset = seq_len(K), max_split_rhat = rhat, setNames(as.data.frame(ranks), TH)),
          file.path(OUT, "mcmc_sbc_ranks.csv"), row.names = FALSE)
in95 <- ranks >= 0.025 & ranks <= 0.975; in90 <- ranks >= 0.05 & ranks <= 0.95
summ <- data.frame(parameter = c(TH, "all"), n_datasets = K,
                   coverage_95 = c(colMeans(in95), mean(in95)), coverage_90 = c(colMeans(in90), mean(in90)),
                   mean_rank = c(colMeans(ranks), mean(ranks)))
write.csv(summ, file.path(OUT, "mcmc_sbc_summary.csv"), row.names = FALSE)
print(summ[summ$parameter == "all", ])
cat("pooled coverage over", length(ranks), "ranks; max split R-hat", round(max(rhat), 3), "\n")
