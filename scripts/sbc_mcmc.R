#!/usr/bin/env Rscript
# =============================================================================
# sbc_mcmc.R — small direct simulation-based calibration of the reference MCMC: reruns
# run_mcmc_reference.R on K synthetic datasets and ranks the true parameters in each
# posterior. Cross-validation only; nothing here feeds the production path.
#
# True parameters are the first K random prior draws (prior_sampling.R) whose simulation passes
# the equilibrium check; they and their noise-free outputs are saved as theta_true.csv and
# x_true.csv. Observations follow the MCMC likelihood: seroprevalence
# ~ Beta(p * phi, (1 - p) * phi), first arrests ~ negative binomial(size, mu), community PWID
# ~ log-normal, with the settings of the targets file. Each dataset gets one chain with the
# production sampler settings, --n-warmup warmup proposals and --max-evals evaluations.
# Usage: Rscript scripts/sbc_mcmc.R [--root .] [--n-datasets 8] [--max-evals 45000] [--cores 8]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) { i <- match(name, args); if (is.na(i)) default else args[i + 1] }
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE); setwd(ROOT)
K <- as.integer(arg("--n-datasets", "8")); MAX_EVALS <- as.integer(arg("--max-evals", "45000"))
N_WARMUP <- as.integer(arg("--n-warmup", "15000"))
CORES <- as.integer(arg("--cores", "8")); SEED <- as.integer(arg("--seed", "20261201"))
FIT <- arg("--fit", "output/calibration/r_20261030_agepop/fit.rds")
OUT <- "output/calibration/mcmc_sbc"
old <- setwd("scripts/model"); setup_env <- new.env(parent = globalenv())
sys.source("setup.R", setup_env); setwd(old)
Rcpp::sourceCpp("scripts/model/sim.cpp")
for (f in c("targets.R", "model_metrics.R", "equilibrium.R", "likelihood.R", "prior_sampling.R"))
  source(file.path("scripts/calibration", f))
TH <- paste0("theta", 1:N_THETA)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
base <- targets <- read_targets("data/calibration_targets.json")  # equilibrium_data() reads targets
check_fit_current(readRDS(FIT), setup_env$params, equilibrium_data(setup_env$data))

# Noise-free 2015 summaries (6 seroprevalences, 6 first arrests, community PWID), or NULL when the
# draw fails the equilibrium check (outside the MCMC target's support).
simulate <- function(theta) {
  r <- simulate_theta(theta, setup_env$params, equilibrium_data(setup_env$data))
  if (is.null(r) || is.null(r$eq) || !r$eq$pass) return(NULL)
  c(r$cs$p_sero, r$cs$first_arrests, r$cs$community_pwid)
}
ps <- prior_spec()
cand <- sample_prior_random(4L * K, ps$mean, ps$sd, SEED, ps$lo, ps$hi)
sims <- parallel::mclapply(seq_len(nrow(cand)), function(i) simulate(cand[i, ]), mc.cores = CORES)
ok <- which(!vapply(sims, is.null, logical(1)))
stopifnot(length(ok) >= K)  # too many equilibrium failures: raise the 4 * K candidates
ok <- ok[seq_len(K)]
theta <- cand[ok, , drop = FALSE]; colnames(theta) <- TH; x <- do.call(rbind, sims[ok])
write.csv(theta, file.path(OUT, "theta_true.csv"), row.names = FALSE)
write.csv(x, file.path(OUT, "x_true.csv"), row.names = FALSE)
cat(length(ok), "of", max(ok), "candidate prior draws pass the equilibrium check\n")
set.seed(SEED)
for (k in seq_len(K)) {
  phi <- base$beta_precision
  tg <- base; tg$seroprevalence <- rbeta(6, x[k, 1:6] * phi, (1 - x[k, 1:6]) * phi)
  tg$first_arrests <- rnbinom(6, size = base$first_arrests_nb_size, mu = x[k, 7:12])
  tg$community_pwid <- x[k, 13] * exp(rnorm(1, 0, base$community_pwid_log_sd))
  jsonlite::write_json(tg, file.path(OUT, sprintf("targets_ds%02d.json", k)), auto_unbox = TRUE, digits = NA,
                       pretty = TRUE)
}

run_one <- function(k) {
  d <- file.path(OUT, sprintf("ds%02d", k))
  st <- system2("Rscript", c("scripts/run_mcmc_reference.R", "--root", ".", "--fit", FIT,
                             "--targets", file.path(OUT, sprintf("targets_ds%02d.json", k)), "--out-dir", d,
                             "--n-chains", 1, "--cores", 1, "--n-warmup", N_WARMUP, "--max-evals", MAX_EVALS,
                             "--seed", SEED + k),
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
}, numeric(N_THETA)))
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
