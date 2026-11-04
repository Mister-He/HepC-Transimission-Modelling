# =============================================================================
# tests/unit/test_mcmc.R — the adaptive random-walk Metropolis sampler (arwm_budget, mcmc.R) recovers
# known targets: a correlated Gaussian and a half-normal truncated by -Inf outside the support
# (REVIEW_REPORT.md R-14). Run: Rscript tests/unit/test_mcmc.R
# =============================================================================

args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))
source(file.path(root, "scripts", "calibration", "mcmc.R"))

# Correlated 2-D Gaussian: mean (1, -2), SDs (1, 3), correlation 0.8.
mu <- c(1, -2); S <- matrix(c(1, 0.8 * 3, 0.8 * 3, 9), 2); Si <- solve(S)
lp_gauss <- function(th) list(lp = -0.5 * drop(t(th - mu) %*% Si %*% (th - mu)), ode = 0L)
res <- arwm_budget(lp_gauss, function() rnorm(2, 0, 3), diag(0.1, 2), diag(1e-8, 2),
                   n_warmup = 3000L, max_evals = 43000L, adapt_start = 500L, seed = 7L)
x <- res$samples
check_identical(res$status, "ok", "sampler finishes within its evaluation budget")
check_identical(res$total_evals, 43000L, "every proposal costs exactly one evaluation")
check_true(all(abs(colMeans(x) - mu) < 0.15 * sqrt(diag(S))), "Gaussian mean recovered (within 0.15 SD)")
check_true(all(abs(apply(x, 2, sd) / sqrt(diag(S)) - 1) < 0.1), "Gaussian SDs recovered (within 10%)")
check_true(abs(cor(x)[1, 2] - 0.8) < 0.05, "Gaussian correlation recovered")
acc <- mean(res$info[, "accepted"])
check_true(acc > 0.15 && acc < 0.6, "adapted acceptance rate is in a sensible range")

# Half-normal: support theta >= 0 enforced by lp = -Inf, as the calibration does for bounds and
# equilibrium failures. Mean sqrt(2/pi), SD sqrt(1 - 2/pi).
lp_half <- function(th) list(lp = if (th < 0) -Inf else -0.5 * th^2, ode = 1L)
res2 <- arwm_budget(lp_half, function() abs(rnorm(1)), matrix(0.25), matrix(1e-8),
                    n_warmup = 2000L, max_evals = 42000L, adapt_start = 500L, seed = 11L)
y <- res2$samples[, 1]
check_true(all(y >= 0), "no retained draw leaves the support")
check_true(abs(mean(y) - sqrt(2 / pi)) < 0.03 && abs(sd(y) - sqrt(1 - 2 / pi)) < 0.03,
           "truncated (half-normal) target recovered")

n_fail <- test_report("tests/unit/test_mcmc.R")
if (n_fail > 0) quit(status = 1)
