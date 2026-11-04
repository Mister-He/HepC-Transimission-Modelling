# =============================================================================
# tests/integration/test_independent_rhs.R — checks sim.cpp against an independent R derivation of the
# full right-hand side (written from the model contract, not from sim.cpp), the liver/excess split of
# HCV deaths, and the whole-step count of run_sim (REVIEW_REPORT.md R-01, R-14, R-15).
# =============================================================================

args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))
old <- setwd(file.path(root, "scripts", "model"))
e <- new.env(parent = globalenv()); sys.source("setup.R", e); setwd(old)
Rcpp::sourceCpp(file.path(root, "scripts", "model", "sim.cpp"))
sys.source(file.path(root, "scripts", "calibration", "model_metrics.R"), envir = environment())

# State y[i, h, k, s]: age i (fastest), HCV state h (u, a, c, t, r), stage k (NC, CC, DC, HCC),
# stratum s (D, J, F, X) -- the column-major reading of the C++ index s*120 + (k-1)*30 + h*6 + i.
indep_rhs <- function(y, P) {
  Y <- array(y, c(6, 5, 4, 4)); dY <- array(0, dim(Y))
  U <- 1; A <- 2; C <- 3; T <- 4; R <- 5; Ds <- 1; Js <- 2; Fs <- 3; Xs <- 4
  inj <- Y[, , , Ds] + Y[, , , Fs]
  active <- apply(inj, 1, sum); infectious <- apply(inj[, c(A, C), , drop = FALSE], 1, sum)
  gam <- P$q * as.vector(P$C_contact %*% ifelse(active > 0, infectious / active, 0))
  g3 <- P$rho * c(P$r3_NC_CC, P$r3_CC_DC, P$r3_CC_HCC, P$r3_DC_HCC) + 1 - P$rho
  pc <- g3 * c(P$p_NC_CC, P$p_CC_DC, P$p_CC_HCC, P$p_DC_HCC)  # chronic: 1->2, 2->3, 2->4, 3->4
  pr <- c(0, P$phi_CC_DC, P$phi_CC_HCC, P$phi_DC_HCC) * pc      # recovered: no progression from NC
  alpha <- c(P$alpha_NC, P$alpha_NC, P$alpha_DC, P$alpha_HCC)
  l1 <- P$c_first * P$lambda3; l2 <- P$lambda2; l3 <- P$lambda3
  for (s in 1:4) for (k in 1:4) for (h in 1:5) {
    smr <- if (h %in% c(A, C, T)) P$smr_viraemic else P$smr_nonviraemic
    liver <- if (k == 3) P$mu_DC * (if (h == R) P$psi_DC else 1) else
      if (k == 4) P$mu_HCC * (if (h == R) P$psi_HCC else 1) else 0
    dY[, h, k, s] <- dY[, h, k, s] - (P$mu * smr + liver) * Y[, h, k, s]
    fD <- if (s == 1) l1 * Y[, h, k, Ds] else 0
    if (s == 1) {  # stratum moves, once per (k, h)
      fF <- l3 * Y[, h, k, Fs]; fJ <- l2 * Y[, h, k, Js]
      dY[, h, k, Ds] <- dY[, h, k, Ds] - fD
      dY[, h, k, Fs] <- dY[, h, k, Fs] - fF + P$pi_recid * fJ
      dY[, h, k, Js] <- dY[, h, k, Js] + fD + fF - fJ
      dY[, h, k, Xs] <- dY[, h, k, Xs] + (1 - P$pi_recid) * fJ
    }
  }
  dY[, U, 1, Ds] <- dY[, U, 1, Ds] + P$beta
  tau_ok <- seq_len(6) >= P$tau_min_age & seq_len(6) <= P$tau_max_age
  for (s in 1:4) for (k in 1:4) {
    if (s %in% c(Ds, Fs)) {
      iu <- gam * Y[, U, k, s]; ir <- gam * Y[, R, k, s]
      dY[, U, k, s] <- dY[, U, k, s] - iu; dY[, R, k, s] <- dY[, R, k, s] - ir
      dY[, A, k, s] <- dY[, A, k, s] + iu + ir
    }
    a <- Y[, A, k, s] / P$iota1
    dY[, A, k, s] <- dY[, A, k, s] - a
    dY[, R, k, s] <- dY[, R, k, s] + P$kappa * a; dY[, C, k, s] <- dY[, C, k, s] + (1 - P$kappa) * a
    st <- P$tau_by_stratum[(s - 1) * 4 + k] * tau_ok * Y[, C, k, s]; fin <- Y[, T, k, s] / P$iota2
    dY[, C, k, s] <- dY[, C, k, s] - st + (1 - alpha[k]) * fin
    dY[, T, k, s] <- dY[, T, k, s] + st - fin
    dY[, R, k, s] <- dY[, R, k, s] + alpha[k] * fin
  }
  for (s in 1:4) for (m in list(c(1, 2, 1), c(2, 3, 2), c(2, 4, 3), c(3, 4, 4))) {
    fc <- pc[m[3]] * Y[, C, m[1], s]; fr <- pr[m[3]] * Y[, R, m[1], s]
    dY[, C, m[1], s] <- dY[, C, m[1], s] - fc; dY[, C, m[2], s] <- dY[, C, m[2], s] + fc
    dY[, R, m[1], s] <- dY[, R, m[1], s] - fr; dY[, R, m[2], s] <- dY[, R, m[2], s] + fr
  }
  as.vector(dY)
}

# One RK4 day followed by the split Euler ageing step (rate 0.1/yr), from a random state with treatment
# by stratum and stage, an age window and the RBV-free DC branch switched on.
set.seed(20261103)
P <- e$params
P$C_contact <- P$C_contact * exp(rnorm(6, 0, 0.5)); P$beta <- P$beta * 1.3; P$c_first <- 2.1
P$tau_by_stratum <- runif(16, 0, 5); P$tau_min_age <- 2L; P$tau_max_age <- 5L; P$alpha_DC <- P$alpha_DC_neg
y <- runif(480, 0, 50)
d <- e$data; d$y0 <- y; d$t_end <- 1 / 365
cpp <- run_sim(P, d)[2, -1]
h <- 1 / 365; f <- function(x) indep_rhs(x, P)
k1 <- f(y); k2 <- f(y + h / 2 * k1); k3 <- f(y + h / 2 * k2); k4 <- f(y + h * k3)
Y <- array(pmax(y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 0), c(6, 5, 4, 4))
moved <- 0.1 * h * Y[1:5, , , , drop = FALSE]
Y[1:5, , , ] <- Y[1:5, , , ] - moved; Y[2:6, , , ] <- Y[2:6, , , ] + moved
check_true(max(abs(cpp - as.vector(Y)) / pmax(abs(cpp), 1e-8)) < 1e-12,
           "sim.cpp one-day step equals an independently derived RHS (RK4 + split ageing)")
no_entry <- modifyList(P, list(beta = rep(0, 6)))
check_equal(sum(f(y)), sum(P$beta) + sum(indep_rhs(y, no_entry)), tol = 1e-10,
            "independent RHS: total change = entry - deaths (mass balance)")

# WHO outcome = liver deaths; the secondary outcome adds the viraemic extrahepatic excess.
r <- rate_components(y, 0, P)
excess <- sum(array(y, c(6, 5, 4, 4))[, 2:4, , ] * P$mu * (P$smr_viraemic - P$smr_nonviraemic))
check_equal(r[["hcv_mortality"]], r[["mortality"]] + excess, tol = 1e-12,
            "liver + extrahepatic HCV deaths split exactly into liver deaths and the viraemic excess")
check_true(r[["mortality"]] > 0 && r[["mortality"]] < r[["hcv_mortality"]],
           "liver deaths (WHO outcome) are positive and below liver + extrahepatic deaths")

# Whole-step count: 0.3 / 0.1 is 2.9999999999999996 in floating point; the run must still reach 0.3.
short <- e$data; short$dt <- 0.1; short$t_end <- 0.3
check_equal(tail(run_sim(e$params, short)[, 1], 1), 0.3, tol = 1e-12,
            "run_sim keeps the last step when t_end / dt is an integer up to rounding")
short$t_end <- 0.25
check_equal(tail(run_sim(e$params, short)[, 1], 1), 0.2, tol = 1e-12,
            "a genuinely fractional t_end / dt still truncates to whole steps")

n_fail <- test_report("tests/integration/test_independent_rhs.R")
if (n_fail > 0) quit(status = 1)
