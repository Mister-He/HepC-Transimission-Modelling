# =============================================================================
# tests/unit/test_likelihood.R — unit tests for the 7-parameter likelihood, prior and parameterisation
# (likelihood.R). Run: Rscript tests/unit/test_likelihood.R
# =============================================================================

.args0 <- commandArgs(FALSE)
.file_arg <- sub("^--file=", "", .args0[grepl("^--file=", .args0)])
root <- normalizePath(file.path(dirname(normalizePath(.file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))
source(file.path(root, "scripts", "calibration", "targets.R"))
targets <- read_targets(file.path(root, "data", "calibration_targets.json"))
source(file.path(root, "scripts", "calibration", "likelihood.R"))

# Synthetic base parameters: the unit test exercises the parameterisation
# (exp transforms, bounds, likelihoods), not the fitted values.
base_params <- list(
  C_contact = matrix(1, nrow = 6, ncol = 6),
  beta = c(10, 300, 100, 50, 20, 1), lambda3 = rep(0.5, 6), c_first = 1
)

check_equal(N_THETA, 7, "exactly 7 fitted parameters (six contact scales and c)")
check_true(all(CONTACT_LO > 0 & CONTACT_LO < CONTACT_HI) && C_LO > 0 && C_LO < C_HI,
           "contact and c bounds are valid")

# theta = 0  =>  all scales = 1, c = 1  =>  params unchanged
pm0 <- build_params(rep(0, N_THETA), base_params)
check_equal(pm0$C_contact, base_params$C_contact,
            "theta = 0 leaves C_contact unchanged")
check_equal(pm0$c_first, 1, "theta = 0 gives c = 1")

# exp transform: a single unit shift scales one row by e; theta7 sets c and nothing else
th1 <- rep(0, N_THETA); th1[1] <- 1
pm1 <- build_params(th1, base_params)
check_equal(pm1$C_contact[1, ], exp(1) * base_params$C_contact[1, ],
            "theta[1] = 1 scales contact row 1 by e")
check_equal(pm1$C_contact[2, ], base_params$C_contact[2, ],
            "other contact rows unchanged")
th7 <- rep(0, N_THETA); th7[7] <- log(0.3)
pm7 <- build_params(th7, base_params)
check_equal(pm7$c_first, 0.3, "theta[7] = log c sets the first-arrest multiplier")
check_equal(pm7$beta, base_params$beta, "entry stays at the CNB counts (no scaling)")
check_equal(pm7$lambda3, base_params$lambda3, "theta[7] leaves the re-arrest rate unchanged")

# bounds
check_true(within_bounds(rep(0, N_THETA)), "zero theta is within bounds")
check_true(!within_bounds(rep(50, N_THETA)), "large theta is out of bounds")
check_true(!within_bounds(c(rep(0, 6), log(0.49))) && within_bounds(c(rep(0, 6), log(0.5))),
           "c must be at least 0.5")

# Fixed prior (not centred on any fit): log contact scales Normal(0, 1.5), log c Normal(0, 0.5).
ps <- prior_spec()
check_equal(ps$mean, rep(0, N_THETA), "prior means are 0 (no data-dependent centre)")
check_equal(ps$sd, c(rep(1.5, 6), 0.5), "prior SDs: 1.5 for contact scales, 0.5 for log c")
check_equal(nlp(c(rep(0, 6), 1)), 2, "the prior term is 0.5 * sum((theta / sd)^2)")
check_equal(nlp(c(1.5, rep(0, 6))), 0.5, "a contact scale one prior SD from 0 adds 0.5")

# Beta likelihood with fixed precision phi: the observed value maximises it, and it matches dbeta.
y <- targets$seroprevalence; phi <- targets$beta_precision
check_true(nll_prev(y) < nll_prev(pmin(y + 0.05, 0.99)) && nll_prev(y) < nll_prev(y - 0.05),
           "the Beta NLL is lower at the observed seroprevalence than 0.05 away")
check_equal(nll_prev(rep(0.3, 6)), -sum(dbeta(y, 0.3 * phi, 0.7 * phi, log = TRUE)),
            "nll_prev is the Beta(p * phi, (1 - p) * phi) negative log density")
check_true(!exists("nll_pop") && !exists("pattern_penalty") && !exists("prior_centre"),
           "no population likelihood, pattern penalty or fit-centred prior")

# build_params rejects the wrong number of parameters
check_error(build_params(rep(0, 12), base_params),
            "build_params rejects length-12 theta")

# Nelder-Mead restarts never end worse and count every evaluation in the history.
source(file.path(root, "scripts", "calibration", "calibrate_nm.R"))
rosen <- function(x) sum(100 * (x[-1] - x[-N_THETA]^2)^2 + (1 - x[-N_THETA])^2)
nm0 <- run_nm_start(rep(-1, N_THETA), rosen, maxit = 500)
nm3 <- run_nm_start(rep(-1, N_THETA), rosen, maxit = 500, restarts = 3)
check_true(nm3$objective <= nm0$objective && nm3$fn_evaluations == nrow(nm3$history),
           "Nelder-Mead restarts do not worsen the fit and log every evaluation")

n_fail <- test_report("tests/unit/test_likelihood.R")
if (n_fail > 0) quit(status = 1)
