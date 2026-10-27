# =============================================================================
# tests/unit/test_likelihood.R — unit tests for the 13-parameter likelihood, prior and parameterisation
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

check_equal(N_THETA, 13, "exactly 13 fitted parameters (six contact scales, six entry scales and c)")
check_true(all(CONTACT_LO > 0 & CONTACT_LO < CONTACT_HI) && C_LO > 0 && C_LO < C_HI &&
             all(ENTRY_LO == 1 & ENTRY_HI > ENTRY_LO),
           "contact, entry and c bounds are valid; entry scales are at least 1")

# theta = 0  =>  all scales = 1, c = 1  =>  params unchanged
pm0 <- build_params(rep(0, N_THETA), base_params)
check_equal(pm0$C_contact, base_params$C_contact,
            "theta = 0 leaves C_contact unchanged")
check_equal(pm0$c_first, 1, "theta = 0 gives c = 1")
check_equal(pm0$beta, base_params$beta, "theta = 0 keeps entry at the CNB counts")

# exp transform: a single unit shift scales one row by e; theta8 scales entry of age group 2 only;
# theta13 sets c and nothing else
th1 <- rep(0, N_THETA); th1[1] <- 1
pm1 <- build_params(th1, base_params)
check_equal(pm1$C_contact[1, ], exp(1) * base_params$C_contact[1, ],
            "theta[1] = 1 scales contact row 1 by e")
check_equal(pm1$C_contact[2, ], base_params$C_contact[2, ],
            "other contact rows unchanged")
th8 <- rep(0, N_THETA); th8[8] <- log(2)
pm8 <- build_params(th8, base_params)
check_equal(pm8$beta, base_params$beta * c(1, 2, 1, 1, 1, 1), "theta[8] = log 2 doubles entry at age 20-29 only")
check_equal(pm8$C_contact, base_params$C_contact, "entry scales leave the contact matrix unchanged")
th13 <- rep(0, N_THETA); th13[13] <- log(0.6)
pm13 <- build_params(th13, base_params)
check_equal(pm13$c_first, 0.6, "theta[13] = log c sets the first-arrest multiplier")
check_equal(pm13$beta, base_params$beta, "theta[13] leaves entry unchanged")
check_equal(pm13$lambda3, base_params$lambda3, "theta[13] leaves the re-arrest rate unchanged")

# bounds
check_true(within_bounds(rep(0, N_THETA)), "zero theta is within bounds")
check_true(!within_bounds(rep(50, N_THETA)), "large theta is out of bounds")
check_true(!within_bounds(c(rep(0, 12), log(0.49))) && within_bounds(c(rep(0, 12), log(0.5))),
           "c must be at least 0.5")
check_true(!within_bounds(c(rep(0, 6), -0.01, rep(0, 6))), "an entry scale below 1 is out of bounds")

# Fixed prior (not centred on any fit): log contact scales Normal(0, 1.5), log entry scales
# Normal(0, 1) truncated at 0 (half-normal), log c Normal(0, 0.5).
ps <- prior_spec()
check_equal(ps$mean, rep(0, N_THETA), "prior means are 0 (no data-dependent centre)")
check_equal(ps$sd, c(rep(1.5, 6), rep(1, 6), 0.5), "prior SDs: 1.5 contact, 1 entry, 0.5 log c")
check_equal(nlp(c(rep(0, 12), 1)), 2, "the prior term is 0.5 * sum((theta / sd)^2)")
check_equal(nlp(c(1.5, rep(0, 12))), 0.5, "a contact scale one prior SD from 0 adds 0.5")

# Beta likelihood with fixed precision phi: the observed value maximises it, and it matches dbeta.
y <- targets$seroprevalence; phi <- targets$beta_precision
check_true(nll_prev(y) < nll_prev(pmin(y + 0.05, 0.99)) && nll_prev(y) < nll_prev(y - 0.05),
           "the Beta NLL is lower at the observed seroprevalence than 0.05 away")
check_equal(nll_prev(rep(0.3, 6)), -sum(dbeta(y, 0.3 * phi, 0.7 * phi, log = TRUE)),
            "nll_prev is the Beta(p * phi, (1 - p) * phi) negative log density")
check_true(!exists("nll_pop") && !exists("pattern_penalty") && !exists("prior_centre"),
           "no population likelihood, pattern penalty or fit-centred prior")

# First arrests: negative binomial around the model flow; community PWID: log-normal.
fa <- targets$first_arrests; size <- targets$first_arrests_nb_size
check_equal(nll_first_arrests(fa * 1.2), -sum(dnbinom(fa, size = size, mu = fa * 1.2, log = TRUE)),
            "nll_first_arrests is the negative binomial negative log density")
check_true(nll_first_arrests(fa) < nll_first_arrests(fa * 1.3) && nll_first_arrests(fa) < nll_first_arrests(fa * 0.7),
           "the first-arrest NLL is lowest near the observed counts")
m <- targets$community_pwid; sdl <- targets$community_pwid_log_sd
check_equal(nll_community(m), -dnorm(0, 0, sdl, log = TRUE), "nll_community is the log-normal density at the estimate")
check_equal(nll_community(m * exp(sdl)) - nll_community(m), 0.5, "one log SD away adds 0.5")
check_equal(nll_targets(list(p_sero = targets$seroprevalence, first_arrests = fa, community_pwid = m)),
            nll_prev(targets$seroprevalence) + nll_first_arrests(fa) + nll_community(m),
            "nll_targets sums the three terms")

# build_params rejects the wrong number of parameters
check_error(build_params(rep(0, 7), base_params),
            "build_params rejects length-7 theta")

# Nelder-Mead restarts never end worse and count every evaluation in the history.
source(file.path(root, "scripts", "calibration", "calibrate_nm.R"))
rosen <- function(x) sum(100 * (x[-1] - x[-N_THETA]^2)^2 + (1 - x[-N_THETA])^2)
nm0 <- run_nm_start(rep(-1, N_THETA), rosen, maxit = 500)
nm3 <- run_nm_start(rep(-1, N_THETA), rosen, maxit = 500, restarts = 3)
check_true(nm3$objective <= nm0$objective && nm3$fn_evaluations == nrow(nm3$history),
           "Nelder-Mead restarts do not worsen the fit and log every evaluation")

n_fail <- test_report("tests/unit/test_likelihood.R")
if (n_fail > 0) quit(status = 1)
