library(data.table)
library(tictoc)
library(dplyr)
library(ggplot2)

# =========================================#
# 1. Preparation
# =========================================#
cat("======Loading libraries and data...======\n")
# Load required functions and data
source("functions/utilis_new.r")
strategies <- read.csv("settings/Final_strategies.csv")

# Change annual proportions to monthly proportions
strategies$pD <- strategies$pD / 12
strategies$pJ1 <- strategies$pJ1 / 12
strategies$pJ2 <- strategies$pJ2 / 12
strategies$pJ3 <- strategies$pJ3 / 12
strategies$pX <- strategies$pX / 12

# Load observed HCV prevalence by age group
obs_data <- list(
    pos = c(55, 145, 183, 164, 212, 299, 222, 190, 133),
    tot = c(307, 797, 829, 633, 598, 642, 481, 439, 366)
)

# Set up to calculate logLL, which can be used for MCMC
nsteps_permonth <- 1 # 4 # 30
intervention <- strategies[1, ]
configuration <- create_configuration(ny = 600, nsteps_permonth = nsteps_permonth) # 250
parameters <- set_parameters(nsteps_permonth = configuration$nsteps_permonth)
initialstates <- make_initial_states(configuration, equilibrium = FALSE, parameters)

# Get this from some personal guess
xparam_init_scaling <- log(c(9.0, 0.2, 0.18, 0.13, 0.44, 0.3, 0.08, 0.02, 0.01))
xparam_init_arrest <- log(c(0.067, 0.21, 0.21, 0.15, 0.2, 0.31, 0.14, 0.4, 0.025) / nsteps_permonth / 12)
xparam_init <- c(xparam_init_scaling, xparam_init_arrest)

## MCMC
# logLL fun for adaptive MCMC
LB_scaling <- log(c(8.5, 0.1, 0.1, 0.1, 0.35, 0.1, 0.06, 0.018, 0.008))
UB_scaling <- log(c(9.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.022, 0.012))

LB_arrest <- xparam_init_arrest - abs(xparam_init_arrest) * 0.1
UB_arrest <- xparam_init_arrest + abs(xparam_init_arrest) * 0.1

LB <- c(LB_scaling, LB_arrest)
UB <- c(UB_scaling, UB_arrest)

labs <- c(paste("<", 20, sep = ""), paste(seq(20, 50, by = 5), seq(20 + 5 - 1, 55 - 1, by = 5), sep = "-"), paste(55, "+", sep = ""))

# =========================================#
# 2. Markov Chain Monte Carlo (MCMC)
# =========================================#
library(mcmc)
library(parallel)
library(coda)
library(future)
library(future.apply)

# Parameters
set.seed(114514)
n_chains  <- 4
n_burnin  <- 3e+4
n_samples <- 2e+4
n_batch_length <- 100
n_batches <- (n_burnin + n_samples) / n_batch_length

# Create some sensible inits for each parallel chain
inits <- replicate(n_chains, runif(length(LB), min = LB, max = UB), simplify = FALSE)

# 2.1 Use grid search to find a good scale for pilot
target_accept <- 0.25
scales_to_try <- c(0.01, 0.02, 0.03, 0.04, 0.05)
accept_rates  <- numeric(length(scales_to_try))

cat("======Start evaluating scalar scales...======\n")
for (i in seq_along(scales_to_try)) {
    test_pilot <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 200, blen = 1, scale = scales_to_try[i])
    accept_rates[i] <- test_pilot$accept
    cat("Scale =", scales_to_try[i], ", Accept rate =", accept_rates[i], "\n")
}

# Choose the best scale
best_scale_idx <- which.min(abs(accept_rates - target_accept))
best_scale <- scales_to_try[best_scale_idx]
cat("Best scale:", best_scale, "with accept rate:", accept_rates[best_scale_idx], "\n")
# Best scale: 0.05 with accept rate: 0.249
# Turns out to be 0.05, but can try again, might vary by seeds
# best_scale <- 0.03

# Run pilot with the best scale
cat("======Run pilot with the best scale...======\n")
pilot <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 1e3, blen = 1, scale = best_scale)
cat("Final pilot acceptance rate:", pilot$accept, "\n")

# Double check if the accept rate is still sensible
if (pilot$accept > 0.5) {
    cat("Warning: Acceptance rate is still too high\n")
    # Lower scale
    pilot <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 1e3, blen = 1, scale = best_scale * 0.1)
    cat("Reduced scale pilot acceptance rate:", pilot$accept, "\n")
} else if (pilot$accept < 0.2) {
    cat("Warning: Acceptance rate is too low\n")
    # Increase scale
    pilot <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 1e3, blen = 1, scale = best_scale * 2)
    cat("Increased scale pilot acceptance rate:", pilot$accept, "\n")
}

# 2.2 Get covariance matrix to further scaling the previous best scale
cat("======Start calculating covariance matrix...======\n")
# Calculate covariance matrix, using a larger regularization
S <- cov(pilot$batch[501:1e3, ]) + 1e-6 * diag(length(xparam_init))
L <- chol(S)   
d <- length(xparam_init)
c.opt <- 2.38 / sqrt(d)
scale_mat <- c.opt * t(L)
# sds <- sqrt(diag(S))
# d <- length(xparam_init)
# cat("Standard deviations:", sds, "\n")

# # Check for any unusually large standard deviations
# if (any(sds > 100)) {
#     cat("Warning: Detected unusually large standard deviations, using fixed scale\n")
#     scale_vec <- rep(best_scale, length(xparam_init))
# } else {
#     c.opt <- 2.38 / sqrt(d)
#     scale_vec <- c.opt * sds
#     # Limit the range of scale_vec
#     scale_vec <- pmax(scale_vec, 0.001)
#     scale_vec <- pmin(scale_vec, 10)
# }
# cat("Final scale vector:", scale_vec, "\n")
# # 0.01932146 0.09241020 0.10644060 0.11404735 0.04434506 0.08024176 0.12374342 0.04834530 0.09606886

# Test final scale_vec
final_test <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 1000, blen = 1, scale = scale_mat)
cat("Final test acceptance rate:", final_test$accept, "\n")
# Final test acceptance rate: 0.292

# 2.3 Run MCMC in parallel with the determined scale_vec
# scale_vec <- c(0.01932146, 0.09241020, 0.10644060, 0.11404735, 0.04434506, 
#                0.08024176, 0.12374342, 0.04834530, 0.09606886)

cat("======Start parallel MCMC...======\n")
time1 <- Sys.time()
run_chain <- function(seed, init) {
  set.seed(seed)
  fit <- metrop(fun_logLL_MCMC, 
                initial = init, 
                nbatch = n_batches,
                blen = n_batch_length,
                scale = scale_mat)
  return(fit)
}

plan(multisession, workers = n_chains)
fits <- future_lapply(1:n_chains, function(k) run_chain(k, inits[[k]]), future.seed = TRUE)
saveRDS(fits, file = "mcmc_raw_results.rds")

# Print acceptance rates for each chain
for (k in 1:n_chains) {
    cat("Chain", k, "acceptance rate:", fits[[k]]$accept, "\n")
}

mcmc_list <- mcmc.list(lapply(fits, function(f) mcmc(f$batch[-c(1:n_burnin/n_batch_length),])))
cat("======Parallel MCMC completes!...======\n")

# =========================================#
# 3. Convergence Diagnostics
# =========================================#
rhat <- gelman.diag(mcmc_list)  # Gelman-Rubin Diagnostics
ess <- effectiveSize(mcmc_list) # Effective Sample Size

cat("======Gelman-Rubin Diagnostics======\n")
print(rhat$psrf)

cat("======Effective Sample Size======\n")
print(ess)

cat("======Bundle Results======\n")
results <- list(
    samples = mcmc_list,
    Rhat    = rhat,
    ESS     = ess
)

# Save to RDS
saveRDS(results, file = "mcmc_results.rds")

time2 <- Sys.time()
cat("======Done======\n")
cat("It takes", difftime(time2, time1, units = "mins"), "minutes in total.\n")
