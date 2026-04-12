# tests/testthat/test-sim.R
# Unit tests for the HCV PWID compartmental model
# Run via: testthat::test_file("tests/testthat/test-sim.R")

library(testthat)
library(Rcpp)
sourceCpp("sim.cpp")
source("setup.R")

# short data fixture for fast tests
data_short        <- data
data_short$t_end  <- 2.0
data_short$dt     <- 1 / 52

# =============================================================================
test_that("idx() helper produces unique indices for all valid inputs", {
  indices <- c()
  for (s in 0:3) for (k in 1:4) for (h in 0:3) for (i in 0:8)
    indices <- c(indices, idx(s, k, h, i))
  expect_equal(length(unique(indices)), 576)
  expect_equal(min(indices), 1L)
  expect_equal(max(indices), 576L)
})

# =============================================================================
test_that("run_sim() returns a numeric matrix with correct dimensions", {
  out <- run_sim(params_s1, data_short)
  n_steps <- as.integer((data_short$t_end - data_short$t_start) /
                          data_short$dt) + 1L
  expect_true(is.matrix(out))
  expect_equal(nrow(out), n_steps)
  expect_equal(ncol(out), 577L)   # 576 compartments + time column
})

# =============================================================================
test_that("no compartment goes negative during a 2-year no-treatment run", {
  out <- run_sim(params_s1, data_short)
  expect_true(all(out[, -1] >= -1e-9),
    info = "Negative compartment values detected (no-treatment scenario)")
})

# =============================================================================
test_that("no NaN or Inf values in output", {
  out <- run_sim(params_s3, data_short)
  expect_true(all(is.finite(out)),
    info = "NaN or Inf detected in simulation output")
})

# =============================================================================
test_that("total population is finite and positive at every time step", {
  out        <- run_sim(params_s3, data_short)
  total_pop  <- rowSums(out[, -1])
  expect_true(all(is.finite(total_pop)))
  expect_true(all(total_pop > 0))
})

# =============================================================================
test_that("time column is monotonically increasing", {
  out <- run_sim(params_s1, data_short)
  diffs <- diff(out[, 1])
  expect_true(all(diffs > 0),
    info = "Time column is not strictly increasing")
})

# =============================================================================
test_that("first time step equals t_start", {
  out <- run_sim(params_s1, data_short)
  expect_equal(out[1, 1], data_short$t_start, tolerance = 1e-10)
})

# =============================================================================
test_that("treatment scenario 4 produces lower chronic burden than scenario 1", {
  out_no_tx <- run_sim(params_s1, data_short)
  out_tx    <- run_sim(params_s4, data_short)

  # sum all chronic (h=2) compartments at final time step
  # chronic columns: every 3rd state in each (stratum, stage, age) block
  chronic_cols <- which(
    sapply(seq_len(576), function(flat) {
      # recover h from flat index: h = ((flat-1) %% 36) %/% 9
      h <- ((flat - 1L) %% 36L) %/% 9L
      h == 2L
    })
  ) + 1L   # +1 for the time column offset

  chronic_no_tx <- sum(out_no_tx[nrow(out_no_tx), chronic_cols])
  chronic_tx    <- sum(out_tx   [nrow(out_tx),    chronic_cols])

  expect_lt(chronic_tx, chronic_no_tx,
    label = "Treatment (s4) should reduce chronic burden vs no treatment (s1)")
})

# =============================================================================
test_that("zero entry rate produces declining total population", {
  p_no_entry        <- params_s1
  p_no_entry$beta   <- rep(0, 9)
  out               <- run_sim(p_no_entry, data_short)
  total_pop         <- rowSums(out[, -1])
  expect_lt(total_pop[nrow(out)], total_pop[1],
    label = "With no entry, total PWID population should decline over time")
})

# =============================================================================
test_that("higher tau leads to more treated individuals at equilibrium", {
  p_low  <- modifyList(params_s1, list(tau = c(0.1, 0.1, 0.0, 0.0)))
  p_high <- modifyList(params_s1, list(tau = c(0.9, 0.9, 0.0, 0.0)))

  out_low  <- run_sim(p_low,  data_short)
  out_high <- run_sim(p_high, data_short)

  # treated (h=3) compartment sum at final step
  treated_cols <- which(
    sapply(seq_len(576), function(flat) {
      h <- ((flat - 1L) %% 36L) %/% 9L
      h == 3L
    })
  ) + 1L

  tx_low  <- sum(out_low [nrow(out_low),  treated_cols])
  tx_high <- sum(out_high[nrow(out_high), treated_cols])

  expect_gt(tx_high, tx_low,
    label = "Higher tau should produce more treated individuals")
})

# =============================================================================
test_that("rho = 1.0 (all GT3) gives higher progression than rho = 0.0", {
  p_all_gt3    <- modifyList(params_s1, list(rho = 1.0))
  p_no_gt3     <- modifyList(params_s1, list(rho = 0.0))

  out_gt3      <- run_sim(p_all_gt3, data_short)
  out_other    <- run_sim(p_no_gt3,  data_short)

  # compare DC (stage 3) chronic compartment sum at final step
  dc_chronic_cols <- which(
    sapply(seq_len(576), function(flat) {
      k <- ((flat - 1L) %/% 36L) %% 4L + 1L
      h <- ((flat - 1L) %% 36L) %/% 9L
      k == 3L && h == 2L
    })
  ) + 1L

  dc_gt3  <- sum(out_gt3  [nrow(out_gt3),   dc_chronic_cols])
  dc_other<- sum(out_other[nrow(out_other),  dc_chronic_cols])

  expect_gt(dc_gt3, dc_other,
    label = "All-GT3 cohort should have higher DC burden than all-other-GT cohort")
})

# =============================================================================
test_that("SVR parameters are within (0, 1)", {
  svr_params <- c("alpha_NC", "alpha_DC_pos", "alpha_DC_neg", "alpha_HCC")
  for (pname in svr_params) {
    val <- params[[pname]]
    expect_gt(val, 0, label = paste(pname, "> 0"))
    expect_lt(val, 1, label = paste(pname, "< 1"))
  }
})

# =============================================================================
test_that("initial conditions vector has correct length", {
  expect_equal(length(data$y0), 576L)
})

# =============================================================================
test_that("simulation is reproducible across two identical calls", {
  out1 <- run_sim(params_s3, data_short)
  out2 <- run_sim(params_s3, data_short)
  expect_equal(out1, out2)
})
