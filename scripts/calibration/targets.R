# =============================================================================
# targets.R — loads the 2015 age-specific calibration targets from JSON and
# derives the binomial prevalence counts. `entry` (SingStat/CNB new drug abusers, 2015) is
# one year's entry into D, observed with log-normal error entry_sigma_log.
# =============================================================================

read_targets <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package jsonlite is required to read calibration_targets.json")
  x <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  stopifnot(length(x$age_groups) == 6L,
            length(x$seroprevalence) == 6L,
            length(x$population) == 6L,
            length(x$entry) == 6L, x$entry_sigma_log > 0)
  x$x_prev <- round(x$seroprevalence * x$population)
  x$p_prev_binom <- x$x_prev / x$population
  x
}
