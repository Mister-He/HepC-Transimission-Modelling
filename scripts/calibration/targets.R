# =============================================================================
# targets.R — loads the 2015 age-specific calibration targets from JSON and
# derives the binomial prevalence counts.
# =============================================================================

read_targets <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package jsonlite is required to read calibration_targets.json")
  x <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  stopifnot(length(x$age_groups) == 6L,
            length(x$seroprevalence) == 6L,
            length(x$population) == 6L)
  x$x_prev <- round(x$seroprevalence * x$population)
  x$p_prev_binom <- x$x_prev / x$population
  x
}
