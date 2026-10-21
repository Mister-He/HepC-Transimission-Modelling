# =============================================================================
# targets.R — loads the 2015 age-specific calibration targets from JSON: the six observed
# seroprevalences and the fixed precision phi of their Beta likelihood (no counts).
# =============================================================================

read_targets <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package jsonlite is required to read calibration_targets.json")
  x <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  stopifnot(length(x$age_groups) == 6L, length(x$seroprevalence) == 6L,
            all(x$seroprevalence > 0 & x$seroprevalence < 1),
            length(x$beta_precision) == 1L, x$beta_precision > 0)
  x
}
