# =============================================================================
# targets.R — loads the 2015 calibration targets from JSON: the six observed prison
# seroprevalences (Beta likelihood, fixed precision phi), the CNB 2015 new drug abusers by age
# (first arrests, negative binomial size) and the community PWID estimate (log-normal SD).
# =============================================================================

read_targets <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package jsonlite is required to read calibration_targets.json")
  x <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  stopifnot(length(x$age_groups) == 6L, length(x$seroprevalence) == 6L,
            all(x$seroprevalence > 0 & x$seroprevalence < 1),
            length(x$beta_precision) == 1L, x$beta_precision > 0,
            length(x$first_arrests) == 6L, all(x$first_arrests >= 0),
            length(x$first_arrests_nb_size) == 1L, x$first_arrests_nb_size > 0,
            length(x$community_pwid) == 1L, x$community_pwid > 0,
            length(x$community_pwid_log_sd) == 1L, x$community_pwid_log_sd > 0)
  x
}
