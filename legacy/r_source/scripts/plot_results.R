# =============================================================================
# plot_results.R — CLI wrapper for the publication-quality calibration figures
#
# The implementation lives in src/calibration/plot_results.R (sourced by
# run_calibration.R); this wrapper forwards command-line arguments to it.
#
# Usage:
#   Rscript scripts/plot_results.R --fit=<fit.rds> --out-dir=<dir>
# =============================================================================

args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- if (length(file_arg) && nzchar(file_arg) && file.exists(file_arg)) {
  normalizePath(file.path(dirname(normalizePath(file_arg)), ".."))
} else {
  getwd()
}
sys.source(file.path(root, "src", "calibration", "plot_results.R"),
           envir = globalenv())
