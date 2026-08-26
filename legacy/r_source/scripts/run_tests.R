# =============================================================================
# run_tests.R — run the unit and/or integration test suites
#
# Usage:
#   Rscript scripts/run_tests.R [--scope unit|integration|all] [--root <root>]
#
# Each tests/**/test_*.R file is executed in its own R session so that
# failures are isolated and reported per file; the runner exits non-zero if
# any test fails. Integration tests compile src/sim.cpp (Rcpp), so they
# require Rcpp and RcppArmadillo.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(name, default = NULL) {
  pos <- match(name, args)
  if (is.na(pos)) return(default)
  args[pos + 1]
}

ROOT  <- normalizePath(arg_val("--root", getwd()), mustWork = TRUE)
SCOPE <- arg_val("--scope", "all")
stopifnot(SCOPE %in% c("unit", "integration", "all"))

dirs <- switch(SCOPE,
  unit        = "tests/unit",
  integration = "tests/integration",
  all         = c("tests/unit", "tests/integration")
)

files <- unlist(lapply(dirs, function(d) {
  p <- file.path(ROOT, d)
  if (!dir.exists(p)) return(character(0))
  list.files(p, pattern = "^test_.*\\.R$", full.names = TRUE)
}))
files <- sort(unique(files))
if (length(files) == 0) {
  stop("No test files found under ", paste(dirs, collapse = ", "))
}

cat("Running", length(files), "test files (scope:", SCOPE, ")\n")
rscript <- file.path(R.home("bin"), "Rscript")
total_fail <- 0L

for (f in files) {
  cat("\n== ", f, " ==\n", sep = "")
  status <- system2(rscript, shQuote(f))
  if (status != 0) {
    cat("!! Test file exited with status", status, "\n")
    total_fail <- total_fail + 1L
  } else {
    cat("OK\n")
  }
}

if (total_fail > 0) {
  stop(sprintf("%d test file(s) failed.", total_fail))
}
cat(sprintf("\nAll %d test file(s) passed.\n", length(files)))
