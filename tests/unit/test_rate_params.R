# =============================================================================
# tests/unit/test_rate_params.R — checks scripts/model/update_rate_params.R on copies of
# setup.R: the stored rate fits reproduce setup.R exactly, and a missing entry fails
# without writing.
# =============================================================================
args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests/helper.R"))

script <- file.path(root, "scripts/model/update_rate_params.R")
setup <- file.path(root, "scripts/model/setup.R")
fits <- file.path(root, "output", c("release_rate_fit.rds", "rearrest_rate_fit.rds"))
run <- function(path) system2("Rscript", c(script, "--root", root, "--setup", path), stdout = FALSE, stderr = FALSE)

if (all(file.exists(fits))) {
  tmp <- tempfile(fileext = ".R"); file.copy(setup, tmp)
  check_identical(run(tmp), 0L, "rate handoff runs on a copy of setup.R")
  check_identical(readLines(tmp), readLines(setup), "stored rate fits reproduce the setup.R rate entries exactly")

  broken <- tempfile(fileext = ".R")
  writeLines(sub("^(\\s*)lambda3(\\s*=)", "\\1lambda3_missing\\2", readLines(setup)), broken)
  before <- readLines(broken)
  check_true(run(broken) != 0L, "rate handoff fails when an expected setup entry is missing")
  check_identical(readLines(broken), before, "a failed rate handoff writes nothing")
} else cat("SKIP: rate-fit RDS files are not present\n")

n_fail <- test_report("tests/unit/test_rate_params.R")
if (n_fail > 0) quit(status = 1)
