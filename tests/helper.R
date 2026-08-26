# =============================================================================
# tests/helper.R — dependency-free assertion helpers for test_*.R files
#
# Each test file sources this helper, runs `check(...)` assertions, and ends
# with:
#   n_fail <- test_report(basename(<this file>))
#   if (n_fail > 0) quit(status = 1)
# =============================================================================

if (!exists(".__TEST__", envir = globalenv())) {
  .__TEST__ <- new.env(parent = emptyenv())
  .__TEST__$checks <- 0L
  .__TEST__$failures <- 0L
  .__TEST__$current_file <- NULL
}

check <- function(cond, label) {
  .__TEST__$checks <- .__TEST__$checks + 1L
  if (isTRUE(cond)) {
    cat(sprintf("  [PASS] %s\n", label))
  } else {
    .__TEST__$failures <- .__TEST__$failures + 1L
    cat(sprintf("  [FAIL] %s\n", label))
  }
  invisible(isTRUE(cond))
}

check_true <- function(x, label) check(isTRUE(x), label)

check_identical <- function(x, y, label) check(identical(x, y), label)

check_equal <- function(x, y, label, tol = 1e-8) {
  check(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)), label)
}

check_error <- function(expr, label) {
  ok <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  check(ok, label)
}

test_report <- function(file = "tests") {
  cat(sprintf("\n[SUMMARY] %s: %d checks, %d failures\n",
              file, .__TEST__$checks, .__TEST__$failures))
  fails <- .__TEST__$failures
  .__TEST__$checks <- 0L
  .__TEST__$failures <- 0L
  fails
}

# Locate the repository root from the currently executing test file.
test_root <- function() {
  args0 <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
  if (length(file_arg) && nzchar(file_arg) && file.exists(file_arg)) {
    normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
  } else {
    normalizePath(getwd())
  }
}
