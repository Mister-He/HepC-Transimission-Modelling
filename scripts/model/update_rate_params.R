#!/usr/bin/env Rscript
# =============================================================================
# update_rate_params.R — copies the Stan rate-fit summaries into scripts/model/setup.R,
# replacing the manual copy step after scripts/model/preprocess.R.
#
# Both fits are list(summary, draws [iteration, chain, variable]) from preprocess.R.
# lambda2 and lambda3 = posterior mean of exp(log_lambda) from the release and re-arrest draws;
# pi_recid = posterior mean of inv_logit(pi_raw). First arrest is not estimated here
# (lambda1 = c_first * lambda3, c_first calibrated). Values are written to 3 significant figures, the precision of the
# existing entries. Every value is validated (six finite positive ages, pi in (0, 1)) and
# every target line must occur exactly once before the file is rewritten; otherwise
# nothing is written. Only the value on each line changes; comments are kept.
# Usage: Rscript scripts/model/update_rate_params.R [--root .] [--setup scripts/model/setup.R]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) { i <- match(name, args); if (is.na(i)) default else args[i + 1] }
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
SETUP <- normalizePath(arg("--setup", file.path(ROOT, "scripts/model/setup.R")), mustWork = TRUE)
RELEASE <- arg("--release-fit", file.path(ROOT, "output/release_rate_fit.rds"))
REARREST <- arg("--rearrest-fit", file.path(ROOT, "output/rearrest_rate_fit.rds"))

as_draws <- function(fit) { d <- fit[[2]]; matrix(d, ncol = dim(d)[3], dimnames = list(NULL, dimnames(d)[[3]])) }
rel_draws <- as_draws(readRDS(RELEASE))
draws <- as_draws(readRDS(REARREST))
cols <- function(prefix, draws) {
  j <- grep(paste0("^", prefix, "\\["), colnames(draws))
  if (!length(j)) stop("rate fit has no ", prefix, " draws")
  draws[, j, drop = FALSE]
}
if (!"pi_raw" %in% colnames(draws)) stop("rearrest fit has no pi_raw draws")

new <- list(
  lambda2 = unname(colMeans(exp(cols("log_lambda", rel_draws)))),
  lambda3 = unname(colMeans(exp(cols("log_lambda", draws)))),
  pi_recid = mean(plogis(draws[, "pi_raw"]))
)
for (nm in c("lambda2", "lambda3"))
  if (length(new[[nm]]) != 6 || any(!is.finite(new[[nm]])) || any(new[[nm]] <= 0))
    stop(nm, " must be six finite positive values")
if (!is.finite(new$pi_recid) || new$pi_recid <= 0 || new$pi_recid >= 1) stop("pi_recid must lie in (0, 1)")

lines <- readLines(SETUP)
fmt <- function(v) if (length(v) == 1) format(signif(v, 3)) else
  paste0("c(", paste(as.character(signif(v, 3)), collapse = ", "), ")")
for (nm in names(new)) {
  pat <- if (nm == "pi_recid") "^(\\s*pi_recid\\s*=\\s*)[0-9.eE+-]+" else
    paste0("^(\\s*", nm, "\\s*=\\s*)c\\([^)]*\\)")
  hit <- grep(pat, lines)
  if (length(hit) != 1) stop("expected exactly one '", nm, "' entry in ", SETUP, ", found ", length(hit))
  lines[hit] <- sub(pat, paste0("\\1", fmt(new[[nm]])), lines[hit])
}
writeLines(lines, SETUP)
for (nm in names(new)) cat(nm, "=", fmt(new[[nm]]), "\n")
