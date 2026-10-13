#!/usr/bin/env Rscript
# =============================================================================
# update_rate_params.R — copies the Stan rate-fit summaries into scripts/model/setup.R,
# replacing the manual copy step after scripts/model/preprocess.R.
#
# lambda2 = exp(posterior mean of the release fit's log_lambda) (preprocess.R convention);
# lambda3 = posterior mean of exp(log_lambda), pi_recid = posterior mean of
# inv_logit(pi_raw) and c_composite = posterior mean of c_composite, all from the
# re-arrest fit draws. Values are written to 3 significant figures, the precision of the
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

release <- readRDS(RELEASE)
rel_sum <- as.data.frame(release$summary("log_lambda"))
rearrest <- readRDS(REARREST)  # list(summary, draws [iteration, chain, variable])
d <- rearrest[[2]]
draws <- matrix(d, ncol = dim(d)[3], dimnames = list(NULL, dimnames(d)[[3]]))
cols <- function(prefix) {
  j <- grep(paste0("^", prefix, "\\["), colnames(draws))
  if (prefix == "c_composite" && !length(j)) j <- grep("^c\\[", colnames(draws))  # fits saved before the rename
  if (!length(j)) stop("rearrest fit has no ", prefix, " draws")
  draws[, j, drop = FALSE]
}
if (!"pi_raw" %in% colnames(draws)) stop("rearrest fit has no pi_raw draws")

new <- list(
  lambda2 = exp(rel_sum$mean[order(as.integer(sub(".*\\[(\\d+)\\]", "\\1", rel_sum$variable)))]),
  lambda3 = unname(colMeans(exp(cols("log_lambda")))),
  c_composite = unname(colMeans(cols("c_composite"))),
  pi_recid = mean(plogis(draws[, "pi_raw"]))
)
for (nm in c("lambda2", "lambda3", "c_composite"))
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
