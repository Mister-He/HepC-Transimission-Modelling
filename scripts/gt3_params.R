#!/usr/bin/env Rscript
# =============================================================================
# gt3_params.R — writes the --param-file inputs for the genotype 3 share sensitivity
# (author decision 2026-11-04): rho 0.10 and 0.50 instead of 0.78, with the SVR composites
# moved accordingly (rho_overrides() in scripts/model/genotype.R). No refit.
# Run: Rscript scripts/gt3_params.R --root .
# Output: output/scenarios_2026_gt3/rho010/ and rho050/ (gt3_params.rds, gt3_params.csv);
# then run scripts/scenarios_2026.R with --param-file (docs/training.md §5).
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) { i <- match(name, args); if (is.na(i)) default else args[i + 1] }
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
FIT <- arg("--fit", "output/calibration/r_20261103_agepop/fit.rds")
source(file.path(ROOT, "scripts/model/genotype.R"))

bp <- readRDS(file.path(ROOT, FIT))$base_params
for (v in c(0.10, 0.50)) {
  dir <- file.path(ROOT, sprintf("output/scenarios_2026_gt3/rho%03d", round(100 * v)))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ov <- rho_overrides(bp, v)
  saveRDS(ov, file.path(dir, "gt3_params.rds"))
  write.csv(as.data.frame(ov), file.path(dir, "gt3_params.csv"), row.names = FALSE)
}
