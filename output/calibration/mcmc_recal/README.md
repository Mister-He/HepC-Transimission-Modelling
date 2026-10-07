# MCMC recalibration output guide

This folder contains the corrected-target adaptive-MCMC posterior samples,
parameter summaries, diagnostics, the combined NPE covariance used to initialize
four chains, and run configuration. The generated `MCMC_STRICT_GATE.txt` is
authoritative for acceptance: every parameter must have pooled ESS above
4,000 and R-hat below 1.001. Samples from a failed gate are diagnostic only.
This run passed, with minimum pooled ESS 10,302 and maximum R-hat 1.000552.
Reruns overwrite these generated files.
