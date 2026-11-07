# =============================================================================
# genotype.R — SVR12 by genotype and the composite rates that depend on the genotype 3 share.
# setup.R stores composites alpha = rho * SVR_GT3 + (1 - rho) * SVR_other (EASL 2020; Xiao 2025);
# GT_SVR holds the per-genotype values behind them. Used by scripts/tornado.R and
# scripts/gt3_params.R. Progression-rate weighting by rho happens inside sim.cpp.
# =============================================================================
GT_SVR <- list(a3NC = 0.976, aoNC = 0.988,   # NC/CC, RBV-free
               a3DCm = 0.45, aoDCm = 0.90,   # DC, RBV-free (the regimen the scenarios use)
               a3DCp = 0.87, aoDCp = 0.92)   # DC with RBV; also used for HCC

# Overrides that set the genotype 3 share to v and move every SVR composite with it.
rho_overrides <- function(bp, v) list(rho = v,
  alpha_NC = bp$alpha_NC + (v - bp$rho) * (GT_SVR$a3NC - GT_SVR$aoNC),
  alpha_DC_neg = bp$alpha_DC_neg + (v - bp$rho) * (GT_SVR$a3DCm - GT_SVR$aoDCm),
  alpha_DC_pos = bp$alpha_DC_pos + (v - bp$rho) * (GT_SVR$a3DCp - GT_SVR$aoDCp),
  alpha_HCC = bp$alpha_HCC + (v - bp$rho) * (GT_SVR$a3DCp - GT_SVR$aoDCp))

# Overrides that set one genotype's SVR (name in GT_SVR) to v; the composite moves by that
# genotype's weight.
svr_overrides <- function(bp, name, v) {
  w <- if (startsWith(name, "a3")) bp$rho else 1 - bp$rho; d <- w * (v - GT_SVR[[name]])
  switch(substr(name, 3, 5), NC = list(alpha_NC = bp$alpha_NC + d), DCm = list(alpha_DC_neg = bp$alpha_DC_neg + d),
         DCp = list(alpha_DC_pos = bp$alpha_DC_pos + d, alpha_HCC = bp$alpha_HCC + d))
}
