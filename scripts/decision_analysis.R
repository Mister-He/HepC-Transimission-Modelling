#!/usr/bin/env Rscript
# =============================================================================
# decision_analysis.R — multi-criteria strategy selection (as in thesis Chapter 1) from the
# per-draw cumulative outcomes written by scenarios_2026.R (scenario_cumulative_draws.csv).
# Strategy sets: the 17 main scenarios (s17) and the 17 x 4 J-coverage combinations (s68: the
# J-coverage grid plus the main run). Criteria, per horizon (2026-2030 and 2026-end year):
#   benefit: infections averted and liver deaths averted vs no treatment, WHO 2030 joint pass;
#   cost:    discounted programme budget (testing, consultation, DAA) and treatment courses.
# Rules: (a) most benefit within a budget cap; (b) lowest cost per outcome averted; (c) TOPSIS
# (0-1 min-max normalisation, cost criteria reverse-coded, equal weights) among Pareto-efficient
# strategies, with a weight sensitivity. Main results use posterior medians (and P(joint pass));
# uncertainty repeats Pareto, TOPSIS and the budget rule within each posterior draw.
# Usage: Rscript scripts/decision_analysis.R [--root .] [--out-dir output/decision_analysis]
# =============================================================================

CRITERIA <- c("infections_averted", "liver_deaths_averted", "who_joint_pass",
              "programme_cost", "treatment_courses")
MAXIMISE <- c(TRUE, TRUE, TRUE, FALSE, FALSE)

# Pareto efficiency: row i is dominated if another row is at least as good on every criterion
# and strictly better on one. M: numeric matrix (rows strategies); maximise: logical per column.
pareto_efficient <- function(M, maximise) {
  G <- sweep(as.matrix(M), 2, ifelse(maximise, 1, -1), `*`)  # larger is better in every column
  vapply(seq_len(nrow(G)), function(i) {
    gi <- matrix(G[i, ], nrow(G), ncol(G), byrow = TRUE)
    !any(rowSums(G >= gi) == ncol(G) & rowSums(G > gi) > 0)
  }, logical(1))
}

# Min-max normalisation to 0-1 with cost criteria reverse-coded (1 = best). A constant column is 1.
normalise01 <- function(M, maximise) {
  M <- as.matrix(M)
  for (j in seq_len(ncol(M))) {
    r <- range(M[, j]); x <- if (diff(r) > 0) (M[, j] - r[1]) / diff(r) else rep(1, nrow(M))
    M[, j] <- if (maximise[j]) x else 1 - x
  }
  M
}

# TOPSIS closeness on a 0-1 benefit-oriented matrix: weighted distances to the ideal (all 1) and
# the anti-ideal (all 0); weights are rescaled to sum to 1.
topsis <- function(N, w) {
  w <- w / sum(w)
  d_plus <- sqrt(colSums((w * t(1 - N))^2)); d_minus <- sqrt(colSums((w * t(N))^2))
  ifelse(d_plus + d_minus > 0, d_minus / (d_plus + d_minus), 0)
}

# TOPSIS among the Pareto-efficient rows only; returns closeness (NA off the Pareto set).
topsis_pareto <- function(M, maximise, w) {
  N <- normalise01(M, maximise); eff <- pareto_efficient(M, maximise)
  out <- rep(NA_real_, nrow(N)); out[eff] <- topsis(N[eff, , drop = FALSE], w); out
}

# Cost-effectiveness frontier for one benefit: strongly dominated (another strategy costs no more
# and gains no less, one strictly), extendedly dominated (a mix of its neighbours would beat it,
# i.e. ICERs fall along the frontier), else on the frontier with its ICER vs the previous one.
efficiency_frontier <- function(id, cost, effect) {
  n <- length(id); status <- rep("frontier", n)
  for (i in seq_len(n)) if (any(cost <= cost[i] & effect >= effect[i] &
                                (cost < cost[i] | effect > effect[i]))) status[i] <- "dominated"
  repeat {
    f <- which(status == "frontier"); f <- f[order(cost[f], -effect[f])]
    if (length(f) < 3) break
    icer <- diff(cost[f]) / diff(effect[f])
    bad <- which(head(icer, -1) > tail(icer, -1))
    if (!length(bad)) break
    status[f[bad[1] + 1]] <- "extendedly dominated"
  }
  f <- which(status == "frontier"); f <- f[order(cost[f])]
  icer <- rep(NA_real_, n); icer[f[-1]] <- diff(cost[f]) / diff(effect[f])
  data.frame(strategy = id, cost = cost, effect = effect, status = status, icer_vs_previous = icer,
             stringsAsFactors = FALSE)
}

# Rule (a): the strategy with the most benefit whose cost is within the cap (ties: lower cost).
best_within_budget <- function(id, cost, benefit, cap) {
  ok <- which(cost <= cap)
  if (!length(ok)) return(NA_character_)
  id[ok[order(-benefit[ok], cost[ok])[1]]]
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  arg <- function(name, default = NULL) { i <- match(name, args); if (is.na(i)) default else args[i + 1] }
  ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
  OUT <- file.path(ROOT, arg("--out-dir", "output/decision_analysis"))
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  N_WEIGHTS <- as.integer(arg("--weight-draws", 1000)); set.seed(as.integer(arg("--seed", 20261112)))

  coverage_of <- function(dir) {
    v <- read.table(file.path(dir, "draw_validity.txt"), col.names = c("k", "v"), as.is = TRUE)
    as.numeric(v$v[v$k == "j_annual_start_probability"])
  }
  read_run <- function(dir) {
    d <- read.csv(file.path(dir, "scenario_cumulative_draws.csv"), check.names = FALSE)
    d$coverage <- coverage_of(dir); d
  }
  main_dir <- file.path(ROOT, "output/scenarios_2026")
  jcov_dirs <- file.path(ROOT, "output/scenarios_2026_jcov", c("j000", "j050", "j090"))
  main <- read_run(main_dir)
  grid <- do.call(rbind, c(lapply(jcov_dirs, read_run), list(main)))
  stopifnot(length(unique(grid$coverage)) == 4)
  horizons <- sort(as.integer(unique(sub(".*_2026_", "", grep("^programme_cost_sgd_2026_", names(main), value = TRUE)))))
  stopifnot(2030L %in% horizons)

  # Criteria matrix of one run table at horizon h, one row per draw x strategy.
  criteria_at <- function(d, h, label) {
    col <- function(k) d[[paste0(k, "_2026_", h)]]
    data.frame(draw = d$draw, strategy = label(d),
               infections_averted = col("infections_averted_vs_no_treatment"),
               liver_deaths_averted = col("liver_deaths_averted_vs_no_treatment"),
               who_joint_pass = as.numeric(as.logical(d$who_joint_pass)),
               programme_cost = col("programme_cost_sgd"), treatment_courses = col("treatment_courses"),
               cpa_infection = col("cost_per_infection_averted_vs_no_treatment_sgd"),
               cpa_liver_death = col("cost_per_liver_death_averted_vs_no_treatment_sgd"),
               stringsAsFactors = FALSE)
  }
  sets <- list(
    s17 = list(data = main, label = function(d) d$scenario_id),
    s68 = list(data = grid, label = function(d) sprintf("%s @ J %s%%", d$scenario_id,
                                                        format(100 * d$coverage, drop0trailing = TRUE))))
  q <- function(x, p) unname(quantile(x, p, na.rm = TRUE, names = FALSE))
  w_equal <- setNames(rep(1, length(CRITERIA)), CRITERIA)
  weight_sets <- c(list(equal = w_equal),
    setNames(lapply(CRITERIA, function(k) replace(w_equal, k, 2)), paste0("double_", CRITERIA)),
    setNames(lapply(CRITERIA, function(k) replace(w_equal, k, 0)), paste0("zero_", CRITERIA)),
    list(benefits_only = c(1, 1, 1, 0, 0), costs_only = c(0, 0, 0, 1, 1)))

  summary_rows <- list(); frontier_rows <- list(); budget_rows <- list(); weight_rows <- list()
  for (set in names(sets)) for (h in horizons) {
    cr <- criteria_at(sets[[set]]$data, h, sets[[set]]$label); hz <- paste0("2026-", h)
    ids <- unique(cr$strategy); by_draw <- split(cr, cr$draw)
    stopifnot(all(vapply(by_draw, function(x) identical(x$strategy, ids), logical(1))))
    # Posterior summary per strategy; the main decision matrix uses medians and P(joint pass).
    s <- data.frame(set = set, horizon = hz, strategy = ids, stringsAsFactors = FALSE)
    for (k in c(CRITERIA[-3], "cpa_infection", "cpa_liver_death")) {
      v <- split(cr[[k]], factor(cr$strategy, ids))
      s[[k]] <- vapply(v, median, numeric(1), na.rm = TRUE)
      s[[paste0(k, "_lo")]] <- vapply(v, q, numeric(1), .025); s[[paste0(k, "_hi")]] <- vapply(v, q, numeric(1), .975)
    }
    s$who_joint_pass <- as.numeric(tapply(cr$who_joint_pass, factor(cr$strategy, ids), mean))
    M <- as.matrix(s[CRITERIA])
    s$pareto_efficient <- pareto_efficient(M, MAXIMISE)
    s$topsis_closeness <- topsis_pareto(M, MAXIMISE, w_equal)
    s$topsis_rank <- rank(-s$topsis_closeness, na.last = "keep", ties.method = "min")
    s$rank_cpa_infection <- rank(s$cpa_infection, na.last = "keep", ties.method = "min")
    s$rank_cpa_liver_death <- rank(s$cpa_liver_death, na.last = "keep", ties.method = "min")
    # Within-draw uncertainty: Pareto set and TOPSIS winner with the draw's own values.
    per_draw <- lapply(by_draw, function(x) {
      Md <- as.matrix(x[CRITERIA]); cl <- topsis_pareto(Md, MAXIMISE, w_equal)
      list(eff = pareto_efficient(Md, MAXIMISE), top = ids[which.max(cl)])
    })
    s$p_pareto_efficient <- rowMeans(vapply(per_draw, `[[`, logical(length(ids)), "eff"))
    tops <- vapply(per_draw, `[[`, character(1), "top")
    s$p_topsis_top <- as.numeric(table(factor(tops, ids))) / length(tops)
    # Random weights: Dirichlet(1) = normalised exponentials (topsis() rescales the weights).
    W <- matrix(rexp(N_WEIGHTS * length(CRITERIA)), ncol = length(CRITERIA))
    rtop <- apply(W, 1, function(w) ids[which.max(topsis_pareto(M, MAXIMISE, w))])
    s$p_top_random_weights <- as.numeric(table(factor(rtop, ids))) / N_WEIGHTS
    summary_rows[[length(summary_rows) + 1]] <- s
    # Efficiency frontiers (medians), cost vs each benefit.
    for (b in c("infections_averted", "liver_deaths_averted"))
      frontier_rows[[length(frontier_rows) + 1]] <- cbind(set = set, horizon = hz,
        benefit = b, efficiency_frontier(ids, s$programme_cost, s[[b]]))
    # Rule (a): budget caps spread over the range of median costs.
    caps <- pretty(range(s$programme_cost), n = 8); caps <- caps[caps >= min(s$programme_cost)]
    for (b in c("infections_averted", "liver_deaths_averted")) for (cap in caps) {
      pick <- best_within_budget(ids, s$programme_cost, s[[b]], cap)
      picks <- vapply(by_draw, function(x) best_within_budget(ids, x$programme_cost, x[[b]], cap), character(1))
      budget_rows[[length(budget_rows) + 1]] <- data.frame(set = set, horizon = hz,
        benefit = b, budget_cap_sgd = cap, best_strategy = pick,
        p_best_strategy = mean(picks == pick, na.rm = TRUE), stringsAsFactors = FALSE)
    }
    # TOPSIS weight sensitivity on the median matrix: fixed weight sets (random weights above).
    for (nm in names(weight_sets)) {
      cl <- topsis_pareto(M, MAXIMISE, weight_sets[[nm]])
      weight_rows[[length(weight_rows) + 1]] <- data.frame(set = set, horizon = hz,
        weights = nm, weight_vector = paste(weight_sets[[nm]], collapse = "/"),
        top_strategy = ids[which.max(cl)], top_closeness = max(cl, na.rm = TRUE), stringsAsFactors = FALSE)
    }
  }
  summ <- do.call(rbind, summary_rows)
  write.csv(summ, file.path(OUT, "decision_criteria_summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, frontier_rows), file.path(OUT, "efficiency_frontier.csv"), row.names = FALSE)
  write.csv(do.call(rbind, budget_rows), file.path(OUT, "budget_rule.csv"), row.names = FALSE)
  write.csv(do.call(rbind, weight_rows), file.path(OUT, "topsis_weight_sensitivity.csv"), row.names = FALSE)
  writeLines(c(paste("sources:", paste(sub(paste0(ROOT, "/"), "", c(main_dir, jcov_dirs), fixed = TRUE), collapse = "; ")),
               paste("posterior draws:", length(unique(main$draw))), paste("horizons: 2026-", horizons, sep = ""),
               paste("random weight draws:", N_WEIGHTS)), file.path(OUT, "decision_source.txt"))
  top <- subset(summ, topsis_rank == 1, c(set, horizon, strategy, topsis_closeness, p_topsis_top))
  print(top, row.names = FALSE)
  cat("Wrote decision analysis to", OUT, "\n")
}
