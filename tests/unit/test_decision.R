# =============================================================================
# tests/unit/test_decision.R — unit tests for the decision-analysis rules in scripts/decision_analysis.R
# (Pareto efficiency, 0-1 normalisation, TOPSIS, efficiency frontier, budget rule).
# =============================================================================

args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests/helper.R"))
source(file.path(root, "scripts/decision_analysis.R"))  # functions only; main block runs under Rscript

# Benefit (max) and cost (min): A cheap/low, B dominated by C, C dear/high, D equals A on both.
M <- cbind(benefit = c(1, 2, 3, 1), cost = c(1, 5, 4, 1))
check_identical(pareto_efficient(M, c(TRUE, FALSE)), c(TRUE, FALSE, TRUE, TRUE),
                "Pareto: strictly worse row is dominated; ties are not")
N <- normalise01(M, c(TRUE, FALSE))
check_equal(N[, "benefit"], c(0, .5, 1, 0), "normalise01: benefit scaled 0-1")
check_equal(N[, "cost"], c(1, 0, .25, 1), "normalise01: cost reverse-coded")
check_equal(normalise01(cbind(c(2, 2)), TRUE), c(1, 1), "normalise01: constant column is 1")
cl <- topsis(rbind(c(1, 1), c(0, 0), c(.5, .5)), c(1, 1))
check_equal(cl, c(1, 0, .5), "TOPSIS: ideal 1, anti-ideal 0, midpoint 0.5")
check_equal(topsis(rbind(c(1, 0), c(0, 1)), c(3, 1))[1] > .5, TRUE, "TOPSIS: weights favour the heavier criterion")
tp <- topsis_pareto(M, c(TRUE, FALSE), c(1, 1))
check_true(is.na(tp[2]) && all(!is.na(tp[-2])), "TOPSIS is computed only on the Pareto set")

# Frontier: O (0,0), A (10,1), B (30,2) is extendedly dominated by A-C mix, C (40,4), D (50,3) dominated.
fr <- efficiency_frontier(c("O", "A", "B", "C", "D"), c(0, 10, 30, 40, 50), c(0, 1, 2, 4, 3))
check_identical(fr$status, c("frontier", "frontier", "extendedly dominated", "frontier", "dominated"),
                "frontier: strong and extended dominance flagged")
check_equal(fr$icer_vs_previous[c(2, 4)], c(10, 10), "frontier: ICERs vs previous frontier strategy")
check_identical(best_within_budget(c("O", "A", "C"), c(0, 10, 40), c(0, 1, 4), 20), "A",
                "budget rule: most benefit within the cap")
check_true(is.na(best_within_budget("A", 10, 1, 5)), "budget rule: NA when nothing is affordable")

n_fail <- test_report("test_decision.R")
if (n_fail > 0) quit(status = 1)
