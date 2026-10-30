# =============================================================================
# tests/unit/test_targets.R — unit tests for loading the 2015 calibration targets and the CNB data.
# =============================================================================

args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests/helper.R"))
source(file.path(root, "scripts/calibration/targets.R"))
targets <- read_targets(file.path(root, "data/calibration_targets.json"))

check_equal(targets$target_year, 2015, "targets describe the 2015 snapshot")
check_equal(length(targets$age_groups), 6, "six age groups are loaded")
check_true(all(targets$seroprevalence > 0 & targets$seroprevalence < 1),
           "six seroprevalence targets in (0, 1)")
check_true(targets$beta_precision > 0, "the Beta likelihood has a positive fixed precision")
check_true(is.null(targets$population) && is.null(targets$x_prev),
           "no prison population or screened-count fields")
sg <- read.csv(file.path(root, "data/singstat_new_drug_abusers.csv"))
dr <- read.csv(file.path(root, "data/singstat_drug_abusers_by_drug.csv"))
new <- dr[dr$year == 2015 & dr$abuser_group == "new", ]
pwid_share <- sum(new$count[new$main_drug %in% c("Heroin", "Methamphetamine")]) / new$count[new$main_drug == "Total"]
check_equal(new$count[new$main_drug == "Total"], sum(sg$new_drug_abusers[sg$year == 2015]),
            "drug-type total equals the age-split total of 2015 new drug abusers")
check_equal(pwid_share, 1127 / 1309, tol = 1e-12, "heroin + methamphetamine share of 2015 new drug abusers")
check_equal(targets$first_arrests, round(sg$new_drug_abusers[sg$year == 2015] * pwid_share),
            "first-arrest target = CNB 2015 new drug abusers x PWID share, rounded")
setup_env <- new.env(); old <- setwd(file.path(root, "scripts/model")); sys.source("setup.R", setup_env); setwd(old)
check_equal(setup_env$params$beta, sg$new_drug_abusers[sg$year == 2015] * pwid_share, tol = 1e-12,
            "entry into D (setup.R beta) = CNB 2015 new drug abusers x PWID share, unrounded")
check_equal(targets$community_pwid, 11000, "community PWID target = Teo 2019 estimate")
check_equal(targets$community_pwid_log_sd, (log(17000) - log(6500)) / (2 * 1.959964), tol = 1e-3,
            "community log SD comes from the 6,500-17,000 interval")
check_true(!"release_sero_2017_target" %in% names(targets),
           "the former 47 percent target is absent")
check_true(is.null(targets$entry) && is.null(targets$population_sigma_log),
           "no entry or population likelihood settings")
check_true(file.exists(file.path(root, "data/singstat_new_drug_abusers.csv")),
           "SingStat new drug abuser data (Table S2 source) are versioned in the repository")

n_fail <- test_report("tests/unit/test_targets.R")
if (n_fail > 0) quit(status = 1)
