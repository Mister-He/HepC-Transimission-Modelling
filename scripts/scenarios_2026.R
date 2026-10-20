#!/usr/bin/env Rscript
# =============================================================================
# scenarios_2026.R — 2026–2030 policy analysis (J treatment from 2019, F/X screening
# quotas from 2026): runs every scenario per posterior draw and reports medians and 95% UIs.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) { i <- match(name, args); if (is.na(i)) default else args[i+1] }
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
FIT <- normalizePath(arg("--fit", file.path(ROOT, "output/calibration/r_20261020_agepop/fit.rds")), mustWork = TRUE)
OUT <- file.path(ROOT, arg("--out-dir", "output/scenarios_2026"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
config <- jsonlite::fromJSON(file.path(ROOT, "data/scenario_policy_config.json"))
stopifnot(config$policy_start_year == 2026, config$j_start_year == 2019,
          config$end_year == 2030,
          identical(as.numeric(config$annual_screening_totals),c(600,1000,2000,5000,8000)),
          config$discount_base_year==2026,config$discount_rate==0.03)
old <- setwd(file.path(ROOT, "scripts/model")); setup_env <- new.env(parent = globalenv())
sys.source("setup.R", setup_env); setwd(old)
Rcpp::sourceCpp(file.path(ROOT, "scripts/model/sim.cpp"))
source(file.path(ROOT, "scripts/calibration/targets.R"))
source(file.path(ROOT, "scripts/calibration/model_metrics.R"))
source(file.path(ROOT, "scripts/calibration/equilibrium.R"))
source(file.path(ROOT, "scripts/calibration/likelihood.R"))
targets <- read_targets(file.path(ROOT, "data/calibration_targets.json"))
fit <- readRDS(FIT)
# --posterior-files: comma-separated CSVs with columns theta1..theta12; each file is one
# sampling group. Default: the reference MCMC posterior.
posterior_files <- strsplit(arg("--posterior-files", "output/calibration/mcmc_reference/mcmc_draws.csv"), ",")[[1]]
posterior_files <- ifelse(grepl("^/", posterior_files), posterior_files, file.path(ROOT, posterior_files))
posterior_files <- normalizePath(posterior_files, mustWork = TRUE)
# Equal numbers of draws are taken from every posterior file.
DRAWS_PER_SEED <- as.integer(arg("--draws-per-seed", 2000)); N_CORES <- as.integer(arg("--cores", 8))
SAMPLE_SEED <- as.integer(arg("--sample-seed", 20261005)); set.seed(SAMPLE_SEED)
theta_cols <- paste0("theta", 1:12)
posterior_all <- lapply(posterior_files, function(f) read.csv(f)[, theta_cols, drop = FALSE])
draw_table <- do.call(rbind, lapply(seq_along(posterior_all), function(k) {
  d <- posterior_all[[k]]; d <- d[sample.int(nrow(d), min(DRAWS_PER_SEED, nrow(d))), , drop = FALSE]
  data.frame(seed = k, draw = seq_len(nrow(d)), d)
}))

idx0 <- function(s, k, h, i) s * 120 + (k - 1) * 30 + h * 6 + i + 2L
rate_components <- function(state, time, pm) {
  y <- as.numeric(state); idxs <- function(s, k, h, age) idx0(s, k, h, age) - 1L
  incidence <- 0; fx_incidence <- 0
  for (age in 0:5) {
    lambda <- 0
    for (source_age in 0:5) {
      infectious <- 0; active <- 0
      for (s in c(0, 2)) for (k in 1:4) {
        infectious <- infectious + y[idxs(s, k, 1, source_age)] + y[idxs(s, k, 2, source_age)]
        active <- active + sum(y[sapply(0:4, function(h) idxs(s, k, h, source_age))])
      }
      if (active > 0) lambda <- lambda + pm$C_contact[age + 1, source_age + 1] * infectious / active
    }
    m <- pm$m_min + (pm$m_max - pm$m_min) / (1 + exp(-(time - pm$m_t0) / pm$m_tau))
    gam <- pm$q * m * lambda
    for (s in c(0, 2)) for (k in 1:4) {
      event_rate <- gam * (y[idxs(s, k, 0, age)] + y[idxs(s, k, 4, age)])
      incidence <- incidence + event_rate
      if (s == 2) fx_incidence <- fx_incidence + event_rate
    }
  }
  mortality <- 0; fx_mortality <- 0
  for (s in 0:3) for (k in 1:4) for (h in 0:4) for (age in 0:5) {
    disease_rate <- if (k == 3) pm$mu_DC * ifelse(h %in% c(3, 4), pm$psi_DC, 1) else
      if (k == 4) pm$mu_HCC * ifelse(h %in% c(3, 4), pm$psi_HCC, 1) else 0
    death_rate <- disease_rate * y[idxs(s, k, h, age)]
    mortality <- mortality + death_rate
    if (s >= 2) fx_mortality <- fx_mortality + death_rate
  }
  c(incidence = incidence, mortality = mortality,
    fx_incidence = fx_incidence, fx_mortality = fx_mortality)
}


strata <- c("D", "J", "F", "X")
idx_state <- function(s, h) as.vector(sapply(1:4, function(k) sapply(0:5, function(i) idx0(s,k,h,i)-1L)))
sum_state <- function(y, s, h=0:4) sum(y[unlist(lapply(h, function(x) idx_state(s,x)))])
state_from <- function(z) as.numeric(z[nrow(z),2:481])
set_tau <- function(pm, s, rate) { pm$tau_by_stratum[(s*4+1):(s*4+4)] <- rate; pm }
tau_j <- -log1p(-config$j_annual_start_probability)
tau_fx <- -log1p(-config$fx_annual_start_probability)
run_year <- function(y0, year, pm, cap=rep(Inf,4)) {
  d <- setup_env$data
  d$y0 <- y0; d$t_start <- 150+year-2016; d$t_end <- d$t_start+1
  d$treatment_capacity <- cap
  z <- run_sim(pm,d)
  if (abs(z[nrow(z),1]-d$t_end)>1e-6) stop("Annual simulation lost endpoint")
  z
}
# Entry-event proxy: antibody tests on antibody-positive prison entrants.
# This avoids recharging the full J stock each year; re-entries can recur.
j_positive_entries <- function(z, pm) {
  x <- vapply(seq_len(nrow(z)), function(i) {
    y <- z[i,2:481]; total <- 0
    for (age in 0:5) for (k in 1:4) for (h in c(1,2,4)) {
      total <- total + pm$lambda1[age+1]*y[idx0(0,k,h,age)-1] +
        pm$lambda3[age+1]*y[idx0(2,k,h,age)-1]
    }
    total
  }, numeric(1))
  sum((head(x,-1)+tail(x,-1))*diff(z[,1])/2)
}
# Screening quotas are one fixed year-start snapshot. No queue is tracked.
screen_plan <- function(y, s, quota) {
  n <- sum_state(y,s,c(0,1,2,4))
  ab <- sum_state(y,s,c(1,2,4)); rna <- sum_state(y,s,c(1,2)); chronic <- sum_state(y,s,2)
  tested <- if (is.infinite(quota)) n else min(quota,n)
  p_ab <- if (n>0) ab/n else 0
  p_rna <- if (n>0) rna/n else 0
  p_c <- if (n>0) chronic/n else 0
  list(eligible=n, chronic_stock=chronic, screening_quota=quota, screened=tested,
       antibody_tests=tested, rna_tests=tested*p_ab, rna_positive=tested*p_rna,
       chronic_identified=tested*p_c, treatment_capacity=tested*p_c,
       p_ab=p_ab,p_rna=p_rna,p_c=p_c)
}
stopifnot(screen_plan(rep(0,480),2,Inf)$treatment_capacity==0,
          screen_plan(rep(0,480),3,100)$antibody_tests==0)
process_row <- function(id, year, s, plan, z, j_ab=0) {
  events <- z[nrow(z),482:493]
  if (s==1) {
    ab <- j_ab; rna <- events[2]; identified <- events[2]
    quota <- NA_real_; eligible <- NA_real_; chronic <- NA_real_
    capacity <- Inf; screened <- ab; rnapos <- events[2]
  } else {
    ab <- plan$antibody_tests; rna <- plan$rna_tests
    identified <- plan$chronic_identified; quota <- plan$screening_quota
    eligible <- plan$eligible; chronic <- plan$chronic_stock
    capacity <- plan$treatment_capacity; screened <- plan$screened
    rnapos <- plan$rna_positive
  }
  starts <- events[s+1]
  if (is.finite(capacity) && starts>capacity+1e-6) stop("Treatment capacity exceeded")
  data.frame(scenario_id=id,year=year,stratum=strata[s+1],
    screening_quota=quota,screening_eligible=eligible,actual_screened=screened,
    antibody_tests=ab,rna_tests=rna,rna_positive=rnapos,
    chronic_stock_start=chronic,chronic_identified=identified,
    treatment_capacity=capacity,treatment_starts=starts,
    treatment_completions=events[4+s+1],treatment_cures=events[8+s+1],
    on_treatment_end=sum_state(state_from(z),s,3),
    unused_treatment_capacity=if (is.finite(capacity)) max(0,capacity-starts) else NA_real_,
    identified_chronic_start_ratio=if (identified>0) starts/identified else NA_real_,
    screened_chronic_stock_fraction=if (is.finite(chronic) && chronic>0) identified/chronic else NA_real_,
    screening_quota_use=if (is.finite(quota) && quota>0) screened/quota else NA_real_,
    stringsAsFactors=FALSE)
}
price <- unlist(config$unit_cost_sgd)
price["sof_vel_course"] <- config$sof_vel_usd_per_week*config$sof_vel_weeks*config$usd_to_sgd
stopifnot(config$sof_vel_usd_per_week*config$sof_vel_weeks==30924,
          abs(price["sof_vel_course"]-40405.2984)<1e-6)
add_cost <- function(d) {
  d$consultation_visits <- d$treatment_starts*config$consultations_per_start
  d$cost_antibody_sgd <- d$antibody_tests*price["antibody"]
  d$cost_rna_sgd <- d$rna_tests*price["rna"]
  d$cost_consultation_sgd <- d$consultation_visits*price["consultation"]
  d$cost_daa_sgd <- d$treatment_starts*price["sof_vel_course"]
  d$cost_total_sgd <- rowSums(d[,c("cost_antibody_sgd","cost_rna_sgd",
    "cost_consultation_sgd","cost_daa_sgd")])
  if (any(d$year>=config$policy_start_year)) {
    d$discount_factor <- (1+config$discount_rate)^-(d$year-config$discount_base_year)
    for (item in c("antibody","rna","consultation","daa","total"))
      d[[paste0("discounted_cost_",item,"_sgd")]] <-
        d[[paste0("cost_",item,"_sgd")]]*d$discount_factor
  }
  d
}

scenario_table <- expand.grid(screening_total=config$annual_screening_totals,
                              f_share=config$f_allocation_fractions,
                              KEEP.OUT.ATTRS=FALSE)
scenario_table$scenario_id <- sprintf("S%d_F%02d",scenario_table$screening_total,
                                      round(100*scenario_table$f_share))
scenario_table$scenario_type <- "finite"
scenario_table$f_screening_quota <- scenario_table$screening_total*scenario_table$f_share
scenario_table$x_screening_quota <- scenario_table$screening_total*(1-scenario_table$f_share)
scenario_table <- rbind(data.frame(screening_total=0,f_share=NA_real_,scenario_id="J_only",
                                   scenario_type="control",f_screening_quota=0,x_screening_quota=0),
                        scenario_table,
                        data.frame(screening_total=Inf,f_share=NA_real_,scenario_id="SInf",
                                   scenario_type="ideal",f_screening_quota=Inf,x_screening_quota=Inf))
scenario_table$description <- ifelse(scenario_table$scenario_type=="control","J background only",
  ifelse(scenario_table$scenario_type=="ideal","Unconstrained annual F/X screening",
    paste0("Annual screening ",scenario_table$screening_total,", F/X ",
           scenario_table$f_share,"/",1-scenario_table$f_share)))
scenario_table$scenario_label <- ifelse(scenario_table$scenario_type=="control","J only",
  ifelse(scenario_table$scenario_type=="ideal","No screening cap",
    sprintf("%s / F %s%%",format(scenario_table$screening_total,trim=TRUE),
      round(100*scenario_table$f_share))))
stopifnot(nrow(scenario_table)==17,sum(scenario_table$scenario_type=="finite")==15,
          !any(scenario_table$screening_total==800))

# One full policy analysis for one posterior parameter vector. Returns NULL if the
# 2015 snapshot is not at equilibrium.
evaluate_draw <- function(theta) {
base <- build_params(theta, fit$base_params)
snapshot_data <- setup_env$data; snapshot_data$t_start <- 0; snapshot_data$t_end <- 150
snapshot <- run_sim(base, snapshot_data)
if (!check_equilibrium(snapshot)$pass) return(NULL)
policy_params <- base
# Costed regimen is SOF/VEL alone. Use the model's RBV-free DC SVR branch.
policy_params$alpha_DC <- policy_params$alpha_DC_neg
policy_params$tau <- rep(1, 4)
policy_params$tau_min_age <- 1L; policy_params$tau_max_age <- 6L
policy_params$tau_by_stratum <- rep(0, 16)
pre_pm <- policy_params
j_pm <- set_tau(policy_params, 1, tau_j)
fx_pm <- set_tau(set_tau(j_pm, 2, tau_fx), 3, tau_fx)
# Annual incidence is an instantaneous rate at the calendar-year endpoint.
rate_row <- function(y, time, year, id) {
  r <- rate_components(y,time,j_pm)
  data.frame(scenario_id=id, year=year, incidence=r["incidence"], mortality=r["mortality"],
    fx_incidence=r["fx_incidence"], fx_mortality=r["fx_mortality"],
    viremic=sum(sapply(0:3,function(s) sum_state(y,s,1:3))),
    total=sum(y), fx_viremic=sum(sapply(2:3,function(s) sum_state(y,s,1:3))),
    fx_population=sum(sapply(2:3,function(s) sum_state(y,s))),
    stringsAsFactors=FALSE)
}
# Same calibrated 2015 state and same 2019--2025 J-treated history for every policy.
initial <- state_from(snapshot)
base_rates <- rate_row(initial,150,2015,"common_history")
history <- list(); hist_j_starts <- numeric(); y <- initial
for (year in 2016:2025) {
  before <- y
  z <- run_year(y,year,if (year<2019) pre_pm else j_pm)
  y <- state_from(z)
  history[[as.character(year)]] <- rate_row(y,151+year-2016,year,"common_history")
  if (year>=2019) hist_j_starts <- c(hist_j_starts,z[nrow(z),483])
}
common_2026_state <- y
history_trajectory <- do.call(rbind,c(list(base_rates),history))
stopifnot(length(hist_j_starts)==7,all(hist_j_starts>0))
run_policy <- function(row) {
  id <- row$scenario_id; y <- common_2026_state
  rows <- list(); trajectories <- list()
  for (year in config$policy_start_year:config$end_year) {
    f <- screen_plan(y,2,row$f_screening_quota)
    x <- screen_plan(y,3,row$x_screening_quota)
    cap <- c(Inf,Inf,f$treatment_capacity,x$treatment_capacity)
    pm <- if (row$scenario_type=="control") j_pm else fx_pm
    z <- run_year(y,year,pm,cap)
    if (year==2030) stopifnot(abs(z[nrow(z),1]-165)<1e-8)
    y <- state_from(z)
    jab <- j_positive_entries(z,pm)
    rows[[paste0(year,"J")]] <- process_row(id,year,1,NULL,z,jab)
    rows[[paste0(year,"F")]] <- process_row(id,year,2,f,z)
    rows[[paste0(year,"X")]] <- process_row(id,year,3,x,z)
    trajectories[[as.character(year)]] <- rate_row(y,151+year-2016,year,id)
  }
  list(process=add_cost(do.call(rbind,rows)),trajectory=do.call(rbind,trajectories),end_state=y)
}
# Counterfactual with no treatment in any stratum, 2015-2030 (shows the effect of J treatment from 2019).
no_tx <- list(base_rates); y <- initial
for (year in 2016:2030) {
  y <- state_from(run_year(y,year,pre_pm))
  no_tx[[as.character(year)]] <- rate_row(y,151+year-2016,year,"no_treatment")
}
no_tx <- do.call(rbind,no_tx); no_tx$scenario_id <- "no_treatment"; rownames(no_tx) <- NULL
for (metric in c("incidence","mortality","fx_incidence","fx_mortality"))
  no_tx[[paste0(metric,"_ratio_2015")]] <- no_tx[[metric]]/base_rates[[metric]]
runs <- lapply(seq_len(nrow(scenario_table)),function(i) run_policy(scenario_table[i,]))
process <- do.call(rbind,lapply(runs,`[[`,"process"))
trajectories <- do.call(rbind,lapply(seq_along(runs),function(i) {
  h <- history_trajectory; h$scenario_id <- scenario_table$scenario_id[i]
  rbind(h,runs[[i]]$trajectory)
}))
rownames(process)<-NULL;rownames(trajectories)<-NULL
stopifnot(all(is.finite(process$cost_total_sgd)),all(process$cost_total_sgd>=0),
          nrow(process)==nrow(scenario_table)*5*3,
          all(table(process$scenario_id,process$year,process$stratum)==1),
          all(abs(process$cost_total_sgd-rowSums(process[,c("cost_antibody_sgd",
            "cost_rna_sgd","cost_consultation_sgd","cost_daa_sgd")]))<1e-7),
          all(abs(process$cost_daa_sgd-process$treatment_starts*price["sof_vel_course"])<1e-7),
          all(abs(process$consultation_visits-process$treatment_starts)<1e-7),
          all(process$treatment_starts[process$stratum %in% c("F","X")] <=
                process$treatment_capacity[process$stratum %in% c("F","X")]+1e-6),
          all(process$treatment_starts[process$scenario_id=="J_only" &
                                         process$stratum %in% c("F","X")]<1e-8))
b <- base_rates[1,]
safe_ratio <- function(x,denominator) if (is.finite(denominator) && denominator>0)
  x/denominator else rep(NA_real_,length(x))
for (metric in c("incidence","mortality","fx_incidence","fx_mortality"))
  trajectories[[paste0(metric,"_ratio_2015")]] <- safe_ratio(trajectories[[metric]],as.numeric(b[[metric]]))
known_pass <- function(ratio,threshold) ifelse(is.finite(ratio),ratio<=threshold,NA)
trajectories$who_incidence_pass <- known_pass(trajectories$incidence_ratio_2015,.2)
trajectories$who_mortality_pass <- known_pass(trajectories$mortality_ratio_2015,.35)
trajectories$who_joint_pass <- ifelse(is.na(trajectories$who_incidence_pass) |
  is.na(trajectories$who_mortality_pass),NA,
  trajectories$who_incidence_pass & trajectories$who_mortality_pass)
trajectories$fx_incidence_pass <- known_pass(trajectories$fx_incidence_ratio_2015,.2)
trajectories$fx_mortality_pass <- known_pass(trajectories$fx_mortality_ratio_2015,.35)
trajectories$fx_joint_pass <- ifelse(is.na(trajectories$fx_incidence_pass) |
  is.na(trajectories$fx_mortality_pass),NA,
  trajectories$fx_incidence_pass & trajectories$fx_mortality_pass)
endpoint <- subset(trajectories,year==2030)
endpoint <- endpoint[match(scenario_table$scenario_id,endpoint$scenario_id),]
fxsum <- aggregate(cbind(actual_screened,antibody_tests,rna_tests,rna_positive,
  chronic_identified,treatment_capacity,treatment_starts,treatment_completions,
  treatment_cures,discounted_cost_total_sgd)~scenario_id,subset(process,stratum %in% c("F","X")),sum)
names(fxsum)[-1] <- paste0("fx_",names(fxsum)[-1])
endpoint <- merge(endpoint,fxsum,by="scenario_id",sort=FALSE)
jtotal <- aggregate(discounted_cost_total_sgd~scenario_id,subset(process,stratum=="J"),sum)
names(jtotal)[2] <- "j_discounted_cost_2026_2030_sgd"
endpoint <- merge(endpoint,jtotal,by="scenario_id",sort=FALSE)
endpoint <- endpoint[match(scenario_table$scenario_id,endpoint$scenario_id),]
endpoint$total_discounted_cost_2026_2030_sgd <- endpoint$j_discounted_cost_2026_2030_sgd+
  endpoint$fx_discounted_cost_total_sgd
control_cost <- endpoint$total_discounted_cost_2026_2030_sgd[endpoint$scenario_id=="J_only"]
endpoint$incremental_discounted_cost_vs_j_only_sgd <- endpoint$total_discounted_cost_2026_2030_sgd-control_cost
endpoint$incidence_reduction <- 1-endpoint$incidence_ratio_2015
endpoint$mortality_reduction <- 1-endpoint$mortality_ratio_2015
endpoint$fx_incidence_reduction <- 1-endpoint$fx_incidence_ratio_2015
endpoint$fx_mortality_reduction <- 1-endpoint$fx_mortality_ratio_2015
# Model outputs count expected courses, not distinct people.
treatment <- aggregate(cbind(treatment_starts,treatment_completions,treatment_cures)~scenario_id+stratum,
                       process,sum)
annual <- process
annual <- annual[order(annual$scenario_id,annual$year,annual$stratum),]
for (item in c("actual_screened","antibody_tests","rna_tests","rna_positive",
               "chronic_identified","treatment_starts","treatment_completions",
               "treatment_cures","consultation_visits","discounted_cost_antibody_sgd",
               "discounted_cost_rna_sgd","discounted_cost_consultation_sgd",
               "discounted_cost_daa_sgd","discounted_cost_total_sgd")) {
  annual[[paste0("cumulative_",item)]] <- ave(annual[[item]],
    interaction(annual$scenario_id,annual$stratum,drop=TRUE),FUN=cumsum)
}
for (id in scenario_table$scenario_id) {
  ix <- annual$scenario_id==id
  stopifnot(abs(sum(annual$discounted_cost_total_sgd[ix])-
    endpoint$total_discounted_cost_2026_2030_sgd[endpoint$scenario_id==id])<1e-5)
}
stopifnot(all(process$year>=2026 & process$year<=2030),
  all(is.finite(process$discounted_cost_total_sgd)),
  all(abs(process$discount_factor-(1.03)^-(process$year-2026))<1e-12))
item_cols <- c("discounted_cost_daa_sgd","discounted_cost_antibody_sgd",
  "discounted_cost_rna_sgd","discounted_cost_consultation_sgd")
item_totals <- aggregate(process[item_cols],list(scenario_id=process$scenario_id),sum)
stopifnot(all(abs(rowSums(item_totals[-1])-
  endpoint$total_discounted_cost_2026_2030_sgd[match(item_totals$scenario_id,endpoint$scenario_id)])<1e-5))
endpoint$total_treatment_starts <- as.numeric(tapply(process$treatment_starts,process$scenario_id,sum)[endpoint$scenario_id])
stratum_totals <- aggregate(cbind(discounted_cost_total_sgd,treatment_starts)~scenario_id+stratum,process,sum)
list(no_tx=no_tx,process=annual,trajectories=trajectories,endpoint=endpoint,treatment=treatment,
     item_totals=item_totals,stratum_totals=stratum_totals)
}

# ---- Aggregate over posterior draws: median and 95% UI (2.5th-97.5th percentile) ----
q_ui <- function(x,p) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_)
  if (all(x==x[1])) return(x[1]); unname(quantile(x,p,type=7,names=FALSE)) }
agg_draws <- function(dfs,keys) {
  d <- dfs[[1]]; n <- nrow(d)
  for (x in dfs) stopifnot(nrow(x)==n,identical(x[keys],d[keys]))
  out <- d[,keys,drop=FALSE]
  for (col in setdiff(names(d),keys)) {
    if (!is.numeric(d[[col]]) && !is.logical(d[[col]])) next
    m <- matrix(vapply(dfs,function(x) as.numeric(x[[col]]),numeric(n)),nrow=n)
    if (is.logical(d[[col]])) { out[[paste0(col,"_prob")]] <- rowMeans(m,na.rm=TRUE); next }
    s <- apply(m,1,function(v) c(median(v,na.rm=TRUE),q_ui(v,.025),q_ui(v,.975)))
    out[[col]] <- s[1,]; out[[paste0(col,"_lo")]] <- s[2,]; out[[paste0(col,"_hi")]] <- s[3,]
  }
  rownames(out) <- NULL; out
}
cat("Running",nrow(scenario_table),"policies x",nrow(draw_table),"posterior draws (",
    paste(table(draw_table$seed),collapse="+"),"per seed ) on",N_CORES,"cores...\n")
results <- parallel::mclapply(seq_len(nrow(draw_table)),function(i)
  tryCatch(evaluate_draw(as.numeric(draw_table[i,theta_cols])),error=function(e) e),
  mc.cores=N_CORES,mc.preschedule=TRUE)
errors <- vapply(results,function(r) inherits(r,"error"),logical(1))
if (any(errors)) stop("Draw failed: ",conditionMessage(results[[which(errors)[1]]]))
valid <- !vapply(results,is.null,logical(1))
cat("Valid draws:",sum(valid),"of",length(valid),"(equilibrium failures:",sum(!valid),")\n")
draw_seed <- draw_table$seed[valid]; results <- results[valid]
MIN_VALID_SHARE <- as.numeric(arg("--min-valid-share", 0.9))  # the NPE cross-validation policy run passes 0
stopifnot(length(results)>=MIN_VALID_SHARE*length(valid))
writeLines(c(paste("draws_requested", length(valid)), paste("draws_valid", sum(valid))), file.path(OUT,"draw_validity.txt"))
collect <- function(name) lapply(results,`[[`,name)
process <- agg_draws(collect("process"),c("scenario_id","year","stratum"))
trajectories <- agg_draws(collect("trajectories"),c("scenario_id","year"))
no_treatment <- agg_draws(collect("no_tx"),c("scenario_id","year"))
write.csv(no_treatment,file.path(OUT,"no_treatment_trajectory.csv"),row.names=FALSE)
endpoint_draws <- collect("endpoint")
endpoint <- agg_draws(endpoint_draws,"scenario_id")
# Draw-level 2030 endpoints for key outcomes (used for Monte Carlo error in the benchmark).
key_cols <- c("scenario_id","incidence_reduction","mortality_reduction","who_incidence_pass",
  "who_mortality_pass","who_joint_pass","total_discounted_cost_2026_2030_sgd","fx_treatment_starts")
write.csv(do.call(rbind,lapply(seq_along(endpoint_draws),function(i)
  cbind(draw=i,seed=draw_seed[i],endpoint_draws[[i]][,intersect(key_cols,names(endpoint_draws[[i]]))]))),
  file.path(OUT,"scenario_2030_draws.csv"),row.names=FALSE)
treatment <- agg_draws(collect("treatment"),c("scenario_id","stratum"))
item_totals <- agg_draws(collect("item_totals"),"scenario_id")
stratum_totals <- agg_draws(collect("stratum_totals"),c("scenario_id","stratum"))
# Pass flags use the median reductions shown in the figures; *_prob columns give the share of draws passing.
known_pass <- function(ratio,threshold) ifelse(is.finite(ratio),ratio<=threshold,NA)
join_pass <- function(a,b) ifelse(is.na(a)|is.na(b),NA,a & b)
for (d in c("trajectories","endpoint")) {
  x <- get(d)
  x$who_incidence_pass <- known_pass(x$incidence_ratio_2015,.2)
  x$who_mortality_pass <- known_pass(x$mortality_ratio_2015,.35)
  x$who_joint_pass <- join_pass(x$who_incidence_pass,x$who_mortality_pass)
  x$fx_incidence_pass <- known_pass(x$fx_incidence_ratio_2015,.2)
  x$fx_mortality_pass <- known_pass(x$fx_mortality_ratio_2015,.35)
  x$fx_joint_pass <- join_pass(x$fx_incidence_pass,x$fx_mortality_pass)
  assign(d,x)
}
endpoint <- merge(endpoint,scenario_table,by="scenario_id",sort=FALSE)
endpoint <- endpoint[match(scenario_table$scenario_id,endpoint$scenario_id),]
rownames(endpoint) <- NULL
endpoint$target_status <- ifelse(is.na(endpoint$who_joint_pass),"Unknown",
  ifelse(endpoint$who_joint_pass,"PASS","FAIL"))
endpoint$n_posterior_draws <- length(results)
endpoint$discount_rate <- config$discount_rate
endpoint$discount_base_year <- config$discount_base_year
endpoint$price_year <- config$price_year
endpoint$price_basis <- config$price_basis
# Per-file medians show whether the posterior files agree.
by_seed <- do.call(rbind,lapply(sort(unique(draw_seed)),function(k) {
  e <- agg_draws(endpoint_draws[draw_seed==k],"scenario_id")
  e <- e[,c("scenario_id","incidence_reduction","incidence_reduction_lo","incidence_reduction_hi",
            "mortality_reduction","mortality_reduction_lo","mortality_reduction_hi",
            "total_discounted_cost_2026_2030_sgd","total_discounted_cost_2026_2030_sgd_lo",
            "total_discounted_cost_2026_2030_sgd_hi","who_joint_pass_prob")]
  data.frame(seed=k,n_draws=sum(draw_seed==k),e)
}))
# Formatted "median (95% UI: lo-hi)" text for tables and the report.
ui_text <- function(col,scale=1,digits=1,prefix="",suffix="") {
  f <- function(v) paste0(prefix,formatC(v*scale,format="f",digits=digits,big.mark=","),suffix)
  sprintf("%s (95%% UI: %s\u2013%s)",f(endpoint[[col]]),f(endpoint[[paste0(col,"_lo")]]),f(endpoint[[paste0(col,"_hi")]]))
}
endpoint$incidence_reduction_text <- ui_text("incidence_reduction",100,1,suffix="%")
endpoint$mortality_reduction_text <- ui_text("mortality_reduction",100,1,suffix="%")
endpoint$total_cost_text <- ui_text("total_discounted_cost_2026_2030_sgd",1e-6,2,prefix="S$",suffix="m")
endpoint$incremental_cost_text <- ui_text("incremental_discounted_cost_vs_j_only_sgd",1e-6,2,prefix="S$",suffix="m")
endpoint$fx_starts_text <- ui_text("fx_treatment_starts",1,0)
endpoint$j_starts_text <- local({
  # J starts: use per-draw sums for the interval.
  jd <- agg_draws(lapply(collect("process"),function(p) {
    a <- aggregate(treatment_starts~scenario_id,subset(p,stratum=="J"),sum)
    a[match(scenario_table$scenario_id,a$scenario_id),]}),"scenario_id")
  f <- function(v) formatC(v,format="f",digits=0,big.mark=",")
  sprintf("%s (95%% UI: %s\u2013%s)",f(jd$treatment_starts),f(jd$treatment_starts_lo),f(jd$treatment_starts_hi))})
endpoint$joint_pass_probability <- endpoint$who_joint_pass_prob
annual <- process
stopifnot(nrow(process)==nrow(scenario_table)*5*3)
writeLines(paste0("Posterior draws: ",DRAWS_PER_SEED," random draws from each posterior file (valid: ",
  paste(paste0("file ",sort(unique(draw_seed)),": ",as.vector(table(draw_seed))),collapse=", "),
  "); sampling seed ",SAMPLE_SEED,". Medians and 95% UIs (2.5th-97.5th percentiles) across draws.\nSources: ",
  paste(posterior_files,collapse="; ")),file.path(OUT,"calibration_source.txt"))
write.csv(scenario_table,file.path(OUT,"scenario_definitions.csv"),row.names=FALSE)
write.csv(trajectories,file.path(OUT,"scenario_trajectories.csv"),row.names=FALSE)
write.csv(treatment,file.path(OUT,"scenario_treatment_summary.csv"),row.names=FALSE)
write.csv(annual,file.path(OUT,"scenario_annual_services_costs.csv"),row.names=FALSE)
write.csv(subset(endpoint,target_status=="PASS"),file.path(OUT,"target_scenario_treatment_summary.csv"),row.names=FALSE)
# Sensitivity on unit costs alone; model events and epidemiology are held fixed.
for (mult in config$cost_sensitivity_multiplier) for (suffix in c("","_lo","_hi")) {
  field <- paste0("cost_x",gsub("\\.","_",as.character(mult)),"_sgd",suffix)
  endpoint[[field]] <- endpoint[[paste0("total_discounted_cost_2026_2030_sgd",suffix)]]*mult
}
write.csv(endpoint,file.path(OUT,"scenario_2030_summary.csv"),row.names=FALSE)
if (requireNamespace("writexl",quietly=TRUE)) writexl::write_xlsx(
  list(Definitions=scenario_table,Endpoint=endpoint,BySeed=by_seed,Annual=annual,Treatment=treatment),
  file.path(OUT,"scenario_results.xlsx"))
# Medians of per-draw totals, for the cost figures in scripts/figures.py.
write.csv(stratum_totals,file.path(OUT,"scenario_stratum_totals.csv"),row.names=FALSE)
write.csv(item_totals,file.path(OUT,"scenario_item_totals.csv"),row.names=FALSE)

if (identical(arg("--no-figures", "false"), "true")) {
  cat("Wrote", nrow(scenario_table), "policies to", OUT, "(figures skipped)\n"); quit(save = "no")
}
library(ggplot2); library(tidyr); library(patchwork)
# Reuse the repository's Nature-style conventions: classic theme, 10 pt base,
# 0.7 pt trajectory lines, viridis policy palette and 300 dpi PNG export.
nature_theme <- function() theme_classic(base_size=10)+
  theme(plot.title=element_text(face="bold"),strip.background=element_blank(),
        strip.text=element_text(face="bold"),legend.position="bottom")
stratum_colours <- c(F="#D55E00",J="#0072B2",X="#009E73")
item_colours <- c(`SOF/VEL`="#0072B2",`HCV antibody`="#D55E00",
                  `HCV RNA`="#009E73",Consultation="#CC79A7")
item_order <- names(item_colours)
finite_levels <- scenario_table$scenario_id[scenario_table$scenario_type=="finite"]
policy_order <- scenario_table$scenario_id
policy_labels <- setNames(scenario_table$scenario_label,scenario_table$scenario_id)
policy_labels[c("J_only","SInf")] <- c("J only [reference]","No cap [reference]")
policy_label_factor <- function(ids) factor(policy_labels[as.character(ids)],
  levels=rev(unname(policy_labels[policy_order])))
finite <- subset(process,stratum %in% c("F","X") & scenario_id %in% finite_levels)
finite$scenario_label <- policy_label_factor(finite$scenario_id)
finite$share_label <- factor(paste0("F ",round(100*scenario_table$f_share[
  match(finite$scenario_id,scenario_table$scenario_id)]),"% / X ",
  round(100*(1-scenario_table$f_share[match(finite$scenario_id,scenario_table$scenario_id)])),"%"),
  levels=c("F 25% / X 75%","F 50% / X 50%","F 75% / X 25%"))
finite$capacity_label <- factor(paste0("S = ",
  scenario_table$screening_total[match(finite$scenario_id,scenario_table$scenario_id)]),
  levels=paste0("S = ",config$annual_screening_totals))

# Match the legacy A/B layout: relative trajectories at left and 2030
# indicator-specific attainment bars at right. The underlying measures remain
# annualized event flows divided by their own 2015 values, not person-time rates.
tr <- trajectories
meta_ix <- match(tr$scenario_id,scenario_table$scenario_id)
tr$screening_line <- ifelse(tr$scenario_id=="J_only","J only",
  ifelse(tr$scenario_id=="SInf","No screening cap",
    paste0(scenario_table$screening_total[meta_ix],"/year")))
line_order <- c(paste0(config$annual_screening_totals,"/year"),"J only","No screening cap")
tr$screening_line <- factor(tr$screening_line,levels=line_order)
tr$allocation_line <- ifelse(scenario_table$scenario_type[meta_ix]=="finite",
  paste0("F ",round(100*scenario_table$f_share[meta_ix]),"%"),"Reference")
tr$allocation_line <- factor(tr$allocation_line,
  levels=c("F 25%","F 50%","F 75%","Reference"))
line_colours <- c(setNames(viridisLite::viridis(5,option="D",begin=.1,end=.85),
  paste0(config$annual_screening_totals,"/year")),`J only`="grey35",`No screening cap`="#000000")
tr_long <- pivot_longer(tr,cols=c(incidence_ratio_2015,mortality_ratio_2015),
  names_to="indicator",values_to="baseline_ratio")
tr_long$indicator <- factor(tr_long$indicator,
  levels=c("incidence_ratio_2015","mortality_ratio_2015"),
  labels=c("Incidence","Mortality"))
trajectory_limits <- data.frame(indicator=factor(c("Incidence","Mortality"),
  levels=levels(tr_long$indicator)),threshold=c(.2,.35))
common_lines <- subset(tr_long,scenario_id=="J_only" & year<=2025)
policy_lines <- subset(tr_long,year>=2025)
p_traj <- ggplot()+
  geom_hline(yintercept=1,linetype=3,colour="grey35",linewidth=.4)+
  geom_hline(data=trajectory_limits,aes(yintercept=threshold),
    linetype=2,colour="black",linewidth=.45)+
  geom_vline(xintercept=2019,linetype=3,colour="grey65",linewidth=.35)+
  geom_vline(xintercept=2026,linetype=2,colour="grey65",linewidth=.35)+
  geom_line(data=common_lines,aes(year,baseline_ratio,group=indicator),
    colour="grey25",linewidth=.8)+
  geom_line(data=policy_lines,aes(year,baseline_ratio,colour=screening_line,
    linetype=allocation_line,group=scenario_id),linewidth=.7)+
  facet_wrap(~indicator,nrow=1)+
  scale_colour_manual(values=line_colours,drop=FALSE)+
  scale_linetype_manual(values=c(`F 25%`="solid",`F 50%`="22",
    `F 75%`="13",Reference="solid"),drop=FALSE)+
  scale_x_continuous(breaks=c(2015,2019,2026,2030),limits=c(2015,2030))+
  scale_y_continuous(labels=scales::percent,limits=c(0,1.2),
    breaks=seq(0,1.2,.2))+
  labs(x="Calendar year",y="Annual events relative to 2015",
    colour="Annual F/X screening",linetype="F allocation",
    title="A  WHO impact indicators, 2015–2030",
    subtitle="Annualized events / 2015; dotted 2019: J starts; dashed 2026: F/X starts")+
  nature_theme()+theme(legend.position="bottom",legend.title.position="left")

# A single 2030 attainment data frame drives both the heatmap and pie filter.
attainment <- endpoint[,c("scenario_id","scenario_type","screening_total","f_share",
 "incidence_reduction","mortality_reduction","incidence_reduction_lo","incidence_reduction_hi",
 "mortality_reduction_lo","mortality_reduction_hi","who_incidence_pass",
 "who_mortality_pass","who_joint_pass","target_status")]
attainment$scenario_label <- policy_label_factor(attainment$scenario_id)
heat <- rbind(
  data.frame(attainment,indicator="Incidence",value=attainment$incidence_reduction,
             value_lo=attainment$incidence_reduction_lo,value_hi=attainment$incidence_reduction_hi,
             pass=attainment$who_incidence_pass),
  data.frame(attainment,indicator="Mortality",value=attainment$mortality_reduction,
             value_lo=attainment$mortality_reduction_lo,value_hi=attainment$mortality_reduction_hi,
             pass=attainment$who_mortality_pass),
  data.frame(attainment,indicator="Joint",value=pmin(attainment$incidence_reduction/.8,
             attainment$mortality_reduction/.65),value_lo=NA_real_,value_hi=NA_real_,
             pass=attainment$who_joint_pass))
heat$indicator <- factor(heat$indicator,levels=c("Incidence","Mortality","Joint"))
heat$status <- ifelse(is.na(heat$pass),"Unknown",ifelse(heat$pass,"PASS","FAIL"))
heat$cell_label <- ifelse(heat$indicator=="Joint",heat$status,
  ifelse(is.finite(heat$value),sprintf("%.1f%%",100*heat$value),"Unknown"))
p <- ggplot(heat,aes(indicator,scenario_label,fill=status))+
  geom_tile(colour="white",linewidth=.5)+
  geom_text(aes(label=cell_label),colour="white",size=3.3)+
  scale_fill_manual(values=c(PASS="#009E73",FAIL="#D55E00",Unknown="grey55"),drop=FALSE)+
  labs(x=NULL,y=NULL,fill="2030 target",title="2030 WHO-style target attainment",
    subtitle="Versus 2015: incidence reduction ≥80%, mortality reduction ≥65%; joint requires both")+
  nature_theme()+theme(axis.text.x=element_text(face="bold"),
    panel.grid=element_blank(),plot.margin=margin(7,20,7,7))
ggsave(file.path(OUT,"scenario_target_heatmap.png"),p,width=12,height=10,dpi=300)

target_bars <- subset(heat,indicator %in% c("Incidence","Mortality"))
target_bars$indicator <- droplevels(target_bars$indicator)
bar_limits <- data.frame(indicator=factor(c("Incidence","Mortality"),
  levels=levels(target_bars$indicator)),threshold=c(.8,.65))
p_bars <- ggplot(target_bars,aes(scenario_label,value,fill=status))+
  geom_col(width=.75)+
  geom_errorbar(aes(ymin=value_lo,ymax=value_hi),width=.35,linewidth=.35,colour="grey20")+
  geom_hline(data=bar_limits,aes(yintercept=threshold),
    linetype=2,colour="black",linewidth=.55)+
  coord_flip()+facet_wrap(~indicator,nrow=1)+
  scale_fill_manual(values=c(PASS="#009E73",FAIL="#D55E00",Unknown="grey55"),drop=FALSE)+
  scale_y_continuous(labels=scales::percent,limits=c(0,1),
    breaks=c(0,.25,.5,.75),expand=expansion(mult=c(0,.02)))+
  labs(x=NULL,y="Reduction versus 2015",fill="2030 target",
    title="B  2030 WHO target attainment",
    subtitle="Bars: median; whiskers: 95% UI. Dashed thresholds: incidence 80%; mortality 65%")+
  nature_theme()+theme(legend.position="bottom",legend.title.position="left")
p <- (p_traj+p_bars+plot_layout(widths=c(1.1,1.9),guides="collect")) &
  theme(legend.position="bottom")
ggsave(file.path(OUT,"scenario_target_trajectories.png"),p,width=22,height=10,dpi=300)

# Two-panel budget and course-count figure, with a single F/J/X legend.
stratum_totals$scenario_label <- policy_label_factor(stratum_totals$scenario_id)
stratum_totals$stratum <- factor(stratum_totals$stratum,levels=c("F","J","X"))
budget_totals <- data.frame(scenario_id=endpoint$scenario_id,
  discounted_cost_total_sgd=endpoint$total_discounted_cost_2026_2030_sgd)
budget_totals$scenario_label <- policy_label_factor(budget_totals$scenario_id)
budget_totals$cost_m <- budget_totals$discounted_cost_total_sgd/1e6
budget_totals$lo_m <- endpoint$total_discounted_cost_2026_2030_sgd_lo/1e6
budget_totals$hi_m <- endpoint$total_discounted_cost_2026_2030_sgd_hi/1e6
start_totals <- data.frame(scenario_id=endpoint$scenario_id,
  treatment_starts=endpoint$total_treatment_starts,
  lo=endpoint$total_treatment_starts_lo,hi=endpoint$total_treatment_starts_hi)
start_totals$scenario_label <- policy_label_factor(start_totals$scenario_id)
p_cost <- ggplot(stratum_totals,aes(scenario_label,discounted_cost_total_sgd/1e6,fill=stratum))+
  geom_col(width=.7)+coord_flip(clip="off")+
  geom_errorbar(data=budget_totals,aes(scenario_label,ymin=lo_m,ymax=hi_m),inherit.aes=FALSE,
    width=.35,linewidth=.35)+
  geom_text(data=budget_totals,aes(scenario_label,hi_m,
    label=sprintf("%.1f",cost_m)),inherit.aes=FALSE,hjust=-.15,size=2.7)+
  scale_fill_manual(values=stratum_colours,breaks=c("F","J","X"),drop=FALSE)+
  scale_y_continuous(expand=expansion(mult=c(0,.19)))+
  labs(x=NULL,y="Discounted cost (2026 SGD millions)",fill="Stratum",
       title="A  2026–2030 discounted costs",
       subtitle="3% annual discount to 2026; label = median total, whiskers = 95% UI")+nature_theme()
p_starts <- ggplot(stratum_totals,aes(scenario_label,treatment_starts,fill=stratum))+
  geom_col(width=.7)+coord_flip(clip="off")+
  geom_errorbar(data=start_totals,aes(scenario_label,ymin=lo,ymax=hi),inherit.aes=FALSE,
    width=.35,linewidth=.35)+
  geom_text(data=start_totals,aes(scenario_label,hi,
    label=scales::comma(round(treatment_starts))),inherit.aes=FALSE,
    hjust=-.1,size=2.7)+
  scale_fill_manual(values=stratum_colours,breaks=c("F","J","X"),drop=FALSE)+
  scale_y_continuous(expand=expansion(mult=c(0,.19)),labels=scales::comma)+
  labs(x=NULL,y="Treatment starts (undiscounted courses)",fill="Stratum",
       title="B  2026–2030 starts",
       subtitle="Service counts are not discounted")+nature_theme()
p <- (p_cost+p_starts+plot_layout(guides="collect")) &
  theme(legend.position="bottom",legend.title.position="left")
ggsave(file.path(OUT,"scenario_cumulative_costs.png"),p,width=16,height=10,dpi=300)

# Annual F/X figure: column = screening slots, solid = actual starts,
# dashed = year-start screening-implied treatment capacity. The grid labels
# identify the split; the two quantities are never placed on one axis.
start_long <- pivot_longer(finite,cols=c(treatment_capacity,treatment_starts),
  names_to="measure",values_to="courses")
start_long$measure <- factor(start_long$measure,
  levels=c("treatment_starts","treatment_capacity"),
  labels=c("Actual starts","Screen-derived treatment quota"))
p <- ggplot(start_long,aes(year,courses,colour=stratum,linetype=measure,
  group=interaction(stratum,measure)))+
  geom_ribbon(data=finite,aes(year,ymin=treatment_starts_lo,ymax=treatment_starts_hi,fill=stratum,
    group=stratum),inherit.aes=FALSE,alpha=.18)+
  scale_fill_manual(values=stratum_colours[c("F","X")],guide="none")+
  geom_line(linewidth=.7)+geom_point(
  data=subset(start_long,measure=="Actual starts"),shape=21,fill="white",size=1.6,stroke=.6)+
  facet_grid(share_label~capacity_label,scales="free_y")+
  scale_colour_manual(values=stratum_colours[c("F","X")],drop=FALSE)+
  scale_linetype_manual(values=c("Actual starts"="solid",
    "Screen-derived treatment quota"="dashed"))+
  scale_x_continuous(breaks=c(2026,2028,2030))+
  labs(x="Calendar year",y="Treatment starts or quota (courses)",colour="Stratum",
    linetype="Measure",title="Annual F/X starts under screening quotas",
    subtitle="Lines: medians; shaded band: 95% UI of actual starts")+
  nature_theme()+theme(strip.text.y=element_text(angle=0),legend.position="bottom")
ggsave(file.path(OUT,"scenario_treatment_starts.png"),p,width=18,height=10,dpi=300)

finite$cumulative_starts <- ave(finite$treatment_starts,
  interaction(finite$scenario_id,finite$stratum,drop=TRUE),FUN=cumsum)
p <- ggplot(finite,aes(year,cumulative_starts,colour=stratum))+geom_line(linewidth=.7)+
  facet_grid(share_label~capacity_label,scales="free_y")+
  scale_colour_manual(values=stratum_colours[c("F","X")],drop=FALSE)+
  scale_x_continuous(breaks=c(2026,2028,2030))+
  labs(x="Calendar year",y="Cumulative treatment starts (courses)",colour="Stratum",
       title="Cumulative F/X treatment starts, 2026–2030")+
  nature_theme()+theme(strip.text.y=element_text(angle=0))
ggsave(file.path(OUT,"scenario_cumulative_treatment.png"),p,width=18,height=10,dpi=300)

util <- subset(finite,is.finite(treatment_capacity)&treatment_capacity>0)
p <- ggplot(util,aes(year,treatment_starts/treatment_capacity,colour=stratum))+
  geom_line(linewidth=.7)+facet_grid(share_label~capacity_label)+
  scale_colour_manual(values=stratum_colours[c("F","X")],drop=FALSE)+
  scale_x_continuous(breaks=c(2026,2028,2030))+
  scale_y_continuous(labels=scales::percent,limits=c(0,1.01))+
  labs(x="Calendar year",y="Screen-derived treatment quota used",colour="Stratum",
       title="Annual F/X treatment quota utilization")+
  nature_theme()+theme(strip.text.y=element_text(angle=0))
ggsave(file.path(OUT,"scenario_capacity_utilization.png"),p,width=18,height=10,dpi=300)

# The cascade is aggregated across F/X; all 15 finite policies are displayed.
cascade <- aggregate(cbind(actual_screened,rna_tests,rna_positive,chronic_identified,
  treatment_starts)~scenario_id+year,finite,sum)
cascade$share_label <- factor(paste0("F ",round(100*scenario_table$f_share[
  match(cascade$scenario_id,scenario_table$scenario_id)]),"% / X ",
  round(100*(1-scenario_table$f_share[match(cascade$scenario_id,scenario_table$scenario_id)])),"%"),
  levels=levels(finite$share_label))
cascade$capacity_label <- factor(paste0("S = ",scenario_table$screening_total[
  match(cascade$scenario_id,scenario_table$scenario_id)]),levels=levels(finite$capacity_label))
cascade_long <- pivot_longer(cascade,cols=c(actual_screened,rna_tests,rna_positive,
  chronic_identified,treatment_starts),names_to="service",values_to="count")
p <- ggplot(cascade_long,aes(year,count,colour=service))+geom_line(linewidth=.7)+
  facet_grid(share_label~capacity_label,scales="free_y")+
  scale_colour_viridis_d(option="D",end=.9)+
  scale_x_continuous(breaks=c(2026,2028,2030))+
  labs(x="Calendar year",y="Expected F+X events",colour="Service",
       title="F/X screening-to-treatment cascade, all finite policies")+
  nature_theme()+theme(strip.text.y=element_text(angle=0))
ggsave(file.path(OUT,"scenario_care_cascade.png"),p,width=18,height=10,dpi=300)

# Annual cost chart uses discounted values for all 17 policies.
annual_items <- pivot_longer(process,cols=c(discounted_cost_antibody_sgd,
  discounted_cost_rna_sgd,discounted_cost_consultation_sgd,
  discounted_cost_daa_sgd),names_to="item",values_to="discounted_sgd")
annual_items$item <- factor(annual_items$item,
  levels=c("discounted_cost_daa_sgd","discounted_cost_antibody_sgd",
    "discounted_cost_rna_sgd","discounted_cost_consultation_sgd"),
  labels=item_order)
annual_items <- aggregate(discounted_sgd~scenario_id+year+item,annual_items,sum)
annual_items$scenario_label <- policy_label_factor(annual_items$scenario_id)
p <- ggplot(annual_items,aes(year,discounted_sgd/1e6,fill=item))+
  geom_col(width=.75)+facet_wrap(~scenario_label,ncol=5,scales="free_y")+
  scale_fill_manual(values=item_colours,drop=FALSE)+
  scale_x_continuous(breaks=c(2026,2028,2030))+
  labs(x="Calendar year",y="Discounted cost (2026 SGD millions)",fill="Cost item",
    title="Annual F/J/X cost components, 2026–2030",
    subtitle="3% annual discount to 2026")+
  nature_theme()
ggsave(file.path(OUT,"scenario_annual_costs.png"),p,width=18,height=12,dpi=300)

# Cost composition uses exactly the endpoint joint PASS status used above.
# Pie shares normalise the per-item medians; the title shows the median total.
item_totals <- item_totals[,c("scenario_id","discounted_cost_daa_sgd","discounted_cost_antibody_sgd",
  "discounted_cost_rna_sgd","discounted_cost_consultation_sgd")]
item_long <- pivot_longer(item_totals,cols=-scenario_id,names_to="item",values_to="discounted_sgd")
item_long$item <- factor(item_long$item,
  levels=c("discounted_cost_daa_sgd","discounted_cost_antibody_sgd",
    "discounted_cost_rna_sgd","discounted_cost_consultation_sgd"),labels=item_order)
pass_ids <- attainment$scenario_id[attainment$target_status=="PASS"]
finite_pass_ids <- intersect(finite_levels,pass_ids)
pie_ids <- c("J_only",policy_order[policy_order %in% pass_ids & policy_order!="J_only"])
stopifnot(setequal(setdiff(pie_ids,"J_only"),setdiff(pass_ids,"J_only")))
pie_title <- function(id) {
  row <- scenario_table[match(id,scenario_table$scenario_id),]
  name <- if (id=="J_only") "Status quo (J only)" else if (id=="SInf")
    "No screening cap" else sprintf("%s/year · F %s%%",row$screening_total,
      round(100*row$f_share))
  sprintf("%s\nS$%.2fm",name,endpoint$total_discounted_cost_2026_2030_sgd[match(id,endpoint$scenario_id)]/1e6)
}
pie_titles <- setNames(vapply(pie_ids,pie_title,character(1)),pie_ids)
slice_order <- c("HCV RNA","HCV antibody","Consultation","SOF/VEL")
short_names <- c(`HCV RNA`="RNA",`HCV antibody`="Ab",Consultation="Visit")
outside_y <- c(`HCV RNA`=-1.42,`HCV antibody`=-1.22,Consultation=-1.02)
pie_wedges <- list();pie_leaders <- list();pie_outside <- list();pie_inside <- list()
for (id in pie_ids) {
  d <- subset(item_long,scenario_id==id)
  d <- d[match(slice_order,as.character(d$item)),]
  total <- sum(d$discounted_sgd)
  fractions <- d$discounted_sgd/total
  stopifnot(is.finite(total),total>0,abs(sum(fractions)-1)<1e-12)
  angle <- -pi/2; cx <- -.28; radius <- .72
  for (j in seq_along(slice_order)) {
    item <- slice_order[j];fraction <- fractions[j]
    next_angle <- angle+2*pi*fraction
    theta <- seq(angle,next_angle,length.out=max(4L,ceiling(720*fraction)))
    pie_wedges[[length(pie_wedges)+1]] <- data.frame(scenario_id=id,
      panel=unname(pie_titles[id]),item=item,
      x=c(cx,cx+radius*cos(theta),cx),y=c(0,radius*sin(theta),0))
    mid <- (angle+next_angle)/2
    if (item=="SOF/VEL") {
      pie_inside[[length(pie_inside)+1]] <- data.frame(panel=unname(pie_titles[id]),
        x=cx-.15,y=.17,label=sprintf("%.1f%%",100*fraction))
    } else {
      pie_leaders[[length(pie_leaders)+1]] <- data.frame(panel=unname(pie_titles[id]),
        x=.51,y=outside_y[item],xend=cx+radius*cos(mid),
        yend=radius*sin(mid))
      pie_outside[[length(pie_outside)+1]] <- data.frame(panel=unname(pie_titles[id]),
        x=.57,y=outside_y[item],label=sprintf("%s %.2f%%",
          short_names[item],100*fraction))
    }
    angle <- next_angle
  }
}
pie_wedges <- do.call(rbind,pie_wedges)
pie_leaders <- do.call(rbind,pie_leaders)
pie_outside <- do.call(rbind,pie_outside)
pie_inside <- do.call(rbind,pie_inside)
for (data_name in c("pie_wedges","pie_leaders","pie_outside","pie_inside"))
  assign(data_name,transform(get(data_name),panel=factor(panel,levels=unname(pie_titles))))
p <- ggplot(pie_wedges,aes(x,y,fill=item,group=interaction(scenario_id,item)))+
  geom_polygon(colour="white",linewidth=.15)+
  geom_segment(data=pie_leaders,aes(x=x,y=y,xend=xend,yend=yend),
    inherit.aes=FALSE,colour="grey30",linewidth=.35,
    arrow=grid::arrow(length=grid::unit(.06,"in"),ends="last"))+
  geom_text(data=pie_outside,aes(x=x,y=y,label=label),inherit.aes=FALSE,
    hjust=0,size=3)+
  geom_text(data=pie_inside,aes(x=x,y=y,label=label),inherit.aes=FALSE,
    colour="white",fontface="bold",size=3.5)+
  facet_wrap(~panel,nrow=1)+coord_fixed(xlim=c(-1.05,1.57),ylim=c(-1.53,.82),clip="off")+
  scale_fill_manual(values=item_colours,breaks=item_order,drop=FALSE)+
  labs(fill="Cost item",title="2026–2030 discounted cost composition",
    subtitle="Status quo reference and scenarios meeting both 2030 targets · 2026 SGD · 3% annual discount")+
  theme_void(base_size=10)+theme(plot.title=element_text(face="bold"),
    strip.text=element_text(face="bold",size=10),strip.background=element_blank(),
    legend.position="bottom",legend.title.position="left",
    panel.spacing=grid::unit(.1,"lines"),
    plot.background=element_rect(fill="white",colour=NA),
    panel.background=element_rect(fill="white",colour=NA),
    legend.background=element_rect(fill="white",colour=NA))
ggsave(file.path(OUT,"scenario_cost_composition_pies.png"),p,width=20,height=5.8,dpi=300,
  device=ragg::agg_png,bg="white")
heat_joint_ids <- as.character(heat$scenario_id[heat$indicator=="Joint" & heat$status=="PASS"])
stopifnot(setequal(finite_pass_ids,intersect(finite_levels,heat_joint_ids)))
cat("Wrote",nrow(scenario_table),"policies to",OUT,"\n")

