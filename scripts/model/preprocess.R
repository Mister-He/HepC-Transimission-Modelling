# =============================================================================
# preprocess.R — cleans Singapore Prison Service records (data/sps.csv) and fits
# release and re-arrest rates with Stan; run from the repository root.
# =============================================================================

library(dplyr)
library(lubridate)

# Clear the environment
rm(list = ls())

gc()

# Read SPS data
sps = read.csv('data/sps.csv')

# ================ Data preparation ================#
# Data cleaning and preprocessing
data = sps %>% 
    # Rename the columns for better readability
    rename(
        "Screening_result" = "Result.from.screening.exercise..Dec.2014.Feb2016.",
        "HepC_database" = "Hep.C.Database..as.of.Dec.2017.",
        "Date_diagnosed" = "Date.Diagnosed..recorded.in.Hep.C.Database"
    ) %>%
    mutate(
        Screening_result = case_when(
            Screening_result == "Borderline Reactive" | Screening_result == "Reactive" ~ "Positive",
            Screening_result == "Non-Reactive" ~ "Negative",
            TRUE ~ "Not Known"
        ),
        HepC_database = case_when(
            HepC_database == " Hep C Positive" ~ "Positive",
            TRUE ~ "Negative"
        )
    ) %>%
    # Focus on those who were arrested after 2014
    mutate(
        ADM_DT = as.Date(ADM_DT, format = "%Y/%m/%d"),
        SENTENCE_DT = as.Date(SENTENCE_DT, format = "%Y/%m/%d"),
        PERM_REL_DATE = as.Date(PERM_REL_DATE, format = "%Y/%m/%d")
    )

# Keep everyone ever held for drug consumption. The modelled PWID population is the screened
# subset (known HCV screening result; unscreened people are assumed not to be PWID): it gives
# the 2015 prison counts and first offenders below. Release and re-arrest rates use all
# consumption records, because screening (Dec 2014-Feb 2016) selects people in prison then,
# whose earlier releases all ended in re-arrest (re-arrest share 0.99 in that subset).
data = data %>% 
  left_join(
    data %>% 
      group_by(IDENTIFER) %>% 
      reframe(
        counts = sum(CONSUMPTION)
      ) %>% 
      ungroup(),
    by = "IDENTIFER"
  ) %>% 
  filter(counts > 0) %>%
  dplyr::select(-counts)

# Merge consecutive incarcerations for the same individual
clean_flag = T
while(clean_flag) {
    temp_data = data %>%
        arrange(IDENTIFER, ADM_DT) %>%
        group_by(IDENTIFER) %>%
        mutate(
            gap = as.numeric(ADM_DT - lag(PERM_REL_DATE)),
            new_incarceration = ifelse(is.na(gap) | gap > 0, 1, 0),
            incarceration_id = cumsum(new_incarceration)
        )
    
    # If there's no negative gap, 
    # it means all consecutive incarcerations have been merged, and we can stop the loop
    if (nrow(filter(temp_data, gap < 0)) == 0) clean_flag = F
    
    data = temp_data %>% 
       group_by(IDENTIFER, incarceration_id) %>%
           reframe(
               IDENTIFER = first(IDENTIFER),
               DRC_IND = max(DRC_IND),
               NATIONALITY = first(NATIONALITY),
               ADM_DT = min(ADM_DT),
               GENDER_CODE = first(GENDER_CODE),
               RACE = first(RACE),
               YEAR_OF_BIRTH = first(YEAR_OF_BIRTH),
               PERM_REL_DATE = max(PERM_REL_DATE),
               Screening_result = first(Screening_result),
               HepC_database = first(HepC_database),
               Date_diagnosed = first(Date_diagnosed)
           ) %>%
           ungroup() %>%
           select(-incarceration_id)
}

# Generate indicator variable for having HepC or not within each incarceration
data = data %>%
    mutate(HEPC_IND = ifelse(
        Screening_result == "Positive", 
        "Positive",
        # Assuming all time positive for 'PRE EXISTING' cases
        ifelse(
            Date_diagnosed == 'PRE EXISTING',
            "Positive",
            ifelse(
                Date_diagnosed %in% c('N.A.', 'No Record'),
                HepC_database,
                # Assume already HCV positive for those admission rate within 3 months
                ifelse(
                    abs(as.numeric(ADM_DT - as.Date(Date_diagnosed, format = "%m/%d/%Y"))) <= 90,
                    "Positive",
                    "Negative"
                )
            )
        )
    ))

# Assign start, end date and time in jail for each incarceration
## Add 1 day to include the release date to avoid 0 day incarceration
data = data %>% 
    mutate(
        start = ADM_DT,
        end = PERM_REL_DATE,
        t_injail = 1 + ifelse(!is.na(end), as.numeric(end-start), as.numeric(as.Date("2017-12-31") - start))
    )

# Assign release or not as event at 2017-12-31
data = data %>% mutate(event = ifelse(is.na(end), 0, 1))

# Assign age group index for each incarceration when admitted
data = data %>% 
    mutate(adm_age = year(start) - YEAR_OF_BIRTH,
           adm_age_grp = case_when(
               adm_age < 20 ~ "<20",
               adm_age >= 20 & adm_age < 30 ~ "20-29",
               adm_age >= 30 & adm_age < 40 ~ "30-39",
               adm_age >= 40 & adm_age < 50 ~ "40-49",
               adm_age >= 50 & adm_age < 60 ~ "50-59",
               adm_age >= 60 ~ "60+",
               T ~ "OUTOFRANGE"
           ))

# Ensure censored record has to show up at last
# Detect are there any record being NA under 'end'
# If there are, it means the record is censored, and we need to make sure it shows up at last for each individual
data = lapply(
    unique(data$IDENTIFER),
    function(id){
        temp = data %>% filter(IDENTIFER == id) %>% arrange(start)
        if(any(is.na(temp$end))){
            loc = which(is.na(temp$end))
            temp = temp[c(1:loc), ]
        }
        temp
    }
) %>% do.call('rbind', .)

# Save the cleaned data for future analysis
saveRDS(data, 'output/cleaned_data.rds')

# Total population counts for calibration targets (use 2015 snapshot)
time_end = as.Date("2015-12-31")
data %>%
    filter(
      start <= time_end, ,
      is.na(end) | end >= time_end,
      Screening_result != "Not Known"
    ) %>%
    mutate(
        c_age = year(time_end) - YEAR_OF_BIRTH,
        c_age_grp = case_when(
            c_age < 20 ~ "<20",
            c_age >= 20 & c_age < 30 ~ "20-29",
            c_age >= 30 & c_age < 40 ~ "30-39",
            c_age >= 40 & c_age < 50 ~ "40-49",
            c_age >= 50 & c_age < 60 ~ "50-59",
            c_age >= 60 ~ "60+",
            T ~ "OUTOFRANGE"
        )
    ) %>%
    group_by(c_age_grp) %>%
    reframe(total = n()) %>% 
    pull

# HCV prevalence target for calibration
data %>%
    filter(
      start <= time_end, ,
      is.na(end) | end >= time_end,
      Screening_result != "Not Known"
    ) %>%
    group_by(IDENTIFER) %>%
    slice(n()) %>%
    ungroup() %>%
    group_by(adm_age_grp) %>%
    reframe(positive = sum(HEPC_IND == "Positive"), total = n()) %>%
    mutate(percent = positive / total) %>%
    pull(percent)


# ================ Release rate estimation ================#
library(dplyr)
library(lubridate)
library(cmdstanr)

# Clean the environment
rm(list = ls())
gc()

# Load the cleaned data
data = readRDS('output/cleaned_data.rds')
data$t_injail = data$t_injail / 365.25 # convert to years

## Lexis splitting
## applied to each incarceration due to age grp specific hazard
split_by_age = function(id, t_tot, age_s, obs_yn, n_age, ageint) {
    out <- c()
    idx <- sum(ageint <= age_s)
    while (t_tot > 0 && idx <= n_age) {
        if (idx == n_age) {
            t0 <- t_tot
            event <- obs_yn # observed
        } else {
            t0 <- min(t_tot, ageint[idx + 1] - max(age_s, ageint[idx]))
            event <- ifelse(abs(t_tot - t0) < 1e-3, obs_yn, 0)
        }
        out <- rbind(out, data.frame(id = id, age = idx, t = t0, event = event))
        idx <- idx + 1
        t_tot <- t_tot - t0
    }
    out
}

## Age group breaks
ageint = c(0, 2:6)*10
n_age = length(ageint)

release_dat = lapply(1:nrow(data), function(i) {
    if (data$t_injail[i] > 0) {
        split_by_age(id = data$IDENTIFER[i], 
                     t_tot = data$t_injail[i], 
                     age_s = data$adm_age[i], 
                     obs_yn = data$event[i], 
                     n_age = n_age, 
                     ageint = ageint)
    }
}) %>% do.call("rbind", .)

## HMC input parameters
set.seed(114514)
N_BURNIN = 1e3
N_ITER = 5e3
N_CHAINS = 4
N_AGE_GRP = length(unique(data$adm_age_grp))

data_list = list(
    n_age = n_age,
    N_seg = nrow(release_dat),
    t = release_dat$t,
    age = release_dat$age,
    event = release_dat$event
)

model = cmdstan_model('scripts/model/release_rate.stan')

## HMC
release_fit = model$sample(
    data = data_list,
    chains = N_CHAINS,
    iter_warmup = N_BURNIN,
    iter_sampling = N_ITER - N_BURNIN,
    parallel_chains = N_CHAINS,
    save_cmdstan_config = T
)

## Release rate estimation results
### Including Rhat and ESS for convergence diagnostics
release_fit$summary()
saveRDS(release_fit, 'output/release_rate_fit.rds')

# =============== Rearrest rate estimation ================#
library(dplyr)
library(lubridate)
library(cmdstanr)

# Clear the environment
rm(list = ls())
gc()

# Load the cleaned data
data = readRDS('output/cleaned_data.rds')
data$t_injail = data$t_injail / 365.25 # convert to years

# Derive rearrest at risk intervals for further HMC estimation
rearrest_dat = lapply(
    unique(data$IDENTIFER), 
    function(id){
        temp = data %>% filter(IDENTIFER == id) %>% arrange(start)
        n.record = nrow(temp)
        
        # If multiple records
        if(n.record > 1){
            out1 = lapply(
                2:n.record, 
                function(i){
                    delta = as.numeric(temp$start[i] - temp$end[i-1]) / 365.25
                    age = temp$adm_age[i-1] + temp$t_injail[i-1]
                    data.frame(
                        id = id,
                        idx = i - 1,
                        t = delta,
                        age_rel = age,
                        censor = 0,
                        rearrest = 1
                    )
                } 
            ) %>% do.call('rbind', .)
        }

        # Deal with last/first record
        temp = temp %>% slice(n.record)
        last_obs = F
        if(temp$end < as.Date('2017-12-31') & !is.na(temp$end)){
            duration = as.numeric(as.Date('2017-12-31') - temp$end) / 365.25
            out0 = data.frame(
                id = id,
                idx = n.record,
                t = duration,
                age_rel = temp$adm_age + temp$t_injail,
                censor = 1,
                rearrest = 0
            )
            last_obs = T
        }

        # Aggregation
        if(n.record > 1 & last_obs) rbind(out1, out0)
        else if(n.record > 1) out1
        else if(last_obs) out0
    }
) %>% do.call('rbind', .)

# Lexis splitting
split_by_age <- function(id, t_tot, age_s, obs_yn, n_age, ageint) {
    out <- c()
    idx <- sum(ageint <= age_s)
    while (t_tot > 0 && idx <= n_age) {
        if (idx == n_age) {
            t0 <- t_tot
            event <- obs_yn # observed
        } else {
            t0 <- min(t_tot, ageint[idx + 1] - max(age_s, ageint[idx]))
            event <- ifelse(abs(t_tot - t0) < 1e-3, obs_yn, 0)
        }
        out <- rbind(out, data.frame(id = id, age = idx, t = t0, event = event))
        idx <- idx + 1
        t_tot <- t_tot - t0
    }
    out
}

## Age group breaks
ageint <- c(0, 2:6) * 10
n_age <- length(ageint)

rearrest_dat_lexis = lapply(
    1:nrow(rearrest_dat),
    function(i){
        temp0 = rearrest_dat %>% slice(i)
        temp1 = split_by_age(
            id = temp0$id,
            t_tot = temp0$t,
            age_s = temp0$age_rel,
            obs_yn = temp0$rearrest,
            n_age = n_age,
            ageint = ageint
        )
        expand_age_grp = rep(0, n_age)
        expand_age_grp[round(temp1$age)] = temp1$t
        output = c(
            range(temp1$age),
            expand_age_grp,
            temp0$rearrest
        )
        return(output)
    }
) %>% do.call('rbind', .)

# Other data for equailibrium calibration: the screened PWID population (see the filter above)
screened = data %>% group_by(IDENTIFER) %>% filter(any(Screening_result != "Not Known")) %>% ungroup()
J_pop = screened %>% 
    filter(start <= as.Date('2015-12-31'), end >= as.Date('2015-12-31') | is.na(end)) %>% 
    group_by(IDENTIFER) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(current_age = 2015 - YEAR_OF_BIRTH) %>%
    mutate(agegrp = cut(current_age, breaks = c(ageint, Inf), labels = seq_along(ageint), right = FALSE)) %>%
    group_by(agegrp) %>%
    reframe(n = n()) %>%
    pull(n)

# Label the first incarceration as first offender
screened <- screened %>%
    group_by(IDENTIFER) %>%
    arrange(start) %>%
    mutate(first_time = ifelse(row_number() == 1, 1, 0)) %>%
    ungroup()

# Overall first-arrest age grp distribution
f_all_dist = screened %>% 
    filter(first_time == 1) %>%
    mutate(start = year(start) - YEAR_OF_BIRTH) %>% 
    mutate(agegrp = cut(start, breaks = c(ageint, Inf), labels = seq_along(ageint), right = FALSE)) %>%
    group_by(agegrp) %>%
    reframe(n = n()) %>%
    pull(n)

# first-arrest age grp distribution by year (2014-2016)
f_yr_dist = lapply(
    2014:2016,
    function(yr){
        screened %>% 
            filter(first_time == 1, year(start) == yr) %>%
            mutate(start = year(start) - YEAR_OF_BIRTH) %>% 
            mutate(agegrp = cut(start, breaks = c(ageint, Inf), labels = seq_along(ageint), right = FALSE)) %>%
            group_by(agegrp) %>%
            reframe(n = n()) %>%
            pull(n)
    }
) %>% do.call('rbind', .) %>% as.matrix

f_obs = nrow(f_yr_dist)

# Death rate and ageing rate: the ODE's values, so this steady-state balance matches sim.cpp
# (background mortality mu * omega; 1/10 per year between 10-year bands).
setup_env = new.env(); old_wd = setwd('scripts/model'); sys.source('setup.R', setup_env); setwd(old_wd)
mu = setup_env$params$mu * setup_env$params$omega
a_age = 0.1

# Fitted release rate from previous step
release_fit = readRDS('output/release_rate_fit.rds')
lambda1 = exp(release_fit$summary('log_lambda')$mean)

## HMC input parameters
set.seed(114514)
N_BURNIN <- 1e3
N_ITER <- 5e3
N_CHAINS <- 4
N_AGE_GRP <- length(unique(data$adm_age_grp))

data_list <- list(
    n_age = n_age,
    release = lambda1,
    F0 = f_yr_dist,
    f_obs = f_obs,
    J = J_pop,
    mu = mu,
    a = a_age,
    N = nrow(rearrest_dat_lexis),
    t = rearrest_dat_lexis[ , 2 + 1:n_age],
    age_r = rearrest_dat_lexis[ , 1:2],
    event = rearrest_dat_lexis[ , 2 + n_age + 1]
)

model <- cmdstan_model("scripts/model/rearrest_rate.stan")

## HMC
rearrest_fit <- model$sample(
    data = data_list,
    chains = N_CHAINS,
    iter_warmup = N_BURNIN,
    iter_sampling = N_ITER - N_BURNIN,
    parallel_chains = N_CHAINS,
    save_cmdstan_config = T
)

## Rearrest rate estimation results
### Including Rhat and ESS for convergence diagnostics
# rearrest_fit$summary()
list(
    rearrest_fit$summary(variables = c(
        "pi_raw", "log_lambda", "f_pop"
    )),
    rearrest_fit$draws(variables = c(
        "pi_raw", "log_lambda", "f_pop"
    ))
) %>% 
saveRDS("output/rearrest_rate_fit.rds")

rearrest_fit$summary(variables = c(
    "pi", "lambda", "f_pop"
))

# # Draw traceplots using bayesplot
# library(bayesplot)
# label_parsed <- function(labels) {
#     labels$variable <- gsub("pi_raw", "pi[raw]", labels$variable)
#     labels$variable <- gsub("log_lambda", "log(lambda)", labels$variable)
#     labels$variable <- gsub("f_pop", "f[pop]", labels$variable)
#     return(labels)
# }
# # Traceplots for pi_raw
# mcmc_trace(rearrest_fit$draws(variables = c(
#     "pi_raw", "log_lambda", "f_pop"
# )))

# Copy the fitted rates into setup.R (replaces manual copying):
#   Rscript scripts/model/update_rate_params.R
