# =============================================================================
# preprocess.R — cleans Singapore Prison Service records (data/sps.csv) and fits
# release and re-arrest rates with Stan; run from the repository root.
# =============================================================================

library(dplyr)
library(lubridate)
library(cmdstanr)

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

# Keep everyone ever held for drug consumption. Release and re-arrest rates use all
# consumption records: screening (Dec 2014-Feb 2016) selects people in prison then, whose
# earlier releases all ended in re-arrest (re-arrest share 0.99 in that subset). The records
# cannot identify true first offenders, so no first-arrest quantity is estimated here.
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

# Keep each person's records up to and including the first censored one (end = NA), in
# order of first appearance and admission date.
data = data %>%
    mutate(.ord = match(IDENTIFER, unique(IDENTIFER))) %>%
    arrange(.ord, start) %>%
    group_by(IDENTIFER) %>%
    filter(cumsum(lag(is.na(end), default = FALSE)) == 0) %>%
    ungroup() %>%
    select(-.ord) %>%
    as.data.frame()

# Save the cleaned data for future analysis
saveRDS(data, 'output/cleaned_data.rds')

# 2015 prison snapshot of the screened PWID, aggregated by age group (not calibration targets;
# the manuscript reports the total against the modelled prison population).
time_end = as.Date("2015-12-31")
in_prison_2015 = data %>%
    filter(start <= time_end, is.na(end) | end >= time_end, Screening_result != "Not Known")
counts_2015 = in_prison_2015 %>%
    mutate(age_group = cut(year(time_end) - YEAR_OF_BIRTH, c(0, 2:6 * 10, Inf), right = FALSE,
                           labels = c("<20", "20-29", "30-39", "40-49", "50-59", "60+"))) %>%
    count(age_group, name = "n_in_prison")
# Seroprevalence by age at admission, from each person's last record in the snapshot.
sero_2015 = in_prison_2015 %>%
    group_by(IDENTIFER) %>% slice(n()) %>% ungroup() %>%
    group_by(age_group = adm_age_grp) %>%
    reframe(positive = sum(HEPC_IND == "Positive"), screened = n())
write.csv(full_join(counts_2015, sero_2015, by = "age_group"), "output/prison_2015_counts.csv", row.names = FALSE)


# ================ Release rate estimation ================#
data$t_injail = data$t_injail / 365.25 # convert to years

## Lexis splitting: cut each spell (length t_tot from age age_s) at the 10-year age bands.
## Bands are walked for all spells at once; the event (obs) goes on the piece that ends within
## 1e-3 years of the spell end, and always on the open 60+ band.
ageint = c(0, 2:6)*10
n_age = length(ageint)
split_by_age = function(t_tot, age_s, obs) {
    idx <- findInterval(age_s, ageint); pieces <- list()
    repeat {
        act <- which(t_tot > 0 & idx <= n_age)
        if (!length(act)) break
        last <- idx[act] == n_age
        t0 <- ifelse(last, t_tot[act],
                     pmin(t_tot[act], ageint[pmin(idx[act] + 1, n_age)] - pmax(age_s[act], ageint[idx[act]])))
        event <- ifelse(last, obs[act], ifelse(abs(t_tot[act] - t0) < 1e-3, obs[act], 0))
        pieces[[length(pieces) + 1]] <- data.frame(row = act, step = length(pieces) + 1, age = idx[act],
                                                   t = t0, event = event)
        idx[act] <- idx[act] + 1; t_tot[act] <- t_tot[act] - t0
    }
    out <- do.call(rbind, pieces)
    out[order(out$row, out$step), ]
}

release_dat = split_by_age(data$t_injail, data$adm_age, data$event)
release_dat = data.frame(id = data$IDENTIFER[release_dat$row], age = release_dat$age,
                         t = release_dat$t, event = release_dat$event)

## HMC input parameters
set.seed(114514)
N_BURNIN = 1e3
N_ITER = 5e3
N_CHAINS = 4

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

## Release rate estimation results (summary and draws, as for the re-arrest fit)
release_fit$summary()
list(release_fit$summary(variables = "log_lambda"), release_fit$draws(variables = "log_lambda")) %>%
saveRDS('output/release_rate_fit.rds')

# =============== Rearrest rate estimation ================#
# At-risk intervals after each release: a gap to the next admission (re-arrest), or for the
# last record, the time from release to 2017-12-31 (censored).
cutoff = as.Date('2017-12-31')
by_person = data %>% group_by(IDENTIFER) %>% mutate(i = row_number(), n.record = n()) %>% ungroup()
gaps = by_person %>%
    group_by(IDENTIFER) %>%
    mutate(t = as.numeric(start - lag(end)) / 365.25, age_rel = lag(adm_age) + lag(t_injail)) %>%
    ungroup() %>% filter(i > 1) %>%
    transmute(id = IDENTIFER, idx = i - 1, t, age_rel, censor = 0, rearrest = 1)
tails = by_person %>%
    filter(i == n.record, !is.na(end), end < cutoff) %>%
    transmute(id = IDENTIFER, idx = n.record, t = as.numeric(cutoff - end) / 365.25,
              age_rel = adm_age + t_injail, censor = 1, rearrest = 0)
rearrest_dat = bind_rows(gaps, tails) %>%
    arrange(match(id, unique(data$IDENTIFER)), idx) %>%
    as.data.frame()

pieces = split_by_age(rearrest_dat$t, rearrest_dat$age_rel, rearrest_dat$rearrest)
stopifnot(all(seq_len(nrow(rearrest_dat)) %in% pieces$row))  # every interval has positive length
t_mat = matrix(0, nrow(rearrest_dat), n_age); t_mat[cbind(pieces$row, round(pieces$age))] = pieces$t
rearrest_dat_lexis = cbind(tapply(pieces$age, pieces$row, min), tapply(pieces$age, pieces$row, max),
                           t_mat, rearrest_dat$rearrest)

## HMC input parameters
set.seed(114514)
N_BURNIN <- 1e3
N_ITER <- 5e3
N_CHAINS <- 4

data_list <- list(
    n_age = n_age,
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
        "pi_raw", "log_lambda"
    )),
    rearrest_fit$draws(variables = c(
        "pi_raw", "log_lambda"
    ))
) %>% 
saveRDS("output/rearrest_rate_fit.rds")

rearrest_fit$summary(variables = c(
    "pi", "lambda"
))

# # Draw traceplots using bayesplot
# library(bayesplot)
# label_parsed <- function(labels) {
#     labels$variable <- gsub("pi_raw", "pi[raw]", labels$variable)
#     labels$variable <- gsub("log_lambda", "log(lambda)", labels$variable)
#     return(labels)
# }
# # Traceplots for pi_raw
# mcmc_trace(rearrest_fit$draws(variables = c(
#     "pi_raw", "log_lambda"
# )))

# Copy the fitted rates into setup.R (replaces manual copying):
#   Rscript scripts/model/update_rate_params.R
