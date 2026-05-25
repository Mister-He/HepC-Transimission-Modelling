# jail release and re-arrest rates
# utility functions---------------
split_by_age <- function(id, t_tot, age_s, obs_yn) {
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


estim <- function(dat) {
  lapply(1:ncol(dat), function(i) {
    c(mean(dat[, i]), quantile(dat[, i], c(0.025, 0.975)))
  }) %>%
    do.call("rbind", .) %>%
    as.matrix()
}

set_panels <- function(nrow = 2, ncol = 1, margin0 = c(1, 0.5, 1, 1), margin1 = c(2, 4.5, 1, 0), byrow = T) {
  if (byrow) {
    par(mfrow = c(nrow, ncol), mar = margin1, oma = margin0)
  } else {
    par(mfcol = c(nrow, ncol), mar = margin1, oma = margin0)
  }
}

add_polygon <- function(dat, x_idx = NA, clr = "indianred2") {
  if (mean(is.na(x_idx) > 0)) x_idx <- 1:nrow(dat)
  polygon(c(x_idx, rev(x_idx)),
    c(dat[, 1], rev(dat[, 2])),
    col = scales::alpha(clr, 0.3), border = NA
  )
}


library(dplyr)
library(cmdstanr)
library(survival)
library(lubridate)

jail_raw <- read.csv("~/Downloads/moh.csv")
jail <- jail_raw %>%
  mutate(
    adm_dt = as.Date(ADM_DT, tryFormats = "%m/%d/%Y"),
    rel_dt = as.Date(PERM_REL_DATE, tryFormats = "%m/%d/%Y"),
    hepC = as.numeric(Result.from.screening.exercise..Dec.2014.Feb2016. %in% c("Borderline Reactive", "Reactive"))
  ) %>%
  mutate(
    end_dt = ifelse(is.na(PERM_REL_DATE), "2017-12-31", as.character(rel_dt)) %>% as.Date(),
    injailatend = is.na(PERM_REL_DATE) %>% as.numeric(),
    adm_yr = year(adm_dt)
  ) %>%
  mutate(
    t_injail = as.numeric(end_dt - adm_dt) / 365,
    adm_age = adm_yr - YEAR_OF_BIRTH
  )
jail$injailatend[jail$rel_dt > as.Date("2017-12-31")] <- 1
ageint <- 3:11 * 5
n_age <- length(ageint)


age_labels = c("15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55+")

jail_consumption = jail#jail[jail$CONSUMPTION ==1,]
counts_by_age_year = lapply(1:nrow(jail_consumption), function(i){
  row = jail_consumption[i,]
  
  # Calendar years this person spans (from admission to end date)
  start_yr = year(row$adm_dt)
  end_yr   = year(row$end_dt)
  
  lapply(start_yr:end_yr, function(yr){
    # Age during this calendar year
    age_in_yr = yr - row$YEAR_OF_BIRTH
    agegp = cut(age_in_yr, breaks = c(ageint, Inf), 
                labels = age_labels, right = FALSE)
    
    data.frame(
      id       = row$IDENTIFER,
      year     = yr,
      age_grp  = as.character(agegp),
      hepC     = row$hepC
    )
  }) %>% do.call('rbind', .)
}) %>% do.call('rbind', .) %>%
  filter(!is.na(age_grp))

# Deduplicate — one row per person per year (in case of multiple admissions)
counts_by_age_year = counts_by_age_year %>%
  group_by(id, year) %>%
  slice(1) %>%
  ungroup()

# Summarise
summary_table = counts_by_age_year %>%
  group_by(year, age_grp) %>%
  summarise(
    hcv_cases  = sum(hepC, na.rm = TRUE),
    total_pop  = n(),
    .groups = "drop"
  ) %>%
  arrange(year, age_grp)

print(summary_table)

part_jail = jail %>% 
  filter(injailatend == 1) %>%
  group_by(AgeGroup_alt) %>%
  summarise(
    hcv_cases = sum(hepC, na.rm = TRUE),
    total_pop = n(),
    .groups = "drop"
  ) %>%
  mutate(prev = hcv_cases / total_pop * 100) %>%
  arrange(AgeGroup_alt)  

sum(part_jail$hcv_cases) / sum(part_jail$total_pop) * 100
