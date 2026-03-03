#jail release and re-arrest rates
#utility functions---------------
split_by_age = function(id, t_tot, age_s, obs_yn){
  out = c()
  idx = sum(ageint<=age_s)
  while(t_tot>0 && idx <= n_age){
    if(idx == n_age){
      t0 = t_tot; event = obs_yn #observed
    }else{
      t0 = min(t_tot, ageint[idx+1] - max(age_s, ageint[idx]))
      event = ifelse(abs(t_tot-t0)<1e-3, obs_yn, 0)
    }
    out = rbind(out, data.frame(id = id, age = idx, t = t0, event = event))
    idx = idx + 1
    t_tot = t_tot - t0
  }
  out
}


estim = function(dat){
  lapply(1:ncol(dat), function(i){
    c(mean(dat[,i]), quantile(dat[,i],c(0.025,0.975)))
  })%>%do.call('rbind',.)%>%as.matrix
}

set_panels=function(nrow=2, ncol=1, margin0=c(1, 0.5, 1, 1), margin1=c(2, 4.5, 1, 0), byrow=T){
  if(byrow){
    par(mfrow=c(nrow,ncol), mar = margin1, oma = margin0)
  }else{
    par(mfcol=c(nrow,ncol), mar = margin1, oma = margin0)
  }
}

add_polygon=function(dat, x_idx = NA, clr='indianred2'){
  if(mean(is.na(x_idx)>0)) x_idx = 1:nrow(dat)
  polygon(c(x_idx, rev(x_idx)), 
          c(dat[,1], rev(dat[,2])), 
          col=scales::alpha(clr,0.3),border=NA)
}


library(dplyr)
library(cmdstanr)
library(survival)
library(lubridate)

jail_raw = read.csv('~/Downloads/moh.csv')
jail = jail_raw %>%
  mutate(adm_dt = as.Date(ADM_DT,tryFormats = '%m/%d/%Y'), 
         rel_dt = as.Date(PERM_REL_DATE,tryFormats = '%m/%d/%Y'), 
         hepC = as.numeric(Result.from.screening.exercise..Dec.2014.Feb2016. %in% c('Borderline Reactive', 'Reactive')))%>%
  mutate(end_dt = ifelse(is.na(PERM_REL_DATE), '2017-12-31', as.character(rel_dt))%>%as.Date(), 
         injailatend = is.na(PERM_REL_DATE)%>%as.numeric, 
         adm_yr = year(adm_dt))%>%
  mutate(t_injail = as.numeric(end_dt - adm_dt)/365, 
         adm_age = adm_yr - YEAR_OF_BIRTH)
jail$injailatend[jail$rel_dt>as.Date('2017-12-31')] = 1
ageint = 3:11*5; n_age = length(ageint)



#release model------------
stan_code = '
data {
  int n_age;
  int<lower=1> N_seg;                      // total number of segments
  vector<lower=0>[N_seg] t;                // segment durations (years)
  array[N_seg] int age;  // age group during segment
  array[N_seg] int event;      // 1 = released at end of segment
}

parameters {
  real l_mu;
  real<lower=0> l_sigma;
  vector[n_age] log_lambda;
}

transformed parameters {
  vector<lower=0>[n_age] lambda = exp(log_lambda);
}

model {
  l_sigma ~ normal(0, 1);
  log_lambda[1] ~ normal(0, 10);
  for(i in 2:n_age){
    log_lambda[i] - log_lambda[i - 1] ~ normal(0, l_sigma);
  }

  for (j in 1:N_seg) {
      target += log(lambda[age[j]]) * event[j]  - lambda[age[j]] * t[j];
  }
}
'

release_dat = lapply(1:nrow(jail), function(i){
  if(jail$t_injail[i]>0){
    split_by_age(jail$IDENTIFER[i], jail$t_injail[i], jail$adm_age[i],  1 - jail$injailatend[i])
  }
})%>%do.call('rbind',.)

data_list = list(
  n_age = n_age,
  N_seg = nrow(release_dat),
  t = release_dat$t,
  age = release_dat$age,
  event = release_dat$event
)

n.sim=5e3
model=cmdstan_model(stan_file = write_stan_file(stan_code))
set.seed(123); n.chains=3; fit=model$sample(data = data_list, chains = n.chains,
                                            iter_sampling = n.sim, iter_warmup = n.sim, 
                                            parallel_chains = n.chains, save_cmdstan_config = TRUE)
out_release=fit$draws(format = 'matrix')
plot.ts(out_release[,paste0('lambda[', 1:n_age, ']')])
lambda1 = out_release[,paste0('lambda[', 1:n_age, ']')]%>%estim%>%round(3)


#survival fit plot
# set_panels(3,3)
# for(i in 1:n_age){
#   km = survfit(Surv(t, event) ~ 1, data = release_dat%>%filter(age == i))
#   plot(km, conf.int = T, ylab = paste0('Group ',i))
#   t_grid = seq(0, 5, 0.05)
#   lines(t_grid, exp(-lambda[i,1]*t_grid), col='red')
#   add_polygon(cbind(exp(-lambda1[i,2]*t_grid), exp(-lambda[i,3]*t_grid)), t_grid)
# }



#re-arrest model-----------
stan_rearrest = '
data {
  int n_age;
  int<lower=1> N;                      // total number of ppl
  matrix<lower=0>[N, n_age] t;                // durations (years)
  array[N, 2] int age_r;             // corresponding age group index range
  array[N] int event_c; // 1 = re-arrested at end of segment
  
  //uncensored part
  /*int<lower=1> N1; //number of segments for uncensored part
  vector<lower=0>[N1] t1; //duration (years)
  array[N1] int age1;  // age group during segment
  array[N1] int event_o;*/      // 1 = re-arrested at end of segment
  
  //equilibrium calibration
  array[n_age] int J; // observed number of ppl in Jail in 2015
  int f_obs; //observed number of years with first arrest distribution
  array[f_obs, n_age] int F0; // observed first offender into Jail by year (2014-2016)
  vector[n_age] mu; //death rate
  vector[n_age] lambda1; //release rate
  int tot_in; //total inflow
}

parameters {
  real<lower=0> l_sigma;
  vector[n_age] log_lambda;
  real pi_raw; 
  //real<lower = 0, upper = 1> pi; 
  //equilibrium calibration
  vector<lower=0>[n_age] j_pop;
  real<lower=0> c; //scaling for first-arrest rates
}

transformed parameters {
  real<lower = 0, upper = 1> pi = inv_logit(pi_raw); 
  vector<lower=0>[n_age] lambda = exp(log_lambda);
  vector<lower=0>[N] hazard;
  vector<lower=0>[N] p_arrest;
  //equilibrium calibration
  vector[n_age] r_pop;
  vector[n_age] f_pop;
  vector[n_age] i_pop;
  real i_sum;
  
  for(i in 1:N){
    hazard[i] = dot_product(t[i, age_r[i,1]:age_r[i,2]], lambda[age_r[i,1]:age_r[i,2]]);
    p_arrest[i] = pi * (1 - exp(-hazard[i]));
  }
  
  //equilibrium calibration
  real a = 0.2; 
  r_pop[1] = fmax(1e-6, pi * lambda1[1] * j_pop[1]  / (lambda[1] + mu[1] + a));
  f_pop[1] = fmax(1e-6, (lambda1[1] + mu[1] + a) * j_pop[1] - lambda[1] * r_pop[1]);
  for(i in 2:(n_age-1)){
    r_pop[i] = fmax(1e-6, (pi * lambda1[i] * j_pop[i] + a * r_pop[i - 1])) / (lambda[i] + mu[i] + a);
    f_pop[i] = fmax(1e-6, (lambda1[i] + mu[i] + a) * j_pop[i] - lambda[i] * r_pop[i] - a * j_pop[i - 1]);
  }
  r_pop[n_age] = fmax(1e-6, (pi * lambda1[n_age] * j_pop[n_age] + a * r_pop[n_age - 1])) / (lambda[n_age] + mu[n_age]);
  f_pop[n_age] = fmax(1e-6, (lambda1[n_age] + mu[n_age]) * j_pop[n_age] - lambda[n_age] * r_pop[n_age] - a * j_pop[n_age - 1]);
  for(i in 1:n_age){
    i_pop[i] = f_pop[i] / fmax(1e-6, c * lambda[i]);
  }
  i_sum = fmax(1e-6, dot_product(to_row_vector(mu), i_pop)) + sum(f_pop);
}

model {
  pi_raw ~ normal(0, 1);
  //pi ~ beta(2, 2);
  l_sigma ~ normal(0, 0.1);
  log_lambda[1] ~ normal(0, 10);
  for(i in 2:n_age){
    log_lambda[i] - log_lambda[i - 1] ~ normal(0, l_sigma);
  }
  
  //censored part
  event_c ~ binomial(1, p_arrest);
  
  //uncensored part
  /*for (j in 1:N1) {
      target += log(lambda[age1[j]]) * event_o[j]  - lambda[age1[j]] * t1[j];
  }*/
  
  //equilibrium part
  sum(J) ~ poisson(sum(j_pop));
  J ~ multinomial(j_pop/sum(j_pop));
  for(i in 1:f_obs){
    F0[i,] ~ multinomial(f_pop/sum(f_pop));
  }
  c ~ normal(1, 1);
  tot_in ~ poisson(i_sum);
}
'

#rearrest data derivation
rearrest_dat0 = lapply(jail$IDENTIFER, function(id){
  dat0 = jail%>%
    filter(IDENTIFER == id)%>%
    arrange(adm_dt)
  n.record = nrow(dat0)
  #with explicit rearrest interval
  if(n.record>1){
    out1 = lapply(2:n.record, function(i){
      delta = as.numeric(dat0$adm_dt[i] - dat0$rel_dt[i-1])/365
      age = dat0$adm_age[i-1] + dat0$t_injail[i-1]
      data.frame(id = id, idx = i-1, t = delta, age_rel = age, censor = 0, rearrest = 1)
    })%>%do.call('rbind',.)
  }
  dat0 = dat0 %>%slice(n.record)
  #release before 2015 --> assuming rearrest before 2016 as they need to have test results
  #last release record after 2015-01-01 --> assuming no re-arrest thereafter
  last_obs = F
  if(dat0$injailatend == 0){
    dur = release = 0
    if(dat0$rel_dt <= as.Date('2014-12-31')){
      #ppl released before data collection
      end = as.Date('2015-12-31')
      dur = as.numeric(end - dat0$rel_dt)/365 #longest possible in community time
      release = 0
    }else if(dat0$end_date < as.Date('2017-12-31') && dat0$t_injail > 0){
      end = as.Date('2017-12-31')
      dur = as.numeric(end - dat0$rel_dt)/365 # censored time (during this period no re-arrest event)
      release = 1
    }
    if(dur+release > 0){
      out0 = data.frame(id = id, idx = n.record, t = dur, 
                       age_rel = dat0$adm_age + dat0$t_injail, censor = 1, rearrest = 1 - release)
      # dat1 = split_by_age(id, dur, dat0$adm_age + dat0$t_injail,  1 - release)
      # t_all = rep(0, n_age)
      # n1 = nrow(dat1)
      # t_all[dat1$age%>%round] = dat1$t; age_range = range(dat1$age)
      # c(id, t_all, age_range, 1 - release)
      last_obs = T
    }
  }
  if(n.record>1 & last_obs){
    rbind(out1, out0)
  }else if(n.record>1){
    out1
  }else if(last_obs){
    out0
  }
})%>%do.call('rbind', .)

#split by age group
rearrest_censor = lapply(which(rearrest_dat0$censor==1), function(i){
  dat0 = rearrest_dat0%>%slice(i)
  dat1 = split_by_age(dat0$id, dat0$t, dat0$age_rel,  1 - dat0$censor)
  t_all = rep(0, n_age)
  n1 = nrow(dat1)
  t_all[dat1$age%>%round] = dat1$t; age_range = range(dat1$age)
  c(id, t_all, age_range, dat0$rearrest)
})%>%do.call('rbind', .)%>%as.matrix

rearrest_obs = lapply(which(rearrest_dat0$censor==0), function(i){
  dat0 = rearrest_dat0%>%slice(i)
  split_by_age(dat0$id, dat0$t, dat0$age_rel,  dat0$rearrest)
})%>%do.call('rbind', .)

#other data (for equilibrium calibration, not necessary atm)
#jail age structure in 2015
j_pop = jail%>%
  filter(adm_dt<=as.Date('2015-12-31'), 
         #rel_dt>=as.Date('2014-12-31'),
         t_injail>0)%>%
  group_by(IDENTIFER)%>%
  slice(1)%>%
  ungroup()%>%
  mutate(current_age = 2015 - YEAR_OF_BIRTH)%>%
  mutate(agegp = cut(current_age, breaks = c(ageint, Inf), labels = seq_along(ageint)))%>%
  group_by(agegp)%>%
  summarise(n = n())%>%pull(n)
#overall first-arrest age composition
first_jail = jail%>%
  filter(t_injail>0)%>%
  arrange(IDENTIFER, adm_dt)%>%
  group_by(IDENTIFER)%>%
  slice(1)%>%
  ungroup()
fadm_dist = first_jail%>%
  mutate(adm_age = year(adm_dt) - YEAR_OF_BIRTH)%>%
  mutate(agegp = cut(adm_age, breaks = c(ageint, Inf), labels = seq_along(ageint)))%>%
  group_by(agegp)%>%
  summarise(n = n())%>%pull(n)
#first-arrest age composition by admission year 
f0_dist = lapply(2014:2016, function(t){
  first_jail%>%
    filter(t_injail>0, year(adm_dt)==t)%>%
    mutate(adm_age = year(adm_dt) - YEAR_OF_BIRTH)%>%
    mutate(agegp = cut(adm_age, breaks = c(ageint, Inf), labels = seq_along(ageint)))%>%
    group_by(agegp)%>%
    summarise(n = n())%>%pull(n)
})%>%do.call('rbind', .)%>%as.matrix
f_obs = nrow(f0_dist)
#death rate
mu = c(0.001267,0.0003,0.0003,0.0004,0.0005,0.0007,0.0014,0.0023,0.0161)
#number of new PWID per year
tot_in = 1861

#data for stan code
data_list1 = list(
  n_age = n_age,
  #data for equilibrium derivation
  lambda1 = lambda1[,1],
  F0 = f0_dist, f_obs = f_obs,
  J = j_pop,
  mu = mu,
  tot_in = tot_in,
  #data for survival analysis
  #uncensored part
  # N1 = nrow(rearrest_obs),
  # t1 = rearrest_obs$t,
  # age1 = rearrest_obs$age,
  # event_o = rearrest_obs$event,
  #censored part (only know y/n for a specific interval)
  N = nrow(rearrest_censor), 
  t = rearrest_censor[,1:n_age + 1],
  age_r = rearrest_censor[,1:2 + n_age + 1]%>%round(),
  event_c = rearrest_censor[,n_age + 4]
)



n.sim=5e3
model=cmdstan_model(stan_file = write_stan_file(stan_rearrest))
set.seed(123); n.chains=3; fit=model$sample(data = data_list1, chains = n.chains,
                                            iter_sampling = n.sim, iter_warmup = n.sim,
                                            parallel_chains = n.chains, save_cmdstan_config = TRUE)
out_re=fit$draws(format = 'matrix')
plot.ts(out_re[,paste0('lambda[', 1:n_age, ']')])
lambda2 = out_re[,paste0('lambda[', 1:n_age, ']')]%>%estim

plot.ts(out_re[,'pi'])
# plot.ts(out_re[,'c'])
pi = out_re[,paste0('pi')]%>%estim
plot.ts(out_re[,paste0('j_pop[', 1:n_age, ']')])
out_re[,paste0('j_pop[', 1:n_age, ']')]%>%estim
out_re[,paste0('f_pop[', 1:n_age, ']')]%>%estim


prob = out_re[,paste0('p_arrest[', 1:data_list1$N, ']')]%>%estim
mean((prob[,1]>0.5) == data_list1$event)


save(out_release, out_re, j_pop, fadm_dist, f0_dist, file = 'jail_sum.RData')


#equilibrium solver----------
# a = 1/5 #aging
# solve_R <- function(J, lambda1, lambda2, mu, pi, a) {
#   R <- rep(0, n_age)
#   for (i in 1:n_age) {
#     R[i] <- (pi * lambda1[i] * J[i] + ifelse(i>1, a * R[i - 1], 0)) / (lambda2[i] + mu[i] + a * (i<n_age))
#   }
#   R
# }
# 
# F_from_R = function(R, J, lambda1, lambda2, mu, a){
#   lapply(1:n_age, function(i){
#     (lambda1[i] + mu[i] + a * (i<n_age)) * J[i] - lambda2[i] * R[i] - ifelse(i>1, a * J[i-1], 0)
#   })%>%unlist
# }
# 
# 
# total_inflow_from_c <- function(c, F0, lambda, mu, a) {
#   #c: scaling factor --  lambda0 = lambda2 * c
#   I <- F0 / (c * lambda)
#   
#   lapply(1:n_age, function(i){
#     (mu[i] + a) * I[i] + F0[i] - ifelse(i>1, a * I[i - 1], 0)
#   })%>%unlist%>%sum
#   
# }
# 
# pi=0.3
# 
# R <- solve_R(j_pop, lambda1[,1], lambda2[,1], mu, pi, a)
# F_from_R(R, j_pop, lambda1[,1], lambda2[,1], mu, a)
# 
# root_fn <- function(c) {
#   total_inflow_from_c(c, f2015_dist, lambda2[,1], mu, a) - tot_in
# }
# 
# uniroot(root_fn, lower = 1e-6, upper = 1e3)$root








