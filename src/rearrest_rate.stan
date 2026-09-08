functions {
  real softplus(real x) {
    return log1p_exp(x);  // smooth, always positive
  }
}

data {
  int<lower=1> n_age;
  int<lower=1> N;                    // total number of ppl
  matrix<lower=0>[N, n_age] t;       // durations (years)
  array[N, 2] int<lower=1, upper=n_age> age_r; // corresponding age group index range
  array[N] int event;                // 1 = re-arrested at end of segment
  
  // equilibrium calibration
  array[n_age] int J;                // observed number of ppl in Jail in 2015
  int f_obs;                         // observed number of years with first arrest distribution
  array[f_obs, n_age] int F0;        // observed first offender into Jail by year (2014-2016)
  vector[n_age] mu;                  // death rate
  vector[n_age] release;             // release rate
  array[n_age] int tot_in;           // age-specific inflow population size
}

transformed data {
  array[n_age] int F0_total;
  for(k in 1:n_age){
    F0_total[k] = 0;
    for(i in 1:f_obs){
      F0_total[k] += F0[i,k];
    }
    F0_total[k] = F0_total[k] %/% f_obs;  // average over years
  }
}

parameters {
  real<lower=0> l_sigma;
  vector[n_age] log_lambda;
  real pi_raw; 
  
  //equilibrium calibration
  real<lower=0> j_total;  // total jail population
  simplex[n_age] j_props; // age distribution of jail population
  vector<lower=0>[n_age] c; // age-specific scaling for first-arrest rates
}

transformed parameters {
  real<lower = 0, upper = 1> pi = inv_logit(pi_raw); 
  vector<lower=0>[n_age] lambda = exp(log_lambda);
  vector<lower=0>[N] hazard;
  // vector<lower=0, upper=1>[N] p_arrest;

  // equilibrium calibration
  vector[n_age] r_pop;
  vector[n_age] f_pop;
  vector[n_age] i_pop;
  vector[n_age] j_pop = j_total * to_vector(j_props);
  
  for(i in 1:N){
    hazard[i] = dot_product(t[i, age_r[i,1]:age_r[i,2]], lambda[age_r[i,1]:age_r[i,2]]);
    // p_arrest[i] = pi * (1 - exp(-hazard[i]));
  }

  // equilibrium calibration
  real a = 0.2; 
  r_pop[1] = softplus(pi * release[1] * j_pop[1]  / (lambda[1] + mu[1] + a));
  f_pop[1] = softplus((release[1] + mu[1] + a) * j_pop[1] - lambda[1] * r_pop[1]);
  for(i in 2:(n_age-1)){
    r_pop[i] = softplus((pi * release[i] * j_pop[i] + a * r_pop[i - 1])) / (lambda[i] + mu[i] + a);
    f_pop[i] = softplus((release[i] + mu[i] + a) * j_pop[i] - lambda[i] * r_pop[i] - a * j_pop[i - 1]);
  }
  r_pop[n_age] = softplus((pi * release[n_age] * j_pop[n_age] + a * r_pop[n_age - 1])) / (lambda[n_age] + mu[n_age]);
  f_pop[n_age] = softplus((release[n_age] + mu[n_age]) * j_pop[n_age] - lambda[n_age] * r_pop[n_age] - a * j_pop[n_age - 1]);
  for(i in 1:n_age){
    i_pop[i] = softplus(f_pop[i] / c[i] / lambda[i]);
  }

  real<lower=0> f_total = sum(f_pop);
  simplex[n_age] f_props = f_pop / f_total;
}

model {
  pi_raw ~ normal(0, 1);
  l_sigma ~ exponential(1);
  log_lambda[1] ~ normal(0, 2);
  for(i in 2:n_age){
    log_lambda[i] - log_lambda[i - 1] ~ normal(0, l_sigma);
  }
  
  // Log-likelihood function
  for (j in 1:N) {
      target += (log(pi) + log(lambda[age_r[j,2]]) - hazard[j]) * event[j] + log(1 - pi + pi * exp(-hazard[j])) * (1 - event[j]);
  }
  // event ~ bernoulli(p_arrest);
  
  // equilibrium part
  sum(J) ~ poisson(sum(j_pop));
  J ~ multinomial(j_props);
  F0_total ~ multinomial(f_props);

  c ~ normal(1, 1);
  tot_in ~ poisson(i_pop);
}
