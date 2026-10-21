// =============================================================================
// rearrest_rate.stan — cure-type survival model of re-arrest after release (probability of
// re-arrest and age-specific hazard; used by preprocess.R).
// =============================================================================

data {
  int<lower=1> n_age;
  int<lower=1> N;                    // total number of ppl
  matrix<lower=0>[N, n_age] t;       // durations (years)
  array[N, 2] int<lower=1, upper=n_age> age_r; // corresponding age group index range
  array[N] int event;                // 1 = re-arrested at end of segment
}

parameters {
  real<lower=0> l_sigma;
  vector[n_age] log_lambda;
  real pi_raw; 
}

transformed parameters {
  real<lower = 0, upper = 1> pi = inv_logit(pi_raw); 
  vector<lower=0>[n_age] lambda = exp(log_lambda);
  vector<lower=0>[N] hazard;

  for(i in 1:N){
    hazard[i] = dot_product(t[i, age_r[i,1]:age_r[i,2]], lambda[age_r[i,1]:age_r[i,2]]);
  }
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
}
