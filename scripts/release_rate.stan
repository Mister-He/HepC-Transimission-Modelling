data {
  int<lower=1> n_age;
  int<lower=1> N_seg;          // total number of segments
  vector<lower=0>[N_seg] t;    // segment durations (years)
  array[N_seg] int<lower=1, upper=n_age> age; // age group during segment
  array[N_seg] int event;      // 1 = released at end of segment
}

parameters {
  real<lower=0> l_sigma;
  vector[n_age] log_lambda;
}

transformed parameters {
  vector<lower=0>[n_age] lambda = exp(log_lambda);
}

model {
  l_sigma ~ exponential(1);
  log_lambda[1] ~ normal(0, 10);
  for(i in 2:n_age){
    log_lambda[i] - log_lambda[i - 1] ~ normal(0, l_sigma);
  }

  for (j in 1:N_seg) {
      target += log(lambda[age[j]]) * event[j]  - lambda[age[j]] * t[j];
  }
}
