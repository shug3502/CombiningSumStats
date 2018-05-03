data {
  int<lower=1> n;       // Number of obs
  real<lower=0> y[n];    // observations
}

parameters {
  real<lower=(0),upper=2> log_theta;
}

transformed parameters {
  real<lower=10^(0),upper=10^(2)> theta;
  theta = 10^log_theta;
}

model {
  // Priors
  log_theta ~ uniform(0, 2);
  // Likelihood
  y ~ uniform(0,theta);
}

generated quantities {
  real y_sim[n];
  for (i in 1:n){
    y_sim[i] = uniform_rng(0,theta);
  }
}
