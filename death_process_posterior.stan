  data {
    int<lower=1> N; //initial condition
    int<lower=1> n; //number of time pts
    int<lower=1> repeats; //number of repeated observations
    int y[repeats,n]; //observations
    real z[repeats]; //extra noise
    real t[n]; //observation times
  }
  parameters {
    real<lower=-4, upper=4> log_k; //decay param
    real<lower=-4, upper=4> log_sigma; //noise
  }
  transformed parameters {
    real<lower=10^(-4),upper=10^4> k;
    real<lower=10^(-4),upper=10^4> sigma;
    k = 10^log_k;
    sigma = 10^log_sigma;
  }
  model {
//    priors
    log_k ~ uniform(-4,4);
    log_sigma ~ uniform(-4,4);
    for (i in 1:repeats){
      y[i,1] ~ binomial(N,exp(-k*t[1]));
      for (j in 2:n){
        y[i,j] ~ binomial(y[i,j-1],exp(-k*(t[j]-t[j-1])));
      }
      z[i] ~ normal(0,sigma); //sigma not sigma^2 i think 
    }
  }
  generated quantities {
    //currently only simulate a single example, not repeats
    int y_sim[n];
    real z_sim;
    y_sim[1] = binomial_rng(N,exp(-k*t[1]));
    for (j in 2:n){
      y_sim[j] = binomial_rng(y_sim[j-1],exp(-k*(t[j]-t[j-1])));
    }
    z_sim = normal_rng(0,sigma); 
  }


