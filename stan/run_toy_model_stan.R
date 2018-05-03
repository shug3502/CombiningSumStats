library(rstan)

n <- 100
theta_true <- 10
data <- 10*runif(n)

estimates <- stan(file = 'toy_model.stan',
                  data = list (
                    n=n,
                    y=data
                  ),
                  seed = 42,
                  chains = 1,
                  iter = 2000,
                  warmup = 1000,
                  control = list(adapt_delta = 0.99)
)

e <- extract(estimates,pars=c("log_theta","theta","y_sim","lp__"),permuted=TRUE)

write.csv(e$log_theta,file='toy_model_stan.csv',row.names=FALSE)

