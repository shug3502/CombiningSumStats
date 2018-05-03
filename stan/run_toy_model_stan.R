library(rstan)
options(mc.cores = parallel::detectCores())

n <- 10
theta_true <- 10
data <- theta_true*runif(n)

estimates <- stan(file = 'toy_model.stan',
                  data = list (
                    n=n,
                    y=data
                  ),
                  seed = 42,
                  chains = 4,
                  iter = 2000,
                  warmup = 1000,
                  control = list(adapt_delta = 0.99)
)

e <- extract(estimates,pars=c("log_theta","theta","y_sim","lp__"),permuted=TRUE)
posterior <- as.array(estimates)
write.csv(e$log_theta,file='toy_model_stan.csv',row.names=FALSE)

summary(estimates,pars=c('log_theta','theta','lp__'))
library(bayesplot)
mcmc_areas(posterior,pars='log_theta')
mcmc_hist(posterior,pars='log_theta')
mcmc_dens_overlay(posterior,pars='log_theta')
mcmc_pairs(posterior,pars=c('log_theta','theta','lp__'))
